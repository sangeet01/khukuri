"""
shard.py
--------
One shard = one memory-mapped float32 matrix on disk.
A shard is append-only while hot, then sealed (frozen) at run end.
Frozen shards are never modified.

Shard file layout:
  <name>.mm   -- float32 memmap, shape (capacity, dim)
  <name>.meta -- JSON: {"n_entries": int, "dim": int, "ids": [...]}
"""

import json
import numpy as np
from pathlib import Path
from typing import Optional


SHARD_CAPACITY = 50_000   # max entries per shard (resize if needed)


class MemoryShard:
    """
    A single on-disk shard.
    Hot shards accept writes. Frozen shards are read-only.
    """

    def __init__(self, path: Path, dim: int, frozen: bool = False):
        self.path   = path
        self.dim    = dim
        self.frozen = frozen
        self._mm: Optional[np.memmap] = None
        self._ids: list[str] = []
        self._n = 0

        meta_path = path.with_suffix(".meta")
        if path.exists() and meta_path.exists():
            self._load()
        elif not frozen:
            self._create()

    # -- write -----------------------------------------------------------------

    def add(self, entry_id: str, vec: np.ndarray) -> int:
        """
        Append a vector. Returns its row index.
        Raises if shard is frozen.
        """
        if self.frozen:
            raise RuntimeError(f"Shard {self.path.name} is frozen.")
        if self._n >= SHARD_CAPACITY:
            raise RuntimeError(f"Shard {self.path.name} full ({SHARD_CAPACITY} entries).")

        row = self._n
        self._mm[row] = vec.astype(np.float32)
        self._mm.flush()
        self._ids.append(entry_id)
        self._n += 1
        self._save_meta()
        return row

    def seal(self) -> None:
        """Freeze this shard. Called at end of run."""
        self.frozen = True
        self._save_meta()

    # -- read ------------------------------------------------------------------

    def query(
        self,
        query_vec: np.ndarray,
        candidate_rows: Optional[list[int]] = None,
        k: int = 10,
    ) -> list[tuple[str, float]]:
        """
        Cosine similarity search within this shard.
        candidate_rows: if given, only score those rows (post metadata filter).
        Returns list of (entry_id, score) sorted descending.
        """
        if self._n == 0:
            return []

        vecs = np.array(self._mm[:self._n], dtype=np.float32)

        if candidate_rows is not None:
            # filter to valid rows only
            rows = [r for r in candidate_rows if 0 <= r < self._n]
            if not rows:
                return []
            row_arr = np.array(rows, dtype=np.int32)
            vecs_sub = vecs[row_arr]
            scores = vecs_sub @ query_vec
            # map back to entry ids
            top_k = min(k, len(rows))
            if top_k == 0:
                return []
            idx = np.argpartition(scores, -top_k)[-top_k:]
            idx = idx[np.argsort(scores[idx])[::-1]]
            return [(self._ids[row_arr[i]], float(scores[i])) for i in idx]
        else:
            scores = vecs @ query_vec   # (n,)
            top_k = min(k, self._n)
            idx = np.argpartition(scores, -top_k)[-top_k:]
            idx = idx[np.argsort(scores[idx])[::-1]]
            return [(self._ids[i], float(scores[i])) for i in idx]

    def get_vector(self, row: int) -> np.ndarray:
        return np.array(self._mm[row], dtype=np.float32)

    @property
    def n_entries(self) -> int:
        return self._n

    @property
    def name(self) -> str:
        return self.path.stem

    # -- internal --------------------------------------------------------------

    def _create(self) -> None:
        self._mm = np.memmap(
            str(self.path), dtype=np.float32, mode="w+",
            shape=(SHARD_CAPACITY, self.dim)
        )
        self._ids = []
        self._n   = 0
        self._save_meta()

    def _load(self) -> None:
        meta = json.loads(self.path.with_suffix(".meta").read_text())
        self._n    = meta["n_entries"]
        self._dim  = meta["dim"]
        self._ids  = meta["ids"]
        mode = "r" if self.frozen else "r+"
        self._mm   = np.memmap(
            str(self.path), dtype=np.float32, mode=mode,
            shape=(SHARD_CAPACITY, self._dim)
        )

    def _save_meta(self) -> None:
        meta = {
            "n_entries": self._n,
            "dim":       self.dim,
            "frozen":    self.frozen,
            "ids":       self._ids,
        }
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps(meta))
        tmp.replace(self.path.with_suffix(".meta"))
