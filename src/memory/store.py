"""
store.py
--------
KhukuriMemory -- the unified memory store for Khukuri.

One KhukuriMemory instance per named store (e.g. "threats", "compounds", "literature").
Each store lives in its own directory under the project root.

Directory layout:
    <project_root>/<store_name>/
        manifest.json       <- metadata + salience for all entries (all runs)
        hot.mm              <- current run's vectors (being written to)
        hot.meta
        frozen/
            run_001.mm      <- sealed at end of run 1
            run_001.meta
            run_002.mm
            ...

Autonomous loop integration:
    On each iteration:  mem.add(...)  -> writes to hot shard
    On iteration end:   mem.commit_run()  -> seals hot -> opens new hot
    Loop runs forever:  memory just grows, query always spans all shards
"""

import time
import shutil
from pathlib import Path
from typing import Any, Optional

import numpy as np

from .encoder import NumenEncoder
from .manifest import Manifest
from .shard import MemoryShard
from .query import QueryEngine, QueryResult


class KhukuriMemory:
    """
    Tiered, persistent, salience-weighted memory store.

    Usage:
        mem = KhukuriMemory("threats", project_root="/path/to/project")

        # add (writes to hot shard only)
        mem.add("threat_042", text="mecA horizontal transfer", meta={"type": "HGT", "score": 0.9})

        # query (fans across all shards)
        results = mem.query("MRSA resistance beta-lactam", k=10, filters={"type": "HGT"})

        # end of autonomous iteration
        mem.commit_run()

        # feedback from the loop
        mem.reinforce("threat_042", positive=True)
    """

    def __init__(
        self,
        store_name: str,
        project_root: str | Path,
        dim: int = 4096,
        encoder: Optional[NumenEncoder] = None,
        salience_weight: float = 0.3,
    ):
        self.store_name = store_name
        self.dim        = dim

        # directory for this store
        self.root = Path(project_root) / store_name
        self.root.mkdir(parents=True, exist_ok=True)
        (self.root / "frozen").mkdir(exist_ok=True)

        self.encoder  = encoder or NumenEncoder(dim=dim)
        self.manifest = Manifest(self.root / "manifest.json")

        # load or create hot shard
        self._hot: MemoryShard = self._open_hot()

        # load frozen shards (read-only)
        self._frozen: list[MemoryShard] = self._load_frozen()

        # query engine
        self._qe = QueryEngine(self.manifest, salience_weight)

        # run counter
        self._run_id = self._infer_run_id()

    # -- write -----------------------------------------------------------------

    def add(
        self,
        entry_id: str,
        text: str,
        meta: Optional[dict[str, Any]] = None,
    ) -> None:
        """
        Encode text and add to hot shard.
        Silently skips if entry_id already exists (idempotent).
        """
        if self.manifest.get(entry_id) is not None:
            return   # already indexed

        vec = self.encoder.encode(text)
        row = self._hot.add(entry_id, vec)
        self.manifest.register(
            entry_id  = entry_id,
            shard_name= self._hot.name,
            row_idx   = row,
            text_key  = text[:120],   # first 120 chars as key preview
            meta      = meta or {},
        )

    def add_batch(
        self,
        entries: list[dict],
    ) -> None:
        """
        Batch add. Each entry: {"id": str, "text": str, "meta": dict}
        Skips duplicates.
        """
        new = [e for e in entries if self.manifest.get(e["id"]) is None]
        if not new:
            return
        texts = [e["text"] for e in new]
        vecs  = self.encoder.encode_batch(texts)
        for i, entry in enumerate(new):
            row = self._hot.add(entry["id"], vecs[i])
            self.manifest.register(
                entry_id  = entry["id"],
                shard_name= self._hot.name,
                row_idx   = row,
                text_key  = entry["text"][:120],
                meta      = entry.get("meta", {}),
            )

    def reinforce(self, entry_id: str, positive: bool = True) -> None:
        """
        Feedback from the autonomous loop.
        Call after you know whether a recalled entry was useful.
        """
        self.manifest.reinforce(entry_id, positive)

    # -- read ------------------------------------------------------------------

    def query(
        self,
        text: str,
        k: int = 10,
        filters: Optional[dict[str, Any]] = None,
    ) -> list[QueryResult]:
        """
        Retrieve top-k entries similar to text, optionally filtered by metadata.

        filters examples:
            {"type": "HGT"}
            {"strain": "MRSA252", "min_score": 0.7}
            {"tag": ["quinolone", "beta-lactam"]}
        """
        qvec = self.encoder.encode(text)
        all_shards = [self._hot] + self._frozen
        return self._qe.query(all_shards, qvec, k=k, filters=filters)

    def query_vec(
        self,
        vec: np.ndarray,
        k: int = 10,
        filters: Optional[dict[str, Any]] = None,
    ) -> list[QueryResult]:
        """Query with a pre-encoded vector (skip encoding step)."""
        all_shards = [self._hot] + self._frozen
        return self._qe.query(all_shards, vec, k=k, filters=filters)

    def get_meta(self, entry_id: str) -> Optional[dict]:
        return self.manifest.meta(entry_id)

    # -- lifecycle -------------------------------------------------------------

    def commit_run(self) -> str:
        """
        Seal the hot shard and open a fresh one.
        Call at the end of each autonomous iteration.
        Returns the name of the sealed shard.
        """
        sealed_name = f"run_{self._run_id:04d}"
        sealed_path = self.root / "frozen" / sealed_name

        # seal current hot
        self._hot.seal()

        # update manifest: rename all "hot" shard refs -> sealed_name
        for eid in self.manifest.all_ids():
            entry = self.manifest.get(eid)
            if entry and entry.get("shard") == "hot":
                entry["shard"] = sealed_name
        self.manifest._save()

        # move hot -> frozen directory
        hot_mm   = self.root / "hot.mm"
        hot_meta = self.root / "hot.meta"
        if hot_mm.exists():
            shutil.move(str(hot_mm),   str(sealed_path.with_suffix(".mm")))
        if hot_meta.exists():
            shutil.move(str(hot_meta), str(sealed_path.with_suffix(".meta")))

        # reload as frozen
        sealed_shard = MemoryShard(
            sealed_path.with_suffix(".mm"), self.dim, frozen=True
        )
        self._frozen.insert(0, sealed_shard)   # most recent first

        # open new hot
        self._run_id += 1
        self._hot = self._open_hot()

        return sealed_name

    def status(self) -> dict:
        """Summary of memory state -- useful for the autonomous loop's logging."""
        return {
            "store":          self.store_name,
            "run_id":         self._run_id,
            "hot_entries":    self._hot.n_entries,
            "frozen_shards":  len(self._frozen),
            "frozen_entries": sum(s.n_entries for s in self._frozen),
            "manifest":       self.manifest.summary(),
        }

    # -- internal --------------------------------------------------------------

    def _open_hot(self) -> MemoryShard:
        return MemoryShard(self.root / "hot.mm", self.dim, frozen=False)

    def _load_frozen(self) -> list[MemoryShard]:
        frozen_dir = self.root / "frozen"
        shards = []
        for mm_path in sorted(frozen_dir.glob("*.mm"), reverse=True):
            shards.append(MemoryShard(mm_path, self.dim, frozen=True))
        return shards

    def _infer_run_id(self) -> int:
        frozen_dir = self.root / "frozen"
        existing = list(frozen_dir.glob("run_*.mm"))
        if not existing:
            return 1
        nums = []
        for p in existing:
            try:
                nums.append(int(p.stem.split("_")[1]))
            except (IndexError, ValueError):
                pass
        return max(nums) + 1 if nums else 1
