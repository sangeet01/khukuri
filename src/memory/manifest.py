"""
manifest.py
-----------
The manifest is the only mutable file touched during a live run.
It tracks:
  - which shard each entry lives in + its row index
  - all metadata (score, strain, timestamps, etc.)
  - salience weight per entry (updated on recall + outcome feedback)

Format on disk: single JSON file, ~1KB per 10 entries.
Never large -- it holds metadata, not vectors.
"""

import json
import time
from pathlib import Path
from typing import Any, Optional


class Manifest:
    """
    Lightweight metadata + salience store.

    salience = recency_weight * win_rate * confidence
    Decays toward 0 if an entry is recalled but never reinforced.
    Strengthens toward 1 if recalled and confirmed as useful.
    """

    SALIENCE_INIT   = 0.5
    SALIENCE_BOOST  = 0.12   # on positive reinforcement
    SALIENCE_DECAY  = 0.04   # on negative / no reinforcement
    SALIENCE_MIN    = 0.05
    SALIENCE_MAX    = 1.0

    def __init__(self, path: Path):
        self.path = path
        self._data: dict[str, dict] = {}   # id -> entry
        if path.exists():
            self._load()

    # -- write -----------------------------------------------------------------

    def register(
        self,
        entry_id: str,
        shard_name: str,
        row_idx: int,
        text_key: str,
        meta: dict[str, Any],
    ) -> None:
        """Register a new entry. Called by MemoryShard on add()."""
        self._data[entry_id] = {
            "shard":     shard_name,
            "row":       row_idx,
            "text_key":  text_key,
            "meta":      meta,
            "salience":  self.SALIENCE_INIT,
            "recalls":   0,
            "wins":      0,
            "added_at":  time.time(),
            "last_seen": None,
        }
        self._save()

    def reinforce(self, entry_id: str, positive: bool = True) -> None:
        """
        Feedback signal from the autonomous loop.
        positive=True  -> compound worked / threat was real -> boost salience
        positive=False -> false positive / irrelevant -> decay salience
        """
        if entry_id not in self._data:
            return
        e = self._data[entry_id]
        e["last_seen"] = time.time()
        e["recalls"]  += 1
        if positive:
            e["wins"]    += 1
            e["salience"] = min(
                e["salience"] + self.SALIENCE_BOOST, self.SALIENCE_MAX
            )
        else:
            e["salience"] = max(
                e["salience"] - self.SALIENCE_DECAY, self.SALIENCE_MIN
            )
        self._save()

    def touch(self, entry_id: str) -> None:
        """Mark recalled (neutral -- no outcome yet)."""
        if entry_id in self._data:
            self._data[entry_id]["last_seen"] = time.time()
            self._data[entry_id]["recalls"]  += 1
            self._save()

    # -- read ------------------------------------------------------------------

    def get(self, entry_id: str) -> Optional[dict]:
        return self._data.get(entry_id)

    def all_ids(self) -> list[str]:
        return list(self._data.keys())

    def salience(self, entry_id: str) -> float:
        e = self._data.get(entry_id)
        return e["salience"] if e else 0.0

    def meta(self, entry_id: str) -> dict:
        e = self._data.get(entry_id)
        return e["meta"] if e else {}

    def filter_ids(self, predicates: dict[str, Any]) -> list[str]:
        """
        Return entry IDs whose metadata satisfies all predicates.

        Supports:
            exact match:    {"strain": "MRSA252"}
            min threshold:  {"min_score": 0.7}   -> meta["score"] >= 0.7
            max threshold:  {"max_score": 0.9}
            list member:    {"tag": ["quinolone", "beta-lactam"]}
        """
        results = []
        for eid, entry in self._data.items():
            m = entry["meta"]
            if self._matches(m, predicates):
                results.append(eid)
        return results

    def shard_rows(self, shard_name: str) -> dict[str, int]:
        """Return {entry_id: row_idx} for all entries in a given shard."""
        return {
            eid: e["row"]
            for eid, e in self._data.items()
            if e["shard"] == shard_name
        }

    def shards_present(self) -> list[str]:
        return list({e["shard"] for e in self._data.values()})

    def summary(self) -> dict:
        n = len(self._data)
        if n == 0:
            return {"entries": 0}
        saliences = [e["salience"] for e in self._data.values()]
        return {
            "entries":        n,
            "shards":         len(self.shards_present()),
            "mean_salience":  round(sum(saliences) / n, 3),
            "top_salience":   round(max(saliences), 3),
        }

    # -- internal --------------------------------------------------------------

    @staticmethod
    def _matches(meta: dict, predicates: dict) -> bool:
        for key, val in predicates.items():
            if key.startswith("min_"):
                field = key[4:]
                if meta.get(field, -1e9) < val:
                    return False
            elif key.startswith("max_"):
                field = key[4:]
                if meta.get(field, 1e9) > val:
                    return False
            elif isinstance(val, list):
                if meta.get(key) not in val:
                    return False
            else:
                if meta.get(key) != val:
                    return False
        return True

    def _save(self) -> None:
        tmp = self.path.with_suffix(".tmp")
        with open(tmp, "w") as f:
            json.dump(self._data, f, indent=2)
        tmp.replace(self.path)   # atomic write

    def _load(self) -> None:
        with open(self.path) as f:
            self._data = json.load(f)
