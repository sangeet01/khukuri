"""
query.py
--------
QueryEngine fans a query across all shards (hot + frozen),
merges results, applies metadata filters, reranks by salience.

This is the read path. All three memory types (Threat, Compound, Literature)
use this engine -- they just pass different filter predicates.
"""

import numpy as np
from typing import Any, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .shard import MemoryShard
    from .manifest import Manifest


class QueryResult:
    """A single result from a memory query."""

    __slots__ = ("entry_id", "score", "salience", "final_score", "meta")

    def __init__(
        self,
        entry_id: str,
        score: float,
        salience: float,
        meta: dict,
        salience_weight: float = 0.3,
    ):
        self.entry_id    = entry_id
        self.score       = score          # raw cosine similarity
        self.salience    = salience       # learned relevance weight
        self.meta        = meta
        # combined: similarity drives most of it, salience breaks ties
        self.final_score = (1 - salience_weight) * score + salience_weight * salience

    def __repr__(self) -> str:
        return (
            f"QueryResult({self.entry_id!r}, "
            f"sim={self.score:.3f}, sal={self.salience:.3f}, "
            f"final={self.final_score:.3f})"
        )


class QueryEngine:
    """
    Multi-shard query with metadata filtering and salience reranking.

    Usage:
        engine = QueryEngine(manifest, salience_weight=0.3)
        results = engine.query(
            shards,
            encoder.encode("some query text"),
            k=10,
            filters={"strain": "MRSA252", "min_score": 0.7},
        )
    """

    def __init__(self, manifest: "Manifest", salience_weight: float = 0.3):
        self.manifest        = manifest
        self.salience_weight = salience_weight

    def query(
        self,
        shards: list["MemoryShard"],
        query_vec: np.ndarray,
        k: int = 10,
        filters: Optional[dict[str, Any]] = None,
    ) -> list[QueryResult]:
        """
        Fan query across all shards, merge, filter, rerank.

        filters: dict of metadata predicates (see Manifest.filter_ids)
        Returns top-k QueryResult objects sorted by final_score descending.
        """
        # -- Step 1: metadata pre-filter (cheap, no vectors) ------------------
        if filters:
            candidate_ids = set(self.manifest.filter_ids(filters))
            if not candidate_ids:
                return []
        else:
            candidate_ids = None   # search everything

        # -- Step 2: per-shard vector search ----------------------------------
        # fetch more per shard than needed so merging doesn't cut good results
        per_shard_k = min(k * 3, k + 20)
        raw: list[tuple[str, float]] = []

        for shard in shards:
            if shard.n_entries == 0:
                continue

            # map candidate_ids to rows in this shard (if filtered)
            if candidate_ids is not None:
                shard_rows_map = self.manifest.shard_rows(shard.name)
                rows = [
                    shard_rows_map[eid]
                    for eid in candidate_ids
                    if eid in shard_rows_map
                ]
                if not rows:
                    continue
            else:
                rows = None

            hits = shard.query(query_vec, candidate_rows=rows, k=per_shard_k)
            raw.extend(hits)

        if not raw:
            return []

        # -- Step 3: deduplicate and take global top-k by raw score -----------
        seen = {}
        for eid, sc in raw:
            if eid not in seen or sc > seen[eid]:
                seen[eid] = sc

        # top-k by raw cosine before reranking
        top_raw = sorted(seen.items(), key=lambda x: x[1], reverse=True)[:k * 2]

        # -- Step 4: salience rerank ------------------------------------------
        results = []
        for eid, sc in top_raw:
            sal  = self.manifest.salience(eid)
            meta = self.manifest.meta(eid)
            results.append(QueryResult(eid, sc, sal, meta, self.salience_weight))

        results.sort(key=lambda r: r.final_score, reverse=True)

        # -- Step 5: mark recalled in manifest --------------------------------
        for r in results[:k]:
            self.manifest.touch(r.entry_id)

        return results[:k]
