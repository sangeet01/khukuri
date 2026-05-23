"""
integrations.py
---------------
Drop-in replacements for the three LimitNumen memory hacks in Khukuri.

Old (patch):
    from src.integrations import ThreatIndex, LiteratureIndex, CompoundMemory

New (this file):
    from src.memory.integrations import ThreatIndex, LiteratureIndex, CompoundMemory

Interface is backward-compatible where it matters.
All three wrap KhukuriMemory with domain-specific add/query methods.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from .store import KhukuriMemory
from .query import QueryResult


# -- ThreatIndex --------------------------------------------------------------

class ThreatIndex:
    """
    Indexes resistance threats (mutations + HGT events).
    Used by PINCER's ThreatAwareFitnessFunction for per-drug top-k retrieval.

    Old behaviour: build() indexes everything upfront, frozen.
    New behaviour: add threats incrementally as PINCER discovers them.
                   Persists across runs -- each run accumulates more threats.
    """

    def __init__(self, project_root: str | Path, dim: int = 4096):
        self._mem = KhukuriMemory(
            "threats", project_root, dim=dim
        )

    def add_threat(
        self,
        threat_id: str,
        description: str,
        threat_type: str = "mutation",   # "mutation" | "HGT"
        gene: str = "",
        strain: str = "",
        severity: float = 0.5,
        **kwargs,
    ) -> None:
        self._mem.add(
            threat_id,
            text=description,
            meta={
                "type":     threat_type,
                "gene":     gene,
                "strain":   strain,
                "severity": severity,
                **kwargs,
            },
        )

    def build(self, threats: list[dict]) -> None:
        """
        Bulk load. Compatible with old ThreatIndex.build(viable_threats) call.
        Each threat dict: {"id", "description", "type", "gene", "strain", ...}
        """
        entries = [
            {
                "id":   t.get("id", t.get("threat_id", str(i))),
                "text": t.get("description", t.get("text", str(t))),
                "meta": {k: v for k, v in t.items() if k not in ("id", "description", "text")},
            }
            for i, t in enumerate(threats)
        ]
        self._mem.add_batch(entries)

    def get_worst_threats(
        self,
        smiles: str,
        k: int = 20,
        threat_type: Optional[str] = None,
        strain: Optional[str] = None,
    ) -> list[QueryResult]:
        """
        Retrieve top-k most relevant threats for a given drug (SMILES).
        Salience ensures frequently-encountered threats surface first.
        """
        filters: dict[str, Any] = {}
        if threat_type:
            filters["type"] = threat_type
        if strain:
            filters["strain"] = strain

        return self._mem.query(smiles, k=k, filters=filters or None)

    def reinforce(self, threat_id: str, was_worst_case: bool) -> None:
        """Call after evolution: did this threat end up being the binding constraint?"""
        self._mem.reinforce(threat_id, positive=was_worst_case)

    def commit_run(self) -> None:
        self._mem.commit_run()

    def status(self) -> dict:
        return self._mem.status()


# -- LiteratureIndex ----------------------------------------------------------

class LiteratureIndex:
    """
    Indexes scientific literature for target discovery.
    Replaces the LimitNumen-backed literature mining hack.

    Adds per-target namespacing: papers are indexed with their target gene/protein
    as metadata, enabling filtered retrieval per target.
    """

    def __init__(self, project_root: str | Path, dim: int = 4096):
        self._mem = KhukuriMemory(
            "literature", project_root, dim=dim
        )

    def index_paper(
        self,
        paper_id: str,
        text: str,
        target: str = "",
        source: str = "pubmed",
        year: Optional[int] = None,
        relevance: float = 0.5,
    ) -> None:
        self._mem.add(
            paper_id,
            text=text,
            meta={
                "target":    target,
                "source":    source,
                "year":      year,
                "relevance": relevance,
            },
        )

    def index_pubmed_query(
        self,
        query_text: str,
        papers: list[dict],
        target: str = "",
    ) -> None:
        """
        Compatible with old LiteratureIndex.index_pubmed_query() call.
        papers: list of {"id", "text", "year", ...}
        """
        entries = [
            {
                "id":   p.get("id", p.get("pmid", str(i))),
                "text": p.get("text", p.get("abstract", "")),
                "meta": {
                    "target":    target,
                    "source":    "pubmed",
                    "year":      p.get("year"),
                    "relevance": p.get("relevance", 0.5),
                },
            }
            for i, p in enumerate(papers)
        ]
        self._mem.add_batch(entries)

    def search(
        self,
        query: str,
        k: int = 10,
        target: Optional[str] = None,
        min_year: Optional[int] = None,
    ) -> list[QueryResult]:
        filters: dict[str, Any] = {}
        if target:
            filters["target"] = target
        if min_year:
            filters["min_year"] = min_year
        return self._mem.query(query, k=k, filters=filters or None)

    def extract_targets(self, results: list[QueryResult]) -> list[str]:
        """Extract unique target names from query results."""
        targets = []
        seen = set()
        for r in results:
            t = r.meta.get("target", "")
            if t and t not in seen:
                targets.append(t)
                seen.add(t)
        return targets

    def commit_run(self) -> None:
        self._mem.commit_run()

    def status(self) -> dict:
        return self._mem.status()


# -- CompoundMemory -----------------------------------------------------------

class CompoundMemory:
    """
    Persistent cross-session memory for drug candidates.
    Replaces LimitNumen compound memory hack.

    Key improvement over the old version:
    - Filter by strain, min_score, tags before vector search
    - Salience weights past winners above past failures
    - commit_run() seals current run's compounds as a frozen shard
      -> next run's seeds are informed by all past runs automatically
    """

    def __init__(self, project_root: str | Path, dim: int = 4096):
        self._mem = KhukuriMemory(
            "compounds", project_root, dim=dim
        )

    def remember(
        self,
        smiles: str,
        compound_id: Optional[str] = None,
        pincer_score: float = 0.0,
        admet: Optional[dict] = None,
        strain: str = "",
        tags: Optional[list[str]] = None,
        **kwargs,
    ) -> str:
        """
        Store a compound. Returns the compound_id used.
        compound_id defaults to a hash of SMILES if not given.
        """
        if compound_id is None:
            compound_id = f"cpd_{abs(hash(smiles)) % 1_000_000:06d}"
        self._mem.add(
            compound_id,
            text=smiles,
            meta={
                "smiles":  smiles,
                "score":   pincer_score,
                "strain":  strain,
                "tags":    tags or [],
                "admet":   admet or {},
                **kwargs,
            },
        )
        return compound_id

    def recall(
        self,
        smiles: str,
        k: int = 5,
        strain: Optional[str] = None,
        min_score: Optional[float] = None,
        tags: Optional[list[str]] = None,
    ) -> list[QueryResult]:
        """
        Retrieve k past compounds similar to smiles.
        Salience surfaces past winners -- compounds that led to good outcomes
        in previous runs will rank above structurally similar but poor performers.
        """
        filters: dict[str, Any] = {}
        if strain:
            filters["strain"] = strain
        if min_score is not None:
            filters["min_score"] = min_score
        if tags:
            filters["tags"] = tags
        return self._mem.query(smiles, k=k, filters=filters or None)

    def reinforce(self, compound_id: str, outcome_positive: bool) -> None:
        """
        Feedback: did this compound advance (positive) or fail (negative)?
        Affects salience -> shapes future seed selection.
        """
        self._mem.reinforce(compound_id, positive=outcome_positive)

    def commit_run(self) -> str:
        """
        Seal this run's compounds.
        Next run's recall() will include them automatically.
        """
        return self._mem.commit_run()

    def top_candidates(
        self,
        k: int = 10,
        strain: Optional[str] = None,
        min_score: float = 0.0,
    ) -> list[QueryResult]:
        """Surface the highest-salience candidates across all runs."""
        # use a generic high-activity query to surface broadly
        return self.recall(
            "antibiotic drug candidate antimicrobial",
            k=k,
            strain=strain,
            min_score=min_score if min_score > 0 else None,
        )

    def status(self) -> dict:
        return self._mem.status()
