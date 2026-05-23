"""
khukuri/src/memory/__init__.py
------------------------------
KhukuriMemory -- tiered, persistent, salience-weighted memory system.

Replaces the three LimitNumen memory hacks in Khukuri:

    Old:
        from src.integrations import ThreatIndex, LiteratureIndex, CompoundMemory

    New:
        from src.memory import ThreatIndex, LiteratureIndex, CompoundMemory

The encoder (LimitNumen / NumenEncoder) is still used internally -- it is
one swappable component, not load-bearing architecture.
"""

from .store import KhukuriMemory
from .integrations import ThreatIndex, LiteratureIndex, CompoundMemory
from .query import QueryResult

__all__ = [
    "KhukuriMemory",
    "ThreatIndex",
    "LiteratureIndex",
    "CompoundMemory",
    "QueryResult",
]
