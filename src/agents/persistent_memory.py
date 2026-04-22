"""
Persistent memory for Khukuri agents.

Saves/loads agent state (failure patterns, meeting logs, learned constraints,
world-state snapshots) to disk so knowledge accumulates across sessions.
"""

import json
import logging
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger('khukuri')

DEFAULT_MEMORY_DIR = Path.home() / ".khukuri" / "memory"


class PersistentMemory:
    """
    Key-value store backed by JSON files.

    Usage:
        mem = PersistentMemory()
        mem.save("failures", {"database": [...], "patterns": {...}})
        data = mem.load("failures")
        mem.append("meeting_log", meeting)
    """

    def __init__(self, memory_dir: Optional[str] = None, project: str = "default"):
        self.memory_dir = Path(memory_dir or DEFAULT_MEMORY_DIR) / project
        self.memory_dir.mkdir(parents=True, exist_ok=True)
        self.project = project
        logger.info(f"PersistentMemory initialised at {self.memory_dir}")

    def save(self, key: str, data: Any) -> bool:
        """Persist `data` under `key`. Overwrites existing value."""
        path = self._path(key)
        try:
            tmp = path.with_suffix(".tmp")
            tmp.write_text(
                json.dumps(data, indent=2, default=_serialise),
                encoding="utf-8",
            )
            tmp.replace(path)
            logger.debug(f"Memory saved: {key}")
            return True
        except Exception as exc:
            logger.error(f"Memory save failed for '{key}': {exc}")
            return False

    def load(self, key: str, default: Any = None) -> Any:
        """Load data for `key`. Returns `default` if not found."""
        path = self._path(key)
        if not path.exists():
            return default if default is not None else {}
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            logger.error(f"Memory load failed for '{key}': {exc}")
            return default if default is not None else {}

    def append(self, key: str, item: Any) -> bool:
        """Append `item` to a list stored under `key`."""
        existing = self.load(key, default=[])
        if not isinstance(existing, list):
            existing = [existing]
        existing.append(item)
        if len(existing) > 1000:
            existing = existing[-1000:]
        return self.save(key, existing)

    def delete(self, key: str) -> bool:
        """Delete a memory key."""
        path = self._path(key)
        if path.exists():
            path.unlink()
            return True
        return False

    def list_keys(self) -> List[str]:
        """Return all stored keys for this project."""
        return [p.stem for p in self.memory_dir.glob("*.json")]

    def clear_all(self) -> bool:
        """Wipe all memory for this project (use with care)."""
        try:
            shutil.rmtree(self.memory_dir)
            self.memory_dir.mkdir(parents=True, exist_ok=True)
            logger.warning(f"Memory cleared for project '{self.project}'")
            return True
        except Exception as exc:
            logger.error(f"Memory clear failed: {exc}")
            return False

    def record_session(self, summary: Dict[str, Any]) -> bool:
        """Append a session summary to the audit log."""
        entry = {
            "timestamp": datetime.now().isoformat(),
            **summary,
        }
        return self.append("session_log", entry)

    def get_stats(self) -> Dict[str, Any]:
        """Return lightweight stats about what's stored."""
        stats: Dict[str, Any] = {"project": self.project, "keys": {}}
        for key in self.list_keys():
            data = self.load(key)
            if isinstance(data, list):
                stats["keys"][key] = f"{len(data)} items"
            elif isinstance(data, dict):
                stats["keys"][key] = f"{len(data)} fields"
            else:
                stats["keys"][key] = "1 value"
        return stats

    def _path(self, key: str) -> Path:
        safe = key.replace("/", "_").replace("\\", "_")
        return self.memory_dir / f"{safe}.json"


class MemoryMixin:
    """
    Mixin that gives an agent class automatic load/save of its state dict.
    """

    _memory_key: str = "agent_state"

    def _state_to_dict(self) -> Dict[str, Any]:
        """Override to return the dict you want to persist."""
        return {}

    def _load_from_dict(self, data: Dict[str, Any]) -> None:
        """Override to restore state from a saved dict."""

    def save_memory(self, memory: "PersistentMemory") -> bool:
        return memory.save(self._memory_key, self._state_to_dict())

    def load_memory(self, memory: "PersistentMemory") -> bool:
        data = memory.load(self._memory_key)
        if data:
            try:
                self._load_from_dict(data)
                logger.info(f"{self.__class__.__name__}: memory restored from disk")
                return True
            except Exception as exc:
                logger.warning(f"Memory restore failed: {exc}")
        return False


def _serialise(obj: Any) -> Any:
    """JSON default serialiser for non-standard types."""
    if hasattr(obj, "isoformat"):
        return obj.isoformat()
    if hasattr(obj, "__dict__"):
        return obj.__dict__
    return str(obj)
