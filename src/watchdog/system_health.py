"""
Khukuri System Health Watchdog — Watchdog 1/3.

Monitors that all Khukuri modules are alive, responsive, and uncorrupted.
Runs as a lightweight check before any discovery run, or on a schedule.

Checks:
    - Core imports (all modules loadable)
    - linearscript backend (SCRIPT round-trip valid)
    - KeyBox engine (NibbleEngine initialises, C/NumPy backend detected)
    - World model integrity (JSON files readable, state consistent)
    - Numen indexes (vectors non-empty, dimensions consistent)
    - PINCER engine (map_threats produces non-empty threat cluster)
    - Resistance database (CARD/built-in data accessible)
    - Molecule generator (produces valid SMILES)
    - Memory dir (readable/writable, not full)

Returns a HealthReport with status per component and overall system status.
Logs warnings for degraded components, errors for failed ones.
"""

import json
import logging
import os
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger('khukuri.watchdog')


# ---------------------------------------------------------------------------
# Health report
# ---------------------------------------------------------------------------

@dataclass
class ComponentHealth:
    name: str
    status: str          # "ok" | "degraded" | "failed"
    message: str = ""
    latency_ms: float = 0.0
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class HealthReport:
    timestamp: str
    overall: str         # "healthy" | "degraded" | "critical"
    components: List[ComponentHealth] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {
            "timestamp": self.timestamp,
            "overall": self.overall,
            "warnings": self.warnings,
            "errors": self.errors,
            "components": [
                {
                    "name": c.name,
                    "status": c.status,
                    "message": c.message,
                    "latency_ms": round(c.latency_ms, 2),
                    "details": c.details,
                }
                for c in self.components
            ],
        }

    def print_summary(self):
        icon = {"healthy": "[OK]", "degraded": "[WARN]", "critical": "[FAIL]"}
        print(f"\n{'='*55}")
        print(f"  Khukuri System Health  {icon.get(self.overall, '?')} {self.overall.upper()}")
        print(f"  {self.timestamp}")
        print(f"{'='*55}")
        for c in self.components:
            sym = {"ok": "OK", "degraded": "WARN", "failed": "FAIL"}.get(c.status, "?")
            print(f"  {sym:6s} {c.name:30s} {c.status:8s}  {c.message}")
        if self.warnings:
            print(f"\n  Warnings:")
            for w in self.warnings:
                print(f"    ! {w}")
        if self.errors:
            print(f"\n  Errors:")
            for e in self.errors:
                print(f"    X {e}")
        print(f"{'='*55}\n")


# ---------------------------------------------------------------------------
# Individual checks
# ---------------------------------------------------------------------------

def _check(name: str, fn) -> ComponentHealth:
    """Run a check function, catch exceptions, return ComponentHealth."""
    t0 = time.time()
    try:
        result = fn()
        latency = (time.time() - t0) * 1000
        if isinstance(result, dict):
            status = result.get("status", "ok")
            message = result.get("message", "")
            details = result.get("details", {})
        else:
            status, message, details = "ok", str(result) if result else "", {}
        return ComponentHealth(name, status, message, latency, details)
    except Exception as exc:
        latency = (time.time() - t0) * 1000
        return ComponentHealth(name, "failed", str(exc)[:120], latency)


def check_core_imports() -> Dict:
    import importlib
    failed = []
    for mod in [
        "src.core", "src.resistance", "src.molecule_design",
        "src.world_model", "src.integrations", "src.workflows",
    ]:
        try:
            importlib.import_module(mod)
        except Exception as e:
            failed.append(f"{mod}: {e}")
    if failed:
        return {"status": "failed", "message": f"{len(failed)} modules failed",
                "details": {"failed": failed}}
    return {"status": "ok", "message": "all modules importable"}


def check_linearscript() -> Dict:
    try:
        from script.rdkit_bridge import SCRIPTFromMol, MolFromSCRIPT
        from script.validator import is_valid_SCRIPT
        from rdkit import Chem
        mol = Chem.MolFromSmiles("c1ccc2ncccc2c1")
        s = SCRIPTFromMol(mol)
        back = MolFromSCRIPT(s)
        if back is None:
            return {"status": "degraded", "message": "round-trip failed"}
        valid = is_valid_SCRIPT(s)
        return {
            "status": "ok" if valid else "degraded",
            "message": f"round-trip {'valid' if valid else 'invalid'}",
            "details": {"script_str": s},
        }
    except ImportError:
        return {"status": "degraded",
                "message": "linearscript not installed (pip install linearscript>=3.0.2)"}


def check_keybox() -> Dict:
    try:
        import sys
        for p in [
            os.path.join(os.path.dirname(__file__), '..', '..', '..', 'keybox'),
            os.path.join(os.path.expanduser('~'), 'keybox'),
        ]:
            if os.path.isdir(os.path.join(p, 'box')):
                sys.path.insert(0, p)
                break
        from box.key.nibble_bridge import NibbleEngine
        engine = NibbleEngine(dim_x=10, dim_y=10, dim_z=10, resolution=1.0)
        return {
            "status": "ok",
            "message": f"NibbleEngine initialised (mode={engine.mode})",
            "details": {"mode": engine.mode},
        }
    except ImportError:
        return {"status": "degraded",
                "message": "KeyBox not found — clone sibling repo"}


def check_world_model(memory_dir: Optional[str] = None) -> Dict:
    try:
        from src.world_model.manager import WorldModelManager
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            wm = WorldModelManager(project="health_check",
                                   memory_dir=memory_dir or tmp)
            wm.record_compound("HC_001", smiles="c1ccccc1", pincer_score=0.5)
            wm._save()
            wm2 = WorldModelManager(project="health_check",
                                    memory_dir=memory_dir or tmp)
            ok = len(wm2.state.compounds) == 1
        return {
            "status": "ok" if ok else "degraded",
            "message": "save/restore cycle passed" if ok else "restore failed",
        }
    except Exception as e:
        return {"status": "failed", "message": str(e)}


def check_numen(memory_dir: Optional[str] = None) -> Dict:
    try:
        from src.integrations.numen_index import NumenRetriever, CompoundMemory
        import tempfile
        nr = NumenRetriever(dim=128)
        nr.add("t1", "MRSA mecA beta-lactam resistance")
        nr.add("t2", "vanA glycopeptide enterococcus")
        results = nr.search("beta-lactam", top_k=1)
        ok = len(results) > 0 and results[0]["id"] == "t1"
        return {
            "status": "ok" if ok else "degraded",
            "message": f"retrieval {'correct' if ok else 'incorrect'}",
            "details": {"top_result": results[0]["id"] if results else None},
        }
    except Exception as e:
        return {"status": "failed", "message": str(e)}


def check_pincer() -> Dict:
    try:
        from src.resistance.pincer_engine import MutationSpaceMapper
        mapper = MutationSpaceMapper()
        threats = mapper.map_binding_pocket("AMILVC", list(range(6)))
        n = len(threats)
        return {
            "status": "ok" if n > 0 else "degraded",
            "message": f"{n} threats mapped",
            "details": {"n_threats": n},
        }
    except Exception as e:
        return {"status": "failed", "message": str(e)}


def check_molecule_generator() -> Dict:
    try:
        from src.molecule_design.generator import MoleculeGenerator
        gen = MoleculeGenerator()
        mols = gen.generate_molecules(n_molecules=3,
                                      target_properties={"mw": (150, 500)})
        n = len(mols)
        return {
            "status": "ok" if n >= 3 else "degraded",
            "message": f"{n}/3 molecules generated",
            "details": {"backend": "SCRIPT" if gen.use_script else "SMILES"},
        }
    except Exception as e:
        return {"status": "failed", "message": str(e)}


def check_memory_dir(memory_dir: Optional[str] = None) -> Dict:
    path = Path(memory_dir or Path.home() / ".khukuri" / "memory")
    try:
        path.mkdir(parents=True, exist_ok=True)
        test_file = path / ".health_check"
        test_file.write_text("ok")
        test_file.unlink()
        # Check disk space
        stat = os.statvfs(str(path))
        free_gb = (stat.f_bavail * stat.f_frsize) / (1024 ** 3)
        status = "ok" if free_gb > 0.5 else "degraded"
        return {
            "status": status,
            "message": f"{free_gb:.1f}GB free at {path}",
            "details": {"path": str(path), "free_gb": round(free_gb, 2)},
        }
    except Exception as e:
        return {"status": "failed", "message": str(e)}


def check_resistance_db() -> Dict:
    try:
        from src.resistance.hgt_mapper import HGTKnowledgeGraphBuilder
        from src.world_model.knowledge_graph import KnowledgeGraph
        kg = KnowledgeGraph()
        builder = HGTKnowledgeGraphBuilder(kg)
        n = builder.load_builtin_mrsa(selective_pressure="high")
        return {
            "status": "ok" if n > 0 else "degraded",
            "message": f"{n} MRSA alleles in built-in DB",
            "details": {"n_alleles": n},
        }
    except Exception as e:
        return {"status": "failed", "message": str(e)}


# ---------------------------------------------------------------------------
# SystemHealthWatchdog
# ---------------------------------------------------------------------------

class SystemHealthWatchdog:
    """
    Khukuri system health watchdog.

    Runs all component checks and produces a HealthReport.
    Can be run manually or called at the start of any discovery run.

    Usage:
        watchdog = SystemHealthWatchdog()
        report = watchdog.check()
        report.print_summary()

        # Raise if critical
        watchdog.check(raise_on_critical=True)

        # Save report
        watchdog.save_report(report)
    """

    CHECKS = [
        ("Core imports",        check_core_imports),
        ("linearscript",        check_linearscript),
        ("KeyBox engine",       check_keybox),
        ("World model",         check_world_model),
        ("Numen indexes",       check_numen),
        ("PINCER engine",       check_pincer),
        ("Molecule generator",  check_molecule_generator),
        ("Resistance database", check_resistance_db),
        ("Memory directory",    check_memory_dir),
    ]

    def __init__(
        self,
        memory_dir: Optional[str] = None,
        report_dir: Optional[str] = None,
    ):
        self.memory_dir = memory_dir
        self.report_dir = Path(
            report_dir or Path.home() / ".khukuri" / "watchdog"
        )
        self.report_dir.mkdir(parents=True, exist_ok=True)

    def check(self, raise_on_critical: bool = False) -> HealthReport:
        """Run all checks. Returns HealthReport."""
        logger.info("SystemHealthWatchdog: running checks...")
        results = []

        for name, fn in self.CHECKS:
            health = _check(name, fn)
            results.append(health)
            if health.status == "failed":
                logger.error(f"  FAIL {name}: {health.message}")
            elif health.status == "degraded":
                logger.warning(f"  WARN {name}: {health.message}")
            else:
                logger.info(f"  OK   {name}: {health.message}")

        # Determine overall status
        statuses = [c.status for c in results]
        if "failed" in statuses:
            overall = "critical"
        elif "degraded" in statuses:
            overall = "degraded"
        else:
            overall = "healthy"

        warnings = [f"{c.name}: {c.message}"
                    for c in results if c.status == "degraded"]
        errors   = [f"{c.name}: {c.message}"
                    for c in results if c.status == "failed"]

        report = HealthReport(
            timestamp=datetime.now().isoformat(),
            overall=overall,
            components=results,
            warnings=warnings,
            errors=errors,
        )

        logger.info(f"SystemHealthWatchdog: {overall.upper()} "
                    f"({len(errors)} errors, {len(warnings)} warnings)")

        if raise_on_critical and overall == "critical":
            raise RuntimeError(
                f"Khukuri system health check CRITICAL: {errors}"
            )

        return report

    def save_report(self, report: HealthReport, filename: Optional[str] = None):
        """Save health report to disk."""
        fname = filename or f"health_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        path = self.report_dir / fname
        path.write_text(json.dumps(report.to_dict(), indent=2))
        logger.info(f"Health report saved: {path}")
        return path

    def load_last_report(self) -> Optional[HealthReport]:
        """Load the most recent saved health report."""
        reports = sorted(self.report_dir.glob("health_*.json"))
        if not reports:
            return None
        data = json.loads(reports[-1].read_text())
        report = HealthReport(
            timestamp=data["timestamp"],
            overall=data["overall"],
            warnings=data.get("warnings", []),
            errors=data.get("errors", []),
        )
        for c in data.get("components", []):
            report.components.append(ComponentHealth(
                name=c["name"],
                status=c["status"],
                message=c["message"],
                latency_ms=c.get("latency_ms", 0),
                details=c.get("details", {}),
            ))
        return report
