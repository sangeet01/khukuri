"""
Khukuri Drug Resistance Sentinel — Watchdog 3/3.

The most important watchdog. Continuously monitors currently deployed
antibiotics against the live resistance landscape and alerts when a drug
is losing ground.

This is the system that should have existed before that toddler's AST report.

What it does:
    - Tracks a panel of monitored antibiotics (SMILES + clinical metadata)
    - Periodically re-scores each drug against the current threat cluster
      using PINCER's minimax fitness
    - Detects when a drug's resistance robustness drops below a threshold
    - Identifies WHICH mutation or HGT event is responsible for the decline
    - Recommends PINCER-generated replacement candidates from CompoundMemory
    - Logs a full alert with the resistance mechanism, gene, and timeline

Alert levels:
    GREEN  — drug holds across all known threats (score > alert_threshold)
    YELLOW — drug weakening against specific threat cluster (early warning)
    RED    — drug failing — resistance mechanism identified, replacement needed
    BLACK  — drug completely compromised — do not use (score near 0)

Designed for:
    - Hospital AMR surveillance
    - Pharmaceutical pipeline monitoring
    - Research tracking of resistance emergence

All alert history persists to disk via PersistentMemory.
"""

import json
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger('khukuri.watchdog.sentinel')


# ---------------------------------------------------------------------------
# Alert levels
# ---------------------------------------------------------------------------

GREEN  = "GREEN"
YELLOW = "YELLOW"
RED    = "RED"
BLACK  = "BLACK"

ALERT_COLORS = {
    GREEN:  "[OK]",
    YELLOW: "[WARN]",
    RED:    "[FAIL]",
    BLACK:  "[CRIT]",
}

DEFAULT_THRESHOLDS = {
    GREEN:  0.65,    # above this = green
    YELLOW: 0.45,    # above this = yellow
    RED:    0.25,    # above this = red
    BLACK:  0.0,     # below yellow threshold = black
}


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class MonitoredDrug:
    """A drug under continuous surveillance."""
    drug_id: str
    name: str
    smiles: str
    drug_class: str                    # beta-lactam, glycopeptide, etc.
    target_pathogen: str               # S. aureus MRSA, K. pneumoniae, etc.
    clinical_status: str = "deployed"  # deployed | pipeline | withdrawn
    notes: str = ""
    baseline_score: Optional[float] = None
    last_score: Optional[float] = None
    last_alert: Optional[str] = None
    alert_history: List[Dict] = field(default_factory=list)


@dataclass
class ResistanceAlert:
    """An alert that a drug is losing ground."""
    timestamp: str
    drug_id: str
    drug_name: str
    alert_level: str                   # GREEN / YELLOW / RED / BLACK
    current_score: float
    previous_score: Optional[float]
    score_delta: Optional[float]
    responsible_threats: List[Dict]    # top threats causing the decline
    resistance_mechanism: str          # mutation | HGT | both
    specific_genes: List[str]          # e.g. ["mecA", "vanA"]
    specific_mutations: List[str]      # e.g. ["A5N", "V10L"]
    recommendation: str
    replacement_candidates: List[str]  # SMILES from CompoundMemory

    def to_dict(self) -> Dict:
        return {
            "timestamp": self.timestamp,
            "drug_id": self.drug_id,
            "drug_name": self.drug_name,
            "alert_level": self.alert_level,
            "current_score": self.current_score,
            "previous_score": self.previous_score,
            "score_delta": self.score_delta,
            "responsible_threats": self.responsible_threats[:5],
            "resistance_mechanism": self.resistance_mechanism,
            "specific_genes": self.specific_genes,
            "specific_mutations": self.specific_mutations,
            "recommendation": self.recommendation,
            "replacement_candidates": self.replacement_candidates[:3],
        }

    def print_alert(self):
        icon = ALERT_COLORS.get(self.alert_level, "?")
        print(f"\n{'='*60}")
        print(f"  {icon} RESISTANCE ALERT - {self.alert_level}")
        print(f"  Drug:      {self.drug_name} ({self.drug_id})")
        print(f"  Time:      {self.timestamp}")
        print(f"{'='*60}")
        print(f"  Score:     {self.current_score:.4f}", end="")
        if self.score_delta is not None:
            print(f"  (Delta {self.score_delta:+.4f})", end="")
        print()
        print(f"  Mechanism: {self.resistance_mechanism}")
        if self.specific_genes:
            print(f"  Genes:     {', '.join(self.specific_genes)}")
        if self.specific_mutations:
            print(f"  Mutations: {', '.join(self.specific_mutations[:5])}")
        if self.responsible_threats:
            print(f"  Top threat: {self.responsible_threats[0].get('id','?')}"
                  f"  fitness={self.responsible_threats[0].get('fitness',0):.3f}")
        print(f"\n  ALERT: {self.recommendation}")
        if self.replacement_candidates:
            print(f"\n  Replacement candidates:")
            for s in self.replacement_candidates[:3]:
                print(f"    {s}")
        print(f"{'='*60}\n")


# ---------------------------------------------------------------------------
# Built-in MRSA antibiotic panel
# ---------------------------------------------------------------------------

MRSA_ANTIBIOTIC_PANEL = [
    MonitoredDrug(
        drug_id="vancomycin",
        name="Vancomycin",
        smiles="CC1C(O)C(Cl)c2cc3cc(O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O[C@@H]4O[C@H](CO[C@@H]5O[C@H](C)[C@@H](O)[C@H](N)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)cc(O)c3c(O)c21",
        drug_class="glycopeptide",
        target_pathogen="S. aureus MRSA",
        clinical_status="deployed",
        notes="Last-resort antibiotic for MRSA. VanA/VanB resistance emerging.",
    ),
    MonitoredDrug(
        drug_id="daptomycin",
        name="Daptomycin",
        smiles="CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N",
        drug_class="lipopeptide",
        target_pathogen="S. aureus MRSA",
        clinical_status="deployed",
        notes="Cell membrane active. Resistance via MprF mutations.",
    ),
    MonitoredDrug(
        drug_id="linezolid",
        name="Linezolid",
        smiles="CC(=O)N[C@@H]1CN(c2ccc(N3CCOCC3)c(F)c2)C(=O)O1",
        drug_class="oxazolidinone",
        target_pathogen="S. aureus MRSA",
        clinical_status="deployed",
        notes="Protein synthesis inhibitor. cfr gene mediates resistance.",
    ),
    MonitoredDrug(
        drug_id="ceftaroline",
        name="Ceftaroline",
        smiles="CC1(C(=O)N2C(C(=O)O)=C(CSc3nnc(NC(=O)C(=NOC(C)(C)C(=O)O)c4csc(N)n4)s3)CS2)SCC1=O",
        drug_class="cephalosporin",
        target_pathogen="S. aureus MRSA",
        clinical_status="deployed",
        notes="5th gen cephalosporin active against MRSA. PBP2a affinity.",
    ),
    MonitoredDrug(
        drug_id="tedizolid",
        name="Tedizolid",
        smiles="Cc1ccc(-c2cnc(F)cn2)cc1[C@@H]1CN(c2ccc(N3N=NC(CO)=C3)c(F)c2)C(=O)O1",
        drug_class="oxazolidinone",
        target_pathogen="S. aureus MRSA",
        clinical_status="deployed",
        notes="Next-gen oxazolidinone. Active against linezolid-resistant strains.",
    ),
    MonitoredDrug(
        drug_id="meropenem",
        name="Meropenem",
        smiles="C[C@@H]1[C@H]2CC(=C(N2C1=O)C(=O)O)SC3CC(N)C3",
        drug_class="carbapenem",
        target_pathogen="K. pneumoniae",
        clinical_status="deployed",
        notes="Broad spectrum. KPC carbapenemases conferring resistance.",
    ),
    MonitoredDrug(
        drug_id="colistin",
        name="Colistin",
        smiles="CCC(C)CC(C)CCCCC(=O)NC(CCN)C(=O)NC(C(C)O)C(=O)NC(CCN)C(=O)N1CCCC1C(=O)NC(CCN)C(=O)NC(CCN)C(=O)NC(C(C)O)C(=O)NC(CCN)C(=O)O",
        drug_class="polymyxin",
        target_pathogen="K. pneumoniae",
        clinical_status="deployed",
        notes="Last resort for carbapenem-resistant GNB. mcr genes spreading.",
    ),
]


# ---------------------------------------------------------------------------
# Drug Resistance Sentinel
# ---------------------------------------------------------------------------

class DrugResistanceSentinel:
    """
    Continuous sentinel monitoring deployed antibiotics against the
    live resistance landscape.

    The system that should have existed before that toddler's AST report.

    Usage:
        sentinel = DrugResistanceSentinel(project="mrsa_surveillance")

        # Load threat landscape
        sentinel.load_threats(viable_threats, hgt_threats)

        # Add drugs to monitor (or use built-in MRSA panel)
        sentinel.load_mrsa_panel()
        sentinel.add_drug(MonitoredDrug(...))

        # Run surveillance
        alerts = sentinel.run_surveillance()
        sentinel.print_dashboard()

        # Subscribe to real-time alerts
        sentinel.on_alert(lambda a: send_email(a))
        sentinel.watch(interval_hours=24)   # background monitoring
    """

    def __init__(
        self,
        project: str = "default",
        memory_dir: Optional[str] = None,
        alert_thresholds: Optional[Dict[str, float]] = None,
        fitness_fn=None,
        on_alert: Optional[Callable] = None,
    ):
        self.project = project
        self.memory_dir = Path(
            memory_dir or Path.home() / ".khukuri" / "sentinel"
        ) / project
        self.memory_dir.mkdir(parents=True, exist_ok=True)

        self.thresholds = alert_thresholds or DEFAULT_THRESHOLDS
        self.fitness_fn = fitness_fn          # ThreatAwareFitnessFunction
        self._on_alert = on_alert             # user callback

        self._drugs: Dict[str, MonitoredDrug] = {}
        self._threats: List = []
        self._hgt_threats: List = []
        self._all_threats: List = []
        self._alert_log: List[ResistanceAlert] = []

        self._restore()
        logger.info(
            f"DrugResistanceSentinel: project='{project}' | "
            f"monitoring {len(self._drugs)} drugs"
        )

    # ------------------------------------------------------------------
    # Setup
    # ------------------------------------------------------------------

    def load_mrsa_panel(self):
        """Load the built-in MRSA antibiotic surveillance panel."""
        for drug in MRSA_ANTIBIOTIC_PANEL:
            self._drugs[drug.drug_id] = drug
        logger.info(f"Sentinel: loaded {len(MRSA_ANTIBIOTIC_PANEL)} MRSA drugs")

    def add_drug(self, drug: MonitoredDrug):
        """Add a drug to the surveillance panel."""
        self._drugs[drug.drug_id] = drug
        logger.info(f"Sentinel: added {drug.name} to surveillance")

    def load_threats(self, threats: List, hgt_threats: Optional[List] = None):
        """Load the current threat landscape."""
        self._threats = threats
        self._hgt_threats = hgt_threats or []
        self._all_threats = threats + self._hgt_threats
        logger.info(
            f"Sentinel: {len(self._threats)} mutational + "
            f"{len(self._hgt_threats)} HGT threats loaded"
        )

    def load_threats_from_pincer(self, pincer_engine):
        """Load threats directly from a PincerEngine instance."""
        self._all_threats = getattr(pincer_engine, 'viable_threats', [])
        hgt_cluster = getattr(pincer_engine, '_hgt_cluster', None)
        self._hgt_threats = hgt_cluster.alleles if hgt_cluster else []
        self._threats = [
            t for t in self._all_threats
            if t not in self._hgt_threats
        ]
        logger.info(
            f"Sentinel: loaded {len(self._threats)} mutational + "
            f"{len(self._hgt_threats)} HGT threats from PincerEngine"
        )

    def set_fitness_fn(self, fitness_fn):
        """Set or replace the fitness function."""
        self.fitness_fn = fitness_fn

    def on_alert(self, callback: Callable):
        """Register a callback for alerts: callback(ResistanceAlert)."""
        self._on_alert = callback

    # ------------------------------------------------------------------
    # Surveillance
    # ------------------------------------------------------------------

    def run_surveillance(
        self,
        drug_ids: Optional[List[str]] = None,
        top_k_threats: int = 30,
    ) -> List[ResistanceAlert]:
        """
        Score all monitored drugs against the current threat landscape.
        Returns list of alerts for drugs that changed status.

        Args:
            drug_ids:      specific drugs to check (None = all)
            top_k_threats: threats to evaluate per drug (via Numen if available)
        """
        if not self._all_threats:
            logger.warning("Sentinel: no threats loaded — call load_threats() first")
            return []

        if not self.fitness_fn:
            logger.warning("Sentinel: no fitness function — using default surrogate")
            from src.resistance.threat_fitness import ThreatAwareFitnessFunction
            self.fitness_fn = ThreatAwareFitnessFunction()

        drugs_to_check = {
            k: v for k, v in self._drugs.items()
            if drug_ids is None or k in drug_ids
        }

        alerts = []
        for drug_id, drug in drugs_to_check.items():
            alert = self._score_drug(drug, top_k_threats)
            if alert:
                alerts.append(alert)
                self._alert_log.append(alert)
                if self._on_alert:
                    self._on_alert(alert)

        self._save()
        logger.info(
            f"Sentinel: surveillance complete — "
            f"{len(drugs_to_check)} drugs checked, {len(alerts)} alerts"
        )
        return alerts

    def _score_drug(self, drug: MonitoredDrug, top_k: int) -> Optional[ResistanceAlert]:
        """Score one drug and generate alert if status changed."""
        try:
            # Get relevant threats
            threats = self._get_relevant_threats(drug, top_k)
            if not threats:
                return None

            # Score
            score = self.fitness_fn(drug.smiles, threats)

            # Determine alert level
            alert_level = self._score_to_level(score)

            # Determine if alert needed
            prev_score = drug.last_score
            prev_level = drug.last_alert
            score_delta = (score - prev_score) if prev_score is not None else None

            # Always alert on first check, RED/BLACK, or significant drop
            needs_alert = (
                prev_score is None
                or alert_level in (RED, BLACK)
                or (score_delta is not None and score_delta < -0.05)
                or alert_level != prev_level
            )

            # Update drug state
            if drug.baseline_score is None:
                drug.baseline_score = score
            drug.last_score = score
            drug.last_alert = alert_level

            if not needs_alert:
                return None

            # Identify responsible threats
            responsible = self._identify_responsible_threats(
                drug.smiles, threats, score
            )

            # Separate HGT vs mutational
            hgt_ids = {
                getattr(t, 'gene_id', '') for t in self._hgt_threats
            }
            specific_genes = list({
                getattr(t, 'gene_id', '')
                for t in responsible
                if getattr(t, 'gene_id', '') in hgt_ids
                and getattr(t, 'gene_id', '')
            })
            specific_mutations = []
            for t in responsible:
                for m in getattr(t, 'mutations', []):
                    wt = getattr(m, 'wild_type', '')
                    pos = getattr(m, 'position', '')
                    mut = getattr(m, 'mutant', '')
                    if wt and mut and mut != 'HGT':
                        specific_mutations.append(f"{wt}{pos}{mut}")

            mechanism = "none"
            if specific_genes and specific_mutations:
                mechanism = "both (mutation + HGT)"
            elif specific_genes:
                mechanism = f"HGT ({', '.join(specific_genes[:3])})"
            elif specific_mutations:
                mechanism = f"mutation ({', '.join(specific_mutations[:3])})"

            recommendation = self._build_recommendation(
                drug, score, alert_level, specific_genes, specific_mutations
            )

            replacements = self._find_replacements(drug, score)

            alert = ResistanceAlert(
                timestamp=datetime.now().isoformat(),
                drug_id=drug.drug_id,
                drug_name=drug.name,
                alert_level=alert_level,
                current_score=score,
                previous_score=prev_score,
                score_delta=score_delta,
                responsible_threats=[
                    {
                        "id": f"{getattr(t,'gene_id','')}_{i}",
                        "fitness": getattr(t, 'fitness', 0),
                        "type": "HGT" if getattr(t,'gene_id','') in hgt_ids else "mutation",
                    }
                    for i, t in enumerate(responsible[:5])
                ],
                resistance_mechanism=mechanism,
                specific_genes=specific_genes,
                specific_mutations=specific_mutations[:10],
                recommendation=recommendation,
                replacement_candidates=replacements,
            )

            drug.alert_history.append(alert.to_dict())
            return alert

        except Exception as exc:
            logger.error(f"Sentinel: scoring failed for {drug.drug_id}: {exc}")
            return None

    def _get_relevant_threats(self, drug: MonitoredDrug, top_k: int) -> List:
        """Get most relevant threats for this drug."""
        try:
            from src.integrations.numen_index import ThreatIndex
            ti = ThreatIndex(dim=256)
            ti.build(self._all_threats)
            return ti.get_worst_threats(drug.smiles, k=min(top_k, len(self._all_threats)))
        except Exception:
            return self._all_threats[:top_k]

    def _identify_responsible_threats(
        self, smiles: str, threats: List, drug_score: float
    ) -> List:
        """Identify which threats are most responsible for a low score."""
        scored = []
        for t in threats:
            try:
                individual = self.fitness_fn(smiles, [t])
                scored.append((individual, t))
            except Exception:
                pass
        scored.sort(key=lambda x: x[0])
        return [t for _, t in scored[:10]]

    def _score_to_level(self, score: float) -> str:
        if score >= self.thresholds[GREEN]:
            return GREEN
        elif score >= self.thresholds[YELLOW]:
            return YELLOW
        elif score >= self.thresholds[RED]:
            return RED
        else:
            return BLACK

    def _build_recommendation(
        self,
        drug: MonitoredDrug,
        score: float,
        level: str,
        genes: List[str],
        mutations: List[str],
    ) -> str:
        if level == GREEN:
            return f"{drug.name} holding — continue monitoring"
        elif level == YELLOW:
            parts = []
            if genes:
                parts.append(f"emerging HGT resistance via {', '.join(genes[:2])}")
            if mutations:
                parts.append(f"mutations {', '.join(mutations[:2])}")
            cause = "; ".join(parts) if parts else "unknown mechanism"
            return (
                f"{drug.name} weakening ({cause}). "
                f"Increase surveillance frequency. Consider combination therapy."
            )
        elif level == RED:
            return (
                f"{drug.name} failing (score={score:.3f}). "
                f"Resistance confirmed. Run PINCER for replacement candidates. "
                f"Consider clinical review."
            )
        else:  # BLACK
            return (
                f"{drug.name} COMPROMISED (score={score:.3f}). "
                f"Do not rely on this drug. "
                f"Immediate PINCER-guided replacement required."
            )

    def _find_replacements(self, drug: MonitoredDrug, score: float) -> List[str]:
        """Find replacement candidates from CompoundMemory."""
        if score > self.thresholds[YELLOW]:
            return []
        try:
            from src.integrations.numen_index import CompoundMemory
            cm = CompoundMemory()
            path = str(self.memory_dir.parent / "numen" / self.project / "compounds.json")
            cm.load(path)
            if len(cm) == 0:
                return []
            similar = cm.recall(drug.smiles, k=5)
            return [
                s["smiles"] for s in similar
                if s.get("pincer_score", 0) > score
                and s.get("smiles") != drug.smiles
            ][:3]
        except Exception:
            return []

    # ------------------------------------------------------------------
    # Dashboard
    # ------------------------------------------------------------------

    def print_dashboard(self):
        """Print surveillance dashboard for all monitored drugs."""
        print(f"\n{'='*65}")
        print(f"  KHUKURI DRUG RESISTANCE SENTINEL DASHBOARD")
        print(f"  Project: {self.project}  |  {datetime.now().strftime('%Y-%m-%d %H:%M')}")
        print(f"  Threats: {len(self._all_threats)} "
              f"({len(self._hgt_threats)} HGT + {len(self._threats)} mutational)")
        print(f"{'='*65}")
        print(f"  {'Drug':20s} {'Class':15s} {'Score':8s} {'Status':8s} {'Delta':8s}")
        print(f"  {'-'*60}")

        for drug in sorted(
            self._drugs.values(),
            key=lambda d: d.last_score or 0.0,
        ):
            icon = ALERT_COLORS.get(drug.last_alert or GREEN, "?")
            score_str = f"{drug.last_score:.4f}" if drug.last_score else "  N/A  "
            delta_str = ""
            if drug.last_score and drug.baseline_score:
                delta = drug.last_score - drug.baseline_score
                delta_str = f"{delta:+.4f}"
            print(
                f"  {icon:6s} {drug.name:18s} {drug.drug_class:15s} "
                f"{score_str:8s} {(drug.last_alert or 'UNSEEN'):8s} {delta_str}"
            )

        print(f"{'='*65}")
        if self._alert_log:
            print(f"  Recent alerts: {len(self._alert_log)}")
            for alert in self._alert_log[-3:]:
                icon = ALERT_COLORS.get(alert.alert_level, "?")
                print(f"    {icon} {alert.drug_name}: {alert.recommendation[:60]}")
        print(f"{'='*65}\n")

    def get_status(self) -> Dict[str, Any]:
        """Return current status of all monitored drugs."""
        return {
            drug_id: {
                "name": drug.name,
                "score": drug.last_score,
                "alert": drug.last_alert,
                "baseline": drug.baseline_score,
                "delta": (
                    (drug.last_score - drug.baseline_score)
                    if drug.last_score and drug.baseline_score else None
                ),
            }
            for drug_id, drug in self._drugs.items()
        }

    def get_critical_drugs(self) -> List[MonitoredDrug]:
        """Return drugs currently at RED or BLACK alert level."""
        return [
            d for d in self._drugs.values()
            if d.last_alert in (RED, BLACK)
        ]

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def _save(self):
        try:
            state = {
                "drugs": {
                    k: {
                        "drug_id": v.drug_id,
                        "name": v.name,
                        "smiles": v.smiles,
                        "drug_class": v.drug_class,
                        "target_pathogen": v.target_pathogen,
                        "clinical_status": v.clinical_status,
                        "baseline_score": v.baseline_score,
                        "last_score": v.last_score,
                        "last_alert": v.last_alert,
                        "alert_history": v.alert_history[-20:],
                    }
                    for k, v in self._drugs.items()
                },
                "alert_count": len(self._alert_log),
            }
            path = self.memory_dir / "sentinel_state.json"
            path.write_text(json.dumps(state, indent=2, default=str))
        except Exception as exc:
            logger.error(f"Sentinel: save failed: {exc}")

    def _restore(self):
        path = self.memory_dir / "sentinel_state.json"
        if not path.exists():
            return
        try:
            state = json.loads(path.read_text())
            for drug_id, data in state.get("drugs", {}).items():
                if drug_id in self._drugs:
                    drug = self._drugs[drug_id]
                    drug.baseline_score = data.get("baseline_score")
                    drug.last_score = data.get("last_score")
                    drug.last_alert = data.get("last_alert")
                    drug.alert_history = data.get("alert_history", [])
            logger.info(f"Sentinel: state restored from {path}")
        except Exception as exc:
            logger.warning(f"Sentinel: restore failed: {exc}")
