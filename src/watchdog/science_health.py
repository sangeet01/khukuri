"""
Khukuri Scientific Health Watchdog — Watchdog 2/3.

Monitors the scientific quality of PINCER evolutionary runs.
Detects when evolution has stalled, drifted, or converged prematurely
so you know whether to trust the apex candidate.

Checks:
    - Score progression    — is the apex score actually improving?
    - Population diversity — has the population converged too early?
    - Dead zone growth     — are dead zones filling (exploration happening)?
    - Threat coverage      — are all threat types being evaluated?
    - Fitness landscape    — is the fitness function differentiating candidates?
    - Stagnation detection — N generations with no improvement
    - HGT threat ratio     — are horizontal threats represented?
    - Apex stability       — does the apex change or flip erratically?

Produces a ScienceHealthReport with per-metric status and actionable
recommendations — e.g. "increase population_size", "add scaffold diversity",
"dead zone threshold too low".

Integrates with PincerEngine via callback — attach once, monitors every run.
"""

import logging
import statistics
from collections import deque
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger('khukuri.watchdog.science')


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class GenerationRecord:
    generation: int
    apex_score: float
    apex_smiles: str
    population_scores: List[float] = field(default_factory=list)
    dead_zones: int = 0
    explored_count: int = 0
    n_threats_evaluated: int = 0


@dataclass
class MetricHealth:
    name: str
    status: str           # "ok" | "warning" | "critical"
    value: Any
    message: str
    recommendation: str = ""


@dataclass
class ScienceHealthReport:
    timestamp: str
    run_id: str
    overall: str          # "healthy" | "warning" | "critical"
    metrics: List[MetricHealth] = field(default_factory=list)
    recommendations: List[str] = field(default_factory=list)
    generation_data: List[Dict] = field(default_factory=list)

    def print_summary(self):
        icon = {"healthy": "[OK]", "warning": "[WARN]", "critical": "[FAIL]"}
        print(f"\n{'='*60}")
        print(f"  PINCER Scientific Health  {icon.get(self.overall,'?')} {self.overall.upper()}")
        print(f"  Run: {self.run_id}  |  {self.timestamp}")
        print(f"{'='*60}")
        for m in self.metrics:
            sym = {"ok": "OK", "warning": "WARN", "critical": "FAIL"}.get(m.status, "?")
            print(f"  {sym:6s} {m.name:35s} {str(m.value):10s}  {m.message}")
        if self.recommendations:
            print(f"\n  Recommendations:")
            for r in self.recommendations:
                print(f"    -> {r}")
        print(f"{'='*60}\n")

    def is_trustworthy(self) -> bool:
        """Returns True if the apex candidate is scientifically trustworthy."""
        return self.overall != "critical"


# ---------------------------------------------------------------------------
# Scientific health checks
# ---------------------------------------------------------------------------

def check_score_progression(records: List[GenerationRecord]) -> MetricHealth:
    """Is the apex score actually improving across generations?"""
    if len(records) < 3:
        return MetricHealth("Score progression", "ok", "N/A",
                            "too few generations to assess")

    scores = [r.apex_score for r in records]
    first_half = scores[:len(scores)//2]
    second_half = scores[len(scores)//2:]

    improvement = max(second_half) - max(first_half)
    total_improvement = scores[-1] - scores[0]

    if total_improvement < 0.001:
        return MetricHealth(
            "Score progression", "critical",
            f"{total_improvement:+.4f}",
            "no improvement across all generations",
            "increase population_size or n_generations; check fitness function",
        )
    elif improvement < 0:
        return MetricHealth(
            "Score progression", "warning",
            f"{total_improvement:+.4f}",
            "improvement in first half only — possible premature convergence",
            "increase mutation rate or scaffold diversity in seed population",
        )
    else:
        return MetricHealth(
            "Score progression", "ok",
            f"{total_improvement:+.4f}",
            f"consistent improvement ({scores[0]:.4f} → {scores[-1]:.4f})",
        )


def check_population_diversity(records: List[GenerationRecord]) -> MetricHealth:
    """Has the population converged too early (low diversity)?"""
    recent = [r for r in records if r.population_scores]
    if not recent:
        return MetricHealth("Population diversity", "ok", "N/A",
                            "no population score data")

    last = recent[-1]
    if len(last.population_scores) < 2:
        return MetricHealth("Population diversity", "ok", "N/A",
                            "single candidate")

    std = statistics.stdev(last.population_scores)
    mean = statistics.mean(last.population_scores)
    cv = std / mean if mean > 0 else 0   # coefficient of variation

    if cv < 0.01:
        return MetricHealth(
            "Population diversity", "critical",
            f"CV={cv:.4f}",
            "population fully converged — all candidates identical",
            "add more scaffold diversity to seed pool; increase mutation operators",
        )
    elif cv < 0.05:
        return MetricHealth(
            "Population diversity", "warning",
            f"CV={cv:.4f}",
            "low diversity — convergence risk",
            "consider reseeding with diverse scaffolds mid-run",
        )
    else:
        return MetricHealth(
            "Population diversity", "ok",
            f"CV={cv:.4f}",
            f"healthy diversity (std={std:.4f})",
        )


def check_dead_zone_growth(records: List[GenerationRecord]) -> MetricHealth:
    """Are dead zones filling — is the algorithm actually exploring?"""
    if len(records) < 5:
        return MetricHealth("Dead zone growth", "ok", "N/A",
                            "too few generations")

    early_dz = records[2].dead_zones
    late_dz = records[-1].dead_zones
    growth = late_dz - early_dz

    if late_dz == 0:
        return MetricHealth(
            "Dead zone growth", "warning",
            "0 dead zones",
            "no dead zones recorded — exploration may be shallow",
            "lower dead_zone_threshold to capture more of chemical space",
        )
    elif growth == 0 and late_dz > 0:
        return MetricHealth(
            "Dead zone growth", "warning",
            f"{late_dz} static",
            "dead zones not growing after early exploration",
            "check if population is stuck in a local minimum",
        )
    else:
        return MetricHealth(
            "Dead zone growth", "ok",
            f"+{growth} zones",
            f"exploration active ({early_dz} → {late_dz} dead zones)",
        )


def check_stagnation(
    records: List[GenerationRecord],
    patience: int = 5,
    min_delta: float = 0.001,
) -> MetricHealth:
    """How many consecutive generations with no meaningful improvement?"""
    if len(records) < patience:
        return MetricHealth("Stagnation", "ok", "N/A", "too few generations")

    stagnant = 0
    for i in range(len(records) - 1, 0, -1):
        delta = records[i].apex_score - records[i-1].apex_score
        if abs(delta) < min_delta:
            stagnant += 1
        else:
            break

    if stagnant >= patience:
        return MetricHealth(
            "Stagnation", "critical",
            f"{stagnant} gens",
            f"no improvement for {stagnant} consecutive generations",
            "evolution stalled — consider early stopping or reseeding",
        )
    elif stagnant >= patience // 2:
        return MetricHealth(
            "Stagnation", "warning",
            f"{stagnant} gens",
            f"slowing — {stagnant} generations without improvement",
            "monitor closely; may need more generations",
        )
    else:
        return MetricHealth(
            "Stagnation", "ok",
            f"{stagnant} gens",
            "evolution actively improving",
        )


def check_fitness_differentiation(records: List[GenerationRecord]) -> MetricHealth:
    """Is the fitness function actually differentiating candidates?"""
    scored = [r for r in records if r.population_scores]
    if not scored:
        return MetricHealth("Fitness differentiation", "ok", "N/A",
                            "no population data")

    all_scores = []
    for r in scored[-3:]:    # last 3 generations
        all_scores.extend(r.population_scores)

    if not all_scores:
        return MetricHealth("Fitness differentiation", "ok", "N/A", "no data")

    score_range = max(all_scores) - min(all_scores)

    if score_range < 0.01:
        return MetricHealth(
            "Fitness differentiation", "critical",
            f"range={score_range:.4f}",
            "fitness function not differentiating — all candidates same score",
            "check ThreatAwareFitnessFunction weights; verify KeyBox/surrogate is active",
        )
    elif score_range < 0.05:
        return MetricHealth(
            "Fitness differentiation", "warning",
            f"range={score_range:.4f}",
            "low differentiation — weak selection pressure",
            "increase threat_weight in ThreatAwareFitnessFunction",
        )
    else:
        return MetricHealth(
            "Fitness differentiation", "ok",
            f"range={score_range:.4f}",
            "strong selection pressure",
        )


def check_apex_stability(records: List[GenerationRecord]) -> MetricHealth:
    """Is the apex candidate stable or flipping erratically?"""
    if len(records) < 4:
        return MetricHealth("Apex stability", "ok", "N/A", "too few generations")

    recent = records[-5:]
    smiles_set = set(r.apex_smiles for r in recent)
    flips = len(smiles_set)

    if flips == len(recent):
        return MetricHealth(
            "Apex stability", "warning",
            f"{flips} unique",
            "apex changing every generation — landscape too rugged",
            "increase population_size for smoother fitness landscape",
        )
    elif flips <= 2:
        return MetricHealth(
            "Apex stability", "ok",
            f"{flips} unique",
            "apex stable across recent generations",
        )
    else:
        return MetricHealth(
            "Apex stability", "ok",
            f"{flips} unique",
            "normal apex refinement",
        )


def check_hgt_coverage(
    records: List[GenerationRecord],
    total_threats: int,
    hgt_threats: int,
) -> MetricHealth:
    """Are horizontal gene transfer threats represented in evaluation?"""
    if total_threats == 0:
        return MetricHealth("HGT coverage", "ok", "N/A", "no threats loaded")

    hgt_ratio = hgt_threats / total_threats

    if hgt_threats == 0:
        return MetricHealth(
            "HGT coverage", "warning",
            "0 HGT threats",
            "no horizontal threats in threat cluster",
            "call engine.map_hgt_threats(knowledge_graph) for dual red team",
        )
    elif hgt_ratio < 0.02:
        return MetricHealth(
            "HGT coverage", "warning",
            f"{hgt_ratio:.1%}",
            "very low HGT proportion",
            "add more resistance genes to KnowledgeGraph",
        )
    else:
        return MetricHealth(
            "HGT coverage", "ok",
            f"{hgt_ratio:.1%}",
            f"{hgt_threats}/{total_threats} threats from HGT",
        )


def check_explored_coverage(records: List[GenerationRecord]) -> MetricHealth:
    """How much chemical space has been explored?"""
    if not records:
        return MetricHealth("Exploration coverage", "ok", "N/A", "no data")

    last = records[-1]
    explored = last.explored_count

    if explored < 10:
        return MetricHealth(
            "Exploration coverage", "warning",
            f"{explored} candidates",
            "very narrow exploration",
            "increase population_size or n_generations",
        )
    elif explored < 50:
        return MetricHealth(
            "Exploration coverage", "ok",
            f"{explored} candidates",
            "moderate exploration",
        )
    else:
        return MetricHealth(
            "Exploration coverage", "ok",
            f"{explored} candidates",
            "broad exploration",
        )


# ---------------------------------------------------------------------------
# Scientific Health Watchdog
# ---------------------------------------------------------------------------

class ScientificHealthWatchdog:
    """
    Monitors the scientific quality of PINCER evolutionary runs.

    Attach to PincerEngine via make_callback() to monitor in real-time,
    or call analyse() post-run with generation records.

    Usage — attach to PINCER:
        watchdog = ScientificHealthWatchdog(run_id="mrsa_run1")
        callback = watchdog.make_callback()
        apex = pincer.evolve(seeds, fitness_fn=fitness_fn, callback=callback)
        report = watchdog.finalise(pincer.get_results())
        report.print_summary()

    Usage — post-run analysis:
        report = watchdog.analyse(generation_records)
        report.print_summary()
        if not report.is_trustworthy():
            logger.warning("Apex candidate not trustworthy — rerun with larger population")
    """

    def __init__(
        self,
        run_id: Optional[str] = None,
        stagnation_patience: int = 5,
        min_score_delta: float = 0.001,
        on_warning: Optional[Callable] = None,
        on_critical: Optional[Callable] = None,
    ):
        self.run_id = run_id or datetime.now().strftime("run_%Y%m%d_%H%M%S")
        self.stagnation_patience = stagnation_patience
        self.min_score_delta = min_score_delta
        self.on_warning = on_warning or (lambda msg: logger.warning(f"WARN {msg}"))
        self.on_critical = on_critical or (lambda msg: logger.error(f"FAIL {msg}"))

        self._records: List[GenerationRecord] = []
        self._total_threats = 0
        self._hgt_threats = 0
        self._last_report: Optional[ScienceHealthReport] = None

    def set_threat_counts(self, total: int, hgt: int):
        """Call before evolve() with threat cluster sizes."""
        self._total_threats = total
        self._hgt_threats = hgt

    def make_callback(self) -> Callable:
        """
        Returns a callback for PincerEngine.evolve().
        Signature: callback(generation, apex_candidate)
        """
        def callback(generation: int, apex):
            record = GenerationRecord(
                generation=generation,
                apex_score=apex.minimax_score,
                apex_smiles=apex.smiles,
                dead_zones=getattr(apex, 'dead_zones', 0),
                explored_count=getattr(apex, 'explored_count', generation),
            )
            self._records.append(record)

            # Live stagnation check every 5 generations
            if generation > 0 and generation % 5 == 0:
                stag = check_stagnation(
                    self._records,
                    patience=self.stagnation_patience,
                    min_delta=self.min_score_delta,
                )
                if stag.status == "critical":
                    self.on_critical(
                        f"Gen {generation}: {stag.message} — {stag.recommendation}"
                    )
                elif stag.status == "warning":
                    self.on_warning(
                        f"Gen {generation}: {stag.message}"
                    )

        return callback

    def record_generation(
        self,
        generation: int,
        apex_score: float,
        apex_smiles: str,
        population_scores: Optional[List[float]] = None,
        dead_zones: int = 0,
        explored_count: int = 0,
    ):
        """Manually record a generation (if not using make_callback)."""
        self._records.append(GenerationRecord(
            generation=generation,
            apex_score=apex_score,
            apex_smiles=apex_smiles,
            population_scores=population_scores or [],
            dead_zones=dead_zones,
            explored_count=explored_count,
        ))

    def finalise(self, pincer_results: Optional[Dict] = None) -> ScienceHealthReport:
        """
        Run all scientific checks and return final ScienceHealthReport.
        Call after evolve() completes.
        """
        if pincer_results:
            self._total_threats = pincer_results.get(
                'viable_threats_count', self._total_threats
            )
            self._hgt_threats = pincer_results.get(
                'hgt_threats_count', self._hgt_threats
            )
            # Sync dead zone and explored counts from final result
            if self._records and pincer_results.get('dead_zones_count') is not None:
                self._records[-1].dead_zones = pincer_results['dead_zones_count']
            if self._records and pincer_results.get('explored_count') is not None:
                self._records[-1].explored_count = pincer_results['explored_count']

        return self.analyse(self._records)

    def analyse(self, records: Optional[List[GenerationRecord]] = None) -> ScienceHealthReport:
        """Run all checks on provided (or stored) records."""
        records = records or self._records

        checks = [
            check_score_progression(records),
            check_population_diversity(records),
            check_dead_zone_growth(records),
            check_stagnation(records, self.stagnation_patience,
                             self.min_score_delta),
            check_fitness_differentiation(records),
            check_apex_stability(records),
            check_hgt_coverage(records, self._total_threats, self._hgt_threats),
            check_explored_coverage(records),
        ]

        statuses = [m.status for m in checks]
        if "critical" in statuses:
            overall = "critical"
        elif "warning" in statuses:
            overall = "warning"
        else:
            overall = "healthy"

        recommendations = [
            m.recommendation for m in checks
            if m.recommendation and m.status in ("warning", "critical")
        ]

        report = ScienceHealthReport(
            timestamp=datetime.now().isoformat(),
            run_id=self.run_id,
            overall=overall,
            metrics=checks,
            recommendations=recommendations,
            generation_data=[
                {
                    "gen": r.generation,
                    "score": r.apex_score,
                    "dead_zones": r.dead_zones,
                }
                for r in records
            ],
        )

        self._last_report = report
        return report

    def get_score_curve(self) -> List[Tuple[int, float]]:
        """Return (generation, apex_score) pairs for plotting."""
        return [(r.generation, r.apex_score) for r in self._records]

    def reset(self):
        """Reset for a new run."""
        self._records = []
        self._total_threats = 0
        self._hgt_threats = 0
        self._last_report = None
