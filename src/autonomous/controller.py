"""
Khukuri Autonomous Loop Controller.

Two operating modes:

    SUPERVISED   — human approves each major decision before execution.
                   Good for research runs, new pathogens, uncertain terrain.

    LIGHTS_OUT   — fully autonomous, runs on schedule, self-monitors.
                   Pulls human back in only when a trigger threshold is crossed.

Trigger conditions that escalate from LIGHTS_OUT → human-in-loop:
    1. DrugResistanceSentinel flags RED or BLACK on any monitored drug
    2. Novel resistance gene detected not in KnowledgeGraph
    3. PINCER apex score drops >0.1 across 3 consecutive runs (landscape shift)
    4. New HGT transfer event found in literature (Numen detects it)
    5. HypothesisEngine confidence > 0.9 (needs wet lab confirmation)
    6. SystemHealthWatchdog reports CRITICAL

Human notification via: console, file log, email (if configured), webhook.

Cron integration:
    # Run lights-out every 24 hours
    0 2 * * * python -m khukuri.autonomous --mode lights_out --project mrsa

    # Run supervised (prompts for approval at decision points)
    python -m khukuri.autonomous --mode supervised --project mrsa
"""

import json
import logging
import os
import time
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

logger = logging.getLogger('khukuri.autonomous')


# ---------------------------------------------------------------------------
# Operating modes
# ---------------------------------------------------------------------------

class Mode(Enum):
    SUPERVISED  = "supervised"
    LIGHTS_OUT  = "lights_out"


# ---------------------------------------------------------------------------
# Escalation triggers
# ---------------------------------------------------------------------------

class TriggerReason(Enum):
    DRUG_RESISTANCE_RED    = "drug_resistance_red"
    DRUG_RESISTANCE_BLACK  = "drug_resistance_black"
    NOVEL_RESISTANCE_GENE  = "novel_resistance_gene"
    LANDSCAPE_SHIFT        = "landscape_shift"
    NEW_HGT_EVENT          = "new_hgt_event"
    HIGH_CONFIDENCE_HYP    = "high_confidence_hypothesis"
    SYSTEM_CRITICAL        = "system_critical"
    MANUAL_OVERRIDE        = "manual_override"


@dataclass
class EscalationEvent:
    """An event requiring human attention."""
    timestamp: str
    trigger: TriggerReason
    severity: str           # "info" | "warning" | "urgent" | "critical"
    message: str
    data: Dict[str, Any] = field(default_factory=dict)
    acknowledged: bool = False
    action_taken: str = ""

    def to_dict(self) -> Dict:
        return {
            "timestamp": self.timestamp,
            "trigger": self.trigger.value,
            "severity": self.severity,
            "message": self.message,
            "data": self.data,
            "acknowledged": self.acknowledged,
            "action_taken": self.action_taken,
        }


# ---------------------------------------------------------------------------
# Decision point (supervised mode)
# ---------------------------------------------------------------------------

@dataclass
class DecisionPoint:
    """A point where supervised mode waits for human approval."""
    question: str
    context: Dict[str, Any]
    options: List[str]
    default: str
    consequence: str        # what happens if approved vs rejected


# ---------------------------------------------------------------------------
# LLM Reasoning layer
# ---------------------------------------------------------------------------

class AutonomousReasoner:
    """
    LLM-powered reasoning over Khukuri's world state.

    Reads: sentinel dashboard, world model summary, hypothesis engine,
           recent PINCER results, open alerts.

    Decides: what to run next, whether to escalate, what seeds to use,
             which targets to prioritise.
    """

    SYSTEM_PROMPT = """You are the autonomous reasoning core of Khukuri, 
a drug discovery platform for antimicrobial resistance.

You receive a structured summary of the current world state and must decide:
1. Should a new PINCER discovery run be triggered? Why?
2. Which pathogen/target should be prioritised?
3. Are there any escalation triggers that need human attention?
4. What seed compounds should the next run use?
5. What hypotheses should be tested?

Be concise, scientific, and decisive. Return valid JSON only."""

    def __init__(self, provider=None):
        self.provider = provider
        self._decision_history: List[Dict] = []

    def reason(self, world_state: Dict, sentinel_status: Dict,
               open_hypotheses: List[Dict]) -> Dict:
        """
        Reason over current state and return a decision dict.

        Returns:
            {
                "run_discovery": bool,
                "target_pathogen": str,
                "priority_targets": [...],
                "seed_focus": str,          # scaffold focus for generator
                "escalate": bool,
                "escalation_reason": str,
                "reasoning": str,
                "next_actions": [...]
            }
        """
        if self.provider is None:
            return self._rule_based_decision(world_state, sentinel_status,
                                              open_hypotheses)
        prompt = self._build_prompt(world_state, sentinel_status,
                                    open_hypotheses)
        try:
            result = self.provider.complete_json(
                [{"role": "user", "content": prompt}],
                system=self.SYSTEM_PROMPT,
                temperature=0.2,
                max_tokens=800,
            )
            self._decision_history.append({
                "timestamp": datetime.now().isoformat(),
                "decision": result,
            })
            return result
        except Exception as exc:
            logger.warning(f"LLM reasoning failed: {exc} — using rule-based fallback")
            return self._rule_based_decision(world_state, sentinel_status,
                                              open_hypotheses)

    def _build_prompt(self, world_state, sentinel_status,
                      open_hypotheses) -> str:
        critical_drugs = [
            k for k, v in sentinel_status.items()
            if v.get("alert") in ("RED", "BLACK")
        ]
        return f"""
WORLD STATE SUMMARY:
{json.dumps(world_state, indent=2, default=str)}

SENTINEL STATUS:
{json.dumps(sentinel_status, indent=2, default=str)}

CRITICAL DRUGS: {critical_drugs}

OPEN HYPOTHESES (top 5):
{json.dumps(open_hypotheses[:5], indent=2, default=str)}

Based on this state, decide what Khukuri should do next.
Return JSON with keys: run_discovery, target_pathogen, priority_targets,
seed_focus, escalate, escalation_reason, reasoning, next_actions.
"""

    def _rule_based_decision(self, world_state, sentinel_status,
                              open_hypotheses) -> Dict:
        """Fallback when no LLM available — pure rule-based decisions."""
        critical = [
            k for k, v in sentinel_status.items()
            if v.get("alert") in ("RED", "BLACK")
        ]
        escalate = len(critical) > 0
        n_compounds = world_state.get("state", {}).get("compounds", 0)

        return {
            "run_discovery": True,
            "target_pathogen": "S. aureus MRSA",
            "priority_targets": ["PBP2a", "FtsZ"],
            "seed_focus": "quinolone" if n_compounds < 10 else "oxazolidinone",
            "escalate": escalate,
            "escalation_reason": (
                f"Drugs at RED/BLACK: {critical}" if escalate else ""
            ),
            "reasoning": (
                "Rule-based: prioritise MRSA, rotate scaffold diversity, "
                f"escalate={'yes' if escalate else 'no'}"
            ),
            "next_actions": [
                "Run PINCER with dual red team",
                "Update HGT threat cluster",
                "Test open hypotheses",
            ],
        }


# ---------------------------------------------------------------------------
# Main autonomous loop controller
# ---------------------------------------------------------------------------

class AutonomousLoopController:
    """
    Khukuri's autonomous loop controller.

    Orchestrates the full discovery cycle with configurable human oversight.

    Usage:
        # Supervised
        controller = AutonomousLoopController(
            mode=Mode.SUPERVISED,
            project="mrsa_campaign",
        )
        controller.run_cycle()

        # Lights-out (cron)
        controller = AutonomousLoopController(
            mode=Mode.LIGHTS_OUT,
            project="mrsa_campaign",
            on_escalation=send_email_alert,
        )
        controller.run_cycle()

        # Set mode at runtime
        controller.set_mode(Mode.LIGHTS_OUT)
    """

    # Thresholds for lights-out → human escalation
    ESCALATION_THRESHOLDS = {
        "apex_score_drop":        0.10,   # drop > 0.1 across 3 runs
        "hypothesis_confidence":   0.90,   # hypothesis needs wet lab
        "consecutive_stagnant_runs": 3,    # runs with no improvement
    }

    def __init__(
        self,
        mode: Mode = Mode.SUPERVISED,
        project: str = "default",
        memory_dir: Optional[str] = None,
        pathogen: str = "S. aureus MRSA",
        provider=None,
        on_escalation: Optional[Callable] = None,
        on_decision_needed: Optional[Callable] = None,
        pdb_path: Optional[str] = None,
    ):
        self.mode = mode
        self.project = project
        self.pathogen = pathogen
        self.on_escalation = on_escalation or self._default_escalation_handler
        self.on_decision_needed = on_decision_needed or self._console_decision
        self.pdb_path = pdb_path

        # Core components
        from ..world_model.manager import WorldModelManager
        from ..integrations.numen_index import KhukuriNumen
        from ..watchdog import SystemHealthWatchdog, DrugResistanceSentinel
        from ..watchdog.science_health import ScientificHealthWatchdog

        self.wm = WorldModelManager(project=project, memory_dir=memory_dir)
        self.numen = KhukuriNumen(project=project, memory_dir=memory_dir)
        self.health_watchdog = SystemHealthWatchdog(memory_dir=memory_dir)
        self.science_watchdog = ScientificHealthWatchdog()
        self.sentinel = DrugResistanceSentinel(project=project,
                                               memory_dir=memory_dir)
        self.reasoner = AutonomousReasoner(provider=provider)

        # State
        self._escalation_log: List[EscalationEvent] = []
        self._run_history: List[Dict] = []
        self._consecutive_stagnant = 0
        self._last_apex_scores: List[float] = []

        # Persistent state dir
        self._state_dir = Path(
            memory_dir or Path.home() / ".khukuri" / "memory"
        ) / project / "autonomous"
        self._state_dir.mkdir(parents=True, exist_ok=True)

        self._restore_state()

        logger.info(
            f"AutonomousLoopController: mode={mode.value} "
            f"project={project} pathogen={pathogen}"
        )

    # ==================================================================
    # Primary interface
    # ==================================================================

    def run_cycle(self) -> Dict[str, Any]:
        """
        Run one complete discovery cycle.

        In SUPERVISED mode: pauses at decision points for human approval.
        In LIGHTS_OUT mode: runs fully autonomously, escalates if needed.

        Returns cycle summary dict.
        """
        cycle_start = datetime.now()
        logger.info(
            f"{'='*50}\n"
            f"Khukuri Autonomous Cycle — {self.mode.value.upper()}\n"
            f"Project: {self.project} | {cycle_start.strftime('%Y-%m-%d %H:%M')}\n"
            f"{'='*50}"
        )

        # ---- Phase 1: System health check ---------------------------
        health = self._run_health_check()
        if health.overall == "critical":
            self._escalate(
                TriggerReason.SYSTEM_CRITICAL,
                "critical",
                f"System health CRITICAL: {health.errors}",
                {"health_report": health.to_dict()},
            )
            if self.mode == Mode.LIGHTS_OUT:
                logger.error("Lights-out: aborting cycle due to critical health")
                return {"status": "aborted", "reason": "system_critical"}

        # ---- Phase 2: Sentinel surveillance -------------------------
        sentinel_alerts = self._run_sentinel()
        sentinel_status = self.sentinel.get_status()

        # ---- Phase 3: LLM reasoning ---------------------------------
        decision = self.reasoner.reason(
            world_state=self.wm.get_summary(),
            sentinel_status=sentinel_status,
            open_hypotheses=self.wm.get_actionable_hypotheses(),
        )
        logger.info(f"Reasoner decision: {decision.get('reasoning', '')[:100]}")

        # ---- Phase 4: Check escalation triggers ---------------------
        escalations = self._check_escalation_triggers(
            decision, sentinel_alerts, health
        )

        # ---- Phase 5: Supervised approval (if needed) ---------------
        if self.mode == Mode.SUPERVISED or escalations:
            approved = self._request_approval(decision, escalations)
            if not approved:
                logger.info("Cycle paused — awaiting human approval")
                return {
                    "status": "paused",
                    "reason": "awaiting_approval",
                    "decision": decision,
                    "escalations": [e.to_dict() for e in escalations],
                }

        # ---- Phase 6: Discovery run ---------------------------------
        if decision.get("run_discovery", True):
            run_result = self._run_discovery(decision)
        else:
            run_result = {"status": "skipped", "reason": decision.get("reasoning")}

        # ---- Phase 7: Post-run analysis -----------------------------
        self._post_run_analysis(run_result, decision)

        # ---- Compile summary ----------------------------------------
        summary = {
            "cycle_id": cycle_start.strftime("%Y%m%d_%H%M%S"),
            "mode": self.mode.value,
            "project": self.project,
            "pathogen": self.pathogen,
            "timestamp": cycle_start.isoformat(),
            "duration_s": (datetime.now() - cycle_start).seconds,
            "health": health.overall,
            "sentinel_alerts": len(sentinel_alerts),
            "escalations": [e.to_dict() for e in escalations],
            "decision": decision,
            "run_result": run_result,
            "status": "completed",
        }

        self._save_cycle(summary)
        self._print_cycle_summary(summary)
        return summary

    def set_mode(self, mode: Mode):
        """Switch operating mode at runtime."""
        old = self.mode
        self.mode = mode
        logger.info(f"Mode changed: {old.value} → {mode.value}")

    def add_escalation_handler(self, handler: Callable):
        """Register a custom escalation handler (e.g. email, Slack, webhook)."""
        original = self.on_escalation
        def combined(event):
            original(event)
            handler(event)
        self.on_escalation = combined

    # ==================================================================
    # Execution phases
    # ==================================================================

    def _run_health_check(self):
        logger.info("Phase 1: System health check")
        report = self.health_watchdog.check()
        logger.info(f"Health: {report.overall}")
        return report

    def _run_sentinel(self) -> List:
        logger.info("Phase 2: Sentinel surveillance")
        self.sentinel.load_mrsa_panel()
        try:
            from ..resistance.pincer_engine import MutationSpaceMapper
            from ..resistance.hgt_mapper import HGTMapper
            from ..world_model.knowledge_graph import KnowledgeGraph
            from ..resistance.threat_fitness import ThreatAwareFitnessFunction

            mapper = MutationSpaceMapper()
            threats = mapper.map_binding_pocket('AMILVCFYWHDEKRST', list(range(12)))
            kg = KnowledgeGraph()
            hgt = HGTMapper(kg)
            cluster = hgt.map_threats('S. aureus', selective_pressure='high')
            self.sentinel.load_threats(threats, cluster.to_viable_threats())
            self.sentinel.set_fitness_fn(ThreatAwareFitnessFunction())
            alerts = self.sentinel.run_surveillance()
            logger.info(f"Sentinel: {len(alerts)} alerts")
            return alerts
        except Exception as exc:
            logger.warning(f"Sentinel run failed: {exc}")
            return []

    def _run_discovery(self, decision: Dict) -> Dict:
        logger.info("Phase 6: Discovery run")
        try:
            from ..workflows.amr_discovery import AMRDiscoveryWorkflow
            from ..molecule_design.generator import MoleculeGenerator

            wf = AMRDiscoveryWorkflow(
                project=self.project,
                pdb_path=self.pdb_path,
            )
            # Use Numen compound memory to enrich seeds
            gen = MoleculeGenerator()
            scaffold_focus = decision.get("seed_focus", "quinolone")
            seeds = gen.seed_for_pincer(n_seeds=10, scaffold_focus=scaffold_focus)

            # Enrich with past winners
            for smi in seeds[:3]:
                similar = self.numen.recall_similar(smi, k=2)
                for s in similar:
                    prev = s.get("smiles", "")
                    if prev and prev not in seeds:
                        seeds.append(prev)

            pathogen = decision.get("target_pathogen", self.pathogen)
            result = wf.run_discovery(
                pathogen=pathogen,
                priority="critical",
                n_iterations=2,
            )

            # Track apex score for landscape shift detection
            iterations = result.get("iterations", [])
            for it in iterations:
                pincer = it.get("pincer_results", {})
                apex = pincer.get("apex_drug", {})
                score = apex.get("minimax_score")
                if score is not None:
                    self._last_apex_scores.append(score)

            return result
        except Exception as exc:
            logger.error(f"Discovery run failed: {exc}")
            return {"status": "failed", "error": str(exc)}

    def _post_run_analysis(self, run_result: Dict, decision: Dict):
        """Update state, check for landscape shifts."""
        # Check stagnation
        if len(self._last_apex_scores) >= 3:
            recent = self._last_apex_scores[-3:]
            drop = max(recent) - min(recent)
            if drop > self.ESCALATION_THRESHOLDS["apex_score_drop"]:
                self._consecutive_stagnant += 1
                if self._consecutive_stagnant >= \
                        self.ESCALATION_THRESHOLDS["consecutive_stagnant_runs"]:
                    self._escalate(
                        TriggerReason.LANDSCAPE_SHIFT,
                        "warning",
                        f"Apex score dropped >{drop:.3f} across last 3 runs "
                        f"— resistance landscape may have shifted",
                        {"recent_scores": recent},
                    )
            else:
                self._consecutive_stagnant = 0

        # Check for novel resistance genes in literature
        self._check_novel_resistance()

        # Save state
        self._save_state()

    def _check_novel_resistance(self):
        """Use Numen to detect new HGT events in recent literature."""
        try:
            # Search for very recent papers on novel resistance
            results = self.numen.search_literature(
                "novel resistance gene horizontal transfer MRSA 2025 2026",
                top_k=5,
            )
            known_genes = {"mecA", "vanA", "vanB", "blaZ", "norA",
                           "tetM", "qnrS", "aac6"}
            for r in results:
                meta = r.get("meta", {})
                title = meta.get("title", "").lower()
                # Simple heuristic — novel gene mentioned not in known set
                for word in title.split():
                    if (len(word) >= 4 and
                            word not in known_genes and
                            any(c.isdigit() for c in word) and
                            word.isalnum()):
                        self._escalate(
                            TriggerReason.NEW_HGT_EVENT,
                            "info",
                            f"Possible novel resistance gene in literature: '{word}'",
                            {"paper": meta},
                        )
                        break
        except Exception:
            pass

    # ==================================================================
    # Escalation system
    # ==================================================================

    def _check_escalation_triggers(
        self,
        decision: Dict,
        sentinel_alerts: List,
        health,
    ) -> List[EscalationEvent]:
        """Check all escalation conditions. Return list of triggered events."""
        events = []

        # 1. RED/BLACK drugs
        for alert in sentinel_alerts:
            if alert.alert_level in ("RED", "BLACK"):
                events.append(EscalationEvent(
                    timestamp=datetime.now().isoformat(),
                    trigger=TriggerReason.DRUG_RESISTANCE_RED
                    if alert.alert_level == "RED"
                    else TriggerReason.DRUG_RESISTANCE_BLACK,
                    severity="urgent" if alert.alert_level == "RED" else "critical",
                    message=(
                        f"{alert.drug_name} at {alert.alert_level}: "
                        f"{alert.resistance_mechanism}"
                    ),
                    data=alert.to_dict(),
                ))

        # 2. High-confidence hypotheses
        for hyp in self.wm.get_actionable_hypotheses(min_confidence=0.9):
            events.append(EscalationEvent(
                timestamp=datetime.now().isoformat(),
                trigger=TriggerReason.HIGH_CONFIDENCE_HYP,
                severity="info",
                message=f"High-confidence hypothesis needs wet lab: "
                        f"{hyp.get('hypothesis',{}).get('hypothesis','')[:80]}",
                data={"hypothesis": hyp},
            ))

        # 3. LLM-flagged escalation
        if decision.get("escalate"):
            events.append(EscalationEvent(
                timestamp=datetime.now().isoformat(),
                trigger=TriggerReason.MANUAL_OVERRIDE,
                severity="warning",
                message=decision.get("escalation_reason", "LLM flagged escalation"),
                data={"decision": decision},
            ))

        # Fire handlers
        for event in events:
            self._escalation_log.append(event)
            self.on_escalation(event)

        return events

    def _escalate(self, trigger: TriggerReason, severity: str,
                  message: str, data: Dict = None):
        """Create and fire an escalation event."""
        event = EscalationEvent(
            timestamp=datetime.now().isoformat(),
            trigger=trigger,
            severity=severity,
            message=message,
            data=data or {},
        )
        self._escalation_log.append(event)
        self.on_escalation(event)

    # ==================================================================
    # Human-in-loop (supervised mode)
    # ==================================================================

    def _request_approval(
        self,
        decision: Dict,
        escalations: List[EscalationEvent],
    ) -> bool:
        """
        Request human approval before proceeding.
        In lights-out mode, only called when escalation triggers fire.
        """
        if escalations:
            # Always pause for critical escalations regardless of mode
            critical = [e for e in escalations
                        if e.severity in ("critical", "urgent")]
            if critical:
                dp = DecisionPoint(
                    question="Critical escalation detected. Proceed with discovery run?",
                    context={
                        "escalations": [e.to_dict() for e in critical],
                        "decision": decision,
                    },
                    options=["approve", "pause", "abort"],
                    default="pause",
                    consequence="approve=run PINCER | pause=hold for review | abort=stop",
                )
                response = self.on_decision_needed(dp)
                return response == "approve"

        if self.mode == Mode.SUPERVISED:
            dp = DecisionPoint(
                question=f"Run discovery for {decision.get('target_pathogen')}?",
                context={"decision": decision},
                options=["approve", "skip", "modify"],
                default="approve",
                consequence="approve=run | skip=skip this cycle | modify=edit params",
            )
            response = self.on_decision_needed(dp)
            return response in ("approve", "y", "yes", "")

        # Lights-out with no critical escalations = auto-approve
        return True

    # ==================================================================
    # Notification handlers
    # ==================================================================

    @staticmethod
    def _default_escalation_handler(event: EscalationEvent):
        """Default: log to console with clear formatting."""
        icons = {
            "info": "ℹ️",
            "warning": "⚠️",
            "urgent": "🔴",
            "critical": "⛔",
        }
        icon = icons.get(event.severity, "?")
        print(f"\n{icon} KHUKURI ESCALATION [{event.severity.upper()}]")
        print(f"   Trigger:  {event.trigger.value}")
        print(f"   Time:     {event.timestamp}")
        print(f"   Message:  {event.message}")
        print()

    @staticmethod
    def _console_decision(dp: DecisionPoint) -> str:
        """Console prompt for supervised mode decisions."""
        print(f"\n{'='*55}")
        print(f"  DECISION REQUIRED")
        print(f"{'='*55}")
        print(f"  {dp.question}")
        print(f"  Options: {' | '.join(dp.options)}")
        print(f"  Default: {dp.default}")
        print(f"  Consequence: {dp.consequence}")
        print(f"{'='*55}")
        try:
            response = input(f"  Your choice [{dp.default}]: ").strip().lower()
            return response if response else dp.default
        except (EOFError, KeyboardInterrupt):
            # Non-interactive (cron) — use default
            return dp.default

    # ==================================================================
    # Persistence
    # ==================================================================

    def _save_cycle(self, summary: Dict):
        path = self._state_dir / f"cycle_{summary['cycle_id']}.json"
        path.write_text(json.dumps(summary, indent=2, default=str))

    def _save_state(self):
        state = {
            "consecutive_stagnant": self._consecutive_stagnant,
            "last_apex_scores": self._last_apex_scores[-20:],
            "escalation_count": len(self._escalation_log),
        }
        (self._state_dir / "controller_state.json").write_text(
            json.dumps(state, indent=2)
        )

    def _restore_state(self):
        path = self._state_dir / "controller_state.json"
        if path.exists():
            try:
                state = json.loads(path.read_text())
                self._consecutive_stagnant = state.get("consecutive_stagnant", 0)
                self._last_apex_scores = state.get("last_apex_scores", [])
                logger.info(
                    f"Controller state restored "
                    f"(stagnant_runs={self._consecutive_stagnant})"
                )
            except Exception:
                pass

    def _print_cycle_summary(self, summary: Dict):
        status_icon = {"completed": "✅", "paused": "⏸️",
                       "aborted": "❌"}.get(summary["status"], "?")
        print(f"\n{'='*55}")
        print(f"  {status_icon} Cycle {summary['cycle_id']} — {summary['status'].upper()}")
        print(f"  Mode: {summary['mode']}  |  Duration: {summary['duration_s']}s")
        print(f"  Health: {summary['health']}  |  Alerts: {summary['sentinel_alerts']}")
        print(f"  Escalations: {len(summary['escalations'])}")
        if summary.get("decision", {}).get("reasoning"):
            print(f"  Reasoning: {summary['decision']['reasoning'][:80]}")
        print(f"{'='*55}\n")
