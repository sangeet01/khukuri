"""
WorldModelManager — persistent, unified interface for Khukuri's world model.

Wires together all five world model components and adds persistence so
knowledge accumulates across sessions rather than resetting each run:

    WorldStateTracker   — compounds, targets, strains, assays, hypotheses
    KnowledgeGraph      — relationship network (also hosts HGT data)
    HypothesisEngine    — generate, test, prioritise hypotheses
    LearningLoop        — Bayesian active learning over compound space
    KosmosEngine        — structured reasoning over world state

Persistence layout (~/.khukuri/memory/<project>/):
    world_state.json    — full WorldStateTracker snapshot
    knowledge_graph.json — node/edge data (serialised)
    learning_loop.json  — observations + model performance
    session_log.json    — lightweight per-run audit trail

Usage in AMRDiscoveryWorkflow:
    # Construction (restores previous session automatically)
    wm = WorldModelManager(project="mrsa_run1")

    # Before a run
    wm.start_run(pathogen="S. aureus MRSA252", context={...})

    # During a run — record everything
    wm.record_compound("COMP_001", smiles="...", admet={...}, pincer_score=0.72)
    wm.record_target("PBP2a", druggability=0.9, essentiality=0.85)
    wm.record_pincer_result(apex_smiles="...", score=0.63, threats=82)
    wm.generate_and_store_hypotheses()

    # After a run
    summary = wm.end_run(results={...})

    # Query accumulated knowledge
    top = wm.get_top_compounds(n=10)
    hyps = wm.get_actionable_hypotheses()
    insight = wm.reason("Which targets have the best druggability profile?")
"""

import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from .state_tracker import WorldStateTracker
from .knowledge_graph import KnowledgeGraph
from .hypothesis_engine import HypothesisEngine
from .learning_loop import LearningLoop
from .kosmos_engine import KosmosEngine

logger = logging.getLogger('khukuri')

DEFAULT_MEMORY_DIR = Path.home() / ".khukuri" / "memory"


class WorldModelManager:
    """
    Persistent, unified interface for Khukuri's world model.

    All five components are wired together here. State is saved to disk
    after every significant update and restored automatically on init.
    """

    def __init__(
        self,
        project: str = "default",
        memory_dir: Optional[str] = None,
        auto_restore: bool = True,
    ):
        self.project = project
        self.memory_dir = Path(memory_dir or DEFAULT_MEMORY_DIR) / project
        self.memory_dir.mkdir(parents=True, exist_ok=True)

        # Instantiate all five components
        self.state     = WorldStateTracker()
        self.graph     = KnowledgeGraph()
        self.hypotheses = HypothesisEngine(self.state, self.graph)
        self.learning  = LearningLoop()
        self.kosmos    = KosmosEngine(self.state, self.graph)

        self._run_count = 0
        self._current_run: Optional[Dict] = None

        if auto_restore:
            self._restore()

        logger.info(
            f"WorldModelManager: project='{project}' | "
            f"compounds={len(self.state.compounds)} "
            f"targets={len(self.state.targets)} "
            f"hypotheses={len(self.state.hypotheses)} "
            f"observations={len(self.learning.observations)}"
        )

    # ==================================================================
    # Run lifecycle
    # ==================================================================

    def start_run(self, pathogen: str, context: Optional[Dict] = None) -> Dict:
        """
        Open a new discovery run. Records context, increments run counter.
        Returns run metadata dict.
        """
        self._run_count += 1
        self._current_run = {
            "run_id": f"{self.project}_run{self._run_count:03d}",
            "pathogen": pathogen,
            "started_at": datetime.now().isoformat(),
            "context": context or {},
        }
        self.kosmos.observe("run_started", self._current_run)
        self.learning.next_iteration()
        logger.info(f"WorldModel: run {self._current_run['run_id']} started for {pathogen}")
        return self._current_run

    def end_run(self, results: Optional[Dict] = None) -> Dict:
        """
        Close current run, generate hypotheses, save everything to disk.
        Returns a summary of what was learned.
        """
        if self._current_run is None:
            logger.warning("end_run called with no active run")
            return {}

        # Generate new hypotheses from accumulated state
        new_hyps = self.hypotheses.generate_hypotheses()
        for hyp in new_hyps:
            self.state.add_hypothesis(
                hyp["hypothesis"],
                hyp.get("evidence_required", []),
                hyp.get("confidence", 0.5),
            )

        summary = {
            **self._current_run,
            "ended_at": datetime.now().isoformat(),
            "results_summary": {
                "compounds_added": len(self.state.compounds),
                "targets_known": len(self.state.targets),
                "hypotheses_total": len(self.state.hypotheses),
                "hypotheses_new": len(new_hyps),
                "learning_iteration": self.learning.iteration,
            },
            "run_results": results or {},
        }

        self.kosmos.observe("run_ended", summary)
        self._save()

        # Append to session log
        self._append_session_log(summary)

        self._current_run = None
        logger.info(f"WorldModel: run ended — {summary['results_summary']}")
        return summary

    # ==================================================================
    # Recording — called during a run
    # ==================================================================

    def record_compound(
        self,
        compound_id: str,
        smiles: str,
        script_str: Optional[str] = None,
        admet: Optional[Dict] = None,
        pincer_score: Optional[float] = None,
        docking_score: Optional[float] = None,
        targets: Optional[List[str]] = None,
    ):
        """Register a compound and update knowledge graph."""
        data = {
            "smiles": smiles,
            "script": script_str,
            "admet": admet or {},
            "pincer_score": pincer_score,
            "docking_score": docking_score,
        }
        self.state.update_compound(compound_id, data)
        self.graph.add_compound(compound_id, data)

        # Wire to targets in graph
        for target_id in (targets or []):
            if docking_score is not None:
                self.graph.add_binding(compound_id, target_id, docking_score)

        # Feed learning loop
        features = self._compound_features(smiles, admet or {})
        activity = pincer_score or (
            max(0.0, -docking_score / 10.0) if docking_score else 0.0
        )
        self.learning.add_observation(compound_id, features, {"activity": activity})

    def record_target(
        self,
        target_id: str,
        druggability: float = 0.5,
        essentiality: float = 0.5,
        sequence: Optional[str] = None,
        organism: Optional[str] = None,
    ):
        """Register a target protein."""
        data = {
            "druggability": druggability,
            "essentiality": essentiality,
            "sequence": sequence,
            "organism": organism,
        }
        self.state.update_target(target_id, data)
        self.graph.add_target(target_id, data)

    def record_strain(
        self,
        strain_id: str,
        resistance_profile: Optional[Dict] = None,
        mic_data: Optional[Dict] = None,
    ):
        """Register a bacterial strain."""
        data = {
            "resistance_profile": resistance_profile or {},
            "mic_data": mic_data or {},
        }
        self.state.update_strain(strain_id, data)
        self.graph.add_strain(strain_id, data)

    def record_pincer_result(
        self,
        apex_smiles: str,
        minimax_score: float,
        n_threats: int,
        dead_zones: int,
        generations: int,
        hgt_threats: int = 0,
    ):
        """Record a PINCER evolutionary run result."""
        compound_id = f"PINCER_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        self.record_compound(
            compound_id,
            smiles=apex_smiles,
            pincer_score=minimax_score,
        )
        self.state.add_experiment(
            "pincer_evolution",
            {
                "n_threats": n_threats,
                "hgt_threats": hgt_threats,
                "generations": generations,
            },
            {
                "apex_smiles": apex_smiles,
                "minimax_score": minimax_score,
                "dead_zones": dead_zones,
            },
        )
        self.kosmos.observe("pincer_result", {
            "apex": apex_smiles,
            "score": minimax_score,
            "threats": n_threats,
        })
        logger.info(
            f"WorldModel: PINCER result recorded "
            f"(apex={apex_smiles[:30]}, score={minimax_score:.4f})"
        )

    def record_assay(
        self,
        compound_id: str,
        strain_id: str,
        mic: Optional[float] = None,
        activity: Optional[float] = None,
    ):
        """Record a biological assay result and test relevant hypotheses."""
        assay_id = f"ASSAY_{compound_id}_{strain_id}"
        results = {"mic": mic, "activity": activity}
        self.state.add_assay_result(assay_id, compound_id, strain_id, results)

        # Auto-test any open hypotheses involving this compound
        self._auto_test_hypotheses(compound_id, results)

    def generate_and_store_hypotheses(self, context: Optional[Dict] = None) -> List[Dict]:
        """Generate new hypotheses from current state and store them."""
        new_hyps = self.hypotheses.generate_hypotheses(context)
        for hyp in new_hyps:
            self.state.add_hypothesis(
                hyp["hypothesis"],
                hyp.get("evidence_required", []),
                hyp.get("confidence", 0.5),
            )
        logger.info(f"WorldModel: {len(new_hyps)} new hypotheses generated")
        return new_hyps

    # ==================================================================
    # Querying accumulated knowledge
    # ==================================================================

    def get_top_compounds(self, n: int = 10, metric: str = "pincer_score") -> List[Dict]:
        """Return top n compounds by metric."""
        compounds = list(self.state.compounds.values())
        scored = [
            c for c in compounds
            if c.get(metric) is not None
        ]
        scored.sort(key=lambda c: c.get(metric, 0.0), reverse=True)
        return scored[:n]

    def get_actionable_hypotheses(self, min_confidence: float = 0.6) -> List[Dict]:
        """Return proposed hypotheses above confidence threshold, prioritised."""
        prioritised = self.hypotheses.prioritize_hypotheses()
        return [
            p for p in prioritised
            if p["hypothesis"].get("confidence", 0.0) >= min_confidence
        ]

    def get_best_next_candidates(
        self,
        candidates: List[Dict],
        n: int = 5,
    ) -> List[Dict]:
        """Use the learning loop to select the most promising candidates."""
        return self.learning.prioritize_next_experiments(candidates, n_select=n)

    def reason(self, question: str, context: Optional[Dict] = None) -> Dict:
        """Structured reasoning over accumulated world state."""
        return self.kosmos.reason(question, context)

    def get_summary(self) -> Dict:
        """Return a complete summary of world model state."""
        state_summary = self.state.get_state_summary()
        graph_stats = self.graph.get_statistics()
        learning_curve = self.learning.get_learning_curve()

        top_compounds = self.get_top_compounds(n=5)

        return {
            "project": self.project,
            "memory_dir": str(self.memory_dir),
            "state": state_summary,
            "graph": graph_stats,
            "learning": {
                "iteration": self.learning.iteration,
                "n_observations": len(self.learning.observations),
                "curve": learning_curve,
            },
            "top_compounds": [
                {
                    "id": c.get("compound_id"),
                    "smiles": c.get("smiles", "")[:50],
                    "pincer_score": c.get("pincer_score"),
                }
                for c in top_compounds
            ],
            "runs_completed": self._run_count,
        }

    # ==================================================================
    # Persistence
    # ==================================================================

    def _save(self):
        """Save all components to disk."""
        try:
            # WorldStateTracker
            self._write_json("world_state.json", self.state.to_dict())

            # KnowledgeGraph — serialise nodes + edges
            self._write_json("knowledge_graph.json", self._serialise_graph())

            # LearningLoop
            self._write_json("learning_loop.json", {
                "observations": self.learning.observations,
                "model_performance": self.learning.model_performance,
                "iteration": self.learning.iteration,
            })

            logger.info(f"WorldModel: state saved to {self.memory_dir}")
        except Exception as exc:
            logger.error(f"WorldModel: save failed — {exc}")

    def _restore(self):
        """Restore all components from disk if data exists."""
        restored = []

        # WorldStateTracker
        ws = self._read_json("world_state.json")
        if ws:
            self.state.from_dict(ws)
            restored.append("world_state")

        # KnowledgeGraph
        kg = self._read_json("knowledge_graph.json")
        if kg:
            self._deserialise_graph(kg)
            restored.append("knowledge_graph")

        # LearningLoop
        ll = self._read_json("learning_loop.json")
        if ll:
            self.learning.observations = ll.get("observations", [])
            self.learning.model_performance = ll.get("model_performance", {})
            self.learning.iteration = ll.get("iteration", 0)
            restored.append("learning_loop")

        if restored:
            logger.info(f"WorldModel: restored [{', '.join(restored)}]")

    def save(self):
        """Public save — call this explicitly if needed."""
        self._save()

    # ==================================================================
    # Internal helpers
    # ==================================================================

    def _auto_test_hypotheses(self, compound_id: str, results: Dict):
        """Test any open hypotheses involving compound_id."""
        for i, hyp in enumerate(self.state.hypotheses):
            if (
                hyp.get("status") == "proposed"
                and compound_id in hyp.get("hypothesis", "")
            ):
                evidence = {}
                if results.get("mic") is not None:
                    evidence["mic_data"] = results["mic"]
                if results.get("activity") is not None:
                    evidence["activity"] = results["activity"]
                if evidence:
                    self.hypotheses.test_hypothesis(i, evidence)

    @staticmethod
    def _compound_features(smiles: str, admet: Dict) -> Dict:
        """Extract numerical features for learning loop."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return {
                    "mw": Descriptors.MolWt(mol),
                    "logp": Descriptors.MolLogP(mol),
                    "hbd": Descriptors.NumHDonors(mol),
                    "hba": Descriptors.NumHAcceptors(mol),
                    **{k: v for k, v in admet.items()
                       if isinstance(v, (int, float))},
                }
        except Exception:
            pass
        return admet

    def _serialise_graph(self) -> Dict:
        """Serialise KnowledgeGraph to a JSON-safe dict."""
        nodes = [
            {"id": n, **{k: _json_safe(v) for k, v in d.items()}}
            for n, d in self.graph.graph.nodes(data=True)
        ]
        edges = [
            {
                "source": u,
                "target": v,
                **{k: _json_safe(val) for k, val in d.items()},
            }
            for u, v, d in self.graph.graph.edges(data=True)
        ]
        return {"nodes": nodes, "edges": edges}

    def _deserialise_graph(self, data: Dict):
        """Restore KnowledgeGraph from serialised dict."""
        for node in data.get("nodes", []):
            node = dict(node)   # copy so we can mutate
            node_id = node.pop("id")
            node.pop("node_type", None)    # KnowledgeGraph.add_* sets this
            node_type_raw = node.pop("node_type_raw", None)
            node_type = node.get("node_subtype", node_type_raw or "unknown")

            if "compound" in str(node_type):
                self.graph.add_compound(node_id, node)
            elif "target" in str(node_type):
                self.graph.add_target(node_id, node)
            elif "gene" in str(node_type):
                self.graph.add_gene(node_id, node)
            elif "strain" in str(node_type) or "species" in str(node_type):
                self.graph.add_strain(node_id, node)
            else:
                self.graph.graph.add_node(node_id, **node)

        for edge in data.get("edges", []):
            edge = dict(edge)
            src = edge.pop("source")
            tgt = edge.pop("target")
            rel = edge.get("relationship_type", "related")
            self.graph.add_relationship(src, tgt, rel, edge)

    def _write_json(self, filename: str, data: Any):
        path = self.memory_dir / filename
        tmp = path.with_suffix(".tmp")
        tmp.write_text(
            json.dumps(data, indent=2, default=_json_serialise),
            encoding="utf-8",
        )
        tmp.replace(path)

    def _read_json(self, filename: str) -> Optional[Any]:
        path = self.memory_dir / filename
        if not path.exists():
            return None
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            logger.warning(f"WorldModel: failed to read {filename}: {exc}")
            return None

    def _append_session_log(self, summary: Dict):
        log_path = self.memory_dir / "session_log.json"
        entries = []
        if log_path.exists():
            try:
                entries = json.loads(log_path.read_text())
            except Exception:
                entries = []
        entries.append({
            "timestamp": datetime.now().isoformat(),
            "run_id": summary.get("run_id"),
            "pathogen": summary.get("pathogen"),
            "results_summary": summary.get("results_summary", {}),
        })
        entries = entries[-200:]   # cap at 200 entries
        log_path.write_text(
            json.dumps(entries, indent=2, default=_json_serialise),
            encoding="utf-8",
        )


# ---------------------------------------------------------------------------
# JSON helpers
# ---------------------------------------------------------------------------

def _json_safe(val: Any) -> Any:
    if isinstance(val, np.ndarray):
        return val.tolist()
    if isinstance(val, (np.integer,)):
        return int(val)
    if isinstance(val, (np.floating,)):
        return float(val)
    return val


def _json_serialise(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if hasattr(obj, "isoformat"):
        return obj.isoformat()
    return str(obj)
