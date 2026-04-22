"""
KhukuriOrchestrator — full-featured multi-agent orchestration system.

Integrates all 6 new capabilities in one unified interface.
"""

import logging
from typing import Any, Callable, Dict, List, Optional

logger = logging.getLogger('khukuri')

from .llm_provider import LLMProvider, create_provider
from .persistent_memory import PersistentMemory
from .parallel_executor import ParallelExecutor, Task
from .domain_agents import (
    ChemistAgent, BiologistAgent, ToxicologistAgent,
    MicrobiologyAgent, GenomicsAgent, CheminformaticsAgent, ClinicalAgent,
    DomainAgent,
)
from .peer_debate import PeerDebate, PanelDiscussion, SocraticDialog, DebateRecord
from .pi_agent import PIAgent
from .meeting_system import MeetingSystem
from .failure_analyzer import FailureAnalyzer
from .competitive_teams import CompetitiveTeams
from .explainer import ExplainableAI


class KhukuriOrchestrator:
    """
    Full-featured orchestrator with all 6 enhancements.

    Quick start:
        orch = KhukuriOrchestrator()
        results = orch.run_discovery("Staphylococcus aureus")
    """

    DEFAULT_PANEL = [
        "chemist", "biologist", "toxicologist",
        "microbiologist", "genomics",
    ]

    def __init__(
        self,
        provider: Optional[LLMProvider] = None,
        openai_client=None,
        provider_name: str = "auto",
        api_key: Optional[str] = None,
        model: Optional[str] = None,
        memory_dir: Optional[str] = None,
        project: str = "default",
        max_workers: int = 6,
        on_progress: Optional[Callable] = None,
    ):
        # 1. LLM Provider
        if provider is not None:
            self.provider = provider
        elif openai_client is not None:
            self.provider = create_provider("openai", openai_client=openai_client)
        else:
            self.provider = create_provider(
                provider_name, api_key=api_key, model=model
            )
        logger.info(f"KhukuriOrchestrator using: {self.provider.name}")

        # 2. Persistent memory
        self.memory = PersistentMemory(memory_dir=memory_dir, project=project)
        self.project = project

        # 3. Parallel executor
        self.executor = ParallelExecutor(
            max_workers=max_workers, on_progress=on_progress
        )

        # 4. Domain specialist agents
        self.agents: Dict[str, DomainAgent] = {
            "chemist": ChemistAgent(self.provider),
            "biologist": BiologistAgent(self.provider),
            "toxicologist": ToxicologistAgent(self.provider),
            "microbiologist": MicrobiologyAgent(self.provider),
            "genomics": GenomicsAgent(self.provider),
            "cheminformatics": CheminformaticsAgent(self.provider),
            "clinical": ClinicalAgent(self.provider),
        }

        # 5. Peer debate
        self._debate_pairs = {
            ("chemist", "toxicologist"): "safety vs potency tradeoffs",
            ("biologist", "genomics"): "target selection and resistance",
        }

        # 6. World model
        from ..world_model import WorldStateTracker, KnowledgeGraph, KosmosEngine
        self.world_state = WorldStateTracker()
        self.knowledge_graph = KnowledgeGraph()
        self.kosmos = KosmosEngine(self.world_state, self.knowledge_graph)
        
        # Load existing world state if available
        saved_state = self.memory.load("world_state")
        if saved_state:
            self.world_state.from_dict(saved_state)
            logger.info("KhukuriOrchestrator: restored world state from persistent memory")

        # Legacy subsystems
        _oa_client = openai_client
        self.lead_agent = PIAgent(openai_client=_oa_client)
        self.meeting_system = MeetingSystem(self.lead_agent, _oa_client)
        self.failure_analyzer = FailureAnalyzer(openai_client=_oa_client)
        self.competitive_teams = CompetitiveTeams(openai_client=_oa_client)
        self.explainer = ExplainableAI(openai_client=_oa_client)

        self.available_tools = {
            'NetworkAnalyzer': 'PPI network analysis and target discovery',
            'MoleculeGenerator': 'AI-powered small molecule generation',
            'ADMETPredictor': 'Drug-likeness, toxicity, pharmacokinetics',
            'ResistancePredictor': 'Resistance evolution and multi-target strategies',
        }

    def run_peer_debate(
        self,
        agent_a: str,
        agent_b: str,
        topic: str,
        context: Dict[str, Any],
        rounds: int = 2,
    ) -> DebateRecord:
        """Run a structured debate between two named domain agents."""
        a = self.agents.get(agent_a)
        b = self.agents.get(agent_b)
        if not a or not b:
            raise ValueError(
                f"Unknown agents: {agent_a}, {agent_b}. "
                f"Available: {list(self.agents.keys())}"
            )
        debate = PeerDebate(a, b, moderator_provider=self.provider)
        record = debate.run(topic, context, rounds=rounds)
        logger.info(
            f"Peer debate ({agent_a} vs {agent_b}): "
            f"consensus={'yes' if record.consensus else 'no'}"
        )
        return record

    def run_panel(
        self,
        topic: str,
        context: Dict[str, Any],
        agent_names: Optional[List[str]] = None,
        max_rounds: int = 2,
    ) -> DebateRecord:
        """Run a round-robin panel discussion with N agents."""
        names = agent_names or self.DEFAULT_PANEL
        agents = [self.agents[n] for n in names if n in self.agents]
        panel = PanelDiscussion(agents, moderator_provider=self.provider)
        return panel.run(topic, context, max_rounds=max_rounds)

    def run_socratic(
        self,
        questioner: str,
        responder: str,
        claim: str,
        context: Dict[str, Any],
        n_questions: int = 3,
    ) -> DebateRecord:
        """Probe a claim with Socratic questioning between two agents."""
        q = self.agents[questioner]
        r = self.agents[responder]
        dialog = SocraticDialog(q, r, n_questions=n_questions)
        return dialog.run(claim, context)

    def run_discovery(self, species_name: str, n_iterations: int = 3, **kwargs) -> Dict[str, Any]:
        """
        Run iterative discovery using agents + world model with persistence.
        """
        logger.info(f"Starting iterative discovery for: {species_name} ({n_iterations} rounds)")
        
        # 1. Observe initial state
        self.kosmos.observe("discovery_start", {"species": species_name, **kwargs})
        
        for i in range(n_iterations):
            logger.info(f"--- Iteration {i+1}/{n_iterations} ---")
            
            # 2. Reason about current progress and strategy
            strategy = self.kosmos.reason(
                f"Iteration {i+1}: What is the optimal next step for {species_name}?",
                {"iteration": i+1, "world_summary": self.world_state.get_state_summary()}
            )
            
            # 3. High-level planning via expert panel
            panel_result = self.run_panel(
                topic=f"Strategic pivot for iteration {i+1} on {species_name}",
                context={"strategy": strategy, "iteration": i+1, **kwargs},
                max_rounds=2
            )
            
            # 4. Record consensus as a new hypothesis
            self.world_state.add_hypothesis(
                hypothesis=panel_result.final_recommendation,
                evidence=[f"Iterative consensus from panel in round {i+1}"],
                confidence=0.7 + (i * 0.1) # Increasing confidence over rounds
            )
            
            # 5. Peer debate on specific technical trade-offs (e.g. Chemist vs Toxicologist)
            if i % 2 == 0: # Every other round, do a deep dive
                tradeoff_record = self.run_peer_debate(
                    "chemist", "toxicologist",
                    topic=f"Safety vs Potency for {species_name} candidates",
                    context={"current_strategy": strategy}
                )
                self.knowledge_graph.add_relationship(
                    "chemist", "toxicologist", "debated",
                    {"result": tradeoff_record.consensus, "iteration": i+1}
                )

            # 6. Tool Execution Phase (The "Doing")
            logger.info(f"Iteration {i+1}: Executing autonomous tool layer...")
            tool_results = self._execute_tool_layer(strategy, i+1)
            
            # 7. Observe the new data from the tools
            self.kosmos.observe(f"tool_execution_round_{i+1}", tool_results)

            # 8. Periodic persistence
            self.memory.save("world_state", self.world_state.to_dict())
            logger.info(f"Iteration {i+1} complete. World state persisted.")

        return {
            "species": species_name,
            "iterations_completed": n_iterations,
            "world_state": self.world_state.get_state_summary(),
            "final_hypotheses": self.world_state.get_active_hypotheses()
        }

    def _execute_tool_layer(self, strategy: Dict[str, Any], iteration: int) -> Dict[str, Any]:
        """
        Autonomous execution of computational tools based on agent strategy.
        """
        results = {"iteration": iteration, "executed": []}
        
        # Determine which tools to run based on strategy/iteration
        # (In a full v2.0, this would be a more complex mapping)
        
        if iteration == 1:
            # Round 1: Target Identification
            logger.info("Auto-triggering Target Discovery...")
            from ..target_discovery.network_analyzer import NetworkAnalyzer
            from ..target_discovery.target_ranker import TargetRanker
            
            analyzer = NetworkAnalyzer()
            ranker = TargetRanker(analyzer)
            # Use real network analyzer logic with mock data for now
            mock_ppi = {("inhA", "katG"): 0.9, ("katG", "ahpC"): 0.7, ("inhA", "fabD"): 0.8}
            network = analyzer.build_ppi_network(mock_ppi)
            targets = ranker.rank_targets()
            
            for t in targets[:3]:
                self.world_state.update_target(t['protein'], {
                    "druggability": t.get('importance_score', 0.5),
                    "rank": t.get('rank')
                })
            results["executed"].append("NetworkAnalyzer + TargetRanker")
            
        elif iteration == 2:
            # Round 2: Lead Generation
            logger.info("Auto-triggering MoleculeGenerator...")
            from ..molecule_design.generator import MoleculeGenerator
            from rdkit import Chem
            generator = MoleculeGenerator()
            mols = generator.generate_library(max_compounds=10)
            for m in mols:
                smiles = Chem.MolToSmiles(m) if hasattr(m, 'GetNumAtoms') else str(m)
                self.world_state.update_compound(smiles, {
                    "status": "generated", 
                    "iteration": iteration,
                    "atoms": m.GetNumAtoms() if hasattr(m, 'GetNumAtoms') else 0
                })
            results["executed"].append("MoleculeGenerator")
            
        return results

    def get_memory_stats(self) -> Dict[str, Any]:
        return self.memory.get_stats()
