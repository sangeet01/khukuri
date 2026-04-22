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

        # 6. World model (placeholder)
        self._world_state = {}

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

    def get_memory_stats(self) -> Dict[str, Any]:
        return self.memory.get_stats()
