"""AI agent system"""

# Core
from .base_agent import BaseAgent
from .virtual_lab import VirtualLab
from .pi_agent import PIAgent
from .meeting_system import MeetingSystem
from .failure_analyzer import FailureAnalyzer
from .competitive_teams import CompetitiveTeams
from .explainer import ExplainableAI

# New: multi-LLM provider
from .llm_provider import (
    LLMProvider,
    AnthropicProvider,
    OpenAIProvider,
    GeminiProvider,
    OllamaProvider,
    FallbackProvider,
    create_provider,
)

# New: persistent memory
from .persistent_memory import PersistentMemory, MemoryMixin

# New: parallel execution
from .parallel_executor import ParallelExecutor, AsyncExecutor, Task, TaskResult

# New: domain specialist agents
from .domain_agents import (
    DomainAgent,
    ChemistAgent,
    BiologistAgent,
    ToxicologistAgent,
    MicrobiologyAgent,
    GenomicsAgent,
    CheminformaticsAgent,
    ClinicalAgent,
    create_agent,
)

# New: peer debate
from .peer_debate import PeerDebate, PanelDiscussion, SocraticDialog, DebateRecord

# Orchestrators
from .hybrid_orchestrator import HybridOrchestrator   # legacy — unchanged
from .orchestrator import KhukuriOrchestrator          # new full-featured

__all__ = [
    # Core
    'BaseAgent', 'VirtualLab', 'PIAgent', 'MeetingSystem',
    'FailureAnalyzer', 'CompetitiveTeams', 'ExplainableAI',
    # LLM providers
    'LLMProvider', 'AnthropicProvider', 'OpenAIProvider',
    'GeminiProvider', 'OllamaProvider', 'FallbackProvider', 'create_provider',
    # Memory
    'PersistentMemory', 'MemoryMixin',
    # Execution
    'ParallelExecutor', 'AsyncExecutor', 'Task', 'TaskResult',
    # Domain agents
    'DomainAgent', 'ChemistAgent', 'BiologistAgent', 'ToxicologistAgent',
    'MicrobiologyAgent', 'GenomicsAgent', 'CheminformaticsAgent', 'ClinicalAgent',
    'create_agent',
    # Debate
    'PeerDebate', 'PanelDiscussion', 'SocraticDialog', 'DebateRecord',
    # Orchestrators
    'HybridOrchestrator', 'KhukuriOrchestrator',
]
