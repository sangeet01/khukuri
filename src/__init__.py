"""Khukuri Virtual Lab - AMR-Focused Drug Discovery Platform"""

__version__ = "3.0.0"
__author__ = "Sangeet Sharma"

from .core import MolecularScorer, Validator, setup_logger
from .target_discovery import NetworkAnalyzer, TargetRanker
from .molecule_design import MoleculeGenerator, PropertyOptimizer
from .docking import VinaWrapper, BindingSiteDetector
from .admet import calculate_drug_likeness, predict_admet

# AMR-specific modules
from .bioknowledge import ResistanceDatabase, PathogenDatabase, TargetProteinDB, DatabaseUpdater
from .genomics import ResistanceGenomicsAnalyzer, MutationTracker, OmicsIntegrator
from .microbiology import MICAnalyzer, AssayTracker, StrainManager

# World model — persistent, unified
from .world_model import (
    WorldStateTracker, KnowledgeGraph, LearningLoop,
    KosmosEngine, HypothesisEngine, WorldModelManager,
)

# Resistance — PINCER dual red team + fitness
from .resistance import (
    PincerEngine, MutationSpaceMapper, ViableMutant,
    ThreatAwareFitnessFunction,
    HGTMapper, TransferThreatCluster, make_pincer_with_hgt,
)

# Integrations — KeyBox + Numen
from .integrations import (
    KhukuriKeyBox, plug_keybox_into_pincer, create_mrsa_keybox, KEYBOX_AVAILABLE,
    NumenRetriever, ThreatIndex, LiteratureIndex, CompoundMemory, KhukuriNumen,
)

# Watchdog — continuous monitoring
from .watchdog import (
    SystemHealthWatchdog, ScientificHealthWatchdog, DrugResistanceSentinel,
    HealthReport, ComponentHealth, ScienceHealthReport, ResistanceAlert,
)

# Agents
from .agents.amr_agents import (MicrobiologyAgent, GenomicsAgent, CheminformaticsAgent,
                                 ResistanceCriticAgent, LiteratureAgent)

# Autonomous — scheduled, lights-out loops
from .autonomous import AutonomousLoopController, PINCERScheduler, Mode

__all__ = [
    # Core
    'MolecularScorer', 'Validator', 'setup_logger',
    # Target Discovery
    'NetworkAnalyzer', 'TargetRanker',
    # Molecule Design
    'MoleculeGenerator', 'PropertyOptimizer',
    # Docking
    'VinaWrapper', 'BindingSiteDetector',
    # ADMET
    'calculate_drug_likeness', 'predict_admet',
    # Bio-Knowledge
    'ResistanceDatabase', 'PathogenDatabase', 'TargetProteinDB', 'DatabaseUpdater',
    # Genomics
    'ResistanceGenomicsAnalyzer', 'MutationTracker', 'OmicsIntegrator',
    # Microbiology
    'MICAnalyzer', 'AssayTracker', 'StrainManager',
    # World Model
    'WorldStateTracker', 'KnowledgeGraph', 'LearningLoop',
    'KosmosEngine', 'HypothesisEngine', 'WorldModelManager',
    # Resistance / PINCER
    'PincerEngine', 'MutationSpaceMapper', 'ViableMutant',
    'ThreatAwareFitnessFunction',
    'HGTMapper', 'TransferThreatCluster', 'make_pincer_with_hgt',
    # Integrations
    'KhukuriKeyBox', 'plug_keybox_into_pincer', 'create_mrsa_keybox', 'KEYBOX_AVAILABLE',
    'NumenRetriever', 'ThreatIndex', 'LiteratureIndex', 'CompoundMemory', 'KhukuriNumen',
    # AMR Agents
    'MicrobiologyAgent', 'GenomicsAgent', 'CheminformaticsAgent',
    'ResistanceCriticAgent', 'LiteratureAgent',
    # Watchdog
    'SystemHealthWatchdog', 'ScientificHealthWatchdog', 'DrugResistanceSentinel',
    'HealthReport', 'ComponentHealth', 'ScienceHealthReport', 'ResistanceAlert',
    # Autonomous
    'AutonomousLoopController', 'PINCERScheduler', 'Mode',
]
