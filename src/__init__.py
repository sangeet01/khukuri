"""Khukuri Virtual Lab - AMR-Focused Drug Discovery Platform"""

__version__ = "3.0.0"
__author__ = "Khukuri Research Team"

from .core import MolecularScorer, Validator, setup_logger
from .target_discovery import NetworkAnalyzer, TargetRanker
from .molecule_design import MoleculeGenerator, PropertyOptimizer
from .docking import VinaWrapper, BindingSiteDetector
from .admet import calculate_drug_likeness, predict_admet

# PINCER counter-evolution engine
from .resistance import PincerEngine, MutationSpaceMapper, MarkovMutationMatrix

# AMR-specific modules
from .bioknowledge import ResistanceDatabase, PathogenDatabase, TargetProteinDB, DatabaseUpdater
from .genomics import ResistanceGenomicsAnalyzer, MutationTracker, OmicsIntegrator
from .microbiology import MICAnalyzer, AssayTracker, StrainManager
from .world_model import WorldStateTracker, KnowledgeGraph, LearningLoop, KosmosEngine, HypothesisEngine
from .agents.amr_agents import (MicrobiologyAgent, GenomicsAgent, CheminformaticsAgent,
                                 ResistanceCriticAgent, LiteratureAgent)

__all__ = [
    # Core
    'MolecularScorer',
    'Validator',
    'setup_logger',
    # Target Discovery
    'NetworkAnalyzer',
    'TargetRanker',
    # Molecule Design
    'MoleculeGenerator',
    'PropertyOptimizer',
    # Docking
    'VinaWrapper',
    'BindingSiteDetector',
    # ADMET
    'calculate_drug_likeness',
    'predict_admet',
    # Bio-Knowledge
    'ResistanceDatabase',
    'PathogenDatabase',
    'TargetProteinDB',
    'DatabaseUpdater',
    # Genomics
    'ResistanceGenomicsAnalyzer',
    'MutationTracker',
    'OmicsIntegrator',
    # Microbiology
    'MICAnalyzer',
    'AssayTracker',
    'StrainManager',
    # World Model
    'WorldStateTracker',
    'KnowledgeGraph',
    'LearningLoop',
    'KosmosEngine',
    'HypothesisEngine',
    # AMR Agents
    'MicrobiologyAgent',
    'GenomicsAgent',
    'CheminformaticsAgent',
    'ResistanceCriticAgent',
    'LiteratureAgent',
    # PINCER Engine
    'PincerEngine',
    'MutationSpaceMapper',
    'MarkovMutationMatrix',
]
