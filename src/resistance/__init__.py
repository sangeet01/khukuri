"""Resistance prediction, prevention, and counter-evolution (PINCER)"""

from .predictor import ResistancePredictor
from .multi_target import MultiTargetDesigner
from .evolution_simulator import EvolutionSimulator
from .mechanisms import load_resistance_mechanisms
from .pincer_engine import (
    PincerEngine,
    MutationSpaceMapper,
    MarkovMutationMatrix,
    MolecularMutator,
    DeadZoneIndex,
    ViableMutant,
    AminoAcidMutation,
    PincerCandidate,
    encode_sequence,
    decode_sequence,
    bitwise_complement,
    hamming_distance,
)

__all__ = [
    'ResistancePredictor',
    'MultiTargetDesigner',
    'EvolutionSimulator',
    'load_resistance_mechanisms',
    # PINCER engine
    'PincerEngine',
    'MutationSpaceMapper',
    'MarkovMutationMatrix',
    'MolecularMutator',
    'DeadZoneIndex',
    'ViableMutant',
    'AminoAcidMutation',
    'PincerCandidate',
    'encode_sequence',
    'decode_sequence',
    'bitwise_complement',
    'hamming_distance',
]
