"""Resistance prediction and prevention"""

from .predictor import ResistancePredictor
from .multi_target import MultiTargetDesigner
from .evolution_simulator import EvolutionSimulator
from .mechanisms import load_resistance_mechanisms
from .pincer_engine import (
    PincerEngine, 
    MarkovMutationMatrix, 
    MutationSpaceMapper,
    encode_sequence,
    decode_sequence,
    bitwise_complement,
    hamming_distance
)
from .threat_fitness import ThreatAwareFitnessFunction
from .hgt_mapper import (
    HGTMapper,
    ResistanceGeneAllele,
    TransferThreatCluster,
    DualRedTeamMixin,
    make_pincer_with_hgt
)

__all__ = [
    'ResistancePredictor', 
    'MultiTargetDesigner', 
    'EvolutionSimulator', 
    'load_resistance_mechanisms',
    'PincerEngine',
    'MarkovMutationMatrix',
    'MutationSpaceMapper',
    'encode_sequence',
    'decode_sequence',
    'bitwise_complement',
    'hamming_distance',
    'ThreatAwareFitnessFunction',
    'HGTMapper',
    'ResistanceGeneAllele',
    'TransferThreatCluster',
    'DualRedTeamMixin',
    'make_pincer_with_hgt'
]
