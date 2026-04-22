"""Simulate resistance evolution and run PINCER counter-evolution"""

import logging
import numpy as np
from .pincer_engine import PincerEngine
from .threat_fitness import ThreatAwareFitnessFunction

logger = logging.getLogger('khukuri')


class EvolutionSimulator:
    """Simulate resistance evolution and run PINCER counter-evolution"""
    
    def __init__(self, population_size=100, n_generations=50):
        self.pincer = PincerEngine(
            population_size=population_size, 
            n_generations=n_generations
        )
        self.fitness_fn = ThreatAwareFitnessFunction()
    
    def simulate(self, target_combo, generations=100, mutation_rate=1e-7):
        """Standard resistance evolution simulation"""
        num_targets = len(target_combo.get('targets', [1]))
        
        # Single target: resistance emerges early
        single_gen = int(generations * 0.3)
        
        # Multi-target: resistance emerges much later
        multi_gen = int(generations * (0.3 + 0.6 * (num_targets - 1) / 3))
        
        return {
            'single_target': {
                'resistance_generation': single_gen,
                'resistance_frequency': 0.8
            },
            'multi_target': {
                'resistance_generation': multi_gen,
                'resistance_frequency': 0.1
            },
            'delay_factor': multi_gen / single_gen if single_gen > 0 else 1
        }

    def run_pincer(self, wild_type_seq, active_site_indices, seed_smiles, known_mutations=None):
        """
        Run the PINCER counter-evolution engine.
        
        1. Red Team: Map viable mutation space
        2. Blue Team: Evolve Skeleton Key drug via minimax
        """
        logger.info("Starting PINCER Counter-Evolution run")
        
        # Phase 1: Red Team
        threats = self.pincer.map_threats(
            wild_type_seq, active_site_indices, known_mutations
        )
        
        # Phase 2: Blue Team
        apex = self.pincer.evolve(seed_smiles, fitness_fn=self.fitness_fn)
        
        results = self.pincer.get_results()
        logger.info(f"PINCER run complete. Apex minimax score: {results['apex_drug']['minimax_score']:.4f}")
        
        return results
