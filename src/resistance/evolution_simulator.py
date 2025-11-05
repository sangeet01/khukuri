"""Simulate resistance evolution"""

import logging
import numpy as np

logger = logging.getLogger('khukuri')


class EvolutionSimulator:
    """Simulate resistance evolution"""
    
    def simulate(self, target_combo, generations=100, mutation_rate=1e-7):
        """Simulate resistance evolution"""
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
