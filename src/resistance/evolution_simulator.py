"""
Resistance evolution simulation powered by the PINCER engine.

Replaces the previous placeholder with actual co-evolutionary simulation
using the Darwin-Godel minimax loop.
"""

import logging
import numpy as np

from .pincer_engine import (
    PincerEngine,
    MutationSpaceMapper,
    MarkovMutationMatrix,
)

logger = logging.getLogger('khukuri')


class EvolutionSimulator:
    """
    Simulate resistance evolution using the PINCER algorithm.

    Red Team: enumerates viable bacterial mutations via MutationSpaceMapper.
    Blue Team: evolves drug candidates via Darwin-Godel minimax loop.

    Can operate in two modes:
    - quick: fast heuristic estimation (no docking, good for screening)
    - full:  runs the complete PINCER evolutionary loop with optional
             external fitness function (docking-based scoring)
    """

    def __init__(self, population_size=50, n_generations=20):
        self.pincer = PincerEngine(
            population_size=population_size,
            n_generations=n_generations,
        )
        self.markov = MarkovMutationMatrix()

    def simulate(self, target_combo, generations=100, mutation_rate=1e-7):
        """
        Quick simulation: estimate resistance emergence timeline.

        Uses the Markov mutation matrix and viable mutation space size
        to compute a biologically grounded estimate of when resistance
        emerges for single-target vs multi-target strategies.

        Args:
            target_combo: dict with 'targets' key listing target proteins
            generations: number of bacterial generations to simulate
            mutation_rate: per-site mutation rate

        Returns:
            dict with resistance timeline for single vs multi-target
        """
        num_targets = len(target_combo.get('targets', [1]))

        # Drake's Rule: mutations per generation
        mu = self.markov.mutations_per_generation

        # Single target: resistance = 1 mutation away
        # Expected generations = 1 / (mu * p_viable)
        # p_viable ~ fraction of mutations that are non-lethal (~0.1)
        p_viable = 0.1
        single_gen = max(1, int(1.0 / (mu * p_viable * mutation_rate * 1e7)))
        single_gen = min(single_gen, int(generations * 0.3))

        # Multi-target: resistance requires independent mutations at each target
        # Probability compounds: (mu * p_viable)^num_targets
        multi_factor = 1.0
        for _ in range(num_targets - 1):
            multi_factor *= (mu * p_viable * 0.01)
            multi_factor = max(multi_factor, 1e-10)

        multi_gen = min(
            int(single_gen / max(multi_factor, 1e-10)),
            generations
        )

        # Resistance frequency at endpoint
        single_freq = min(0.95, 1.0 - np.exp(-generations / max(single_gen, 1)))
        multi_freq = min(0.95, 1.0 - np.exp(-generations / max(multi_gen, 1)))

        return {
            'single_target': {
                'resistance_generation': single_gen,
                'resistance_frequency': round(single_freq, 4),
            },
            'multi_target': {
                'resistance_generation': multi_gen,
                'resistance_frequency': round(multi_freq, 4),
            },
            'delay_factor': round(multi_gen / single_gen, 2) if single_gen > 0 else 1,
            'mutations_per_generation': round(mu, 6),
        }

    def run_pincer(self, wild_type_seq, active_site_indices,
                   seed_smiles, known_mutations=None,
                   fitness_fn=None, callback=None):
        """
        Full PINCER run: map threats, then evolve Skeleton Key drug.

        Args:
            wild_type_seq: amino acid sequence of target protein
            active_site_indices: binding pocket residue indices
            seed_smiles: initial drug SMILES list
            known_mutations: optional clinical mutation data
            fitness_fn: optional docking-based fitness function
            callback: optional per-generation callback

        Returns:
            dict with apex drug, top candidates, threat analysis, history
        """
        logger.info("PINCER: Phase 1 -- Mapping viable threat landscape")
        threats = self.pincer.map_threats(
            wild_type_seq, active_site_indices, known_mutations
        )

        logger.info(
            f"PINCER: {len(threats)} viable threats mapped. "
            f"Phase 2 -- Evolving Skeleton Key"
        )
        apex = self.pincer.evolve(
            seed_smiles,
            fitness_fn=fitness_fn,
            callback=callback,
        )

        results = self.pincer.get_results()
        results['threat_analysis'] = {
            'total_viable': len(threats),
            'dead_zones': len(self.pincer.mapper.dead_zones),
            'single_mutants': sum(
                1 for t in threats if len(t.mutations) == 1
            ),
            'double_mutants': sum(
                1 for t in threats if len(t.mutations) == 2
            ),
        }

        logger.info(
            f"PINCER complete: apex={apex.smiles} "
            f"minimax={apex.minimax_score:.4f}"
        )
        return results
