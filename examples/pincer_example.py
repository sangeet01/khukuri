"""
PINCER Algorithm -- Standalone Example

Demonstrates the complete PINCER counter-evolution pipeline:
1. Red Team: Map the viable mutation space of a bacterial binding pocket
2. Blue Team: Evolve a Skeleton Key drug via Darwin-Godel minimax loop
3. Report: Top candidates with worst-case binding analysis

Usage:
    python examples/pincer_example.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.resistance import (
    PincerEngine,
    MarkovMutationMatrix,
    MutationSpaceMapper,
    encode_sequence,
    decode_sequence,
    bitwise_complement,
    hamming_distance,
)


def demo_2bit_encoding():
    """Demonstrate Durbin's 2-bit DNA encoding and O(1) complement."""
    print("=" * 60)
    print("PHASE 0: 2-Bit DNA Encoding (Durbin)")
    print("=" * 60)

    seq = "ATCGATCG"
    bits = encode_sequence(seq)
    complement = bitwise_complement(bits)

    print(f"  Original:   {seq}")
    print(f"  Bits:       {list(bits)}")
    print(f"  Complement: {decode_sequence(complement)}")
    print(f"  (Computed via single bitwise NOT operation)")
    print()


def demo_markov_matrix():
    """Demonstrate the Markov mutation matrix with transition bias."""
    print("=" * 60)
    print("PHASE 0.5: Markov Mutation Matrix")
    print("=" * 60)

    markov = MarkovMutationMatrix()
    print(f"  Mutation rate: {markov.mutation_rate_per_bp:.1e} per bp per gen")
    print(f"  Genome length: {markov.genome_length:,} bp")
    print(f"  Expected mutations/gen: {markov.mutations_per_generation:.4f}")
    print()

    print("  Substitution probability matrix (A=0, C=1, G=2, T=3):")
    labels = ['A', 'C', 'G', 'T']
    print(f"       {'    '.join(labels)}")
    for i, row in enumerate(markov.matrix):
        vals = '  '.join(f"{v:.3f}" for v in row)
        print(f"  {labels[i]}:  {vals}")

    print()
    print("  Note: A<->G (transition) has higher probability than A<->C (transversion)")
    print()


def demo_pincer():
    """Run the full PINCER pipeline."""

    # --- Configuration ---
    wild_type_seq = "AMILVCFYWHDEKRSTGNPQ"
    active_site = list(range(len(wild_type_seq)))

    # Known clinical mutations (e.g., from rpoB resistance data)
    known_mutations = {
        5: [{'mutant': 'L', 'fitness_cost': 0.15}],
        10: [{'mutant': 'N', 'fitness_cost': 0.10}],
    }

    seed_drugs = [
        "c1ccccc1",      # benzene scaffold
        "c1ccncc1",      # pyridine scaffold
        "c1ccoc1",       # furan scaffold
        "c1ccc(F)cc1",   # fluorobenzene
        "c1ccc(Cl)cc1",  # chlorobenzene
    ]

    # --- Phase 1: Red Team ---
    print("=" * 60)
    print("PHASE 1: RED TEAM -- Mapping Viable Mutation Space")
    print("=" * 60)

    pincer = PincerEngine(population_size=50, n_generations=10)
    threats = pincer.map_threats(
        wild_type_seq, active_site, known_mutations
    )

    singles = sum(1 for t in threats if len(t.mutations) == 1)
    doubles = sum(1 for t in threats if len(t.mutations) == 2)

    print(f"  Wild-type sequence:  {wild_type_seq}")
    print(f"  Active site size:    {len(active_site)} residues")
    print(f"  Total viable threats: {len(threats)}")
    print(f"    Single mutants:    {singles}")
    print(f"    Double mutants:    {doubles}")
    print(f"  Dead zones:          {len(pincer.mapper.dead_zones)}")
    print()

    if threats:
        print("  Top 5 highest-fitness threats:")
        top = sorted(threats, key=lambda t: t.fitness, reverse=True)[:5]
        for i, t in enumerate(top):
            muts = ', '.join(
                f"{m.wild_type}{m.position}{m.mutant}" for m in t.mutations
            )
            print(f"    {i+1}. {muts}  fitness={t.fitness:.3f}  dG={t.folding_energy:.2f}")
    print()

    # --- Phase 2: Blue Team ---
    print("=" * 60)
    print("PHASE 2: BLUE TEAM -- Darwin-Godel Evolutionary Loop")
    print("=" * 60)

    def log_callback(gen, apex):
        if gen % 3 == 0 or gen == pincer.n_generations - 1:
            print(f"  Gen {gen:3d}: minimax={apex.minimax_score:.4f}  {apex.smiles}")

    apex = pincer.evolve(seed_drugs, callback=log_callback)
    print()

    # --- Results ---
    print("=" * 60)
    print("RESULTS: PINCER Counter-Evolution Analysis")
    print("=" * 60)

    results = pincer.get_results()

    print(f"  Apex Skeleton Key:")
    print(f"    SMILES:       {results['apex_drug']['smiles']}")
    print(f"    Minimax score: {results['apex_drug']['minimax_score']:.4f}")
    print(f"    Worst mutant:  {results['apex_drug']['worst_mutant']}")
    print()

    print(f"  Top 5 candidates:")
    for i, c in enumerate(results['top_5']):
        print(f"    {i+1}. {c['smiles']:30s} minimax={c['minimax_score']:.4f}")
    print()

    print(f"  Statistics:")
    print(f"    Viable threats:  {results['viable_threats_count']}")
    print(f"    Dead zones:      {results['dead_zones_count']}")
    print(f"    Total explored:  {results['explored_count']}")
    print(f"    Generations:     {results['generations']}")
    print()
    print("  The Pincer snaps shut.")


if __name__ == '__main__':
    demo_2bit_encoding()
    demo_markov_matrix()
    demo_pincer()
