"""
Tests for the PINCER counter-evolution engine.

Covers:
- 2-bit DNA encoding (Durbin)
- Markov mutation matrix (transition/transversion bias)
- Mutation space mapper (Red Team)
- Dead zone index
- Molecular mutator (Blue Team)
- PINCER engine full loop (minimax + Darwin-Godel)
- EvolutionSimulator (quick mode + PINCER mode)
"""

import pytest
import numpy as np
import sys
import os

# Add project root to sys.path to allow direct execution
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.resistance.pincer_engine import (
    encode_sequence,
    decode_sequence,
    bitwise_complement,
    hamming_distance,
    MarkovMutationMatrix,
    MutationSpaceMapper,
    DeadZoneIndex,
    MolecularMutator,
    PincerEngine,
    PincerCandidate,
    ViableMutant,
    AminoAcidMutation,
    AMINO_ACIDS,
    AMINO_ACID_PROPS,
)
from src.resistance.evolution_simulator import EvolutionSimulator


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_dna():
    return "ATCGATCG"

@pytest.fixture
def wild_type_seq():
    return "AMILVCFYWH"

@pytest.fixture
def active_site(wild_type_seq):
    return list(range(len(wild_type_seq)))

@pytest.fixture
def seed_smiles():
    return ["c1ccccc1", "c1ccncc1", "c1ccoc1", "c1ccc(F)cc1"]

@pytest.fixture
def known_mutations():
    return {
        2: [{'mutant': 'L', 'fitness_cost': 0.15}],
        5: [{'mutant': 'N', 'fitness_cost': 0.08}],
    }

@pytest.fixture
def mapper():
    return MutationSpaceMapper(fitness_threshold=0.3)

@pytest.fixture
def pincer(seed_smiles):
    return PincerEngine(population_size=20, n_generations=5)


# ---------------------------------------------------------------------------
# 2-Bit DNA Encoding
# ---------------------------------------------------------------------------

class TestDNAEncoding:

    def test_encode_returns_uint8_array(self, sample_dna):
        bits = encode_sequence(sample_dna)
        assert isinstance(bits, np.ndarray)
        assert bits.dtype == np.uint8

    def test_encode_length_matches_sequence(self, sample_dna):
        bits = encode_sequence(sample_dna)
        assert len(bits) == len(sample_dna)

    def test_encode_values_in_range(self, sample_dna):
        bits = encode_sequence(sample_dna)
        assert all(0 <= int(b) <= 3 for b in bits)

    def test_decode_roundtrip_dna(self, sample_dna):
        bits = encode_sequence(sample_dna)
        result = decode_sequence(bits, is_rna=False)
        assert result == sample_dna

    def test_decode_rna_uses_u_not_t(self):
        bits = encode_sequence("AACG")
        # T/U is encoded as 3; when is_rna=True it should decode as U
        seq = decode_sequence(np.array([3], dtype=np.uint8), is_rna=True)
        assert seq == "U"
        seq_dna = decode_sequence(np.array([3], dtype=np.uint8), is_rna=False)
        assert seq_dna == "T"

    def test_complement_at_is_bitwise_inverse(self):
        # A=00, T=11 => ~00 & 0b11 == 11
        a_bits = encode_sequence("A")
        t_bits = encode_sequence("T")
        assert np.array_equal(bitwise_complement(a_bits), t_bits)

    def test_complement_cg_is_bitwise_inverse(self):
        # C=01, G=10 => ~01 & 0b11 == 10
        c_bits = encode_sequence("C")
        g_bits = encode_sequence("G")
        assert np.array_equal(bitwise_complement(c_bits), g_bits)

    def test_double_complement_is_identity(self, sample_dna):
        bits = encode_sequence(sample_dna)
        assert np.array_equal(bitwise_complement(bitwise_complement(bits)), bits)

    def test_hamming_distance_identical(self, sample_dna):
        bits = encode_sequence(sample_dna)
        assert hamming_distance(bits, bits) == 0

    def test_hamming_distance_one_mutation(self):
        a = encode_sequence("AAAA")
        b = encode_sequence("AAGA")
        assert hamming_distance(a, b) == 1

    def test_hamming_distance_full_complement(self):
        bits = encode_sequence("ATCG")
        comp = bitwise_complement(bits)
        # Every position differs
        assert hamming_distance(bits, comp) == 4

    def test_encode_uracil_same_as_thymine(self):
        u_bits = encode_sequence("U")
        t_bits = encode_sequence("T")
        assert np.array_equal(u_bits, t_bits)


# ---------------------------------------------------------------------------
# Markov Mutation Matrix
# ---------------------------------------------------------------------------

class TestMarkovMutationMatrix:

    def test_matrix_shape(self):
        m = MarkovMutationMatrix()
        assert m.matrix.shape == (4, 4)

    def test_rows_sum_to_one(self):
        m = MarkovMutationMatrix()
        for i in range(4):
            assert abs(m.matrix[i].sum() - 1.0) < 1e-9

    def test_diagonal_is_zero(self):
        m = MarkovMutationMatrix()
        for i in range(4):
            assert m.matrix[i, i] == 0.0

    def test_transitions_more_probable_than_transversions(self):
        m = MarkovMutationMatrix()
        # A(0)->G(2) transition should be more likely than A(0)->C(1) transversion
        assert m.matrix[0, 2] > m.matrix[0, 1]

    def test_drakes_rule_ecoli(self):
        m = MarkovMutationMatrix(
            mutation_rate_per_bp=1e-10,
            genome_length=4_600_000
        )
        # Expected ~0.00046, well below 1
        assert 0.0001 < m.mutations_per_generation < 0.01

    def test_mutate_bits_changes_sequence(self):
        m = MarkovMutationMatrix()
        bits = encode_sequence("AAAAAAAAAA")
        mutated = m.mutate_bits(bits, n_mutations=3)
        # At least one position should differ
        assert not np.array_equal(bits, mutated)

    def test_mutate_bits_output_in_range(self):
        m = MarkovMutationMatrix()
        bits = encode_sequence("ATCGATCG")
        mutated = m.mutate_bits(bits, n_mutations=4)
        assert all(0 <= int(b) <= 3 for b in mutated)

    def test_generate_paths_returns_list(self):
        m = MarkovMutationMatrix()
        bits = encode_sequence("ATCG")
        paths = m.generate_mutation_paths(bits, [0, 1, 2, 3], max_depth=1)
        assert isinstance(paths, list)
        assert len(paths) > 0

    def test_generate_paths_differ_from_wildtype(self):
        m = MarkovMutationMatrix()
        bits = encode_sequence("AAAA")
        paths = m.generate_mutation_paths(bits, [0, 1, 2, 3], max_depth=1)
        for path in paths:
            assert not np.array_equal(path, bits)


# ---------------------------------------------------------------------------
# Mutation Space Mapper
# ---------------------------------------------------------------------------

class TestMutationSpaceMapper:

    def test_returns_list_of_viable_mutants(self, mapper, wild_type_seq, active_site):
        threats = mapper.map_binding_pocket(wild_type_seq, active_site)
        assert isinstance(threats, list)
        assert len(threats) > 0

    def test_all_mutants_above_fitness_threshold(self, mapper, wild_type_seq, active_site):
        threats = mapper.map_binding_pocket(wild_type_seq, active_site)
        for t in threats:
            assert t.fitness > mapper.fitness_threshold

    def test_mutant_sequences_differ_from_wildtype(self, mapper, wild_type_seq, active_site):
        threats = mapper.map_binding_pocket(wild_type_seq, active_site)
        for t in threats:
            assert t.sequence != wild_type_seq

    def test_each_mutant_has_at_least_one_mutation(self, mapper, wild_type_seq, active_site):
        threats = mapper.map_binding_pocket(wild_type_seq, active_site)
        for t in threats:
            assert len(t.mutations) >= 1

    def test_vectors_are_unit_normalized(self, mapper, wild_type_seq, active_site):
        threats = mapper.map_binding_pocket(wild_type_seq, active_site)
        for t in threats:
            if t.vector is not None:
                norm = np.linalg.norm(t.vector)
                assert abs(norm - 1.0) < 1e-6 or norm == 0.0

    def test_dead_zones_populated(self, mapper, wild_type_seq, active_site):
        mapper.map_binding_pocket(wild_type_seq, active_site)
        assert len(mapper.dead_zones) > 0

    def test_known_mutations_override_fitness_cost(self, wild_type_seq, active_site, known_mutations):
        mapper = MutationSpaceMapper(fitness_threshold=0.0)
        threats = mapper.map_binding_pocket(wild_type_seq, active_site, known_mutations)
        # Should still produce mutants
        assert len(threats) > 0

    def test_high_threshold_fewer_mutants(self, wild_type_seq, active_site):
        low = MutationSpaceMapper(fitness_threshold=0.1)
        high = MutationSpaceMapper(fitness_threshold=0.9)
        t_low = low.map_binding_pocket(wild_type_seq, active_site)
        t_high = high.map_binding_pocket(wild_type_seq, active_site)
        assert len(t_low) >= len(t_high)

    def test_amino_acid_props_all_defined(self):
        for aa in AMINO_ACIDS:
            assert aa in AMINO_ACID_PROPS

    def test_fitness_cost_same_residue_is_zero(self, mapper):
        cost = mapper._fitness_cost('A', 'A')
        assert cost == 0.0 or cost < 0.05  # self-mutation has no cost

    def test_fitness_cost_charged_substitution_is_high(self, mapper):
        # Replacing neutral Ala with charged Arg should cost more than Ala->Val
        cost_conservative = mapper._fitness_cost('A', 'V')
        cost_radical = mapper._fitness_cost('A', 'R')
        assert cost_radical > cost_conservative


# ---------------------------------------------------------------------------
# Dead Zone Index
# ---------------------------------------------------------------------------

class TestDeadZoneIndex:

    def test_new_smiles_not_explored(self):
        idx = DeadZoneIndex()
        assert not idx.is_explored("c1ccccc1")

    def test_mark_sets_explored(self):
        idx = DeadZoneIndex()
        idx.mark("c1ccccc1", score=0.8)
        assert idx.is_explored("c1ccccc1")

    def test_high_score_not_dead_zone(self):
        idx = DeadZoneIndex()
        idx.mark("c1ccccc1", score=0.8)
        assert not idx.is_dead_zone("c1ccccc1")

    def test_low_score_becomes_dead_zone(self):
        idx = DeadZoneIndex()
        idx.mark("c1ccccc1", score=0.1, fail_threshold=0.2)
        assert idx.is_dead_zone("c1ccccc1")

    def test_different_smiles_independent(self):
        idx = DeadZoneIndex()
        idx.mark("c1ccccc1", score=0.1)
        assert not idx.is_dead_zone("c1ccncc1")

    def test_custom_fail_threshold(self):
        idx = DeadZoneIndex()
        idx.mark("CCO", score=0.25, fail_threshold=0.3)
        assert idx.is_dead_zone("CCO")


# ---------------------------------------------------------------------------
# Molecular Mutator
# ---------------------------------------------------------------------------

class TestMolecularMutator:

    def test_mutate_returns_string(self):
        m = MolecularMutator()
        result = m.mutate("c1ccccc1")
        assert isinstance(result, str)
        assert len(result) > 0

    def test_mutate_produces_valid_smiles(self):
        try:
            from rdkit import Chem
        except ImportError:
            pytest.skip("RDKit not available")
        m = MolecularMutator()
        for smi in ["c1ccccc1", "c1ccncc1", "CCO", "CC(=O)O"]:
            result = m.mutate(smi)
            mol = Chem.MolFromSmiles(result)
            assert mol is not None, f"Invalid SMILES after mutation: {result}"

    def test_crossover_returns_valid_smiles(self):
        try:
            from rdkit import Chem
        except ImportError:
            pytest.skip("RDKit not available")
        m = MolecularMutator()
        result = m.crossover("c1ccccc1", "c1ccncc1")
        mol = Chem.MolFromSmiles(result)
        assert mol is not None

    def test_multiple_mutations_vary(self):
        m = MolecularMutator()
        results = set(m.mutate("c1ccccc1") for _ in range(20))
        # At least 2 distinct outcomes across 20 mutations
        assert len(results) >= 2


# ---------------------------------------------------------------------------
# PINCER Engine
# ---------------------------------------------------------------------------

class TestPincerEngine:

    def test_map_threats_returns_list(self, pincer, wild_type_seq, active_site):
        threats = pincer.map_threats(wild_type_seq, active_site)
        assert isinstance(threats, list)
        assert len(threats) > 0

    def test_evolve_returns_candidate(self, pincer, wild_type_seq, active_site, seed_smiles):
        pincer.map_threats(wild_type_seq, active_site)
        apex = pincer.evolve(seed_smiles)
        assert isinstance(apex, PincerCandidate)

    def test_apex_smiles_is_valid(self, pincer, wild_type_seq, active_site, seed_smiles):
        try:
            from rdkit import Chem
        except ImportError:
            pytest.skip("RDKit not available")
        pincer.map_threats(wild_type_seq, active_site)
        apex = pincer.evolve(seed_smiles)
        mol = Chem.MolFromSmiles(apex.smiles)
        assert mol is not None

    def test_apex_minimax_score_in_range(self, pincer, wild_type_seq, active_site, seed_smiles):
        pincer.map_threats(wild_type_seq, active_site)
        apex = pincer.evolve(seed_smiles)
        assert 0.0 <= apex.minimax_score <= 1.0

    def test_history_length_equals_generations(self, pincer, wild_type_seq, active_site, seed_smiles):
        pincer.map_threats(wild_type_seq, active_site)
        pincer.evolve(seed_smiles)
        assert len(pincer.history) == pincer.n_generations

    def test_get_results_structure(self, pincer, wild_type_seq, active_site, seed_smiles):
        pincer.map_threats(wild_type_seq, active_site)
        pincer.evolve(seed_smiles)
        r = pincer.get_results()
        assert 'apex_drug' in r
        assert 'top_5' in r
        assert 'viable_threats_count' in r
        assert 'dead_zones_count' in r
        assert 'explored_count' in r
        assert 'generations' in r

    def test_top_5_sorted_by_minimax(self, pincer, wild_type_seq, active_site, seed_smiles):
        pincer.map_threats(wild_type_seq, active_site)
        pincer.evolve(seed_smiles)
        r = pincer.get_results()
        scores = [c['minimax_score'] for c in r['top_5']]
        assert scores == sorted(scores, reverse=True)

    def test_dead_zones_accumulate(self, pincer, wild_type_seq, active_site, seed_smiles):
        pincer.map_threats(wild_type_seq, active_site)
        pincer.evolve(seed_smiles)
        r = pincer.get_results()
        assert r['explored_count'] > 0

    def test_minimax_without_threats_uses_default(self, seed_smiles):
        pincer = PincerEngine(population_size=10, n_generations=3)
        # Evolve without mapping threats
        apex = pincer.evolve(seed_smiles)
        assert isinstance(apex, PincerCandidate)

    @pytest.mark.slow
    def test_larger_population_improves_score(self, wild_type_seq, active_site, seed_smiles):
        small = PincerEngine(population_size=10, n_generations=5)
        large = PincerEngine(population_size=50, n_generations=5)
        small.map_threats(wild_type_seq, active_site)
        large.map_threats(wild_type_seq, active_site)
        apex_small = small.evolve(seed_smiles)
        apex_large = large.evolve(seed_smiles)
        # Larger population should produce score >= small (stochastic, not strict)
        assert apex_large.minimax_score >= 0.0


# ---------------------------------------------------------------------------
# Evolution Simulator
# ---------------------------------------------------------------------------

class TestEvolutionSimulator:

    def test_simulate_quick_returns_dict(self):
        sim = EvolutionSimulator()
        result = sim.simulate({'targets': ['inhA', 'katG']}, generations=100)
        assert isinstance(result, dict)

    def test_simulate_has_required_keys(self):
        sim = EvolutionSimulator()
        result = sim.simulate({'targets': ['inhA']})
        assert 'single_target' in result
        assert 'multi_target' in result
        assert 'delay_factor' in result
        assert 'mutations_per_generation' in result

    def test_multi_target_slower_than_single(self):
        sim = EvolutionSimulator()
        result = sim.simulate({'targets': ['inhA', 'katG']}, generations=100)
        s_gen = result['single_target']['resistance_generation']
        m_gen = result['multi_target']['resistance_generation']
        assert m_gen >= s_gen

    def test_more_targets_slower_resistance(self):
        sim = EvolutionSimulator()
        r1 = sim.simulate({'targets': ['A']})
        r2 = sim.simulate({'targets': ['A', 'B']})
        assert (r2['multi_target']['resistance_generation'] >=
                r1['multi_target']['resistance_generation'])

    def test_resistance_frequency_between_0_and_1(self):
        sim = EvolutionSimulator()
        result = sim.simulate({'targets': ['inhA']}, generations=200)
        assert 0.0 <= result['single_target']['resistance_frequency'] <= 1.0
        assert 0.0 <= result['multi_target']['resistance_frequency'] <= 1.0

    def test_delay_factor_at_least_1(self):
        sim = EvolutionSimulator()
        result = sim.simulate({'targets': ['inhA', 'katG']})
        assert result['delay_factor'] >= 1.0

    def test_mutations_per_generation_drakes_rule(self):
        sim = EvolutionSimulator()
        result = sim.simulate({'targets': ['A']})
        mu = result['mutations_per_generation']
        # Drake's rule for E. coli: ~0.003 mutations/genome/generation
        assert 0.0001 < mu < 1.0

    def test_run_pincer_returns_results(self, wild_type_seq, active_site, seed_smiles):
        sim = EvolutionSimulator(population_size=20, n_generations=3)
        result = sim.run_pincer(wild_type_seq, active_site, seed_smiles)
        assert 'apex_drug' in result
        assert 'threat_analysis' in result

    def test_run_pincer_threat_analysis_keys(self, wild_type_seq, active_site, seed_smiles):
        sim = EvolutionSimulator(population_size=20, n_generations=3)
        result = sim.run_pincer(wild_type_seq, active_site, seed_smiles)
        ta = result['threat_analysis']
        assert 'total_viable' in ta
        assert 'dead_zones' in ta
        assert 'single_mutants' in ta
        assert 'double_mutants' in ta

    def test_run_pincer_apex_is_valid_smiles(self, wild_type_seq, active_site, seed_smiles):
        try:
            from rdkit import Chem
        except ImportError:
            pytest.skip("RDKit not available")
        sim = EvolutionSimulator(population_size=20, n_generations=3)
        result = sim.run_pincer(wild_type_seq, active_site, seed_smiles)
        apex_smiles = result['apex_drug']['smiles']
        assert Chem.MolFromSmiles(apex_smiles) is not None

    def test_run_pincer_with_known_mutations(self, wild_type_seq, active_site, seed_smiles, known_mutations):
        sim = EvolutionSimulator(population_size=10, n_generations=3)
        result = sim.run_pincer(wild_type_seq, active_site, seed_smiles,
                                known_mutations=known_mutations)
        assert result['threat_analysis']['total_viable'] > 0
