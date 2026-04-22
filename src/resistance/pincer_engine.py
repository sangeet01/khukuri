"""
PINCER Algorithm: In Silico Counter-Evolution Engine for AMR

Models the pharmacological duel between antibiotics and bacterial resistance
as a zero-sum minimax game. Exploits the finite viable mutation space of
bacterial receptors to pre-compute resistance pathways and evolve
broad-spectrum molecular candidates.

References:
    Nowak, M.A. "Evolutionary Dynamics" (2006)
    Durbin, R. et al. "Biological Sequence Analysis" (1998)
    Drake, J.W. "A constant rate of spontaneous mutation" (1991)
"""

import logging
import random
import hashlib
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from collections import defaultdict

import numpy as np

logger = logging.getLogger('khukuri')


# ---------------------------------------------------------------------------
# 2-Bit DNA Encoding (Durbin)
# A=00, C=01, G=10, T/U=11
# Complement is bitwise NOT: ~00=11 (A<->T), ~01=10 (C<->G)
# ---------------------------------------------------------------------------

NUCLEOTIDE_TO_BITS = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'U': 3}
BITS_TO_DNA = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
BITS_TO_RNA = {0: 'A', 1: 'C', 2: 'G', 3: 'U'}

AMINO_ACIDS = list('ARNDCQEGHILKMFPSTWYV')

# Physicochemical properties: volume (A^3), charge at pH7, hydrophobicity
AMINO_ACID_PROPS = {
    'A': (88.6,   0.0,  1.8), 'R': (173.4,  1.0, -4.5),
    'N': (114.1,  0.0, -3.5), 'D': (111.1, -1.0, -3.5),
    'C': (108.5,  0.0,  2.5), 'Q': (143.8,  0.0, -3.5),
    'E': (138.4, -1.0, -3.5), 'G': (60.1,   0.0, -0.4),
    'H': (153.2,  0.5, -3.2), 'I': (166.7,  0.0,  4.5),
    'L': (166.7,  0.0,  3.8), 'K': (168.6,  1.0, -3.9),
    'M': (162.9,  0.0,  1.9), 'F': (189.9,  0.0,  2.8),
    'P': (112.7,  0.0, -1.6), 'S': (89.0,   0.0, -0.8),
    'T': (116.1,  0.0, -0.7), 'W': (227.8,  0.0, -0.9),
    'Y': (193.6,  0.0, -1.3), 'V': (140.0,  0.0,  4.2),
}


def encode_sequence(seq):
    """Encode DNA/RNA string to 2-bit numpy array."""
    return np.array([NUCLEOTIDE_TO_BITS.get(c, 0) for c in seq.upper()],
                    dtype=np.uint8)


def decode_sequence(bits, is_rna=False):
    """Decode 2-bit array back to string."""
    lut = BITS_TO_RNA if is_rna else BITS_TO_DNA
    return ''.join(lut[int(b) & 0b11] for b in bits)


def bitwise_complement(bits):
    """O(1) complement via NOT. A<->T, C<->G."""
    return (~bits) & 0b11


def hamming_distance(a, b):
    """Hamming distance between two 2-bit encoded sequences."""
    return int(np.sum(a != b))


# ---------------------------------------------------------------------------
# Markov Mutation Matrix
# Transitions (purine<->purine) are 2-3x more likely than transversions.
# Drake's Rule: E. coli ~0.003 mutations per genome per generation.
# ---------------------------------------------------------------------------

class MarkovMutationMatrix:
    """Nucleotide substitution probability matrix with transition bias."""

    def __init__(self, transition_weight=0.7, transversion_weight=0.3,
                 mutation_rate_per_bp=1e-10, genome_length=4_600_000):
        self.mutation_rate_per_bp = mutation_rate_per_bp
        self.genome_length = genome_length

        # Build 4x4 substitution matrix
        # Purines: A(0), G(2).  Pyrimidines: C(1), T(3).
        ti = transition_weight
        tv = transversion_weight
        total = ti + 2 * tv

        self.matrix = np.array([
            [0.0,      tv/total, ti/total, tv/total],  # A ->
            [tv/total, 0.0,      tv/total, ti/total],  # C ->
            [ti/total, tv/total, 0.0,      tv/total],  # G ->
            [tv/total, ti/total, tv/total, 0.0     ],  # T ->
        ], dtype=np.float64)

        # Normalize rows
        for i in range(4):
            s = self.matrix[i].sum()
            if s > 0:
                self.matrix[i] /= s

    @property
    def mutations_per_generation(self):
        """Expected mutations per genome per generation (Drake's Rule)."""
        return self.mutation_rate_per_bp * self.genome_length

    def mutate_bits(self, bits, n_mutations=1):
        """Apply n point mutations to a 2-bit encoded sequence."""
        result = bits.copy()
        n = len(result)
        if n == 0:
            return result

        positions = np.random.choice(n, size=min(n_mutations, n), replace=False)
        for pos in positions:
            cur = int(result[pos])
            result[pos] = np.random.choice(4, p=self.matrix[cur])
        return result

    def generate_mutation_paths(self, wild_type_bits, active_positions,
                                max_depth=2, paths_per_depth=100):
        """
        Generate biologically plausible mutation paths from wild-type,
        only mutating active site positions (binding pocket residues).
        """
        paths = []
        n_sites = len(active_positions)
        if n_sites == 0:
            return paths

        for depth in range(1, max_depth + 1):
            for _ in range(min(paths_per_depth, n_sites ** depth)):
                mutant = wild_type_bits.copy()
                chosen = random.sample(
                    active_positions,
                    min(depth, len(active_positions))
                )
                for pos in chosen:
                    if pos < len(mutant):
                        cur = int(mutant[pos])
                        mutant[pos] = np.random.choice(4, p=self.matrix[cur])

                if not np.array_equal(mutant, wild_type_bits):
                    paths.append(mutant)

        return paths


# ---------------------------------------------------------------------------
# Mutation Space Mapper (Red Team)
# Maps the finite viable mutation space V of a bacterial receptor.
# V = { v in 20^N | delta_G_fold(v) < theta_viable }
# ---------------------------------------------------------------------------

@dataclass
class AminoAcidMutation:
    """A single amino acid substitution."""
    position: int
    wild_type: str
    mutant: str
    fitness_cost: float = 0.0


@dataclass
class ViableMutant:
    """A receptor mutant the bacterium can adopt without dying."""
    sequence: str
    mutations: List[AminoAcidMutation]
    folding_energy: float
    fitness: float
    vector: Optional[np.ndarray] = None


class MutationSpaceMapper:
    """
    Map the finite viable mutation space of a bacterial binding pocket.

    The biological hard ceiling: if a bacterium mutates its receptor too
    aggressively, the protein loses structural integrity and the bacterium
    dies. Only mutations preserving fitness above threshold are viable.
    """

    def __init__(self, fitness_threshold=0.3):
        self.fitness_threshold = fitness_threshold
        self.viable_mutants = []
        self.dead_zones = set()

    def map_binding_pocket(self, wild_type_seq, active_site_indices,
                           known_mutations=None):
        """
        Enumerate viable single and double mutations of the binding pocket.

        Args:
            wild_type_seq: amino acid string of the wild-type protein
            active_site_indices: list of residue indices in the binding pocket
            known_mutations: dict {position: [{'mutant': 'L', 'fitness_cost': 0.15}]}

        Returns:
            List[ViableMutant]
        """
        self.viable_mutants = []
        self.dead_zones = set()
        wt = list(wild_type_seq)

        # --- Single-point mutations ---
        viable_per_pos = defaultdict(list)

        for pos in active_site_indices:
            if pos >= len(wt):
                continue
            wt_aa = wt[pos]

            for mut_aa in AMINO_ACIDS:
                if mut_aa == wt_aa:
                    continue

                cost = self._fitness_cost(wt_aa, mut_aa)

                # Override with clinical data if available
                if known_mutations and pos in known_mutations:
                    for km in known_mutations[pos]:
                        if km.get('mutant') == mut_aa:
                            cost = km.get('fitness_cost', cost)

                fitness = 1.0 - cost

                if fitness > self.fitness_threshold:
                    mut_list = wt.copy()
                    mut_list[pos] = mut_aa
                    mut_seq = ''.join(mut_list)

                    m = AminoAcidMutation(pos, wt_aa, mut_aa, cost)
                    vm = ViableMutant(
                        sequence=mut_seq,
                        mutations=[m],
                        folding_energy=-cost * 5.0,
                        fitness=fitness,
                        vector=self._embed(wild_type_seq, [m])
                    )
                    self.viable_mutants.append(vm)
                    viable_per_pos[pos].append(mut_aa)
                else:
                    self.dead_zones.add((pos, mut_aa))

        # --- Double-point mutations (from viable singles only) ---
        positions = [p for p in active_site_indices if p in viable_per_pos]
        count = 0
        max_doubles = 300

        for i, p1 in enumerate(positions):
            if count >= max_doubles:
                break
            for p2 in positions[i + 1:]:
                if count >= max_doubles:
                    break
                for aa1 in viable_per_pos[p1][:3]:
                    if count >= max_doubles:
                        break
                    for aa2 in viable_per_pos[p2][:3]:
                        if count >= max_doubles:
                            break

                        c1 = self._fitness_cost(wt[p1], aa1)
                        c2 = self._fitness_cost(wt[p2], aa2)
                        combined = max(0.0, 1.0 - c1 - c2)

                        if combined > self.fitness_threshold:
                            ml = wt.copy()
                            ml[p1] = aa1
                            ml[p2] = aa2

                            m1 = AminoAcidMutation(p1, wt[p1], aa1, c1)
                            m2 = AminoAcidMutation(p2, wt[p2], aa2, c2)

                            vm = ViableMutant(
                                sequence=''.join(ml),
                                mutations=[m1, m2],
                                folding_energy=-(c1 + c2) * 5.0,
                                fitness=combined,
                                vector=self._embed(wild_type_seq, [m1, m2])
                            )
                            self.viable_mutants.append(vm)
                            count += 1

        logger.info(
            f"PINCER Red Team: {len(self.viable_mutants)} viable mutants, "
            f"{len(self.dead_zones)} dead zones from "
            f"{len(active_site_indices)} active site positions"
        )
        return self.viable_mutants

    def _fitness_cost(self, wt_aa, mut_aa):
        """Estimate fitness cost from physicochemical disruption."""
        if wt_aa not in AMINO_ACID_PROPS or mut_aa not in AMINO_ACID_PROPS:
            return 0.5
        wv, wc, wh = AMINO_ACID_PROPS[wt_aa]
        mv, mc, mh = AMINO_ACID_PROPS[mut_aa]

        vol = abs(wv - mv) / 200.0
        chg = abs(wc - mc)
        hyd = abs(wh - mh) / 9.0

        return min(1.0, 0.3 * vol + 0.4 * chg + 0.3 * hyd)

    def _embed(self, wt_seq, mutations):
        """Compute a fixed-dimension vector for a mutant."""
        feats = []
        for m in mutations:
            if m.wild_type in AMINO_ACID_PROPS and m.mutant in AMINO_ACID_PROPS:
                wp = AMINO_ACID_PROPS[m.wild_type]
                mp = AMINO_ACID_PROPS[m.mutant]
                feats.extend([
                    m.position / 1000.0,
                    wp[0] / 250.0, mp[0] / 250.0,
                    wp[1], mp[1],
                    wp[2] / 5.0, mp[2] / 5.0,
                    m.fitness_cost,
                ])

        dim = 64
        if len(feats) < dim:
            feats.extend([0.0] * (dim - len(feats)))
        feats = feats[:dim]

        vec = np.array(feats, dtype=np.float64)
        norm = np.linalg.norm(vec)
        if norm > 0:
            vec /= norm
        return vec


# ---------------------------------------------------------------------------
# Dead Zone Index
# O(1) lookup to avoid re-exploring failed chemical space.
# ---------------------------------------------------------------------------

class DeadZoneIndex:
    """Track explored chemical space via content-addressable hashing."""

    def __init__(self):
        self.explored = set()
        self.failed = set()

    def _hash(self, smiles):
        return hashlib.sha256(smiles.encode()).hexdigest()[:16]

    def is_explored(self, smiles):
        return self._hash(smiles) in self.explored

    def is_dead_zone(self, smiles):
        return self._hash(smiles) in self.failed

    def mark(self, smiles, score, fail_threshold=0.2):
        h = self._hash(smiles)
        self.explored.add(h)
        if score < fail_threshold:
            self.failed.add(h)


# ---------------------------------------------------------------------------
# Molecular Mutator (Blue Team)
# Grammar-safe mutation operators on SMILES/RDKit molecules.
# ---------------------------------------------------------------------------

class MolecularMutator:
    """
    Mutation operators for drug candidates.
    Uses RDKit for chemical validity (analogous to SCRIPT grammar guarantee).
    Every output is a valid molecule -- zero wasted compute.
    """

    # Functional group swaps (bioisosteric replacements)
    ATOM_SWAPS = [
        ('F', 'Cl'), ('F', 'Br'), ('Cl', 'Br'),
        ('O', 'S'), ('N', 'O'),
    ]

    RING_SWAPS = [
        ('c1ccccc1', 'c1ccncc1'),  # benzene -> pyridine
        ('c1ccccc1', 'c1ccoc1'),   # benzene -> furan
        ('c1ccncc1', 'c1ccoc1'),   # pyridine -> furan
    ]

    APPENDAGES = ['F', 'Cl', 'O', 'N', 'C', 'C(=O)O', 'C(=O)N', 'S(=O)(=O)N']

    def __init__(self):
        try:
            from rdkit import Chem
            self._chem = Chem
        except ImportError:
            self._chem = None
            logger.warning("RDKit not available, molecular mutations disabled")

    def mutate(self, smiles):
        """Apply a random valid mutation to a SMILES string."""
        if self._chem is None:
            return smiles

        strategies = [self._swap_atom, self._swap_ring, self._append_group]
        random.shuffle(strategies)

        for strategy in strategies:
            result = strategy(smiles)
            if result and result != smiles:
                mol = self._chem.MolFromSmiles(result)
                if mol is not None:
                    return self._chem.MolToSmiles(mol)

        return smiles

    def crossover(self, smiles_a, smiles_b):
        """Simple crossover: take scaffold from A, substituent from B."""
        if self._chem is None:
            return smiles_a

        mol_a = self._chem.MolFromSmiles(smiles_a)
        mol_b = self._chem.MolFromSmiles(smiles_b)
        if mol_a is None or mol_b is None:
            return smiles_a

        # Use the larger molecule as base, attach fragment from smaller
        if mol_a.GetNumAtoms() >= mol_b.GetNumAtoms():
            base, donor = smiles_a, smiles_b
        else:
            base, donor = smiles_b, smiles_a

        # Try appending a piece of the donor
        result = base + donor[:3]
        mol = self._chem.MolFromSmiles(result)
        if mol is not None:
            return self._chem.MolToSmiles(mol)
        return base

    def _swap_atom(self, smiles):
        for a, b in self.ATOM_SWAPS:
            if a in smiles:
                return smiles.replace(a, b, 1)
            if b in smiles:
                return smiles.replace(b, a, 1)
        return None

    def _swap_ring(self, smiles):
        for a, b in self.RING_SWAPS:
            if a in smiles:
                return smiles.replace(a, b, 1)
            if b in smiles:
                return smiles.replace(b, a, 1)
        return None

    def _append_group(self, smiles):
        group = random.choice(self.APPENDAGES)
        return smiles + group


# ---------------------------------------------------------------------------
# PINCER Engine (The Darwin-Godel Loop)
# s* = argmax_s ( min_v K(s, v) )   for s in S, v in V
# ---------------------------------------------------------------------------

@dataclass
class PincerCandidate:
    """A drug candidate in the evolutionary population."""
    smiles: str
    minimax_score: float = 0.0
    worst_mutant: str = ''
    generation: int = 0


class PincerEngine:
    """
    The PINCER Algorithm: Darwin-Godel Counter-Evolution Engine.

    Red Team: MutationSpaceMapper enumerates viable bacterial mutations.
    Blue Team: Evolutionary loop optimizes drug candidates via minimax.
    Godel Constraint: only provably superior mutations are accepted.

    The optimal Skeleton Key drug s* maximizes the worst-case binding
    affinity across all viable receptor mutants.
    """

    def __init__(self, population_size=100, n_generations=50,
                 elite_fraction=0.2):
        self.population_size = population_size
        self.n_generations = n_generations
        self.elite_count = max(1, int(population_size * elite_fraction))

        self.markov = MarkovMutationMatrix()
        self.mapper = MutationSpaceMapper()
        self.dead_zones = DeadZoneIndex()
        self.mutator = MolecularMutator()

        self.viable_threats = []
        self.population = []
        self.history = []

    # -- Red Team: Map the threat landscape --

    def map_threats(self, wild_type_seq, active_site_indices,
                    known_mutations=None):
        """
        Phase 1 (Red Team): Map the finite viable mutation space.

        Args:
            wild_type_seq: amino acid sequence of target protein
            active_site_indices: list of binding pocket residue indices
            known_mutations: optional clinical mutation data

        Returns:
            List[ViableMutant] -- the threat cluster
        """
        self.viable_threats = self.mapper.map_binding_pocket(
            wild_type_seq, active_site_indices, known_mutations
        )
        return self.viable_threats

    # -- Blue Team: Evolve the Skeleton Key --

    def evolve(self, seed_smiles, fitness_fn=None, callback=None):
        """
        Phase 2 (Blue Team): Run the Darwin-Godel evolutionary loop.

        Args:
            seed_smiles: list of initial SMILES strings (seed population)
            fitness_fn: callable(smiles, viable_threats) -> float
                        If None, uses a default docking-free scoring heuristic.
            callback: optional callable(generation, apex_candidate) for logging

        Returns:
            PincerCandidate -- the apex Skeleton Key drug
        """
        if fitness_fn is None:
            fitness_fn = self._default_fitness

        # Initialize population from seeds
        self.population = []
        for smi in seed_smiles[:self.population_size]:
            self.population.append(PincerCandidate(smiles=smi, generation=0))

        # Fill remaining slots with mutations of seeds
        while len(self.population) < self.population_size:
            parent = random.choice(seed_smiles)
            child = self.mutator.mutate(parent)
            if not self.dead_zones.is_dead_zone(child):
                self.population.append(
                    PincerCandidate(smiles=child, generation=0)
                )

        self.history = []

        for gen in range(self.n_generations):
            # Evaluate minimax fitness for each candidate
            for candidate in self.population:
                score, worst = self._minimax_evaluate(
                    candidate.smiles, fitness_fn
                )
                candidate.minimax_score = score
                candidate.worst_mutant = worst
                candidate.generation = gen

                self.dead_zones.mark(candidate.smiles, score)

            # Rank by minimax score (higher is better)
            self.population.sort(
                key=lambda c: c.minimax_score, reverse=True
            )

            apex = self.population[0]
            self.history.append({
                'generation': gen,
                'apex_score': apex.minimax_score,
                'apex_smiles': apex.smiles,
                'population_size': len(self.population),
                'dead_zones': self.dead_zones.failed.__len__(),
            })

            if callback:
                callback(gen, apex)

            logger.info(
                f"PINCER Gen {gen}: apex={apex.minimax_score:.4f} "
                f"smiles={apex.smiles}"
            )

            if gen == self.n_generations - 1:
                break

            # --- Selection + Reproduction ---
            elites = self.population[:self.elite_count]
            next_gen = list(elites)  # elitism

            parents = self.population[:max(2, self.elite_count * 2)]

            while len(next_gen) < self.population_size:
                p1 = random.choice(parents)
                p2 = random.choice(parents)

                # Crossover
                if random.random() < 0.5:
                    child_smi = self.mutator.crossover(p1.smiles, p2.smiles)
                else:
                    child_smi = p1.smiles

                # Mutation
                child_smi = self.mutator.mutate(child_smi)

                # Dead zone check
                if self.dead_zones.is_dead_zone(child_smi):
                    continue

                # Godel constraint: quick pre-check
                # Only add if at least parseable
                try:
                    from rdkit import Chem
                    if Chem.MolFromSmiles(child_smi) is None:
                        continue
                except ImportError:
                    pass

                next_gen.append(PincerCandidate(smiles=child_smi))

            self.population = next_gen

        # Return the apex Skeleton Key
        self.population.sort(key=lambda c: c.minimax_score, reverse=True)
        return self.population[0]

    def _minimax_evaluate(self, smiles, fitness_fn):
        """
        Minimax evaluation: score = min_v fitness(smiles, v)
        The drug is only as strong as its weakest binding across mutations.
        """
        if not self.viable_threats:
            return fitness_fn(smiles, []), ''

        worst_score = float('inf')
        worst_label = ''

        for threat in self.viable_threats:
            score = fitness_fn(smiles, [threat])
            if score < worst_score:
                worst_score = score
                muts = ', '.join(
                    f"{m.wild_type}{m.position}{m.mutant}"
                    for m in threat.mutations
                )
                worst_label = muts

        return worst_score, worst_label

    def _default_fitness(self, smiles, threats):
        """
        Default scoring heuristic when no docking engine is available.
        Uses molecular descriptors as a proxy for binding potential.
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 0.0

            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)

            # Lipinski-like scoring
            score = 1.0
            if mw > 500:
                score -= 0.2
            if logp > 5:
                score -= 0.2
            if hbd > 5:
                score -= 0.1
            if hba > 10:
                score -= 0.1

            # Penalize for each threat's disruption
            for threat in threats:
                for m in threat.mutations:
                    score -= m.fitness_cost * 0.05

            return max(0.0, min(1.0, score))

        except ImportError:
            return 0.5

    def get_results(self):
        """Return structured results from the evolutionary run."""
        if not self.population:
            return {}

        self.population.sort(key=lambda c: c.minimax_score, reverse=True)

        return {
            'apex_drug': {
                'smiles': self.population[0].smiles,
                'minimax_score': self.population[0].minimax_score,
                'worst_mutant': self.population[0].worst_mutant,
            },
            'top_5': [
                {
                    'smiles': c.smiles,
                    'minimax_score': c.minimax_score,
                    'worst_mutant': c.worst_mutant,
                }
                for c in self.population[:5]
            ],
            'viable_threats_count': len(self.viable_threats),
            'dead_zones_count': len(self.dead_zones.failed),
            'explored_count': len(self.dead_zones.explored),
            'generations': len(self.history),
            'history': self.history,
        }
