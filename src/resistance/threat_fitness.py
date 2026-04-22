"""
ThreatAwareFitnessFunction — bridges PINCER's red and blue teams.

The core problem: _default_fitness scores all seed drugs ~0.965 because
Lipinski penalties are small and threat penalization (fitness_cost * 0.05)
is noise. The evolutionary loop has nothing to climb.

This module fixes that by:

1. Embedding each drug candidate into the same 64-dim physicochemical space
   as the threat vectors (ViableMutant.vector).

2. For each viable threat, computing a *binding disruption score* — how much
   that mutation would degrade the drug's binding character.

3. Minimax fitness = min over threats of (1 - disruption). A drug that
   survives every mutation scores near 1.0. A drug wiped out by one mutation
   scores near 0.0. This creates sharp differentiation.

4. A secondary Lipinski/drug-likeness gate is applied after minimax so
   the search stays in drug-like chemical space.

Design notes:
- The drug embedding mirrors MutationSpaceMapper._embed() so both vectors
  live in the same feature space and dot-product similarity is meaningful.
- No docking required. This is a pure physicochemical surrogate intended
  as a drop-in until KeyBox is available.
- The function is stateless and picklable — safe for ParallelExecutor.
"""

import logging
from typing import List, Optional

import numpy as np

logger = logging.getLogger('khukuri')

# Physicochemical properties imported from pincer_engine to stay in sync
from .pincer_engine import AMINO_ACID_PROPS, ViableMutant


# ---------------------------------------------------------------------------
# Drug embedding
# ---------------------------------------------------------------------------

# Molecular descriptor targets per amino acid property dimension.
# We map each drug descriptor to the same axes used by _embed():
#   axis 0: position-normalised (skip for drug — no position)
#   axis 1-2: volume proxies  → MW, TPSA
#   axis 3-4: charge proxies  → HBD, HBA (net charge character)
#   axis 5-6: hydrophobicity  → LogP, aromatic fraction
#   axis 7:   fitness cost proxy → rule-of-five violations / 5

def embed_drug(smiles: str, dim: int = 64) -> Optional[np.ndarray]:
    """
    Embed a drug SMILES into the same 64-dim physicochemical space
    as ViableMutant.vector so cosine similarity is meaningful.

    Returns None if the molecule is invalid.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mw    = Descriptors.MolWt(mol)
        logp  = Descriptors.MolLogP(mol)
        tpsa  = Descriptors.TPSA(mol)
        hbd   = Descriptors.NumHDonors(mol)
        hba   = Descriptors.NumHAcceptors(mol)
        rotb  = rdMolDescriptors.CalcNumRotatableBonds(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        arom  = rdMolDescriptors.CalcNumAromaticRings(mol)
        n_atoms = mol.GetNumAtoms()
        n_heavy = mol.GetNumHeavyAtoms()

        # Fraction of atoms that are polar (N or O)
        polar = sum(
            1 for a in mol.GetAtoms()
            if a.GetAtomicNum() in (7, 8)
        ) / max(n_heavy, 1)

        # Aromatic fraction
        arom_frac = arom / max(rings, 1) if rings > 0 else 0.0

        # Rule-of-five violation count (0–4)
        ro5_violations = sum([
            mw > 500,
            logp > 5,
            hbd > 5,
            hba > 10,
        ])

        # Build base feature vector (mirrors threat vector axes)
        # Normalise to [0, 1] or [-1, 1] to match threat embedding scale
        feats = [
            # Volume proxies (axes 1-2 in threat vector)
            min(mw / 600.0, 1.0),
            min(tpsa / 150.0, 1.0),
            # Charge proxies (axes 3-4)
            min(hbd / 5.0, 1.0),
            min(hba / 10.0, 1.0),
            # Hydrophobicity (axes 5-6)
            np.clip(logp / 5.0, -1.0, 1.0),
            arom_frac,
            # Fitness proxy (axis 7)
            ro5_violations / 4.0,
            # Extra descriptors to fill 64 dims
            min(rotb / 10.0, 1.0),
            min(rings / 6.0, 1.0),
            min(n_heavy / 40.0, 1.0),
            polar,
        ]

        # Tile the base features to fill dim (repeating with slight noise
        # so the vector has non-trivial structure across all 64 dims)
        rng = np.random.default_rng(hash(smiles) % (2**31))
        base = np.array(feats, dtype=np.float64)
        tiles = int(np.ceil(dim / len(feats)))
        full = np.tile(base, tiles)[:dim]
        # Add tiny deterministic perturbation so identical scaffolds
        # with different substituents get distinct vectors
        full += rng.normal(0, 0.01, dim)

        norm = np.linalg.norm(full)
        if norm > 0:
            full /= norm
        return full

    except Exception as exc:
        logger.debug(f"embed_drug failed for {smiles}: {exc}")
        return None


# ---------------------------------------------------------------------------
# Binding disruption model
# ---------------------------------------------------------------------------

def binding_disruption(drug_vec: np.ndarray, threat: ViableMutant) -> float:
    """
    Estimate how much a receptor mutation disrupts binding of this drug.

    Uses cosine similarity between the drug embedding and the threat vector.
    High similarity → the drug's physicochemical character matches the
    mutation's profile → MORE disruption (the mutation exploits drug features).
    Low similarity → the drug is orthogonal to this mutation → less disruption.

    Returns a disruption score in [0, 1].
    """
    if threat.vector is None or drug_vec is None:
        # Fallback: use raw fitness cost of mutations
        avg_cost = np.mean([m.fitness_cost for m in threat.mutations]) if threat.mutations else 0.5
        return avg_cost * 0.5

    tv = threat.vector
    # Align dimensions
    min_dim = min(len(drug_vec), len(tv))
    d = drug_vec[:min_dim]
    t = tv[:min_dim]

    cosine = float(np.dot(d, t))   # both normalised, so this is cosine sim

    # Transform: high cosine (drug looks like mutation target) → high disruption
    # Scale by threat fitness so high-fitness mutations (dangerous ones) matter more
    disruption = ((cosine + 1.0) / 2.0) * threat.fitness

    return float(np.clip(disruption, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Drug-likeness gate
# ---------------------------------------------------------------------------

def drug_likeness_score(smiles: str) -> float:
    """
    Return a [0, 1] drug-likeness score.
    1.0 = perfect Lipinski compliance + good QED-like properties.
    0.0 = completely non-drug-like.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, QED

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0

        qed = QED.qed(mol)           # [0, 1] — already normalised
        mw  = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)

        # Hard penalties
        score = qed
        if mw > 600:
            score *= 0.5
        if mw < 100:
            score *= 0.7
        if logp > 6:
            score *= 0.6
        if logp < -3:
            score *= 0.8

        return float(np.clip(score, 0.0, 1.0))

    except Exception:
        return 0.5


# ---------------------------------------------------------------------------
# Main fitness function
# ---------------------------------------------------------------------------

class ThreatAwareFitnessFunction:
    """
    Callable fitness function for PincerEngine.evolve().

    Usage:
        fitness_fn = ThreatAwareFitnessFunction(
            threat_weight=0.7,
            drug_likeness_weight=0.3,
            disruption_penalty=0.85,
        )
        apex = pincer.evolve(seed_smiles, fitness_fn=fitness_fn)

    Args:
        threat_weight: weight given to minimax binding survival (vs drug-likeness)
        drug_likeness_weight: weight given to QED / Lipinski score
        disruption_penalty: how harshly a single bad mutation penalises the score.
                            1.0 = pure minimax (worst threat dominates fully).
                            0.5 = average across threats (softer).
    """

    def __init__(
        self,
        threat_weight: float = 0.7,
        drug_likeness_weight: float = 0.3,
        disruption_penalty: float = 0.85,
    ):
        self.threat_weight = threat_weight
        self.dl_weight = drug_likeness_weight
        self.disruption_penalty = disruption_penalty
        self._cache: dict = {}      # smiles → (drug_vec, dl_score)

    def __call__(self, smiles: str, threats: List[ViableMutant]) -> float:
        """
        Compute fitness for `smiles` against `threats`.

        Compatible with PincerEngine's signature:
            fitness_fn(smiles, viable_threats) -> float
        """
        # --- Drug embedding (cached) ---
        if smiles not in self._cache:
            drug_vec = embed_drug(smiles)
            dl = drug_likeness_score(smiles)
            self._cache[smiles] = (drug_vec, dl)

        drug_vec, dl_score = self._cache[smiles]

        if drug_vec is None:
            return 0.0      # invalid molecule

        # --- Threat survival score ---
        if not threats:
            threat_score = 1.0
        else:
            # Compute disruption for each threat
            disruptions = [binding_disruption(drug_vec, t) for t in threats]

            # Worst-case disruption (minimax)
            worst = max(disruptions)

            # Survival = 1 - disruption, penalised by disruption_penalty
            # disruption_penalty < 1 softens the minimax for early generations
            threat_score = 1.0 - (worst * self.disruption_penalty)
            threat_score = float(np.clip(threat_score, 0.0, 1.0))

        # --- Combined score ---
        score = (
            self.threat_weight * threat_score
            + self.dl_weight * dl_score
        )
        return float(np.clip(score, 0.0, 1.0))

    def explain(self, smiles: str, threats: List[ViableMutant]) -> dict:
        """
        Return a breakdown of the fitness score for inspection / debugging.
        """
        drug_vec, dl_score = self._cache.get(smiles, (embed_drug(smiles), drug_likeness_score(smiles)))

        if drug_vec is None:
            return {"error": "invalid molecule", "score": 0.0}

        disruptions = {
            ', '.join(f"{m.wild_type}{m.position}{m.mutant}" for m in t.mutations):
            binding_disruption(drug_vec, t)
            for t in threats[:10]      # top 10 for readability
        }

        worst_threat = max(disruptions, key=disruptions.get) if disruptions else "none"
        worst_disruption = disruptions.get(worst_threat, 0.0)
        threat_score = 1.0 - (worst_disruption * self.disruption_penalty)

        return {
            "smiles": smiles,
            "total_score": self(smiles, threats),
            "threat_score": float(np.clip(threat_score, 0.0, 1.0)),
            "drug_likeness": dl_score,
            "worst_threat": worst_threat,
            "worst_disruption": worst_disruption,
            "threat_breakdown": disruptions,
            "weights": {
                "threat": self.threat_weight,
                "drug_likeness": self.dl_weight,
                "disruption_penalty": self.disruption_penalty,
            },
        }

    def plug_in_keybox(self, keybox_fn):
        """
        When KeyBox is ready, call this to replace the surrogate
        binding_disruption with real voxel physics scores.

        keybox_fn signature: (smiles: str, threat: ViableMutant) -> float
        The returned float should be a binding affinity in [0, 1]
        (higher = better binding).
        """
        import functools

        original_call = self.__class__.__call__

        @functools.wraps(original_call)
        def keybox_aware(self_inner, smiles, threats):
            if not threats:
                return drug_likeness_score(smiles)

            affinities = [keybox_fn(smiles, t) for t in threats]
            worst_affinity = min(affinities)          # minimax
            threat_score = worst_affinity

            _, dl_score = self_inner._cache.get(
                smiles,
                (None, drug_likeness_score(smiles))
            )
            return float(np.clip(
                self_inner.threat_weight * threat_score
                + self_inner.dl_weight * dl_score,
                0.0, 1.0
            ))

        self.__class__.__call__ = keybox_aware
        logger.info("ThreatAwareFitnessFunction: KeyBox backend plugged in")
