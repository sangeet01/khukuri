"""
KeyBox → Khukuri bridge.

Connects KeyBox's NibbleEngine (11-channel voxel docking) to PINCER's
ThreatAwareFitnessFunction, replacing the physicochemical surrogate with
real binding physics.

Also provides KhukuriKeyBox — a unified interface for the AMRDiscoveryWorkflow
to call docking without knowing KeyBox internals.

Usage:
    from src.integrations.keybox_bridge import KhukuriKeyBox, plug_keybox_into_pincer

    # 1. Initialise
    kb = KhukuriKeyBox(pdb_path="PBP2a.pdb", binding_site_center=(x, y, z))

    # 2. Plug into PINCER fitness function
    fitness_fn = ThreatAwareFitnessFunction()
    plug_keybox_into_pincer(fitness_fn, kb)

    # 3. Run PINCER — now uses real voxel physics
    apex = pincer.evolve(seeds, fitness_fn=fitness_fn)

    # 4. Score any SMILES directly
    score = kb.score(smiles)
"""

import logging
import sys
import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger('khukuri')


# ---------------------------------------------------------------------------
# KeyBox import (graceful — KeyBox may be a sibling repo or installed pkg)
# ---------------------------------------------------------------------------

def _import_keybox():
    """Try multiple import strategies for KeyBox."""
    # 1. Already installed as package
    try:
        from box.key.nibble_bridge import NibbleEngine
        from box.key.engine import KeyBoxSystem
        from box.key.models import Molecule, MoleculeType
        return NibbleEngine, KeyBoxSystem, Molecule, MoleculeType
    except ImportError:
        pass

    # 2. Sibling repo on the filesystem
    for candidate in [
        os.path.join(os.path.dirname(__file__), '..', '..', '..', 'keybox'),
        os.path.join(os.path.expanduser('~'), 'keybox'),
        os.path.join(os.path.expanduser('~'), 'claude', 'keybox'),
    ]:
        box_path = os.path.join(candidate, 'box')
        if os.path.isdir(box_path):
            if candidate not in sys.path:
                sys.path.insert(0, candidate)
            try:
                from box.key.nibble_bridge import NibbleEngine
                from box.key.engine import KeyBoxSystem
                from box.key.models import Molecule, MoleculeType
                logger.info(f"KeyBox loaded from {candidate}")
                return NibbleEngine, KeyBoxSystem, Molecule, MoleculeType
            except ImportError:
                pass

    return None, None, None, None


NibbleEngine, KeyBoxSystem, Molecule, MoleculeType = _import_keybox()
KEYBOX_AVAILABLE = NibbleEngine is not None

if KEYBOX_AVAILABLE:
    logger.info(f"KeyBox bridge: NibbleEngine available")
else:
    logger.warning("KeyBox bridge: KeyBox not found — install or place sibling to khukuri/")


# ---------------------------------------------------------------------------
# Molecule builder (SMILES → KeyBox Molecule)
# ---------------------------------------------------------------------------

def _smiles_to_keybox_molecule(smiles: str, molecule_type=None):
    """
    Convert a SMILES string to a KeyBox Molecule object with
    3D coordinates, partial charges and hydrophobicity.
    """
    if not KEYBOX_AVAILABLE:
        return None

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
        from rdkit.Chem import rdPartialCharges

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
            # Fallback: distance geometry
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.UFFOptimizeMolecule(mol)
        mol = Chem.RemoveHs(mol)

        conf = mol.GetConformer()
        coords = [
            list(conf.GetAtomPosition(i))
            for i in range(mol.GetNumAtoms())
        ]

        # Gasteiger partial charges
        rdPartialCharges.ComputeGasteigerCharges(mol)
        charges = [
            float(mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge'))
            for i in range(mol.GetNumAtoms())
        ]
        # Replace NaN charges with 0
        charges = [0.0 if (c != c) else c for c in charges]

        # Per-atom hydrophobicity proxy (carbon=1, heteroatom=0)
        hydro = [
            1.0 if mol.GetAtomWithIdx(i).GetAtomicNum() == 6 else 0.0
            for i in range(mol.GetNumAtoms())
        ]

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)

        mol_type = molecule_type or MoleculeType.DRUG

        return Molecule(
            smiles=smiles,
            molecule_type=mol_type,
            coords=coords,
            charges=charges,
            hydrophobicity=hydro,
            molecular_weight=mw,
            logp=logp,
        )

    except Exception as exc:
        logger.warning(f"SMILES→Molecule failed for {smiles}: {exc}")
        return None


# ---------------------------------------------------------------------------
# KhukuriKeyBox — unified docking interface for Khukuri
# ---------------------------------------------------------------------------

class KhukuriKeyBox:
    """
    Unified KeyBox interface for Khukuri.

    Wraps NibbleEngine to provide:
        score(smiles) -> float          binding affinity in [0, 1]
        score_batch(smiles_list)        parallel scoring
        score_vs_threat(smiles, threat) for ThreatAwareFitnessFunction

    Args:
        pdb_path: path to receptor PDB file (e.g. PBP2a_MRSA.pdb)
        binding_site_center: (x, y, z) center of binding pocket in Angstroms
        pocket_range: radius around center to include (default 10.0 Å)
        blur_radius: Gaussian blur for atom projection (default 1.5 Å)
        dim: voxel grid dimension (default 30)
        resolution: Angstroms per voxel (default 0.8)
        score_range: (min, max) raw affinity range for normalisation
    """

    def __init__(
        self,
        pdb_path: Optional[str] = None,
        binding_site_center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
        pocket_range: float = 10.0,
        blur_radius: float = 1.5,
        dim: int = 30,
        resolution: float = 0.8,
        score_range: Tuple[float, float] = (-500.0, 500.0),
    ):
        self.pdb_path = pdb_path
        self.center = binding_site_center
        self.pocket_range = pocket_range
        self.blur_radius = blur_radius
        self.score_range = score_range
        self._cache: Dict[str, float] = {}

        if not KEYBOX_AVAILABLE:
            logger.warning("KhukuriKeyBox: KeyBox unavailable, scores will be 0.5")
            self._engine = None
            self._pocket_loaded = False
            return

        self._engine = NibbleEngine(dim_x=dim, dim_y=dim, dim_z=dim,
                                    resolution=resolution)
        logger.info(f"KhukuriKeyBox: NibbleEngine mode={self._engine.mode}")

        self._pocket_loaded = False
        if pdb_path and os.path.exists(pdb_path):
            self._load_pocket()
        elif pdb_path:
            logger.warning(f"KhukuriKeyBox: PDB not found at {pdb_path}")

    def _load_pocket(self):
        """Load and voxelise the binding pocket from PDB."""
        try:
            # Auto-detect binding site if no center provided
            if self.center == (0.0, 0.0, 0.0):
                n_atoms = self._engine.load_full_protein(
                    self.pdb_path, self.blur_radius
                )
                if n_atoms > 0:
                    sites = self._engine.find_binding_sites(n_sites=1)
                    if sites:
                        self.center = sites[0]['center']
                        logger.info(
                            f"KhukuriKeyBox: active site auto-detected at "
                            f"{self.center}  score={sites[0]['score']:.4f}  "
                            f"(from {n_atoms} atoms)"
                        )
                        self._pocket_loaded = True
                        self._detected_sites = sites
                        return
            # Manual center provided — load just the pocket region
            n_atoms = self._engine.load_pdb_pocket(
                self.pdb_path, self.center,
                self.pocket_range, self.blur_radius
            )
            self._pocket_loaded = n_atoms > 0
            logger.info(
                f"KhukuriKeyBox: pocket loaded "
                f"({n_atoms} atoms, mode={self._engine.mode})"
            )
        except Exception as exc:
            logger.error(f"KhukuriKeyBox: pocket load failed: {exc}")
            self._pocket_loaded = False

    def detect_binding_sites(self, n_sites=3) -> list:
        """
        Run full inverse-vector active site detection on the loaded protein.
        Returns top n_sites candidates sorted by druggability score.
        """
        if not KEYBOX_AVAILABLE or self._engine is None:
            return []
        if not self._pocket_loaded:
            logger.warning("KhukuriKeyBox: load a protein first")
            return []
        sites = self._engine.find_binding_sites(n_sites=n_sites)
        for s in sites:
            logger.info(
                f"  Site {s['rank']}: center={s['center']}  "
                f"score={s['score']:.4f}  "
                f"hba={s['channel_profile']['hba']:.3f}  "
                f"hbd={s['channel_profile']['hbd']:.3f}  "
                f"lipo={s['channel_profile']['lipo']:.3f}"
            )
        return sites

    def load_pocket(self, pdb_path: str,
                    center: Tuple[float, float, float],
                    pocket_range: float = 10.0):
        """Load a new pocket (call this to swap targets)."""
        self.pdb_path = pdb_path
        self.center = center
        self.pocket_range = pocket_range
        self._cache.clear()
        self._load_pocket()

    def score(self, smiles: str) -> float:
        """
        Score a SMILES string against the loaded pocket.
        Returns affinity in [0, 1] — higher = better binding.
        """
        if smiles in self._cache:
            return self._cache[smiles]

        if not KEYBOX_AVAILABLE or self._engine is None:
            return 0.5

        if not self._pocket_loaded:
            logger.warning("KhukuriKeyBox: no pocket loaded, returning 0.5")
            return 0.5

        try:
            # Build a minimal KeyBoxSystem to hold the molecule
            drug_system = _MinimalDrugSystem(smiles)
            if drug_system.molecule is None:
                return 0.0

            raw = self._engine.project_molecule(drug_system, self.center)
            normalised = self._normalise(raw)
            self._cache[smiles] = normalised
            return normalised

        except Exception as exc:
            logger.warning(f"KhukuriKeyBox.score failed for {smiles}: {exc}")
            return 0.0

    def score_batch(self, smiles_list: List[str]) -> List[float]:
        """Score a list of SMILES strings."""
        return [self.score(s) for s in smiles_list]

    def score_vs_threat(self, smiles: str, threat) -> float:
        """
        Score a drug against a specific resistance threat.

        For mutational threats: score against the mutant pocket if available,
        otherwise use base pocket score penalised by threat fitness cost.
        For HGT threats: penalise by transfer probability × fitness.
        """
        base_score = self.score(smiles)

        # Check if threat has mutation info (ViableMutant-like)
        mutations = getattr(threat, 'mutations', [])
        fitness_cost = 0.0
        for m in mutations:
            fitness_cost += getattr(m, 'fitness_cost', 0.0)
        avg_cost = fitness_cost / max(len(mutations), 1)

        # HGT threat: additional penalty from transfer probability
        transfer_prob = getattr(threat, 'transfer_probability', 0.0)
        hgt_penalty = transfer_prob * 0.3   # scaled — high prob = bigger hit

        # Combined: base affinity penalised by resistance mechanisms
        combined = base_score * (1.0 - avg_cost * 0.5) * (1.0 - hgt_penalty)
        return float(np.clip(combined, 0.0, 1.0))

    def get_top_binders(
        self,
        smiles_list: List[str],
        top_n: int = 10,
    ) -> List[Dict]:
        """Score all candidates and return top_n sorted by affinity."""
        scored = [
            {"smiles": s, "affinity": self.score(s)}
            for s in smiles_list
        ]
        scored.sort(key=lambda x: x["affinity"], reverse=True)
        return scored[:top_n]

    def _normalise(self, raw: float) -> float:
        """Normalise raw affinity score to [0, 1]."""
        lo, hi = self.score_range
        if hi == lo:
            return 0.5
        return float(np.clip((raw - lo) / (hi - lo), 0.0, 1.0))

    @property
    def is_ready(self) -> bool:
        return KEYBOX_AVAILABLE and self._pocket_loaded


# ---------------------------------------------------------------------------
# Minimal drug system shim (KeyBoxSystem expects .apis + .excipients)
# ---------------------------------------------------------------------------

class _MinimalDrugSystem:
    """Minimal stub so NibbleEngine.project_molecule() gets what it needs."""

    def __init__(self, smiles: str):
        self.molecule = _smiles_to_keybox_molecule(smiles)
        self.apis = [self.molecule] if self.molecule else []
        self.excipients = []


# ---------------------------------------------------------------------------
# Plug KeyBox into ThreatAwareFitnessFunction
# ---------------------------------------------------------------------------

def plug_keybox_into_pincer(fitness_fn, keybox: KhukuriKeyBox) -> None:
    """
    Replace ThreatAwareFitnessFunction's surrogate with real KeyBox physics.

    After calling this, fitness_fn(smiles, threats) uses voxel docking
    scores instead of physicochemical similarity.

    Args:
        fitness_fn: a ThreatAwareFitnessFunction instance
        keybox: an initialised KhukuriKeyBox with pocket loaded
    """
    if not keybox.is_ready:
        logger.warning(
            "plug_keybox_into_pincer: KeyBox not ready "
            "(no pocket loaded or KeyBox unavailable). "
            "Fitness function unchanged."
        )
        return

    # Use the existing plug_in_keybox hook on ThreatAwareFitnessFunction
    def keybox_fn(smiles: str, threat) -> float:
        return keybox.score_vs_threat(smiles, threat)

    fitness_fn.plug_in_keybox(keybox_fn)
    logger.info("plug_keybox_into_pincer: KeyBox physics active in PINCER")


# ---------------------------------------------------------------------------
# Wire into AMRDiscoveryWorkflow
# ---------------------------------------------------------------------------

def create_mrsa_keybox(
    pdb_path: Optional[str] = None,
    binding_site_center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> KhukuriKeyBox:
    """
    Factory for MRSA-specific KeyBox instance.

    PBP2a binding site center defaults are approximate —
    provide your own from the PDB structure (e.g. 2OLV).
    """
    # Known approximate PBP2a active site from PDB 2OLV
    default_center = (23.5, 18.2, 31.7)

    return KhukuriKeyBox(
        pdb_path=pdb_path,
        binding_site_center=binding_site_center or default_center,
        pocket_range=12.0,
        blur_radius=1.5,
        dim=30,
        resolution=0.8,
        score_range=(-300.0, 300.0),
    )
