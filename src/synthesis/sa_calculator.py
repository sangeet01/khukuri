"""Synthetic accessibility calculator"""

import logging
from rdkit.Chem import rdMolDescriptors

logger = logging.getLogger('khukuri')


def calculate_sa_score(mol):
    """Calculate synthetic accessibility score (1-10, lower is better)"""
    try:
        sa_score = rdMolDescriptors.BertzCT(mol) / 100
        rings = mol.GetRingInfo().NumRings()
        heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])
        
        complexity_penalty = rings * 0.2 + heteroatoms * 0.1
        final_sa = min(10, sa_score + complexity_penalty)
        
        return max(1, final_sa)
    except:
        rings = mol.GetRingInfo().NumRings()
        heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])
        return min(10, 1 + rings * 0.5 + heteroatoms * 0.3)
