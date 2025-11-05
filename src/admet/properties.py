"""Physicochemical property calculations"""

import logging
from rdkit.Chem import Descriptors
import numpy as np

logger = logging.getLogger('khukuri')


def calculate_properties(mol):
    """Calculate all physicochemical properties"""
    try:
        return {
            'mw': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'hbd': Descriptors.NumHDonors(mol),
            'hba': Descriptors.NumHAcceptors(mol),
            'tpsa': Descriptors.TPSA(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'aromatic_rings': Descriptors.NumAromaticRings(mol),
            'heavy_atoms': mol.GetNumHeavyAtoms(),
            'formal_charge': Descriptors.FormalCharge(mol),
            'molar_refractivity': Descriptors.MolMR(mol)
        }
    except Exception as e:
        logger.error(f"Property calculation failed: {e}")
        return {}


def predict_solubility(mol):
    """Predict aqueous solubility (LogS)"""
    try:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        rotbonds = Descriptors.NumRotatableBonds(mol)
        logs = 0.16 - 0.63*np.log10(mw) - 0.0062*mw + 0.066*rotbonds - 0.74*logp
        return logs
    except:
        return None


def calculate_bbb_permeability(mol):
    """Calculate Blood-Brain Barrier permeability"""
    try:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        
        bbb_score = 1.0
        if mw > 450: bbb_score -= 0.3
        if logp < 1 or logp > 3: bbb_score -= 0.2
        if tpsa > 90: bbb_score -= 0.3
        if hbd > 5: bbb_score -= 0.2
        
        return max(0, bbb_score)
    except:
        return None
