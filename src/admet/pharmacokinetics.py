"""Pharmacokinetic predictions"""

import logging
from rdkit.Chem import Descriptors

logger = logging.getLogger('khukuri')


def predict_bioavailability(mol):
    """Predict oral bioavailability"""
    try:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        rotbonds = Descriptors.NumRotatableBonds(mol)
        
        ba_score = 1.0
        if mw > 500: ba_score -= 0.2
        if logp < 0 or logp > 5: ba_score -= 0.2
        if tpsa > 140: ba_score -= 0.3
        if hbd > 5: ba_score -= 0.2
        if rotbonds > 10: ba_score -= 0.1
        
        return max(0, ba_score)
    except:
        return 0.5


def predict_half_life(mol):
    """Predict plasma half-life (hours)"""
    try:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        half_life = 2 + (logp * 3) + (mw / 200)
        return max(0.5, min(24, half_life))
    except:
        return 4.0


def predict_clearance(mol):
    """Predict plasma clearance (mL/min/kg)"""
    try:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        clearance = 20 - (logp * 2) + (mw / 100)
        return max(1, min(50, clearance))
    except:
        return 15.0


def predict_protein_binding(mol):
    """Predict protein binding percentage"""
    try:
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        binding = 50 + (logp * 15) - (tpsa * 0.2)
        return max(10, min(99, binding))
    except:
        return 70.0


def predict_volume_distribution(mol):
    """Predict volume of distribution (L/kg)"""
    try:
        logp = Descriptors.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        vd = 1 + (logp * 2) + (mw / 500)
        return max(0.1, min(20, vd))
    except:
        return 1.0


def predict_dose(mol, binding_affinity):
    """Predict therapeutic dose (mg)"""
    if binding_affinity is None:
        return None
    try:
        mw = Descriptors.MolWt(mol)
        ba = predict_bioavailability(mol)
        potency_factor = abs(binding_affinity) / 10
        dose = (100 / potency_factor) * (mw / 300) / max(0.1, ba)
        return max(1, min(1000, dose))
    except:
        return None


def predict_admet(mol):
    """Comprehensive ADMET prediction"""
    try:
        return {
            'bioavailability': predict_bioavailability(mol),
            'half_life': predict_half_life(mol),
            'clearance': predict_clearance(mol),
            'protein_binding': predict_protein_binding(mol),
            'volume_distribution': predict_volume_distribution(mol)
        }
    except Exception as e:
        logger.error(f"ADMET prediction failed: {e}")
        return {}
