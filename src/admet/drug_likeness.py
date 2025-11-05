"""Drug-likeness calculations"""

import logging
from rdkit.Chem import Descriptors, QED

logger = logging.getLogger('khukuri')


def lipinski_ro5(mol):
    """Calculate Lipinski Rule of 5 violations"""
    try:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if hbd > 5: violations += 1
        if hba > 10: violations += 1
        
        return violations
    except Exception as e:
        logger.error(f"Lipinski calculation failed: {e}")
        return 5


def calculate_qed(mol):
    """Calculate Quantitative Estimate of Drug-likeness"""
    try:
        return QED.qed(mol)
    except:
        # Fallback calculation
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        mw_score = 1 - abs(mw - 350) / 150 if abs(mw - 350) < 150 else 0
        logp_score = 1 - abs(logp - 2.5) / 2.5 if abs(logp - 2.5) < 2.5 else 0
        hbd_score = 1 - hbd / 5 if hbd <= 5 else 0
        hba_score = 1 - hba / 10 if hba <= 10 else 0
        
        return (mw_score + logp_score + hbd_score + hba_score) / 4


def calculate_drug_likeness(mol):
    """Comprehensive drug-likeness score"""
    try:
        lipinski_violations = lipinski_ro5(mol)
        qed_score = calculate_qed(mol)
        
        mw = Descriptors.MolWt(mol)
        mw_score = max(0, 1 - abs(mw - 350) / 200)
        
        logp = Descriptors.MolLogP(mol)
        logp_score = max(0, 1 - abs(logp - 2.5) / 3)
        
        drug_likeness = (
            (1 - lipinski_violations / 5) * 0.4 +
            qed_score * 0.3 +
            mw_score * 0.15 +
            logp_score * 0.15
        )
        
        return max(0, min(1, drug_likeness))
    except Exception as e:
        logger.error(f"Drug-likeness calculation failed: {e}")
        return 0.0
