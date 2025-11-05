"""Toxicity prediction"""

import logging
from rdkit import Chem

logger = logging.getLogger('khukuri')


def check_structural_alerts(mol):
    """Check for structural toxicity alerts"""
    alerts = []
    toxic_patterns = {
        'nitro_group': '[N+](=O)[O-]',
        'alpha_beta_unsaturated_carbonyl': 'C=C-C=O',
        'halogenated_alkene': '[Cl,Br,I]C=C',
        'hydrazine': 'N-N',
        'sulfonyl_phosphonyl': '[S,P](=O)(=O)',
        'epoxide': 'C1OC1',
        'aziridine': 'C1NC1',
        'quinone': 'C1=CC(=O)C=CC1=O'
    }
    
    for alert_name, pattern in toxic_patterns.items():
        try:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                alerts.append(alert_name)
        except:
            continue
    
    return alerts


def predict_toxicity(mol):
    """Predict toxicity risk score (0-1, higher is riskier)"""
    try:
        alerts = check_structural_alerts(mol)
        alert_score = min(1.0, len(alerts) / 5.0)
        
        # Additional heuristics
        from rdkit.Chem import Descriptors
        logp = Descriptors.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        
        # Very lipophilic compounds tend to be more toxic
        lipophilicity_penalty = 0
        if logp > 5:
            lipophilicity_penalty = 0.2
        
        # Very large compounds may have toxicity issues
        size_penalty = 0
        if mw > 600:
            size_penalty = 0.1
        
        total_toxicity = min(1.0, alert_score + lipophilicity_penalty + size_penalty)
        
        return {
            'toxicity_score': total_toxicity,
            'structural_alerts': alerts,
            'risk_level': 'high' if total_toxicity > 0.6 else 'medium' if total_toxicity > 0.3 else 'low'
        }
    except Exception as e:
        logger.error(f"Toxicity prediction failed: {e}")
        return {'toxicity_score': 0.5, 'structural_alerts': [], 'risk_level': 'unknown'}


def predict_cyp_inhibition(mol):
    """Predict CYP450 inhibition risk"""
    try:
        from rdkit.Chem import Descriptors
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        
        risk_score = 0
        if mw > 400: risk_score += 1
        if logp > 3: risk_score += 1
        if tpsa < 60: risk_score += 1
        
        return min(5, risk_score)
    except:
        return 2
