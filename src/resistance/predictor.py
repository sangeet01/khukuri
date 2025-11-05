"""Resistance prediction"""

import logging
from .mechanisms import load_resistance_mechanisms

logger = logging.getLogger('khukuri')


class ResistancePredictor:
    """Predict resistance likelihood"""
    
    def __init__(self):
        self.mechanisms = load_resistance_mechanisms()
    
    def predict_likelihood(self, mol, pathogen_type, target_protein):
        """Predict resistance development likelihood"""
        resistance_score = 0.5  # Base probability
        
        # Check if target has known resistance mutations
        if pathogen_type in self.mechanisms:
            if target_protein in self.mechanisms[pathogen_type]:
                resistance_score += 0.2
        
        # Molecular complexity reduces resistance
        from rdkit.Chem import Descriptors
        complexity = Descriptors.BertzCT(mol)
        if complexity > 500:
            resistance_score -= 0.1
        
        resistance_score = max(0.1, min(0.9, resistance_score))
        
        return {
            'resistance_probability': resistance_score,
            'risk_level': 'high' if resistance_score > 0.7 else 'medium' if resistance_score > 0.4 else 'low',
            'known_mechanisms': self.mechanisms.get(pathogen_type, {}).get(target_protein, [])
        }
