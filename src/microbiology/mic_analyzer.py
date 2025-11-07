"""MIC (Minimum Inhibitory Concentration) analysis"""

import logging
from typing import Dict, List, Optional
import numpy as np

logger = logging.getLogger('khukuri')


class MICAnalyzer:
    """Analyze MIC data and predict resistance"""
    
    def __init__(self):
        self.breakpoints = self._initialize_breakpoints()
        self.mic_data = {}
    
    def _initialize_breakpoints(self) -> Dict:
        """CLSI/EUCAST breakpoints (μg/mL)"""
        return {
            'rifampicin': {'susceptible': 1.0, 'resistant': 4.0},
            'isoniazid': {'susceptible': 0.2, 'resistant': 1.0},
            'ciprofloxacin': {'susceptible': 1.0, 'resistant': 4.0},
            'vancomycin': {'susceptible': 4.0, 'resistant': 16.0},
            'meropenem': {'susceptible': 2.0, 'resistant': 8.0},
            'colistin': {'susceptible': 2.0, 'resistant': 4.0},
            'linezolid': {'susceptible': 4.0, 'resistant': 8.0}
        }
    
    def add_mic_result(self, compound_id: str, strain_id: str, drug: str, mic_value: float):
        """Add MIC test result"""
        key = f"{compound_id}_{strain_id}_{drug}"
        self.mic_data[key] = {
            'compound_id': compound_id,
            'strain_id': strain_id,
            'drug': drug,
            'mic': mic_value,
            'interpretation': self._interpret_mic(drug, mic_value)
        }
        logger.info(f"Added MIC result: {key} = {mic_value} ug/mL")
    
    def _interpret_mic(self, drug: str, mic_value: float) -> str:
        """Interpret MIC value using breakpoints"""
        if drug not in self.breakpoints:
            return 'unknown'
        
        bp = self.breakpoints[drug]
        if mic_value <= bp['susceptible']:
            return 'susceptible'
        elif mic_value >= bp['resistant']:
            return 'resistant'
        else:
            return 'intermediate'
    
    def get_mic_profile(self, strain_id: str) -> Dict:
        """Get complete MIC profile for strain"""
        profile = {
            'strain_id': strain_id,
            'results': [],
            'resistance_count': 0,
            'susceptible_count': 0
        }
        
        for key, data in self.mic_data.items():
            if data['strain_id'] == strain_id:
                profile['results'].append(data)
                if data['interpretation'] == 'resistant':
                    profile['resistance_count'] += 1
                elif data['interpretation'] == 'susceptible':
                    profile['susceptible_count'] += 1
        
        return profile
    
    def compare_compounds(self, compound_ids: List[str], strain_id: str) -> List[Dict]:
        """Compare MIC values across compounds"""
        results = []
        
        for compound_id in compound_ids:
            compound_mics = [data for key, data in self.mic_data.items()
                           if data['compound_id'] == compound_id and data['strain_id'] == strain_id]
            
            if compound_mics:
                avg_mic = np.mean([d['mic'] for d in compound_mics])
                results.append({
                    'compound_id': compound_id,
                    'avg_mic': avg_mic,
                    'n_tests': len(compound_mics),
                    'interpretations': [d['interpretation'] for d in compound_mics]
                })
        
        return sorted(results, key=lambda x: x['avg_mic'])
    
    def predict_cross_resistance(self, strain_id: str, known_drug: str, target_drug: str) -> Dict:
        """Predict cross-resistance based on mechanism"""
        cross_resistance_map = {
            'rifampicin': ['rifabutin', 'rifapentine'],
            'isoniazid': ['ethionamide'],
            'ciprofloxacin': ['levofloxacin', 'moxifloxacin'],
            'meropenem': ['imipenem', 'ertapenem']
        }
        
        prediction = {
            'known_drug': known_drug,
            'target_drug': target_drug,
            'cross_resistance_likely': False,
            'confidence': 0.0
        }
        
        # Check if strain is resistant to known drug
        profile = self.get_mic_profile(strain_id)
        known_resistant = any(r['drug'] == known_drug and r['interpretation'] == 'resistant' 
                             for r in profile['results'])
        
        if known_resistant and target_drug in cross_resistance_map.get(known_drug, []):
            prediction['cross_resistance_likely'] = True
            prediction['confidence'] = 0.7
        
        return prediction
    
    def calculate_resistance_index(self, strain_id: str) -> float:
        """Calculate overall resistance index (0-1)"""
        profile = self.get_mic_profile(strain_id)
        total = len(profile['results'])
        
        if total == 0:
            return 0.0
        
        return profile['resistance_count'] / total
