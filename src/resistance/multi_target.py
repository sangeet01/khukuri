"""Multi-target strategies"""

import logging
import networkx as nx

logger = logging.getLogger('khukuri')


class MultiTargetDesigner:
    """Design multi-target strategies"""
    
    def find_combinations(self, ranked_targets, network, max_combinations=10):
        """Find optimal multi-target combinations"""
        combinations = []
        top_targets = [t['protein'] for t in ranked_targets[:20]]
        
        for i, target1 in enumerate(top_targets):
            for target2 in top_targets[i+1:]:
                # Prefer non-interacting targets (different pathways)
                if not network.has_edge(target1, target2):
                    score = ranked_targets[i]['importance_score'] + ranked_targets[i+1]['importance_score']
                    combinations.append({
                        'targets': [target1, target2],
                        'combined_score': score,
                        'synergy_potential': 'high'
                    })
        
        combinations.sort(key=lambda x: x['combined_score'], reverse=True)
        return combinations[:max_combinations]
    
    def calculate_resistance_reduction(self, num_targets, single_target_resistance=1e-7):
        """Calculate resistance reduction for multi-target approach"""
        multi_target_resistance = single_target_resistance ** num_targets
        reduction_factor = single_target_resistance / multi_target_resistance
        
        return {
            'single_target_probability': single_target_resistance,
            'multi_target_probability': multi_target_resistance,
            'resistance_reduction_factor': reduction_factor
        }
