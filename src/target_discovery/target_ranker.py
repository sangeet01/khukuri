"""Target ranking algorithms"""

import logging
import numpy as np

logger = logging.getLogger('khukuri')


class TargetRanker:
    """Rank drug targets by importance"""
    
    def __init__(self, network_analyzer):
        self.network_analyzer = network_analyzer
    
    def rank_targets(self, disease_nodes=None):
        """Rank targets based on network centrality"""
        centrality_metrics = self.network_analyzer.calculate_centrality()
        
        if not centrality_metrics:
            logger.warning("No centrality metrics available")
            return []
        
        # Combine centrality scores
        combined_scores = {}
        for node in centrality_metrics['degree'].keys():
            score = (
                centrality_metrics['degree'].get(node, 0) * 0.3 +
                centrality_metrics['betweenness'].get(node, 0) * 0.4 +
                centrality_metrics['closeness'].get(node, 0) * 0.2 +
                centrality_metrics['eigenvector'].get(node, 0) * 0.1
            )
            combined_scores[node] = score
        
        # Sort by score
        ranked = sorted(combined_scores.items(), key=lambda x: x[1], reverse=True)
        
        return [
            {
                'protein': protein,
                'importance_score': score,
                'rank': i + 1,
                'centrality': {
                    'degree': centrality_metrics['degree'].get(protein, 0),
                    'betweenness': centrality_metrics['betweenness'].get(protein, 0)
                }
            }
            for i, (protein, score) in enumerate(ranked)
        ]
    
    def filter_druggable_targets(self, ranked_targets, min_score=0.1):
        """Filter for druggable targets"""
        return [t for t in ranked_targets if t['importance_score'] >= min_score]
