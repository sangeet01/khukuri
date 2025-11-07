"""Learning loop for continuous improvement"""

import logging
from typing import Dict, List, Optional
import numpy as np

logger = logging.getLogger('khukuri')


class LearningLoop:
    """Bayesian optimization and active learning for compound prioritization"""
    
    def __init__(self):
        self.observations = []
        self.model_performance = {}
        self.iteration = 0
    
    def add_observation(self, compound_id: str, features: Dict, outcome: Dict):
        """Add experimental observation"""
        self.observations.append({
            'compound_id': compound_id,
            'features': features,
            'outcome': outcome,
            'iteration': self.iteration
        })
        logger.info(f"Added observation for {compound_id}")
    
    def prioritize_next_experiments(self, candidates: List[Dict], n_select: int = 5) -> List[Dict]:
        """Prioritize next compounds to test using acquisition function"""
        if not self.observations:
            # Random selection for first iteration
            return self._random_selection(candidates, n_select)
        
        # Calculate acquisition scores
        scored_candidates = []
        for candidate in candidates:
            score = self._calculate_acquisition_score(candidate)
            scored_candidates.append({
                **candidate,
                'acquisition_score': score
            })
        
        # Select top candidates
        scored_candidates.sort(key=lambda x: x['acquisition_score'], reverse=True)
        selected = scored_candidates[:n_select]
        
        logger.info(f"Selected {len(selected)} candidates for next iteration")
        return selected
    
    def _calculate_acquisition_score(self, candidate: Dict) -> float:
        """Calculate acquisition score (exploration vs exploitation)"""
        # Expected improvement calculation
        predicted_value = self._predict_outcome(candidate)
        uncertainty = self._estimate_uncertainty(candidate)
        
        # Balance exploration and exploitation
        best_observed = self._get_best_observed_value()
        improvement = max(0, predicted_value - best_observed)
        
        # UCB-style acquisition
        acquisition_score = predicted_value + 2.0 * uncertainty
        
        return acquisition_score
    
    def _predict_outcome(self, candidate: Dict) -> float:
        """Predict outcome for candidate"""
        if not self.observations:
            return 0.5
        
        # Simple similarity-based prediction
        features = candidate.get('features', {})
        similarities = []
        outcomes = []
        
        for obs in self.observations:
            sim = self._calculate_similarity(features, obs['features'])
            similarities.append(sim)
            outcomes.append(obs['outcome'].get('activity', 0.0))
        
        # Weighted average by similarity
        similarities = np.array(similarities)
        outcomes = np.array(outcomes)
        
        if similarities.sum() > 0:
            weights = similarities / similarities.sum()
            prediction = np.dot(weights, outcomes)
        else:
            prediction = np.mean(outcomes)
        
        return float(prediction)
    
    def _estimate_uncertainty(self, candidate: Dict) -> float:
        """Estimate prediction uncertainty"""
        if not self.observations:
            return 1.0
        
        # Uncertainty based on distance to observed compounds
        features = candidate.get('features', {})
        min_distance = float('inf')
        
        for obs in self.observations:
            distance = 1.0 - self._calculate_similarity(features, obs['features'])
            min_distance = min(min_distance, distance)
        
        # Higher uncertainty for compounds far from observations
        uncertainty = min(min_distance, 1.0)
        return uncertainty
    
    def _calculate_similarity(self, features1: Dict, features2: Dict) -> float:
        """Calculate feature similarity"""
        # Simple Jaccard-like similarity
        common_keys = set(features1.keys()) & set(features2.keys())
        if not common_keys:
            return 0.0
        
        similarities = []
        for key in common_keys:
            val1 = features1[key]
            val2 = features2[key]
            
            if isinstance(val1, (int, float)) and isinstance(val2, (int, float)):
                # Numerical similarity
                max_val = max(abs(val1), abs(val2), 1.0)
                sim = 1.0 - abs(val1 - val2) / max_val
                similarities.append(sim)
        
        return np.mean(similarities) if similarities else 0.0
    
    def _get_best_observed_value(self) -> float:
        """Get best observed outcome"""
        if not self.observations:
            return 0.0
        
        values = [obs['outcome'].get('activity', 0.0) for obs in self.observations]
        return max(values)
    
    def _random_selection(self, candidates: List[Dict], n_select: int) -> List[Dict]:
        """Random selection for initial iteration"""
        indices = np.random.choice(len(candidates), min(n_select, len(candidates)), replace=False)
        return [candidates[i] for i in indices]
    
    def update_model(self):
        """Update predictive model with new observations"""
        if len(self.observations) < 5:
            logger.info("Insufficient data for model update")
            return
        
        # Calculate model performance
        recent_obs = self.observations[-10:]
        predictions = []
        actuals = []
        
        for obs in recent_obs:
            pred = self._predict_outcome({'features': obs['features']})
            actual = obs['outcome'].get('activity', 0.0)
            predictions.append(pred)
            actuals.append(actual)
        
        # Calculate MAE
        mae = np.mean(np.abs(np.array(predictions) - np.array(actuals)))
        
        self.model_performance[self.iteration] = {
            'mae': float(mae),
            'n_observations': len(self.observations)
        }
        
        logger.info(f"Model updated - MAE: {mae:.3f}")
    
    def next_iteration(self):
        """Move to next iteration"""
        self.iteration += 1
        self.update_model()
        logger.info(f"Starting iteration {self.iteration}")
    
    def get_learning_curve(self) -> Dict:
        """Get learning curve data"""
        return {
            'iterations': list(self.model_performance.keys()),
            'performance': [v['mae'] for v in self.model_performance.values()],
            'n_observations': [v['n_observations'] for v in self.model_performance.values()]
        }
