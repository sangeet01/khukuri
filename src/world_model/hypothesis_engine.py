"""Hypothesis generation and testing engine"""

import logging
from typing import Dict, List, Optional
import numpy as np

logger = logging.getLogger('khukuri')


class HypothesisEngine:
    """Generate and test scientific hypotheses"""
    
    def __init__(self, state_tracker, knowledge_graph):
        self.state_tracker = state_tracker
        self.knowledge_graph = knowledge_graph
        self.hypothesis_templates = self._initialize_templates()
    
    def _initialize_templates(self) -> List[Dict]:
        """Hypothesis templates"""
        return [
            {
                'template': 'Compound {compound} binds target {target} with high affinity',
                'evidence_required': ['docking_score', 'binding_assay'],
                'confidence_threshold': 0.7
            },
            {
                'template': 'Compound {compound} shows activity against {pathogen}',
                'evidence_required': ['mic_data', 'growth_inhibition'],
                'confidence_threshold': 0.8
            },
            {
                'template': 'Target {target} is essential for {pathogen} survival',
                'evidence_required': ['knockout_study', 'literature'],
                'confidence_threshold': 0.75
            },
            {
                'template': 'Mutation {mutation} in {gene} confers resistance to {drug}',
                'evidence_required': ['genomic_data', 'phenotype_data'],
                'confidence_threshold': 0.85
            }
        ]
    
    def generate_hypotheses(self, context: Optional[Dict] = None) -> List[Dict]:
        """Generate hypotheses from current state"""
        hypotheses = []
        
        # Get current state
        compounds = self.state_tracker.query_state('compound')
        targets = self.state_tracker.query_state('target')
        
        # Generate compound-target hypotheses
        for compound in compounds[:5]:
            for target in targets[:3]:
                hyp = self._generate_binding_hypothesis(compound, target)
                if hyp:
                    hypotheses.append(hyp)
        
        # Generate resistance hypotheses
        strains = self.state_tracker.query_state('strain')
        for strain in strains[:3]:
            hyp = self._generate_resistance_hypothesis(strain)
            if hyp:
                hypotheses.append(hyp)
        
        logger.info(f"Generated {len(hypotheses)} hypotheses")
        return hypotheses
    
    def _generate_binding_hypothesis(self, compound: Dict, target: Dict) -> Optional[Dict]:
        """Generate compound-target binding hypothesis"""
        compound_id = compound.get('compound_id')
        target_id = target.get('target_id')
        
        # Check if already tested
        existing = [h for h in self.state_tracker.hypotheses 
                   if compound_id in h.get('hypothesis', '') and target_id in h.get('hypothesis', '')]
        
        if existing:
            return None
        
        hypothesis = f"Compound {compound_id} binds target {target_id} with high affinity"
        
        # Estimate prior confidence
        druggability = target.get('druggability', 0.5)
        drug_likeness = compound.get('drug_likeness', 0.5)
        confidence = (druggability + drug_likeness) / 2
        
        return {
            'hypothesis': hypothesis,
            'type': 'binding',
            'entities': {'compound': compound_id, 'target': target_id},
            'confidence': confidence,
            'evidence_required': ['docking_score', 'binding_assay'],
            'status': 'proposed'
        }
    
    def _generate_resistance_hypothesis(self, strain: Dict) -> Optional[Dict]:
        """Generate resistance mechanism hypothesis"""
        strain_id = strain.get('strain_id')
        
        if 'resistance_profile' not in strain:
            return None
        
        hypothesis = f"Strain {strain_id} exhibits resistance via efflux pump upregulation"
        
        return {
            'hypothesis': hypothesis,
            'type': 'resistance_mechanism',
            'entities': {'strain': strain_id},
            'confidence': 0.6,
            'evidence_required': ['expression_data', 'efflux_assay'],
            'status': 'proposed'
        }
    
    def test_hypothesis(self, hypothesis_idx: int, evidence: Dict) -> Dict:
        """Test hypothesis with evidence"""
        if hypothesis_idx >= len(self.state_tracker.hypotheses):
            return {'error': 'Hypothesis not found'}
        
        hypothesis = self.state_tracker.hypotheses[hypothesis_idx]
        
        # Check evidence
        required = hypothesis.get('evidence_required', [])
        provided = list(evidence.keys())
        
        evidence_score = len(set(required) & set(provided)) / len(required) if required else 0
        
        # Evaluate evidence quality
        quality_scores = []
        for ev_type, ev_data in evidence.items():
            if ev_type == 'docking_score':
                # Good docking score < -7 kcal/mol
                score = 1.0 if ev_data < -7 else 0.5
            elif ev_type == 'mic_data':
                # Low MIC is good
                score = 1.0 if ev_data < 1.0 else 0.5
            else:
                score = 0.7  # Default
            
            quality_scores.append(score)
        
        avg_quality = np.mean(quality_scores) if quality_scores else 0.5
        
        # Final confidence
        final_confidence = (evidence_score + avg_quality) / 2
        
        # Determine status
        if final_confidence > 0.7:
            status = 'validated'
        elif final_confidence < 0.3:
            status = 'rejected'
        else:
            status = 'inconclusive'
        
        # Update hypothesis
        self.state_tracker.update_hypothesis_status(hypothesis_idx, status, {
            'evidence': evidence,
            'confidence': final_confidence,
            'evidence_score': evidence_score,
            'quality_score': avg_quality
        })
        
        logger.info(f"Hypothesis {hypothesis_idx} tested: {status} (confidence: {final_confidence:.2f})")
        
        return {
            'hypothesis_idx': hypothesis_idx,
            'status': status,
            'confidence': final_confidence,
            'evidence_completeness': evidence_score,
            'evidence_quality': avg_quality
        }
    
    def prioritize_hypotheses(self) -> List[Dict]:
        """Prioritize hypotheses for testing"""
        active = self.state_tracker.get_active_hypotheses()
        
        # Score by confidence and feasibility
        scored = []
        for i, hyp in enumerate(active):
            score = hyp.get('confidence', 0.5)
            
            # Boost if evidence is easy to obtain
            required = hyp.get('evidence_required', [])
            if 'docking_score' in required:
                score += 0.1  # Easy to compute
            
            scored.append({
                'index': i,
                'hypothesis': hyp,
                'priority_score': score
            })
        
        # Sort by priority
        scored.sort(key=lambda x: x['priority_score'], reverse=True)
        
        return scored
    
    def suggest_experiments(self, hypothesis_idx: int) -> List[Dict]:
        """Suggest experiments to test hypothesis"""
        if hypothesis_idx >= len(self.state_tracker.hypotheses):
            return []
        
        hypothesis = self.state_tracker.hypotheses[hypothesis_idx]
        required = hypothesis.get('evidence_required', [])
        
        experiments = []
        
        for ev_type in required:
            if ev_type == 'docking_score':
                experiments.append({
                    'type': 'molecular_docking',
                    'parameters': {
                        'compound': hypothesis['entities'].get('compound'),
                        'target': hypothesis['entities'].get('target')
                    },
                    'cost': 'low',
                    'time': 'hours'
                })
            elif ev_type == 'mic_data':
                experiments.append({
                    'type': 'mic_assay',
                    'parameters': {
                        'compound': hypothesis['entities'].get('compound'),
                        'strain': hypothesis['entities'].get('strain')
                    },
                    'cost': 'medium',
                    'time': 'days'
                })
            elif ev_type == 'binding_assay':
                experiments.append({
                    'type': 'spr_binding',
                    'parameters': {
                        'compound': hypothesis['entities'].get('compound'),
                        'target': hypothesis['entities'].get('target')
                    },
                    'cost': 'high',
                    'time': 'weeks'
                })
        
        return experiments
