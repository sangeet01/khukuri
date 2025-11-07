"""Kosmos-style world model engine with memory and code execution"""

import logging
from typing import Dict, List, Any, Optional
from datetime import datetime
import json

logger = logging.getLogger('khukuri')


class KosmosEngine:
    """Kosmos-style engine: memory, code execution, structured reasoning"""
    
    def __init__(self, state_tracker, knowledge_graph):
        self.state_tracker = state_tracker
        self.knowledge_graph = knowledge_graph
        self.memory = []
        self.code_history = []
        self.observations = []
    
    def observe(self, observation_type: str, data: Dict):
        """Record observation (assay result, literature finding, etc.)"""
        obs = {
            'timestamp': datetime.now().isoformat(),
            'type': observation_type,
            'data': data
        }
        self.observations.append(obs)
        self.memory.append(f"Observed {observation_type}: {data}")
        logger.info(f"Recorded observation: {observation_type}")
    
    def reason(self, question: str, context: Optional[Dict] = None) -> Dict:
        """Structured reasoning over world state"""
        # Query relevant state
        relevant_compounds = self.state_tracker.query_state('compound', context)
        relevant_targets = self.state_tracker.query_state('target', context)
        
        # Find patterns in knowledge graph
        patterns = self._find_patterns()
        
        # Generate reasoning
        reasoning = {
            'question': question,
            'relevant_entities': {
                'compounds': len(relevant_compounds),
                'targets': len(relevant_targets)
            },
            'patterns': patterns,
            'conclusion': self._synthesize_conclusion(question, patterns)
        }
        
        self.memory.append(f"Reasoned about: {question}")
        return reasoning
    
    def _find_patterns(self) -> List[Dict]:
        """Find patterns in knowledge graph"""
        patterns = []
        
        # Multi-target compounds
        multi_target = self.knowledge_graph.find_multi_target_compounds(min_targets=2)
        if multi_target:
            patterns.append({
                'type': 'multi_target',
                'count': len(multi_target),
                'significance': 'high'
            })
        
        # Resistance correlations
        if len(self.observations) > 5:
            patterns.append({
                'type': 'resistance_trend',
                'count': len(self.observations),
                'significance': 'medium'
            })
        
        return patterns
    
    def _synthesize_conclusion(self, question: str, patterns: List[Dict]) -> str:
        """Synthesize conclusion from patterns"""
        if not patterns:
            return "Insufficient data for conclusion"
        
        if any(p['type'] == 'multi_target' for p in patterns):
            return "Multi-target strategy recommended based on compound analysis"
        
        return "Continue data collection for better insights"
    
    def execute_analysis(self, analysis_type: str, parameters: Dict) -> Dict:
        """Execute analysis code"""
        code_record = {
            'timestamp': datetime.now().isoformat(),
            'type': analysis_type,
            'parameters': parameters,
            'status': 'pending'
        }
        
        try:
            if analysis_type == 'resistance_prediction':
                result = self._predict_resistance(parameters)
            elif analysis_type == 'target_prioritization':
                result = self._prioritize_targets(parameters)
            elif analysis_type == 'compound_optimization':
                result = self._optimize_compound(parameters)
            else:
                result = {'error': 'Unknown analysis type'}
            
            code_record['status'] = 'completed'
            code_record['result'] = result
            
        except Exception as e:
            code_record['status'] = 'failed'
            code_record['error'] = str(e)
            logger.error(f"Analysis failed: {e}")
        
        self.code_history.append(code_record)
        return code_record
    
    def _predict_resistance(self, params: Dict) -> Dict:
        """Predict resistance for compound"""
        compound_id = params.get('compound_id')
        compound_state = self.state_tracker.get_compound_state(compound_id)
        
        if not compound_state:
            return {'error': 'Compound not found'}
        
        # Simple prediction based on state
        resistance_score = 0.5
        if 'assay_results' in compound_state:
            resistance_score = 0.3  # Lower if tested
        
        return {
            'compound_id': compound_id,
            'resistance_score': resistance_score,
            'confidence': 0.7
        }
    
    def _prioritize_targets(self, params: Dict) -> List[Dict]:
        """Prioritize targets based on druggability and resistance"""
        targets = self.state_tracker.query_state('target')
        
        prioritized = sorted(
            targets,
            key=lambda t: t.get('druggability', 0),
            reverse=True
        )
        
        return prioritized[:5]
    
    def _optimize_compound(self, params: Dict) -> Dict:
        """Suggest compound optimizations"""
        compound_id = params.get('compound_id')
        
        return {
            'compound_id': compound_id,
            'suggestions': [
                'Increase lipophilicity',
                'Add polar group for solubility',
                'Optimize for target binding'
            ]
        }
    
    def generate_report(self) -> Dict:
        """Generate comprehensive report (Kosmos-style)"""
        report = {
            'timestamp': datetime.now().isoformat(),
            'summary': self.state_tracker.get_state_summary(),
            'observations': len(self.observations),
            'analyses_run': len(self.code_history),
            'memory_size': len(self.memory),
            'key_findings': self._extract_key_findings(),
            'recommendations': self._generate_recommendations()
        }
        
        return report
    
    def _extract_key_findings(self) -> List[str]:
        """Extract key findings from memory"""
        findings = []
        
        # Check for multi-target compounds
        multi_target = self.knowledge_graph.find_multi_target_compounds()
        if multi_target:
            findings.append(f"Identified {len(multi_target)} multi-target compounds")
        
        # Check for validated hypotheses
        validated = self.state_tracker.get_validated_hypotheses()
        if validated:
            findings.append(f"Validated {len(validated)} hypotheses")
        
        return findings
    
    def _generate_recommendations(self) -> List[str]:
        """Generate actionable recommendations"""
        recommendations = []
        
        # Based on state
        summary = self.state_tracker.get_state_summary()
        
        if summary['compounds'] < 10:
            recommendations.append("Generate more compound candidates")
        
        if summary['assays'] < 5:
            recommendations.append("Conduct more biological assays")
        
        if summary['hypotheses']['proposed'] > 0:
            recommendations.append("Test pending hypotheses")
        
        return recommendations
