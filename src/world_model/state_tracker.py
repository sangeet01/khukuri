"""Track world state for AMR discovery (Kosmos-style)"""

import logging
from typing import Dict, List, Optional, Any
from datetime import datetime
import json

logger = logging.getLogger('khukuri')


class WorldStateTracker:
    """Track complete state of AMR discovery process"""
    
    def __init__(self):
        self.compounds = {}
        self.targets = {}
        self.strains = {}
        self.assays = {}
        self.hypotheses = []
        self.experiments = []
        self.state_history = []
    
    def update_compound(self, compound_id: str, data: Dict):
        """Update compound state"""
        if compound_id not in self.compounds:
            self.compounds[compound_id] = {'compound_id': compound_id, 'history': []}
        
        self.compounds[compound_id].update(data)
        self.compounds[compound_id]['history'].append({
            'timestamp': datetime.now().isoformat(),
            'update': data
        })
        
        self._record_state_change('compound_update', compound_id, data)
    
    def update_target(self, target_id: str, data: Dict):
        """Update target state"""
        if target_id not in self.targets:
            self.targets[target_id] = {'target_id': target_id, 'history': []}
        
        self.targets[target_id].update(data)
        self.targets[target_id]['history'].append({
            'timestamp': datetime.now().isoformat(),
            'update': data
        })
        
        self._record_state_change('target_update', target_id, data)
    
    def update_strain(self, strain_id: str, data: Dict):
        """Update strain state"""
        if strain_id not in self.strains:
            self.strains[strain_id] = {'strain_id': strain_id, 'history': []}
        
        self.strains[strain_id].update(data)
        self.strains[strain_id]['history'].append({
            'timestamp': datetime.now().isoformat(),
            'update': data
        })
        
        self._record_state_change('strain_update', strain_id, data)
    
    def add_assay_result(self, assay_id: str, compound_id: str, strain_id: str, results: Dict):
        """Add assay result to world state"""
        self.assays[assay_id] = {
            'assay_id': assay_id,
            'compound_id': compound_id,
            'strain_id': strain_id,
            'results': results,
            'timestamp': datetime.now().isoformat()
        }
        
        # Update compound with assay result
        if compound_id in self.compounds:
            if 'assay_results' not in self.compounds[compound_id]:
                self.compounds[compound_id]['assay_results'] = []
            self.compounds[compound_id]['assay_results'].append(assay_id)
        
        self._record_state_change('assay_result', assay_id, results)
    
    def add_hypothesis(self, hypothesis: str, evidence: List[str], confidence: float):
        """Add scientific hypothesis"""
        hyp = {
            'hypothesis': hypothesis,
            'evidence': evidence,
            'confidence': confidence,
            'timestamp': datetime.now().isoformat(),
            'status': 'proposed'
        }
        self.hypotheses.append(hyp)
        logger.info(f"Added hypothesis: {hypothesis[:50]}...")
    
    def update_hypothesis_status(self, hypothesis_idx: int, status: str, results: Optional[Dict] = None):
        """Update hypothesis status (validated/rejected)"""
        if 0 <= hypothesis_idx < len(self.hypotheses):
            self.hypotheses[hypothesis_idx]['status'] = status
            if results:
                self.hypotheses[hypothesis_idx]['results'] = results
            logger.info(f"Updated hypothesis {hypothesis_idx} status: {status}")
    
    def add_experiment(self, experiment_type: str, parameters: Dict, results: Optional[Dict] = None):
        """Add experiment to history"""
        exp = {
            'experiment_type': experiment_type,
            'parameters': parameters,
            'results': results,
            'timestamp': datetime.now().isoformat()
        }
        self.experiments.append(exp)
        self._record_state_change('experiment', experiment_type, parameters)
    
    def get_compound_state(self, compound_id: str) -> Optional[Dict]:
        """Get complete compound state"""
        return self.compounds.get(compound_id)
    
    def get_target_state(self, target_id: str) -> Optional[Dict]:
        """Get complete target state"""
        return self.targets.get(target_id)
    
    def get_active_hypotheses(self) -> List[Dict]:
        """Get all active hypotheses"""
        return [h for h in self.hypotheses if h['status'] == 'proposed']
    
    def get_validated_hypotheses(self) -> List[Dict]:
        """Get validated hypotheses"""
        return [h for h in self.hypotheses if h['status'] == 'validated']
    
    def query_state(self, entity_type: str, filters: Optional[Dict] = None) -> List[Dict]:
        """Query world state"""
        if entity_type == 'compound':
            entities = self.compounds.values()
        elif entity_type == 'target':
            entities = self.targets.values()
        elif entity_type == 'strain':
            entities = self.strains.values()
        elif entity_type == 'assay':
            entities = self.assays.values()
        else:
            return []
        
        if not filters:
            return list(entities)
        
        # Apply filters
        results = []
        for entity in entities:
            match = all(entity.get(k) == v for k, v in filters.items())
            if match:
                results.append(entity)
        
        return results
    
    def _record_state_change(self, change_type: str, entity_id: str, data: Any):
        """Record state change in history"""
        self.state_history.append({
            'timestamp': datetime.now().isoformat(),
            'change_type': change_type,
            'entity_id': entity_id,
            'data': data
        })
    
    def get_state_summary(self) -> Dict:
        """Get summary of current world state"""
        return {
            'compounds': len(self.compounds),
            'targets': len(self.targets),
            'strains': len(self.strains),
            'assays': len(self.assays),
            'hypotheses': {
                'total': len(self.hypotheses),
                'proposed': len([h for h in self.hypotheses if h['status'] == 'proposed']),
                'validated': len([h for h in self.hypotheses if h['status'] == 'validated']),
                'rejected': len([h for h in self.hypotheses if h['status'] == 'rejected'])
            },
            'experiments': len(self.experiments),
            'state_changes': len(self.state_history)
        }
    
    def export_state(self, filepath: str):
        """Export world state to file"""
        state = {
            'compounds': self.compounds,
            'targets': self.targets,
            'strains': self.strains,
            'assays': self.assays,
            'hypotheses': self.hypotheses,
            'experiments': self.experiments
        }
        
        with open(filepath, 'w') as f:
            json.dump(state, f, indent=2)
        
        logger.info(f"Exported world state to {filepath}")
