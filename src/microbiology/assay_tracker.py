"""Track biological assay results"""

import logging
from typing import Dict, List, Optional
from datetime import datetime

logger = logging.getLogger('khukuri')


class AssayTracker:
    """Track and manage biological assay results"""
    
    def __init__(self):
        self.assays = {}
        self.assay_counter = 0
    
    def add_assay(self, assay_type: str, compound_id: str, strain_id: str, 
                  results: Dict, metadata: Optional[Dict] = None) -> str:
        """Add assay result"""
        self.assay_counter += 1
        assay_id = f"ASSAY_{self.assay_counter:06d}"
        
        self.assays[assay_id] = {
            'assay_id': assay_id,
            'assay_type': assay_type,
            'compound_id': compound_id,
            'strain_id': strain_id,
            'results': results,
            'metadata': metadata or {},
            'timestamp': datetime.now().isoformat()
        }
        
        logger.info(f"Added {assay_type} assay: {assay_id}")
        return assay_id
    
    def get_assay(self, assay_id: str) -> Optional[Dict]:
        """Get assay by ID"""
        return self.assays.get(assay_id)
    
    def get_compound_assays(self, compound_id: str) -> List[Dict]:
        """Get all assays for compound"""
        return [assay for assay in self.assays.values() 
                if assay['compound_id'] == compound_id]
    
    def get_strain_assays(self, strain_id: str) -> List[Dict]:
        """Get all assays for strain"""
        return [assay for assay in self.assays.values() 
                if assay['strain_id'] == strain_id]
    
    def get_assays_by_type(self, assay_type: str) -> List[Dict]:
        """Get assays by type"""
        return [assay for assay in self.assays.values() 
                if assay['assay_type'] == assay_type]
    
    def summarize_compound_activity(self, compound_id: str) -> Dict:
        """Summarize activity across all assays"""
        assays = self.get_compound_assays(compound_id)
        
        summary = {
            'compound_id': compound_id,
            'total_assays': len(assays),
            'assay_types': {},
            'active_count': 0,
            'inactive_count': 0
        }
        
        for assay in assays:
            assay_type = assay['assay_type']
            summary['assay_types'][assay_type] = summary['assay_types'].get(assay_type, 0) + 1
            
            # Determine activity
            if assay['results'].get('active', False):
                summary['active_count'] += 1
            else:
                summary['inactive_count'] += 1
        
        if summary['total_assays'] > 0:
            summary['activity_rate'] = summary['active_count'] / summary['total_assays']
        else:
            summary['activity_rate'] = 0.0
        
        return summary
