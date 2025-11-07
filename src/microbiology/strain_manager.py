"""Manage bacterial strain information"""

import logging
from typing import Dict, List, Optional

logger = logging.getLogger('khukuri')


class StrainManager:
    """Manage bacterial strain collection"""
    
    def __init__(self):
        self.strains = {}
        self._initialize_reference_strains()
    
    def _initialize_reference_strains(self):
        """Initialize common reference strains"""
        reference_strains = {
            'ATCC_25923': {
                'species': 'Staphylococcus aureus',
                'type': 'reference',
                'resistance_profile': 'susceptible',
                'description': 'Quality control strain'
            },
            'ATCC_BAA-1556': {
                'species': 'Staphylococcus aureus',
                'type': 'reference',
                'resistance_profile': 'MRSA',
                'description': 'Methicillin-resistant reference'
            },
            'ATCC_25922': {
                'species': 'Escherichia coli',
                'type': 'reference',
                'resistance_profile': 'susceptible',
                'description': 'Quality control strain'
            },
            'ATCC_BAA-2452': {
                'species': 'Klebsiella pneumoniae',
                'type': 'reference',
                'resistance_profile': 'NDM-1',
                'description': 'Carbapenem-resistant reference'
            }
        }
        
        for strain_id, info in reference_strains.items():
            self.add_strain(strain_id, info)
    
    def add_strain(self, strain_id: str, strain_info: Dict):
        """Add strain to collection"""
        self.strains[strain_id] = {
            'strain_id': strain_id,
            **strain_info
        }
        logger.info(f"Added strain: {strain_id}")
    
    def get_strain(self, strain_id: str) -> Optional[Dict]:
        """Get strain information"""
        return self.strains.get(strain_id)
    
    def get_strains_by_species(self, species: str) -> List[Dict]:
        """Get all strains of species"""
        return [strain for strain in self.strains.values() 
                if species.lower() in strain.get('species', '').lower()]
    
    def get_resistant_strains(self) -> List[Dict]:
        """Get all resistant strains"""
        return [strain for strain in self.strains.values() 
                if strain.get('resistance_profile', 'susceptible') != 'susceptible']
    
    def get_reference_strains(self) -> List[Dict]:
        """Get reference strains"""
        return [strain for strain in self.strains.values() 
                if strain.get('type') == 'reference']
    
    def update_strain(self, strain_id: str, updates: Dict):
        """Update strain information"""
        if strain_id in self.strains:
            self.strains[strain_id].update(updates)
            logger.info(f"Updated strain: {strain_id}")
        else:
            logger.warning(f"Strain not found: {strain_id}")
