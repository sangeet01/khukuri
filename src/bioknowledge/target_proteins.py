"""Target protein database and analysis"""

import logging
from typing import Dict, List, Optional

logger = logging.getLogger('khukuri')


class TargetProteinDB:
    """Database of validated antibiotic targets"""
    
    def __init__(self):
        self.targets = self._initialize_targets()
    
    def _initialize_targets(self) -> Dict:
        """Initialize validated antibiotic targets"""
        return {
            'InhA': {
                'function': 'Enoyl-ACP reductase',
                'pathway': 'Fatty acid synthesis',
                'organism': 'M. tuberculosis',
                'druggability': 0.9,
                'known_inhibitors': ['isoniazid', 'ethionamide'],
                'binding_site': 'NAD+ binding pocket',
                'pdb_ids': ['1ENY', '2X22', '4TZK']
            },
            'DNA gyrase': {
                'function': 'DNA topoisomerase',
                'pathway': 'DNA replication',
                'organism': 'Multiple',
                'druggability': 0.85,
                'known_inhibitors': ['fluoroquinolones', 'novobiocin'],
                'binding_site': 'ATP binding site',
                'pdb_ids': ['1KZN', '2XCT']
            },
            'PBP2a': {
                'function': 'Penicillin-binding protein',
                'pathway': 'Cell wall synthesis',
                'organism': 'S. aureus',
                'druggability': 0.75,
                'known_inhibitors': ['ceftaroline', 'ceftobiprole'],
                'binding_site': 'Active site serine',
                'pdb_ids': ['1VQQ', '3ZG0']
            },
            'MurA': {
                'function': 'UDP-N-acetylglucosamine enolpyruvyl transferase',
                'pathway': 'Peptidoglycan synthesis',
                'organism': 'Multiple',
                'druggability': 0.8,
                'known_inhibitors': ['fosfomycin'],
                'binding_site': 'Substrate binding pocket',
                'pdb_ids': ['1UAE', '1EYN']
            },
            'FtsZ': {
                'function': 'Cell division protein',
                'pathway': 'Cytokinesis',
                'organism': 'Multiple',
                'druggability': 0.7,
                'known_inhibitors': ['PC190723', 'TXA709'],
                'binding_site': 'GTP binding site',
                'pdb_ids': ['2VAP', '3VOA']
            },
            'NDM-1': {
                'function': 'Metallo-beta-lactamase',
                'pathway': 'Antibiotic resistance',
                'organism': 'Gram-negative',
                'druggability': 0.65,
                'known_inhibitors': ['aspergillomarasmine A'],
                'binding_site': 'Zinc binding site',
                'pdb_ids': ['3SPU', '4EXS']
            }
        }
    
    def get_target_info(self, target_name: str) -> Optional[Dict]:
        """Get target protein information"""
        return self.targets.get(target_name)
    
    def get_druggable_targets(self, min_score: float = 0.7) -> List[str]:
        """Get targets above druggability threshold"""
        return [name for name, info in self.targets.items() 
                if info.get('druggability', 0) >= min_score]
    
    def get_targets_by_pathway(self, pathway: str) -> List[str]:
        """Get targets in specific pathway"""
        return [name for name, info in self.targets.items() 
                if pathway.lower() in info.get('pathway', '').lower()]
    
    def get_pdb_structures(self, target_name: str) -> List[str]:
        """Get PDB IDs for target"""
        info = self.targets.get(target_name, {})
        return info.get('pdb_ids', [])
    
    def add_target(self, name: str, target_info: Dict):
        """Add new target protein"""
        self.targets[name] = target_info
        logger.info(f"Added target protein: {name}")
