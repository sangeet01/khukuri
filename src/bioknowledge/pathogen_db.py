"""Pathogen taxonomy and target protein database"""

import logging
from typing import Dict, List, Optional

logger = logging.getLogger('khukuri')


class PathogenDatabase:
    """Pathogen taxonomy and essential proteins"""
    
    def __init__(self):
        self.pathogens = self._initialize_pathogens()
        self.priority_pathogens = self._initialize_priority_list()
    
    def _initialize_pathogens(self) -> Dict:
        """Initialize pathogen database with WHO priority pathogens"""
        return {
            'Mycobacterium tuberculosis': {
                'taxonomy': {'kingdom': 'Bacteria', 'phylum': 'Actinobacteria', 'class': 'Actinomycetes'},
                'priority': 'critical',
                'essential_genes': ['inhA', 'katG', 'rpoB', 'embB', 'pncA'],
                'targets': ['InhA', 'KatG', 'RNA polymerase', 'Arabinosyltransferase'],
                'resistance_profile': ['rifampicin', 'isoniazid', 'pyrazinamide', 'ethambutol']
            },
            'Staphylococcus aureus': {
                'taxonomy': {'kingdom': 'Bacteria', 'phylum': 'Firmicutes', 'class': 'Bacilli'},
                'priority': 'high',
                'essential_genes': ['mecA', 'pbp2a', 'femA', 'murA'],
                'targets': ['PBP2a', 'FemA', 'MurA', 'DNA gyrase'],
                'resistance_profile': ['methicillin', 'vancomycin']
            },
            'Escherichia coli': {
                'taxonomy': {'kingdom': 'Bacteria', 'phylum': 'Proteobacteria', 'class': 'Gammaproteobacteria'},
                'priority': 'high',
                'essential_genes': ['ftsZ', 'murA', 'dnaA', 'gyrA'],
                'targets': ['FtsZ', 'MurA', 'DNA gyrase', 'Topoisomerase IV'],
                'resistance_profile': ['beta-lactam', 'fluoroquinolone', 'carbapenem']
            },
            'Pseudomonas aeruginosa': {
                'taxonomy': {'kingdom': 'Bacteria', 'phylum': 'Proteobacteria', 'class': 'Gammaproteobacteria'},
                'priority': 'critical',
                'essential_genes': ['oprD', 'mexAB', 'ampC', 'gyrA'],
                'targets': ['OprD', 'MexAB-OprM', 'AmpC', 'DNA gyrase'],
                'resistance_profile': ['carbapenem', 'fluoroquinolone', 'aminoglycoside']
            },
            'Acinetobacter baumannii': {
                'taxonomy': {'kingdom': 'Bacteria', 'phylum': 'Proteobacteria', 'class': 'Gammaproteobacteria'},
                'priority': 'critical',
                'essential_genes': ['blaOXA', 'adeABC', 'ompA'],
                'targets': ['OXA carbapenemase', 'AdeABC efflux', 'OmpA'],
                'resistance_profile': ['carbapenem', 'colistin']
            },
            'Klebsiella pneumoniae': {
                'taxonomy': {'kingdom': 'Bacteria', 'phylum': 'Proteobacteria', 'class': 'Gammaproteobacteria'},
                'priority': 'critical',
                'essential_genes': ['blaNDM', 'blaKPC', 'ompK36'],
                'targets': ['NDM-1', 'KPC', 'OmpK36'],
                'resistance_profile': ['carbapenem', 'colistin']
            }
        }
    
    def _initialize_priority_list(self) -> Dict:
        """WHO priority pathogen list"""
        return {
            'critical': ['Acinetobacter baumannii', 'Pseudomonas aeruginosa', 'Klebsiella pneumoniae'],
            'high': ['Staphylococcus aureus', 'Helicobacter pylori', 'Salmonella spp'],
            'medium': ['Streptococcus pneumoniae', 'Haemophilus influenzae']
        }
    
    def get_pathogen_info(self, pathogen_name: str) -> Optional[Dict]:
        """Get pathogen information"""
        return self.pathogens.get(pathogen_name)
    
    def get_essential_genes(self, pathogen_name: str) -> List[str]:
        """Get essential genes for pathogen"""
        info = self.pathogens.get(pathogen_name, {})
        return info.get('essential_genes', [])
    
    def get_target_proteins(self, pathogen_name: str) -> List[str]:
        """Get target proteins for pathogen"""
        info = self.pathogens.get(pathogen_name, {})
        return info.get('targets', [])
    
    def get_priority_pathogens(self, priority_level: str = 'critical') -> List[str]:
        """Get pathogens by priority level"""
        return self.priority_pathogens.get(priority_level, [])
    
    def get_resistance_profile(self, pathogen_name: str) -> List[str]:
        """Get known resistance profile"""
        info = self.pathogens.get(pathogen_name, {})
        return info.get('resistance_profile', [])
