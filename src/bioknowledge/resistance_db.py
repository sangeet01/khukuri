"""Resistance gene and mechanism database (CARD, ResFinder, MEGARes integration)"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger('khukuri')


class ResistanceDatabase:
    """Curated AMR gene and mechanism database"""
    
    def __init__(self, db_path: Optional[str] = None):
        self.db_path = db_path or Path(__file__).parent.parent.parent / "config" / "resistance_db.json"
        self.genes = {}
        self.mechanisms = {}
        self.load_database()
    
    def load_database(self):
        """Load resistance database"""
        try:
            if Path(self.db_path).exists():
                with open(self.db_path) as f:
                    data = json.load(f)
                    self.genes = data.get('genes', {})
                    self.mechanisms = data.get('mechanisms', {})
            else:
                self._initialize_default_db()
        except Exception as e:
            logger.warning(f"Failed to load resistance DB: {e}, using defaults")
            self._initialize_default_db()
    
    def _initialize_default_db(self):
        """Initialize with common AMR genes"""
        self.genes = {
            'mecA': {'type': 'beta_lactam', 'organism': 'S. aureus', 'mechanism': 'target_alteration'},
            'vanA': {'type': 'glycopeptide', 'organism': 'Enterococcus', 'mechanism': 'target_alteration'},
            'blaCTX-M': {'type': 'beta_lactam', 'organism': 'E. coli', 'mechanism': 'antibiotic_inactivation'},
            'blaNDM': {'type': 'carbapenem', 'organism': 'Gram-negative', 'mechanism': 'antibiotic_inactivation'},
            'tetM': {'type': 'tetracycline', 'organism': 'Multiple', 'mechanism': 'efflux_pump'},
            'ermB': {'type': 'macrolide', 'organism': 'Multiple', 'mechanism': 'target_modification'},
            'aac(6\')': {'type': 'aminoglycoside', 'organism': 'Multiple', 'mechanism': 'antibiotic_inactivation'}
        }
        
        self.mechanisms = {
            'target_alteration': ['PBP modification', 'ribosome modification', 'DNA gyrase mutation'],
            'antibiotic_inactivation': ['beta-lactamase', 'aminoglycoside-modifying enzyme', 'chloramphenicol acetyltransferase'],
            'efflux_pump': ['MexAB-OprM', 'AcrAB-TolC', 'Tet efflux'],
            'reduced_permeability': ['porin loss', 'LPS modification']
        }
    
    def query_gene(self, gene_name: str) -> Optional[Dict]:
        """Query resistance gene information"""
        return self.genes.get(gene_name)
    
    def query_mechanism(self, mechanism: str) -> List[str]:
        """Query resistance mechanisms"""
        return self.mechanisms.get(mechanism, [])
    
    def get_genes_by_organism(self, organism: str) -> List[str]:
        """Get resistance genes for specific organism"""
        return [gene for gene, info in self.genes.items() 
                if organism.lower() in info.get('organism', '').lower()]
    
    def get_genes_by_drug_class(self, drug_class: str) -> List[str]:
        """Get resistance genes for drug class"""
        return [gene for gene, info in self.genes.items() 
                if drug_class.lower() in info.get('type', '').lower()]
    
    def add_gene(self, gene_name: str, gene_info: Dict):
        """Add new resistance gene"""
        self.genes[gene_name] = gene_info
        logger.info(f"Added resistance gene: {gene_name}")
    
    def save_database(self):
        """Save database to file"""
        try:
            with open(self.db_path, 'w') as f:
                json.dump({'genes': self.genes, 'mechanisms': self.mechanisms}, f, indent=2)
            logger.info(f"Saved resistance database to {self.db_path}")
        except Exception as e:
            logger.error(f"Failed to save database: {e}")
