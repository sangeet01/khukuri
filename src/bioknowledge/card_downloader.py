"""Auto-download and parse CARD resistance database"""

import logging
import requests
import json
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger('khukuri')


class CARDDownloader:
    """Download and parse CARD AMR database"""
    
    def __init__(self, data_dir: Optional[str] = None):
        self.data_dir = Path(data_dir) if data_dir else Path(__file__).parent.parent.parent / "data" / "resistance"
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.db_file = self.data_dir / "card_database.json"
    
    def download_card(self) -> bool:
        """Download CARD database from API"""
        try:
            logger.info("Downloading CARD database...")
            
            # CARD prevalence data API
            url = "https://card.mcmaster.ca/latest/prevalence"
            response = requests.get(url, timeout=120)
            
            if response.status_code == 200:
                data = response.json()
                
                with open(self.db_file, 'w') as f:
                    json.dump(data, f, indent=2)
                
                logger.info(f"✓ Downloaded CARD database")
                return True
            else:
                logger.warning(f"CARD API returned {response.status_code}, using fallback")
                return False
                
        except Exception as e:
            logger.warning(f"CARD download failed: {e}, using fallback")
            return False
    
    def parse_card_to_resistance_db(self) -> Dict:
        """Parse CARD JSON to Khukuri resistance database format"""
        try:
            if not self.db_file.exists():
                logger.info("CARD database not found, downloading...")
                self.download_card()
            
            if self.db_file.exists():
                with open(self.db_file) as f:
                    card_data = json.load(f)
                
                resistance_genes = {}
                
                # Parse CARD structure
                for entry_id, entry_data in card_data.items():
                    if not isinstance(entry_data, dict):
                        continue
                    
                    gene_name = entry_data.get('model_name', entry_id)
                    
                    resistance_genes[gene_name] = {
                        'type': entry_data.get('drug_class', 'unknown'),
                        'organism': entry_data.get('organism', 'Multiple'),
                        'mechanism': entry_data.get('resistance_mechanism', 'unknown'),
                        'source': 'CARD'
                    }
                
                logger.info(f"Parsed {len(resistance_genes)} genes from CARD")
                return resistance_genes
            
            return self._get_fallback_db()
            
        except Exception as e:
            logger.error(f"CARD parsing failed: {e}")
            return self._get_fallback_db()
    
    def _get_fallback_db(self) -> Dict:
        """Comprehensive fallback database"""
        logger.info("Using curated fallback resistance database")
        return {
            # Beta-lactam resistance
            'mecA': {'type': 'beta_lactam', 'organism': 'S. aureus', 'mechanism': 'target_alteration', 'source': 'fallback'},
            'blaCTX-M': {'type': 'beta_lactam', 'organism': 'E. coli', 'mechanism': 'antibiotic_inactivation', 'source': 'fallback'},
            'blaNDM': {'type': 'carbapenem', 'organism': 'Gram-negative', 'mechanism': 'antibiotic_inactivation', 'source': 'fallback'},
            'blaKPC': {'type': 'carbapenem', 'organism': 'K. pneumoniae', 'mechanism': 'antibiotic_inactivation', 'source': 'fallback'},
            'blaOXA': {'type': 'carbapenem', 'organism': 'Acinetobacter', 'mechanism': 'antibiotic_inactivation', 'source': 'fallback'},
            
            # Glycopeptide resistance
            'vanA': {'type': 'glycopeptide', 'organism': 'Enterococcus', 'mechanism': 'target_alteration', 'source': 'fallback'},
            'vanB': {'type': 'glycopeptide', 'organism': 'Enterococcus', 'mechanism': 'target_alteration', 'source': 'fallback'},
            
            # TB resistance
            'inhA': {'type': 'isoniazid', 'organism': 'M. tuberculosis', 'mechanism': 'target_mutation', 'source': 'fallback'},
            'katG': {'type': 'isoniazid', 'organism': 'M. tuberculosis', 'mechanism': 'target_mutation', 'source': 'fallback'},
            'rpoB': {'type': 'rifampicin', 'organism': 'M. tuberculosis', 'mechanism': 'target_mutation', 'source': 'fallback'},
            'embB': {'type': 'ethambutol', 'organism': 'M. tuberculosis', 'mechanism': 'target_mutation', 'source': 'fallback'},
            'gyrA': {'type': 'fluoroquinolone', 'organism': 'M. tuberculosis', 'mechanism': 'target_mutation', 'source': 'fallback'},
            
            # Aminoglycoside resistance
            'aac(6\')': {'type': 'aminoglycoside', 'organism': 'Multiple', 'mechanism': 'antibiotic_inactivation', 'source': 'fallback'},
            'aph(3\')': {'type': 'aminoglycoside', 'organism': 'Multiple', 'mechanism': 'antibiotic_inactivation', 'source': 'fallback'},
            
            # Tetracycline resistance
            'tetM': {'type': 'tetracycline', 'organism': 'Multiple', 'mechanism': 'efflux_pump', 'source': 'fallback'},
            'tetA': {'type': 'tetracycline', 'organism': 'Multiple', 'mechanism': 'efflux_pump', 'source': 'fallback'},
            
            # Macrolide resistance
            'ermB': {'type': 'macrolide', 'organism': 'Multiple', 'mechanism': 'target_modification', 'source': 'fallback'},
            'mphA': {'type': 'macrolide', 'organism': 'Multiple', 'mechanism': 'antibiotic_inactivation', 'source': 'fallback'},
            
            # Fluoroquinolone resistance
            'qnrA': {'type': 'fluoroquinolone', 'organism': 'Multiple', 'mechanism': 'target_protection', 'source': 'fallback'},
            'qnrB': {'type': 'fluoroquinolone', 'organism': 'Multiple', 'mechanism': 'target_protection', 'source': 'fallback'},
            
            # Colistin resistance
            'mcr-1': {'type': 'colistin', 'organism': 'E. coli', 'mechanism': 'target_modification', 'source': 'fallback'},
            'mcr-2': {'type': 'colistin', 'organism': 'E. coli', 'mechanism': 'target_modification', 'source': 'fallback'}
        }
    
    def update_resistance_db_file(self, output_path: Optional[str] = None):
        """Update resistance_db.json with CARD data"""
        genes = self.parse_card_to_resistance_db()
        
        output_file = Path(output_path) if output_path else Path(__file__).parent.parent.parent / "config" / "resistance_db.json"
        
        db_data = {
            'genes': genes,
            'mechanisms': {
                'target_alteration': ['PBP modification', 'ribosome modification', 'DNA gyrase mutation'],
                'target_mutation': ['point mutation', 'deletion', 'insertion'],
                'antibiotic_inactivation': ['beta-lactamase', 'aminoglycoside-modifying enzyme', 'chloramphenicol acetyltransferase'],
                'efflux_pump': ['MexAB-OprM', 'AcrAB-TolC', 'Tet efflux'],
                'reduced_permeability': ['porin loss', 'LPS modification'],
                'target_protection': ['ribosomal protection', 'DNA gyrase protection'],
                'target_modification': ['methylation', 'phosphorylation']
            },
            'source': 'CARD' if self.db_file.exists() else 'fallback',
            'gene_count': len(genes)
        }
        
        with open(output_file, 'w') as f:
            json.dump(db_data, f, indent=2)
        
        logger.info(f"Updated {output_file} with {len(genes)} genes")
        return output_file
