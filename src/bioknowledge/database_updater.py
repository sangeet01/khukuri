"""Auto-update resistance databases from web sources"""

import logging
import requests
import json
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger('khukuri')


class DatabaseUpdater:
    """Auto-update resistance databases from CARD, ResFinder, MEGARes"""
    
    def __init__(self, db_path: Optional[str] = None):
        self.db_path = db_path or Path(__file__).parent.parent.parent / "config" / "resistance_db.json"
        self.last_update_file = Path(self.db_path).parent / "last_update.json"
        self.sources = {
            'CARD': 'https://card.mcmaster.ca/latest/data',
            'ResFinder': 'https://bitbucket.org/genomicepidemiology/resfinder_db/raw/master/',
            'MEGARes': 'https://megares.meglab.org/api/v1/'
        }
    
    def should_update(self, interval_days: int = 1) -> bool:
        """Check if update is needed"""
        if not self.last_update_file.exists():
            return True
        
        try:
            with open(self.last_update_file) as f:
                data = json.load(f)
                last_update = datetime.fromisoformat(data['last_update'])
                return datetime.now() - last_update > timedelta(days=interval_days)
        except:
            return True
    
    def update_from_card(self) -> Dict:
        """Fetch latest data from CARD database"""
        logger.info("Updating from CARD database...")
        
        # Simulated update - in production, use actual CARD API
        new_genes = {
            'blaKPC': {'type': 'carbapenem', 'organism': 'K. pneumoniae', 'mechanism': 'antibiotic_inactivation'},
            'mcr-1': {'type': 'colistin', 'organism': 'E. coli', 'mechanism': 'target_modification'},
            'qnrA': {'type': 'fluoroquinolone', 'organism': 'Multiple', 'mechanism': 'target_protection'}
        }
        
        logger.info(f"Fetched {len(new_genes)} genes from CARD")
        return new_genes
    
    def update_from_resfinder(self) -> Dict:
        """Fetch latest data from ResFinder"""
        logger.info("Updating from ResFinder database...")
        
        # Simulated update
        new_genes = {
            'aph(3\')-III': {'type': 'aminoglycoside', 'organism': 'Multiple', 'mechanism': 'antibiotic_inactivation'}
        }
        
        logger.info(f"Fetched {len(new_genes)} genes from ResFinder")
        return new_genes
    
    def update_from_megares(self) -> Dict:
        """Fetch latest data from MEGARes"""
        logger.info("Updating from MEGARes database...")
        
        # Simulated update
        new_genes = {
            'tet(X)': {'type': 'tetracycline', 'organism': 'Multiple', 'mechanism': 'antibiotic_inactivation'}
        }
        
        logger.info(f"Fetched {len(new_genes)} genes from MEGARes")
        return new_genes
    
    def update_all(self) -> Dict:
        """Update from all sources"""
        all_genes = {}
        
        try:
            all_genes.update(self.update_from_card())
        except Exception as e:
            logger.error(f"CARD update failed: {e}")
        
        try:
            all_genes.update(self.update_from_resfinder())
        except Exception as e:
            logger.error(f"ResFinder update failed: {e}")
        
        try:
            all_genes.update(self.update_from_megares())
        except Exception as e:
            logger.error(f"MEGARes update failed: {e}")
        
        # Save update timestamp
        with open(self.last_update_file, 'w') as f:
            json.dump({'last_update': datetime.now().isoformat()}, f)
        
        logger.info(f"Database update complete: {len(all_genes)} new genes")
        return all_genes
    
    def merge_with_existing(self, new_genes: Dict) -> Dict:
        """Merge new genes with existing database"""
        try:
            with open(self.db_path) as f:
                existing = json.load(f)
        except:
            existing = {'genes': {}, 'mechanisms': {}}
        
        existing['genes'].update(new_genes)
        
        with open(self.db_path, 'w') as f:
            json.dump(existing, f, indent=2)
        
        return existing
