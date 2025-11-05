"""Fetch data from external APIs"""

import logging
import requests
from time import sleep

logger = logging.getLogger('khukuri')


class DataFetcher:
    """Fetch omics and clinical data"""
    
    def __init__(self, timeout=10, max_retries=3):
        self.timeout = timeout
        self.max_retries = max_retries
    
    def fetch_string_ppi(self, species_name):
        """Fetch PPI data from STRING database"""
        try:
            species_id = self._get_string_species_id(species_name)
            if not species_id:
                logger.error(f"Species not found: {species_name}")
                return {}
            
            url = f"https://string-db.org/api/json/network?identifiers={species_name}&species={species_id}"
            response = requests.get(url, timeout=self.timeout)
            
            if response.status_code != 200:
                logger.error(f"STRING API error: {response.status_code}")
                return {}
            
            ppi_data = response.json()
            interactions = {}
            
            for edge in ppi_data:
                try:
                    score = float(edge.get('score', 0))
                    if score > 0.4:
                        interactions[(edge['preferredName_A'], edge['preferredName_B'])] = score
                except (ValueError, TypeError, KeyError):
                    continue
            
            logger.info(f"Fetched {len(interactions)} PPI interactions")
            return interactions
            
        except Exception as e:
            logger.error(f"STRING fetch failed: {e}")
            return {}
    
    def _get_string_species_id(self, species_name):
        """Get STRING species ID"""
        try:
            url = f"https://string-db.org/api/json/resolve?identifier={species_name.replace(' ', '%20')}"
            response = requests.get(url, timeout=self.timeout)
            if response.status_code == 200 and response.json():
                return response.json()[0]['taxonId']
        except:
            pass
        return None
    
    def fetch_kegg_pathways(self, organism_code='hsa'):
        """Fetch pathways from KEGG"""
        try:
            url = f"https://rest.kegg.jp/list/pathway/{organism_code}"
            response = requests.get(url, timeout=self.timeout)
            
            if response.status_code != 200:
                return {}
            
            pathways = {}
            for line in response.text.strip().split('\n'):
                if '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        pathways[parts[0].strip()] = parts[1].strip()
            
            logger.info(f"Fetched {len(pathways)} KEGG pathways")
            return pathways
            
        except Exception as e:
            logger.error(f"KEGG fetch failed: {e}")
            return {}
