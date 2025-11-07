"""Dynamic PDB structure downloader - agents decide what to download"""

import logging
import requests
from pathlib import Path
from typing import Dict, List, Optional
from Bio.PDB import PDBList

logger = logging.getLogger('khukuri')


class StructureDownloader:
    """Download PDB structures on-demand based on agent decisions"""
    
    def __init__(self, data_dir: Optional[str] = None):
        self.data_dir = Path(data_dir) if data_dir else Path(__file__).parent.parent.parent / "data" / "structures"
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.pdbl = PDBList()
    
    def search_pdb_by_protein(self, protein_name: str, organism: str = None) -> List[Dict]:
        """Search PDB for protein structures"""
        try:
            search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
            
            query = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_words",
                        "value": protein_name
                    }
                },
                "return_type": "entry",
                "request_options": {
                    "results_content_type": ["experimental"],
                    "sort": [{"sort_by": "score", "direction": "desc"}]
                }
            }
            
            if organism:
                query["query"] = {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        query["query"],
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "rcsb_entity_source_organism.scientific_name",
                                "operator": "contains_words",
                                "value": organism
                            }
                        }
                    ]
                }
            
            response = requests.post(search_url, json=query, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                results = []
                
                for entry in data.get('result_set', [])[:5]:  # Top 5 results
                    pdb_id = entry.get('identifier')
                    if pdb_id:
                        results.append({
                            'pdb_id': pdb_id,
                            'score': entry.get('score', 0)
                        })
                
                logger.info(f"Found {len(results)} structures for {protein_name}")
                return results
            
            return []
            
        except Exception as e:
            logger.error(f"PDB search failed: {e}")
            return []
    
    def download_structure(self, pdb_id: str) -> Optional[Path]:
        """Download single PDB structure"""
        try:
            pdb_file = self.data_dir / f"{pdb_id.upper()}.pdb"
            
            # Check if already exists
            if pdb_file.exists():
                logger.info(f"Structure {pdb_id} already exists")
                return pdb_file
            
            logger.info(f"Downloading {pdb_id}...")
            
            # Download using BioPython
            filename = self.pdbl.retrieve_pdb_file(
                pdb_id, 
                pdir=str(self.data_dir), 
                file_format='pdb',
                overwrite=False
            )
            
            # Rename to clean format
            if filename:
                old_path = Path(filename)
                if old_path.exists():
                    old_path.rename(pdb_file)
                    logger.info(f"✓ Downloaded {pdb_id}")
                    return pdb_file
            
            return None
            
        except Exception as e:
            logger.error(f"Failed to download {pdb_id}: {e}")
            return None
    
    def download_for_target(self, target_name: str, organism: str = None, pdb_id: str = None) -> Optional[Path]:
        """Download structure for target - search if PDB ID not provided"""
        if pdb_id:
            return self.download_structure(pdb_id)
        
        # Search PDB for target
        results = self.search_pdb_by_protein(target_name, organism)
        
        if results:
            best_match = results[0]
            logger.info(f"Best match for {target_name}: {best_match['pdb_id']}")
            return self.download_structure(best_match['pdb_id'])
        
        logger.warning(f"No structures found for {target_name}")
        return None
    
    def download_batch(self, targets: List[Dict]) -> Dict[str, Path]:
        """Download multiple structures from agent decisions
        
        Args:
            targets: List of dicts with 'name', optional 'pdb_id', 'organism'
        
        Returns:
            Dict mapping target names to downloaded file paths
        """
        results = {}
        
        for target in targets:
            name = target.get('name')
            pdb_id = target.get('pdb_id')
            organism = target.get('organism')
            
            if not name:
                continue
            
            path = self.download_for_target(name, organism, pdb_id)
            if path:
                results[name] = path
        
        logger.info(f"Downloaded {len(results)}/{len(targets)} structures")
        return results
    
    def get_available_structures(self) -> List[str]:
        """List downloaded structures"""
        return [f.stem for f in self.data_dir.glob("*.pdb")]
    
    def get_structure_path(self, pdb_id: str) -> Optional[Path]:
        """Get path to structure, download if needed"""
        pdb_file = self.data_dir / f"{pdb_id.upper()}.pdb"
        
        if pdb_file.exists():
            return pdb_file
        
        return self.download_structure(pdb_id)
