"""Automated data updater - fetches latest data from web daily"""

import json
import requests
from datetime import datetime, timedelta
from pathlib import Path
from .logger import setup_logger

logger = setup_logger(__name__)


class DataUpdater:
    """Automatically update local databases from web sources"""
    
    def __init__(self, data_dir="data"):
        self.data_dir = Path(data_dir)
        self.metadata_file = self.data_dir / "update_metadata.json"
        self.sources = {
            "resistance_mutations": self._fetch_resistance_data,
            "drug_targets": self._fetch_drug_targets,
            "clinical_trials": self._fetch_clinical_trials,
            "pdb_structures": self._fetch_pdb_updates
        }
    
    def should_update(self, dataset_name):
        """Check if dataset needs update (>24 hours old)"""
        if not self.metadata_file.exists():
            return True
        
        try:
            with open(self.metadata_file) as f:
                metadata = json.load(f)
            
            last_update = metadata.get(dataset_name, {}).get("last_update")
            if not last_update:
                return True
            
            last_time = datetime.fromisoformat(last_update)
            return datetime.now() - last_time > timedelta(hours=24)
        except Exception as e:
            logger.warning(f"Error reading metadata: {e}")
            return True
    
    def update_all(self, force=False):
        """Update all datasets if needed"""
        logger.info("Starting data update check...")
        results = {}
        
        for dataset_name, fetch_func in self.sources.items():
            if force or self.should_update(dataset_name):
                logger.info(f"Updating {dataset_name}...")
                try:
                    data = fetch_func()
                    self._save_dataset(dataset_name, data)
                    results[dataset_name] = "updated"
                except Exception as e:
                    logger.error(f"Failed to update {dataset_name}: {e}")
                    results[dataset_name] = f"failed: {e}"
            else:
                logger.info(f"{dataset_name} is up to date")
                results[dataset_name] = "skipped"
        
        return results
    
    def _fetch_resistance_data(self):
        """Fetch resistance mutations from CARD database"""
        logger.info("Fetching resistance data from CARD...")
        
        # CARD API endpoint
        url = "https://card.mcmaster.ca/latest/data"
        
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                # Parse and structure data
                raw_data = response.json()
                structured = self._parse_card_data(raw_data)
                return structured
        except Exception as e:
            logger.warning(f"CARD fetch failed: {e}, using fallback")
        
        # Fallback: fetch from PubMed recent papers
        return self._fetch_resistance_from_pubmed()
    
    def _fetch_drug_targets(self):
        """Fetch validated drug targets from DrugBank/ChEMBL"""
        logger.info("Fetching drug targets...")
        
        targets = []
        
        # ChEMBL API
        try:
            url = "https://www.ebi.ac.uk/chembl/api/data/target.json"
            params = {"limit": 100, "target_type": "SINGLE PROTEIN"}
            response = requests.get(url, params=params, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                for target in data.get("targets", []):
                    targets.append({
                        "id": target.get("target_chembl_id"),
                        "name": target.get("pref_name"),
                        "organism": target.get("organism"),
                        "type": target.get("target_type")
                    })
        except Exception as e:
            logger.warning(f"ChEMBL fetch failed: {e}")
        
        return {"targets": targets, "count": len(targets)}
    
    def _fetch_clinical_trials(self):
        """Fetch recent clinical trials from ClinicalTrials.gov"""
        logger.info("Fetching clinical trials...")
        
        try:
            url = "https://clinicaltrials.gov/api/query/study_fields"
            params = {
                "expr": "drug discovery OR antimicrobial",
                "fields": "NCTId,BriefTitle,Condition,InterventionName,Phase",
                "min_rnk": 1,
                "max_rnk": 100,
                "fmt": "json"
            }
            
            response = requests.get(url, params=params, timeout=30)
            if response.status_code == 200:
                data = response.json()
                trials = data.get("StudyFieldsResponse", {}).get("StudyFields", [])
                return {"trials": trials, "count": len(trials)}
        except Exception as e:
            logger.warning(f"ClinicalTrials fetch failed: {e}")
        
        return {"trials": [], "count": 0}
    
    def _fetch_pdb_updates(self):
        """Fetch recent PDB structure updates"""
        logger.info("Fetching PDB updates...")
        
        try:
            # Get recent entries from last week
            url = "https://www.rcsb.org/pdb/rest/getCurrent"
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                pdb_ids = response.text.strip().split("\n")[:50]  # Top 50
                return {"pdb_ids": pdb_ids, "count": len(pdb_ids)}
        except Exception as e:
            logger.warning(f"PDB fetch failed: {e}")
        
        return {"pdb_ids": [], "count": 0}
    
    def _fetch_resistance_from_pubmed(self):
        """Fallback: fetch resistance data from recent PubMed papers"""
        logger.info("Fetching resistance data from PubMed...")
        
        try:
            # Search recent papers on antimicrobial resistance
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {
                "db": "pubmed",
                "term": "antimicrobial resistance mutations",
                "retmax": 20,
                "retmode": "json",
                "sort": "date"
            }
            
            response = requests.get(search_url, params=params, timeout=30)
            if response.status_code == 200:
                data = response.json()
                pmids = data.get("esearchresult", {}).get("idlist", [])
                
                return {
                    "source": "pubmed",
                    "pmids": pmids,
                    "count": len(pmids),
                    "note": "Manual curation needed"
                }
        except Exception as e:
            logger.error(f"PubMed fetch failed: {e}")
        
        return {"source": "none", "count": 0}
    
    def _parse_card_data(self, raw_data):
        """Parse CARD database format"""
        # Simplified parser - adapt to actual CARD format
        mutations = {}
        
        for entry in raw_data.get("data", []):
            gene = entry.get("gene_name")
            if gene:
                mutations[gene] = {
                    "mutations": entry.get("mutations", []),
                    "resistance": entry.get("resistance_type"),
                    "organism": entry.get("organism")
                }
        
        return mutations
    
    def _save_dataset(self, dataset_name, data):
        """Save dataset to JSON file"""
        output_file = self.data_dir / f"{dataset_name}.json"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        # Update metadata
        self._update_metadata(dataset_name)
        logger.info(f"Saved {dataset_name} to {output_file}")
    
    def _update_metadata(self, dataset_name):
        """Update metadata file with last update time"""
        metadata = {}
        if self.metadata_file.exists():
            with open(self.metadata_file) as f:
                metadata = json.load(f)
        
        metadata[dataset_name] = {
            "last_update": datetime.now().isoformat(),
            "status": "success"
        }
        
        with open(self.metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)


def update_databases(force=False):
    """Convenience function to update all databases"""
    updater = DataUpdater()
    return updater.update_all(force=force)


if __name__ == "__main__":
    # Run update
    results = update_databases(force=True)
    print("\nUpdate Results:")
    for dataset, status in results.items():
        print(f"  {dataset}: {status}")
