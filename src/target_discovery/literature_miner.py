"""Mine scientific literature"""

import logging
from Bio import Entrez, Medline

logger = logging.getLogger('khukuri')


class LiteratureMiner:
    """Mine PubMed for target validation"""
    
    def __init__(self, email='khukuri@example.com'):
        Entrez.email = email
    
    def search_pubmed(self, query, max_results=50):
        """Search PubMed"""
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            record = Entrez.read(handle)
            handle.close()
            return record.get("IdList", [])
        except Exception as e:
            logger.error(f"PubMed search failed: {e}")
            return []
    
    def fetch_abstracts(self, pmid_list):
        """Fetch abstracts for PMIDs"""
        if not pmid_list:
            return []
        
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid_list, rettype="medline", retmode="text")
            records = Medline.parse(handle)
            
            abstracts = []
            for record in records:
                if 'AB' in record:
                    abstracts.append({
                        'pmid': record.get('PMID', ''),
                        'title': record.get('TI', ''),
                        'abstract': record.get('AB', ''),
                        'year': record.get('DP', '')
                    })
            handle.close()
            return abstracts
            
        except Exception as e:
            logger.error(f"Abstract fetch failed: {e}")
            return []
    
    def extract_evidence(self, query, max_results=50):
        """Extract evidence from literature"""
        pmids = self.search_pubmed(query, max_results)
        abstracts = self.fetch_abstracts(pmids)
        logger.info(f"Found {len(abstracts)} relevant papers")
        return abstracts
