"""Track and analyze resistance mutations"""

import logging
from typing import Dict, List
from collections import defaultdict

logger = logging.getLogger('khukuri')


class MutationTracker:
    """Track resistance mutations across strains"""
    
    def __init__(self):
        self.mutations = defaultdict(list)
        self.strain_profiles = {}
    
    def add_strain(self, strain_id: str, mutations: Dict[str, List[str]]):
        """Add strain mutation profile"""
        self.strain_profiles[strain_id] = mutations
        
        for gene, muts in mutations.items():
            for mut in muts:
                self.mutations[gene].append({
                    'strain': strain_id,
                    'mutation': mut
                })
        
        logger.info(f"Added strain {strain_id} with {sum(len(m) for m in mutations.values())} mutations")
    
    def get_mutation_frequency(self, gene: str, mutation: str) -> float:
        """Calculate mutation frequency across strains"""
        if not self.strain_profiles:
            return 0.0
        
        count = sum(1 for strain_id, muts in self.strain_profiles.items() 
                   if mutation in muts.get(gene, []))
        
        return count / len(self.strain_profiles)
    
    def find_co_occurring_mutations(self, gene: str, mutation: str) -> List[Dict]:
        """Find mutations that co-occur with target mutation"""
        co_occurring = defaultdict(int)
        
        # Find strains with target mutation
        target_strains = [strain_id for strain_id, muts in self.strain_profiles.items()
                         if mutation in muts.get(gene, [])]
        
        # Count co-occurring mutations
        for strain_id in target_strains:
            for g, muts in self.strain_profiles[strain_id].items():
                for m in muts:
                    if g != gene or m != mutation:
                        co_occurring[f"{g}:{m}"] += 1
        
        # Calculate frequencies
        results = []
        for mut_key, count in co_occurring.items():
            g, m = mut_key.split(':')
            results.append({
                'gene': g,
                'mutation': m,
                'frequency': count / len(target_strains),
                'count': count
            })
        
        return sorted(results, key=lambda x: x['frequency'], reverse=True)
    
    def get_strain_profile(self, strain_id: str) -> Dict:
        """Get complete mutation profile for strain"""
        return self.strain_profiles.get(strain_id, {})
    
    def compare_strains(self, strain1: str, strain2: str) -> Dict:
        """Compare mutation profiles between strains"""
        profile1 = self.strain_profiles.get(strain1, {})
        profile2 = self.strain_profiles.get(strain2, {})
        
        all_genes = set(profile1.keys()) | set(profile2.keys())
        
        comparison = {
            'shared_mutations': [],
            'unique_to_strain1': [],
            'unique_to_strain2': [],
            'similarity_score': 0.0
        }
        
        total_mutations = 0
        shared_count = 0
        
        for gene in all_genes:
            muts1 = set(profile1.get(gene, []))
            muts2 = set(profile2.get(gene, []))
            
            shared = muts1 & muts2
            unique1 = muts1 - muts2
            unique2 = muts2 - muts1
            
            comparison['shared_mutations'].extend([{'gene': gene, 'mutation': m} for m in shared])
            comparison['unique_to_strain1'].extend([{'gene': gene, 'mutation': m} for m in unique1])
            comparison['unique_to_strain2'].extend([{'gene': gene, 'mutation': m} for m in unique2])
            
            total_mutations += len(muts1 | muts2)
            shared_count += len(shared)
        
        if total_mutations > 0:
            comparison['similarity_score'] = shared_count / total_mutations
        
        return comparison
