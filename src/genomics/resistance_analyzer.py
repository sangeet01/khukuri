"""Genomic resistance pattern analysis"""

import logging
from typing import Dict, List, Optional
import numpy as np

logger = logging.getLogger('khukuri')


class ResistanceGenomicsAnalyzer:
    """Analyze genomic patterns of resistance"""
    
    def __init__(self):
        self.known_mutations = self._initialize_mutations()
    
    def _initialize_mutations(self) -> Dict:
        """Known resistance mutations"""
        return {
            'rpoB': {
                'S531L': {'drug': 'rifampicin', 'frequency': 0.65, 'fitness_cost': 0.15},
                'H526Y': {'drug': 'rifampicin', 'frequency': 0.15, 'fitness_cost': 0.10},
                'D516V': {'drug': 'rifampicin', 'frequency': 0.10, 'fitness_cost': 0.12}
            },
            'gyrA': {
                'S83L': {'drug': 'fluoroquinolone', 'frequency': 0.45, 'fitness_cost': 0.08},
                'D87N': {'drug': 'fluoroquinolone', 'frequency': 0.30, 'fitness_cost': 0.10}
            },
            'katG': {
                'S315T': {'drug': 'isoniazid', 'frequency': 0.70, 'fitness_cost': 0.05}
            }
        }
    
    def analyze_mutation_profile(self, gene: str, mutations: List[str]) -> Dict:
        """Analyze mutation profile for resistance"""
        results = {
            'gene': gene,
            'mutations': [],
            'predicted_resistance': [],
            'fitness_cost': 0.0,
            'confidence': 0.0
        }
        
        if gene not in self.known_mutations:
            results['confidence'] = 0.1
            return results
        
        gene_mutations = self.known_mutations[gene]
        total_frequency = 0.0
        total_cost = 0.0
        
        for mutation in mutations:
            if mutation in gene_mutations:
                mut_info = gene_mutations[mutation]
                results['mutations'].append({
                    'mutation': mutation,
                    'drug': mut_info['drug'],
                    'frequency': mut_info['frequency']
                })
                results['predicted_resistance'].append(mut_info['drug'])
                total_frequency += mut_info['frequency']
                total_cost += mut_info['fitness_cost']
        
        if results['mutations']:
            results['fitness_cost'] = total_cost / len(results['mutations'])
            results['confidence'] = min(total_frequency / len(mutations), 1.0)
        
        return results
    
    def predict_resistance_evolution(self, current_mutations: Dict, drug_pressure: str) -> Dict:
        """Predict likely resistance evolution under drug pressure"""
        predictions = {
            'likely_mutations': [],
            'timeline_generations': 0,
            'probability': 0.0
        }
        
        # Find mutations associated with drug
        for gene, mutations in self.known_mutations.items():
            for mutation, info in mutations.items():
                if info['drug'] == drug_pressure and mutation not in current_mutations.get(gene, []):
                    predictions['likely_mutations'].append({
                        'gene': gene,
                        'mutation': mutation,
                        'frequency': info['frequency'],
                        'fitness_cost': info['fitness_cost']
                    })
        
        if predictions['likely_mutations']:
            # Estimate timeline based on frequency and fitness cost
            avg_freq = np.mean([m['frequency'] for m in predictions['likely_mutations']])
            avg_cost = np.mean([m['fitness_cost'] for m in predictions['likely_mutations']])
            
            predictions['probability'] = avg_freq * (1 - avg_cost)
            predictions['timeline_generations'] = int(100 / (avg_freq * 100 + 1))
        
        return predictions
    
    def identify_compensatory_mutations(self, primary_mutation: str, gene: str) -> List[Dict]:
        """Identify potential compensatory mutations"""
        # Simplified model - in production, use ML or database
        compensatory = []
        
        if gene == 'rpoB' and primary_mutation in ['S531L', 'H526Y']:
            compensatory.append({
                'gene': 'rpoC',
                'mutation': 'I491V',
                'effect': 'Restores fitness',
                'evidence': 'literature'
            })
        
        return compensatory
