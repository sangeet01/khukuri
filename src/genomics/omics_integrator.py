"""Integrate multi-omics data for AMR analysis"""

import logging
from typing import Dict, List, Optional
import numpy as np

logger = logging.getLogger('khukuri')


class OmicsIntegrator:
    """Integrate genomics, transcriptomics, and metabolomics data"""
    
    def __init__(self):
        self.genomic_data = {}
        self.transcriptomic_data = {}
        self.metabolomic_data = {}
    
    def add_genomic_data(self, sample_id: str, data: Dict):
        """Add genomic data (mutations, SNPs)"""
        self.genomic_data[sample_id] = data
        logger.info(f"Added genomic data for {sample_id}")
    
    def add_transcriptomic_data(self, sample_id: str, expression_data: Dict[str, float]):
        """Add gene expression data"""
        self.transcriptomic_data[sample_id] = expression_data
        logger.info(f"Added transcriptomic data for {sample_id}")
    
    def add_metabolomic_data(self, sample_id: str, metabolite_data: Dict[str, float]):
        """Add metabolite concentration data"""
        self.metabolomic_data[sample_id] = metabolite_data
        logger.info(f"Added metabolomic data for {sample_id}")
    
    def integrate_resistance_profile(self, sample_id: str) -> Dict:
        """Integrate multi-omics for resistance profile"""
        profile = {
            'sample_id': sample_id,
            'resistance_score': 0.0,
            'mechanisms': [],
            'confidence': 0.0
        }
        
        scores = []
        
        # Genomic contribution
        if sample_id in self.genomic_data:
            genomic_score = self._analyze_genomic_resistance(self.genomic_data[sample_id])
            scores.append(genomic_score)
            profile['mechanisms'].append('genetic_mutations')
        
        # Transcriptomic contribution
        if sample_id in self.transcriptomic_data:
            transcriptomic_score = self._analyze_expression_resistance(self.transcriptomic_data[sample_id])
            scores.append(transcriptomic_score)
            if transcriptomic_score > 0.5:
                profile['mechanisms'].append('efflux_upregulation')
        
        # Metabolomic contribution
        if sample_id in self.metabolomic_data:
            metabolomic_score = self._analyze_metabolic_resistance(self.metabolomic_data[sample_id])
            scores.append(metabolomic_score)
        
        if scores:
            profile['resistance_score'] = np.mean(scores)
            profile['confidence'] = len(scores) / 3.0  # Max confidence with all 3 omics
        
        return profile
    
    def _analyze_genomic_resistance(self, genomic_data: Dict) -> float:
        """Analyze genomic data for resistance markers"""
        resistance_genes = genomic_data.get('resistance_genes', [])
        mutations = genomic_data.get('mutations', {})
        
        score = 0.0
        if resistance_genes:
            score += min(len(resistance_genes) * 0.2, 0.6)
        
        if mutations:
            high_impact = sum(1 for m in mutations.values() if m.get('impact') == 'high')
            score += min(high_impact * 0.15, 0.4)
        
        return min(score, 1.0)
    
    def _analyze_expression_resistance(self, expression_data: Dict[str, float]) -> float:
        """Analyze expression data for resistance patterns"""
        # Efflux pump genes
        efflux_genes = ['mexA', 'mexB', 'acrA', 'acrB', 'tolC']
        efflux_expression = [expression_data.get(gene, 0.0) for gene in efflux_genes]
        
        # Target modification genes
        modification_genes = ['ermB', 'ermC', 'mphA']
        modification_expression = [expression_data.get(gene, 0.0) for gene in modification_genes]
        
        score = 0.0
        if efflux_expression:
            avg_efflux = np.mean(efflux_expression)
            if avg_efflux > 2.0:  # 2-fold upregulation
                score += 0.4
        
        if modification_expression:
            avg_mod = np.mean(modification_expression)
            if avg_mod > 2.0:
                score += 0.3
        
        return min(score, 1.0)
    
    def _analyze_metabolic_resistance(self, metabolite_data: Dict[str, float]) -> float:
        """Analyze metabolite data for resistance indicators"""
        # Simplified model - look for metabolic shifts
        score = 0.0
        
        # Check for altered cell wall precursors
        if 'UDP-NAG' in metabolite_data and metabolite_data['UDP-NAG'] > 1.5:
            score += 0.2
        
        # Check for stress metabolites
        if 'ppGpp' in metabolite_data and metabolite_data['ppGpp'] > 2.0:
            score += 0.3
        
        return min(score, 1.0)
    
    def find_correlations(self, gene: str, metabolite: str) -> Dict:
        """Find correlations between gene expression and metabolite levels"""
        expression_values = []
        metabolite_values = []
        
        for sample_id in self.transcriptomic_data.keys():
            if sample_id in self.metabolomic_data:
                expr = self.transcriptomic_data[sample_id].get(gene, 0.0)
                metab = self.metabolomic_data[sample_id].get(metabolite, 0.0)
                expression_values.append(expr)
                metabolite_values.append(metab)
        
        correlation = 0.0
        if len(expression_values) > 2:
            correlation = np.corrcoef(expression_values, metabolite_values)[0, 1]
        
        return {
            'gene': gene,
            'metabolite': metabolite,
            'correlation': float(correlation),
            'n_samples': len(expression_values)
        }
