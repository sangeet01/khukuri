"""Specialized agents for AMR discovery"""

import logging
from .base_agent import BaseAgent

logger = logging.getLogger('khukuri')


class MicrobiologyAgent(BaseAgent):
    """Agent specialized in microbiology and pathogen analysis"""
    
    def __init__(self, openai_client=None):
        super().__init__(
            role="Microbiologist",
            expertise="pathogen biology, AMR mechanisms, MIC analysis",
            openai_client=openai_client
        )
    
    def analyze(self, data, question):
        """Analyze from microbiology perspective"""
        context = f"""As a microbiologist specializing in AMR, analyze:
        
Pathogen: {data.get('pathogen', 'N/A')}
Strain: {data.get('strain_id', 'N/A')}
Resistance Profile: {data.get('resistance_profile', 'N/A')}
MIC Data: {data.get('mic_data', 'N/A')}

Question: {question}

Provide assessment of resistance mechanisms and compound efficacy."""
        
        if self.openai_client:
            return self._ai_analyze(data, context)
        
        return {
            'role': self.role,
            'assessment': 'Resistance mechanisms identified',
            'concerns': ['Monitor efflux pump activity', 'Check for target mutations'],
            'recommendations': ['Test against resistant strains', 'Evaluate cross-resistance']
        }


class GenomicsAgent(BaseAgent):
    """Agent specialized in genomics and resistance patterns"""
    
    def __init__(self, openai_client=None):
        super().__init__(
            role="Genomics Specialist",
            expertise="resistance genomics, mutation analysis, evolutionary patterns",
            openai_client=openai_client
        )
    
    def analyze(self, data, question):
        """Analyze from genomics perspective"""
        context = f"""As a genomics specialist, analyze:
        
Resistance Genes: {data.get('resistance_genes', [])}
Mutations: {data.get('mutations', {})}
Evolutionary Pressure: {data.get('drug_pressure', 'N/A')}

Question: {question}

Assess genetic basis of resistance and evolution potential."""
        
        if self.openai_client:
            return self._ai_analyze(data, context)
        
        return {
            'role': self.role,
            'assessment': 'Genetic resistance markers identified',
            'concerns': ['High mutation frequency', 'Compensatory mutations possible'],
            'recommendations': ['Monitor resistance evolution', 'Consider multi-target approach']
        }


class CheminformaticsAgent(BaseAgent):
    """Agent specialized in cheminformatics and molecular design"""
    
    def __init__(self, openai_client=None):
        super().__init__(
            role="Cheminformatics Specialist",
            expertise="molecular design, QSAR, docking, ADMET",
            openai_client=openai_client
        )
    
    def analyze(self, data, question):
        """Analyze from cheminformatics perspective"""
        context = f"""As a cheminformatics specialist, analyze:
        
Compound: {data.get('smiles', 'N/A')}
Target: {data.get('target', 'N/A')}
Binding Affinity: {data.get('binding_affinity', 'N/A')}
ADMET: {data.get('admet', {})}

Question: {question}

Evaluate molecular properties and optimization potential."""
        
        if self.openai_client:
            return self._ai_analyze(data, context)
        
        return {
            'role': self.role,
            'assessment': 'Compound shows drug-like properties',
            'concerns': ['Optimize binding affinity', 'Improve metabolic stability'],
            'recommendations': ['Explore bioisosteres', 'Enhance target selectivity']
        }


class ResistanceCriticAgent(BaseAgent):
    """Agent specialized in resistance prediction and criticism"""
    
    def __init__(self, openai_client=None):
        super().__init__(
            role="Resistance Critic",
            expertise="resistance prediction, multi-target strategies, evolutionary analysis",
            openai_client=openai_client
        )
    
    def analyze(self, data, question):
        """Critically analyze resistance potential"""
        context = f"""As a resistance prediction expert, critically evaluate:
        
Compound: {data.get('compound_id', 'N/A')}
Target: {data.get('target', 'N/A')}
Mechanism: {data.get('mechanism', 'N/A')}
Known Resistance: {data.get('known_resistance', [])}

Question: {question}

Identify resistance risks and suggest mitigation strategies."""
        
        if self.openai_client:
            return self._ai_analyze(data, context)
        
        return {
            'role': self.role,
            'assessment': 'Moderate resistance risk',
            'concerns': ['Single target vulnerability', 'Known resistance mechanisms exist'],
            'recommendations': ['Develop multi-target strategy', 'Design resistance-refractory compounds']
        }


class LiteratureAgent(BaseAgent):
    """Agent specialized in literature synthesis and evidence gathering"""
    
    def __init__(self, openai_client=None):
        super().__init__(
            role="Literature Analyst",
            expertise="scientific literature, evidence synthesis, hypothesis generation",
            openai_client=openai_client
        )
    
    def analyze(self, data, question):
        """Synthesize literature evidence"""
        context = f"""As a literature analyst, synthesize evidence for:
        
Topic: {data.get('topic', 'N/A')}
Target: {data.get('target', 'N/A')}
Pathogen: {data.get('pathogen', 'N/A')}

Question: {question}

Provide evidence-based assessment and identify knowledge gaps."""
        
        if self.openai_client:
            return self._ai_analyze(data, context)
        
        return {
            'role': self.role,
            'assessment': 'Literature supports target validity',
            'evidence': ['Multiple studies confirm target essentiality', 'Clinical data available'],
            'recommendations': ['Review recent resistance reports', 'Identify novel mechanisms']
        }
