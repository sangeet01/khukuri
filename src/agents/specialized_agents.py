"""Specialized AI agents for drug discovery"""

from .base_agent import BaseAgent


class ChemistAgent(BaseAgent):
    """Agent specialized in medicinal chemistry"""
    
    def __init__(self, api_key=None):
        super().__init__(
            name="Chemist",
            role="medicinal_chemistry",
            api_key=api_key
        )
    
    def analyze(self, data):
        """Analyze from chemistry perspective"""
        prompt = f"""As a medicinal chemist, analyze this molecule:
        SMILES: {data.get('smiles', 'N/A')}
        Properties: {data.get('properties', {})}
        
        Provide brief assessment of drug-likeness and synthesis."""
        
        return self._call_llm(prompt, fallback={
            "assessment": "Drug-like properties acceptable",
            "concerns": ["Verify synthetic accessibility"],
            "recommendations": ["Consider bioisosteric replacements"]
        })


class BiologistAgent(BaseAgent):
    """Agent specialized in biology and pharmacology"""
    
    def __init__(self, api_key=None):
        super().__init__(
            name="Biologist",
            role="pharmacology",
            api_key=api_key
        )
    
    def analyze(self, data):
        """Analyze from biology perspective"""
        prompt = f"""As a pharmacologist, analyze this drug candidate:
        Target: {data.get('target', 'N/A')}
        Binding: {data.get('binding_affinity', 'N/A')}
        ADMET: {data.get('admet', {})}
        
        Assess pharmacological viability."""
        
        return self._call_llm(prompt, fallback={
            "assessment": "Pharmacologically viable",
            "concerns": ["Verify target selectivity"],
            "recommendations": ["Test in vitro activity"]
        })


class ToxicologistAgent(BaseAgent):
    """Agent specialized in toxicology"""
    
    def __init__(self, api_key=None):
        super().__init__(
            name="Toxicologist",
            role="toxicology",
            api_key=api_key
        )
    
    def analyze(self, data):
        """Analyze toxicity risks"""
        prompt = f"""As a toxicologist, assess safety:
        Molecule: {data.get('smiles', 'N/A')}
        Toxicity alerts: {data.get('toxicity', {})}
        
        Identify safety concerns."""
        
        return self._call_llm(prompt, fallback={
            "assessment": "No major toxicity alerts",
            "concerns": ["Monitor hepatotoxicity"],
            "recommendations": ["Conduct safety studies"]
        })
