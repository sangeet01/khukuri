"""Virtual lab orchestrator"""

import logging
from .base_agent import BaseAgent

logger = logging.getLogger('khukuri')


class VirtualLab:
    """Multi-agent virtual lab"""
    
    def __init__(self, openai_client=None):
        self.openai_client = openai_client
        self.agents = self._initialize_agents()
    
    def _initialize_agents(self):
        """Initialize agent team"""
        return {
            'chemist': BaseAgent('Chemist', 'molecular design and synthesis', self.openai_client),
            'biologist': BaseAgent('Biologist', 'target biology and mechanisms', self.openai_client),
            'pharmacologist': BaseAgent('Pharmacologist', 'ADMET and drug properties', self.openai_client)
        }
    
    def run_meeting(self, data, question):
        """Run multi-agent meeting"""
        results = {}
        for name, agent in self.agents.items():
            results[name] = agent.analyze(data, question)
        
        logger.info(f"Meeting completed with {len(results)} agents")
        return results
