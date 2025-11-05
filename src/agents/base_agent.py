"""Base agent class"""

import logging
import json

logger = logging.getLogger('khukuri')


class BaseAgent:
    """Base class for AI agents"""
    
    def __init__(self, role, expertise, openai_client=None):
        self.role = role
        self.expertise = expertise
        self.openai_client = openai_client
        self.memory = []
    
    def analyze(self, data, question):
        """Analyze data and answer question"""
        if self.openai_client:
            return self._ai_analyze(data, question)
        return self._fallback_analyze(data, question)
    
    def _ai_analyze(self, data, question):
        """AI-powered analysis"""
        try:
            response = self.openai_client.chat.completions.create(
                model="gpt-3.5-turbo",
                messages=[
                    {"role": "system", "content": f"You are a {self.role} with expertise in {self.expertise}."},
                    {"role": "user", "content": f"{question}\n\nData: {json.dumps(data, default=str)}"}
                ],
                temperature=0.7,
                max_tokens=500
            )
            return json.loads(response.choices[0].message.content)
        except:
            return self._fallback_analyze(data, question)
    
    def _fallback_analyze(self, data, question):
        """Fallback analysis without AI"""
        return {
            'role': self.role,
            'analysis': f"Analysis of {question}",
            'recommendations': ['Proceed with validation', 'Monitor results'],
            'confidence': 'medium'
        }
