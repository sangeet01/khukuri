"""Explainable AI - generates human-readable explanations for decisions"""

import logging
import json
from typing import Dict, List, Any, Optional
from .base_agent import BaseAgent

logger = logging.getLogger('khukuri')


class ExplainableAI(BaseAgent):
    """Generates clear explanations for AI decisions"""
    
    def __init__(self, openai_client=None):
        super().__init__(
            role="Explainability Specialist",
            expertise="translating complex AI decisions into clear explanations",
            openai_client=openai_client
        )
        self.explanation_history = []
    
    def explain_ranking(self, compounds: List[Dict[str, Any]],
                       ranking_criteria: Dict[str, float]) -> Dict[str, Any]:
        """
        Explain why compounds were ranked in a certain order
        
        Args:
            compounds: List of compounds with scores
            ranking_criteria: Weights used for ranking
        
        Returns:
            Human-readable explanation
        """
        logger.info(f"Explaining ranking of {len(compounds)} compounds")
        
        if self.openai_client:
            return self._ai_explain_ranking(compounds, ranking_criteria)
        return self._rule_based_ranking_explanation(compounds, ranking_criteria)
    
    def _ai_explain_ranking(self, compounds: List[Dict[str, Any]],
                           criteria: Dict[str, float]) -> Dict[str, Any]:
        """AI-powered ranking explanation"""
        try:
            # Take top 5 for explanation
            top_compounds = compounds[:5]
            
            prompt = f"""You are explaining why drug candidates were ranked in a specific order.

RANKING CRITERIA (weights):
{json.dumps(criteria, indent=2)}

TOP COMPOUNDS:
{json.dumps(top_compounds, indent=2, default=str)}

Explain in simple terms:
1. Why #1 compound is ranked first
2. Key differences between top compounds
3. What makes a compound score well
4. Trade-offs in the ranking

Return JSON:
{{
  "summary": "One sentence summary of ranking logic",
  "top_compound_explanation": {{
    "compound_name": "name",
    "why_ranked_first": "Clear explanation",
    "key_strengths": ["strength1", "strength2"],
    "potential_concerns": ["concern1", "concern2"]
  }},
  "ranking_factors": [
    {{
      "factor": "binding_affinity",
      "importance": "high/medium/low",
      "explanation": "Why this matters"
    }}
  ],
  "trade_offs": ["Trade-off 1", "Trade-off 2"],
  "confidence": "high/medium/low",
  "layman_summary": "Explanation for non-experts"
}}"""
            
            response = self.openai_client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {"role": "system", "content": "You are an Explainability Specialist. Make complex science understandable. Return only valid JSON."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.3,
                max_tokens=1500
            )
            
            content = response.choices[0].message.content.strip()
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0].strip()
            elif "```" in content:
                content = content.split("```")[1].split("```")[0].strip()
            
            explanation = json.loads(content)
            self.explanation_history.append({
                'type': 'ranking',
                'explanation': explanation
            })
            
            return explanation
            
        except Exception as e:
            logger.warning(f"AI explanation failed: {e}, using rule-based")
            return self._rule_based_ranking_explanation(compounds, criteria)
    
    def _rule_based_ranking_explanation(self, compounds: List[Dict[str, Any]],
                                       criteria: Dict[str, float]) -> Dict[str, Any]:
        """Rule-based ranking explanation"""
        
        if not compounds:
            return {'summary': 'No compounds to rank', 'confidence': 'low'}
        
        top = compounds[0]
        top_name = top.get('name', 'Compound #1')
        
        # Find strongest factor
        scores = top.get('scores', {})
        if scores:
            best_factor = max(scores.items(), key=lambda x: x[1])
            best_factor_name, best_score = best_factor
        else:
            best_factor_name, best_score = 'overall_score', top.get('total_score', 0)
        
        return {
            'summary': f"{top_name} ranked first due to superior {best_factor_name}",
            'top_compound_explanation': {
                'compound_name': top_name,
                'why_ranked_first': f"Highest {best_factor_name} score ({best_score:.2f})",
                'key_strengths': [f"Strong {best_factor_name}"],
                'potential_concerns': ['Requires experimental validation']
            },
            'ranking_factors': [
                {
                    'factor': factor,
                    'importance': 'high' if weight > 0.3 else 'medium' if weight > 0.1 else 'low',
                    'explanation': f"Weighted at {weight*100:.0f}% in ranking"
                }
                for factor, weight in criteria.items()
            ],
            'trade_offs': ['Balancing multiple properties'],
            'confidence': 'medium',
            'layman_summary': f"{top_name} is predicted to work best based on computational analysis"
        }
    
    def explain_decision(self, decision: str, context: Dict[str, Any],
                        reasoning: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Explain any AI decision
        
        Args:
            decision: The decision made
            context: Context in which decision was made
            reasoning: Optional reasoning data
        
        Returns:
            Clear explanation
        """
        logger.info(f"Explaining decision: {decision}")
        
        if self.openai_client:
            return self._ai_explain_decision(decision, context, reasoning)
        return self._simple_decision_explanation(decision, context)
    
    def _ai_explain_decision(self, decision: str, context: Dict[str, Any],
                            reasoning: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """AI-powered decision explanation"""
        try:
            prompt = f"""You are explaining an AI decision in drug discovery.

DECISION: {decision}

CONTEXT:
{json.dumps(context, indent=2, default=str)}

REASONING DATA:
{json.dumps(reasoning, indent=2, default=str) if reasoning else "Not provided"}

Explain:
1. What decision was made
2. Why this decision makes sense
3. What data supported it
4. What alternatives were considered
5. Potential risks/limitations

Return JSON:
{{
  "decision_summary": "Clear statement of what was decided",
  "rationale": "Why this decision was made",
  "supporting_evidence": ["Evidence 1", "Evidence 2"],
  "alternatives_considered": ["Alternative 1", "Alternative 2"],
  "risks": ["Risk 1", "Risk 2"],
  "confidence_level": "high/medium/low",
  "next_steps": ["What to do next"],
  "for_experts": "Technical explanation",
  "for_non_experts": "Simple explanation"
}}"""
            
            response = self.openai_client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {"role": "system", "content": "You are an Explainability Specialist. Return only valid JSON."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.3,
                max_tokens=1500
            )
            
            content = response.choices[0].message.content.strip()
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0].strip()
            elif "```" in content:
                content = content.split("```")[1].split("```")[0].strip()
            
            explanation = json.loads(content)
            self.explanation_history.append({
                'type': 'decision',
                'decision': decision,
                'explanation': explanation
            })
            
            return explanation
            
        except Exception as e:
            logger.warning(f"AI explanation failed: {e}, using simple explanation")
            return self._simple_decision_explanation(decision, context)
    
    def _simple_decision_explanation(self, decision: str,
                                    context: Dict[str, Any]) -> Dict[str, Any]:
        """Simple decision explanation"""
        return {
            'decision_summary': decision,
            'rationale': 'Based on available data and computational predictions',
            'supporting_evidence': ['Computational analysis', 'Historical data'],
            'alternatives_considered': ['Other approaches were evaluated'],
            'risks': ['Predictions require experimental validation'],
            'confidence_level': 'medium',
            'next_steps': ['Validate computationally', 'Test experimentally'],
            'for_experts': f"Decision: {decision}. Context: {context}",
            'for_non_experts': f"We decided to {decision} based on computer predictions"
        }
    
    def generate_report(self, results: Dict[str, Any],
                       report_type: str = 'comprehensive') -> str:
        """
        Generate human-readable report
        
        Args:
            results: Discovery results
            report_type: 'comprehensive', 'executive', or 'technical'
        
        Returns:
            Formatted report
        """
        logger.info(f"Generating {report_type} report")
        
        if report_type == 'executive':
            return self._executive_summary(results)
        elif report_type == 'technical':
            return self._technical_report(results)
        else:
            return self._comprehensive_report(results)
    
    def _executive_summary(self, results: Dict[str, Any]) -> str:
        """Executive summary for non-technical stakeholders"""
        
        n_compounds = len(results.get('candidates', []))
        top_compound = results.get('candidates', [{}])[0] if results.get('candidates') else {}
        
        return f"""
EXECUTIVE SUMMARY
=================

PROJECT: {results.get('project_name', 'Drug Discovery')}
DATE: {results.get('date', 'N/A')}

KEY FINDINGS:
-------------
• Identified {n_compounds} promising drug candidates
• Top candidate: {top_compound.get('name', 'N/A')}
• Predicted success rate: {results.get('success_probability', 'N/A')}

RECOMMENDATION:
---------------
{results.get('recommendation', 'Proceed with experimental validation of top candidates')}

NEXT STEPS:
-----------
1. Synthesize top 3-5 candidates
2. Conduct in vitro testing
3. Evaluate safety profile

TIMELINE: {results.get('estimated_timeline', '3-6 months for initial validation')}
BUDGET: {results.get('estimated_cost', 'TBD based on synthesis complexity')}
"""
    
    def _technical_report(self, results: Dict[str, Any]) -> str:
        """Technical report for scientists"""
        
        return f"""
TECHNICAL REPORT
================

METHODOLOGY:
------------
{results.get('methodology', 'Computational drug discovery pipeline')}

RESULTS:
--------
Candidates Generated: {len(results.get('candidates', []))}
Docking Scores: {results.get('docking_summary', 'N/A')}
ADMET Pass Rate: {results.get('admet_pass_rate', 'N/A')}

TOP CANDIDATES:
---------------
{chr(10).join([f"  {i+1}. {c.get('name', 'N/A')}: Score {c.get('total_score', 'N/A')}" for i, c in enumerate(results.get('candidates', [])[:5])])}

COMPUTATIONAL DETAILS:
----------------------
{json.dumps(results.get('computational_details', {}), indent=2)}

VALIDATION PLAN:
----------------
{results.get('validation_plan', 'Standard in vitro/in vivo testing')}
"""
    
    def _comprehensive_report(self, results: Dict[str, Any]) -> str:
        """Comprehensive report with all details"""
        
        exec_summary = self._executive_summary(results)
        tech_report = self._technical_report(results)
        
        return f"{exec_summary}\n\n{tech_report}\n\nFULL RESULTS:\n{json.dumps(results, indent=2, default=str)}"
    
    def get_explanation_history(self) -> List[Dict[str, Any]]:
        """Get history of all explanations"""
        return self.explanation_history
