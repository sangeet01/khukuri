"""Failure Analysis Agent - learns from experimental failures"""

import logging
import json
from typing import Dict, List, Any, Optional
from datetime import datetime
from .base_agent import BaseAgent

logger = logging.getLogger('khukuri')


class FailureAnalyzer(BaseAgent):
    """Analyzes why compounds fail and learns from mistakes"""
    
    def __init__(self, openai_client=None):
        super().__init__(
            role="Failure Analysis Specialist",
            expertise="analyzing experimental failures and extracting learnings",
            openai_client=openai_client
        )
        self.failure_database = []
        self.learned_patterns = {
            'avoid_scaffolds': [],
            'avoid_substituents': [],
            'property_thresholds': {},
            'failure_modes': {}
        }
    
    def analyze_failure(self, compound_data: Dict[str, Any], 
                       experimental_result: Dict[str, Any]) -> Dict[str, Any]:
        """
        Analyze why a compound failed
        
        Args:
            compound_data: Compound structure and predicted properties
            experimental_result: Actual experimental outcomes
        
        Returns:
            Analysis with root causes and recommendations
        """
        logger.info(f"Analyzing failure for compound: {compound_data.get('name', 'unknown')}")
        
        failure_record = {
            'timestamp': datetime.now().isoformat(),
            'compound': compound_data,
            'experimental': experimental_result,
            'analysis': {}
        }
        
        if self.openai_client:
            analysis = self._ai_analyze_failure(compound_data, experimental_result)
        else:
            analysis = self._rule_based_analysis(compound_data, experimental_result)
        
        failure_record['analysis'] = analysis
        self.failure_database.append(failure_record)
        
        # Update learned patterns
        self._update_learned_patterns(analysis)
        
        return analysis
    
    def _ai_analyze_failure(self, compound_data: Dict[str, Any],
                           experimental_result: Dict[str, Any]) -> Dict[str, Any]:
        """AI-powered failure analysis"""
        try:
            prompt = f"""You are a Failure Analysis Specialist analyzing why a drug candidate failed.

COMPOUND DATA:
{json.dumps(compound_data, indent=2, default=str)}

EXPERIMENTAL RESULT:
{json.dumps(experimental_result, indent=2, default=str)}

PREVIOUS FAILURES:
{json.dumps(self.failure_database[-5:], indent=2, default=str) if self.failure_database else "None"}

Analyze:
1. Root cause of failure
2. Which structural features contributed
3. Which property predictions were wrong
4. Patterns across multiple failures
5. Actionable recommendations

Return JSON:
{{
  "root_cause": "Primary reason for failure",
  "contributing_factors": ["Factor 1", "Factor 2"],
  "structural_issues": {{
    "problematic_scaffold": "scaffold name or null",
    "problematic_groups": ["group1", "group2"]
  }},
  "prediction_errors": {{
    "property": "predicted vs actual"
  }},
  "patterns_detected": ["Pattern 1", "Pattern 2"],
  "recommendations": {{
    "avoid": ["What to avoid in future"],
    "prioritize": ["What to prioritize"],
    "modify": ["How to modify approach"]
  }},
  "confidence": "high/medium/low"
}}"""
            
            response = self.openai_client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {"role": "system", "content": "You are a Failure Analysis Specialist. Return only valid JSON."},
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
            
            return json.loads(content)
            
        except Exception as e:
            logger.warning(f"AI failure analysis failed: {e}, using rule-based")
            return self._rule_based_analysis(compound_data, experimental_result)
    
    def _rule_based_analysis(self, compound_data: Dict[str, Any],
                            experimental_result: Dict[str, Any]) -> Dict[str, Any]:
        """Rule-based failure analysis"""
        
        predicted = compound_data.get('predicted_properties', {})
        actual = experimental_result.get('measured_properties', {})
        
        root_cause = "Unknown"
        contributing_factors = []
        
        # Analyze prediction errors
        if 'activity' in predicted and 'activity' in actual:
            pred_activity = predicted['activity']
            actual_activity = actual['activity']
            if abs(pred_activity - actual_activity) > 5:
                root_cause = "Poor activity prediction"
                contributing_factors.append("Docking score unreliable")
        
        # Check property violations
        if actual.get('solubility', 1) < 0.01:
            contributing_factors.append("Poor solubility")
        
        if actual.get('toxicity', 0) > 0.5:
            contributing_factors.append("High toxicity")
        
        if actual.get('stability', 1) < 0.3:
            contributing_factors.append("Chemical instability")
        
        return {
            'root_cause': root_cause,
            'contributing_factors': contributing_factors,
            'structural_issues': {
                'problematic_scaffold': None,
                'problematic_groups': []
            },
            'prediction_errors': {
                'activity': f"Predicted {predicted.get('activity', 'N/A')}, Got {actual.get('activity', 'N/A')}"
            },
            'patterns_detected': [],
            'recommendations': {
                'avoid': ['Similar scaffolds with poor solubility'],
                'prioritize': ['Compounds with better predicted ADMET'],
                'modify': ['Add solubilizing groups', 'Reduce molecular weight']
            },
            'confidence': 'medium'
        }
    
    def _update_learned_patterns(self, analysis: Dict[str, Any]):
        """Update learned patterns from analysis"""
        
        # Track problematic scaffolds
        scaffold = analysis.get('structural_issues', {}).get('problematic_scaffold')
        if scaffold and scaffold not in self.learned_patterns['avoid_scaffolds']:
            self.learned_patterns['avoid_scaffolds'].append(scaffold)
        
        # Track problematic groups
        groups = analysis.get('structural_issues', {}).get('problematic_groups', [])
        for group in groups:
            if group not in self.learned_patterns['avoid_substituents']:
                self.learned_patterns['avoid_substituents'].append(group)
        
        # Track failure modes
        root_cause = analysis.get('root_cause', 'unknown')
        if root_cause not in self.learned_patterns['failure_modes']:
            self.learned_patterns['failure_modes'][root_cause] = 0
        self.learned_patterns['failure_modes'][root_cause] += 1
    
    def get_design_constraints(self) -> Dict[str, Any]:
        """Get design constraints based on learned failures"""
        return {
            'avoid_scaffolds': self.learned_patterns['avoid_scaffolds'],
            'avoid_substituents': self.learned_patterns['avoid_substituents'],
            'common_failure_modes': sorted(
                self.learned_patterns['failure_modes'].items(),
                key=lambda x: x[1],
                reverse=True
            )[:5],
            'total_failures_analyzed': len(self.failure_database)
        }
    
    def generate_report(self) -> str:
        """Generate failure analysis report"""
        
        if not self.failure_database:
            return "No failures analyzed yet."
        
        report = f"""
FAILURE ANALYSIS REPORT
=======================
Total Failures Analyzed: {len(self.failure_database)}

LEARNED PATTERNS:
-----------------
Scaffolds to Avoid: {len(self.learned_patterns['avoid_scaffolds'])}
{chr(10).join(['  - ' + s for s in self.learned_patterns['avoid_scaffolds'][:5]])}

Substituents to Avoid: {len(self.learned_patterns['avoid_substituents'])}
{chr(10).join(['  - ' + s for s in self.learned_patterns['avoid_substituents'][:5]])}

COMMON FAILURE MODES:
---------------------
{chr(10).join([f'  {i+1}. {mode}: {count} occurrences' for i, (mode, count) in enumerate(sorted(self.learned_patterns['failure_modes'].items(), key=lambda x: x[1], reverse=True)[:5])])}

RECENT FAILURES:
----------------
{chr(10).join([f"  - {f['compound'].get('name', 'unknown')}: {f['analysis'].get('root_cause', 'unknown')}" for f in self.failure_database[-5:]])}

RECOMMENDATIONS:
----------------
Based on {len(self.failure_database)} failures:
  1. Screen compounds against learned constraints before synthesis
  2. Focus on scaffolds with no failure history
  3. Prioritize compounds with diverse mechanisms
  4. Validate predictions with multiple methods
"""
        return report
    
    def should_synthesize(self, compound_data: Dict[str, Any]) -> Dict[str, bool]:
        """Check if compound should be synthesized based on learned failures"""
        
        warnings = []
        
        # Check scaffold
        scaffold = compound_data.get('scaffold', '')
        if scaffold in self.learned_patterns['avoid_scaffolds']:
            warnings.append(f"Scaffold '{scaffold}' has failed before")
        
        # Check substituents
        substituents = compound_data.get('substituents', [])
        for sub in substituents:
            if sub in self.learned_patterns['avoid_substituents']:
                warnings.append(f"Substituent '{sub}' has caused failures")
        
        # Check properties against learned thresholds
        props = compound_data.get('predicted_properties', {})
        if props.get('solubility', 1) < 0.01:
            warnings.append("Very low predicted solubility")
        
        should_proceed = len(warnings) == 0
        
        return {
            'should_synthesize': should_proceed,
            'warnings': warnings,
            'confidence': 'high' if len(self.failure_database) > 10 else 'low'
        }
