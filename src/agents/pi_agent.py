"""Lead Scientist Agent - orchestrates drug discovery workflow"""

import logging
import json
from typing import Dict, List, Optional, Any
from .base_agent import BaseAgent

logger = logging.getLogger('khukuri')


class PIAgent(BaseAgent):
    """Lead Scientist agent that orchestrates the discovery process"""
    
    def __init__(self, openai_client=None):
        super().__init__(
            role="Lead Scientist",
            expertise="orchestrating AI-driven antimicrobial drug discovery",
            openai_client=openai_client
        )
        self.scientist_agents = []
        self.meeting_history = []
        self.project_context = {}
    
    def create_team(self, project_description: str) -> List[Dict[str, str]]:
        """Create scientist team based on project needs"""
        logger.info(f"PI creating team for: {project_description}")
        
        if self.openai_client:
            return self._ai_create_team(project_description)
        return self._default_team()
    
    def _ai_create_team(self, project_description: str) -> List[Dict[str, str]]:
        """AI-powered team creation"""
        try:
            prompt = f"""You are a Lead Scientist creating a team for antimicrobial drug discovery.

Project: {project_description}

Create a team of 3 scientist agents. For each agent, specify:
- title: Agent's role name
- expertise: Specific scientific expertise
- goal: What they aim to achieve
- role: Their responsibilities in the project

Return ONLY a JSON array of 3 agents. Example format:
[
  {{
    "title": "Computational Biologist",
    "expertise": "protein structure analysis and molecular docking",
    "goal": "identify optimal binding sites and predict binding affinity",
    "role": "analyze target proteins and guide docking simulations"
  }}
]"""
            
            response = self.openai_client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {"role": "system", "content": "You are a Lead Scientist in antimicrobial drug discovery. Return only valid JSON."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.7,
                max_tokens=800
            )
            
            content = response.choices[0].message.content.strip()
            # Extract JSON if wrapped in markdown
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0].strip()
            elif "```" in content:
                content = content.split("```")[1].split("```")[0].strip()
            
            team = json.loads(content)
            self.scientist_agents = team
            logger.info(f"PI created team: {[a['title'] for a in team]}")
            return team
            
        except Exception as e:
            logger.warning(f"AI team creation failed: {e}, using default team")
            return self._default_team()
    
    def _default_team(self) -> List[Dict[str, str]]:
        """Default team composition"""
        team = [
            {
                "title": "Computational Biologist",
                "expertise": "protein structure analysis and molecular docking",
                "goal": "identify optimal binding sites and predict binding affinity",
                "role": "analyze target proteins and guide docking simulations"
            },
            {
                "title": "Medicinal Chemist",
                "expertise": "small molecule design and ADMET optimization",
                "goal": "design drug-like molecules with favorable properties",
                "role": "generate and optimize molecular candidates"
            },
            {
                "title": "Resistance Specialist",
                "expertise": "antimicrobial resistance mechanisms and evolution",
                "goal": "predict and prevent resistance development",
                "role": "analyze resistance pathways and guide multi-target strategies"
            }
        ]
        self.scientist_agents = team
        return team
    
    def convene_team_meeting(self, agenda: str, context: Dict[str, Any], 
                            agenda_questions: Optional[List[str]] = None,
                            rounds: int = 3) -> Dict[str, Any]:
        """Run collaborative team meeting with all agents"""
        logger.info(f"PI convening team meeting: {agenda}")
        
        meeting_record = {
            'type': 'team_meeting',
            'agenda': agenda,
            'agenda_questions': agenda_questions or [],
            'rounds': rounds,
            'discussions': [],
            'summary': {},
            'recommendations': []
        }
        
        if self.openai_client:
            result = self._ai_team_meeting(agenda, context, agenda_questions, rounds)
        else:
            result = self._fallback_team_meeting(agenda, context)
        
        meeting_record.update(result)
        self.meeting_history.append(meeting_record)
        return result
    
    def _ai_team_meeting(self, agenda: str, context: Dict[str, Any],
                        agenda_questions: Optional[List[str]], rounds: int) -> Dict[str, Any]:
        """AI-powered team meeting"""
        try:
            # Build meeting context
            team_members = ", ".join([a['title'] for a in self.scientist_agents])
            questions_text = "\n".join([f"{i+1}. {q}" for i, q in enumerate(agenda_questions or [])])
            
            prompt = f"""You are the Lead Scientist leading a team meeting.

TEAM MEMBERS: {team_members}

AGENDA: {agenda}

AGENDA QUESTIONS:
{questions_text if questions_text else "None"}

CONTEXT:
{json.dumps(context, indent=2, default=str)}

PREVIOUS MEETINGS:
{json.dumps([m['summary'] for m in self.meeting_history[-3:]], indent=2) if self.meeting_history else "None"}

Conduct a {rounds}-round discussion where:
1. You provide initial thoughts and questions
2. Each team member provides their perspective
3. You synthesize and make decisions

Return a JSON object with:
{{
  "pi_opening": "Your initial thoughts and questions",
  "team_input": {{
    "agent_name": "Their perspective and recommendations"
  }},
  "synthesis": "Your synthesis of the discussion",
  "recommendations": ["Specific actionable recommendations"],
  "answers": {{"question": "answer"}} (if agenda questions provided),
  "next_steps": ["Next actions to take"]
}}"""
            
            response = self.openai_client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {"role": "system", "content": "You are a Lead Scientist in antimicrobial drug discovery. Return only valid JSON."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.7,
                max_tokens=2000
            )
            
            content = response.choices[0].message.content.strip()
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0].strip()
            elif "```" in content:
                content = content.split("```")[1].split("```")[0].strip()
            
            result = json.loads(content)
            logger.info("Team meeting completed successfully")
            return result
            
        except Exception as e:
            logger.warning(f"AI team meeting failed: {e}, using fallback")
            return self._fallback_team_meeting(agenda, context)
    
    def _fallback_team_meeting(self, agenda: str, context: Dict[str, Any]) -> Dict[str, Any]:
        """Fallback team meeting without AI"""
        return {
            'pi_opening': f"Let's discuss: {agenda}",
            'team_input': {
                agent['title']: f"Analysis from {agent['expertise']} perspective"
                for agent in self.scientist_agents
            },
            'synthesis': "Team agrees to proceed with standard workflow",
            'recommendations': [
                "Use network pharmacology for target discovery",
                "Apply molecular docking for binding prediction",
                "Implement resistance analysis"
            ],
            'answers': {},
            'next_steps': ["Execute recommended workflow"]
        }
    
    def design_workflow(self, available_tools: List[str], context: Dict[str, Any]) -> Dict[str, Any]:
        """Design optimal workflow using available tools"""
        logger.info(f"PI designing workflow with tools: {available_tools}")
        
        if self.openai_client:
            return self._ai_design_workflow(available_tools, context)
        return self._default_workflow(available_tools)
    
    def _ai_design_workflow(self, available_tools: List[str], context: Dict[str, Any]) -> Dict[str, Any]:
        """AI-powered workflow design"""
        try:
            prompt = f"""You are the Lead Scientist designing a drug discovery workflow.

AVAILABLE TOOLS:
{json.dumps(available_tools, indent=2)}

PROJECT CONTEXT:
{json.dumps(context, indent=2, default=str)}

TEAM RECOMMENDATIONS:
{json.dumps([m['recommendations'] for m in self.meeting_history[-2:]], indent=2) if self.meeting_history else "None"}

Design an optimal workflow that:
1. Uses available tools efficiently
2. Addresses the project goals
3. Follows team recommendations
4. Minimizes computational time

Return a JSON object:
{{
  "workflow_steps": [
    {{
      "step": 1,
      "tool": "tool_name",
      "purpose": "why this tool",
      "inputs": ["required inputs"],
      "outputs": ["expected outputs"]
    }}
  ],
  "rationale": "Why this workflow design",
  "estimated_time": "time estimate",
  "success_criteria": ["How to measure success"]
}}"""
            
            response = self.openai_client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {"role": "system", "content": "You are a Lead Scientist designing drug discovery workflows. Return only valid JSON."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.6,
                max_tokens=1500
            )
            
            content = response.choices[0].message.content.strip()
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0].strip()
            elif "```" in content:
                content = content.split("```")[1].split("```")[0].strip()
            
            workflow = json.loads(content)
            logger.info(f"Workflow designed with {len(workflow.get('workflow_steps', []))} steps")
            return workflow
            
        except Exception as e:
            logger.warning(f"AI workflow design failed: {e}, using default")
            return self._default_workflow(available_tools)
    
    def _default_workflow(self, available_tools: List[str]) -> Dict[str, Any]:
        """Default workflow design"""
        return {
            'workflow_steps': [
                {
                    'step': 1,
                    'tool': 'NetworkAnalyzer',
                    'purpose': 'Discover and rank drug targets',
                    'inputs': ['species_name'],
                    'outputs': ['ranked_targets', 'ppi_network']
                },
                {
                    'step': 2,
                    'tool': 'MoleculeGenerator',
                    'purpose': 'Design candidate molecules',
                    'inputs': ['target_structure', 'binding_site'],
                    'outputs': ['candidate_molecules']
                },
                {
                    'step': 3,
                    'tool': 'VinaWrapper',
                    'purpose': 'Predict binding affinity',
                    'inputs': ['molecules', 'target_pdb'],
                    'outputs': ['docking_scores', 'binding_poses']
                },
                {
                    'step': 4,
                    'tool': 'ResistancePredictor',
                    'purpose': 'Assess resistance risk',
                    'inputs': ['molecules', 'target_info'],
                    'outputs': ['resistance_scores', 'evolution_pathways']
                },
                {
                    'step': 5,
                    'tool': 'ADMETPredictor',
                    'purpose': 'Filter for drug-likeness',
                    'inputs': ['molecules'],
                    'outputs': ['admet_scores', 'filtered_candidates']
                }
            ],
            'rationale': 'Standard antimicrobial discovery pipeline',
            'estimated_time': '1-2 hours',
            'success_criteria': [
                'At least 5 candidates with binding affinity < -7 kcal/mol',
                'Low resistance risk (< 20%)',
                'Pass Lipinski rules'
            ]
        }
    
    def summarize_meeting(self, meeting_result: Dict[str, Any]) -> str:
        """Create detailed meeting summary for future reference"""
        summary = f"""
MEETING SUMMARY
===============
Agenda: {meeting_result.get('agenda', 'N/A')}

PI Opening: {meeting_result.get('pi_opening', 'N/A')}

Team Input:
{json.dumps(meeting_result.get('team_input', {}), indent=2)}

Synthesis: {meeting_result.get('synthesis', 'N/A')}

Recommendations:
{chr(10).join(['- ' + r for r in meeting_result.get('recommendations', [])])}

Next Steps:
{chr(10).join(['- ' + s for s in meeting_result.get('next_steps', [])])}
"""
        return summary
