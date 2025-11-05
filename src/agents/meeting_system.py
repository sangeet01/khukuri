"""Meeting system for agent collaboration"""

import logging
from typing import Dict, List, Any, Optional
from .base_agent import BaseAgent

logger = logging.getLogger('khukuri')


class MeetingSystem:
    """Manages team and individual meetings between agents"""
    
    def __init__(self, pi_agent, openai_client=None):
        self.pi_agent = pi_agent
        self.openai_client = openai_client
        self.meeting_log = []
    
    def team_meeting(self, agenda: str, context: Dict[str, Any],
                    agenda_questions: Optional[List[str]] = None,
                    agenda_rules: Optional[List[str]] = None,
                    rounds: int = 3) -> Dict[str, Any]:
        """
        Conduct team meeting with all agents
        
        Args:
            agenda: Meeting topic/goal
            context: Relevant data and information
            agenda_questions: Specific questions to answer
            agenda_rules: Rules/constraints to follow
            rounds: Number of discussion rounds
        
        Returns:
            Meeting results with recommendations and decisions
        """
        logger.info(f"Starting team meeting: {agenda}")
        
        meeting_input = {
            'agenda': agenda,
            'context': context,
            'agenda_questions': agenda_questions or [],
            'agenda_rules': agenda_rules or [],
            'rounds': rounds
        }
        
        # PI convenes the meeting
        result = self.pi_agent.convene_team_meeting(
            agenda=agenda,
            context=context,
            agenda_questions=agenda_questions,
            rounds=rounds
        )
        
        # Log the meeting
        self.meeting_log.append({
            'type': 'team_meeting',
            'input': meeting_input,
            'result': result
        })
        
        return result
    
    def individual_meeting(self, agent_title: str, task: str, 
                          context: Dict[str, Any],
                          with_critic: bool = True,
                          rounds: int = 2) -> Dict[str, Any]:
        """
        Conduct individual meeting with one agent
        
        Args:
            agent_title: Which agent to meet with
            task: Specific task for the agent
            context: Relevant data and information
            with_critic: Whether to include critic feedback
            rounds: Number of improvement rounds
        
        Returns:
            Agent's response and final answer
        """
        logger.info(f"Starting individual meeting with {agent_title}: {task}")
        
        # Find the agent
        agent_info = next(
            (a for a in self.pi_agent.scientist_agents if a['title'] == agent_title),
            None
        )
        
        if not agent_info:
            logger.warning(f"Agent {agent_title} not found, using Lead Scientist")
            agent_info = {
                'title': 'Lead Scientist',
                'expertise': self.pi_agent.expertise
            }
        
        # Create agent instance
        agent = BaseAgent(
            role=agent_info['title'],
            expertise=agent_info['expertise'],
            openai_client=self.openai_client
        )
        
        # Initial response
        response = agent.analyze(context, task)
        
        # Critic rounds if requested
        if with_critic and rounds > 0:
            critic = BaseAgent(
                role="Scientific Critic",
                expertise="providing critical feedback for scientific rigor",
                openai_client=self.openai_client
            )
            
            for round_num in range(rounds):
                critique = critic.analyze(
                    {'response': response, 'task': task},
                    "Provide critical feedback to improve this response"
                )
                
                # Agent improves based on critique
                response = agent.analyze(
                    {'previous_response': response, 'critique': critique, 'task': task},
                    "Improve your response based on the critique"
                )
        
        result = {
            'agent': agent_title,
            'task': task,
            'final_response': response,
            'rounds_completed': rounds
        }
        
        # Log the meeting
        self.meeting_log.append({
            'type': 'individual_meeting',
            'agent': agent_title,
            'task': task,
            'result': result
        })
        
        return result
    
    def parallel_meetings(self, agenda: str, context: Dict[str, Any],
                         n_parallel: int = 3,
                         temperature: float = 0.8) -> Dict[str, Any]:
        """
        Run multiple parallel meetings and merge results
        
        Args:
            agenda: Meeting topic
            context: Meeting context
            n_parallel: Number of parallel meetings
            temperature: LLM temperature for creativity
        
        Returns:
            Merged best results from all parallel meetings
        """
        logger.info(f"Running {n_parallel} parallel meetings: {agenda}")
        
        parallel_results = []
        
        # Run parallel meetings
        for i in range(n_parallel):
            result = self.pi_agent.convene_team_meeting(
                agenda=agenda,
                context=context,
                agenda_questions=None,
                rounds=2
            )
            parallel_results.append(result)
        
        # Merge results
        if self.openai_client:
            merged = self._ai_merge_results(parallel_results, agenda)
        else:
            merged = self._simple_merge_results(parallel_results)
        
        # Log
        self.meeting_log.append({
            'type': 'parallel_meetings',
            'n_parallel': n_parallel,
            'agenda': agenda,
            'merged_result': merged
        })
        
        return merged
    
    def _ai_merge_results(self, results: List[Dict[str, Any]], agenda: str) -> Dict[str, Any]:
        """AI-powered merging of parallel meeting results"""
        try:
            import json
            
            prompt = f"""You are the Lead Scientist merging results from {len(results)} parallel meetings.

AGENDA: {agenda}

PARALLEL MEETING RESULTS:
{json.dumps(results, indent=2, default=str)}

Merge these results by:
1. Identifying the best recommendations from each
2. Combining complementary insights
3. Resolving any conflicts
4. Creating a single optimal answer

Return a JSON object:
{{
  "merged_recommendations": ["Best recommendations from all meetings"],
  "key_insights": ["Most important insights"],
  "consensus_points": ["What all meetings agreed on"],
  "unique_contributions": {{"meeting_N": "unique insight"}},
  "final_decision": "The optimal path forward"
}}"""
            
            response = self.openai_client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {"role": "system", "content": "You are a Lead Scientist merging meeting results. Return only valid JSON."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.2,  # Lower temperature for consistency
                max_tokens=1500
            )
            
            content = response.choices[0].message.content.strip()
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0].strip()
            elif "```" in content:
                content = content.split("```")[1].split("```")[0].strip()
            
            return json.loads(content)
            
        except Exception as e:
            logger.warning(f"AI merge failed: {e}, using simple merge")
            return self._simple_merge_results(results)
    
    def _simple_merge_results(self, results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Simple merging without AI"""
        all_recommendations = []
        for result in results:
            all_recommendations.extend(result.get('recommendations', []))
        
        # Remove duplicates while preserving order
        unique_recommendations = []
        seen = set()
        for rec in all_recommendations:
            if rec not in seen:
                unique_recommendations.append(rec)
                seen.add(rec)
        
        return {
            'merged_recommendations': unique_recommendations[:5],  # Top 5
            'key_insights': ['Multiple perspectives considered'],
            'consensus_points': ['Standard workflow recommended'],
            'final_decision': 'Proceed with merged recommendations'
        }
    
    def get_meeting_summary(self) -> str:
        """Get summary of all meetings"""
        summary = f"MEETING LOG ({len(self.meeting_log)} meetings)\n"
        summary += "=" * 50 + "\n\n"
        
        for i, meeting in enumerate(self.meeting_log, 1):
            summary += f"{i}. {meeting['type'].upper()}\n"
            if meeting['type'] == 'team_meeting':
                summary += f"   Agenda: {meeting['input']['agenda']}\n"
                summary += f"   Recommendations: {len(meeting['result'].get('recommendations', []))}\n"
            elif meeting['type'] == 'individual_meeting':
                summary += f"   Agent: {meeting['agent']}\n"
                summary += f"   Task: {meeting['task']}\n"
            summary += "\n"
        
        return summary
