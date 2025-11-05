"""Competitive Agent Teams - multiple teams compete for best solution"""

import logging
import json
from typing import Dict, List, Any, Optional
from .pi_agent import PIAgent
from .meeting_system import MeetingSystem

logger = logging.getLogger('khukuri')


class CompetitiveTeams:
    """Manages multiple competing agent teams"""
    
    def __init__(self, openai_client=None):
        self.openai_client = openai_client
        self.teams = {}
        self.competition_history = []
    
    def create_competing_teams(self, project_description: str,
                              strategies: List[str]) -> Dict[str, Any]:
        """
        Create multiple teams with different strategies
        
        Args:
            project_description: What to work on
            strategies: List of strategy names (e.g., ['conservative', 'aggressive'])
        
        Returns:
            Dictionary of team_id -> team_info
        """
        logger.info(f"Creating {len(strategies)} competing teams")
        
        for strategy in strategies:
            team_id = f"team_{strategy}"
            
            # Create PI for this team
            pi_agent = PIAgent(self.openai_client)
            
            # Customize team based on strategy
            if strategy == 'conservative':
                context = f"{project_description}. Strategy: Prioritize safety and proven approaches."
            elif strategy == 'aggressive':
                context = f"{project_description}. Strategy: Prioritize novelty and high-risk high-reward."
            elif strategy == 'balanced':
                context = f"{project_description}. Strategy: Balance innovation with feasibility."
            else:
                context = f"{project_description}. Strategy: {strategy}"
            
            # Create team
            team_members = pi_agent.create_team(context)
            
            self.teams[team_id] = {
                'strategy': strategy,
                'pi_agent': pi_agent,
                'meeting_system': MeetingSystem(pi_agent, self.openai_client),
                'team_members': team_members,
                'results': []
            }
        
        return {tid: {'strategy': t['strategy'], 'members': t['team_members']} 
                for tid, t in self.teams.items()}
    
    def run_competition(self, task: str, context: Dict[str, Any],
                       evaluation_criteria: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Run competition where all teams solve the same task
        
        Args:
            task: Task description
            context: Task context
            evaluation_criteria: How to judge solutions
        
        Returns:
            Competition results with winner
        """
        logger.info(f"Running competition: {task}")
        
        if not self.teams:
            raise ValueError("No teams created. Call create_competing_teams first.")
        
        # Each team works on the task
        team_solutions = {}
        
        for team_id, team in self.teams.items():
            logger.info(f"Team {team_id} ({team['strategy']}) working on task...")
            
            # Team meeting to solve task
            solution = team['meeting_system'].team_meeting(
                agenda=task,
                context=context,
                rounds=3
            )
            
            team_solutions[team_id] = {
                'strategy': team['strategy'],
                'solution': solution,
                'team_members': [m['title'] for m in team['team_members']]
            }
            
            team['results'].append(solution)
        
        # Evaluate solutions
        evaluation = self._evaluate_solutions(
            team_solutions,
            evaluation_criteria or ['novelty', 'feasibility', 'impact']
        )
        
        # Record competition
        competition_record = {
            'task': task,
            'teams': list(team_solutions.keys()),
            'solutions': team_solutions,
            'evaluation': evaluation,
            'winner': evaluation['winner']
        }
        
        self.competition_history.append(competition_record)
        
        return competition_record
    
    def _evaluate_solutions(self, team_solutions: Dict[str, Any],
                           criteria: List[str]) -> Dict[str, Any]:
        """Evaluate and rank team solutions"""
        
        if self.openai_client:
            return self._ai_evaluate(team_solutions, criteria)
        return self._simple_evaluate(team_solutions)
    
    def _ai_evaluate(self, team_solutions: Dict[str, Any],
                    criteria: List[str]) -> Dict[str, Any]:
        """AI-powered solution evaluation"""
        try:
            prompt = f"""You are judging a competition between drug discovery teams.

EVALUATION CRITERIA:
{chr(10).join([f'- {c}' for c in criteria])}

TEAM SOLUTIONS:
{json.dumps(team_solutions, indent=2, default=str)}

Evaluate each team's solution on each criterion (score 1-10).
Determine the winner and explain why.

Return JSON:
{{
  "scores": {{
    "team_id": {{
      "criterion1": score,
      "criterion2": score,
      "total": total_score
    }}
  }},
  "rankings": ["team_id1", "team_id2", "team_id3"],
  "winner": "team_id",
  "winner_rationale": "Why this team won",
  "strengths_by_team": {{
    "team_id": ["strength1", "strength2"]
  }},
  "best_ideas": ["Best idea from any team"]
}}"""
            
            response = self.openai_client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {"role": "system", "content": "You are an impartial judge evaluating scientific solutions. Return only valid JSON."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.2,
                max_tokens=1500
            )
            
            content = response.choices[0].message.content.strip()
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0].strip()
            elif "```" in content:
                content = content.split("```")[1].split("```")[0].strip()
            
            return json.loads(content)
            
        except Exception as e:
            logger.warning(f"AI evaluation failed: {e}, using simple evaluation")
            return self._simple_evaluate(team_solutions)
    
    def _simple_evaluate(self, team_solutions: Dict[str, Any]) -> Dict[str, Any]:
        """Simple evaluation without AI"""
        
        scores = {}
        for team_id, solution in team_solutions.items():
            # Simple scoring based on number of recommendations
            recs = solution['solution'].get('recommendations', [])
            score = len(recs) * 2  # More recommendations = higher score
            
            scores[team_id] = {
                'recommendations_count': len(recs),
                'total': score
            }
        
        # Rank by score
        rankings = sorted(scores.keys(), key=lambda x: scores[x]['total'], reverse=True)
        winner = rankings[0] if rankings else list(team_solutions.keys())[0]
        
        return {
            'scores': scores,
            'rankings': rankings,
            'winner': winner,
            'winner_rationale': f"Most comprehensive solution with {scores[winner]['recommendations_count']} recommendations",
            'strengths_by_team': {
                tid: [f"{team_solutions[tid]['strategy']} approach"]
                for tid in team_solutions.keys()
            },
            'best_ideas': []
        }
    
    def merge_best_ideas(self) -> Dict[str, Any]:
        """Merge best ideas from all teams"""
        
        if not self.competition_history:
            return {'merged_ideas': [], 'source': 'none'}
        
        latest_competition = self.competition_history[-1]
        
        if self.openai_client:
            return self._ai_merge_ideas(latest_competition)
        return self._simple_merge_ideas(latest_competition)
    
    def _ai_merge_ideas(self, competition: Dict[str, Any]) -> Dict[str, Any]:
        """AI-powered idea merging"""
        try:
            prompt = f"""You are merging the best ideas from competing teams.

COMPETITION RESULTS:
{json.dumps(competition, indent=2, default=str)}

Create a merged solution that:
1. Takes the best ideas from each team
2. Combines complementary approaches
3. Resolves conflicts intelligently
4. Creates a superior hybrid solution

Return JSON:
{{
  "merged_recommendations": ["Best recommendations from all teams"],
  "hybrid_approach": "Description of merged strategy",
  "contributions_by_team": {{
    "team_id": "What this team contributed"
  }},
  "synergies": ["How ideas complement each other"],
  "conflicts_resolved": ["How conflicts were resolved"]
}}"""
            
            response = self.openai_client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {"role": "system", "content": "You are synthesizing ideas from multiple teams. Return only valid JSON."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.4,
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
            return self._simple_merge_ideas(competition)
    
    def _simple_merge_ideas(self, competition: Dict[str, Any]) -> Dict[str, Any]:
        """Simple idea merging"""
        
        all_recommendations = []
        contributions = {}
        
        for team_id, solution in competition['solutions'].items():
            recs = solution['solution'].get('recommendations', [])
            all_recommendations.extend(recs)
            contributions[team_id] = f"Contributed {len(recs)} recommendations"
        
        # Remove duplicates
        unique_recs = list(set(all_recommendations))
        
        return {
            'merged_recommendations': unique_recs,
            'hybrid_approach': 'Combined all team recommendations',
            'contributions_by_team': contributions,
            'synergies': ['Multiple perspectives considered'],
            'conflicts_resolved': []
        }
    
    def get_competition_summary(self) -> str:
        """Get summary of all competitions"""
        
        if not self.competition_history:
            return "No competitions run yet."
        
        summary = f"""
COMPETITION SUMMARY
===================
Total Competitions: {len(self.competition_history)}

TEAM PERFORMANCE:
-----------------
"""
        
        # Count wins per team
        wins = {}
        for comp in self.competition_history:
            winner = comp['winner']
            wins[winner] = wins.get(winner, 0) + 1
        
        for team_id, win_count in sorted(wins.items(), key=lambda x: x[1], reverse=True):
            strategy = self.teams.get(team_id, {}).get('strategy', 'unknown')
            summary += f"  {team_id} ({strategy}): {win_count} wins\n"
        
        summary += f"\nRECENT COMPETITIONS:\n"
        summary += "-" * 50 + "\n"
        
        for comp in self.competition_history[-3:]:
            summary += f"  Task: {comp['task']}\n"
            summary += f"  Winner: {comp['winner']}\n"
            summary += f"  Rationale: {comp['evaluation'].get('winner_rationale', 'N/A')}\n\n"
        
        return summary
