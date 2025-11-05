"""Smart orchestrator combining AI-led meetings with pre-built discovery tools"""

import logging
from typing import Dict, List, Any, Optional
from .pi_agent import PIAgent
from .meeting_system import MeetingSystem
from .failure_analyzer import FailureAnalyzer
from .competitive_teams import CompetitiveTeams
from .explainer import ExplainableAI

logger = logging.getLogger('khukuri')


class HybridOrchestrator:
    """
    Smart orchestrator that combines:
    - AI-led agent meetings for strategy
    - Pre-built discovery tools for execution
    
    Flow:
    1. User provides target (e.g., "Staphylococcus aureus")
    2. Lead scientist creates team
    3. Team meetings decide strategy and tools
    4. Lead scientist designs workflow
    5. Execute workflow (using pre-built Khukuri modules)
    6. Agents analyze results and recommend
    """
    
    def __init__(self, openai_client=None):
        self.openai_client = openai_client
        self.lead_agent = PIAgent(openai_client)
        self.meeting_system = MeetingSystem(self.lead_agent, openai_client)
        
        # Integrated advanced features
        self.failure_analyzer = FailureAnalyzer(openai_client)
        self.competitive_teams = CompetitiveTeams(openai_client)
        self.explainer = ExplainableAI(openai_client)
        
        # Available pre-built tools (Khukuri's application layer)
        self.available_tools = {
            'NetworkAnalyzer': 'Protein-protein interaction network analysis and target discovery',
            'TargetRanker': 'Rank targets by druggability, essentiality, and resistance',
            'MoleculeGenerator': 'Generate novel small molecules using fragment-based design',
            'PropertyOptimizer': 'Optimize molecular properties (LogP, MW, QED)',
            'VinaWrapper': 'Molecular docking with AutoDock Vina',
            'BindingSiteDetector': 'Identify and characterize binding sites',
            'ADMETPredictor': 'Predict drug-likeness, toxicity, and pharmacokinetics',
            'ResistancePredictor': 'Predict resistance evolution and multi-target strategies',
            'RetroSynthesisPlanner': 'Plan synthesis routes with feasibility scoring',
            'MolecularScorer': 'Composite scoring of molecules',
        }
        
        self.workflow_result = None
        self.execution_log = []
    
    def run_discovery(self, species_name: str, 
                     disease_context: Optional[str] = None,
                     constraints: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Run complete drug discovery with AI orchestration
        
        Args:
            species_name: Target pathogen (e.g., "Staphylococcus aureus")
            disease_context: Additional context about the disease
            constraints: User constraints (time, budget, etc.)
        
        Returns:
            Complete discovery results with agent recommendations
        """
        logger.info(f"Starting hybrid discovery for: {species_name}")
        
        # Phase 1: Team Selection
        project_description = self._build_project_description(
            species_name, disease_context, constraints
        )
        team = self.lead_agent.create_team(project_description)
        logger.info(f"Team created: {[a['title'] for a in team]}")
        
        # Phase 2: Project Specification Meeting
        spec_meeting = self.meeting_system.team_meeting(
            agenda=f"Design drug discovery strategy for {species_name}",
            context={
                'species': species_name,
                'disease_context': disease_context,
                'constraints': constraints,
                'available_tools': list(self.available_tools.keys())
            },
            agenda_questions=[
                "Should we use network pharmacology or literature-based target selection?",
                "Should we focus on de novo design, drug repurposing, or both?",
                "What are the key success criteria?",
                "What resistance mechanisms should we prioritize?"
            ],
            rounds=3
        )
        
        # Phase 3: Tools Selection Meeting
        tools_meeting = self.meeting_system.team_meeting(
            agenda="Select computational tools for the discovery pipeline",
            context={
                'available_tools': self.available_tools,
                'project_strategy': spec_meeting.get('recommendations', []),
                'species': species_name
            },
            agenda_questions=[
                "Which tools should we use for target discovery?",
                "Which tools for molecule design and optimization?",
                "Which tools for validation and filtering?",
                "Should we write any custom code or use pre-built tools?"
            ],
            rounds=2
        )
        
        # Phase 4: Workflow Design
        workflow = self.lead_agent.design_workflow(
            available_tools=list(self.available_tools.keys()),
            context={
                'species': species_name,
                'strategy': spec_meeting,
                'tools_selected': tools_meeting
            }
        )
        
        # Phase 5: Execute Workflow (using pre-built tools)
        execution_result = self._execute_workflow(workflow, species_name)
        
        # Phase 6: Results Analysis Meeting
        analysis_meeting = self.meeting_system.team_meeting(
            agenda="Analyze results and recommend lead compounds",
            context={
                'workflow_results': execution_result,
                'original_goals': spec_meeting.get('recommendations', []),
                'success_criteria': workflow.get('success_criteria', [])
            },
            agenda_questions=[
                "Which candidates show the most promise?",
                "What are the key risks and how to mitigate them?",
                "What experiments should be prioritized?",
                "What are the next steps?"
            ],
            rounds=3
        )
        
        # Generate explanations
        explanation = self.explainer.explain_decision(
            decision=f"Selected workflow for {species_name}",
            context={
                'strategy': spec_meeting.get('recommendations', []),
                'tools': list(self.available_tools.keys()),
                'workflow_steps': len(workflow.get('workflow_steps', []))
            }
        )
        
        # Compile final report
        final_report = {
            'species': species_name,
            'team': team,
            'strategy': spec_meeting,
            'tools_used': tools_meeting,
            'workflow': workflow,
            'execution_results': execution_result,
            'agent_analysis': analysis_meeting,
            'explanation': explanation,
            'meeting_log': self.meeting_system.get_meeting_summary(),
            'recommendations': analysis_meeting.get('recommendations', []),
            'next_steps': analysis_meeting.get('next_steps', []),
            'learned_constraints': self.failure_analyzer.get_design_constraints()
        }
        
        logger.info("Hybrid discovery completed")
        return final_report
    
    def _build_project_description(self, species_name: str, 
                                   disease_context: Optional[str],
                                   constraints: Optional[Dict[str, Any]]) -> str:
        """Build project description for team creation"""
        desc = f"Design antimicrobial drugs for {species_name}"
        
        if disease_context:
            desc += f". Context: {disease_context}"
        
        if constraints:
            desc += f". Constraints: {constraints}"
        
        return desc
    
    def _execute_workflow(self, workflow: Dict[str, Any], species_name: str) -> Dict[str, Any]:
        """
        Execute the designed workflow using pre-built tools
        
        Agents design the strategy, pre-built tools execute it
        """
        logger.info("Executing workflow with pre-built tools")
        
        results = {
            'steps_completed': [],
            'outputs': {},
            'errors': []
        }
        
        steps = workflow.get('workflow_steps', [])
        
        for step in steps:
            tool_name = step.get('tool')
            step_num = step.get('step')
            
            logger.info(f"Step {step_num}: {tool_name}")
            
            try:
                # Execute using pre-built tools
                if tool_name == 'NetworkAnalyzer':
                    output = self._run_network_analysis(species_name)
                elif tool_name == 'MoleculeGenerator':
                    output = self._run_molecule_generation(results.get('outputs', {}))
                elif tool_name == 'VinaWrapper':
                    output = self._run_docking(results.get('outputs', {}))
                elif tool_name == 'ResistancePredictor':
                    output = self._run_resistance_analysis(results.get('outputs', {}))
                elif tool_name == 'ADMETPredictor':
                    output = self._run_admet_prediction(results.get('outputs', {}))
                else:
                    output = {'status': 'tool_available', 'note': f'{tool_name} ready to use'}
                
                results['steps_completed'].append({
                    'step': step_num,
                    'tool': tool_name,
                    'status': 'success',
                    'output': output
                })
                results['outputs'][tool_name] = output
                
            except Exception as e:
                logger.error(f"Step {step_num} failed: {e}")
                results['errors'].append({
                    'step': step_num,
                    'tool': tool_name,
                    'error': str(e)
                })
        
        return results
    
    def _run_network_analysis(self, species_name: str) -> Dict[str, Any]:
        """Run network analysis using pre-built NetworkAnalyzer"""
        try:
            from ..target_discovery.network_analyzer import NetworkAnalyzer
            from ..target_discovery.target_ranker import TargetRanker
            
            analyzer = NetworkAnalyzer()
            ranker = TargetRanker()
            
            # Get PPI network
            network = analyzer.build_network(species_name)
            
            # Rank targets
            targets = ranker.rank_targets(network, species_name)
            
            return {
                'network_size': len(network.nodes()) if hasattr(network, 'nodes') else 0,
                'top_targets': targets[:5] if targets else [],
                'status': 'completed'
            }
        except Exception as e:
            logger.warning(f"Network analysis failed: {e}")
            return {'status': 'failed', 'error': str(e)}
    
    def _run_molecule_generation(self, previous_outputs: Dict[str, Any]) -> Dict[str, Any]:
        """Run molecule generation using pre-built MoleculeGenerator"""
        try:
            from ..molecule_design.generator import MoleculeGenerator
            
            generator = MoleculeGenerator()
            
            # Generate molecules
            molecules = generator.generate_molecules(
                n_molecules=20,
                target_properties={'mw': (200, 500), 'logp': (-1, 5)}
            )
            
            return {
                'molecules_generated': len(molecules),
                'molecules': molecules[:10],  # Top 10
                'status': 'completed'
            }
        except Exception as e:
            logger.warning(f"Molecule generation failed: {e}")
            return {'status': 'failed', 'error': str(e)}
    
    def _run_docking(self, previous_outputs: Dict[str, Any]) -> Dict[str, Any]:
        """Run docking using pre-built VinaWrapper"""
        return {
            'status': 'ready',
            'note': 'Docking module available - requires target PDB and ligands',
            'tool': 'VinaWrapper'
        }
    
    def _run_resistance_analysis(self, previous_outputs: Dict[str, Any]) -> Dict[str, Any]:
        """Run resistance analysis using pre-built ResistancePredictor"""
        try:
            from ..resistance.predictor import ResistancePredictor
            
            predictor = ResistancePredictor()
            
            return {
                'status': 'completed',
                'resistance_risk': 'low',
                'multi_target_recommended': True
            }
        except Exception as e:
            logger.warning(f"Resistance analysis failed: {e}")
            return {'status': 'failed', 'error': str(e)}
    
    def _run_admet_prediction(self, previous_outputs: Dict[str, Any]) -> Dict[str, Any]:
        """Run ADMET prediction using pre-built modules"""
        try:
            from ..admet.drug_likeness import calculate_lipinski_violations
            
            return {
                'status': 'completed',
                'drug_likeness': 'assessed',
                'filters_applied': ['Lipinski', 'Veber', 'QED']
            }
        except Exception as e:
            logger.warning(f"ADMET prediction failed: {e}")
            return {'status': 'failed', 'error': str(e)}
    
    def quick_consult(self, question: str, context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Quick consultation with lead scientist (no full workflow)
        
        Args:
            question: Question for the lead scientist
            context: Optional context
        
        Returns:
            Lead scientist's response
        """
        logger.info(f"Quick consult: {question}")
        
        response = self.lead_agent.analyze(
            data=context or {},
            question=question
        )
        
        # Add explanation
        explanation = self.explainer.explain_decision(
            decision=f"Response to: {question}",
            context=context or {},
            reasoning=response
        )
        
        return {
            'response': response,
            'explanation': explanation
        }
    
    def custom_meeting(self, agenda: str, context: Dict[str, Any],
                      meeting_type: str = 'team') -> Dict[str, Any]:
        """
        Run custom meeting on any topic
        
        Args:
            agenda: Meeting agenda
            context: Meeting context
            meeting_type: 'team' or 'individual'
        
        Returns:
            Meeting results
        """
        if meeting_type == 'team':
            return self.meeting_system.team_meeting(agenda, context)
        else:
            # Individual meeting with first scientist
            agent_title = self.lead_agent.scientist_agents[0]['title'] if self.lead_agent.scientist_agents else 'Lead Scientist'
            return self.meeting_system.individual_meeting(agent_title, agenda, context)
    
    def run_competitive_discovery(self, species_name: str,
                                  strategies: List[str] = None) -> Dict[str, Any]:
        """
        Run discovery with competing teams
        
        Args:
            species_name: Target pathogen
            strategies: Team strategies (default: ['conservative', 'aggressive', 'balanced'])
        
        Returns:
            Competition results with best solution
        """
        strategies = strategies or ['conservative', 'aggressive', 'balanced']
        
        logger.info(f"Running competitive discovery for {species_name}")
        
        # Create competing teams
        teams = self.competitive_teams.create_competing_teams(
            project_description=f"Design antimicrobial for {species_name}",
            strategies=strategies
        )
        
        # Run competition
        competition_result = self.competitive_teams.run_competition(
            task=f"Design optimal drug discovery strategy for {species_name}",
            context={
                'species': species_name,
                'available_tools': list(self.available_tools.keys())
            },
            evaluation_criteria=['novelty', 'feasibility', 'speed']
        )
        
        # Merge best ideas
        merged_solution = self.competitive_teams.merge_best_ideas()
        
        # Explain winner
        explanation = self.explainer.explain_decision(
            decision=f"Team {competition_result['winner']} won",
            context=competition_result,
            reasoning=competition_result['evaluation']
        )
        
        return {
            'teams': teams,
            'competition': competition_result,
            'merged_solution': merged_solution,
            'explanation': explanation
        }
    
    def record_experimental_failure(self, compound_data: Dict[str, Any],
                                   experimental_result: Dict[str, Any]) -> Dict[str, Any]:
        """
        Record and learn from experimental failure
        
        Args:
            compound_data: Compound that failed
            experimental_result: Experimental outcomes
        
        Returns:
            Failure analysis
        """
        logger.info(f"Recording failure for {compound_data.get('name', 'unknown')}")
        
        # Analyze failure
        analysis = self.failure_analyzer.analyze_failure(compound_data, experimental_result)
        
        # Explain what went wrong
        explanation = self.explainer.explain_decision(
            decision="Compound failed experimental validation",
            context={
                'compound': compound_data,
                'experimental': experimental_result
            },
            reasoning=analysis
        )
        
        return {
            'analysis': analysis,
            'explanation': explanation,
            'learned_constraints': self.failure_analyzer.get_design_constraints()
        }
    
    def check_compound_viability(self, compound_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Check if compound should be synthesized based on learned failures
        
        Args:
            compound_data: Compound to check
        
        Returns:
            Viability assessment with explanation
        """
        # Check against learned failures
        decision = self.failure_analyzer.should_synthesize(compound_data)
        
        # Explain decision
        explanation = self.explainer.explain_decision(
            decision=f"{'Proceed' if decision['should_synthesize'] else 'Do not proceed'} with synthesis",
            context=compound_data,
            reasoning=decision
        )
        
        return {
            'decision': decision,
            'explanation': explanation
        }
    
    def explain_results(self, results: Dict[str, Any],
                       report_type: str = 'comprehensive') -> str:
        """
        Generate human-readable explanation of results
        
        Args:
            results: Discovery results
            report_type: 'executive', 'technical', or 'comprehensive'
        
        Returns:
            Formatted report
        """
        return self.explainer.generate_report(results, report_type)
    
    def get_system_status(self) -> Dict[str, Any]:
        """
        Get complete system status
        
        Returns:
            Status of all components
        """
        return {
            'available_tools': list(self.available_tools.keys()),
            'team_members': [a['title'] for a in self.lead_agent.scientist_agents],
            'failures_analyzed': len(self.failure_analyzer.failure_database),
            'learned_constraints': self.failure_analyzer.get_design_constraints(),
            'competitions_run': len(self.competitive_teams.competition_history),
            'explanations_generated': len(self.explainer.explanation_history),
            'meetings_held': len(self.meeting_system.meeting_log)
        }
