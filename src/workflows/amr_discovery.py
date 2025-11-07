"""Autonomous AMR-focused drug discovery workflow"""

import logging
from typing import Dict, List, Optional
from ..bioknowledge import ResistanceDatabase, PathogenDatabase, TargetProteinDB, DatabaseUpdater
from ..genomics import ResistanceGenomicsAnalyzer, MutationTracker, OmicsIntegrator
from ..microbiology import MICAnalyzer, AssayTracker, StrainManager
from ..world_model import WorldStateTracker, KnowledgeGraph, LearningLoop, KosmosEngine, HypothesisEngine
from ..agents.amr_agents import (MicrobiologyAgent, GenomicsAgent, CheminformaticsAgent,
                                  ResistanceCriticAgent, LiteratureAgent)
from ..molecule_design import MoleculeGenerator, PropertyOptimizer
from ..docking import VinaWrapper
from ..admet import predict_admet
from ..resistance import ResistancePredictor

logger = logging.getLogger('khukuri')


class AMRDiscoveryWorkflow:
    """Complete AMR-focused autonomous discovery workflow with Kosmos-style capabilities"""
    
    def __init__(self, openai_client=None, auto_update_db: bool = True):
        # Auto-update databases
        if auto_update_db:
            self._update_databases()
        
        # Bio-knowledge layer
        self.resistance_db = ResistanceDatabase()
        self.pathogen_db = PathogenDatabase()
        self.target_db = TargetProteinDB()
        
        # Genomics layer
        self.genomics_analyzer = ResistanceGenomicsAnalyzer()
        self.mutation_tracker = MutationTracker()
        self.omics_integrator = OmicsIntegrator()
        
        # Microbiology layer
        self.mic_analyzer = MICAnalyzer()
        self.assay_tracker = AssayTracker()
        self.strain_manager = StrainManager()
        
        # World model (Kosmos-style)
        self.world_state = WorldStateTracker()
        self.knowledge_graph = KnowledgeGraph()
        self.learning_loop = LearningLoop()
        self.kosmos = KosmosEngine(self.world_state, self.knowledge_graph)
        self.hypothesis_engine = HypothesisEngine(self.world_state, self.knowledge_graph)
        
        # Agents
        self.agents = {
            'microbiology': MicrobiologyAgent(openai_client),
            'genomics': GenomicsAgent(openai_client),
            'cheminformatics': CheminformaticsAgent(openai_client),
            'resistance_critic': ResistanceCriticAgent(openai_client),
            'literature': LiteratureAgent(openai_client)
        }
        
        # Existing modules
        self.molecule_generator = MoleculeGenerator()
        self.property_optimizer = PropertyOptimizer()
        self.resistance_predictor = ResistancePredictor()
    
    def _update_databases(self):
        """Auto-update resistance databases"""
        updater = DatabaseUpdater()
        if updater.should_update(interval_days=1):
            logger.info("Updating resistance databases...")
            new_genes = updater.update_all()
            updater.merge_with_existing(new_genes)
            logger.info(f"Database updated with {len(new_genes)} new genes")
    
    def run_discovery(self, pathogen: str, priority: str = 'critical', 
                     n_compounds: int = 20, n_iterations: int = 3) -> Dict:
        """Run complete AMR discovery workflow"""
        logger.info(f"Starting AMR discovery for {pathogen}")
        
        # 1. Initial observation (Kosmos-style)
        pathogen_info = self.pathogen_db.get_pathogen_info(pathogen)
        self.kosmos.observe('pathogen_selected', pathogen_info)
        
        # 2. Target identification
        targets = self._identify_targets(pathogen, priority)
        logger.info(f"Identified {len(targets)} targets")
        
        # 3. Resistance analysis
        resistance_profile = self._analyze_resistance(pathogen)
        logger.info(f"Resistance profile: {resistance_profile}")
        
        # 4. Generate hypotheses
        hypotheses = self.hypothesis_engine.generate_hypotheses({'pathogen': pathogen})
        for hyp in hypotheses[:5]:  # Top 5
            self.world_state.add_hypothesis(
                hyp['hypothesis'],
                hyp.get('evidence_required', []),
                hyp.get('confidence', 0.5)
            )
        
        # 5. Iterative compound discovery
        results = {
            'pathogen': pathogen,
            'targets': targets,
            'resistance_profile': resistance_profile,
            'hypotheses': len(self.world_state.hypotheses),
            'iterations': []
        }
        
        for iteration in range(n_iterations):
            logger.info(f"Starting iteration {iteration + 1}/{n_iterations}")
            
            # Kosmos reasoning
            reasoning = self.kosmos.reason(
                f"What compounds show promise against {pathogen}?",
                {'pathogen': pathogen}
            )
            
            iter_results = self._run_iteration(
                pathogen=pathogen,
                targets=targets,
                n_compounds=n_compounds
            )
            iter_results['reasoning'] = reasoning
            
            results['iterations'].append(iter_results)
            self.learning_loop.next_iteration()
        
        # 6. Final recommendations and report
        results['recommendations'] = self._generate_recommendations()
        results['world_state_summary'] = self.world_state.get_state_summary()
        results['kosmos_report'] = self.kosmos.generate_report()
        results['validated_hypotheses'] = len(self.world_state.get_validated_hypotheses())
        
        logger.info("AMR discovery workflow completed")
        return results
    
    def _identify_targets(self, pathogen: str, priority: str) -> List[Dict]:
        """Identify and prioritize targets"""
        # Get pathogen info
        pathogen_info = self.pathogen_db.get_pathogen_info(pathogen)
        if not pathogen_info:
            logger.warning(f"Pathogen {pathogen} not in database")
            return []
        
        # Get essential genes and targets
        essential_genes = pathogen_info.get('essential_genes', [])
        target_proteins = pathogen_info.get('targets', [])
        
        targets = []
        for target_name in target_proteins:
            target_info = self.target_db.get_target_info(target_name)
            if target_info:
                targets.append({
                    'name': target_name,
                    'druggability': target_info.get('druggability', 0.5),
                    'pathway': target_info.get('pathway', 'unknown'),
                    'pdb_ids': target_info.get('pdb_ids', [])
                })
                
                # Add to world state
                self.world_state.update_target(target_name, target_info)
                self.knowledge_graph.add_target(target_name, target_info)
        
        # Sort by druggability
        targets.sort(key=lambda x: x['druggability'], reverse=True)
        
        return targets[:5]  # Top 5 targets
    
    def _analyze_resistance(self, pathogen: str) -> Dict:
        """Analyze resistance profile"""
        pathogen_info = self.pathogen_db.get_pathogen_info(pathogen)
        if not pathogen_info:
            return {}
        
        resistance_drugs = pathogen_info.get('resistance_profile', [])
        
        # Get resistance genes
        resistance_genes = []
        for drug in resistance_drugs:
            genes = self.resistance_db.get_genes_by_drug_class(drug)
            resistance_genes.extend(genes)
        
        profile = {
            'resistant_to': resistance_drugs,
            'resistance_genes': list(set(resistance_genes)),
            'priority': pathogen_info.get('priority', 'medium')
        }
        
        return profile
    
    def _run_iteration(self, pathogen: str, targets: List[Dict], n_compounds: int) -> Dict:
        """Run single discovery iteration"""
        # 1. Generate/select compounds
        if self.learning_loop.iteration == 0:
            # Initial generation
            compounds = self._generate_compounds(n_compounds)
        else:
            # Active learning selection
            candidates = self._generate_compounds(n_compounds * 3)
            compounds = self.learning_loop.prioritize_next_experiments(candidates, n_compounds)
        
        # 2. Multi-agent analysis
        analyzed_compounds = []
        for compound in compounds[:10]:  # Analyze top 10
            analysis = self._multi_agent_analysis(compound, targets[0], pathogen)
            compound['agent_analysis'] = analysis
            analyzed_compounds.append(compound)
        
        # 3. Predict resistance
        for compound in analyzed_compounds:
            resistance_pred = self._predict_resistance(compound, pathogen)
            compound['resistance_prediction'] = resistance_pred
        
        # 4. Update world state
        for compound in analyzed_compounds:
            self.world_state.update_compound(compound['compound_id'], compound)
            self.knowledge_graph.add_compound(compound['compound_id'], compound)
        
        return {
            'n_compounds': len(analyzed_compounds),
            'top_compounds': analyzed_compounds[:5],
            'avg_resistance_score': sum(c['resistance_prediction']['resistance_probability'] 
                                       for c in analyzed_compounds) / len(analyzed_compounds)
        }
    
    def _generate_compounds(self, n_compounds: int) -> List[Dict]:
        """Generate compound candidates"""
        molecules = self.molecule_generator.generate_library(max_compounds=n_compounds)
        
        compounds = []
        for i, mol in enumerate(molecules):
            compound_id = f"COMP_{self.learning_loop.iteration:03d}_{i:04d}"
            
            # Calculate features
            from rdkit.Chem import Descriptors
            features = {
                'mw': Descriptors.MolWt(mol),
                'logp': Descriptors.MolLogP(mol),
                'hbd': Descriptors.NumHDonors(mol),
                'hba': Descriptors.NumHAcceptors(mol)
            }
            
            compounds.append({
                'compound_id': compound_id,
                'mol': mol,
                'smiles': mol.GetSmiles() if hasattr(mol, 'GetSmiles') else str(mol),
                'features': features
            })
        
        return compounds
    
    def _multi_agent_analysis(self, compound: Dict, target: Dict, pathogen: str) -> Dict:
        """Multi-agent collaborative analysis"""
        analyses = {}
        
        # Prepare data
        data = {
            'compound_id': compound['compound_id'],
            'smiles': compound.get('smiles'),
            'target': target['name'],
            'pathogen': pathogen,
            'features': compound.get('features', {})
        }
        
        # Each agent analyzes
        for agent_name, agent in self.agents.items():
            question = f"Evaluate this compound for {pathogen} targeting {target['name']}"
            analyses[agent_name] = agent.analyze(data, question)
        
        # Synthesize consensus
        consensus = self._synthesize_consensus(analyses)
        
        return {
            'individual_analyses': analyses,
            'consensus': consensus
        }
    
    def _synthesize_consensus(self, analyses: Dict) -> Dict:
        """Synthesize consensus from agent analyses"""
        all_concerns = []
        all_recommendations = []
        
        for agent_name, analysis in analyses.items():
            all_concerns.extend(analysis.get('concerns', []))
            all_recommendations.extend(analysis.get('recommendations', []))
        
        return {
            'major_concerns': list(set(all_concerns))[:3],
            'top_recommendations': list(set(all_recommendations))[:3],
            'agent_agreement': len(analyses)
        }
    
    def _predict_resistance(self, compound: Dict, pathogen: str) -> Dict:
        """Predict resistance development"""
        mol = compound.get('mol')
        if not mol:
            return {'resistance_probability': 0.5, 'risk_level': 'unknown'}
        
        # Use existing resistance predictor
        prediction = self.resistance_predictor.predict_likelihood(
            mol=mol,
            pathogen_type=pathogen,
            target_protein='unknown'
        )
        
        return prediction
    
    def _generate_recommendations(self) -> List[Dict]:
        """Generate final recommendations"""
        recommendations = []
        
        # Find best compounds
        compounds = self.world_state.query_state('compound')
        if compounds:
            # Sort by resistance score (lower is better)
            sorted_compounds = sorted(
                compounds,
                key=lambda x: x.get('resistance_prediction', {}).get('resistance_probability', 1.0)
            )
            
            recommendations.append({
                'type': 'lead_compound',
                'compound_id': sorted_compounds[0]['compound_id'],
                'reason': 'Lowest resistance probability'
            })
        
        # Multi-target compounds
        multi_target = self.knowledge_graph.find_multi_target_compounds(min_targets=2)
        if multi_target:
            recommendations.append({
                'type': 'multi_target',
                'compound_id': multi_target[0]['compound_id'],
                'reason': f"Targets {multi_target[0]['target_count']} proteins"
            })
        
        return recommendations


def run_amr_discovery(pathogen: str, priority: str = 'critical', 
                     n_compounds: int = 20, n_iterations: int = 3,
                     auto_update_db: bool = True,
                     openai_api_key: Optional[str] = None) -> Dict:
    """Convenience function to run AMR discovery with Kosmos-style capabilities"""
    openai_client = None
    if openai_api_key:
        try:
            from openai import OpenAI
            openai_client = OpenAI(api_key=openai_api_key)
        except ImportError:
            logger.warning("OpenAI not available, using fallback analysis")
    
    workflow = AMRDiscoveryWorkflow(openai_client, auto_update_db=auto_update_db)
    return workflow.run_discovery(pathogen, priority, n_compounds, n_iterations)
