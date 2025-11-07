"""Example: AMR-focused drug discovery workflow"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.workflows import run_amr_discovery, AMRDiscoveryWorkflow
from src.bioknowledge import ResistanceDatabase, PathogenDatabase
from src.microbiology import MICAnalyzer, StrainManager
from src.world_model import WorldStateTracker, KnowledgeGraph
from src.core import setup_logger

# Setup logging
logger = setup_logger('khukuri', log_file='amr_discovery.log')


def example_basic_amr_discovery():
    """Basic AMR discovery workflow"""
    print("=" * 60)
    print("Example 1: Basic AMR Discovery")
    print("=" * 60)
    
    # Run discovery for M. tuberculosis
    results = run_amr_discovery(
        pathogen='Mycobacterium tuberculosis',
        priority='critical',
        n_compounds=10,
        n_iterations=2
    )
    
    print(f"\nPathogen: {results['pathogen']}")
    print(f"Targets identified: {len(results['targets'])}")
    print(f"Iterations completed: {len(results['iterations'])}")
    print(f"\nTop targets:")
    for target in results['targets'][:3]:
        print(f"  - {target['name']} (druggability: {target['druggability']:.2f})")
    
    print(f"\nRecommendations:")
    for rec in results['recommendations']:
        print(f"  - {rec['type']}: {rec['compound_id']} - {rec['reason']}")
    
    return results


def example_resistance_analysis():
    """Analyze resistance patterns"""
    print("\n" + "=" * 60)
    print("Example 2: Resistance Pattern Analysis")
    print("=" * 60)
    
    # Initialize databases
    resistance_db = ResistanceDatabase()
    pathogen_db = PathogenDatabase()
    
    # Query resistance genes
    pathogen = 'Staphylococcus aureus'
    genes = resistance_db.get_genes_by_organism(pathogen)
    
    print(f"\nResistance genes for {pathogen}:")
    for gene in genes:
        info = resistance_db.query_gene(gene)
        print(f"  - {gene}: {info['type']} resistance via {info['mechanism']}")
    
    # Get pathogen info
    pathogen_info = pathogen_db.get_pathogen_info(pathogen)
    print(f"\nEssential targets:")
    for target in pathogen_info['targets']:
        print(f"  - {target}")


def example_mic_tracking():
    """Track MIC data"""
    print("\n" + "=" * 60)
    print("Example 3: MIC Data Tracking")
    print("=" * 60)
    
    mic_analyzer = MICAnalyzer()
    strain_manager = StrainManager()
    
    # Add test results
    mic_analyzer.add_mic_result('COMP_001', 'ATCC_25923', 'rifampicin', 0.5)
    mic_analyzer.add_mic_result('COMP_001', 'ATCC_25923', 'isoniazid', 0.1)
    mic_analyzer.add_mic_result('COMP_002', 'ATCC_25923', 'rifampicin', 2.0)
    
    # Get profile
    profile = mic_analyzer.get_mic_profile('ATCC_25923')
    
    print(f"\nMIC Profile for {profile['strain_id']}:")
    print(f"  Susceptible: {profile['susceptible_count']}")
    print(f"  Resistant: {profile['resistance_count']}")
    print(f"\nResults:")
    for result in profile['results']:
        print(f"  - {result['compound_id']} vs {result['drug']}: "
              f"{result['mic']} μg/mL ({result['interpretation']})")


def example_world_model():
    """Demonstrate world model tracking"""
    print("\n" + "=" * 60)
    print("Example 4: World Model State Tracking")
    print("=" * 60)
    
    world_state = WorldStateTracker()
    knowledge_graph = KnowledgeGraph()
    
    # Add entities
    world_state.update_compound('COMP_001', {
        'smiles': 'CC(C)Cc1ccc(cc1)C(C)C(O)=O',
        'mw': 206.28,
        'activity': 'active'
    })
    
    world_state.update_target('InhA', {
        'function': 'Enoyl-ACP reductase',
        'druggability': 0.9
    })
    
    world_state.add_assay_result('ASSAY_001', 'COMP_001', 'ATCC_25923', {
        'mic': 0.5,
        'interpretation': 'susceptible'
    })
    
    # Add to knowledge graph
    knowledge_graph.add_compound('COMP_001', {'mw': 206.28})
    knowledge_graph.add_target('InhA', {'druggability': 0.9})
    knowledge_graph.add_binding('COMP_001', 'InhA', affinity=-8.5, confidence=0.85)
    
    # Query state
    summary = world_state.get_state_summary()
    print(f"\nWorld State Summary:")
    print(f"  Compounds: {summary['compounds']}")
    print(f"  Targets: {summary['targets']}")
    print(f"  Assays: {summary['assays']}")
    print(f"  Hypotheses: {summary['hypotheses']['total']}")
    
    # Query knowledge graph
    targets = knowledge_graph.query_compound_targets('COMP_001')
    print(f"\nCompound COMP_001 targets:")
    for target in targets:
        print(f"  - {target['target_id']}: {target['affinity']} kcal/mol")
    
    kg_stats = knowledge_graph.get_statistics()
    print(f"\nKnowledge Graph:")
    print(f"  Nodes: {kg_stats['total_nodes']}")
    print(f"  Edges: {kg_stats['total_edges']}")


def example_full_workflow():
    """Complete AMR discovery workflow with all components"""
    print("\n" + "=" * 60)
    print("Example 5: Complete AMR Workflow")
    print("=" * 60)
    
    workflow = AMRDiscoveryWorkflow()
    
    # Run discovery
    results = workflow.run_discovery(
        pathogen='Pseudomonas aeruginosa',
        priority='critical',
        n_compounds=15,
        n_iterations=2
    )
    
    print(f"\nDiscovery Results:")
    print(f"  Pathogen: {results['pathogen']}")
    print(f"  Priority: {results['resistance_profile']['priority']}")
    print(f"  Resistant to: {', '.join(results['resistance_profile']['resistant_to'])}")
    
    print(f"\nIterations:")
    for i, iteration in enumerate(results['iterations'], 1):
        print(f"  Iteration {i}:")
        print(f"    Compounds analyzed: {iteration['n_compounds']}")
        print(f"    Avg resistance score: {iteration['avg_resistance_score']:.3f}")
    
    print(f"\nFinal Recommendations:")
    for rec in results['recommendations']:
        print(f"  - {rec['type'].upper()}: {rec['compound_id']}")
        print(f"    Reason: {rec['reason']}")
    
    # Export world state
    workflow.world_state.export_state('amr_world_state.json')
    print(f"\nWorld state exported to amr_world_state.json")


if __name__ == '__main__':
    print("\n" + "=" * 60)
    print("Khukuri AMR Discovery Examples")
    print("=" * 60)
    
    try:
        # Run examples
        example_basic_amr_discovery()
        example_resistance_analysis()
        example_mic_tracking()
        example_world_model()
        example_full_workflow()
        
        print("\n" + "=" * 60)
        print("All examples completed successfully!")
        print("=" * 60)
        
    except Exception as e:
        logger.error(f"Example failed: {e}", exc_info=True)
        print(f"\nError: {e}")
