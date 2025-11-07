"""Integration tests for AMR workflow"""

import pytest
from src.workflows import AMRDiscoveryWorkflow
from src.bioknowledge import ResistanceDatabase, PathogenDatabase
from src.world_model import WorldStateTracker, KnowledgeGraph


def test_amr_workflow_initialization():
    """Test AMR workflow initialization"""
    workflow = AMRDiscoveryWorkflow()
    
    assert workflow.resistance_db is not None
    assert workflow.pathogen_db is not None
    assert workflow.world_state is not None
    assert workflow.knowledge_graph is not None
    assert len(workflow.agents) == 5


def test_identify_targets():
    """Test target identification"""
    workflow = AMRDiscoveryWorkflow()
    targets = workflow._identify_targets('Mycobacterium tuberculosis', 'critical')
    
    assert len(targets) > 0
    assert all('druggability' in t for t in targets)


def test_analyze_resistance():
    """Test resistance analysis"""
    workflow = AMRDiscoveryWorkflow()
    profile = workflow._analyze_resistance('Staphylococcus aureus')
    
    assert 'resistant_to' in profile
    assert 'resistance_genes' in profile


def test_world_state_integration():
    """Test world state integration"""
    workflow = AMRDiscoveryWorkflow()
    
    # Add compound
    workflow.world_state.update_compound('COMP_001', {'smiles': 'CCO'})
    
    # Add to knowledge graph
    workflow.knowledge_graph.add_compound('COMP_001', {'mw': 46})
    
    # Verify
    assert 'COMP_001' in workflow.world_state.compounds
    assert 'COMP_001' in workflow.knowledge_graph.graph.nodes()


def test_multi_agent_analysis():
    """Test multi-agent analysis"""
    workflow = AMRDiscoveryWorkflow()
    
    compound = {'compound_id': 'COMP_001', 'smiles': 'CCO'}
    target = {'name': 'InhA', 'druggability': 0.9}
    
    analysis = workflow._multi_agent_analysis(compound, target, 'M. tuberculosis')
    
    assert 'individual_analyses' in analysis
    assert 'consensus' in analysis
    assert len(analysis['individual_analyses']) == 5


def test_end_to_end_workflow():
    """Test complete end-to-end workflow"""
    workflow = AMRDiscoveryWorkflow()
    
    results = workflow.run_discovery(
        pathogen='Staphylococcus aureus',
        priority='high',
        n_compounds=5,
        n_iterations=1
    )
    
    assert 'pathogen' in results
    assert 'targets' in results
    assert 'resistance_profile' in results
    assert 'iterations' in results
    assert len(results['iterations']) == 1
