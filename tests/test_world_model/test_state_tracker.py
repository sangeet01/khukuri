"""Tests for WorldStateTracker"""

import pytest
from src.world_model import WorldStateTracker


def test_state_tracker_initialization():
    """Test state tracker initialization"""
    tracker = WorldStateTracker()
    assert len(tracker.compounds) == 0
    assert len(tracker.targets) == 0
    assert len(tracker.hypotheses) == 0


def test_update_compound():
    """Test compound update"""
    tracker = WorldStateTracker()
    tracker.update_compound('COMP_001', {'smiles': 'CCO', 'mw': 46.07})
    
    assert 'COMP_001' in tracker.compounds
    assert tracker.compounds['COMP_001']['smiles'] == 'CCO'
    assert len(tracker.compounds['COMP_001']['history']) == 1


def test_update_target():
    """Test target update"""
    tracker = WorldStateTracker()
    tracker.update_target('InhA', {'function': 'Enoyl-ACP reductase'})
    
    assert 'InhA' in tracker.targets
    assert tracker.targets['InhA']['function'] == 'Enoyl-ACP reductase'


def test_add_assay_result():
    """Test adding assay result"""
    tracker = WorldStateTracker()
    tracker.update_compound('COMP_001', {'smiles': 'CCO'})
    tracker.add_assay_result('ASSAY_001', 'COMP_001', 'STRAIN_001', {'mic': 0.5})
    
    assert 'ASSAY_001' in tracker.assays
    assert 'assay_results' in tracker.compounds['COMP_001']


def test_add_hypothesis():
    """Test adding hypothesis"""
    tracker = WorldStateTracker()
    tracker.add_hypothesis('Compound binds InhA', ['docking', 'assay'], 0.85)
    
    assert len(tracker.hypotheses) == 1
    assert tracker.hypotheses[0]['confidence'] == 0.85
    assert tracker.hypotheses[0]['status'] == 'proposed'


def test_update_hypothesis_status():
    """Test updating hypothesis status"""
    tracker = WorldStateTracker()
    tracker.add_hypothesis('Test hypothesis', ['evidence'], 0.8)
    tracker.update_hypothesis_status(0, 'validated', {'result': 'confirmed'})
    
    assert tracker.hypotheses[0]['status'] == 'validated'
    assert 'results' in tracker.hypotheses[0]


def test_get_state_summary():
    """Test state summary"""
    tracker = WorldStateTracker()
    tracker.update_compound('COMP_001', {'smiles': 'CCO'})
    tracker.update_target('InhA', {'function': 'test'})
    tracker.add_hypothesis('Test', ['evidence'], 0.8)
    
    summary = tracker.get_state_summary()
    assert summary['compounds'] == 1
    assert summary['targets'] == 1
    assert summary['hypotheses']['total'] == 1


def test_query_state():
    """Test state query"""
    tracker = WorldStateTracker()
    tracker.update_compound('COMP_001', {'smiles': 'CCO', 'activity': 'active'})
    tracker.update_compound('COMP_002', {'smiles': 'CC', 'activity': 'inactive'})
    
    results = tracker.query_state('compound', {'activity': 'active'})
    assert len(results) == 1
    assert results[0]['compound_id'] == 'COMP_001'
