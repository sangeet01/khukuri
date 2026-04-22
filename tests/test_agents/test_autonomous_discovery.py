import pytest
import logging
import sys
import os

# Add project root to path for direct execution
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from src.agents import KhukuriOrchestrator, FallbackProvider

def test_autonomous_discovery_loop():
    """Verify the multi-iteration autonomous discovery loop and tool execution."""
    orch = KhukuriOrchestrator(provider=FallbackProvider(), project="test_autonomous")
    
    # Run 2 iterations to trigger both tool layers (TargetRanker and MoleculeGenerator)
    results = orch.run_discovery(
        species_name="Mycobacterium tuberculosis",
        n_iterations=2
    )
    
    # Verify iteration count
    assert results['iterations_completed'] == 2
    
    # Verify Tool Layer 1: Targets identified
    targets = orch.world_state.query_state('target')
    assert len(targets) > 0
    assert any(t['target_id'] in ['inhA', 'katG', 'ahpC'] for t in targets)
    
    # Verify Tool Layer 2: Compounds generated
    compounds = orch.world_state.query_state('compound')
    assert len(compounds) >= 10
    
    # Verify Hypotheses
    assert len(results['final_hypotheses']) >= 2
    
    # Verify Persistence (stats)
    stats = orch.get_memory_stats()
    assert 'world_state' in stats['keys']

if __name__ == "__main__":
    print("Running Autonomous Discovery Loop Tests...")
    test_autonomous_discovery_loop()
    print("Tests Passed!")
