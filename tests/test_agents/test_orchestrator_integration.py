import pytest
import sys
import os
import logging

# Add project root to path for direct execution
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from src.agents import KhukuriOrchestrator, FallbackProvider

# Configure logging for test visibility
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('khukuri')

def test_orchestrator_initialization():
    """Verify that the orchestrator loads all agents and providers correctly."""
    provider = FallbackProvider()
    orch = KhukuriOrchestrator(provider=provider, project="test_init")
    
    assert orch.provider.name == "FallbackProvider"
    assert "chemist" in orch.agents
    assert "biologist" in orch.agents
    assert hasattr(orch, 'world_state')

def test_orchestrator_discovery_flow():
    """Verify a single-round discovery flow with world model integration."""
    orch = KhukuriOrchestrator(provider=FallbackProvider(), project="test_flow")
    
    results = orch.run_discovery(
        species_name="Mycobacterium tuberculosis",
        n_iterations=1,
        priority="critical"
    )
    
    assert results['species'] == "Mycobacterium tuberculosis"
    assert results['iterations_completed'] == 1
    assert 'world_state' in results

if __name__ == "__main__":
    print("Running Orchestrator Integration Tests...")
    test_orchestrator_initialization()
    test_orchestrator_discovery_flow()
    print("Tests Passed!")
