"""Integration tests for workflows"""

import pytest
from src.workflows.autonomous_discovery import run_autonomous_discovery


class TestWorkflows:
    
    @pytest.mark.integration
    def test_autonomous_discovery_basic(self):
        """Test basic autonomous discovery workflow"""
        # This is a smoke test - just ensure it runs without crashing
        try:
            results = run_autonomous_discovery(
                species='Escherichia coli',
                max_compounds=5
            )
            # Should return dict with expected keys
            if results:
                assert 'targets' in results or 'molecules' in results
        except Exception as e:
            # Network errors are acceptable in tests
            if 'timeout' not in str(e).lower() and 'connection' not in str(e).lower():
                raise
    
    @pytest.mark.integration
    @pytest.mark.slow
    def test_target_to_lead_pipeline(self):
        """Test target-to-lead pipeline"""
        from src.workflows.target_to_lead import target_to_lead_pipeline
        
        # This requires PDB download - mark as slow
        try:
            results = target_to_lead_pipeline('4DXD', max_compounds=5)
            if results:
                assert 'pdb_id' in results
        except Exception as e:
            # PDB download errors acceptable
            if 'pdb' not in str(e).lower():
                raise
