"""Pytest configuration"""

import pytest


def pytest_configure(config):
    """Configure pytest markers"""
    config.addinivalue_line("markers", "integration: mark test as integration test")
    config.addinivalue_line("markers", "slow: mark test as slow running")
    config.addinivalue_line("markers", "requires_network: mark test as requiring network")


@pytest.fixture
def sample_smiles():
    """Sample SMILES strings for testing"""
    return [
        'CCO',  # Ethanol
        'CC(=O)Oc1ccccc1C(=O)O',  # Aspirin
        'c1ccccc1',  # Benzene
        'c1ccncc1'  # Pyridine
    ]


@pytest.fixture
def sample_ppi_data():
    """Sample PPI data for testing"""
    return {
        ('ProteinA', 'ProteinB'): 0.8,
        ('ProteinB', 'ProteinC'): 0.6,
        ('ProteinA', 'ProteinC'): 0.5,
        ('ProteinC', 'ProteinD'): 0.7
    }
