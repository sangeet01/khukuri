"""Test network analyzer"""

import pytest
from src.target_discovery.network_analyzer import NetworkAnalyzer


class TestNetworkAnalyzer:
    
    @pytest.fixture
    def analyzer(self):
        return NetworkAnalyzer()
    
    @pytest.fixture
    def sample_ppi_data(self):
        return {
            ('ProteinA', 'ProteinB'): 0.8,
            ('ProteinB', 'ProteinC'): 0.6,
            ('ProteinA', 'ProteinC'): 0.5
        }
    
    def test_build_ppi_network(self, analyzer, sample_ppi_data):
        """Test PPI network building"""
        network = analyzer.build_ppi_network(sample_ppi_data)
        assert network.number_of_nodes() == 3
        assert network.number_of_edges() == 3
    
    def test_calculate_centrality(self, analyzer, sample_ppi_data):
        """Test centrality calculation"""
        analyzer.build_ppi_network(sample_ppi_data)
        centrality = analyzer.calculate_centrality()
        assert 'degree' in centrality
        assert 'betweenness' in centrality
    
    def test_get_network_stats(self, analyzer, sample_ppi_data):
        """Test network statistics"""
        analyzer.build_ppi_network(sample_ppi_data)
        stats = analyzer.get_network_stats()
        assert 'nodes' in stats
        assert stats['nodes'] == 3
