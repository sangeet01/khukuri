"""PPI network analysis"""

import logging
import networkx as nx
from community import community_louvain

logger = logging.getLogger('khukuri')


class NetworkAnalyzer:
    """Analyze protein-protein interaction networks"""
    
    def __init__(self):
        self.network = nx.Graph()
    
    def build_ppi_network(self, ppi_data):
        """Build PPI network from interaction data"""
        self.network = nx.Graph()
        for (protein_a, protein_b), score in ppi_data.items():
            if score > 0.4:
                self.network.add_edge(protein_a, protein_b, weight=score)
        logger.info(f"Built network: {self.network.number_of_nodes()} nodes, {self.network.number_of_edges()} edges")
        return self.network
    
    def calculate_centrality(self):
        """Calculate network centrality metrics"""
        if self.network.number_of_nodes() == 0:
            return {}
        
        return {
            'degree': nx.degree_centrality(self.network),
            'betweenness': nx.betweenness_centrality(self.network),
            'closeness': nx.closeness_centrality(self.network),
            'eigenvector': nx.eigenvector_centrality(self.network, max_iter=1000)
        }
    
    def detect_communities(self):
        """Detect functional modules using Louvain algorithm"""
        if self.network.number_of_nodes() == 0:
            return {}
        return community_louvain.best_partition(self.network)
    
    def get_network_stats(self):
        """Get network statistics"""
        if self.network.number_of_nodes() == 0:
            return {}
        
        return {
            'nodes': self.network.number_of_nodes(),
            'edges': self.network.number_of_edges(),
            'density': nx.density(self.network),
            'avg_clustering': nx.average_clustering(self.network),
            'connected_components': nx.number_connected_components(self.network)
        }
