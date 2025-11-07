"""Knowledge graph for AMR relationships"""

import logging
from typing import Dict, List, Optional, Tuple
import networkx as nx

logger = logging.getLogger('khukuri')


class KnowledgeGraph:
    """Knowledge graph linking compounds, targets, genes, and resistance"""
    
    def __init__(self):
        self.graph = nx.MultiDiGraph()
    
    def add_compound(self, compound_id: str, properties: Dict):
        """Add compound node"""
        self.graph.add_node(compound_id, node_type='compound', **properties)
    
    def add_target(self, target_id: str, properties: Dict):
        """Add target protein node"""
        self.graph.add_node(target_id, node_type='target', **properties)
    
    def add_gene(self, gene_id: str, properties: Dict):
        """Add gene node"""
        self.graph.add_node(gene_id, node_type='gene', **properties)
    
    def add_strain(self, strain_id: str, properties: Dict):
        """Add strain node"""
        self.graph.add_node(strain_id, node_type='strain', **properties)
    
    def add_relationship(self, source: str, target: str, relationship_type: str, 
                        properties: Optional[Dict] = None):
        """Add relationship between entities"""
        self.graph.add_edge(source, target, 
                          relationship_type=relationship_type, 
                          **(properties or {}))
        logger.debug(f"Added relationship: {source} --[{relationship_type}]--> {target}")
    
    def add_binding(self, compound_id: str, target_id: str, affinity: float, 
                   confidence: float = 1.0):
        """Add compound-target binding relationship"""
        self.add_relationship(compound_id, target_id, 'binds_to', {
            'affinity': affinity,
            'confidence': confidence
        })
    
    def add_resistance_gene(self, strain_id: str, gene_id: str, mechanism: str):
        """Add strain-gene resistance relationship"""
        self.add_relationship(strain_id, gene_id, 'has_resistance_gene', {
            'mechanism': mechanism
        })
    
    def add_target_gene(self, target_id: str, gene_id: str):
        """Add target-gene encoding relationship"""
        self.add_relationship(gene_id, target_id, 'encodes')
    
    def query_compound_targets(self, compound_id: str) -> List[Dict]:
        """Get all targets for compound"""
        if compound_id not in self.graph:
            return []
        
        targets = []
        for target in self.graph.successors(compound_id):
            edge_data = self.graph.get_edge_data(compound_id, target)
            if edge_data:
                for key, data in edge_data.items():
                    if data.get('relationship_type') == 'binds_to':
                        targets.append({
                            'target_id': target,
                            'affinity': data.get('affinity'),
                            'confidence': data.get('confidence', 1.0)
                        })
        
        return targets
    
    def query_strain_resistance(self, strain_id: str) -> List[Dict]:
        """Get resistance genes for strain"""
        if strain_id not in self.graph:
            return []
        
        resistance_genes = []
        for gene in self.graph.successors(strain_id):
            edge_data = self.graph.get_edge_data(strain_id, gene)
            if edge_data:
                for key, data in edge_data.items():
                    if data.get('relationship_type') == 'has_resistance_gene':
                        resistance_genes.append({
                            'gene_id': gene,
                            'mechanism': data.get('mechanism')
                        })
        
        return resistance_genes
    
    def find_multi_target_compounds(self, min_targets: int = 2) -> List[Dict]:
        """Find compounds that bind multiple targets"""
        multi_target = []
        
        for node in self.graph.nodes():
            if self.graph.nodes[node].get('node_type') == 'compound':
                targets = self.query_compound_targets(node)
                if len(targets) >= min_targets:
                    multi_target.append({
                        'compound_id': node,
                        'targets': targets,
                        'target_count': len(targets)
                    })
        
        return sorted(multi_target, key=lambda x: x['target_count'], reverse=True)
    
    def find_shortest_path(self, source: str, target: str) -> Optional[List[str]]:
        """Find shortest path between entities"""
        try:
            return nx.shortest_path(self.graph, source, target)
        except nx.NetworkXNoPath:
            return None
    
    def get_neighbors(self, node_id: str, relationship_type: Optional[str] = None) -> List[str]:
        """Get neighboring nodes"""
        if node_id not in self.graph:
            return []
        
        neighbors = []
        for neighbor in self.graph.neighbors(node_id):
            if relationship_type:
                edge_data = self.graph.get_edge_data(node_id, neighbor)
                if edge_data:
                    for key, data in edge_data.items():
                        if data.get('relationship_type') == relationship_type:
                            neighbors.append(neighbor)
                            break
            else:
                neighbors.append(neighbor)
        
        return neighbors
    
    def get_subgraph(self, node_ids: List[str]) -> nx.MultiDiGraph:
        """Extract subgraph containing specified nodes"""
        return self.graph.subgraph(node_ids).copy()
    
    def get_statistics(self) -> Dict:
        """Get graph statistics"""
        node_types = {}
        for node in self.graph.nodes():
            node_type = self.graph.nodes[node].get('node_type', 'unknown')
            node_types[node_type] = node_types.get(node_type, 0) + 1
        
        relationship_types = {}
        for u, v, data in self.graph.edges(data=True):
            rel_type = data.get('relationship_type', 'unknown')
            relationship_types[rel_type] = relationship_types.get(rel_type, 0) + 1
        
        return {
            'total_nodes': self.graph.number_of_nodes(),
            'total_edges': self.graph.number_of_edges(),
            'node_types': node_types,
            'relationship_types': relationship_types
        }
