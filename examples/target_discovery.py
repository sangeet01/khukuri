"""Target discovery workflow example"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.target_discovery import NetworkAnalyzer, TargetRanker, DataFetcher
from src.core import setup_logger

logger = setup_logger(__name__)


def main():
    """Run target discovery for tuberculosis"""
    
    seed_genes = ["inhA", "katG", "rpoB", "embB", "pncA"]
    logger.info(f"Starting target discovery with {len(seed_genes)} seed genes")
    
    fetcher = DataFetcher()
    ppi_data = fetcher.fetch_string_ppi(seed_genes, species=83332)
    
    if not ppi_data:
        logger.warning("No PPI data found, using seed genes only")
        ppi_data = [(g1, g2, 0.5) for i, g1 in enumerate(seed_genes) 
                    for g2 in seed_genes[i+1:]]
    
    analyzer = NetworkAnalyzer()
    network = analyzer.build_ppi_network(ppi_data)
    
    metrics = analyzer.calculate_network_metrics(network)
    logger.info(f"Network: {metrics['num_nodes']} nodes, {metrics['num_edges']} edges")
    
    ranker = TargetRanker()
    ranked_targets = ranker.rank_targets(network, seed_genes)
    
    logger.info("\nTop 10 Drug Targets:")
    for i, target in enumerate(ranked_targets[:10], 1):
        logger.info(f"{i}. {target['gene']}: score={target['score']:.3f}")
    
    return ranked_targets


if __name__ == "__main__":
    targets = main()
    print(f"\n✓ Target discovery completed!")
    print(f"  Identified {len(targets)} potential targets")
