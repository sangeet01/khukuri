"""Autonomous drug discovery workflow"""

import logging
from ..target_discovery import NetworkAnalyzer, TargetRanker, DataFetcher
from ..molecule_design import MoleculeGenerator, PropertyOptimizer
from ..docking import ReceptorPrep, BindingSiteDetector, VinaWrapper
from ..admet import predict_admet
from ..resistance import ResistancePredictor

logger = logging.getLogger('khukuri')


def run_autonomous_discovery(species, target_type='bacterial', max_compounds=50):
    """Run complete autonomous discovery workflow"""
    logger.info(f"Starting autonomous discovery for {species}")
    
    # 1. Target Discovery
    fetcher = DataFetcher()
    ppi_data = fetcher.fetch_string_ppi(species)
    
    analyzer = NetworkAnalyzer()
    network = analyzer.build_ppi_network(ppi_data)
    
    ranker = TargetRanker(analyzer)
    targets = ranker.rank_targets()
    
    if not targets:
        logger.error("No targets identified")
        return None
    
    top_target = targets[0]
    logger.info(f"Top target: {top_target['protein']}")
    
    # 2. Molecule Design
    generator = MoleculeGenerator()
    molecules = generator.generate_library(max_compounds=max_compounds)
    
    optimizer = PropertyOptimizer()
    optimized = optimizer.optimize_admet(molecules, {'min_drug_likeness': 0.5})
    
    logger.info(f"Generated {len(optimized)} optimized molecules")
    
    # 3. Results
    return {
        'targets': targets[:10],
        'molecules': optimized[:20],
        'top_target': top_target,
        'network_stats': analyzer.get_network_stats()
    }
