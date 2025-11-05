"""Autonomous drug discovery workflow example"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.workflows import run_autonomous_discovery
from src.core import setup_logger

logger = setup_logger(__name__)


def main():
    """Run complete autonomous drug discovery workflow"""
    
    logger.info("=" * 60)
    logger.info("Autonomous Drug Discovery Workflow")
    logger.info("=" * 60)
    
    config = {
        "disease": "tuberculosis",
        "target_genes": ["inhA", "katG", "rpoB"],
        "num_candidates": 5,
        "max_iterations": 3,
        "scoring_weights": {
            "binding_affinity": 0.4,
            "drug_likeness": 0.25,
            "synthetic_accessibility": 0.15,
            "resistance_score": 0.2
        }
    }
    
    logger.info(f"Disease: {config['disease']}")
    logger.info(f"Target genes: {', '.join(config['target_genes'])}")
    logger.info(f"Generating {config['num_candidates']} candidates")
    
    logger.info("\nStarting autonomous discovery...")
    results = run_autonomous_discovery(**config)
    
    logger.info("\n" + "=" * 60)
    logger.info("Results Summary")
    logger.info("=" * 60)
    
    logger.info(f"Targets identified: {len(results.get('targets', []))}")
    logger.info(f"Molecules generated: {len(results.get('molecules', []))}")
    logger.info(f"Lead candidates: {len(results.get('leads', []))}")
    
    if results.get('leads'):
        logger.info("\nTop Lead Candidates:")
        for i, lead in enumerate(results['leads'][:3], 1):
            logger.info(f"\n{i}. SMILES: {lead['smiles']}")
            logger.info(f"   Score: {lead['score']:.3f}")
            logger.info(f"   Binding: {lead.get('binding_affinity', 'N/A')}")
            logger.info(f"   Drug-likeness: {lead.get('drug_likeness', 'N/A')}")
    
    return results


if __name__ == "__main__":
    results = main()
    print(f"\n✓ Autonomous workflow completed!")
    print(f"  Generated {len(results.get('leads', []))} lead candidates")
