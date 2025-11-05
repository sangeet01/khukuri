"""Target to lead pipeline"""

import logging
from ..molecule_design import MoleculeGenerator, PropertyOptimizer
from ..docking import ReceptorPrep, BindingSiteDetector, VinaWrapper
from ..admet import predict_admet, calculate_drug_likeness
from ..core import MolecularScorer

logger = logging.getLogger('khukuri')


def target_to_lead_pipeline(pdb_id, max_compounds=100):
    """Complete target to lead pipeline"""
    logger.info(f"Starting target-to-lead for {pdb_id}")
    
    # 1. Prepare receptor
    prep = ReceptorPrep()
    pdb_file = prep.download_pdb(pdb_id)
    if not pdb_file:
        logger.error("PDB download failed")
        return None
    
    clean_file = prep.clean_receptor(pdb_file)
    
    # 2. Detect binding site
    detector = BindingSiteDetector(clean_file)
    binding_site = detector.auto_detect()
    
    # 3. Generate molecules
    generator = MoleculeGenerator()
    molecules = generator.generate_library(max_compounds=max_compounds)
    
    # 4. Optimize
    optimizer = PropertyOptimizer()
    optimized = optimizer.optimize_admet(molecules, {'min_drug_likeness': 0.5})
    
    # 5. Score
    scorer = MolecularScorer()
    scored = []
    for mol in optimized[:20]:
        score, metrics, _ = scorer.calculate_composite_score(mol)
        scored.append({
            'molecule': mol,
            'score': score,
            'metrics': metrics
        })
    
    scored.sort(key=lambda x: x['score'], reverse=True)
    
    logger.info(f"Pipeline complete: {len(scored)} lead candidates")
    
    return {
        'pdb_id': pdb_id,
        'binding_site': binding_site,
        'lead_candidates': scored[:10]
    }
