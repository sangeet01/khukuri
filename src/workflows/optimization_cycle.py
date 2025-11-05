"""Iterative optimization workflow"""

from rdkit import Chem
from ..molecule_design import PropertyOptimizer
from ..admet import predict_admet
from ..core import MolecularScorer, setup_logger

logger = setup_logger(__name__)


def run_optimization_cycle(initial_molecule, target_properties, max_iterations=10):
    """
    Run iterative optimization cycle
    
    Args:
        initial_molecule: Starting molecule (SMILES or Mol)
        target_properties: Dict of target property values
        max_iterations: Maximum optimization iterations
        
    Returns:
        Dict with optimized molecule and history
    """
    logger.info(f"Starting optimization cycle (max {max_iterations} iterations)")
    
    # Parse molecule
    if isinstance(initial_molecule, str):
        mol = Chem.MolFromSmiles(initial_molecule)
    else:
        mol = initial_molecule
    
    if mol is None:
        raise ValueError("Invalid molecule")
    
    # Initialize
    optimizer = PropertyOptimizer()
    scorer = MolecularScorer()
    
    best_mol = mol
    best_score = 0.0
    history = []
    
    # Optimization loop
    for iteration in range(max_iterations):
        logger.info(f"Iteration {iteration + 1}/{max_iterations}")
        
        # Optimize
        optimized = optimizer.optimize_admet(mol, target_properties)
        
        # Score
        score, metrics, _ = scorer.calculate_composite_score(optimized)
        
        # Track history
        history.append({
            "iteration": iteration + 1,
            "smiles": Chem.MolToSmiles(optimized),
            "score": score,
            "metrics": metrics
        })
        
        # Update best
        if score > best_score:
            best_mol = optimized
            best_score = score
            logger.info(f"New best score: {score:.3f}")
        
        # Check convergence
        if len(history) > 2:
            recent_scores = [h["score"] for h in history[-3:]]
            if max(recent_scores) - min(recent_scores) < 0.01:
                logger.info("Converged")
                break
        
        mol = optimized
    
    # Final results
    admet = predict_admet(best_mol)
    
    results = {
        "optimized_molecule": best_mol,
        "smiles": Chem.MolToSmiles(best_mol),
        "score": best_score,
        "admet": admet,
        "iterations": len(history),
        "history": history
    }
    
    logger.info(f"Optimization complete: {len(history)} iterations, score={best_score:.3f}")
    return results
