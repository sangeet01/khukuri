#!/usr/bin/env python
"""Quick smoke test for Khukuri installation"""

import sys
from pathlib import Path

# Add parent directory to path
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

def test_imports():
    """Test all module imports"""
    print("Testing imports...")
    
    try:
        from src.core import MolecularScorer, Validator, setup_logger
        print("[OK] Core modules")
        
        from src.admet import calculate_drug_likeness, predict_admet
        print("[OK] ADMET modules")
        
        from src.target_discovery import NetworkAnalyzer, TargetRanker
        print("[OK] Target discovery modules")
        
        from src.molecule_design import MoleculeGenerator, PropertyOptimizer
        print("[OK] Molecule design modules")
        
        from src.docking import ReceptorPrep, BindingSiteDetector
        print("[OK] Docking modules")
        
        from src.resistance import ResistancePredictor, MultiTargetDesigner, PincerEngine
        print("[OK] Resistance modules (including PINCER)")
        
        from src.synthesis import RetroSynthesizer, RoutePlanner
        print("[OK] Synthesis modules")
        
        from src.agents import BaseAgent, VirtualLab
        print("[OK] Agent modules")
        
        from src.workflows import run_autonomous_discovery
        print("[OK] Workflow modules")
        
        return True
    except ImportError as e:
        print(f"[FAIL] Import failed: {e}")
        return False


def test_basic_functionality():
    """Test basic functionality"""
    print("\nTesting basic functionality...")
    
    try:
        from rdkit import Chem
        from src.core import MolecularScorer
        
        mol = Chem.MolFromSmiles('CCO')
        scorer = MolecularScorer()
        score, _, _ = scorer.calculate_composite_score(mol)
        
        print(f"[OK] Molecular scoring works (score: {score:.3f})")
        
        # Test PINCER core with ThreatAwareFitnessFunction and HGT
        from src.resistance import PincerEngine, ThreatAwareFitnessFunction, make_pincer_with_hgt
        from src.world_model import KnowledgeGraph
        
        kg = KnowledgeGraph()
        pincer = make_pincer_with_hgt(population_size=10, n_generations=2)
        pincer.map_threats("AMILV", [0, 1, 2])
        pincer.map_hgt_threats(kg, target_strain="S. aureus")
        
        pincer.evolve(["c1ccccc1"], fitness_fn=ThreatAwareFitnessFunction())
        print("[OK] PINCER Dual Red Team (Vertical + Horizontal) works")
        
        return True
    except Exception as e:
        print(f"[FAIL] Functionality test failed: {e}")
        return False


def main():
    """Run all tests"""
    print("=" * 60)
    print("Khukuri Quick Test")
    print("=" * 60)
    
    results = []
    results.append(test_imports())
    results.append(test_basic_functionality())
    
    print("\n" + "=" * 60)
    if all(results):
        print("[SUCCESS] All tests passed! Khukuri is ready to use.")
        return 0
    else:
        print("[ERROR] Some tests failed. Check installation.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
