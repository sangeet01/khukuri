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
        print("✓ Core modules")
        
        from src.admet import calculate_drug_likeness, predict_admet
        print("✓ ADMET modules")
        
        from src.target_discovery import NetworkAnalyzer, TargetRanker
        print("✓ Target discovery modules")
        
        from src.molecule_design import MoleculeGenerator, PropertyOptimizer
        print("✓ Molecule design modules")
        
        from src.docking import ReceptorPrep, BindingSiteDetector
        print("✓ Docking modules")
        
        from src.resistance import ResistancePredictor, MultiTargetDesigner
        print("✓ Resistance modules")
        
        from src.synthesis import RetroSynthesizer, RoutePlanner
        print("✓ Synthesis modules")
        
        from src.agents import BaseAgent, VirtualLab
        print("✓ Agent modules")
        
        from src.workflows import run_autonomous_discovery
        print("✓ Workflow modules")
        
        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
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
        
        print(f"✓ Molecular scoring works (score: {score:.3f})")
        return True
    except Exception as e:
        print(f"✗ Functionality test failed: {e}")
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
        print("✓ All tests passed! Khukuri is ready to use.")
        return 0
    else:
        print("✗ Some tests failed. Check installation.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
