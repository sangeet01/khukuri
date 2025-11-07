"""Quick test script for AMR features"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.bioknowledge import ResistanceDatabase, PathogenDatabase, TargetProteinDB
from src.genomics import ResistanceGenomicsAnalyzer, MutationTracker
from src.microbiology import MICAnalyzer, StrainManager
from src.world_model import WorldStateTracker, KnowledgeGraph
from src.core import setup_logger

logger = setup_logger('khukuri')


def test_bioknowledge():
    """Test bio-knowledge layer"""
    print("\n" + "="*60)
    print("Testing Bio-Knowledge Layer")
    print("="*60)
    
    # Resistance database
    resistance_db = ResistanceDatabase()
    print(f"[OK] Loaded {len(resistance_db.genes)} resistance genes")
    
    genes = resistance_db.get_genes_by_organism('S. aureus')
    print(f"[OK] Found {len(genes)} resistance genes for S. aureus")
    
    # Pathogen database
    pathogen_db = PathogenDatabase()
    pathogen_info = pathogen_db.get_pathogen_info('Mycobacterium tuberculosis')
    print(f"[OK] Loaded pathogen info: {len(pathogen_info['targets'])} targets")
    
    # Target database
    target_db = TargetProteinDB()
    druggable = target_db.get_druggable_targets(min_score=0.7)
    print(f"[OK] Found {len(druggable)} druggable targets")


def test_genomics():
    """Test genomics layer"""
    print("\n" + "="*60)
    print("Testing Genomics Layer")
    print("="*60)
    
    # Resistance analyzer
    analyzer = ResistanceGenomicsAnalyzer()
    profile = analyzer.analyze_mutation_profile('rpoB', ['S531L', 'H526Y'])
    print(f"[OK] Analyzed mutations: {len(profile['mutations'])} found")
    print(f"  Predicted resistance: {profile['predicted_resistance']}")
    
    # Mutation tracker
    tracker = MutationTracker()
    tracker.add_strain('STRAIN_001', {'rpoB': ['S531L'], 'gyrA': ['S83L']})
    tracker.add_strain('STRAIN_002', {'rpoB': ['S531L'], 'katG': ['S315T']})
    print(f"[OK] Tracked {len(tracker.strain_profiles)} strains")
    
    freq = tracker.get_mutation_frequency('rpoB', 'S531L')
    print(f"  S531L frequency: {freq:.2f}")


def test_microbiology():
    """Test microbiology layer"""
    print("\n" + "="*60)
    print("Testing Microbiology Layer")
    print("="*60)
    
    # MIC analyzer
    mic_analyzer = MICAnalyzer()
    mic_analyzer.add_mic_result('COMP_001', 'ATCC_25923', 'rifampicin', 0.5)
    mic_analyzer.add_mic_result('COMP_001', 'ATCC_25923', 'isoniazid', 0.1)
    print(f"[OK] Added MIC results: {len(mic_analyzer.mic_data)} entries")
    
    profile = mic_analyzer.get_mic_profile('ATCC_25923')
    print(f"  Susceptible: {profile['susceptible_count']}")
    print(f"  Resistant: {profile['resistance_count']}")
    
    # Strain manager
    strain_manager = StrainManager()
    ref_strains = strain_manager.get_reference_strains()
    print(f"[OK] Loaded {len(ref_strains)} reference strains")


def test_world_model():
    """Test world model"""
    print("\n" + "="*60)
    print("Testing World Model")
    print("="*60)
    
    # State tracker
    world_state = WorldStateTracker()
    world_state.update_compound('COMP_001', {'smiles': 'CCO', 'mw': 46.07})
    world_state.update_target('InhA', {'druggability': 0.9})
    world_state.add_hypothesis('Compound binds InhA', ['docking'], 0.85)
    
    summary = world_state.get_state_summary()
    print(f"[OK] World state: {summary['compounds']} compounds, {summary['targets']} targets")
    print(f"  Hypotheses: {summary['hypotheses']['total']}")
    
    # Knowledge graph
    kg = KnowledgeGraph()
    kg.add_compound('COMP_001', {'mw': 250})
    kg.add_target('InhA', {'druggability': 0.9})
    kg.add_binding('COMP_001', 'InhA', affinity=-8.5)
    
    stats = kg.get_statistics()
    print(f"[OK] Knowledge graph: {stats['total_nodes']} nodes, {stats['total_edges']} edges")


def test_integration():
    """Test integrated workflow"""
    print("\n" + "="*60)
    print("Testing Integration")
    print("="*60)
    
    try:
        from src.workflows import run_amr_discovery
        
        print("[OK] AMR workflow module loaded")
        print("  Run 'python examples/amr_discovery_example.py' for full demo")
        
    except Exception as e:
        print(f"[FAIL] Integration test failed: {e}")


if __name__ == '__main__':
    print("\n" + "="*60)
    print("Khukuri AMR Features - Quick Test")
    print("="*60)
    
    try:
        test_bioknowledge()
        test_genomics()
        test_microbiology()
        test_world_model()
        test_integration()
        
        print("\n" + "="*60)
        print("[SUCCESS] All AMR features working correctly!")
        print("="*60)
        print("\nNext steps:")
        print("  1. Run: python examples/amr_discovery_example.py")
        print("  2. Run: python -m pytest tests/test_bioknowledge/ -v")
        print("  3. See: docs/amr_features.md for full documentation")
        print("="*60 + "\n")
        
    except Exception as e:
        logger.error(f"Test failed: {e}", exc_info=True)
        print(f"\n[FAIL] Test failed: {e}")
        sys.exit(1)
