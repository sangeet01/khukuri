"""Example: Agent-driven structure download and docking"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.docking import StructureDownloader
from src.bioknowledge import ResistanceDatabase
from src.core import setup_logger

logger = setup_logger('khukuri')


def simulate_agent_decision(pathogen: str, resistance_db: ResistanceDatabase):
    """Simulate LLM agents deciding which targets to use"""
    
    print(f"\n[Agent Meeting] Analyzing {pathogen}...")
    
    # Get resistance genes for pathogen
    genes = resistance_db.get_genes_by_organism(pathogen)
    
    print(f"[Genomics Agent] Found {len(genes)} resistance genes")
    print(f"  Genes: {', '.join(genes[:5])}")
    
    # Agents decide targets (simulated)
    if 'tuberculosis' in pathogen.lower():
        targets = [
            {'name': 'InhA', 'pdb_id': '1P44', 'organism': 'Mycobacterium tuberculosis'},
            {'name': 'KatG', 'pdb_id': '2CCA', 'organism': 'Mycobacterium tuberculosis'},
            {'name': 'RpoB', 'pdb_id': '5UH5', 'organism': 'Mycobacterium tuberculosis'}
        ]
        print("[Cheminformatics Agent] Recommending TB targets: InhA, KatG, RpoB")
    
    elif 'aureus' in pathogen.lower():
        targets = [
            {'name': 'PBP2a', 'pdb_id': '1VQQ', 'organism': 'Staphylococcus aureus'},
            {'name': 'DHFR', 'pdb_id': '2W9S', 'organism': 'Staphylococcus aureus'}
        ]
        print("[Cheminformatics Agent] Recommending S. aureus targets: PBP2a, DHFR")
    
    else:
        # Search PDB dynamically
        targets = [
            {'name': 'DNA gyrase', 'organism': pathogen},
            {'name': 'RNA polymerase', 'organism': pathogen}
        ]
        print(f"[Cheminformatics Agent] Will search PDB for {pathogen} targets")
    
    print("[Resistance Critic] Approved target selection")
    
    return targets


def main():
    print("="*60)
    print("Agent-Driven Structure Download Demo")
    print("="*60)
    
    # Initialize
    resistance_db = ResistanceDatabase()
    downloader = StructureDownloader()
    
    # Simulate agent meeting for TB
    pathogen = "Mycobacterium tuberculosis"
    targets = simulate_agent_decision(pathogen, resistance_db)
    
    # Download structures based on agent decisions
    print(f"\n[System] Downloading structures for {len(targets)} targets...")
    structures = downloader.download_batch(targets)
    
    print(f"\n[System] Downloaded {len(structures)} structures:")
    for target_name, path in structures.items():
        print(f"  ✓ {target_name}: {path.name}")
    
    # Show available structures
    available = downloader.get_available_structures()
    print(f"\n[System] Total structures available: {len(available)}")
    
    print("\n" + "="*60)
    print("Demo complete!")
    print("\nIn real workflow:")
    print("  1. Agents analyze pathogen and resistance data")
    print("  2. Agents decide which targets to use")
    print("  3. System downloads structures on-demand")
    print("  4. Docking proceeds with downloaded structures")
    print("="*60)


if __name__ == '__main__':
    main()
