"""Basic usage examples"""

from khukuri.workflows import run_autonomous_discovery, target_to_lead_pipeline
from khukuri.core import setup_logger

# Setup logging
logger = setup_logger('khukuri_example')

# Example 1: Autonomous discovery
print("=" * 60)
print("EXAMPLE 1: Autonomous Target Discovery")
print("=" * 60)

results = run_autonomous_discovery(
    species='Staphylococcus aureus',
    target_type='bacterial',
    max_compounds=50
)

if results:
    print(f"\nTop Target: {results['top_target']['protein']}")
    print(f"Importance Score: {results['top_target']['importance_score']:.3f}")
    print(f"Generated {len(results['molecules'])} optimized molecules")

# Example 2: Target to lead pipeline
print("\n" + "=" * 60)
print("EXAMPLE 2: Target-to-Lead Pipeline")
print("=" * 60)

lead_results = target_to_lead_pipeline(
    pdb_id='4DXD',
    max_compounds=100
)

if lead_results:
    print(f"\nPDB ID: {lead_results['pdb_id']}")
    print(f"Binding Site: {lead_results['binding_site']}")
    print(f"Lead Candidates: {len(lead_results['lead_candidates'])}")
    
    # Show top candidate
    if lead_results['lead_candidates']:
        top = lead_results['lead_candidates'][0]
        print(f"\nTop Candidate Score: {top['score']:.3f}")
        print(f"Drug-likeness: {top['metrics']['drug_likeness']:.3f}")
