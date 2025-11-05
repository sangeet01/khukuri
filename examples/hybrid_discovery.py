"""
Example: Hybrid Discovery System
Combines Stanford Virtual Lab's PI-led meetings with Khukuri's pre-built tools
"""

import os
import sys
import json

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.agents import HybridOrchestrator


def main():
    """Run hybrid discovery example"""
    
    print("=" * 70)
    print("KHUKURI SMART DISCOVERY SYSTEM")
    print("AI-Led Strategy + Pre-Built Tools")
    print("=" * 70)
    print()
    
    # Initialize orchestrator
    # Note: Pass OpenAI client for AI-powered meetings, or None for fallback mode
    openai_client = None  # Set to your OpenAI client if available
    
    orchestrator = HybridOrchestrator(openai_client=openai_client)
    
    # Example 1: Full Discovery Pipeline
    print("EXAMPLE 1: Full Discovery Pipeline")
    print("-" * 70)
    
    species = "Staphylococcus aureus"
    disease_context = "MRSA infection with high resistance rates"
    constraints = {
        'time_limit': '2 hours',
        'focus': 'novel mechanisms',
        'avoid': 'beta-lactams'
    }
    
    print(f"Target: {species}")
    print(f"Context: {disease_context}")
    print(f"Constraints: {constraints}")
    print()
    
    # Run discovery
    results = orchestrator.run_discovery(
        species_name=species,
        disease_context=disease_context,
        constraints=constraints
    )
    
    # Display results
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    
    print("\n1. TEAM CREATED:")
    for agent in results['team']:
        print(f"   - {agent['title']}: {agent['expertise']}")
    
    print("\n2. STRATEGY RECOMMENDATIONS:")
    for rec in results['strategy'].get('recommendations', [])[:3]:
        print(f"   - {rec}")
    
    print("\n3. WORKFLOW DESIGNED:")
    for step in results['workflow'].get('workflow_steps', [])[:5]:
        print(f"   Step {step['step']}: {step['tool']} - {step['purpose']}")
    
    print("\n4. EXECUTION RESULTS:")
    for step in results['execution_results'].get('steps_completed', [])[:5]:
        print(f"   ✓ {step['tool']}: {step['status']}")
    
    print("\n5. AGENT ANALYSIS:")
    for rec in results['agent_analysis'].get('recommendations', [])[:3]:
        print(f"   - {rec}")
    
    print("\n6. NEXT STEPS:")
    for step in results.get('next_steps', [])[:3]:
        print(f"   - {step}")
    
    # Save full report
    output_file = "hybrid_discovery_report.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n✓ Full report saved to: {output_file}")
    
    # Example 2: Quick Consultation
    print("\n\n" + "=" * 70)
    print("EXAMPLE 2: Quick Consultation with Lead Scientist")
    print("-" * 70)
    
    question = "What are the key challenges in designing antibiotics for Gram-negative bacteria?"
    print(f"Question: {question}")
    print()
    
    response = orchestrator.quick_consult(
        question=question,
        context={'pathogen_type': 'Gram-negative'}
    )
    
    print("Lead Scientist Response:")
    print(json.dumps(response, indent=2, default=str))
    
    # Example 3: Custom Meeting
    print("\n\n" + "=" * 70)
    print("EXAMPLE 3: Custom Team Meeting")
    print("-" * 70)
    
    meeting_result = orchestrator.custom_meeting(
        agenda="Evaluate combination therapy strategies for resistant pathogens",
        context={
            'pathogen': 'Pseudomonas aeruginosa',
            'resistance_profile': 'MDR',
            'current_treatments': ['carbapenems', 'colistin']
        },
        meeting_type='team'
    )
    
    print("Meeting Recommendations:")
    for rec in meeting_result.get('recommendations', [])[:3]:
        print(f"   - {rec}")
    
    print("\n" + "=" * 70)
    print("SMART DISCOVERY SYSTEM DEMONSTRATION COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
