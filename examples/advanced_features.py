"""
Example: Advanced Features
Demonstrates Failure Analysis, Competitive Teams, and Explainable AI
"""

import os
import sys
import json

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.agents import FailureAnalyzer, CompetitiveTeams, ExplainableAI


def demo_failure_analysis():
    """Demo: Learning from experimental failures"""
    print("=" * 70)
    print("FEATURE 1: FAILURE ANALYSIS")
    print("=" * 70)
    print()
    
    analyzer = FailureAnalyzer(openai_client=None)
    
    # Simulate compound failures
    failures = [
        {
            'compound': {
                'name': 'Compound_A',
                'scaffold': 'quinolone',
                'predicted_properties': {'activity': 8.5, 'solubility': 0.05}
            },
            'experimental': {
                'measured_properties': {'activity': 3.2, 'solubility': 0.005, 'toxicity': 0.7}
            }
        },
        {
            'compound': {
                'name': 'Compound_B',
                'scaffold': 'quinolone',
                'predicted_properties': {'activity': 7.8, 'solubility': 0.08}
            },
            'experimental': {
                'measured_properties': {'activity': 2.5, 'solubility': 0.01, 'toxicity': 0.6}
            }
        },
        {
            'compound': {
                'name': 'Compound_C',
                'scaffold': 'beta-lactam',
                'predicted_properties': {'activity': 9.0, 'solubility': 0.5}
            },
            'experimental': {
                'measured_properties': {'activity': 8.8, 'solubility': 0.45, 'toxicity': 0.1}
            }
        }
    ]
    
    print("Analyzing experimental failures...\n")
    
    for failure in failures:
        analysis = analyzer.analyze_failure(
            failure['compound'],
            failure['experimental']
        )
        
        print(f"Compound: {failure['compound']['name']}")
        print(f"  Root Cause: {analysis['root_cause']}")
        print(f"  Contributing Factors: {', '.join(analysis['contributing_factors'])}")
        print(f"  Recommendations: {analysis['recommendations']['avoid'][0]}")
        print()
    
    # Get learned constraints
    print("\nLEARNED CONSTRAINTS:")
    print("-" * 70)
    constraints = analyzer.get_design_constraints()
    print(f"Scaffolds to Avoid: {constraints['avoid_scaffolds']}")
    print(f"Common Failures: {constraints['common_failure_modes']}")
    print()
    
    # Check if new compound should be synthesized
    new_compound = {
        'name': 'Compound_D',
        'scaffold': 'quinolone',  # Known to fail
        'predicted_properties': {'activity': 8.0, 'solubility': 0.1}
    }
    
    decision = analyzer.should_synthesize(new_compound)
    print(f"Should synthesize {new_compound['name']}? {decision['should_synthesize']}")
    if decision['warnings']:
        print(f"Warnings: {', '.join(decision['warnings'])}")
    print()
    
    # Generate report
    print("\nFAILURE ANALYSIS REPORT:")
    print("-" * 70)
    print(analyzer.generate_report())


def demo_competitive_teams():
    """Demo: Multiple teams competing"""
    print("\n" + "=" * 70)
    print("FEATURE 2: COMPETITIVE TEAMS")
    print("=" * 70)
    print()
    
    competition = CompetitiveTeams(openai_client=None)
    
    # Create competing teams
    print("Creating competing teams...\n")
    teams = competition.create_competing_teams(
        project_description="Design antibiotic for MRSA",
        strategies=['conservative', 'aggressive', 'balanced']
    )
    
    for team_id, team_info in teams.items():
        print(f"{team_id} ({team_info['strategy']}):")
        for member in team_info['members']:
            print(f"  - {member['title']}: {member['expertise']}")
        print()
    
    # Run competition
    print("Running competition...\n")
    result = competition.run_competition(
        task="Select optimal drug discovery strategy for MRSA",
        context={
            'pathogen': 'Staphylococcus aureus',
            'resistance_profile': 'High',
            'timeline': '6 months'
        },
        evaluation_criteria=['feasibility', 'innovation', 'speed']
    )
    
    print(f"WINNER: {result['winner']}")
    print(f"Rationale: {result['evaluation']['winner_rationale']}")
    print()
    
    print("TEAM RANKINGS:")
    for i, team_id in enumerate(result['evaluation']['rankings'], 1):
        strategy = teams[team_id]['strategy']
        print(f"  {i}. {team_id} ({strategy})")
    print()
    
    # Merge best ideas
    print("Merging best ideas from all teams...\n")
    merged = competition.merge_best_ideas()
    print("MERGED RECOMMENDATIONS:")
    for rec in merged['merged_recommendations'][:5]:
        print(f"  - {rec}")
    print()
    
    # Competition summary
    print("\nCOMPETITION SUMMARY:")
    print("-" * 70)
    print(competition.get_competition_summary())


def demo_explainable_ai():
    """Demo: Explainable AI reports"""
    print("\n" + "=" * 70)
    print("FEATURE 3: EXPLAINABLE AI")
    print("=" * 70)
    print()
    
    explainer = ExplainableAI(openai_client=None)
    
    # Explain compound ranking
    print("Explaining compound rankings...\n")
    
    compounds = [
        {
            'name': 'Compound_X',
            'total_score': 8.5,
            'scores': {
                'binding_affinity': 9.2,
                'drug_likeness': 8.0,
                'resistance_risk': 0.05
            }
        },
        {
            'name': 'Compound_Y',
            'total_score': 7.8,
            'scores': {
                'binding_affinity': 8.5,
                'drug_likeness': 7.5,
                'resistance_risk': 0.15
            }
        },
        {
            'name': 'Compound_Z',
            'total_score': 7.2,
            'scores': {
                'binding_affinity': 7.8,
                'drug_likeness': 7.0,
                'resistance_risk': 0.10
            }
        }
    ]
    
    ranking_criteria = {
        'binding_affinity': 0.4,
        'drug_likeness': 0.35,
        'resistance_risk': 0.25
    }
    
    explanation = explainer.explain_ranking(compounds, ranking_criteria)
    
    print(f"SUMMARY: {explanation['summary']}")
    print()
    
    top = explanation['top_compound_explanation']
    print(f"TOP COMPOUND: {top['compound_name']}")
    print(f"Why Ranked First: {top['why_ranked_first']}")
    print(f"Key Strengths:")
    for strength in top['key_strengths']:
        print(f"  - {strength}")
    print()
    
    print("RANKING FACTORS:")
    for factor in explanation['ranking_factors']:
        print(f"  - {factor['factor']} ({factor['importance']}): {factor['explanation']}")
    print()
    
    print(f"FOR NON-EXPERTS: {explanation['layman_summary']}")
    print()
    
    # Explain a decision
    print("\nExplaining a decision...\n")
    
    decision_explanation = explainer.explain_decision(
        decision="Use network pharmacology for target discovery",
        context={
            'pathogen': 'E. coli',
            'available_data': 'PPI network, resistance genes',
            'timeline': 'urgent'
        }
    )
    
    print(f"DECISION: {decision_explanation['decision_summary']}")
    print(f"RATIONALE: {decision_explanation['rationale']}")
    print()
    print("SUPPORTING EVIDENCE:")
    for evidence in decision_explanation['supporting_evidence']:
        print(f"  - {evidence}")
    print()
    print(f"CONFIDENCE: {decision_explanation['confidence_level']}")
    print()
    
    # Generate reports
    print("\nGenerating reports...\n")
    
    results = {
        'project_name': 'MRSA Antibiotic Discovery',
        'date': '2024-01-15',
        'candidates': compounds,
        'success_probability': '75%',
        'recommendation': 'Proceed with top 3 candidates',
        'estimated_timeline': '4-6 months',
        'estimated_cost': '$50K-$100K'
    }
    
    print("EXECUTIVE SUMMARY:")
    print("-" * 70)
    print(explainer.generate_report(results, report_type='executive'))


def main():
    """Run all demos"""
    
    print("\n")
    print("*" * 70)
    print("KHUKURI ADVANCED FEATURES DEMONSTRATION")
    print("*" * 70)
    print()
    
    # Demo 1: Failure Analysis
    demo_failure_analysis()
    
    # Demo 2: Competitive Teams
    demo_competitive_teams()
    
    # Demo 3: Explainable AI
    demo_explainable_ai()
    
    print("\n" + "=" * 70)
    print("ALL DEMONSTRATIONS COMPLETE")
    print("=" * 70)
    print()
    print("These features enhance Khukuri with:")
    print("  ✓ Learning from experimental failures")
    print("  ✓ Exploring diverse solutions via competition")
    print("  ✓ Transparent, explainable AI decisions")
    print()


if __name__ == "__main__":
    main()
