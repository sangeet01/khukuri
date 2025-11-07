"""Tests for ResistanceGenomicsAnalyzer"""

import pytest
from src.genomics import ResistanceGenomicsAnalyzer


def test_analyzer_initialization():
    """Test analyzer initialization"""
    analyzer = ResistanceGenomicsAnalyzer()
    assert len(analyzer.known_mutations) > 0


def test_analyze_mutation_profile():
    """Test mutation profile analysis"""
    analyzer = ResistanceGenomicsAnalyzer()
    profile = analyzer.analyze_mutation_profile('rpoB', ['S531L', 'H526Y'])
    
    assert profile['gene'] == 'rpoB'
    assert len(profile['mutations']) == 2
    assert 'rifampicin' in profile['predicted_resistance']
    assert profile['confidence'] > 0


def test_predict_resistance_evolution():
    """Test resistance evolution prediction"""
    analyzer = ResistanceGenomicsAnalyzer()
    prediction = analyzer.predict_resistance_evolution(
        current_mutations={'rpoB': ['S531L']},
        drug_pressure='rifampicin'
    )
    
    assert 'likely_mutations' in prediction
    assert 'probability' in prediction
    assert 'timeline_generations' in prediction


def test_identify_compensatory_mutations():
    """Test compensatory mutation identification"""
    analyzer = ResistanceGenomicsAnalyzer()
    compensatory = analyzer.identify_compensatory_mutations('S531L', 'rpoB')
    
    assert isinstance(compensatory, list)
