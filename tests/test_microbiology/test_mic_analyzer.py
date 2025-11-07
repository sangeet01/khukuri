"""Tests for MICAnalyzer"""

import pytest
from src.microbiology import MICAnalyzer


def test_mic_analyzer_initialization():
    """Test MIC analyzer initialization"""
    analyzer = MICAnalyzer()
    assert len(analyzer.breakpoints) > 0


def test_add_mic_result():
    """Test adding MIC result"""
    analyzer = MICAnalyzer()
    analyzer.add_mic_result('COMP_001', 'STRAIN_001', 'rifampicin', 0.5)
    assert len(analyzer.mic_data) == 1


def test_interpret_mic():
    """Test MIC interpretation"""
    analyzer = MICAnalyzer()
    
    # Susceptible
    interp = analyzer._interpret_mic('rifampicin', 0.5)
    assert interp == 'susceptible'
    
    # Resistant
    interp = analyzer._interpret_mic('rifampicin', 8.0)
    assert interp == 'resistant'
    
    # Intermediate
    interp = analyzer._interpret_mic('rifampicin', 2.0)
    assert interp == 'intermediate'


def test_get_mic_profile():
    """Test getting MIC profile"""
    analyzer = MICAnalyzer()
    analyzer.add_mic_result('COMP_001', 'STRAIN_001', 'rifampicin', 0.5)
    analyzer.add_mic_result('COMP_002', 'STRAIN_001', 'isoniazid', 5.0)
    
    profile = analyzer.get_mic_profile('STRAIN_001')
    assert profile['strain_id'] == 'STRAIN_001'
    assert len(profile['results']) == 2
    assert profile['susceptible_count'] == 1
    assert profile['resistance_count'] == 1


def test_calculate_resistance_index():
    """Test resistance index calculation"""
    analyzer = MICAnalyzer()
    analyzer.add_mic_result('COMP_001', 'STRAIN_001', 'rifampicin', 0.5)
    analyzer.add_mic_result('COMP_002', 'STRAIN_001', 'isoniazid', 5.0)
    
    index = analyzer.calculate_resistance_index('STRAIN_001')
    assert 0.0 <= index <= 1.0
    assert index == 0.5  # 1 resistant out of 2
