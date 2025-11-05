"""Test molecular scorer"""

import pytest
from rdkit import Chem
from src.core.molecular_scorer import MolecularScorer


class TestMolecularScorer:
    
    @pytest.fixture
    def scorer(self):
        return MolecularScorer()
    
    @pytest.fixture
    def test_mol(self):
        return Chem.MolFromSmiles('CCO')
    
    def test_calculate_drug_likeness(self, scorer, test_mol):
        """Test drug-likeness calculation"""
        score = scorer.calculate_drug_likeness(test_mol)
        assert 0 <= score <= 1
    
    def test_calculate_sa_score(self, scorer, test_mol):
        """Test SA score calculation"""
        score = scorer.calculate_sa_score(test_mol)
        assert 1 <= score <= 10
    
    def test_assess_toxicity_risk(self, scorer, test_mol):
        """Test toxicity assessment"""
        score = scorer.assess_toxicity_risk(test_mol)
        assert 0 <= score <= 1
    
    def test_evaluate_pharmacokinetics(self, scorer, test_mol):
        """Test PK evaluation"""
        score = scorer.evaluate_pharmacokinetics(test_mol)
        assert 0 <= score <= 1
    
    def test_calculate_composite_score(self, scorer, test_mol):
        """Test composite score"""
        score, metrics, normalized = scorer.calculate_composite_score(test_mol, binding_affinity=-7.5)
        assert 0 <= score <= 1
        assert 'drug_likeness' in metrics
        assert 'binding_affinity' in normalized
