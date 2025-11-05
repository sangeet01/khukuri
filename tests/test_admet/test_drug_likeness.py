"""Test ADMET drug-likeness"""

import pytest
from rdkit import Chem
from src.admet.drug_likeness import lipinski_ro5, calculate_qed, calculate_drug_likeness


class TestDrugLikeness:
    
    @pytest.fixture
    def aspirin(self):
        return Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')
    
    @pytest.fixture
    def ethanol(self):
        return Chem.MolFromSmiles('CCO')
    
    def test_lipinski_ro5_pass(self, ethanol):
        """Test Lipinski for drug-like molecule"""
        violations = lipinski_ro5(ethanol)
        assert violations == 0
    
    def test_lipinski_ro5_violations(self):
        """Test Lipinski violations"""
        # Large molecule
        large_mol = Chem.MolFromSmiles('C' * 100)
        violations = lipinski_ro5(large_mol)
        assert violations > 0
    
    def test_calculate_qed(self, aspirin):
        """Test QED calculation"""
        qed = calculate_qed(aspirin)
        assert 0 <= qed <= 1
    
    def test_calculate_drug_likeness(self, aspirin):
        """Test overall drug-likeness"""
        score = calculate_drug_likeness(aspirin)
        assert 0 <= score <= 1
