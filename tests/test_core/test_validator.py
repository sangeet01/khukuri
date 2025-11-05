"""Test validator module"""

import pytest
from rdkit import Chem
from src.core.validator import Validator


class TestValidator:
    
    def test_validate_smiles_valid(self):
        """Test valid SMILES validation"""
        assert Validator.validate_smiles('CCO') == True
        assert Validator.validate_smiles('c1ccccc1') == True
    
    def test_validate_smiles_invalid(self):
        """Test invalid SMILES validation"""
        assert Validator.validate_smiles('') == False
        assert Validator.validate_smiles(None) == False
        assert Validator.validate_smiles('INVALID') == False
    
    def test_validate_molecule_valid(self):
        """Test valid molecule validation"""
        mol = Chem.MolFromSmiles('CCO')
        assert Validator.validate_molecule(mol) == True
    
    def test_validate_molecule_invalid(self):
        """Test invalid molecule validation"""
        assert Validator.validate_molecule(None) == False
    
    def test_validate_pdb_id_valid(self):
        """Test valid PDB ID"""
        assert Validator.validate_pdb_id('4DXD') == True
        assert Validator.validate_pdb_id('1ABC') == True
    
    def test_validate_pdb_id_invalid(self):
        """Test invalid PDB ID"""
        assert Validator.validate_pdb_id('') == False
        assert Validator.validate_pdb_id('TOOLONG') == False
        assert Validator.validate_pdb_id('12') == False
    
    def test_validate_binding_affinity(self):
        """Test binding affinity validation"""
        assert Validator.validate_binding_affinity(-7.5) == True
        assert Validator.validate_binding_affinity(None) == False
    
    def test_validate_coordinates(self):
        """Test coordinate validation"""
        assert Validator.validate_coordinates([1.0, 2.0, 3.0]) == True
        assert Validator.validate_coordinates([1.0, 2.0]) == False
        assert Validator.validate_coordinates(None) == False
