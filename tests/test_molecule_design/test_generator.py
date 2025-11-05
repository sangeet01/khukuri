"""Test molecule generator"""

import pytest
from src.molecule_design.generator import MoleculeGenerator


class TestMoleculeGenerator:
    
    @pytest.fixture
    def generator(self):
        return MoleculeGenerator()
    
    def test_generate_library(self, generator):
        """Test library generation"""
        mols = generator.generate_library(max_compounds=10)
        assert len(mols) > 0
        assert all(mol is not None for mol in mols)
    
    def test_generate_3d_conformer(self, generator):
        """Test 3D conformer generation"""
        from rdkit import Chem
        mol = Chem.MolFromSmiles('CCO')
        mol_3d = generator.generate_3d_conformer(mol)
        assert mol_3d is not None
