"""Molecule generation"""

import logging
from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger('khukuri')


class MoleculeGenerator:
    """Generate drug-like molecules"""
    
    def __init__(self):
        self.generated = []
    
    def generate_library(self, scaffolds=None, substituents=None, max_compounds=100):
        """Generate molecule library"""
        if scaffolds is None:
            scaffolds = ['c1ccccc1', 'c1ccncc1', 'c1ccoc1']
        if substituents is None:
            substituents = ['F', 'Cl', 'N', 'O', 'C']
        
        molecules = []
        for scaffold in scaffolds:
            for sub in substituents:
                if len(molecules) >= max_compounds:
                    break
                try:
                    smiles = f"{scaffold}{sub}"
                    mol = Chem.MolFromSmiles(smiles)
                    if mol and mol.GetNumAtoms() >= 5:
                        molecules.append(mol)
                except:
                    continue
        
        self.generated = molecules
        logger.info(f"Generated {len(molecules)} molecules")
        return molecules
    
    def generate_3d_conformer(self, mol):
        """Generate 3D conformer"""
        try:
            mol_3d = Chem.AddHs(mol)
            if AllChem.EmbedMolecule(mol_3d, randomSeed=42) == 0:
                AllChem.UFFOptimizeMolecule(mol_3d, maxIters=1000)
                return mol_3d
        except:
            pass
        return None
