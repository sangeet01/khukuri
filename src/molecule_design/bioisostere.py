"""Bioisostere replacement"""

import logging
from rdkit import Chem

logger = logging.getLogger('khukuri')


class BioisostereReplacer:
    """Replace functional groups with bioisosteres"""
    
    def __init__(self):
        self.replacements = {
            'c1ccccc1': ['c1ccncc1', 'c1ccoc1', 'c1ccsc1'],  # Benzene
            'C(=O)O': ['c1nnnn1', 'S(=O)(=O)N'],  # Carboxylic acid
            'C(=O)N': ['c1nncn1', 'CS(=O)(=O)N'],  # Amide
            'CO': ['CF', 'CN'],  # Alcohol
            'C': ['F', 'CF3']  # Methyl
        }
    
    def find_replacements(self, mol):
        """Find bioisostere replacements"""
        smiles = Chem.MolToSmiles(mol)
        variants = []
        
        for original, replacements in self.replacements.items():
            if original in smiles:
                for replacement in replacements:
                    try:
                        new_smiles = smiles.replace(original, replacement, 1)
                        new_mol = Chem.MolFromSmiles(new_smiles)
                        if new_mol:
                            variants.append(new_mol)
                    except:
                        continue
        
        logger.info(f"Generated {len(variants)} bioisostere variants")
        return variants
