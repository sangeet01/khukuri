"""Property optimization"""

import logging
from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger('khukuri')


class PropertyOptimizer:
    """Optimize molecular properties"""
    
    def optimize_structure(self, mol):
        """Optimize 3D structure"""
        try:
            mol_h = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol_h, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol_h)
            return Chem.RemoveHs(mol_h)
        except:
            return mol
    
    def optimize_admet(self, molecules, target_properties):
        """Filter molecules by ADMET properties"""
        from ..admet import calculate_drug_likeness, predict_bioavailability
        
        optimized = []
        for mol in molecules:
            drug_likeness = calculate_drug_likeness(mol)
            bioavailability = predict_bioavailability(mol)
            
            if drug_likeness >= target_properties.get('min_drug_likeness', 0.5):
                if bioavailability >= target_properties.get('min_bioavailability', 0.3):
                    optimized.append(mol)
        
        logger.info(f"Optimized to {len(optimized)} molecules")
        return optimized
