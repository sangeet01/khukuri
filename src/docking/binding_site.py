"""Binding site detection"""

import logging
import numpy as np
from Bio.PDB import PDBParser

logger = logging.getLogger('khukuri')


class BindingSiteDetector:
    """Detect binding sites in protein structures"""
    
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.parser = PDBParser(QUIET=True)
    
    def auto_detect(self):
        """Auto-detect binding site"""
        # Try ligand-based detection first
        ligand_site = self._detect_ligand_site()
        if ligand_site is not None:
            logger.info(f"Detected ligand binding site: {ligand_site}")
            return ligand_site
        
        # Fallback to geometric center
        center = self._geometric_center()
        logger.info(f"Using geometric center: {center}")
        return center
    
    def _detect_ligand_site(self):
        """Detect site from existing ligands"""
        ligand_coords = []
        try:
            with open(self.pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('HETATM') and not line[17:20].strip() in ['HOH', 'WAT']:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ligand_coords.append([x, y, z])
            
            if ligand_coords:
                return np.mean(ligand_coords, axis=0)
        except:
            pass
        return None
    
    def _geometric_center(self):
        """Calculate geometric center"""
        try:
            structure = self.parser.get_structure('protein', self.pdb_file)
            ca_coords = []
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if 'CA' in residue:
                            ca_coords.append(residue['CA'].get_coord())
            
            if ca_coords:
                return np.mean(ca_coords, axis=0)
        except:
            pass
        return np.array([0.0, 0.0, 0.0])
