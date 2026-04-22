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
    
    def get_pocket_residues(self, center, radius=10.0):
        """Get residues within a radius of the specified center"""
        pocket_residues = []
        try:
            structure = self.parser.get_structure('protein', self.pdb_file)
            for model in structure:
                for chain in model:
                    for residue in chain:
                        # Check distance of any atom in residue to center
                        for atom in residue:
                            dist = np.linalg.norm(atom.get_coord() - center)
                            if dist <= radius:
                                # Get 1-letter code
                                res_name = residue.get_resname()
                                res_id = residue.get_id()[1]
                                chain_id = chain.get_id()
                                pocket_residues.append({
                                    'id': res_id,
                                    'chain': chain_id,
                                    'name': res_name,
                                    'code': self._get_one_letter(res_name)
                                })
                                break
            
            logger.info(f"Extracted {len(pocket_residues)} pocket residues")
            return pocket_residues
        except Exception as e:
            logger.error(f"Failed to get pocket residues: {e}")
            return []

    def get_sequence(self):
        """Extract full amino acid sequence from PDB"""
        seq = []
        try:
            structure = self.parser.get_structure('protein', self.pdb_file)
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id()[0] == ' ': # Standard residue
                            seq.append(self._get_one_letter(residue.get_resname()))
            return "".join(seq)
        except:
            return ""

    def _get_one_letter(self, resname):
        """Convert 3-letter residue name to 1-letter code"""
        map_3to1 = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        return map_3to1.get(resname, 'X')

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
