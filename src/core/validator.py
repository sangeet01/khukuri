"""Data validation utilities"""

import logging
from rdkit import Chem
import requests

logger = logging.getLogger('khukuri')


class Validator:
    """Centralized validation for all data types"""
    
    @staticmethod
    def validate_molecule(mol):
        """Validate RDKit molecule object"""
        if mol is None:
            logger.error("Invalid molecule: None")
            return False
        try:
            if mol.GetNumAtoms() < 3:
                logger.warning(f"Molecule has too few atoms: {mol.GetNumAtoms()}")
                return False
            if mol.GetNumHeavyAtoms() < 2:
                logger.warning(f"Molecule has too few heavy atoms: {mol.GetNumHeavyAtoms()}")
                return False
            Chem.SanitizeMol(mol)
            return True
        except Exception as e:
            logger.error(f"Molecule validation failed: {e}")
            return False
    
    @staticmethod
    def validate_smiles(smiles):
        """Validate SMILES string"""
        if not smiles or not isinstance(smiles, str):
            logger.error(f"Invalid SMILES: {smiles}")
            return False
        try:
            mol = Chem.MolFromSmiles(smiles)
            return Validator.validate_molecule(mol)
        except Exception as e:
            logger.error(f"SMILES validation failed for '{smiles}': {e}")
            return False
    
    @staticmethod
    def validate_pdb_id(pdb_id):
        """Validate PDB ID format"""
        if not pdb_id or not isinstance(pdb_id, str):
            logger.error(f"Invalid PDB ID: {pdb_id}")
            return False
        if len(pdb_id) != 4 or not pdb_id.isalnum():
            logger.error(f"PDB ID must be 4 alphanumeric characters: {pdb_id}")
            return False
        return True
    
    @staticmethod
    def validate_pdb_exists(pdb_id):
        """Check if PDB ID exists in RCSB database"""
        try:
            url = f"https://www.rcsb.org/structure/{pdb_id}"
            response = requests.head(url, timeout=10)
            return response.status_code == 200
        except requests.RequestException:
            return False
    
    @staticmethod
    def validate_binding_affinity(affinity):
        """Validate binding affinity value"""
        if affinity is None:
            return False
        try:
            affinity = float(affinity)
            if affinity > 0 or affinity < -20:
                logger.warning(f"Unusual binding affinity: {affinity} kcal/mol")
            return True
        except (ValueError, TypeError):
            logger.error(f"Invalid binding affinity: {affinity}")
            return False
    
    @staticmethod
    def validate_coordinates(coords):
        """Validate 3D coordinates"""
        if coords is None or len(coords) != 3:
            logger.error(f"Invalid coordinates: {coords}")
            return False
        try:
            coords = [float(c) for c in coords]
            if any(abs(c) > 1000 for c in coords):
                logger.warning(f"Unusual coordinate values: {coords}")
            return True
        except (ValueError, TypeError):
            logger.error(f"Coordinate validation failed")
            return False
