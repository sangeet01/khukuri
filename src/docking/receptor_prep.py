"""Receptor preparation"""

import logging
import os
from Bio.PDB import PDBList
from pathlib import Path

logger = logging.getLogger('khukuri')


class ReceptorPrep:
    """Prepare PDB structures for docking"""
    
    def download_pdb(self, pdb_id, output_dir='.'):
        """Download PDB file"""
        try:
            pdbl = PDBList()
            pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir, file_format='pdb', overwrite=True)
            pdb_file = os.path.join(output_dir, f'pdb{pdb_id.lower()}.ent')
            
            if not os.path.exists(pdb_file):
                pdb_file = os.path.join(output_dir, f'{pdb_id.lower()}.pdb')
            
            if os.path.exists(pdb_file):
                logger.info(f"Downloaded {pdb_id}")
                return pdb_file
            
            logger.error(f"PDB file not found: {pdb_id}")
            return None
            
        except Exception as e:
            logger.error(f"PDB download failed: {e}")
            return None
    
    def clean_receptor(self, pdb_file, output_file='receptor_clean.pdb'):
        """Clean PDB file (remove HETATM, water)"""
        try:
            with open(pdb_file, 'r') as f_in, open(output_file, 'w') as f_out:
                atom_count = 0
                for line in f_in:
                    if line.startswith('ATOM') or line.startswith('TER') or line.startswith('END'):
                        f_out.write(line)
                        if line.startswith('ATOM'):
                            atom_count += 1
            
            logger.info(f"Cleaned receptor: {atom_count} atoms")
            return output_file if atom_count > 0 else None
            
        except Exception as e:
            logger.error(f"Receptor cleaning failed: {e}")
            return None
