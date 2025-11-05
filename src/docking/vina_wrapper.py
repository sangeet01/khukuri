"""AutoDock Vina wrapper"""

import logging
import os
import subprocess
from rdkit import Chem

logger = logging.getLogger('khukuri')


class VinaWrapper:
    """Interface to AutoDock Vina"""
    
    def __init__(self, receptor_pdbqt, binding_site, box_size=(40, 40, 40), exhaustiveness=8):
        self.receptor = receptor_pdbqt
        self.center = binding_site
        self.box_size = box_size
        self.exhaustiveness = exhaustiveness
    
    def dock_molecule(self, mol, output_prefix='ligand'):
        """Dock single molecule"""
        try:
            # Save as SDF
            sdf_file = f"{output_prefix}.sdf"
            writer = Chem.SDWriter(sdf_file)
            writer.write(mol)
            writer.close()
            
            # Convert to PDBQT
            pdbqt_file = f"{output_prefix}.pdbqt"
            cmd = f"obabel {sdf_file} -O {pdbqt_file} --partialcharge gasteiger -h"
            subprocess.run(cmd, shell=True, capture_output=True)
            
            if not os.path.exists(pdbqt_file):
                return None
            
            # Run Vina
            output_file = f"{output_prefix}_out.pdbqt"
            vina_cmd = (
                f"vina --receptor {self.receptor} --ligand {pdbqt_file} "
                f"--center_x {self.center[0]:.2f} --center_y {self.center[1]:.2f} --center_z {self.center[2]:.2f} "
                f"--size_x {self.box_size[0]} --size_y {self.box_size[1]} --size_z {self.box_size[2]} "
                f"--out {output_file} --exhaustiveness {self.exhaustiveness}"
            )
            subprocess.run(vina_cmd, shell=True, capture_output=True)
            
            # Extract score
            score = self._extract_score(output_file)
            return {'score': score, 'output_file': output_file}
            
        except Exception as e:
            logger.error(f"Docking failed: {e}")
            return None
    
    def _extract_score(self, output_file):
        """Extract Vina score"""
        if not os.path.exists(output_file):
            return None
        try:
            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith('REMARK VINA RESULT:'):
                        parts = line.split()
                        if len(parts) >= 4:
                            return float(parts[3])
        except:
            pass
        return None
