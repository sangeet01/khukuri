import ctypes
import os
import numpy as np
from typing import Tuple, List, Optional, Dict
import time

# --- C Constants ---
N_CHANNELS = 11
CH_STERIC_DEMAND = 0
CH_ELEC_DEMAND = 1
CH_HBA_DEMAND = 2
CH_HBD_DEMAND = 3
CH_LIPO_DEMAND = 4
CH_AROM_DEMAND = 5
CH_METAL_DEMAND = 6
CH_CATION_DEMAND = 7
CH_ANION_DEMAND = 8
CH_PHOBIC_CORE = 9
CH_SOLVENT_EXPO = 10

# --- C Structure Definitions ---
class NibbleGridStruct(ctypes.Structure):
    _fields_ = [
        ("dim_x", ctypes.c_int),
        ("dim_y", ctypes.c_int),
        ("dim_z", ctypes.c_int),
        ("resolution", ctypes.c_float),
        ("data", ctypes.POINTER(ctypes.c_float))
    ]

# --- NumPy Fallback Implementation ---
class NumPyNibbleEngine:
    """
    A pure NumPy implementation of the Nibble engine logic.
    Provides O(1) matrix operations via NumPy's C-level vectorization.
    Used as a high-performance fallback when nibble.dll cannot be loaded.
    """
    def __init__(self, dim_x=30, dim_y=30, dim_z=30, resolution=0.8):
        self.dim_x = dim_x
        self.dim_y = dim_y
        self.dim_z = dim_z
        self.resolution = resolution
        self.shape = (dim_x, dim_y, dim_z, N_CHANNELS)
        self.data = np.zeros(self.shape, dtype=np.float32)
        self.reset_steric()

    def reset_steric(self):
        self.data[..., CH_STERIC_DEMAND] = 1.0

    def project_atom(self, x, y, z, radius, channels):
        """Gaussian projection into the voxel grid."""
        # Calculate grid coords
        res = self.resolution
        ix, iy, iz = int(x/res), int(y/res), int(z/res)
        rv = int(radius/res) + 1
        
        # Slicing bounds
        x0, x1 = max(0, ix-rv), min(self.dim_x, ix+rv+1)
        y0, y1 = max(0, iy-rv), min(self.dim_y, iy+rv+1)
        z0, z1 = max(0, iz-rv), min(self.dim_z, iz+rv+1)
        
        # Grid coordinates for the sub-patch
        gx, gy, gz = np.meshgrid(
            np.arange(x0, x1) * res,
            np.arange(y0, y1) * res,
            np.arange(z0, z1) * res,
            indexing='ij'
        )
        
        # Distances
        d2 = (gx - x)**2 + (gy - y)**2 + (gz - z)**2
        mask = d2 < (radius * 2)**2
        
        if not np.any(mask):
            return

        # Gaussian weights
        weights = np.exp(-d2[mask] / (2.0 * radius**2))
        
        # Update channels
        for c_idx, val in enumerate(channels):
            if val != 0:
                self.data[x0:x1, y0:y1, z0:z1, c_idx][mask] += val * weights
        
        # Clip to prevent numerical explosion
        self.data[x0:x1, y0:y1, z0:z1, :][mask] = np.clip(
            self.data[x0:x1, y0:y1, z0:z1, :][mask], -200.0, 200.0
        )

    def load_pdb_pocket(self, pdb_path, center, range_val, blur_radius):
        """Parses PDB and projects demands."""
        self.reset_steric()
        atoms_found = 0
        cx, cy, cz = center
        sx, sy, sz = cx - range_val, cy - range_val, cz - range_val
        
        if not os.path.exists(pdb_path):
            return 0
            
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        ax = float(line[30:38])
                        ay = float(line[38:46])
                        az = float(line[46:54])
                        
                        dist2 = (ax-cx)**2 + (ay-cy)**2 + (az-cz)**2
                        if dist2 > range_val**2:
                            continue
                            
                        # B-factor adaptive blur
                        atom_blur = blur_radius
                        try:
                            b_val = float(line[60:66])
                            if b_val > 0:
                                b_sigma = np.sqrt(b_val / 78.957)
                                atom_blur *= (1.0 + min(2.0, b_sigma))
                        except: pass
                        
                        element = line[76:78].strip() or line[12:14].strip()
                        demands = np.zeros(N_CHANNELS)
                        demands[CH_STERIC_DEMAND] = -1.0
                        
                        if 'O' in element:
                            demands[CH_HBA_DEMAND] = 1.0
                            demands[CH_ELEC_DEMAND] = -0.6
                        elif 'N' in element:
                            demands[CH_HBD_DEMAND] = 1.0
                            demands[CH_ELEC_DEMAND] = 0.4
                        elif 'C' in element:
                            demands[CH_LIPO_DEMAND] = 0.7
                        
                        self.project_atom(ax-sx, ay-sy, az-sz, atom_blur, demands)
                        atoms_found += 1
                    except: continue
        return atoms_found

    def compute_affinity(self, other_engine):
        """Calculates binding affinity via dot product."""
        ps = self.data[..., CH_STERIC_DEMAND]
        ds = other_engine.data[..., CH_STERIC_DEMAND]
        
        # Hard clash penalty
        clash_mask = (ds > 0.1) & (ps < 0.05)
        score = -500.0 * np.sum(ds[clash_mask])
        
        # Channel overlap
        # We exclude steric from the sum since it's handled above or by direct mult
        score += np.sum(self.data * other_engine.data)
        return float(score)

    def find_binding_sites(self, n_sites=3, min_distance=5.0):
        """
        Find putative binding sites using inverse vector field analysis.

        The active site is where steric demand is LOW (space exists) but
        all other channel demands are HIGH (chemistry demands interaction).
        No prior knowledge of pocket location needed.

        Inverse score per voxel:
            site_score = (1 - steric) * mean(channels 1-10)

        Args:
            n_sites:      number of candidate sites to return
            min_distance: minimum Angstrom distance between sites (suppresses
                          peaks that are too close — non-maximum suppression)

        Returns:
            List of dicts, each with:
                center:    (x, y, z) in Angstroms
                score:     field-theoretic druggability score
                voxel_idx: (ix, iy, iz) grid indices
                channel_profile: per-channel mean in a 3-voxel radius
        """
        steric = self.data[..., CH_STERIC_DEMAND]          # (x, y, z)
        other  = self.data[..., 1:].mean(axis=-1)          # mean of channels 1-10

        # Normalise both fields to [0, 1]
        def _norm(arr):
            lo, hi = arr.min(), arr.max()
            return (arr - lo) / (hi - lo + 1e-9)

        steric_n = _norm(np.abs(steric))   # abs — steric is often negative
        other_n  = _norm(other)

        # Inverse vector score: empty space WITH chemical demand
        site_score = (1.0 - steric_n) * other_n

        # Smooth to reduce noise
        from scipy.ndimage import gaussian_filter
        site_score = gaussian_filter(site_score, sigma=1.0)

        # Non-maximum suppression — find peaks separated by min_distance
        min_voxels = max(1, int(min_distance / self.resolution))
        sites = []
        score_map = site_score.copy()

        for _ in range(n_sites):
            idx = np.unravel_index(np.argmax(score_map), score_map.shape)
            peak_score = float(score_map[idx])

            if peak_score < 1e-6:
                break

            # Convert voxel index to Angstrom coordinates
            ix, iy, iz = idx
            center = (
                ix * self.resolution,
                iy * self.resolution,
                iz * self.resolution,
            )

            # Per-channel profile in a small sphere around the peak
            r = min(2, min_voxels)
            x0, x1 = max(0, ix-r), min(self.dim_x, ix+r+1)
            y0, y1 = max(0, iy-r), min(self.dim_y, iy+r+1)
            z0, z1 = max(0, iz-r), min(self.dim_z, iz+r+1)
            local = self.data[x0:x1, y0:y1, z0:z1, :]
            channel_profile = {
                'steric':    float(local[..., CH_STERIC_DEMAND].mean()),
                'elec':      float(local[..., CH_ELEC_DEMAND].mean()),
                'hba':       float(local[..., CH_HBA_DEMAND].mean()),
                'hbd':       float(local[..., CH_HBD_DEMAND].mean()),
                'lipo':      float(local[..., CH_LIPO_DEMAND].mean()),
                'arom':      float(local[..., CH_AROM_DEMAND].mean()),
                'metal':     float(local[..., CH_METAL_DEMAND].mean()),
                'cation':    float(local[..., CH_CATION_DEMAND].mean()),
                'anion':     float(local[..., CH_ANION_DEMAND].mean()),
                'phobic':    float(local[..., CH_PHOBIC_CORE].mean()),
                'solvent':   float(local[..., CH_SOLVENT_EXPO].mean()),
            }

            sites.append({
                'center': center,
                'score': peak_score,
                'voxel_idx': idx,
                'channel_profile': channel_profile,
                'rank': len(sites) + 1,
            })

            # Suppress region around this peak for next iteration
            score_map[
                max(0, ix - min_voxels): ix + min_voxels + 1,
                max(0, iy - min_voxels): iy + min_voxels + 1,
                max(0, iz - min_voxels): iz + min_voxels + 1,
            ] = 0.0

        return sites

# --- Native C implementation Bridge ---
class NativeNibbleEngine:
    def __init__(self, lib, dim_x=30, dim_y=30, dim_z=30, resolution=0.8):
        self._lib = lib
        self.dim_x = dim_x
        self.dim_y = dim_y
        self.dim_z = dim_z
        self.resolution = resolution
        self.grid_ptr = self._lib.nibble_create_grid(dim_x, dim_y, dim_z, resolution)

    def __del__(self):
        if hasattr(self, 'grid_ptr') and self.grid_ptr:
            self._lib.nibble_free_grid(self.grid_ptr)

    def load_pdb_pocket(self, pdb_path, center, range_val, blur_radius):
        return self._lib.nibble_load_pdb_pocket(
            self.grid_ptr, pdb_path.encode('utf-8'),
            float(center[0]), float(center[1]), float(center[2]),
            float(range_val), float(blur_radius)
        )

    def project_molecule_data(self, coords, charges, hydrophobicity, start_coords):
        # Create temporary drug grid
        drug_ptr = self._lib.nibble_create_grid(self.dim_x, self.dim_y, self.dim_z, self.resolution)
        for i in range(len(coords)):
            ch_vals = (ctypes.c_float * N_CHANNELS)(0)
            ch_vals[CH_STERIC_DEMAND] = 1.0
            ch_vals[CH_ELEC_DEMAND] = float(charges[i])
            ch_vals[CH_LIPO_DEMAND] = float(hydrophobicity[i])
            
            self._lib.nibble_project_atom(
                drug_ptr,
                float(coords[i][0] - start_coords[0]),
                float(coords[i][1] - start_coords[1]),
                float(coords[i][2] - start_coords[2]),
                1.5, ch_vals
            )
        affinity = self._lib.nibble_compute_affinity(self.grid_ptr, drug_ptr)
        self._lib.nibble_free_grid(drug_ptr)
        return float(affinity)

# --- Main Interface ---
class NibbleEngine:
    """
    Orchestrator for the Nibble docking engine.
    Automatically selects the best available backend (Native C or NumPy).
    """
    _lib = None
    _tried_load = False

    def __init__(self, dim_x=30, dim_y=30, dim_z=30, resolution=0.8):
        self._init_lib()
        self.dim_x = dim_x
        self.dim_y = dim_y
        self.dim_z = dim_z
        self.resolution = resolution
        
        if self._lib:
            self.backend = NativeNibbleEngine(self._lib, dim_x, dim_y, dim_z, resolution)
            self.mode = "C-Native"
        else:
            self.backend = NumPyNibbleEngine(dim_x, dim_y, dim_z, resolution)
            self.mode = "NumPy-Fallback"

    @classmethod
    def _init_lib(cls):
        if cls._tried_load: return
        cls._tried_load = True
        
        lib_name = 'nibble.dll' if os.name == 'nt' else 'libnibble.so'
        lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'nibble', lib_name))
        
        try:
            cls._lib = ctypes.CDLL(lib_path)
            # Setup argtypes
            cls._lib.nibble_create_grid.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_float]
            cls._lib.nibble_create_grid.restype = ctypes.POINTER(NibbleGridStruct)
            cls._lib.nibble_free_grid.argtypes = [ctypes.POINTER(NibbleGridStruct)]
            
            cls._lib.nibble_reset_steric.argtypes = [ctypes.POINTER(NibbleGridStruct)]
            
            cls._lib.nibble_load_pdb_pocket.argtypes = [
                ctypes.POINTER(NibbleGridStruct), ctypes.c_char_p, 
                ctypes.c_float, ctypes.c_float, ctypes.c_float, 
                ctypes.c_float, ctypes.c_float
            ]
            cls._lib.nibble_load_pdb_pocket.restype = ctypes.c_int
            
            cls._lib.nibble_compute_affinity.argtypes = [
                ctypes.POINTER(NibbleGridStruct), ctypes.POINTER(NibbleGridStruct)
            ]
            cls._lib.nibble_compute_affinity.restype = ctypes.c_float
            
            cls._lib.nibble_update_local.argtypes = [
                ctypes.POINTER(NibbleGridStruct), ctypes.c_int, ctypes.c_int, ctypes.c_int,
                ctypes.c_int, ctypes.POINTER(ctypes.c_float)
            ]
            
            cls._lib.nibble_project_atom.argtypes = [
                ctypes.POINTER(NibbleGridStruct), ctypes.c_float, ctypes.c_float, ctypes.c_float,
                ctypes.c_float, ctypes.POINTER(ctypes.c_float)
            ]
        except Exception:
            cls._lib = None

    def load_pdb_pocket(self, pdb_path, center, range_val=10.0, blur_radius=1.5):
        return self.backend.load_pdb_pocket(pdb_path, center, range_val, blur_radius)

    def load_full_protein(self, pdb_path, blur_radius=1.5):
        """
        Load an entire protein PDB without specifying a binding site center.
        Use find_binding_sites() afterwards to detect pockets automatically.

        The grid is sized to fit the protein — dim and resolution set at
        NibbleEngine construction time should be large enough.
        (Tip: use dim=60, resolution=1.0 for a typical ~300 residue protein)
        """
        # Load with a large range centered at origin — covers whole protein
        return self.backend.load_pdb_pocket(
            pdb_path,
            center=(
                self.dim_x * self.backend.resolution / 2,
                self.dim_y * self.backend.resolution / 2,
                self.dim_z * self.backend.resolution / 2,
            ),
            range_val=max(self.dim_x, self.dim_y, self.dim_z) * self.backend.resolution,
            blur_radius=blur_radius,
        )

    def find_binding_sites(self, n_sites=3, min_distance=5.0):
        """
        Detect putative binding sites using inverse vector field analysis.

        No prior knowledge of pocket location required.

        Inverse score = (1 - steric) × mean(other channels)
        Peak voxels of this field = where space exists AND chemistry demands binding.

        Args:
            n_sites:      number of candidate sites (default 3)
            min_distance: minimum Angstrom separation between sites

        Returns:
            List of site dicts sorted by score:
                center:          (x, y, z) Angstroms
                score:           field-theoretic druggability score [0,1]
                voxel_idx:       (ix, iy, iz)
                channel_profile: per-channel mean in pocket neighbourhood
                rank:            1 = best

        Example:
            engine = NibbleEngine(dim_x=60, dim_y=60, dim_z=60, resolution=1.0)
            engine.load_full_protein("PBP2a.pdb")
            sites = engine.find_binding_sites(n_sites=3)
            best = sites[0]
            print(f"Best site: {best['center']}  score={best['score']:.4f}")
        """
        return self.backend.find_binding_sites(n_sites=n_sites,
                                               min_distance=min_distance)

    def project_molecule(self, drug_engine, start_coords):
        """Projects APIs/Excipients and returns affinity."""
        if self.mode == "C-Native":
            # Extract combined data for C
            all_coords = []
            all_charges = []
            all_hydro = []
            for mol in drug_engine.apis + drug_engine.excipients:
                all_coords.extend(mol.coords)
                all_charges.extend(mol.charges)
                all_hydro.extend(mol.hydrophobicity)
            return self.backend.project_molecule_data(all_coords, all_charges, all_hydro, start_coords)
        else:
            # NumPy path
            drug_grid = NumPyNibbleEngine(self.dim_x, self.dim_y, self.dim_z, self.resolution)
            for mol in drug_engine.apis + drug_engine.excipients:
                for i in range(len(mol.coords)):
                    ch_vals = [0.0] * N_CHANNELS
                    ch_vals[CH_STERIC_DEMAND] = 1.0
                    ch_vals[CH_ELEC_DEMAND] = float(mol.charges[i])
                    ch_vals[CH_LIPO_DEMAND] = float(mol.hydrophobicity[i])
                    drug_grid.project_atom(
                        mol.coords[i][0] - start_coords[0],
                        mol.coords[i][1] - start_coords[1],
                        mol.coords[i][2] - start_coords[2],
                        1.5, ch_vals
                    )
            return self.backend.compute_affinity(drug_grid)
