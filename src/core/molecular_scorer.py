"""Unified molecular scoring system"""

import logging
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, QED
import yaml
from pathlib import Path

logger = logging.getLogger('khukuri')


class MolecularScorer:
    """Unified system for molecular scoring and analysis"""
    
    def __init__(self, config_path=None):
        """Initialize with optional config file"""
        self.config = self._load_config(config_path)
    
    def _load_config(self, config_path):
        """Load configuration from YAML or use defaults"""
        default_config = {
            'weights': {
                'binding_affinity': 0.4,
                'drug_likeness': 0.25,
                'synthetic_accessibility': 0.15,
                'toxicity_risk': 0.1,
                'pharmacokinetic_profile': 0.1
            },
            'thresholds': {
                'binding_affinity_cutoff': -5.0,
                'lipinski_violations_max': 1,
                'qed_min': 0.5,
                'sa_max': 6.0,
                'toxicity_alerts_max': 2
            }
        }
        
        if config_path and Path(config_path).exists():
            with open(config_path) as f:
                return yaml.safe_load(f)
        return default_config
    
    def calculate_composite_score(self, mol, binding_affinity=None):
        """Calculate comprehensive composite score"""
        try:
            metrics = {
                'binding_affinity': binding_affinity if binding_affinity else 0,
                'drug_likeness': self.calculate_drug_likeness(mol),
                'synthetic_accessibility': self.calculate_sa_score(mol),
                'toxicity_risk': self.assess_toxicity_risk(mol),
                'pharmacokinetic_profile': self.evaluate_pharmacokinetics(mol)
            }
            
            normalized = self._normalize_metrics(metrics)
            
            composite_score = sum(
                normalized[metric] * self.config['weights'][metric]
                for metric in self.config['weights']
            )
            
            return composite_score, metrics, normalized
            
        except Exception as e:
            logger.error(f"Composite scoring failed: {e}")
            return 0.0, {}, {}
    
    def calculate_drug_likeness(self, mol):
        """Calculate drug-likeness score"""
        try:
            lipinski_violations = self._lipinski_violations(mol)
            qed_score = QED.qed(mol)
            
            mw = Descriptors.MolWt(mol)
            mw_score = max(0, 1 - abs(mw - 350) / 200)
            
            logp = Descriptors.MolLogP(mol)
            logp_score = max(0, 1 - abs(logp - 2.5) / 3)
            
            drug_likeness = (
                (1 - lipinski_violations / 5) * 0.4 +
                qed_score * 0.3 +
                mw_score * 0.15 +
                logp_score * 0.15
            )
            
            return max(0, min(1, drug_likeness))
        except Exception as e:
            logger.error(f"Drug-likeness calculation failed: {e}")
            return 0.0
    
    def calculate_sa_score(self, mol):
        """Calculate synthetic accessibility score (1-10, lower is better)"""
        try:
            from rdkit.Chem import rdMolDescriptors
            sa_score = rdMolDescriptors.BertzCT(mol) / 100
            
            rings = mol.GetRingInfo().NumRings()
            heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])
            
            complexity_penalty = rings * 0.2 + heteroatoms * 0.1
            final_sa = min(10, sa_score + complexity_penalty)
            
            return max(1, final_sa)
        except:
            rings = mol.GetRingInfo().NumRings()
            heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])
            return min(10, 1 + rings * 0.5 + heteroatoms * 0.3)
    
    def assess_toxicity_risk(self, mol):
        """Assess toxicity risk (0-1, higher is riskier)"""
        try:
            alerts = self._check_toxicity_alerts(mol)
            return min(1.0, alerts / 5.0)
        except:
            return 0.5
    
    def evaluate_pharmacokinetics(self, mol):
        """Evaluate pharmacokinetic profile (0-1, higher is better)"""
        try:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            hbd = Descriptors.NumHDonors(mol)
            rotbonds = Descriptors.NumRotatableBonds(mol)
            
            pk_score = 1.0
            if mw > 500: pk_score -= 0.2
            if logp < 0 or logp > 5: pk_score -= 0.2
            if tpsa > 140: pk_score -= 0.3
            if hbd > 5: pk_score -= 0.2
            if rotbonds > 10: pk_score -= 0.1
            
            return max(0, pk_score)
        except:
            return 0.5
    
    def _lipinski_violations(self, mol):
        """Calculate Lipinski Rule of 5 violations"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if hbd > 5: violations += 1
        if hba > 10: violations += 1
        return violations
    
    def _check_toxicity_alerts(self, mol):
        """Check for structural toxicity alerts"""
        alerts = 0
        toxic_patterns = [
            '[N+](=O)[O-]',  # Nitro groups
            'C=C-C=O',       # α,β-unsaturated carbonyl
            '[Cl,Br,I]C=C',  # Halogenated alkenes
            'N-N',           # Hydrazines
            '[S,P](=O)(=O)'  # Sulfonyl/phosphonyl
        ]
        for pattern in toxic_patterns:
            try:
                if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                    alerts += 1
            except:
                continue
        return alerts
    
    def _normalize_metrics(self, metrics):
        """Normalize metrics to 0-1 scale"""
        normalized = {}
        
        # Binding affinity: more negative is better
        if metrics['binding_affinity'] is not None:
            affinity_score = max(0, min(1, (-metrics['binding_affinity'] + 15) / 15))
        else:
            affinity_score = 0
        normalized['binding_affinity'] = affinity_score
        
        normalized['drug_likeness'] = metrics['drug_likeness']
        normalized['synthetic_accessibility'] = 1 - (metrics['synthetic_accessibility'] - 1) / 9
        normalized['toxicity_risk'] = 1 - metrics['toxicity_risk']
        normalized['pharmacokinetic_profile'] = metrics['pharmacokinetic_profile']
        
        return normalized
