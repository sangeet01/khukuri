# API Reference

## Core Module

### MolecularScorer
```python
from src.core import MolecularScorer
scorer = MolecularScorer()
score, metrics, breakdown = scorer.calculate_composite_score(mol)
```

### Validator
```python
from src.core import Validator
validator = Validator()
is_valid = validator.validate_molecule(mol)
```

## Target Discovery

### NetworkAnalyzer
```python
from src.target_discovery import NetworkAnalyzer
analyzer = NetworkAnalyzer()
network = analyzer.build_ppi_network(ppi_data)
```

### TargetRanker
```python
from src.target_discovery import TargetRanker
ranker = TargetRanker()
targets = ranker.rank_targets(network, seed_genes)
```

## Molecule Design

### MoleculeGenerator
```python
from src.molecule_design import MoleculeGenerator
generator = MoleculeGenerator()
molecules = generator.generate_library(scaffold, num_molecules=100)
```

### PropertyOptimizer
```python
from src.molecule_design import PropertyOptimizer
optimizer = PropertyOptimizer()
optimized = optimizer.optimize_admet(mol, target_properties)
```

## Docking

### VinaWrapper
```python
from src.docking import VinaWrapper
vina = VinaWrapper()
results = vina.run_docking(receptor, ligand, center, box_size)
```

## ADMET

```python
from src.admet import calculate_drug_likeness, predict_admet
drug_likeness = calculate_drug_likeness(mol)
admet_profile = predict_admet(mol)
```

## Workflows

```python
from src.workflows import run_autonomous_discovery
results = run_autonomous_discovery(
    disease="tuberculosis",
    target_genes=["inhA", "katG"],
    num_candidates=10
)
```
