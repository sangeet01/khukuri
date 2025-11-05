# User Guide

## Installation

```bash
cd khukuri
pip install -r requirements.txt
python scripts/quick_test.py
```

## Quick Start

### Target Discovery
```python
from src.target_discovery import NetworkAnalyzer, TargetRanker, DataFetcher

fetcher = DataFetcher()
ppi_data = fetcher.fetch_string_ppi(["inhA", "katG"], species=83332)

analyzer = NetworkAnalyzer()
network = analyzer.build_ppi_network(ppi_data)

ranker = TargetRanker()
targets = ranker.rank_targets(network, ["inhA", "katG"])
```

### Molecule Generation
```python
from src.molecule_design import MoleculeGenerator
from rdkit import Chem

generator = MoleculeGenerator()
scaffold = Chem.MolFromSmiles("c1ccccc1")
molecules = generator.generate_library(scaffold, num_molecules=100)
```

### ADMET Prediction
```python
from src.admet import predict_admet
from rdkit import Chem

mol = Chem.MolFromSmiles("CCO")
admet = predict_admet(mol)
print(f"Drug-likeness: {admet['qed']:.3f}")
```

### Molecular Docking
```python
from src.docking import ReceptorPrep, VinaWrapper

prep = ReceptorPrep()
receptor = prep.download_pdb("1HSG")

vina = VinaWrapper()
results = vina.run_docking(receptor, mol, center=(15,15,15), box_size=(20,20,20))
```

### Autonomous Workflow
```python
from src.workflows import run_autonomous_discovery

results = run_autonomous_discovery(
    disease="tuberculosis",
    target_genes=["inhA", "katG"],
    num_candidates=10
)
```

## Configuration

Edit `config/default_config.yaml`:
```yaml
scoring:
  weights:
    binding_affinity: 0.4
    drug_likeness: 0.25
```

## Examples

See `examples/` directory for complete working examples.
