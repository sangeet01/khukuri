# Khukuri Virtual Lab

Production-ready AI-powered drug discovery platform with modular architecture.

## Features

- **Target Discovery**: PPI network analysis, target ranking, literature mining
- **Molecule Design**: AI-powered generation, property optimization, fragment-based design
- **Molecular Docking**: AutoDock Vina integration, binding site detection, pose analysis
- **ADMET Prediction**: Drug-likeness, toxicity, pharmacokinetics
- **Resistance Prediction**: Multi-target strategies, evolution simulation
- **Retrosynthesis**: Route planning, synthetic accessibility scoring
- **Multi-Agent System**: AI agents for autonomous drug discovery workflows

## Quick Start

### Installation

```bash
# Install dependencies
pip install -r requirements.txt

# For development
pip install -r requirements-dev.txt

# Setup AutoDock Vina (Unix/Linux/macOS)
bash scripts/setup_vina.sh
```

### Quick Test

```bash
python scripts/quick_test.py
```

### Run Tests

```bash
python -m pytest tests/ -v
```

## Project Structure

```
khukuri/
├── src/                    # Source code (9 modules, 39 files)
│   ├── core/              # Logging, validation, scoring
│   ├── target_discovery/  # Network analysis, target ranking
│   ├── molecule_design/   # Generation, optimization
│   ├── docking/           # Vina wrapper, pose analysis
│   ├── admet/             # Properties, toxicity, PK/PD
│   ├── resistance/        # Prediction, multi-target
│   ├── synthesis/         # Retrosynthesis, SA scoring
│   ├── agents/            # AI agents, orchestrator
│   └── workflows/         # End-to-end pipelines
├── tests/                 # Test suite (13 files)
├── config/                # Configuration files
├── scripts/               # Automation scripts
└── examples/              # Usage examples
```

## Usage Example

```python
from src.workflows import run_autonomous_discovery

# Run autonomous drug discovery
results = run_autonomous_discovery(
    disease="tuberculosis",
    target_genes=["inhA", "katG"],
    num_candidates=10
)
```

## Dependencies

- RDKit >= 2022.09.1
- NetworkX >= 2.6.0
- BioPython >= 1.79
- NumPy, Pandas, SciPy
- PyYAML >= 6.0
- Requests >= 2.26.0
- python-louvain >= 0.16
- OpenAI >= 1.0.0 (optional)


Authors: Sangeet Sharma and Bipin Gyawali


