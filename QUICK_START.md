# Khukuri v3.0 - Quick Start Guide

<p align="center">
  <img src="assets/logo.png" width="300" alt="Khukuri Virtual Lab Logo">
</p>

Khukuri is a production-grade, autonomous drug discovery platform designed for Antimicrobial Resistance (AMR). It integrates multi-agent orchestration, physics-based docking, and the PINCER counter-evolution engine.

---

## 🚀 What's New in v3.0 (The PINCER Update)

Khukuri now features the **PINCER Counter-Evolution Engine**, a Darwin-Gödel co-evolutionary system for preempting antimicrobial resistance:

- **Darwin-Gödel Engine**: Fully autonomous in-silico drug discovery pipeline that optimizes against predicted future mutations.
- **Minimax Pharmacological Duel**: Models resistance as a zero-sum game ($s^* = \arg\max \min K$), ensuring drugs are effective against the worst-case bacterial threats.
- **2-Bit DNA Hardware Logic**: Optimized binary encoding for Markov-biased mutation prediction and $O(1)$ complementation.
- **Mutation Space Mapping**: Pre-computes the finite cluster of viable bacterial threats (the "Red Team").
- **Skeleton Key Objective**: Evolving molecules that maintain binding affinity across an entire cluster of predicted mutants.

---

## 📦 Installation

```bash
# Clone the repository
git clone https://github.com/your-repo/khukuri.git
cd khukuri

# Install dependencies
pip install -r requirements.txt
```

---

## 🔍 Quick Test (30 seconds)

Verify your entire installation (Core, ADMET, Docking, PINCER, Synthesis, and Agents) using the automated smoke test:

```bash
python scripts/quick_test.py
```

Expected output:
```
[OK] Core modules
[OK] ADMET modules
[OK] Target discovery modules
[OK] Molecule design modules
[OK] Docking modules
[OK] Resistance modules (including PINCER)
...
[SUCCESS] All tests passed! Khukuri is ready to use.
```

---

## 🛠️ Core Modules & Usage

### 1. PINCER Counter-Evolution (New)
Pre-empt resistance by evolving "Skeleton Key" molecules that target the entire viable mutation space of a pocket.

```python
from src.resistance import PincerEngine

# Initialize the engine
pincer = PincerEngine(population_size=100, n_generations=50)

# Phase 1: Map the viable threat landscape (Red Team)
threats = pincer.map_threats("AMILVCFYWHDEKRSTGNPQ", list(range(20)))

# Phase 2: Evolve the Skeleton Key drug (Blue Team)
seed_drugs = ["c1ccccc1", "c1ccncc1", "c1ccoc1"]
apex = pincer.evolve(seed_drugs)
print(f"Apex Skeleton Key: {apex.smiles} (Minimax Score: {apex.minimax_score:.4f})")
```

### 2. Autonomous AMR Discovery (The Macro Loop)
The primary entry point for drug discovery. This workflow automates target selection, compound generation, and PINCER analysis.

```python
from src.workflows import run_amr_discovery

# Discover compounds for M. tuberculosis with PINCER pre-computation
results = run_amr_discovery(
    pathogen='Mycobacterium tuberculosis',
    priority='critical',
    n_compounds=20,
    n_iterations=3
)

print(f"PINCER threats mapped: {results['pincer_analysis']['viable_threats_count']}")
print(f"Top recommendation: {results['recommendations'][0]['compound_id']}")
```

### 3. Molecular Docking & Binding Site Detection
Physics-based validation of drug candidates against bacterial targets using AutoDock Vina.

```python
from src.docking import BindingSiteDetector, VinaWrapper

# Detect pockets in a protein structure
detector = BindingSiteDetector()
pockets = detector.detect_pockets("target.pdb")

# Run docking
vina = VinaWrapper()
result = vina.dock("ligand.sdf", "target.pdbqt", pockets[0]['center'])
print(f"Binding Affinity: {result['affinity']} kcal/mol")
```

### 4. ADMET & Molecular Design
AI-driven molecule generation, fragment-based design, and property optimization.

```python
from src.molecule_design import MoleculeGenerator
from src.admet import predict_admet

# Generate candidates
generator = MoleculeGenerator()
candidates = generator.generate(n_molecules=100)

# Filter by drug-likeness
for smiles in candidates:
    props = predict_admet(smiles)
    if props['drug_likeness'] > 0.8:
        print(f"High-quality candidate: {smiles}")
```

### 5. World Model & Knowledge Graph
Persistent state tracking (Kosmos Engine) and Bayesian learning loop.

```python
from src.world_model import WorldStateTracker, KnowledgeGraph

world = WorldStateTracker()
kg = KnowledgeGraph()

# Add discovery observations
kg.add_binding('COMP_001', 'InhA', affinity=-8.5)
world.update_compound('COMP_001', {'smiles': 'CCO', 'activity': 'active'})
world.add_hypothesis('Compound binds InhA', ['docking'], 0.85)
```

### 6. Resistance Analysis (Genomics & Microbiology)
Analyze clinical resistance mutations and track MIC (Minimum Inhibitory Concentration) data.

```python
from src.bioknowledge import ResistanceDatabase
from src.genomics import ResistanceGenomicsAnalyzer
from src.microbiology import MICAnalyzer

# Query resistance genes
db = ResistanceDatabase()
genes = db.get_genes_by_organism('S. aureus')

# Analyze mutations (e.g. rpoB S531L)
analyzer = ResistanceGenomicsAnalyzer()
profile = analyzer.analyze_mutation_profile('rpoB', ['S531L'])

# Track MIC data
mic = MICAnalyzer()
mic.add_mic_result('COMP_001', 'ATCC_25923', 'rifampicin', 0.5)
```

---

## 🧪 Run Examples & Tests

### Examples
| Example | Description |
| :--- | :--- |
| `python examples/amr_discovery_example.py` | Full end-to-end AMR discovery workflow |
| `python examples/pincer_example.py` | Standalone demonstration of the PINCER engine |
| `python examples/world_model_example.py` | Knowledge graph and state tracking |

### Tests
```bash
# Run the full test suite (Core, ADMET, Genomics, etc.)
python -m pytest tests/ -v

# Run the specific PINCER counter-evolution tests
make test-pincer
```

---

## ⚙️ Configuration

Customize the discovery pipeline by editing `config/amr_config.yaml`:

```yaml
# PINCER Engine Parameters
pincer:
  population_size: 100
  n_generations: 50
  transition_weight: 0.7  # Markov transition bias
  fitness_threshold: 0.3  # Biological hard ceiling for viable mutants

# Discovery Workflow Parameters
amr_discovery:
  default_pathogen: "Mycobacterium tuberculosis"
  n_compounds_per_iteration: 20
  n_iterations: 3
```

---

## 🏗️ Architecture

```
Khukuri Virtual Lab
├── workflows/         # End-to-end pipelines (AMR Discovery)
├── resistance/        # PINCER Engine, Darwin-Godel Loop
├── docking/           # Vina Wrapper, Pocket Detection
├── molecule_design/   # Generation & Optimization
├── world_model/       # Kosmos Engine, Knowledge Graph
├── bioknowledge/      # Pathogen & Resistance Databases
└── agents/            # Multi-agent Orchestration
```

---

## 📜 Version & Status

**Version**: 3.0.0 (PINCER Counter-Evolution Update)
**Status**: Production-ready ✓
**Author**: Khukuri Research Team
