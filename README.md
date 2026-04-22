# Khukuri Virtual Lab v3.0

<p align="center">
  <img src="assets/logo.png" width="400" alt="Khukuri Virtual Lab Logo">
</p>

Production-ready AI-powered drug discovery platform with the PINCER counter-evolution engine.

## Features

### Core Modules
- **Target Discovery**: PPI network analysis, target ranking, literature mining
- **Molecule Design**: AI-powered generation, property optimization, fragment-based design
- **Molecular Docking**: AutoDock Vina integration, binding site detection, pose analysis
- **ADMET Prediction**: Drug-likeness, toxicity, pharmacokinetics
- **Resistance Prediction**: Multi-target strategies, evolution simulation
- **Retrosynthesis**: Route planning, synthetic accessibility scoring

### PINCER Counter-Evolution Engine (NEW)
The PINCER algorithm models the pharmacological duel between antibiotics and
bacterial resistance as a zero-sum minimax game:

- **Markov Mutation Matrix**: 2-bit DNA encoding with transition/transversion bias
  for biologically accurate mutation prediction
- **Mutation Space Mapper (Red Team)**: Enumerates all viable receptor mutations
  within the finite biological hard ceiling
- **Darwin-Godel Loop (Blue Team)**: Evolves drug candidates via constrained
  evolution, optimizing worst-case binding across the full mutation cluster
- **Dead Zone Index**: O(1) content-addressable tracking of explored chemical
  space to prevent redundant computation
- **Minimax Fitness**: Drugs are scored by their weakest binding across all
  viable mutations -- the Skeleton Key objective

### AMR-Specific Features
- **Bio-Knowledge Layer**: CARD/ResFinder/MEGARes integration, pathogen database, target proteins
- **Genomics Analysis**: Resistance mutations, evolution prediction, multi-omics integration
- **Microbiology Tools**: MIC analysis, assay tracking, strain management
- **World Model**: Kosmos-style state tracking, knowledge graph, learning loop
- **Specialized Agents**: Microbiology, Genomics, Cheminformatics, Resistance Critic, Literature

## Quick Start

### Installation

```bash
# Install dependencies
pip install -r requirements.txt

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
|-- src/                    # Source code
|   |-- core/              # Logging, validation, scoring
|   |-- target_discovery/  # Network analysis, target ranking
|   |-- molecule_design/   # Generation, optimization
|   |-- docking/           # Vina wrapper, pose analysis
|   |-- admet/             # Properties, toxicity, PK/PD
|   |-- resistance/        # PINCER engine, prediction, multi-target
|   |   |-- pincer_engine.py      # Core PINCER algorithm
|   |   |-- evolution_simulator.py # PINCER-powered evolution sim
|   |   |-- predictor.py          # Resistance likelihood
|   |   +-- multi_target.py       # Combination strategies
|   |-- synthesis/         # Retrosynthesis, SA scoring
|   |-- agents/            # AI agents, orchestrator
|   |-- world_model/       # Kosmos engine, knowledge graph
|   +-- workflows/         # End-to-end pipelines (AMR discovery)
|-- tests/                 # Test suite
|-- config/                # Configuration files
|-- scripts/               # Automation scripts
+-- examples/              # Usage examples
```

## Usage Examples

### PINCER Counter-Evolution (NEW)

```python
from src.resistance import PincerEngine

# Initialize the engine
pincer = PincerEngine(population_size=100, n_generations=50)

# Phase 1 (Red Team): Map viable mutation space
wild_type = "AMILVCFYWHDEKRSTGNPQ"
active_site = list(range(len(wild_type)))
threats = pincer.map_threats(wild_type, active_site)
print(f"Mapped {len(threats)} viable resistance pathways")

# Phase 2 (Blue Team): Evolve Skeleton Key drug
seed_drugs = ["c1ccccc1", "c1ccncc1", "c1ccoc1"]
apex = pincer.evolve(seed_drugs)
print(f"Skeleton Key: {apex.smiles} (minimax={apex.minimax_score:.4f})")

# Full results
results = pincer.get_results()
print(f"Top 5 candidates: {[c['smiles'] for c in results['top_5']]}")
```

### AMR Discovery (with PINCER)

```python
from src.workflows import run_amr_discovery

# Run AMR-focused discovery -- PINCER runs automatically
results = run_amr_discovery(
    pathogen="Mycobacterium tuberculosis",
    priority="critical",
    n_compounds=20,
    n_iterations=3
)

print(f"Targets: {results['targets']}")
print(f"PINCER threats mapped: {results['pincer_analysis']['viable_threats_count']}")
print(f"Apex drug: {results['pincer_analysis']['apex_drug']}")
print(f"Recommendations: {results['recommendations']}")
```

### Traditional Discovery

```python
from src.workflows import run_autonomous_discovery

# Run autonomous drug discovery
results = run_autonomous_discovery(
    disease="tuberculosis",
    target_genes=["inhA", "katG"],
    num_candidates=10
)
```

### Resistance Analysis

```python
from src.bioknowledge import ResistanceDatabase
from src.genomics import ResistanceGenomicsAnalyzer
from src.microbiology import MICAnalyzer

# Query resistance genes
resistance_db = ResistanceDatabase()
genes = resistance_db.get_genes_by_organism('S. aureus')

# Analyze mutations
analyzer = ResistanceGenomicsAnalyzer()
profile = analyzer.analyze_mutation_profile('rpoB', ['S531L'])

# Track MIC data
mic_analyzer = MICAnalyzer()
mic_analyzer.add_mic_result('COMP_001', 'ATCC_25923', 'rifampicin', 0.5)
```

### Low-Level PINCER Components

```python
from src.resistance import (
    MarkovMutationMatrix,
    MutationSpaceMapper,
    encode_sequence,
    bitwise_complement,
    decode_sequence,
)

# 2-bit DNA encoding
bits = encode_sequence("ATCGATCG")
complement = bitwise_complement(bits)
print(decode_sequence(complement))  # TAGCTAGC

# Markov mutation matrix
markov = MarkovMutationMatrix()
print(f"Mutations per generation: {markov.mutations_per_generation:.4f}")
mutated = markov.mutate_bits(bits, n_mutations=2)
print(f"Mutated: {decode_sequence(mutated)}")

# Map binding pocket mutations
mapper = MutationSpaceMapper(fitness_threshold=0.3)
viable = mapper.map_binding_pocket(
    wild_type_seq="AMILVCFYWH",
    active_site_indices=[0, 1, 2, 3, 4],
)
print(f"Viable mutants: {len(viable)}")
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

## Architecture

The PINCER algorithm operates as the evolutionary core of the Khukuri platform:

```
Target Side (Red Team):                    Solution Side (Blue Team):
  Wild-Type Receptor DNA                     Seed Drug Molecules (SMILES)
    |                                          |
    v                                          v
  Markov Mutation Engine                     Darwin-Godel Evolutionary Loop
  (2-bit, transition bias)                   (grammar-safe molecular mutations)
    |                                          |
    v                                          v
  Predicted Mutation Cluster                 Candidate Population
    |                                          |
    v                                          v
  Folding Energy Filter                      Minimax Fitness Evaluation
  (lethality cliff)                          (worst-case binding across threats)
    |                                          |
    v                                          v
  Viable Threat Cluster -----> Godel Proof Validator <----- Dead Zone Check
    |                              |
    v                              v
  LimitNumen Vector Index      Apex Skeleton Key s*
                                   |
                                   v
                               KHUKURI Orchestrator
                               (ADMET / Synthesis / MATP)
                                   |
                                   v
                               Optimized Lab Candidate
```

## Status

✅ **Production-ready**: 70+ files, ~3,500 lines of code, 0% mocks, 0% duplication


