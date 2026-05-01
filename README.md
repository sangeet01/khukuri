# Khukuri Virtual Lab 

<p align="center">
  <img src="assets/logo.png" width="400" alt="Khukuri Virtual Lab Logo">
</p>

AI-powered drug discovery platform with the PINCER counter-evolution engine.

## Features

### Core Modules
- **Target Discovery**: PPI network analysis, target ranking, **LimitNumen literature mining**
- **Molecule Design**: **linearscript grammar-guided generation**, property optimization
- **Molecular Docking**: **KeyBox 11-channel voxel physics**, **Inverse Vector Field site detection** for docking and other applications
- **ADMET Prediction**: Drug-likeness, toxicity, pharmacokinetics
- **Resistance Prediction**: Multi-target strategies, evolution simulation
- **Retrosynthesis**: Route planning, synthetic accessibility scoring

### PINCER Counter-Evolution Engine
The Pincer algorithm models the pharmacological duel between antibiotics and bacterial resistance as a zero-sum minimax game:

- **Dual-Threat Red Team**: Maps both vertical **Markov Mutation Matrices** (2-bit DNA encoding with transition bias) and **Horizontal Gene Transfer (HGT)** resistance networks.
- **Mutation Space Mapper**: Enumerates all viable receptor mutations within the finite biological hard ceiling.
- **Darwin-Godel Loop (Blue Team)**: Evolves drug candidates via grammar-safe mutations, optimizing worst-case binding across the full mutation cluster.
- **ThreatIndex & Dead Zone**: O(1) LimitNumen vector retrieval for fast threat evaluation and content-addressable space tracking.
- **Minimax Fitness**: Drugs are scored by their weakest binding across all viable mutations -- the **Skeleton Key** objective.

### AMR-Specific Features
- **Bio-Knowledge Layer**: CARD/ResFinder/MEGARes integration, pathogen database, target proteins
- **Genomics Analysis**: Resistance mutations, evolution prediction, multi-omics integration
- **Microbiology Tools**: MIC analysis, assay tracking, strain management
- **Persistent World Model**: Accumulates intelligence (compounds, targets, hypotheses) seamlessly across sessions.
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
|   |-- integrations/      # KeyBox Bridge, LimitNumen Index
|   |-- world_model/       # Kosmos engine, knowledge graph
|   +-- workflows/         # End-to-end pipelines (AMR discovery)
|-- tests/                 # Test suite
|-- config/                # Configuration files
|-- scripts/               # Automation scripts
+-- examples/              # Usage examples
```

## Usage Examples

### PINCER Counter-Evolution

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

### AMR Discovery (with KeyBox, Numen & PINCER)

```python
from src.workflows import AMRDiscoveryWorkflow

# Run AMR-focused discovery -- KeyBox, Numen, and PINCER run automatically
workflow = AMRDiscoveryWorkflow(
    pdb_path="PBP2a_2OLV.pdb", # KeyBox auto-detects the binding site via field physics
    project="mrsa_run_01"      # Persistent world model remembers past runs
)

results = workflow.run_discovery(
    pathogen="Staphylococcus aureus MRSA",
    priority="critical",
    n_compounds=20,
    n_iterations=3
)

print(f"Targets Identified: {len(results['targets'])}")
print(f"Numen Literature Targets: {[t['name'] for t in results['targets'] if t.get('source') == 'numen_literature']}")
print(f"Apex drug: {results['iterations'][-1]['top_compounds'][0]['smiles']}")
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
- linearscript >= 3.0.1
- NumPy, Pandas, SciPy
- PyYAML >= 6.0
- Requests >= 2.26.0
- python-louvain >= 0.16
- OpenAI >= 1.0.0 (optional)

## Architecture

The PINCER algorithm operates as the evolutionary core of the Khukuri platform:

```text
Literature Mining (Numen) ──→ Target Discovery
         ↓
    Red Team (PINCER)                      Solution Side (Blue Team):
      |-- Vertical Mutation Matrix           Seed Drugs (Numen CompoundMemory)
      |-- Horizontal Gene Transfer Mapper      |
         ↓                                     v
    Threat Index (Numen) ──→ Fast retrieval ──→ Darwin-Godel Evolutionary Loop
         ↓                                     (SCRIPT grammar-safe mutations)
    Viable Threat Cluster                      |
         |                                     v
         |                               Candidate Population
         |                                     |
         +------------------------------------ v
                                  Minimax Fitness Evaluation
                                  (KeyBox Voxel Physics Docking)
                                               |
                                               v
     Dead Zone Check <-------------- Apex Skeleton Key s*
                                               |
                                               v
                                   World Model Persistence
                                  (ADMET / Synthesis / MATP)
                                               |
                                               v
                                     Optimized Lab Candidate
```

## Status

✅ **Research-grade**: 70+ files, ~3,500 lines of code, 0% mocks, 0% duplication


#

Khukuri V3.0
