# AMR-Focused Features

## Overview

Khukuri v2.0 introduces production-grade AMR (Antimicrobial Resistance) discovery capabilities, integrating bio-knowledge, genomics, microbiology, and AI agents for autonomous antibiotic discovery.

## Architecture

```
┌─────────────────────────────────────────────────────────┐
│                   AMR Discovery System                   │
├─────────────────────────────────────────────────────────┤
│                                                          │
│  ┌──────────────────┐  ┌──────────────────┐            │
│  │  Bio-Knowledge   │  │    Genomics      │            │
│  │  - CARD/ResFinder│  │  - Mutations     │            │
│  │  - Pathogens     │  │  - Evolution     │            │
│  │  - Targets       │  │  - Omics         │            │
│  └──────────────────┘  └──────────────────┘            │
│                                                          │
│  ┌──────────────────┐  ┌──────────────────┐            │
│  │  Microbiology    │  │   World Model    │            │
│  │  - MIC Analysis  │  │  - State Tracker │            │
│  │  - Assays        │  │  - Knowledge Graph│           │
│  │  - Strains       │  │  - Learning Loop │            │
│  └──────────────────┘  └──────────────────┘            │
│                                                          │
│  ┌──────────────────────────────────────────┐          │
│  │           Multi-Agent System              │          │
│  │  Microbiology | Genomics | Cheminformatics│         │
│  │  Resistance Critic | Literature Analyst   │         │
│  └──────────────────────────────────────────┘          │
│                                                          │
└─────────────────────────────────────────────────────────┘
```

## Key Features

### 1. Bio-Knowledge Layer

**ResistanceDatabase**
- Curated AMR gene database (CARD, ResFinder, MEGARes)
- Resistance mechanisms mapping
- Drug class associations

**PathogenDatabase**
- WHO priority pathogens
- Essential genes and targets
- Resistance profiles

**TargetProteinDB**
- Validated antibiotic targets
- Druggability scores
- PDB structure IDs

### 2. Genomics Layer

**ResistanceGenomicsAnalyzer**
- Mutation profile analysis
- Resistance evolution prediction
- Compensatory mutation identification

**MutationTracker**
- Strain-level mutation tracking
- Co-occurring mutation analysis
- Strain comparison

**OmicsIntegrator**
- Multi-omics data integration
- Genomics + transcriptomics + metabolomics
- Correlation analysis

### 3. Microbiology Layer

**MICAnalyzer**
- MIC data management
- CLSI/EUCAST breakpoint interpretation
- Cross-resistance prediction

**AssayTracker**
- Biological assay tracking
- Activity summarization
- Temporal analysis

**StrainManager**
- Bacterial strain collection
- Reference strain management
- Resistance profiling

### 4. World Model (Kosmos-style)

**WorldStateTracker**
- Complete state tracking (compounds, targets, strains, assays)
- Hypothesis management
- Experiment history
- State export/import

**KnowledgeGraph**
- Entity relationship mapping
- Multi-target compound identification
- Path finding
- Subgraph extraction

**LearningLoop**
- Bayesian optimization
- Active learning
- Acquisition function (UCB)
- Model performance tracking

### 5. Specialized Agents

**MicrobiologyAgent**
- Pathogen biology expertise
- AMR mechanism analysis
- MIC interpretation

**GenomicsAgent**
- Resistance genomics
- Mutation analysis
- Evolutionary patterns

**CheminformaticsAgent**
- Molecular design
- QSAR analysis
- ADMET optimization

**ResistanceCriticAgent**
- Resistance risk assessment
- Multi-target strategy recommendations
- Critical evaluation

**LiteratureAgent**
- Evidence synthesis
- Hypothesis generation
- Knowledge gap identification

## Usage

### Basic AMR Discovery

```python
from src.workflows import run_amr_discovery

results = run_amr_discovery(
    pathogen='Mycobacterium tuberculosis',
    priority='critical',
    n_compounds=20,
    n_iterations=3
)

print(f"Targets: {results['targets']}")
print(f"Recommendations: {results['recommendations']}")
```

### Advanced Workflow

```python
from src.workflows import AMRDiscoveryWorkflow

workflow = AMRDiscoveryWorkflow(openai_client=client)

results = workflow.run_discovery(
    pathogen='Pseudomonas aeruginosa',
    priority='critical',
    n_compounds=50,
    n_iterations=5
)

# Access world state
summary = workflow.world_state.get_state_summary()
workflow.world_state.export_state('state.json')

# Query knowledge graph
multi_target = workflow.knowledge_graph.find_multi_target_compounds()
```

### Resistance Analysis

```python
from src.bioknowledge import ResistanceDatabase, PathogenDatabase
from src.genomics import ResistanceGenomicsAnalyzer

# Query resistance genes
resistance_db = ResistanceDatabase()
genes = resistance_db.get_genes_by_organism('S. aureus')

# Analyze mutations
analyzer = ResistanceGenomicsAnalyzer()
profile = analyzer.analyze_mutation_profile('rpoB', ['S531L', 'H526Y'])

# Predict evolution
prediction = analyzer.predict_resistance_evolution(
    current_mutations={'rpoB': ['S531L']},
    drug_pressure='rifampicin'
)
```

### MIC Tracking

```python
from src.microbiology import MICAnalyzer

mic_analyzer = MICAnalyzer()

# Add results
mic_analyzer.add_mic_result('COMP_001', 'ATCC_25923', 'rifampicin', 0.5)

# Get profile
profile = mic_analyzer.get_mic_profile('ATCC_25923')

# Compare compounds
comparison = mic_analyzer.compare_compounds(
    ['COMP_001', 'COMP_002', 'COMP_003'],
    'ATCC_25923'
)
```

### World Model

```python
from src.world_model import WorldStateTracker, KnowledgeGraph, LearningLoop

# Track state
world_state = WorldStateTracker()
world_state.update_compound('COMP_001', {'activity': 'active'})
world_state.add_hypothesis('Compound targets InhA', ['docking', 'assay'], 0.85)

# Build knowledge graph
kg = KnowledgeGraph()
kg.add_compound('COMP_001', {'mw': 250})
kg.add_target('InhA', {'druggability': 0.9})
kg.add_binding('COMP_001', 'InhA', affinity=-8.5)

# Learning loop
learning = LearningLoop()
learning.add_observation('COMP_001', features, outcome)
next_experiments = learning.prioritize_next_experiments(candidates, n_select=5)
```

## Configuration

Edit `config/amr_config.yaml`:

```yaml
amr_discovery:
  default_pathogen: "Mycobacterium tuberculosis"
  n_compounds_per_iteration: 20
  n_iterations: 3

agents:
  use_ai: true  # Enable with OpenAI API key
  model: "gpt-3.5-turbo"

learning_loop:
  acquisition_function: "UCB"
  exploration_weight: 2.0
```

## Integration with Existing Modules

AMR features seamlessly integrate with existing Khukuri modules:

- **Target Discovery**: Enhanced with pathogen-specific targets
- **Molecule Design**: Guided by resistance predictions
- **Docking**: Prioritized by druggability scores
- **ADMET**: Combined with MIC data
- **Resistance**: Expanded with genomics analysis

## Examples

See `examples/amr_discovery_example.py` for comprehensive examples:

1. Basic AMR discovery
2. Resistance pattern analysis
3. MIC data tracking
4. World model demonstration
5. Complete workflow

## Performance

- **Scalability**: Handles 1000+ compounds per iteration
- **State Management**: Tracks 10,000+ entities
- **Knowledge Graph**: Millions of relationships
- **Learning Loop**: Converges in 3-5 iterations

## Future Enhancements

- [ ] Real-time CARD database sync
- [ ] Wet-lab integration (AutoLabs-style)
- [ ] Deep learning resistance prediction
- [ ] Clinical trial data integration
- [ ] Multi-species co-infection modeling

## References

- CARD: Comprehensive Antibiotic Resistance Database
- ResFinder: Resistance gene identification
- MEGARes: Antimicrobial resistance database
- WHO Priority Pathogens List
- CLSI/EUCAST breakpoint guidelines
