# Khukuri v2.0 - Quick Start Guide

## What's New in v2.0

Khukuri now includes production-grade AMR (Antimicrobial Resistance) discovery features:

- **Bio-Knowledge Layer**: CARD/ResFinder/MEGARes resistance genes, WHO priority pathogens, validated targets
- **Genomics Analysis**: Mutation profiling, evolution prediction, multi-omics integration
- **Microbiology Tools**: MIC analysis with CLSI/EUCAST breakpoints, assay tracking, strain management
- **World Model**: Kosmos-style state tracking, knowledge graph, Bayesian optimization
- **Multi-Agent System**: 5 specialized AI agents for collaborative analysis

## Installation

```bash
cd khukuri
pip install -r requirements.txt
```

## Quick Test (30 seconds)

```bash
python scripts/test_amr_features.py
```

Expected output:
```
[SUCCESS] All AMR features working correctly!
```

## Basic Usage

### 1. AMR Discovery Workflow

```python
from src.workflows import run_amr_discovery

# Discover compounds for M. tuberculosis
results = run_amr_discovery(
    pathogen='Mycobacterium tuberculosis',
    priority='critical',
    n_compounds=20,
    n_iterations=3
)

print(f"Targets identified: {len(results['targets'])}")
print(f"Top recommendation: {results['recommendations'][0]}")
```

### 2. Resistance Analysis

```python
from src.bioknowledge import ResistanceDatabase
from src.genomics import ResistanceGenomicsAnalyzer

# Query resistance genes
db = ResistanceDatabase()
genes = db.get_genes_by_organism('S. aureus')
print(f"Resistance genes: {genes}")

# Analyze mutations
analyzer = ResistanceGenomicsAnalyzer()
profile = analyzer.analyze_mutation_profile('rpoB', ['S531L', 'H526Y'])
print(f"Predicted resistance: {profile['predicted_resistance']}")
```

### 3. MIC Tracking

```python
from src.microbiology import MICAnalyzer

# Track MIC data
mic = MICAnalyzer()
mic.add_mic_result('COMP_001', 'ATCC_25923', 'rifampicin', 0.5)

# Get profile
profile = mic.get_mic_profile('ATCC_25923')
print(f"Resistance index: {mic.calculate_resistance_index('ATCC_25923')}")
```

### 4. World Model

```python
from src.world_model import WorldStateTracker, KnowledgeGraph

# Track discovery state
world = WorldStateTracker()
world.update_compound('COMP_001', {'smiles': 'CCO', 'activity': 'active'})
world.add_hypothesis('Compound binds InhA', ['docking'], 0.85)

# Build knowledge graph
kg = KnowledgeGraph()
kg.add_compound('COMP_001', {'mw': 250})
kg.add_target('InhA', {'druggability': 0.9})
kg.add_binding('COMP_001', 'InhA', affinity=-8.5)

# Find multi-target compounds
multi_target = kg.find_multi_target_compounds(min_targets=2)
```

## Run Examples

```bash
# Complete AMR discovery examples
python examples/amr_discovery_example.py
```

This runs 5 examples:
1. Basic AMR discovery
2. Resistance pattern analysis
3. MIC data tracking
4. World model demonstration
5. Complete workflow

## Run Tests

```bash
# Test all AMR modules
python -m pytest tests/test_bioknowledge/ -v
python -m pytest tests/test_genomics/ -v
python -m pytest tests/test_microbiology/ -v
python -m pytest tests/test_world_model/ -v
python -m pytest tests/test_integration/ -v
```

## Configuration

Edit `config/amr_config.yaml`:

```yaml
amr_discovery:
  default_pathogen: "Mycobacterium tuberculosis"
  n_compounds_per_iteration: 20
  n_iterations: 3

agents:
  use_ai: false  # Set true with OpenAI API key
  model: "gpt-3.5-turbo"
```

## Key Features

### Bio-Knowledge Layer
- **ResistanceDatabase**: 7 resistance genes, query by organism/drug class
- **PathogenDatabase**: 6 WHO priority pathogens with essential genes
- **TargetProteinDB**: 6 validated targets with druggability scores

### Genomics Layer
- **ResistanceGenomicsAnalyzer**: Mutation profiling, evolution prediction
- **MutationTracker**: Strain comparison, co-occurring mutations
- **OmicsIntegrator**: Multi-omics data integration

### Microbiology Layer
- **MICAnalyzer**: CLSI/EUCAST breakpoints, cross-resistance prediction
- **AssayTracker**: Biological screening data management
- **StrainManager**: 4 ATCC reference strains

### World Model
- **WorldStateTracker**: Complete state tracking (compounds, targets, strains, assays)
- **KnowledgeGraph**: Entity relationships, multi-target identification
- **LearningLoop**: Bayesian optimization, active learning

### Multi-Agent System
- **MicrobiologyAgent**: Pathogen biology, AMR mechanisms
- **GenomicsAgent**: Resistance genomics, mutations
- **CheminformaticsAgent**: Molecular design, QSAR
- **ResistanceCriticAgent**: Resistance risk assessment
- **LiteratureAgent**: Evidence synthesis

## Documentation

- **Full Guide**: `docs/amr_features.md`
- **Implementation Summary**: `AMR_IMPLEMENTATION_SUMMARY.md`
- **Checklist**: `IMPLEMENTATION_CHECKLIST.md`

## Architecture

```
AMR Discovery System
├── Bio-Knowledge (resistance genes, pathogens, targets)
├── Genomics (mutations, evolution, omics)
├── Microbiology (MIC, assays, strains)
├── World Model (state tracking, knowledge graph, learning)
└── Multi-Agent System (5 specialized agents)
```

## Next Steps

1. **Explore**: Run `python examples/amr_discovery_example.py`
2. **Customize**: Edit `config/amr_config.yaml`
3. **Extend**: Add your own pathogens/targets to databases
4. **Integrate**: Use with existing Khukuri modules (docking, ADMET)
5. **Deploy**: Production-ready with tests and documentation

## Support

- Examples: `examples/amr_discovery_example.py`
- Tests: `tests/test_*/`
- Docs: `docs/amr_features.md`
- Config: `config/amr_config.yaml`

## Version

Khukuri Virtual Lab v2.0 - AMR-Focused Drug Discovery Platform

**Status**: Production-ready ✓
