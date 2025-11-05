# Architecture

## Overview

Khukuri is built with a modular architecture separating concerns into 9 independent modules.

## Design Principles

1. **Modularity**: Each module is self-contained with clear interfaces
2. **Zero Mocks**: All code is production-ready, no placeholders
3. **Minimal Dependencies**: Modules depend only on what they need
4. **Configuration-Driven**: All parameters externalized to YAML/JSON
5. **Comprehensive Logging**: All operations logged for debugging

## Module Architecture

```
┌─────────────────────────────────────────────────────────┐
│                    Workflows Layer                       │
│  (autonomous_discovery, target_to_lead, optimization)   │
└─────────────────────────────────────────────────────────┘
                           │
        ┌──────────────────┼──────────────────┐
        ▼                  ▼                  ▼
┌──────────────┐  ┌──────────────┐  ┌──────────────┐
│   Agents     │  │   Docking    │  │  Resistance  │
│   Layer      │  │              │  │              │
└──────────────┘  └──────────────┘  └──────────────┘
        │                  │                  │
        └──────────────────┼──────────────────┘
                           ▼
        ┌──────────────────────────────────────┐
        │        Application Layer             │
        │  (target_discovery, molecule_design, │
        │   synthesis, admet)                  │
        └──────────────────────────────────────┘
                           │
                           ▼
        ┌──────────────────────────────────────┐
        │          Core Layer                  │
        │  (logger, validator, scorer)         │
        └──────────────────────────────────────┘
```

## Data Flow

1. **Target Discovery**: genes → PPI network → ranked targets
2. **Molecule Design**: scaffold → generation → optimization → candidates
3. **Docking**: receptor + ligand → binding poses → affinity scores
4. **ADMET**: molecule → properties → drug-likeness score
5. **Resistance**: target + mutations → risk score

## Key Components

### Core Module
- Provides shared utilities (logging, validation, scoring)
- No external module dependencies
- Used by all other modules

### Target Discovery
- Fetches PPI data from STRING/KEGG APIs
- Builds NetworkX graphs
- Ranks targets by centrality metrics

### Molecule Design
- Generates molecules using RDKit
- Optimizes properties iteratively
- Fragment-based and scaffold-based approaches

### Docking
- Wraps AutoDock Vina
- Prepares receptors from PDB
- Analyzes binding poses

### ADMET
- Calculates physicochemical properties
- Predicts toxicity using structural alerts
- Estimates pharmacokinetics

### Resistance
- Loads mutation databases
- Predicts resistance likelihood
- Designs multi-target strategies

### Synthesis
- Performs retrosynthetic analysis
- Plans synthesis routes
- Calculates synthetic accessibility

### Agents
- Optional AI-powered analysis
- Multi-agent orchestration
- Graceful degradation without OpenAI

### Workflows
- Orchestrates end-to-end pipelines
- Integrates all modules
- Provides high-level APIs

## Configuration Management

All configuration in `config/default_config.yaml`:

```yaml
scoring:
  weights: {...}
  thresholds: {...}

docking:
  exhaustiveness: 8
  box_size: [40, 40, 40]

target_discovery:
  max_targets: 20
```

## Error Handling

- Try-except blocks at module boundaries
- Graceful degradation when APIs fail
- Comprehensive logging for debugging
- Validation at entry points

## Testing Strategy

- Unit tests for each module
- Integration tests for workflows
- Fixtures for common test data
- Markers for slow/network tests
