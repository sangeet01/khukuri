
# Khukuri Virtual Lab

<p align="center">
  <img src="assets/logo.png" width="400" alt="Khukuri Virtual Lab Logo"/>
</p>

<p align="center">
  <strong>Closed-loop autonomous drug discovery for antimicrobial resistance.</strong><br/>
  Any bacteria. Any mutated receptor.
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-3.0.0-blue"/>
  <img src="https://img.shields.io/badge/python-3.8%2B-green"/>
  <img src="https://img.shields.io/badge/license-MIT-orange"/>
  <img src="https://img.shields.io/badge/status-research--grade-purple"/>
</p>

---

## Why Khukuri Exists

During hospital visits as a pharmacy undergraduate, I read resistance reports that stayed with me. One in particular — a child, barely a year old, with an infection that had failed every drug in the column. Every single antibiotic. Every test drug. Resistant.

I looked for a computational system that could address this. Something that modelled resistance not as a static profile but as what it actually is — a dynamic, evolving, networked process. I didn't find one.

So I built one.

Khukuri began as a translation of Stanford Zou Group's virtual lab (designed for viral nanobody discovery) into bacterial drug discovery. That translation forced every subsequent decision. Bacteria mutate faster. Resistance transfers horizontally. Targets shift. The problem is harder. After months of continuous research and development, it became something else entirely.

---

## What Khukuri Does

Khukuri is a closed-loop autonomous drug discovery platform for antimicrobial resistance. Given a bacterial target, it:

1. **Identifies and ranks drug targets** via PPI network analysis and Numen-powered literature mining
2. **Maps the complete resistance landscape** — both mutational evolution (vertical) and horizontal gene transfer (horizontal) — using the PINCER engine
3. **Evolves drug candidates** that remain effective across all futures the bacterium can reach, not just its current state
4. **Scores candidates with real binding physics** via KeyBox voxel docking (11-channel field theory, C core)
5. **Filters for drug-likeness** via ADMET prediction and synthesisability scoring
6. **Accumulates knowledge** across every run via a persistent world model — each session starts smarter than the last

The loop closes. Each run feeds the next.

---

## The Stack

Khukuri is built on three companion projects, each solving a problem the previous one created:

| Project | Role | Why it exists |
|---|---|---|
| [linearscript](https://pypi.org/project/linearscript/) | Molecular representation | SMILES is unreliable for evolutionary algorithms. SCRIPT uses Paninian grammar — every output is valid by construction. |
| [KeyBox](https://github.com/sangeet01/keybox) | Voxel physics engine (C) | MD simulation needs a supercomputer. KeyBox is 11-channel field theory on a laptop, fast enough to serve as a fitness function inside an evolutionary loop. |
| [LimitNumen](https://github.com/sangeet01/limitnumen) | Retrieval (no training) | Character n-gram hashing, arbitrary dimension, beats BM25 at 93.9% Recall@100. Powers literature mining, threat retrieval, and compound memory. |


---

## Core Architecture

```
                    KHUKURI CLOSED LOOP
                    
  ┌─────────────────────────────────────────────────────────────┐
  │                                                             │
  │  LimitNumen Literature Mining                               │
  │       ↓                                                     │
  │  Target Discovery (PPI Network + Literature)                │
  │       ↓                                                     │
  │  ┌─────────────────────────────────────────────┐           │
  │  │           PINCER ENGINE                     │           │
  │  │                                             │           │
  │  │  Red Team (Vertical)    Red Team (HGT)      │           │
  │  │  MutationSpaceMapper    HGTMapper            │           │
  │  │  Markov Mutation Matrix KnowledgeGraph       │           │
  │  │         ↓                    ↓              │           │
  │  │    Viable Threat Cluster (432+ threats)     │           │
  │  │         ↓                                   │           │
  │  │  LimitNumen ThreatIndex (fast retrieval)    │           │
  │  │         ↓                                   │           │
  │  │  Blue Team: Darwin-Gödel Evolutionary Loop  │           │
  │  │  linearscript seeds → KeyBox fitness        │           │
  │  │  Minimax: best drug across ALL futures      │           │
  │  │         ↓                                   │           │
  │  │    Apex Skeleton Key s*                     │           │
  │  └─────────────────────────────────────────────┘           │
  │       ↓                                                     │
  │  ADMET + Synthesisability Filtering                         │
  │       ↓                                                     │
  │  World Model (persistent across sessions)                   │
  │  LimitNumen CompoundMemory → seeds next run                 │
  │       ↓                                                     │
  │  Optimised Lead Candidate                                   │
  │       ↑_____________________________________________        │
  │                    (loop closes)                            │
  └─────────────────────────────────────────────────────────────┘
```

---

## Pincer: The Resistance Engine

Pincer models the pharmacological duel between antibiotics and bacterial resistance as a zero-sum minimax game.

**The core insight:** Conventional drug discovery asks "what is the best drug for the current resistance profile?" Pincer asks "what is the best drug across all futures the bacterium can reach?" That is a fundamentally different — and harder — question.

### Red Team: Dual Resistance Mapping

**Vertical threats (mutation):**
- `MutationSpaceMapper` enumerates all viable single and double-point mutations in the binding pocket
- Markov Mutation Matrix with 2-bit DNA encoding and transition/transversion bias
- Fitness-weighted by structural and energetic cost of each mutation
- Dead Zone Index tracks O(1) explored chemical space

**Horizontal threats (HGT):**
- `HGTMapper` queries the `KnowledgeGraph` for resistance genes that could arrive via mobile genetic elements
- Transfer probability weighted by MGE type × phylogenetic distance × antibiotic selective pressure
- Built-in MRSA dataset: mecA, mecC, blaZ, vanA, vanB, aac6-aph2, tetM, norA, qnrS
- Both threat types merge into one unified threat cluster

### Blue Team: Darwin-Gödel Evolution

- Seed population generated by `MoleculeGenerator` using linearscript (100% valid molecules)
- `ThreatAwareFitnessFunction` scores each candidate against worst-case threats
- `LimitNumen ThreatIndex` retrieves top-k relevant threats per drug (6.9× speedup over full loop)
- KeyBox `NibbleEngine` provides real 11-channel voxel physics scoring when a PDB is loaded
- Evolutionary loop terminates at the Apex Skeleton Key — the drug that wins the minimax game

---

## KeyBox Integration

KeyBox is a separate voxel physics engine that replaces AutoDock Vina for in-loop fitness evaluation. It has 11 channels: steric demand, electrostatics, H-bond donor/acceptor, lipophilicity, aromaticity, metal coordination, cation/anion, hydrophobic core, solvent exposure.

**Active site detection (no prior knowledge required):**

```python
from src.integrations import KhukuriKeyBox

# Inverse vector field analysis — finds pockets automatically
kb = KhukuriKeyBox(pdb_path="PBP2a.pdb")  # center=(0,0,0) triggers auto-detect
sites = kb.detect_binding_sites(n_sites=3)
# Returns ranked sites with channel profiles: hba, hbd, lipo, arom, elec...
```

The inverse vector principle: active sites are where steric demand is low (space exists) but all other channels are high (chemistry demands interaction). `score = (1 - steric) × mean(other channels)`.

**Plug into Pincer:**

```python
from src.integrations import plug_keybox_into_pincer
from src.resistance.threat_fitness import ThreatAwareFitnessFunction

fitness_fn = ThreatAwareFitnessFunction()
plug_keybox_into_pincer(fitness_fn, kb)
# Pincer now uses real voxel physics — one call, surrogate swapped out
```

---

## LimitNumen Integration

Three use cases, one retrieval engine:

**1. Threat retrieval (Pincer)**
```python
from src.integrations import ThreatIndex
ti = ThreatIndex(dim=256)
ti.build(viable_threats)                      # index 400+ threats
worst = ti.get_worst_threats(smiles, k=20)    # retrieve top-k per drug
# 6.9x speedup over full loop evaluation
```

**2. Literature mining (target discovery)**
```python
from src.integrations import LiteratureIndex
li = LiteratureIndex()
li.index_pubmed_query("PBP2a MRSA resistance", max_results=50)
results = li.search("mecA horizontal transfer MRSA")
targets = li.extract_targets(results)  # ['mecA', 'pbp2a', 'vanA', ...]
```

**3. Compound memory (cross-session learning)**
```python
from src.integrations import CompoundMemory
cm = CompoundMemory()
cm.remember(smiles, pincer_score=0.81, admet={...})
similar = cm.recall(new_smiles, k=5)  # past winners inform new seeds
```

---

## World Model

A persistent, unified knowledge layer that accumulates across every run:

```python
from src.world_model import WorldModelManager

wm = WorldModelManager(project="mrsa_campaign")
# Restores everything from previous sessions automatically

wm.start_run("S. aureus MRSA252")
wm.record_target("PBP2a", druggability=0.9, essentiality=0.85)
wm.record_compound("COMP_001", smiles="...", pincer_score=0.81)
wm.record_pincer_result(apex_smiles="...", score=0.81, n_threats=432, ...)
summary = wm.end_run()
# State saved. Next run starts with this knowledge.
```

Five components: `WorldStateTracker`, `KnowledgeGraph`, `HypothesisEngine`, `LearningLoop`, `KosmosEngine`. All wired, all persistent.

---

## Quick Start

### Installation

```bash
git clone https://github.com/sangeet01/khukuri.git
cd khukuri
pip install -r requirements.txt

# KeyBox (sibling repo — provides real docking physics)
git clone https://github.com/sangeet01/keybox.git ../keybox

# Quick test
python scripts/quick_test.py
```

### Minimal example

```python
from src.workflows.amr_discovery import AMRDiscoveryWorkflow

workflow = AMRDiscoveryWorkflow(
    project="mrsa_run1",
    pdb_path="PBP2a_2OLV.pdb",   # optional — auto-detects binding site
)
results = workflow.run_discovery(
    pathogen="Staphylococcus aureus MRSA",
    priority="critical",
    n_iterations=3,
)
print(results["recommendations"])
```

### Pincer standalone

```python
from src.resistance.pincer_engine import PincerEngine
from src.resistance.hgt_mapper import make_pincer_with_hgt
from src.resistance.threat_fitness import ThreatAwareFitnessFunction
from src.world_model.knowledge_graph import KnowledgeGraph
from src.molecule_design.generator import MoleculeGenerator

# Dual red team
engine = make_pincer_with_hgt(population_size=100, n_generations=50)
engine.map_threats(wild_type_seq, active_site_indices)
engine.map_hgt_threats(KnowledgeGraph(), target_strain="S. aureus")

# SCRIPT-backed seeds
gen = MoleculeGenerator()
seeds = gen.seed_for_pincer(n_seeds=20, scaffold_focus="quinolone")

# Evolve
fitness_fn = ThreatAwareFitnessFunction()
apex = engine.evolve(seeds, fitness_fn=fitness_fn)
results = engine.get_results()
print(f"Apex: {results['apex_drug']['smiles']}")
print(f"Score: {results['apex_drug']['minimax_score']:.4f}")
print(f"Worst threat: {results['apex_drug']['worst_mutant']}")
```

---

## Project Structure

```
khukuri/
├── src/
│   ├── core/               # Logging, validation, scoring
│   ├── target_discovery/   # PPI network, target ranking, literature mining
│   ├── molecule_design/    # linearscript-backed generation, optimisation
│   ├── docking/            # AutoDock Vina wrapper, binding site detection
│   ├── admet/              # Drug-likeness, toxicity, pharmacokinetics
│   ├── resistance/         # PINCER engine, HGT mapper, threat fitness
│   │   ├── pincer_engine.py        # Core PINCER algorithm
│   │   ├── hgt_mapper.py           # Horizontal gene transfer red team
│   │   ├── threat_fitness.py       # ThreatAwareFitnessFunction
│   │   └── evolution_simulator.py  # Full evolution simulation
│   ├── synthesis/          # Retrosynthesis, SA scoring, route planning
│   ├── bioknowledge/       # CARD/ResFinder/MEGARes, pathogen/target DBs
│   ├── genomics/           # Resistance mutations, multi-omics
│   ├── microbiology/       # MIC analysis, assay tracking, strain management
│   ├── world_model/        # Persistent state, knowledge graph, hypotheses
│   │   └── manager.py              # WorldModelManager — unified interface
│   ├── integrations/       # KeyBox bridge, Numen index
│   │   ├── keybox_bridge.py        # KhukuriKeyBox, auto active site detection
│   │   └── numen_index.py          # ThreatIndex, LiteratureIndex, CompoundMemory
│   ├── agents/             # Multi-agent system, PI agent, peer debate
│   └── workflows/          # AMRDiscoveryWorkflow, autonomous discovery
├── tests/
├── examples/
│   └── pincer_example.py
├── scripts/
│   └── quick_test.py
└── docs/
```

---

## Dependencies

```
# Core
rdkit>=2022.09.1
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0

# Molecular representation
linearscript>=3.0.2

# Network analysis
networkx>=2.6.0
python-louvain>=0.16

# Bioinformatics
biopython>=1.79

# Utilities
requests>=2.26.0
pyyaml>=6.0

# Optional — AI agents
openai>=1.0.0
anthropic>=0.20.0

# Optional — visualisation
matplotlib>=3.4.0
py3Dmol>=2.0.0
```

**External engines (sibling repos):**
- [KeyBox](https://github.com/sangeet01/keybox) — voxel physics docking (C + Python)

---

## Status

✅ Research-grade: 80+ files, ~5,000 lines, 0% mocks, 0% duplication  
✅ PINCER dual red team (mutation + HGT) — implemented and tested  
✅ KeyBox integration — voxel docking + inverse vector active site detection  
✅ LimitNumen integration — threat retrieval, literature mining, compound memory  
✅ World model — persistent across sessions, full knowledge accumulation  
✅ linearscript backend — 100% valid molecule generation  
⏳ Wet lab validation — in progress  
⏳ Clinical dataset benchmarking — planned  

---

## Theoretical Foundations

- **PINCER** — Darwin-Gödel machine applied to pharmacological minimax: self-modifying evolutionary search over a resistance landscape that changes as you search it
- **Dual red team** — resistance as a network process (HGT) + point process (mutation), unified in one minimax objective
- **Inverse vector field** — active site detection from field-theoretic demand complementarity, not geometric cavity detection
- **Network pharmacology** — target essentiality weighted by conservation across HGT donor species; dynamic rather than static target network
- **LimitNumen** — character n-gram hashing with CRC32 and log-saturation; training-free dense retrieval beating BM25 at 93.9% Recall@100

---

## Acknowledgements

This project began with a hospital visit — I read a resistance report of a one-year-old child whose infection had failed every antibiotic in the column. Every drug. Every test. Resistant.

That report is why Khukuri exists.

**Claude (Anthropic)** — engineering partner throughout development: wiring architecture, handling large codebases, bridging theory to implementation. If I am the architect with the vision, Claude is the engineer who closed the gaps.

**Stanford Zou Group** — their virtual lab paper that started the translation.

**The AMR research community** — whose papers I devoured like a caterpillar.

---

## Citation

If you use Khukuri in your research:

```bibtex
@software{khukuri2026,
  author  = {Sharma, Sangeet},
  title   = {Khukuri Virtual Lab: Closed-loop autonomous drug discovery for AMR},
  year    = {2026},
  url     = {https://github.com/sangeet01/khukuri},
  version = {3.0.0}
}
```

---

*Built by a pharmacy undergraduate in Nepal, for the children whose infections run out of options.*

