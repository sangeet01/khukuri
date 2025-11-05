# Smart Discovery System

## Overview

Khukuri's **Smart Discovery System** combines:

1. **AI-Led Meetings**: Agents collaborate on strategy and decisions
2. **Pre-Built Tools**: Robust modules for fast execution

## Architecture

```
USER INPUT ("Staphylococcus aureus")
    ↓
LEAD SCIENTIST (creates team, leads meetings)
    ↓
TEAM MEETINGS (agents debate strategy)
    ↓
WORKFLOW DESIGN (Lead scientist designs optimal workflow)
    ↓
EXECUTION (pre-built discovery tools)
    ↓
ANALYSIS MEETING (agents evaluate results)
    ↓
RECOMMENDATIONS
```

## Key Components

### 1. Lead Scientist (`PIAgent`)

Lead scientist that orchestrates the entire discovery process.

**Capabilities:**
- Creates scientist teams based on project needs
- Convenes and leads team meetings
- Designs workflows using available tools
- Synthesizes agent discussions into decisions

**Example:**
```python
from src.agents import PIAgent

lead = PIAgent(openai_client=client)

# Create team
team = lead.create_team("Design antibiotics for MRSA")

# Design workflow
workflow = lead.design_workflow(
    available_tools=['NetworkAnalyzer', 'VinaWrapper', 'ADMETPredictor'],
    context={'species': 'Staphylococcus aureus'}
)
```

### 2. Meeting System (`MeetingSystem`)

Manages collaborative meetings between agents.

**Meeting Types:**

**Team Meeting**: All agents discuss broad topics
```python
from src.agents import MeetingSystem

meetings = MeetingSystem(lead_agent, openai_client)

result = meetings.team_meeting(
    agenda="Select drug discovery strategy",
    context={'species': 'E. coli', 'resistance': 'high'},
    agenda_questions=[
        "Should we use network pharmacology?",
        "Focus on de novo or repurposing?"
    ],
    rounds=3
)
```

**Individual Meeting**: One agent tackles specific task
```python
result = meetings.individual_meeting(
    agent_title="Medicinal Chemist",
    task="Design molecules with improved LogP",
    context={'current_molecules': molecules},
    with_critic=True,
    rounds=2
)
```

**Parallel Meetings**: Multiple meetings merged for best results
```python
result = meetings.parallel_meetings(
    agenda="Identify optimal targets",
    context={'species': 'P. aeruginosa'},
    n_parallel=3
)
```

### 3. Smart Orchestrator (`HybridOrchestrator`)

Main interface combining AI meetings with pre-built tools.

**Full Discovery Pipeline:**
```python
from src.agents import HybridOrchestrator

orchestrator = HybridOrchestrator(openai_client=client)

results = orchestrator.run_discovery(
    species_name="Staphylococcus aureus",
    disease_context="MRSA with high resistance",
    constraints={'time_limit': '2 hours', 'focus': 'novel mechanisms'}
)

# Results include:
# - Team composition
# - Strategy recommendations
# - Designed workflow
# - Execution results
# - Agent analysis
# - Next steps
```

**Quick Consultation:**
```python
response = orchestrator.quick_consult(
    question="What are key challenges for Gram-negative antibiotics?",
    context={'pathogen_type': 'Gram-negative'}
)
```

**Custom Meeting:**
```python
result = orchestrator.custom_meeting(
    agenda="Evaluate combination therapy strategies",
    context={'pathogen': 'P. aeruginosa', 'resistance': 'MDR'},
    meeting_type='team'
)
```

## Discovery Phases

### Phase 1: Team Selection
Lead scientist creates team based on project needs.

**Default Team:**
- Computational Biologist (structure analysis, docking)
- Medicinal Chemist (molecule design, ADMET)
- Resistance Specialist (resistance prediction, multi-target)

**AI-Powered:** Lead scientist creates custom teams based on specific challenges.

### Phase 2: Project Specification
Team meeting to decide strategy.

**Decisions Made:**
- Network pharmacology vs literature-based targets?
- De novo design vs drug repurposing vs both?
- Key success criteria
- Priority resistance mechanisms

### Phase 3: Tools Selection
Team meeting to select computational tools.

**Available Tools:**
- `NetworkAnalyzer`: PPI network analysis
- `TargetRanker`: Target prioritization
- `MoleculeGenerator`: Fragment-based design
- `PropertyOptimizer`: Property optimization
- `VinaWrapper`: Molecular docking
- `BindingSiteDetector`: Binding site analysis
- `ADMETPredictor`: Drug-likeness prediction
- `ResistancePredictor`: Resistance analysis
- `RetroSynthesisPlanner`: Synthesis planning
- `MolecularScorer`: Composite scoring

### Phase 4: Workflow Design
Lead scientist designs optimal workflow combining selected tools.

**Example Workflow:**
```python
{
  'workflow_steps': [
    {
      'step': 1,
      'tool': 'NetworkAnalyzer',
      'purpose': 'Discover drug targets',
      'inputs': ['species_name'],
      'outputs': ['ranked_targets', 'ppi_network']
    },
    {
      'step': 2,
      'tool': 'MoleculeGenerator',
      'purpose': 'Design candidates',
      'inputs': ['target_structure'],
      'outputs': ['candidate_molecules']
    },
    # ... more steps
  ],
  'rationale': 'Network-based multi-target approach',
  'estimated_time': '1-2 hours'
}
```

### Phase 5: Execution
Execute workflow using pre-built Khukuri tools.

**Hybrid Approach:**
- Agents **design** the workflow
- Pre-built tools **execute** the workflow
- Best of both: flexibility + reliability

### Phase 6: Results Analysis
Team meeting to analyze results and recommend leads.

**Analysis Includes:**
- Which candidates show promise?
- Key risks and mitigation strategies
- Prioritized experiments
- Next steps

## Comparison: Modes

| Aspect | Fixed Pipeline | Smart Discovery |
|--------|---------------|----------------|
| **Strategy** | Pre-defined | Agents decide |
| **Tools** | Pre-built | Pre-built (agents select) |
| **Workflow** | Fixed | Agents design |
| **Flexibility** | Low | High |
| **Reliability** | High | High |
| **Speed** | Fast | Medium |

## Usage Modes

### Mode 1: Full AI (with OpenAI)
```python
import openai

client = openai.OpenAI(api_key="your-key")
orchestrator = HybridOrchestrator(openai_client=client)

# Agents use GPT-4 for meetings and decisions
results = orchestrator.run_discovery("E. coli")
```

### Mode 2: Fallback (without OpenAI)
```python
orchestrator = HybridOrchestrator(openai_client=None)

# Uses default strategies and pre-built tools
results = orchestrator.run_discovery("E. coli")
```

### Mode 3: Hybrid (selective AI)
```python
# Use AI for strategy, fallback for execution
orchestrator = HybridOrchestrator(openai_client=client)

# Custom control over which parts use AI
result = orchestrator.custom_meeting(
    agenda="Specific question",
    context=data,
    meeting_type='team'
)
```

## Benefits

### 1. Flexibility
- Agents adapt strategy to each pathogen
- Custom workflows for unique challenges
- Not locked into fixed pipeline

### 2. Reliability
- Pre-built tools are tested and robust
- Graceful degradation without AI
- Production-ready execution

### 3. Transparency
- Meeting logs show decision process
- Understand why certain tools were chosen
- Audit trail for scientific rigor

### 4. Extensibility
- Easy to add new tools to available_tools
- Agents automatically consider new tools
- Modular architecture

## Example Output

```json
{
  "species": "Staphylococcus aureus",
  "team": [
    {
      "title": "Computational Biologist",
      "expertise": "protein structure and docking"
    },
    {
      "title": "Medicinal Chemist",
      "expertise": "small molecule design"
    },
    {
      "title": "Resistance Specialist",
      "expertise": "resistance mechanisms"
    }
  ],
  "strategy": {
    "recommendations": [
      "Use network pharmacology for multi-target approach",
      "Focus on essential proteins with low resistance",
      "Combine de novo design with repurposing"
    ]
  },
  "workflow": {
    "workflow_steps": [...],
    "estimated_time": "1-2 hours"
  },
  "execution_results": {
    "steps_completed": [...],
    "outputs": {...}
  },
  "agent_analysis": {
    "recommendations": [
      "Prioritize compounds with multi-target activity",
      "Validate top 5 candidates experimentally",
      "Monitor for resistance development"
    ]
  }
}
```

## Best Practices

1. **Provide Context**: More context = better agent decisions
2. **Set Constraints**: Time, budget, focus areas guide workflow
3. **Review Meetings**: Check meeting logs for decision rationale
4. **Iterate**: Use quick_consult for follow-up questions
5. **Validate**: Always validate computational predictions experimentally

## Future Enhancements

- [ ] Agent memory across sessions
- [ ] Learning from experimental results
- [ ] Custom agent creation by users
- [ ] Integration with lab automation
- [ ] Real-time collaboration features

## References

- Khukuri: Network pharmacology + AI agents for antimicrobials
- Smart Discovery: Flexible AI-led drug discovery
