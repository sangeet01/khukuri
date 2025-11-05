# Data Directory

This directory contains automatically updated datasets from public databases.

## Auto-Update System

Data is automatically fetched daily from:
- **CARD** - Resistance mutations
- **ChEMBL** - Drug targets
- **ClinicalTrials.gov** - Clinical trials
- **PDB** - Protein structures
- **PubMed** - Recent research papers

## Setup Auto-Update

### Unix/Linux/macOS
```bash
bash scripts/setup_auto_update.sh
```

### Windows
```bash
scripts\setup_auto_update.bat
```
(Run as Administrator)

## Manual Update

```bash
python scripts/update_data.py
```

## Data Files

- `resistance_mutations.json` - Latest resistance data
- `drug_targets.json` - Validated drug targets
- `clinical_trials.json` - Recent clinical trials
- `pdb_structures.json` - New PDB entries
- `update_metadata.json` - Update timestamps

## Configuration

Edit `config/data_sources.yaml` to:
- Enable/disable sources
- Change update frequency
- Add API keys
- Configure filters

## Data Freshness

Check last update:
```python
from src.core.data_updater import DataUpdater
updater = DataUpdater()
# Check metadata file
```

## Backup

Old data is automatically backed up for 30 days in `data/backups/`
