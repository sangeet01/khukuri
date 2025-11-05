"""Resistance mechanisms database"""

import json
import logging
from pathlib import Path

logger = logging.getLogger('khukuri')


def load_resistance_mechanisms(config_path=None):
    """Load resistance mechanisms from JSON"""
    if config_path is None:
        config_path = Path(__file__).parent.parent.parent / 'config' / 'resistance_db.json'
    
    try:
        with open(config_path) as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Failed to load resistance DB: {e}")
        return {}
