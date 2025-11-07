"""Daily data update script - run via cron/task scheduler"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.core.data_updater import update_databases
from src.core import setup_logger

logger = setup_logger(__name__)


def main():
    """Run daily data update"""
    logger.info("=" * 60)
    logger.info("Daily Data Update")
    logger.info("=" * 60)
    
    try:
        results = update_databases(force=False)
        
        logger.info("\nUpdate Results:")
        for dataset, status in results.items():
            logger.info(f"  {dataset}: {status}")
        
        # Count updates
        updated = sum(1 for s in results.values() if s == "updated")
        logger.info(f"\nTotal datasets updated: {updated}/{len(results)}")
        
        return 0
    except Exception as e:
        logger.error(f"Update failed: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
