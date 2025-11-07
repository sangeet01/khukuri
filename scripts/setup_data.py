"""Setup script - download resistance database and initial structures"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.bioknowledge.card_downloader import CARDDownloader
from src.core import setup_logger

logger = setup_logger('khukuri')


def main():
    print("="*60)
    print("Khukuri Data Setup")
    print("="*60)
    
    # Download CARD database
    print("\n[1/1] Setting up CARD resistance database...")
    card = CARDDownloader()
    
    success = card.download_card()
    db_file = card.update_resistance_db_file()
    
    if success:
        print(f"[OK] CARD database downloaded and parsed")
    else:
        print(f"[WARN] Using fallback resistance database (22 genes)")
    
    print(f"[OK] Database saved to: {db_file}")
    
    print("\n" + "="*60)
    print("Setup complete!")
    print("\nNote: PDB structures will be downloaded automatically")
    print("when agents decide which targets to use.")
    print("\nRun: python examples/amr_discovery_example.py")
    print("="*60)


if __name__ == '__main__':
    main()
