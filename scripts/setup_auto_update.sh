#!/bin/bash
# Setup automatic daily updates via cron

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "Setting up automatic daily data updates..."

# Create cron job
CRON_CMD="0 2 * * * cd $PROJECT_DIR && python scripts/update_data.py >> logs/data_update.log 2>&1"

# Check if cron job already exists
(crontab -l 2>/dev/null | grep -v "update_data.py"; echo "$CRON_CMD") | crontab -

echo "✓ Cron job installed: Daily update at 2:00 AM"
echo "  Command: $CRON_CMD"
echo ""
echo "To view cron jobs: crontab -l"
echo "To remove: crontab -e (then delete the line)"
echo ""
echo "Creating logs directory..."
mkdir -p "$PROJECT_DIR/logs"

echo "✓ Setup complete!"
