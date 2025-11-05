#!/bin/bash
# Setup AutoDock Vina

echo "Setting up AutoDock Vina..."

# Detect OS
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo "Detected Linux"
    sudo apt-get update
    sudo apt-get install -y autodock-vina openbabel
elif [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Detected macOS"
    brew install autodock-vina open-babel
else
    echo "Windows detected - manual installation required"
    echo "Download from: http://vina.scripps.edu"
fi

# Verify installation
if command -v vina &> /dev/null; then
    echo "✓ AutoDock Vina installed successfully"
    vina --version
else
    echo "⚠️  Installation failed or manual setup required"
fi
