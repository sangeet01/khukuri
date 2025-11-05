#!/bin/bash
# Install Khukuri dependencies

echo "Installing Khukuri Virtual Lab dependencies..."

python -m pip install --upgrade pip
pip install -r requirements.txt
pip install pytest pytest-cov black flake8

echo "✓ Dependencies installed!"
echo "Run 'pytest' to verify installation"
