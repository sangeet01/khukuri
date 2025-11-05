#!/bin/bash
# Run Khukuri test suite

echo "Running Khukuri test suite..."

# Run all tests
pytest tests/ -v

# Run with coverage
echo ""
echo "Running with coverage..."
pytest tests/ --cov=src --cov-report=html --cov-report=term

echo ""
echo "✓ Tests complete! Coverage report in htmlcov/index.html"
