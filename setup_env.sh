#!/bin/bash
# Setup script for DARTSGPT environment

echo "DARTSGPT Environment Setup"
echo "=========================="

# Create virtual environment
echo "Creating virtual environment..."
python3 -m venv venv

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install requirements
echo "Installing requirements..."
pip install -r requirements.txt

# Make CLI executable
echo "Making CLI executable..."
chmod +x cli.py
chmod +x bin/dartsgpt

echo ""
echo "✅ Setup complete!"
echo ""
echo "To activate the environment, run:"
echo "  source venv/bin/activate"
echo ""
echo "To test the CLI, run:"
echo "  python -m dartsgpt --help"
echo "  # or"
echo "  ./bin/dartsgpt --help"