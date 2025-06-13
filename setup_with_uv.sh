#!/bin/bash
# Setup script for DARTSGPT with uv

echo "DARTSGPT Setup with uv"
echo "====================="

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "❌ uv is not installed. Please install it first:"
    echo "   curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi

echo "✓ uv is installed"

# Change to dartsgpt directory
cd dartsgpt

# Create virtual environment with uv
echo "Creating virtual environment with uv..."
uv venv

# Sync dependencies
echo "Installing dependencies..."
uv sync

# Make CLI executable
echo "Making CLI executable..."
chmod +x cli.py
chmod +x bin/dartsgpt

echo ""
echo "✅ Setup complete!"
echo ""
echo "To activate the environment, run:"
echo "  source dartsgpt/.venv/bin/activate"
echo ""
echo "To test the CLI, run:"
echo "  cd dartsgpt"
echo "  uv run python -m dartsgpt --help"
echo "  # or"
echo "  uv run dartsgpt --help"