#!/usr/bin/env python
"""Script to run the DARTSGPT CLI."""

import sys
from pathlib import Path

# Add parent to path if needed
sys.path.insert(0, str(Path(__file__).parent.parent))

if __name__ == "__main__":
    from dartsgpt.cli import cli
    cli()