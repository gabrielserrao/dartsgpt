#!/usr/bin/env python
"""Test DARTSGPT CLI functionality without LangChain dependencies."""

import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

print("DARTSGPT CLI Test (Simple)")
print("=" * 50)

# Test configuration
print("\n1. Testing configuration...")
try:
    from dartsgpt.config import get_settings
    settings = get_settings()
    print(f"   ✓ Configuration loaded")
    print(f"   ✓ OpenAI Model: {settings.openai_model}")
    api_key_set = bool(settings.openai_api_key)
    print(f"   ✓ API Key: {'Set' if api_key_set else 'Not set'}")
except Exception as e:
    print(f"   ❌ Error: {e}")

# Test templates
print("\n2. Testing template database...")
try:
    from dartsgpt.knowledge.templates.template_database import TEMPLATES
    print(f"   ✓ Templates loaded: {len(TEMPLATES)} templates")
    # Show first 3 templates
    for i, (name, template) in enumerate(list(TEMPLATES.items())[:3]):
        print(f"   - {name}: {template.description}")
except Exception as e:
    print(f"   ❌ Error: {e}")

# Test physics features
print("\n3. Testing physics features...")
try:
    from dartsgpt.knowledge.physics_features import PHYSICS_KEYWORDS, UNIT_CONVERSIONS
    print(f"   ✓ Physics types: {list(PHYSICS_KEYWORDS.keys())}")
    print(f"   ✓ Unit conversions: {len(UNIT_CONVERSIONS)} conversions")
except Exception as e:
    print(f"   ❌ Error: {e}")

# Test example prompts
print("\n4. Testing example prompts...")
try:
    from dartsgpt.knowledge.example_prompts import EXAMPLE_PROMPTS
    print(f"   ✓ Physics with examples: {list(EXAMPLE_PROMPTS.keys())}")
except Exception as e:
    print(f"   ❌ Error: {e}")

print("\n" + "=" * 50)
print("Test complete!")