#!/usr/bin/env python
"""Simple test for DARTSGPT components without full langchain."""

import os
import sys
from pathlib import Path

# Add dartsgpt to path
sys.path.insert(0, str(Path(__file__).parent))

# Set environment variable to avoid langchain import
os.environ["OPENAI_API_KEY"] = os.getenv("OPENAI_API_KEY", "test-key")

print("DARTSGPT Simple Component Test")
print("=" * 50)

# Test 1: Config loading
print("\n1. Testing Configuration...")
try:
    from dartsgpt.config import get_settings
    settings = get_settings()
    print(f"   ✅ Settings loaded")
    print(f"   ✅ Model: {settings.openai_model}")
    print(f"   ✅ API Key: {'*' * 20}...{settings.openai_api_key[-4:] if settings.openai_api_key else 'Not set'}")
except Exception as e:
    print(f"   ❌ Config error: {e}")

# Test 2: Template database
print("\n2. Testing Template Database...")
try:
    from dartsgpt.knowledge.templates.template_database import TEMPLATES
    print(f"   ✅ Templates loaded: {len(TEMPLATES)} available")
    
    # Show first 3 templates
    for i, (name, template) in enumerate(list(TEMPLATES.items())[:3]):
        print(f"   - {name}: {template.physics_type} physics")
except Exception as e:
    print(f"   ❌ Template error: {e}")

# Test 3: Physics features
print("\n3. Testing Physics Features...")
try:
    from dartsgpt.knowledge.physics_features import PHYSICS_KEYWORDS, UNIT_CONVERSIONS
    print(f"   ✅ Physics keywords loaded: {len(PHYSICS_KEYWORDS)} types")
    print(f"   ✅ Unit conversions loaded: {len(UNIT_CONVERSIONS)} conversions")
    
    # Test a conversion
    ft_to_m = UNIT_CONVERSIONS.get("ft_to_m", 0)
    print(f"   ✅ Test conversion: 100 ft = {100 * ft_to_m:.2f} m")
except Exception as e:
    print(f"   ❌ Physics features error: {e}")

# Test 4: Example prompts
print("\n4. Testing Example Prompts...")
try:
    from dartsgpt.knowledge.example_prompts import EXAMPLE_PROMPTS
    print(f"   ✅ Example prompts loaded: {len(EXAMPLE_PROMPTS)} physics types")
    
    # Show physics types
    print(f"   Available types: {', '.join(EXAMPLE_PROMPTS.keys())}")
except Exception as e:
    print(f"   ❌ Example prompts error: {e}")

# Test 5: CLI commands (without running them)
print("\n5. Testing CLI Structure...")
try:
    from dartsgpt.cli import cli
    commands = []
    for name, cmd in cli.commands.items():
        commands.append(name)
    print(f"   ✅ CLI commands available: {', '.join(commands)}")
except Exception as e:
    print(f"   ❌ CLI error: {e}")

print("\n" + "=" * 50)
print("Simple test complete!")
print("\nTo run the full system, you may need to:")
print("1. Install compatible versions of langchain and pydantic")
print("2. Use: pip install langchain==0.0.350 pydantic==1.10.13")
print("3. Or create a fresh virtual environment")