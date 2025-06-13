#!/usr/bin/env python
"""Test CLI commands without full agent pipeline."""

import sys
from pathlib import Path

# Add dartsgpt to path
sys.path.insert(0, str(Path(__file__).parent))

print("DARTSGPT CLI Test")
print("=" * 50)

# Test list-templates command
print("\n1. Testing list-templates command...")
try:
    from dartsgpt.knowledge.templates.template_database import TEMPLATES
    
    print(f"\nAvailable DARTS Templates ({len(TEMPLATES)} total):\n")
    
    for i, (name, template) in enumerate(list(TEMPLATES.items())[:5]):
        star = "⭐" if template.is_golden else "  "
        print(f"{star} {name:20} - {template.description}")
        print(f"    Physics: {template.physics_type:15} Complexity: {template.complexity}")
    
    print(f"\n... and {len(TEMPLATES) - 5} more templates")
except Exception as e:
    print(f"   ❌ Error: {e}")

# Test show-examples command
print("\n\n2. Testing show-examples command...")
try:
    from dartsgpt.knowledge.example_prompts import EXAMPLE_PROMPTS
    
    physics_type = "compositional"
    examples = EXAMPLE_PROMPTS.get(physics_type, {})
    
    print(f"\nExamples for {physics_type.upper()} physics:\n")
    
    print("Simple prompts:")
    for prompt in examples.get('simple', [])[:2]:
        print(f"  • {prompt}")
    
    print("\nDetailed prompts:")
    for prompt in examples.get('detailed', [])[:2]:
        print(f"  • {prompt}")
except Exception as e:
    print(f"   ❌ Error: {e}")

# Test intent parsing (simplified)
print("\n\n3. Testing basic intent parsing (without LLM)...")
try:
    from dartsgpt.knowledge.physics_features import PHYSICS_KEYWORDS
    
    test_prompts = [
        "Create a CO2 injection model",
        "Build a waterflood simulation",
        "Model a geothermal reservoir"
    ]
    
    for prompt in test_prompts:
        prompt_lower = prompt.lower()
        detected_physics = None
        
        # Simple keyword matching
        for physics_type, keywords in PHYSICS_KEYWORDS.items():
            if any(kw in prompt_lower for kw in keywords):
                detected_physics = physics_type
                break
        
        print(f"\nPrompt: {prompt}")
        print(f"Detected physics: {detected_physics or 'unknown'}")
except Exception as e:
    print(f"   ❌ Error: {e}")

# Test parameter extraction (simplified)
print("\n\n4. Testing basic parameter extraction (without LLM)...")
try:
    import re
    
    test_prompt = "Create a 50x50x10 reservoir with 25% porosity and 100 mD permeability"
    
    # Grid pattern
    grid_pattern = r'(\d+)\s*[xX×]\s*(\d+)\s*[xX×]\s*(\d+)'
    grid_match = re.search(grid_pattern, test_prompt)
    
    # Porosity pattern
    poro_pattern = r'(\d+(?:\.\d+)?)\s*%\s*porosity'
    poro_match = re.search(poro_pattern, test_prompt, re.I)
    
    # Permeability pattern
    perm_pattern = r'(\d+(?:\.\d+)?)\s*(?:md|mD)'
    perm_match = re.search(perm_pattern, test_prompt, re.I)
    
    print(f"\nPrompt: {test_prompt}")
    print("\nExtracted parameters:")
    
    if grid_match:
        print(f"  Grid: {grid_match.group(1)}x{grid_match.group(2)}x{grid_match.group(3)}")
    
    if poro_match:
        print(f"  Porosity: {poro_match.group(1)}%")
    
    if perm_match:
        print(f"  Permeability: {perm_match.group(1)} mD")
except Exception as e:
    print(f"   ❌ Error: {e}")

print("\n\n" + "=" * 50)
print("CLI test complete!")
print("\nTo use the full CLI with agents, run:")
print("  python -m dartsgpt --help")
print("\nOr if installed:")
print("  dartsgpt --help")