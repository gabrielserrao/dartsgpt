#!/usr/bin/env python
"""Quick example runner for DARTSGPT."""

import asyncio
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from cli import generate


async def run_examples():
    """Run several example generations."""
    
    examples = [
        {
            "name": "CO2 Storage",
            "prompt": "Create a CO2 injection model for carbon storage in a saline aquifer",
            "output_dir": "output/co2_storage"
        },
        {
            "name": "Waterflood",
            "prompt": "Create a five-spot waterflood pattern in a 50x50x3 reservoir",
            "output_dir": "output/waterflood"
        },
        {
            "name": "Geothermal",
            "prompt": "Build a geothermal doublet system with 50°C injection into 200°C reservoir",
            "output_dir": "output/geothermal"
        }
    ]
    
    print("DARTSGPT Example Runner")
    print("=" * 50)
    
    for example in examples:
        print(f"\n\nGenerating: {example['name']}")
        print(f"Prompt: {example['prompt']}")
        print("-" * 50)
        
        try:
            await generate(
                prompt=example['prompt'],
                output_dir=example['output_dir'],
                verbose=True,
                no_validation=False
            )
            print(f"\n✅ {example['name']} generation complete!")
            print(f"   Files saved to: {example['output_dir']}/")
        except Exception as e:
            print(f"\n❌ {example['name']} generation failed: {e}")
    
    print("\n" + "=" * 50)
    print("All examples complete!")
    print("\nGenerated files:")
    for example in examples:
        output_path = Path(example['output_dir'])
        if output_path.exists():
            print(f"\n{example['name']}:")
            for file in output_path.glob("*.py"):
                print(f"  - {file}")


if __name__ == "__main__":
    asyncio.run(run_examples())