#!/usr/bin/env python3
"""
Browse and display DARTS example models and their natural language prompts.
"""

import os
import json
from pathlib import Path


def load_example_prompts():
    """Load the example prompts from JSON file."""
    prompts_file = Path(__file__).parent / "example_prompts.json"
    with open(prompts_file, 'r') as f:
        return json.load(f)


def list_physics_types():
    """List all available physics types."""
    examples_dir = Path(__file__).parent
    physics_types = []
    
    for item in examples_dir.iterdir():
        if item.is_dir() and not item.name.startswith('.') and not item.name.startswith('_'):
            if (item / "model.py").exists():
                physics_types.append(item.name)
    
    return sorted(physics_types)


def display_example(physics_type, prompts_data):
    """Display information about a specific physics type."""
    print(f"\n{'='*60}")
    print(f"Physics Type: {physics_type.upper()}")
    print(f"{'='*60}")
    
    # Check if README exists
    readme_path = Path(__file__).parent / physics_type / "README.md"
    if readme_path.exists():
        with open(readme_path, 'r') as f:
            lines = f.readlines()
            # Print first few lines of description
            for line in lines[2:10]:
                if line.strip() and not line.startswith('#'):
                    print(line.strip())
                if "## Model Description" in line:
                    break
    
    # Display prompts
    if physics_type in prompts_data:
        print(f"\n{'Simple Prompts:'}")
        print("-" * 40)
        for prompt in prompts_data[physics_type]["simple"]:
            print(f"• {prompt}")
        
        print(f"\n{'Detailed Prompts:'}")
        print("-" * 40)
        for i, prompt in enumerate(prompts_data[physics_type]["detailed"], 1):
            print(f"\n{i}. {prompt}")
        
        print(f"\n{'Keywords:'}")
        print("-" * 40)
        print(", ".join(prompts_data[physics_type]["keywords"]))


def main():
    """Main function to browse examples."""
    print("DARTS Example Browser")
    print("====================")
    
    # Load prompts data
    prompts_data = load_example_prompts()
    
    # List available physics types
    physics_types = list_physics_types()
    
    while True:
        print("\nAvailable physics types:")
        for i, ptype in enumerate(physics_types, 1):
            print(f"{i}. {ptype}")
        print("0. Exit")
        
        try:
            choice = input("\nSelect a physics type (0-{}): ".format(len(physics_types)))
            choice = int(choice)
            
            if choice == 0:
                print("Exiting...")
                break
            elif 1 <= choice <= len(physics_types):
                display_example(physics_types[choice-1], prompts_data)
                input("\nPress Enter to continue...")
            else:
                print("Invalid choice. Please try again.")
        except (ValueError, KeyboardInterrupt):
            print("\nInvalid input or interrupted. Exiting...")
            break


if __name__ == "__main__":
    main()