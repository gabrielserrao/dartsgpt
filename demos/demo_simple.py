#!/usr/bin/env python
"""Simple demonstration of DARTSGPT capabilities without full LLM pipeline."""

import sys
from pathlib import Path
import json
import re

# Add dartsgpt to path
sys.path.insert(0, str(Path(__file__).parent))

from dartsgpt.knowledge.templates.template_database import TEMPLATES
from dartsgpt.knowledge.physics_features import PHYSICS_KEYWORDS, UNIT_CONVERSIONS
from dartsgpt.knowledge.example_prompts import EXAMPLE_PROMPTS


def simple_intent_classifier(prompt: str) -> dict:
    """Simple keyword-based intent classifier."""
    prompt_lower = prompt.lower()
    
    # Detect physics type
    physics_type = "dead_oil"  # default
    for ptype, keywords in PHYSICS_KEYWORDS.items():
        if any(kw in prompt_lower for kw in keywords):
            physics_type = ptype
            break
    
    # Detect features
    features = []
    if "injection" in prompt_lower:
        features.append("injection")
    if "production" in prompt_lower:
        features.append("production")
    if "thermal" in prompt_lower or "temperature" in prompt_lower:
        features.append("thermal")
    if "five-spot" in prompt_lower or "5-spot" in prompt_lower:
        features.append("five-spot")
    
    # Detect complexity
    if "simple" in prompt_lower or "basic" in prompt_lower:
        complexity = "simple"
    elif "complex" in prompt_lower or "advanced" in prompt_lower:
        complexity = "complex"
    else:
        complexity = "moderate"
    
    return {
        "physics_type": physics_type,
        "features": features,
        "complexity": complexity
    }


def simple_parameter_extractor(prompt: str) -> dict:
    """Simple regex-based parameter extractor."""
    params = {
        "grid_parameters": {},
        "rock_properties": {},
        "well_parameters": {}
    }
    
    # Grid dimensions
    grid_pattern = r'(\d+)\s*[xX×]\s*(\d+)\s*[xX×]\s*(\d+)'
    grid_match = re.search(grid_pattern, prompt)
    if grid_match:
        params["grid_parameters"]["nx"] = int(grid_match.group(1))
        params["grid_parameters"]["ny"] = int(grid_match.group(2))
        params["grid_parameters"]["nz"] = int(grid_match.group(3))
    
    # Porosity
    poro_pattern = r'(\d+(?:\.\d+)?)\s*%\s*porosity'
    poro_match = re.search(poro_pattern, prompt, re.I)
    if poro_match:
        params["rock_properties"]["porosity"] = float(poro_match.group(1)) / 100
    
    # Permeability
    perm_pattern = r'(\d+(?:\.\d+)?)\s*(?:md|mD)'
    perm_match = re.search(perm_pattern, prompt, re.I)
    if perm_match:
        perm_md = float(perm_match.group(1))
        params["rock_properties"]["permeability"] = perm_md  # Keep in mD
    
    # Temperature
    temp_pattern = r'(\d+(?:\.\d+)?)\s*°?[CF]'
    temp_match = re.search(temp_pattern, prompt)
    if temp_match:
        temp_value = float(temp_match.group(1))
        if 'F' in prompt[temp_match.start():temp_match.end()]:
            temp_value = UNIT_CONVERSIONS["f_to_c"](temp_value)
        params["temperature"] = temp_value
    
    return params


def simple_template_selector(physics_type: str, features: list) -> str:
    """Simple template selector based on physics type."""
    # Find templates matching physics type
    candidates = [
        (name, template) for name, template in TEMPLATES.items()
        if template.physics_type == physics_type
    ]
    
    if not candidates:
        # Fallback to dead oil
        candidates = [
            (name, template) for name, template in TEMPLATES.items()
            if template.physics_type == "dead_oil"
        ]
    
    # Prefer golden templates
    golden = [(n, t) for n, t in candidates if t.is_golden]
    if golden:
        return golden[0][0]
    
    # Return first match
    return candidates[0][0] if candidates else "2ph_do"


def generate_simple_code(template_name: str, parameters: dict) -> dict:
    """Generate simple DARTS code based on template and parameters."""
    template = TEMPLATES[template_name]
    
    # Get grid parameters
    grid = parameters.get("grid_parameters", {})
    nx = grid.get("nx", 50)
    ny = grid.get("ny", 50)
    nz = grid.get("nz", 10)
    
    # Get rock properties
    rock = parameters.get("rock_properties", {})
    poro = rock.get("porosity", 0.2)
    perm = rock.get("permeability", 100)  # permeability in mD
    
    model_code = f'''"""Generated DARTS model using template: {template_name}"""

import numpy as np
from darts.models.cicd_model import CICDModel
from darts.physics.{template.physics_type} import {template.physics_type.title().replace("_", "")}


class Model(CICDModel):
    """Generated {template.description}"""
    
    def __init__(self):
        super().__init__()
        self.physics_type = "{template.physics_type}"
        
    def set_physics(self):
        """Set up physics for the model."""
        # TODO: Configure physics based on {template.physics_type}
        pass
        
    def set_reservoir(self):
        """Set up reservoir properties."""
        # Grid dimensions
        self.nx = {nx}
        self.ny = {ny}
        self.nz = {nz}
        
        # Cell dimensions (m)
        self.dx = 10.0
        self.dy = 10.0
        self.dz = 2.0
        
        # Rock properties
        self.porosity = np.ones((self.nx, self.ny, self.nz)) * {poro}
        self.permeability = np.ones((self.nx, self.ny, self.nz)) * {perm}  # mD
        
    def set_wells(self):
        """Set up well locations and perforations."""
        # TODO: Configure wells based on requirements
        pass
        
    def set_initial_conditions(self):
        """Set initial conditions."""
        # TODO: Set initial pressure, saturation, etc.
        pass
        
    def set_well_controls(self):
        """Set well control schedules."""
        # TODO: Configure well controls
        pass
'''

    main_code = f'''"""Main script to run the DARTS model."""

import numpy as np
from model import Model


def main():
    # Create model instance
    model = Model()
    
    # Initialize model
    model.init()
    
    # Run simulation
    model.run()
    
    # Save results
    model.save_results()
    

if __name__ == "__main__":
    main()
'''

    return {
        "model_code": model_code,
        "main_code": main_code,
        "template_used": template_name,
        "parameters": parameters
    }


def demo_generation(prompt: str):
    """Demonstrate the generation process."""
    print(f"\n{'='*60}")
    print(f"PROMPT: {prompt}")
    print(f"{'='*60}")
    
    # Step 1: Intent Classification
    print("\n1. INTENT CLASSIFICATION")
    intent = simple_intent_classifier(prompt)
    print(f"   Physics type: {intent['physics_type']}")
    print(f"   Features: {intent['features']}")
    print(f"   Complexity: {intent['complexity']}")
    
    # Step 2: Parameter Extraction
    print("\n2. PARAMETER EXTRACTION")
    parameters = simple_parameter_extractor(prompt)
    print(f"   Grid: {parameters['grid_parameters']}")
    print(f"   Rock: {parameters['rock_properties']}")
    
    # Step 3: Template Selection
    print("\n3. TEMPLATE SELECTION")
    template_name = simple_template_selector(intent['physics_type'], intent['features'])
    template = TEMPLATES[template_name]
    print(f"   Selected: {template_name}")
    print(f"   Description: {template.description}")
    
    # Step 4: Code Generation
    print("\n4. CODE GENERATION")
    result = generate_simple_code(template_name, parameters)
    print(f"   Model code: {len(result['model_code'])} characters")
    print(f"   Main code: {len(result['main_code'])} characters")
    
    # Save output
    output_dir = Path("demo_output")
    output_dir.mkdir(exist_ok=True)
    
    (output_dir / "model.py").write_text(result["model_code"])
    (output_dir / "main.py").write_text(result["main_code"])
    (output_dir / "metadata.json").write_text(json.dumps({
        "prompt": prompt,
        "intent": intent,
        "parameters": parameters,
        "template": template_name
    }, indent=2))
    
    print(f"\n✅ Files saved to: {output_dir}/")
    print(f"   - model.py")
    print(f"   - main.py")
    print(f"   - metadata.json")


def main():
    """Run demonstrations."""
    print("DARTSGPT Simple Demonstration")
    print("=" * 60)
    print("\nThis demo shows DARTSGPT capabilities using simple rule-based processing")
    print("(without requiring LLM API calls)")
    
    # Demo 1: CO2 Storage
    demo_generation("Create a CO2 injection model for carbon storage in a 100x100x20 reservoir with 25% porosity")
    
    # Demo 2: Waterflood
    demo_generation("Build a waterflood model with 50x50x5 grid and 200 mD permeability")
    
    # Demo 3: Geothermal
    demo_generation("Model a geothermal reservoir at 150°C with injection and production wells")
    
    print(f"\n\n{'='*60}")
    print("Demo complete! Check the 'demo_output' directory for generated files.")
    print("\nNote: This is a simplified demo. The full system uses LLMs for:")
    print("  - More accurate intent classification")
    print("  - Better parameter extraction")
    print("  - Intelligent template selection")
    print("  - Complete code generation with proper physics setup")


if __name__ == "__main__":
    main()