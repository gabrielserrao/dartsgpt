# DARTS Examples Directory

This directory contains representative examples of different physics types available in DARTS, organized to help users understand how to create natural language prompts for DARTSGPT.

## Directory Structure

```
examples/
├── compositional/     # Multi-component, multi-phase flow
├── dead_oil/         # Two-phase black oil models
├── geothermal/       # Thermal reservoir simulations
├── ccs/              # CO2 storage and sequestration
├── poroelastic/      # Coupled flow and geomechanics
├── chemical/         # Reactive transport with chemistry
├── example_catalog.md    # Detailed catalog of all examples
├── example_prompts.json  # Natural language prompts database
├── browse_examples.py    # Interactive example browser
└── README.md            # This file
```

## Quick Start

1. **Browse Examples Interactively:**
   ```bash
   python browse_examples.py
   ```

2. **View Complete Catalog:**
   Open `example_catalog.md` for a comprehensive overview of all examples.

3. **Explore Specific Physics:**
   Each subdirectory contains:
   - `model.py` - The DARTS model definition
   - `main.py` - Script to run the simulation
   - `README.md` - Detailed description and prompt examples

## Physics Types Overview

### 1. Compositional Flow
- Multi-component systems (CO2, methane, water)
- Phase equilibrium calculations
- Component tracking in multiple phases

### 2. Dead Oil (Black Oil)
- Simple two-phase flow (water-oil)
- Waterflooding and oil recovery
- No compositional effects

### 3. Geothermal Systems
- Thermal effects and heat transfer
- Water-steam phase transitions
- Temperature-dependent properties

### 4. Carbon Capture and Storage (CCS)
- CO2 injection into saline aquifers
- Dissolution and trapping mechanisms
- Long-term storage modeling

### 5. Poroelastic Modeling
- Coupled flow and mechanical deformation
- Reservoir compaction and subsidence
- Stress-dependent properties

### 6. Chemical Reactions
- Mineral dissolution/precipitation
- Reactive transport
- Dynamic porosity changes

## Using These Examples

### For Learning DARTS
1. Start with simple examples (dead_oil for basic flow)
2. Study the model structure and physics setup
3. Run the examples and analyze outputs
4. Modify parameters to see effects

### For Creating DARTSGPT Prompts
1. Review the natural language prompts in each README
2. Check `example_prompts.json` for variations
3. Use keywords from successful examples
4. Be specific about physics, components, and grid

### Example Prompt Patterns

**Simple:** "Create a [physics_type] model with [key_features]"

**Detailed:** "Build a [physics_type] model with:
- Components: [list]
- Grid: [dimensions]
- Properties: [values]
- Wells: [configuration]
- Runtime: [duration]"

## Tips for Natural Language Prompts

1. **Be Specific:** Mention the physics type explicitly
2. **Include Key Parameters:** Grid size, components, phases
3. **Specify Goals:** What you want to simulate or study
4. **Use Domain Terms:** Use proper reservoir engineering terminology
5. **Provide Context:** Explain the application if relevant

## Common Patterns

All DARTS models follow a similar structure:
- Model class definition
- Reservoir setup (grid, properties)
- Physics configuration (components, phases)
- Well placement and controls
- Simulation parameters
- Output generation

## Need Help?

- Run `python browse_examples.py` for interactive exploration
- Check individual README files for detailed information
- Refer to `example_catalog.md` for comprehensive documentation
- Use `example_prompts.json` as a reference for natural language