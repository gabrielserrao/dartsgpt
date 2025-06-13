# Dead Oil (Black Oil) Example

This example demonstrates a simple two-phase flow model without gas or compositional effects.

## Model Description

- **Components:** Water, Oil
- **Phases:** Water, Oil
- **Grid:** 1D structured (100×1×1)
- **Wells:** 1 water injector, 1 producer
- **Physics:** Immiscible displacement

## Natural Language Prompt Examples

### Basic Waterflood
"Create a simple waterflood model in a 1D reservoir with 100 cells."

### Oil Production
"I need to simulate oil production using water injection. The reservoir should be 1 km long with 100 mD permeability and 20% porosity."

### Detailed Specification
"Build a dead oil reservoir model with:
- Two-phase flow (water and oil)
- 100 cells in x-direction, 10m each
- Permeability: 300 mD in all directions
- Porosity: 0.2
- Water injection at the left boundary
- Oil production at the right boundary
- Initial oil saturation: 1.0
- Runtime: 300 days"

### Variations

1. **2D Waterflood:**
   "Create a 2D areal waterflood with 50×50 cells and a five-spot pattern"

2. **Heterogeneous Reservoir:**
   "Add layers with different permeabilities: 100 mD, 500 mD, and 50 mD"

3. **Multiple Wells:**
   "Create a line drive pattern with 3 injectors and 3 producers"

4. **Gravity Effects:**
   "Make it a vertical cross-section to study gravity segregation"

## Key Files

- `model.py`: Defines the reservoir model and fluid properties
- `main.py`: Runs the simulation and plots results

## Customization Points

1. **Relative Permeability:** Modify the PhaseRelPerm functions
2. **Fluid Properties:** Change density and viscosity values
3. **Grid:** Adjust dimensions and cell sizes
4. **Initial Conditions:** Change initial saturations
5. **Well Controls:** Switch between rate and pressure control