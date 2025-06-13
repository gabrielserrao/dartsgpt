# Compositional Flow Example

This example demonstrates a two-phase compositional flow model with CO2, methane (C1), and water.

## Model Description

- **Components:** CO2, C1 (methane), H2O
- **Phases:** Gas, Oil
- **Flash:** Constant K-value
- **Grid:** 1D structured (1000×1×1)
- **Wells:** 1 injector, 1 producer

## Natural Language Prompt Examples

### Basic Compositional Model
"Create a compositional flow model with CO2 and methane in a 1000-cell 1D reservoir. Use constant K-values for phase behavior."

### CO2 Injection
"I need to simulate CO2 injection into a reservoir containing methane. The model should have two phases (gas and oil) with water as a third component. Place an injection well at the left boundary and a production well at the right."

### Detailed Specification
"Build a compositional reservoir model with:
- Components: CO2, C1, and H2O
- Phases: gas and oil
- 1D grid with 1000 cells, each 1m × 10m × 10m
- Permeability: 100 mD
- Porosity: 0.3
- Injection well at cell (1,1,1)
- Production well at cell (1000,1,1)
- Runtime: 1000 days"

### Variations

1. **Different Components:**
   "Create a compositional model with nitrogen and CO2 instead of methane"

2. **Different Grid:**
   "Make it a 2D model with 100×50 cells"

3. **Multiple Wells:**
   "Add two more production wells in the middle of the reservoir"

4. **Different Phase Behavior:**
   "Use Peng-Robinson EOS instead of constant K-values"

## Key Files

- `model.py`: Defines the reservoir model, physics, and properties
- `main.py`: Runs the simulation and generates output

## Customization Points

1. **Components:** Change the components list and molecular weights
2. **Phase Behavior:** Replace ConstantK with more complex flash calculations
3. **Grid:** Modify nx, ny, nz for different dimensions
4. **Properties:** Adjust permeability, porosity, and other rock properties
5. **Wells:** Add more wells or change their locations and controls