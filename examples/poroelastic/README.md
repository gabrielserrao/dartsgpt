# Poroelastic Modeling Example

This example demonstrates coupled flow and geomechanics for reservoir deformation analysis.

## Model Description

- **Physics:** Single-phase flow with mechanical deformation
- **Coupling:** Two-way coupled (flow affects stress, stress affects flow)
- **Analytical Solutions:** Mandel, Terzaghi, Bai problems
- **Grid:** Unstructured mesh support
- **Mechanics:** Linear elasticity

## Natural Language Prompt Examples

### Basic Poroelastic Model
"Create a poroelastic model to study reservoir compaction due to depletion."

### Consolidation Problem
"Model the Terzaghi consolidation problem with a 1D column under constant load."

### Detailed Specification
"Build a poroelastic model with:
- Single-phase water flow
- Mechanical deformation with Biot coupling
- Young's modulus: 5 GPa
- Poisson's ratio: 0.25
- Biot coefficient: 0.8
- Initial pressure: 10 MPa
- Applied stress: 1 MPa compression
- Drainage boundary at top
- Runtime: until steady state"

### Variations

1. **Injection-Induced Seismicity:**
   "Model fault reactivation due to fluid injection with poroelastic stresses"

2. **Thermal Effects:**
   "Add thermoporoelastic coupling for cold water injection"

3. **Reservoir-Caprock System:**
   "Create a multi-layer model with different mechanical properties"

4. **Complex Loading:**
   "Apply cyclic loading to study fatigue effects"

## Key Files

- `model.py`: Defines the coupled problem setup
- `main.py`: Runs simulation and compares with analytical solutions
- `reservoir.py`: Custom reservoir class for mechanical problems

## Customization Points

1. **Mechanical Properties:** Vary Young's modulus, Poisson's ratio
2. **Coupling:** Adjust Biot coefficient and compressibility
3. **Boundary Conditions:** Apply different mechanical constraints
4. **Mesh:** Use different element types (hex, tet, wedge)
5. **Constitutive Models:** Extend to plasticity or damage