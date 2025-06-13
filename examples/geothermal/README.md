# Geothermal System Example

This example demonstrates a two-phase geothermal model with water/steam phases and thermal effects.

## Model Description

- **Component:** Water
- **Phases:** Liquid water, Steam
- **Grid:** 1D structured (500×1×1)
- **Thermal:** Heat conduction and convection
- **Wells:** 1 hot water injector, 1 producer

## Natural Language Prompt Examples

### Basic Geothermal Model
"Create a geothermal reservoir model with water and steam phases including temperature effects."

### Hot Water Injection
"Simulate hot water injection at 300°C into a geothermal reservoir. The model should track temperature changes and phase transitions between water and steam."

### Detailed Specification
"Build a geothermal reservoir model with:
- Single component (H2O) in water and steam phases
- 500 cells, each 10m × 10m × 1m
- Permeability: 300 mD
- Porosity: 0.2
- Rock heat capacity: 2200 J/kg/K
- Rock thermal conductivity: 181.44 W/m/K
- Injection temperature: 300°C
- Initial reservoir temperature: 200°C
- Runtime: 1000 days"

### Variations

1. **Enhanced Geothermal System (EGS):**
   "Create an EGS model with a fractured reservoir and cold water injection"

2. **Doublet System:**
   "Model a geothermal doublet with reinjection of produced fluids"

3. **Multi-layer System:**
   "Add caprock and basement layers with different thermal properties"

4. **Natural Convection:**
   "Create a 2D model to study natural convection cells"

## Key Files

- `model.py`: Defines the thermal reservoir model and properties
- `main.py`: Runs the simulation and analyzes results

## Customization Points

1. **Thermal Properties:** Modify heat capacity and thermal conductivity
2. **Phase Behavior:** Adjust steam table properties
3. **Grid:** Create 2D or 3D models for complex flow patterns
4. **Boundary Conditions:** Add heat loss to caprock/basement
5. **Well Design:** Implement complex well geometries