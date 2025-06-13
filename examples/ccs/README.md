# Carbon Capture and Storage (CCS) Example

This example demonstrates CO2 injection into a saline aquifer for carbon sequestration.

## Model Description

- **Components:** CO2, Brine (with dissolved components)
- **Phases:** Gas (CO2-rich), Aqueous (brine)
- **Flash:** Negative flash with cubic EOS
- **Grid:** 2D structured (100×1×40)
- **Physics:** CO2 dissolution, density effects

## Natural Language Prompt Examples

### Basic CCS Model
"Create a CO2 storage model in a saline aquifer with 40 layers."

### CO2 Sequestration
"I need to simulate CO2 injection into a deep saline aquifer. The model should include CO2 dissolution in brine and track the CO2 plume migration."

### Detailed Specification
"Build a CCS model with:
- CO2-brine system with dissolution
- 100×1×40 cells (radial-like grid in x)
- Logarithmic spacing in x-direction
- 40 layers, each 2m thick
- Depth: 1000m to 1080m
- Permeability: 100 mD horizontal, 10 mD vertical
- Porosity: 0.2
- CO2 injection at bottom corner
- Open boundary at far end
- Runtime: 30 years injection + 70 years monitoring"

### Variations

1. **Structural Trapping:**
   "Add an anticline structure to study structural trapping of CO2"

2. **Multiple Wells:**
   "Create a multi-well injection scenario with 3 CO2 injectors"

3. **Caprock Integrity:**
   "Include a low-permeability caprock layer and study CO2 breakthrough"

4. **Residual Trapping:**
   "Model CO2 injection followed by water chase to enhance residual trapping"

## Key Files

- `model.py`: Defines the CCS model with complex thermodynamics
- `main.py`: Runs long-term injection and monitoring

## Customization Points

1. **Thermodynamics:** Switch between different EOS models
2. **Dissolution:** Modify CO2 solubility correlations
3. **Grid:** Create truly radial or unstructured grids
4. **Heterogeneity:** Add realistic permeability distributions
5. **Monitoring:** Track different storage mechanisms