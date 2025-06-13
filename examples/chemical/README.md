# Chemical Dissolution Example

This example demonstrates reactive transport with mineral dissolution using PHREEQC coupling.

## Model Description

- **Chemistry:** PHREEQC geochemical engine
- **Reactions:** Calcite dissolution/precipitation
- **Transport:** Advection with reactions
- **Porosity:** Dynamic porosity changes
- **Kinetics:** Rate-limited reactions

## Natural Language Prompt Examples

### Basic Dissolution Model
"Create a model for calcite dissolution with acidic water injection."

### Reactive Transport
"Simulate CO2-saturated water injection causing mineral dissolution and porosity changes in a carbonate reservoir."

### Detailed Specification
"Build a reactive transport model with:
- Calcite mineral dissolution
- CO2-saturated injection water
- Kinetic dissolution rate: 1e-6 mol/m²/s
- Initial porosity: 0.2
- Initial calcite volume fraction: 0.1
- 100 cells in 1D
- Track pH, Ca²⁺, and porosity evolution
- Runtime: 100 days"

### Variations

1. **Multiple Minerals:**
   "Add dolomite and anhydrite to create a complex carbonate system"

2. **Precipitation:**
   "Model scale formation due to incompatible water mixing"

3. **Acidizing:**
   "Simulate acid injection for well stimulation with wormhole formation"

4. **Biogeochemistry:**
   "Include microbial reactions for enhanced oil recovery"

## Key Files

- `model.py`: Defines the reactive transport model
- `main.py`: Runs simulation and analyzes geochemical changes
- `phreeqc.dat`: PHREEQC database file

## Customization Points

1. **Minerals:** Add or remove minerals from the system
2. **Kinetics:** Modify reaction rates and mechanisms
3. **Database:** Switch between different thermodynamic databases
4. **Species:** Include additional aqueous species
5. **Coupling:** Adjust porosity-permeability relationships