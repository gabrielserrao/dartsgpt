# DARTS Example Catalog

This catalog contains representative examples of different physics types available in DARTS, along with natural language prompts that could be used to generate similar models.

## 1. Compositional Flow

**Location:** `compositional/`

**Description:** Two-phase compositional flow model simulating CO2, C1 (methane), and water components in gas and oil phases. Uses constant K-value flash calculations.

**Key Features:**
- 3 components: CO2, C1, H2O
- 2 phases: gas, oil
- Constant K-value flash
- 1D reservoir (1000 cells)
- Injection and production wells

**Example Natural Language Prompts:**
- "Create a compositional model with CO2 and methane in a 1D reservoir"
- "Simulate CO2 injection into an oil reservoir with methane"
- "Model a two-phase compositional flow with CO2, methane and water"
- "Build a simple CO2 storage simulation with compositional physics"

## 2. Dead Oil (Black Oil)

**Location:** `dead_oil/`

**Description:** Two-phase dead oil model simulating water and oil flow without gas phase or compositional effects.

**Key Features:**
- 2 components: water, oil
- 2 phases: water, oil
- Simple density and viscosity models
- 1D reservoir (100 cells)
- Water injection for oil displacement

**Example Natural Language Prompts:**
- "Create a simple water flooding model"
- "Simulate oil production with water injection"
- "Model two-phase flow without gas"
- "Build a dead oil reservoir simulation"
- "Create a waterflood in a 1D reservoir"

## 3. Geothermal Systems

**Location:** `geothermal/`

**Description:** Two-phase geothermal model with water component in liquid and steam phases, including thermal effects.

**Key Features:**
- 1 component: water
- 2 phases: water, steam
- Thermal properties (enthalpy, heat conduction)
- Temperature-dependent properties
- Heat capacity and thermal conductivity

**Example Natural Language Prompts:**
- "Create a geothermal reservoir model"
- "Simulate hot water injection in a geothermal system"
- "Model steam generation in a thermal reservoir"
- "Build a two-phase geothermal simulation with heat transfer"
- "Create a hot water injection model with thermal effects"

## 4. Carbon Capture and Storage (CCS)

**Location:** `ccs/`

**Description:** CO2 storage model with negative flash calculations for CO2-brine systems, including dissolution and density effects.

**Key Features:**
- CO2-brine system
- Negative flash with cubic EOS
- Density models (Garcia 2001)
- Multi-layer reservoir (100x1x40)
- CO2 dissolution in brine

**Example Natural Language Prompts:**
- "Create a CO2 storage model in a saline aquifer"
- "Simulate CO2 injection for carbon sequestration"
- "Model CO2 dissolution in brine"
- "Build a CCS simulation with multiple layers"
- "Create a carbon storage model with density effects"

## 5. Poroelastic Modeling

**Location:** `poroelastic/`

**Description:** Coupled flow and geomechanics model for single-phase flow with mechanical deformation.

**Key Features:**
- Single-phase flow
- Mechanical deformation
- Poroelastic coupling
- Analytical solutions (Mandel, Terzaghi, Bai)
- Unstructured mesh support

**Example Natural Language Prompts:**
- "Create a poroelastic model with mechanical deformation"
- "Simulate reservoir compaction due to production"
- "Model coupled flow and geomechanics"
- "Build a single-phase flow model with stress effects"
- "Create a consolidation model with pore pressure coupling"

## 6. Chemical Reactions

**Location:** `chemical/`

**Description:** Reactive transport model with mineral dissolution using PHREEQC coupling.

**Key Features:**
- Chemical reactions via PHREEQC
- Mineral dissolution/precipitation
- Porosity changes
- Reactive transport
- Kinetic reactions

**Example Natural Language Prompts:**
- "Create a model with calcite dissolution"
- "Simulate acid injection with mineral reactions"
- "Model reactive transport with porosity changes"
- "Build a chemical flooding simulation"
- "Create a dissolution model with geochemistry"

## Usage Guidelines

1. **Starting Simple:** Begin with the simplest model that captures your physics (e.g., dead oil for basic flow)
2. **Adding Complexity:** Gradually add features like thermal effects, compositional behavior, or mechanics
3. **Parameter Values:** The examples use typical values that can be modified for specific cases
4. **Grid Resolution:** Examples use relatively coarse grids that can be refined for production runs
5. **Well Controls:** Examples show basic well configurations that can be extended

## Common Patterns

All models follow a similar structure:
1. Import necessary modules
2. Define Model class inheriting from appropriate base
3. Set reservoir properties (grid, porosity, permeability)
4. Define physics (components, phases, properties)
5. Set up wells and controls
6. Configure simulation parameters
7. Run simulation in main.py

## Tips for Natural Language Generation

When creating prompts for DARTSGPT:
- Be specific about physics type (compositional, thermal, mechanical)
- Mention key components or phases
- Specify grid dimensions if important
- Include well configuration details
- Mention any special features (reactions, thermal effects, etc.)