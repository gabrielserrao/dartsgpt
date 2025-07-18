"""Example prompts for different physics types."""

EXAMPLE_PROMPTS = {
    "compositional": {
        "simple": [
            "Create a compositional model with CO2 and methane",
            "Build a two-phase compositional flow simulation",
            "Model CO2 injection with compositional physics"
        ],
        "detailed": [
            "Create a compositional reservoir model with CO2, methane, and water components. Use a 1D grid with 1000 cells, constant K-values for phase behavior, and place injection and production wells at the boundaries.",
            "Build a CO2 storage model using compositional flow. The reservoir should have 100 mD permeability, 30% porosity, and track component distribution in gas and oil phases over 1000 days.",
            "Simulate enhanced gas recovery by CO2 injection. Model a 2D reservoir (100x50 cells) with three components (CO2, C1, H2O) and monitor composition changes."
        ]
    },
    
    "dead_oil": {
        "simple": [
            "Create a waterflood model",
            "Simulate oil production with water injection",
            "Build a two-phase flow model without gas"
        ],
        "detailed": [
            "Create a waterflooding simulation in a 1D reservoir. Use 100 cells with 300 mD permeability and 20% porosity. Inject water at the left and produce oil at the right.",
            "Model a five-spot waterflood pattern in a 2D reservoir (50x50 cells). Include relative permeability effects and track oil recovery over 5 years.",
            "Build a layered reservoir model with water injection. Use three layers with different permeabilities (100, 500, 50 mD) and study the sweep efficiency."
        ]
    },
    
    "geothermal": {
        "simple": [
            "Create a geothermal reservoir model",
            "Simulate hot water injection",
            "Model a two-phase thermal system"
        ],
        "detailed": [
            "Build a geothermal reservoir model with water/steam phases. Include heat conduction with rock thermal conductivity of 2 W/m/K and heat capacity of 1000 J/kg/K. Inject water at 50°C.",
            "Create an enhanced geothermal system (EGS) model. Use a fractured reservoir with cold water injection at 50°C into a 200°C reservoir and track thermal breakthrough.",
            "Model a geothermal doublet system with production and reinjection wells. Include thermal effects and optimize well spacing to maximize heat extraction over 30 years."
        ]
    },
    
    "ccs": {
        "simple": [
            "Create a CO2 storage model",
            "Simulate CO2 injection into brine",
            "Model carbon sequestration"
        ],
        "detailed": [
            "Build a CO2 storage model in a saline aquifer with 40 layers. Use negative flash for CO2-brine equilibrium and track CO2 plume migration over 100 years.",
            "Create a CCS simulation with structural trapping. Model CO2 injection into an anticline structure with caprock seal and monitor pressure buildup and CO2 saturation.",
            "Simulate CO2 injection with dissolution and residual trapping. Use a 2D radial grid, inject CO2 for 30 years followed by 70 years of monitoring, and quantify different trapping mechanisms."
        ]
    },
    
    "poroelastic": {
        "simple": [
            "Create a poroelastic model",
            "Simulate reservoir compaction",
            "Model coupled flow and mechanics"
        ],
        "detailed": [
            "Build a poroelastic model for reservoir compaction during depletion. Use Young's modulus of 10 GPa, Poisson's ratio of 0.25, and track surface subsidence.",
            "Create a model for injection-induced seismicity. Include a fault with frictional properties and study how pore pressure changes affect fault stability.",
            "Model the Mandel-Cryer effect in a poroelastic medium. Use appropriate boundary conditions and compare numerical results with analytical solutions."
        ]
    },
    
    "chemical": {
        "simple": [
            "Create a mineral dissolution model",
            "Simulate reactive transport",
            "Model chemical flooding"
        ],
        "detailed": [
            "Build a calcite dissolution model with CO2-saturated water injection. Use PHREEQC for geochemistry and track porosity changes due to mineral dissolution.",
            "Create a reactive transport model for acid stimulation. Inject HCl to dissolve carbonates and model wormhole formation with dynamic porosity-permeability coupling.",
            "Simulate mineral scaling during water injection. Model the mixing of incompatible waters leading to barite precipitation and permeability reduction."
        ]
    }
}