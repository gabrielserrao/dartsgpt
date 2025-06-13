"""Physics and feature keywords for DARTSGPT intent classification."""

# Keywords associated with each physics type
PHYSICS_KEYWORDS = {
    "compositional": [
        "co2", "carbon dioxide", "gas injection", "miscible", "component",
        "compositional", "eos", "equation of state", "phase behavior",
        "vle", "vapor liquid equilibrium", "k-value", "fugacity",
        "gas storage", "carbon storage", "ccs", "ccus", "eor gas"
    ],
    "dead_oil": [
        "water flooding", "waterflood", "immiscible", "dead oil",
        "two phase", "2 phase", "oil water", "simple oil", "black oil simple",
        "secondary recovery", "water injection", "waterdrive"
    ],
    "black_oil": [
        "black oil", "solution gas", "dissolved gas", "three phase",
        "3 phase", "gas oil water", "pvt", "bubble point", "rs", "bo",
        "conventional", "primary recovery", "gas cap", "goc", "owc"
    ],
    "geothermal": [
        "geothermal", "heat", "temperature", "thermal", "enthalpy",
        "steam", "hot water", "heat extraction", "egs", "enhanced geothermal",
        "doublet", "thermal breakthrough", "heat exchanger", "binary cycle"
    ],
    "poroelastic": [
        "stress", "strain", "mechanics", "poroelastic", "geomechanics",
        "subsidence", "compaction", "fault", "fracture", "deformation",
        "coupled", "consolidation", "effective stress", "pore pressure coupling"
    ],
    "specialized": [
        "chemical", "polymer", "surfactant", "foam", "asphaltene",
        "reaction", "kinetic", "biodegradation", "microbial", "alkaline",
        "asp", "chemical eor", "ift", "interfacial tension", "adsorption"
    ]
}

# Keywords associated with specific features
FEATURE_KEYWORDS = {
    "thermal": [
        "temperature", "heat", "thermal", "steam", "hot", "cold",
        "heat transfer", "conduction", "convection"
    ],
    "wells": [
        "well", "injection", "production", "injector", "producer",
        "perforations", "completion", "horizontal", "vertical", "deviated"
    ],
    "aquifer": [
        "aquifer", "water drive", "aquifer support", "analytical aquifer",
        "numerical aquifer", "carter tracy", "fetkovich"
    ],
    "faults": [
        "fault", "sealing", "transmissibility", "barrier", "flow barrier",
        "compartment", "discontinuity"
    ],
    "gravity": [
        "gravity", "segregation", "gravity drainage", "gravity override",
        "density difference", "buoyancy"
    ],
    "capillary": [
        "capillary", "capillary pressure", "imbibition", "drainage",
        "wettability", "contact angle"
    ],
    "relative_permeability": [
        "relative permeability", "kr", "kro", "krw", "krg",
        "saturation function", "corey", "brooks corey"
    ],
    "pvt": [
        "pvt", "fluid properties", "viscosity", "density", "compressibility",
        "fvf", "formation volume factor", "gas oil ratio", "gor"
    ],
    "grid": [
        "grid", "mesh", "cells", "blocks", "cartesian", "corner point",
        "unstructured", "voronoi", "pebi"
    ],
    "heterogeneity": [
        "heterogeneous", "layers", "facies", "permeability distribution",
        "porosity distribution", "anisotropy", "kv/kh"
    ]
}

# Complexity indicators
COMPLEXITY_INDICATORS = {
    "simple": [
        "simple", "basic", "test", "example", "tutorial",
        "homogeneous", "single", "uniform"
    ],
    "moderate": [
        "typical", "standard", "common", "field", "realistic",
        "heterogeneous", "multiple", "layers"
    ],
    "complex": [
        "complex", "advanced", "detailed", "full field",
        "history match", "optimization", "uncertainty",
        "coupled", "multiphysics", "compositional"
    ]
}

# Unit conversion factors
UNIT_CONVERSIONS = {
    # Length
    "ft_to_m": 0.3048,
    "m_to_ft": 3.28084,
    
    # Pressure
    "psi_to_pa": 6894.76,
    "bar_to_pa": 100000,
    "atm_to_pa": 101325,
    
    # Temperature
    "f_to_c": lambda f: (f - 32) * 5/9,
    "c_to_f": lambda c: c * 9/5 + 32,
    "c_to_k": lambda c: c + 273.15,
    
    # Volume
    "bbl_to_m3": 0.158987,
    "ft3_to_m3": 0.0283168,
    
    # Rate
    "stb_d_to_m3_s": 1.84013e-6,
    "mscf_d_to_m3_s": 0.327741,
    
    # Permeability
    "md_to_m2": 9.869233e-16,
    "d_to_m2": 9.869233e-13
}