"""
DARTS Model Generated from User Prompt
=====================================

Original prompt: "Build a ccs model of a reservoir with porosity 0.25, perm 100 md and 2 wells"

Model Configuration:
- Physics Type: compositional (Basic two-phase compositional flow model)
- Template Used: 2ph_comp
- Intent Classification: {'physics_type': 'compositional', 'features': ['Flow type: injection', 'Well types: vertical'], 'complexity': 'simple'}

This model was automatically generated based on your natural language description.
The code includes explanatory comments showing how your requirements were interpreted.
"""

import numpy as np
from darts.models.reservoirmodel import ReservoirModel
from darts.physics.compositional import Compositional
from darts.tools.keyword_file_tools import load_single_keyword


class Model(ReservoirModel):
    """
    Basic two-phase compositional flow model
    
    This model implements compositional physics
    as requested in your prompt. The specific physics type was chosen because:
    Physics type 'compositional' was selected based on keywords in your prompt
    """
    
    def __init__(self):
        """Initialize the DARTS reservoir model."""
        super().__init__()
        self.physics_type = "compositional"
        # This physics type handles Flow type: injection, Well types: vertical features
        
    def set_physics(self):
        """
        Set up the physics engine for the simulation.
        
        Using Compositional physics which provides:
        - Compositional flow with EOS
        - 
        - 
        - 
        - 
        """
        self.physics = Compositional()
        
    def set_reservoir(self):
        """
        Define reservoir geometry and rock properties.
        
        Grid dimensions using default values (not specified in prompt)
        """
        # Grid dimensions
        self.nx = 50  # Number of cells in X direction
        self.ny = 50  # Number of cells in Y direction  
        self.nz = 10  # Number of cells in Z direction (layers)
        
        # Cell dimensions (m)
        self.dx = 10.0  # Cell size in X direction
        self.dy = 10.0  # Cell size in Y direction
        self.dz = 2.0   # Cell thickness
        
        # Total reservoir dimensions: 500m x 500m x 20m
        
        # Rock properties
        # Using default porosity (not specified)
        self.porosity = np.ones((self.nx, self.ny, self.nz)) * 0.2
        
        # Using default permeability (not specified)
        self.permeability = np.ones((self.nx, self.ny, self.nz)) * 9.869233e-14  # m²
        # Note: 9.869233e-14 m² = 100.0 mD
        
        # Initialize the reservoir grid with these properties
        self.init_reservoir()
        
    def set_wells(self):
        """
        Configure well locations and completions.
        
        Using default well configuration
        """
        # Injector well
        # Default injector added
        self.add_well(
            "INJ1",                      # Well name
            welltype='injector',         # Well type
            i=1, j=1,                   # Location at corner of reservoir
            k_top=1, k_bottom=self.nz   # Completed through all layers
        )
        
        # Producer well  
        # Default producer added
        self.add_well(
            "PROD1",                     # Well name
            welltype='producer',         # Well type
            i=self.nx, j=self.ny,       # Location at opposite corner
            k_top=1, k_bottom=self.nz   # Completed through all layers
        )
        
        # Note: This creates a diagonal flow pattern across the reservoir
        
    def set_initial_conditions(self):
        """
        Set initial reservoir conditions.
        
        Using default temperature
        Using default pressure
        """
        self.set_uniform_initial_conditions(
            pressure=200.0,    # Initial pressure in bar
            temperature=50.0      # Initial temperature in °C
        )
        
        # Special initial conditions for compositional
        
    def set_well_controls(self):
        """
        Define well operating conditions and control modes.
        
        Well controls are set based on the simulation objectives:
        
        
        
        """
        # Injector operating at constant rate
        # Using default injection rate
        self.set_well_rate(
            "INJ1", 
            rate=100.0  # m³/day 
        )
        
        # Producer operating at constant bottom hole pressure
        # Using default production pressure
        self.set_well_bhp(
            "PROD1", 
            bhp=180.0   # bar (bottom hole pressure)
        )
