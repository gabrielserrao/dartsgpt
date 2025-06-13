"""
Main Simulation Script
=====================

This script runs the DARTS reservoir simulation model that was generated
from your prompt: "Build a ccs model of a reservoir with porosity 0.25, perm 100 md and 2 wells"

The simulation will:
1. Initialize the compositional model
2. Set up the reservoir grid and properties
3. Configure wells and operating conditions
4. Run the simulation for the specified time period
5. Save results for visualization and analysis
"""

import numpy as np
from model import Model


def main():
    """
    Main function to execute the reservoir simulation.
    
    This simulation addresses your request for:
    Build a ccs model of a reservoir with porosity 0.25, perm 100 md and 2 wells
    """
    
    # Create model instance
    print("Initializing compositional reservoir model...")
    model = Model()
    
    # Configure physics engine
    print("Setting up physics...")
    model.set_physics()
    
    # Build reservoir grid and assign properties
    print("Creating reservoir grid (50x50x10 cells)...")
    model.set_reservoir()
    
    # Add wells to the model
    print("Configuring wells...")
    model.set_wells()
    
    # Set initial conditions
    print("Setting initial conditions...")
    model.set_initial_conditions()
    
    # Configure well controls
    print("Setting well operating conditions...")
    model.set_well_controls()
    
    # Set simulation parameters
    print("Configuring simulation parameters...")
    model.set_simulation_params(
        first_ts=0.001,    # Initial time step (days)
        max_ts=10.0,       # Maximum time step (days)
        final_time=365.0   # Total simulation time (days) - 1 year
    )
    
    # 
    # 
    
    # Run the simulation
    print("\nStarting simulation...")
    print("This may take a few minutes depending on model complexity...")
    model.run()
    
    # Save simulation results
    print("\nSaving results...")
    model.save_results("results")
    print("Results saved to 'results' directory")
    
    # Post-processing suggestions based on your objectives:
    
    
    
    
    

if __name__ == "__main__":
    main()
