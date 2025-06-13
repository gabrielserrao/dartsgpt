"""Multi-agent orchestrator for DARTSGPT using LangGraph."""

import operator
from typing import Annotated, List, Literal, TypedDict, Dict, Any
from pathlib import Path
import json

from langchain_core.messages import BaseMessage, HumanMessage, SystemMessage
from langchain_openai import ChatOpenAI
from langgraph.graph import StateGraph, MessagesState, START, END
from langgraph.types import Command
from langgraph.prebuilt import create_react_agent
from langgraph.checkpoint.memory import MemorySaver

from config import get_settings
from knowledge.templates.template_database import TEMPLATES
from knowledge.physics_features import PHYSICS_KEYWORDS, UNIT_CONVERSIONS
from knowledge.example_prompts import EXAMPLE_PROMPTS


# Define the main state for our orchestrator
class DARTSGPTState(MessagesState):
    """State for the DARTSGPT multi-agent system."""
    current_agent: str
    prompt: str
    intent: Dict[str, Any]
    selected_template: str
    parameters: Dict[str, Any]
    model_code: str
    main_code: str
    validation_result: Dict[str, Any]
    final_output: Dict[str, Any]


# Intent Classifier Agent
def create_intent_classifier(model: ChatOpenAI) -> callable:
    """Create the intent classification agent."""
    
    def intent_classifier(state: DARTSGPTState) -> Command[Literal["template_selector", "parameter_extractor", END]]:
        """Classify the user's intent and extract physics type."""
        
        system_prompt = """You are an expert at understanding DARTS reservoir simulation requests.
        Analyze the user's prompt and extract:
        
        1. Physics type - Choose the MOST appropriate:
           - compositional: CO2, gas injection, EOS, miscible displacement, carbon storage
           - dead_oil: waterflood, immiscible, basic oil-water, steam injection
           - black_oil: PVT, solution gas, bubble point, gas cap, volatile oil
           - geothermal: temperature, heat, steam, thermal, hot water, EGS
           - poroelastic: geomechanics, stress, compaction, subsidence, coupled
           - chemical: surfactant, polymer, ASP, EOR chemicals, foam
        
        2. Key features (list all that apply):
           - Flow type: injection, production, natural flow
           - Well types: vertical, horizontal, multilateral, smart
           - Special: thermal effects, reactions, tracers, monitoring
        
        3. Complexity level:
           - simple: basic physics, few wells, homogeneous
           - moderate: heterogeneous, multiple wells, time-varying
           - complex: coupled physics, complex geology, optimization
        
        Return ONLY a JSON object with fields: physics_type, features (array), complexity"""
        
        messages = [
            SystemMessage(content=system_prompt),
            HumanMessage(content=f"Analyze this prompt: {state['prompt']}")
        ]
        
        response = model.invoke(messages)
        
        # Parse response (in production, use structured output)
        try:
            intent = json.loads(response.content)
        except:
            # Fallback to simple keyword matching
            prompt_lower = state['prompt'].lower()
            physics_type = "dead_oil"
            for ptype, keywords in PHYSICS_KEYWORDS.items():
                if any(kw in prompt_lower for kw in keywords):
                    physics_type = ptype
                    break
            
            intent = {
                "physics_type": physics_type,
                "features": [],
                "complexity": "moderate"
            }
        
        print(f"🔍 Intent classified: {intent}")
        
        return Command(
            goto="template_selector",
            update={
                "intent": intent,
                "current_agent": "template_selector"
            }
        )
    
    return intent_classifier


# Template Selector Agent
def create_template_selector(model: ChatOpenAI) -> callable:
    """Create the template selection agent."""
    
    def template_selector(state: DARTSGPTState) -> Command[Literal["parameter_extractor", END]]:
        """Select the best DARTS template based on intent."""
        
        physics_type = state['intent'].get('physics_type', 'dead_oil')
        features = state['intent'].get('features', [])
        
        # Find matching templates
        candidates = [
            (name, template) for name, template in TEMPLATES.items()
            if template.physics_type == physics_type
        ]
        
        if not candidates:
            candidates = [
                (name, template) for name, template in TEMPLATES.items()
                if template.physics_type == "dead_oil"
            ]
        
        # Prefer golden templates
        golden = [(n, t) for n, t in candidates if t.is_golden]
        if golden:
            selected = golden[0][0]
        else:
            selected = candidates[0][0] if candidates else "2ph_do"
        
        print(f"📋 Template selected: {selected}")
        
        return Command(
            goto="parameter_extractor",
            update={
                "selected_template": selected,
                "current_agent": "parameter_extractor"
            }
        )
    
    return template_selector


# Parameter Extractor Agent
def create_parameter_extractor(model: ChatOpenAI) -> callable:
    """Create the parameter extraction agent."""
    
    def parameter_extractor(state: DARTSGPTState) -> Command[Literal["code_generator", END]]:
        """Extract numerical parameters from the prompt."""
        
        system_prompt = """You are an expert at extracting numerical parameters from reservoir descriptions.
        
        Extract ALL numerical values and their units:
        
        1. Grid/Geometry:
           - Grid dimensions: nx, ny, nz (e.g., "100x100x20" → nx:100, ny:100, nz:20)
           - Cell sizes: dx, dy, dz in meters or feet
           - Domain size: length, width, height
           - Depth: reservoir depth below surface
        
        2. Rock Properties:
           - Porosity: ALWAYS convert to fraction (0-1)
             * "25%" or "25 percent" → 0.25
             * "0.25" → 0.25
           - Permeability: extract value and unit
             * "100 mD" or "100 millidarcy" → 100 (unit: mD)
             * "1e-13 m²" → 1e-13 (unit: m2)
           - Rock compressibility, thermal properties
        
        3. Fluid Properties:
           - Temperature: value and unit (C, F, K)
           - Pressure: value and unit (bar, psi, Pa, atm)
           - Viscosity, density, saturation
        
        4. Well Parameters:
           - Rates: injection/production rates with units (m³/day, bbl/day, kg/s)
           - Pressures: BHP, wellhead pressure
           - Well depth, radius, skin factor
        
        5. Time Parameters:
           - Simulation time, time steps, cycles
        
        Return a structured JSON with sections:
        {
          "grid_parameters": {...},
          "rock_properties": {...},
          "fluid_properties": {...},
          "well_parameters": {...},
          "time_parameters": {...}
        }
        
        Include units in a separate field when applicable."""
        
        messages = [
            SystemMessage(content=system_prompt),
            HumanMessage(content=f"Extract parameters from: {state['prompt']}")
        ]
        
        response = model.invoke(messages)
        
        # Parse response (in production, use structured output)
        try:
            parameters = json.loads(response.content)
            # Post-process porosity if needed
            if "rock_properties" in parameters:
                if "porosity" in parameters["rock_properties"]:
                    poro = parameters["rock_properties"]["porosity"]
                    # If porosity > 1, assume it's a percentage and convert
                    if isinstance(poro, (int, float)) and poro > 1:
                        parameters["rock_properties"]["porosity"] = poro / 100.0
            
            # Merge fluid_properties into main parameters for backward compatibility
            if "fluid_properties" in parameters:
                for key, value in parameters["fluid_properties"].items():
                    parameters[key] = value
        except:
            # Fallback to regex extraction
            import re
            parameters = {"grid_parameters": {}, "rock_properties": {}}
            
            # Grid dimensions
            grid_pattern = r'(\d+)\s*[xX×]\s*(\d+)\s*[xX×]\s*(\d+)'
            grid_match = re.search(grid_pattern, state['prompt'])
            if grid_match:
                parameters["grid_parameters"]["nx"] = int(grid_match.group(1))
                parameters["grid_parameters"]["ny"] = int(grid_match.group(2))
                parameters["grid_parameters"]["nz"] = int(grid_match.group(3))
            
            # Porosity extraction
            poro_pattern = r'(\d+(?:\.\d+)?)\s*%\s*porosity'
            poro_match = re.search(poro_pattern, state['prompt'], re.IGNORECASE)
            if poro_match:
                parameters["rock_properties"]["porosity"] = float(poro_match.group(1)) / 100.0
            
            # Permeability extraction
            perm_pattern = r'(\d+(?:\.\d+)?)\s*m[Dd]'
            perm_match = re.search(perm_pattern, state['prompt'], re.IGNORECASE)
            if perm_match:
                parameters["rock_properties"]["permeability"] = float(perm_match.group(1))
        
        print(f"🔢 Parameters extracted: {parameters}")
        
        return Command(
            goto="code_generator",
            update={
                "parameters": parameters,
                "current_agent": "code_generator"
            }
        )
    
    return parameter_extractor


# Code Generator Agent
def create_code_generator(model: ChatOpenAI) -> callable:
    """Create the code generation agent."""
    
    def code_generator(state: DARTSGPTState) -> Command[Literal["validator", END]]:
        """Generate DARTS model code based on template and parameters."""
        
        template_name = state['selected_template']
        template = TEMPLATES[template_name]
        parameters = state.get('parameters', {})
        prompt = state['prompt']
        intent = state.get('intent', {})
        
        # Get parameters with defaults, handling None values
        grid = parameters.get('grid_parameters', {})
        nx = grid.get('nx') if grid.get('nx') is not None else 50
        ny = grid.get('ny') if grid.get('ny') is not None else 50
        nz = grid.get('nz') if grid.get('nz') is not None else 10
        
        rock = parameters.get('rock_properties', {})
        poro = rock.get('porosity') if rock.get('porosity') is not None else 0.2
        perm_md = rock.get('permeability') if rock.get('permeability') is not None else 100.0  # in mD
        # Convert permeability from mD to m²
        perm = perm_md * 9.869233e-16  # m²
        
        well_params = parameters.get('well_parameters', {})
        inj_rate = well_params.get('injection_rate') if well_params.get('injection_rate') is not None else 100.0
        prod_bhp = well_params.get('production_pressure') if well_params.get('production_pressure') is not None else 180.0
        
        # Check for temperature in multiple places
        temp = parameters.get('temperature')
        if temp is None and 'fluid_properties' in parameters:
            temp = parameters['fluid_properties'].get('temperature')
        if temp is None:
            temp = 50.0
            
        # Check for pressure in multiple places  
        pressure = parameters.get('pressure')
        if pressure is None and 'fluid_properties' in parameters:
            pressure = parameters['fluid_properties'].get('pressure')
        if pressure is None:
            pressure = 200.0
        
        # Create comments explaining the model based on user prompt
        physics_comment = f"Physics type '{template.physics_type}' was selected based on keywords in your prompt"
        if 'co2' in prompt.lower():
            physics_comment += " (CO2 injection implies compositional physics)"
        elif 'geothermal' in prompt.lower():
            physics_comment += " (geothermal keywords detected)"
        elif 'waterflood' in prompt.lower():
            physics_comment += " (waterflood implies dead oil or black oil physics)"
            
        grid_comment = "Grid dimensions "
        if nx == 50 and ny == 50 and nz == 10:
            grid_comment += "using default values (not specified in prompt)"
        else:
            grid_comment += f"extracted from prompt: {nx}x{ny}x{nz}"
            
        # Generate model code with detailed comments
        model_code = f'''"""
DARTS Model Generated from User Prompt
=====================================

Original prompt: "{prompt}"

Model Configuration:
- Physics Type: {template.physics_type} ({template.description})
- Template Used: {template_name}
- Intent Classification: {intent}

This model was automatically generated based on your natural language description.
The code includes explanatory comments showing how your requirements were interpreted.
"""

import numpy as np
from darts.models.reservoirmodel import ReservoirModel
from darts.physics.{template.physics_type} import {template.physics_type.title().replace("_", "")}
from darts.tools.keyword_file_tools import load_single_keyword


class Model(ReservoirModel):
    """
    {template.description}
    
    This model implements {template.physics_type.replace('_', ' ')} physics
    as requested in your prompt. The specific physics type was chosen because:
    {physics_comment}
    """
    
    def __init__(self):
        """Initialize the DARTS reservoir model."""
        super().__init__()
        self.physics_type = "{template.physics_type}"
        # This physics type handles {', '.join(intent.get('features', ['standard flow']))} features
        
    def set_physics(self):
        """
        Set up the physics engine for the simulation.
        
        Using {template.physics_type.title().replace("_", "")} physics which provides:
        - {'Compositional flow with EOS' if 'comp' in template.physics_type else ''}
        - {'Thermal effects' if 'thermal' in template.physics_type or 'geothermal' in template.physics_type else ''}
        - {'Black oil PVT relations' if 'black_oil' in template.physics_type else ''}
        - {'Geomechanics coupling' if 'poroelastic' in template.physics_type else ''}
        - {'Chemical species transport' if 'chemical' in template.physics_type else ''}
        """
        self.physics = {template.physics_type.title().replace("_", "")}()
        
    def set_reservoir(self):
        """
        Define reservoir geometry and rock properties.
        
        {grid_comment}
        """
        # Grid dimensions
        self.nx = {nx}  # Number of cells in X direction
        self.ny = {ny}  # Number of cells in Y direction  
        self.nz = {nz}  # Number of cells in Z direction (layers)
        
        # Cell dimensions (m)
        self.dx = 10.0  # Cell size in X direction
        self.dy = 10.0  # Cell size in Y direction
        self.dz = 2.0   # Cell thickness
        
        # Total reservoir dimensions: {nx*10}m x {ny*10}m x {nz*2}m
        
        # Rock properties
        # {'Porosity extracted from prompt' if poro != 0.2 else 'Using default porosity (not specified)'}
        self.porosity = np.ones((self.nx, self.ny, self.nz)) * {poro}
        
        # {'Permeability extracted from prompt' if perm != 100e-15 else 'Using default permeability (not specified)'}
        self.permeability = np.ones((self.nx, self.ny, self.nz)) * {perm}  # m²
        # Note: {perm} m² = {perm/9.869233e-16:.1f} mD
        
        # Initialize the reservoir grid with these properties
        self.init_reservoir()
        
    def set_wells(self):
        """
        Configure well locations and completions.
        
        {'Wells configured based on injection/production keywords in prompt' if any(kw in prompt.lower() for kw in ['inject', 'produc']) else 'Using default well configuration'}
        """
        # Injector well
        # {'Injection well requested in prompt' if 'inject' in prompt.lower() else 'Default injector added'}
        self.add_well(
            "INJ1",                      # Well name
            welltype='injector',         # Well type
            i=1, j=1,                   # Location at corner of reservoir
            k_top=1, k_bottom=self.nz   # Completed through all layers
        )
        
        # Producer well  
        # {'Production well requested in prompt' if 'produc' in prompt.lower() else 'Default producer added'}
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
        
        {'Temperature specified in prompt: ' + str(temp) + '°C' if temp != 50.0 else 'Using default temperature'}
        {'Pressure specified in prompt: ' + str(pressure) + ' bar' if pressure != 200.0 else 'Using default pressure'}
        """
        self.set_uniform_initial_conditions(
            pressure={pressure},    # Initial pressure in bar
            temperature={temp}      # Initial temperature in °C
        )
        
        # {'Special initial conditions for ' + template.physics_type if template.physics_type in ['compositional', 'chemical'] else ''}
        
    def set_well_controls(self):
        """
        Define well operating conditions and control modes.
        
        Well controls are set based on the simulation objectives:
        {'- Injection for CO2 storage' if 'co2' in prompt.lower() else ''}
        {'- Water injection for pressure support' if 'waterflood' in prompt.lower() else ''}
        {'- Production for reservoir depletion' if 'depletion' in prompt.lower() else ''}
        """
        # Injector operating at constant rate
        # {'Injection rate extracted from prompt' if inj_rate != 100.0 else 'Using default injection rate'}
        self.set_well_rate(
            "INJ1", 
            rate={inj_rate}  # m³/day {'(CO2 at surface conditions)' if 'co2' in prompt.lower() else '(water)' if 'water' in prompt.lower() else ''}
        )
        
        # Producer operating at constant bottom hole pressure
        # {'Production pressure extracted from prompt' if prod_bhp != 180.0 else 'Using default production pressure'}
        self.set_well_bhp(
            "PROD1", 
            bhp={prod_bhp}   # bar (bottom hole pressure)
        )
'''

        # Generate main code with comments
        main_code = f'''"""
Main Simulation Script
=====================

This script runs the DARTS reservoir simulation model that was generated
from your prompt: "{prompt}"

The simulation will:
1. Initialize the {template.physics_type.replace('_', ' ')} model
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
    {prompt}
    """
    
    # Create model instance
    print("Initializing {template.physics_type.replace('_', ' ')} reservoir model...")
    model = Model()
    
    # Configure physics engine
    print("Setting up physics...")
    model.set_physics()
    
    # Build reservoir grid and assign properties
    print("Creating reservoir grid ({nx}x{ny}x{nz} cells)...")
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
    
    # {'Adjust simulation time for CO2 storage (typically longer)' if 'co2' in prompt.lower() and 'storage' in prompt.lower() else ''}
    # {'Consider longer simulation for EOR processes' if 'eor' in prompt.lower() else ''}
    
    # Run the simulation
    print("\\nStarting simulation...")
    print("This may take a few minutes depending on model complexity...")
    model.run()
    
    # Save simulation results
    print("\\nSaving results...")
    model.save_results("results")
    print("Results saved to 'results' directory")
    
    # Post-processing suggestions based on your objectives:
    {'# - Analyze CO2 plume migration and trapping' if 'co2' in prompt.lower() else ''}
    {'# - Evaluate thermal breakthrough times' if 'geothermal' in prompt.lower() else ''}
    {'# - Calculate oil recovery factor' if any(kw in prompt.lower() for kw in ['oil', 'recovery']) else ''}
    {'# - Monitor pressure buildup/depletion' if 'pressure' in prompt.lower() else ''}
    

if __name__ == "__main__":
    main()
'''
        
        print(f"💻 Code generated with detailed comments: {len(model_code)} + {len(main_code)} characters")
        
        return Command(
            goto="validator",
            update={
                "model_code": model_code,
                "main_code": main_code,
                "current_agent": "validator"
            }
        )
    
    return code_generator


# Validator Agent
def create_validator(model: ChatOpenAI) -> callable:
    """Create the validation agent."""
    
    def validator(state: DARTSGPTState) -> Command[Literal[END]]:
        """Validate the generated code."""
        
        validation_result = {
            "syntax_valid": True,
            "imports_valid": True,
            "structure_valid": True,
            "issues": []
        }
        
        # Basic validation checks
        model_code = state['model_code']
        
        # Check for required imports
        if "from darts.models" not in model_code:
            validation_result["imports_valid"] = False
            validation_result["issues"].append("Missing DARTS model import")
        
        # Check for required methods
        required_methods = ["set_physics", "set_reservoir", "set_wells"]
        for method in required_methods:
            if f"def {method}" not in model_code:
                validation_result["structure_valid"] = False
                validation_result["issues"].append(f"Missing required method: {method}")
        
        validation_passed = all([
            validation_result["syntax_valid"],
            validation_result["imports_valid"],
            validation_result["structure_valid"]
        ])
        
        print(f"✅ Validation {'passed' if validation_passed else 'failed'}")
        
        # Prepare final output
        final_output = {
            "success": True,
            "model_code": state['model_code'],
            "main_code": state['main_code'],
            "template_used": state['selected_template'],
            "parameters": state['parameters'],
            "validation": validation_result
        }
        
        return Command(
            goto=END,
            update={
                "validation_result": validation_result,
                "final_output": final_output,
                "current_agent": "completed"
            }
        )
    
    return validator


# Supervisor Agent
def create_supervisor(model: ChatOpenAI) -> callable:
    """Create the supervisor agent that orchestrates the workflow."""
    
    def supervisor(state: DARTSGPTState) -> Command[Literal["intent_classifier", "template_selector", "parameter_extractor", "code_generator", "validator", END]]:
        """Supervise the multi-agent workflow."""
        
        # Initial routing - start with intent classification
        if not state.get('intent'):
            print("\n🎯 Supervisor: Starting with intent classification")
            return Command(goto="intent_classifier", update={"current_agent": "intent_classifier"})
        
        # Check current state and route accordingly
        current = state.get('current_agent', '')
        
        if current == "completed":
            print("\n🎯 Supervisor: Workflow completed successfully!")
            return Command(goto=END)
        
        # Let the workflow continue
        print(f"\n🎯 Supervisor: Current agent is {current}")
        return Command(goto=current)
    
    return supervisor


def build_dartsgpt_graph(model: ChatOpenAI = None) -> StateGraph:
    """Build the complete DARTSGPT multi-agent graph."""
    
    if model is None:
        settings = get_settings()
        model = ChatOpenAI(
            model=settings.openai_model,
            api_key=settings.openai_api_key,
            temperature=0.7
        )
    
    # Initialize the graph
    builder = StateGraph(DARTSGPTState)
    
    # Create agents
    supervisor_agent = create_supervisor(model)
    intent_agent = create_intent_classifier(model)
    template_agent = create_template_selector(model)
    parameter_agent = create_parameter_extractor(model)
    code_agent = create_code_generator(model)
    validator_agent = create_validator(model)
    
    # Add nodes
    builder.add_node("supervisor", supervisor_agent)
    builder.add_node("intent_classifier", intent_agent)
    builder.add_node("template_selector", template_agent)
    builder.add_node("parameter_extractor", parameter_agent)
    builder.add_node("code_generator", code_agent)
    builder.add_node("validator", validator_agent)
    
    # Set entry point
    builder.add_edge(START, "supervisor")
    
    # Add edges from agents back to supervisor
    builder.add_edge("intent_classifier", "template_selector")
    builder.add_edge("template_selector", "parameter_extractor")
    builder.add_edge("parameter_extractor", "code_generator")
    builder.add_edge("code_generator", "validator")
    builder.add_edge("validator", END)
    
    # Compile the graph
    return builder.compile()


def run_dartsgpt_orchestrator(prompt: str, output_dir: str = "output") -> Dict[str, Any]:
    """Run the DARTSGPT orchestrator with a given prompt."""
    
    print(f"\n{'='*60}")
    print(f"DARTSGPT Multi-Agent Orchestrator")
    print(f"{'='*60}")
    print(f"Prompt: {prompt}")
    print(f"{'='*60}\n")
    
    # Build the graph
    graph = build_dartsgpt_graph()
    
    # Initialize state
    initial_state = {
        "messages": [HumanMessage(content=prompt)],
        "prompt": prompt
    }
    
    # Run the graph
    result = graph.invoke(initial_state)
    
    # Extract final output
    final_output = result.get('final_output', {})
    
    if final_output.get('success') and output_dir is not None:
        # Save generated files
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        (output_path / "model.py").write_text(final_output['model_code'])
        (output_path / "main.py").write_text(final_output['main_code'])
        (output_path / "metadata.json").write_text(json.dumps({
            "prompt": prompt,
            "template_used": final_output['template_used'],
            "parameters": final_output['parameters'],
            "validation": final_output['validation']
        }, indent=2))
        
        print(f"\n✅ Files saved to: {output_path}/")
        print(f"   - model.py")
        print(f"   - main.py")
        print(f"   - metadata.json")
    
    return final_output


if __name__ == "__main__":
    # Test the orchestrator
    test_prompt = "Create a CO2 injection model for carbon storage in a 100x100x20 reservoir with 25% porosity"
    result = run_dartsgpt_orchestrator(test_prompt)