"""Code Generator Agent for DARTSGPT.

This agent generates complete DARTS model code by:
- Loading selected template
- Filling template with extracted parameters
- Generating well configurations
- Creating main.py runner
- Adding necessary imports
"""

from typing import Dict, Any, List, Optional, Tuple
import json
import os
from pathlib import Path

from langchain.prompts import PromptTemplate
from langchain.tools import BaseTool, tool
from langchain.schema import HumanMessage, SystemMessage
from langchain_core.output_parsers import JsonOutputParser
from pydantic import BaseModel, Field

from ..agents.base_agent import BaseAgent, AgentInput, AgentOutput
from ..knowledge.templates.template_database import TEMPLATES


class CodeGeneration(BaseModel):
    """Code generation result."""
    model_code: str = Field(description="Generated model.py code")
    main_code: str = Field(description="Generated main.py code")
    imports: List[str] = Field(description="Required imports")
    well_configuration: Dict[str, Any] = Field(description="Well setup details")
    modifications: List[str] = Field(description="Modifications made to template")
    warnings: List[str] = Field(description="Any warnings or notes")
    reasoning: str = Field(description="Explanation of code generation")


class CodeGeneratorAgent(BaseAgent):
    """Agent that generates complete DARTS model code."""
    
    def __init__(self, llm=None):
        """Initialize the Code Generator Agent."""
        super().__init__(llm)
        self.output_parser = JsonOutputParser(pydantic_object=CodeGeneration)
        self.template_path = Path(__file__).parent.parent / "knowledge" / "templates" / "template_files"
    
    def get_tools(self) -> List[BaseTool]:
        """Get tools for code generation."""
        
        @tool
        def load_template_code(template_name: str) -> Dict[str, str]:
            """Load template code files.
            
            Args:
                template_name: Name of the template
                
            Returns:
                Dictionary with model_code and main_code
            """
            if template_name not in TEMPLATES:
                return {"error": f"Template {template_name} not found"}
            
            template_info = TEMPLATES[template_name]
            
            # For now, return a basic template structure
            # In production, this would load actual template files
            model_template = f"""# Generated from template: {template_name}
from darts.models.cicd_model import CICDModel
from darts.physics.{template_info.physics_type} import {template_info.physics_type.title().replace('_', '')}

class Model(CICDModel):
    \"\"\"Generated DARTS model for {template_info.description}\"\"\"
    
    def __init__(self):
        super().__init__()
        
        # Physics
        self.physics_type = "{template_info.physics_type}"
        
    def set_physics(self):
        \"\"\"Set up physics for the model.\"\"\"
        # Physics setup will be filled by code generator
        pass
        
    def set_reservoir(self):
        \"\"\"Set up reservoir properties.\"\"\"
        # Reservoir setup will be filled by code generator
        pass
        
    def set_wells(self):
        \"\"\"Set up well locations and perforations.\"\"\"
        # Well setup will be filled by code generator
        pass
        
    def set_initial_conditions(self):
        \"\"\"Set initial conditions.\"\"\"
        # Initial conditions will be filled by code generator
        pass
        
    def set_well_controls(self):
        \"\"\"Set well control schedules.\"\"\"
        # Well controls will be filled by code generator
        pass
"""

            main_template = f"""# Generated main.py for {template_name}
import numpy as np
from model import Model

def main():
    # Create model instance
    model = Model()
    
    # Initialize model
    model.init()
    
    # Run simulation
    model.run()
    
    # Save results
    model.save_results()
    
if __name__ == "__main__":
    main()
"""
            
            return {
                "model_code": model_template,
                "main_code": main_template
            }
        
        @tool
        def generate_physics_setup(physics_type: str, parameters: Dict[str, Any]) -> str:
            """Generate physics setup code.
            
            Args:
                physics_type: Type of physics
                parameters: Physics parameters
                
            Returns:
                Generated physics setup code
            """
            code_lines = []
            
            if physics_type == "compositional":
                code_lines.append("# Compositional physics setup")
                code_lines.append("self.n_components = 3")
                code_lines.append("self.phases = ['gas', 'oil']")
                code_lines.append("self.components = ['CO2', 'C1', 'H2O']")
                code_lines.append("")
                code_lines.append("# Flash calculation")
                code_lines.append("self.flash = ConstantK(self.components)")
                code_lines.append("self.flash.k_values = [3.0, 1.5, 0.01]")
                
            elif physics_type == "dead_oil":
                code_lines.append("# Dead oil physics setup")
                code_lines.append("self.n_components = 2")
                code_lines.append("self.phases = ['water', 'oil']")
                code_lines.append("self.components = ['H2O', 'oil']")
                code_lines.append("")
                code_lines.append("# Density models")
                code_lines.append("self.density_water = 1000  # kg/m3")
                code_lines.append("self.density_oil = 800    # kg/m3")
                
            elif physics_type == "geothermal":
                code_lines.append("# Geothermal physics setup")
                code_lines.append("self.thermal = True")
                code_lines.append("self.phases = ['water', 'steam']")
                code_lines.append("")
                code_lines.append("# Thermal properties")
                rock_heat_cap = parameters.get("rock_heat_capacity", 1000)
                rock_thermal_cond = parameters.get("rock_thermal_conductivity", 2.0)
                code_lines.append(f"self.rock_heat_capacity = {rock_heat_cap}  # J/kg/K")
                code_lines.append(f"self.rock_thermal_conductivity = {rock_thermal_cond}  # W/m/K")
                
            elif physics_type == "black_oil":
                code_lines.append("# Black oil physics setup")
                code_lines.append("self.n_components = 3")
                code_lines.append("self.phases = ['water', 'oil', 'gas']")
                code_lines.append("")
                code_lines.append("# PVT properties")
                code_lines.append("self.pvt_table = self.load_pvt_data('pvt_data.txt')")
                
            return "\n        ".join(code_lines)
        
        @tool
        def generate_reservoir_setup(grid_params: Dict[str, Any], 
                                   rock_props: Dict[str, Any]) -> str:
            """Generate reservoir setup code.
            
            Args:
                grid_params: Grid parameters
                rock_props: Rock properties
                
            Returns:
                Generated reservoir setup code
            """
            code_lines = []
            
            # Grid setup
            code_lines.append("# Grid dimensions")
            nx = grid_params.get("nx", 50)
            ny = grid_params.get("ny", 50)
            nz = grid_params.get("nz", 10)
            code_lines.append(f"self.nx = {nx}")
            code_lines.append(f"self.ny = {ny}")
            code_lines.append(f"self.nz = {nz}")
            code_lines.append("")
            
            # Cell dimensions
            dx = grid_params.get("dx", 10.0)
            dy = grid_params.get("dy", 10.0)
            dz = grid_params.get("dz", 2.0)
            code_lines.append("# Cell dimensions (m)")
            code_lines.append(f"self.dx = {dx}")
            code_lines.append(f"self.dy = {dy}")
            code_lines.append(f"self.dz = {dz}")
            code_lines.append("")
            
            # Rock properties
            code_lines.append("# Rock properties")
            porosity = rock_props.get("porosity", 0.2)
            permeability = rock_props.get("permeability", 100)  # Default 100 mD
            
            code_lines.append(f"self.porosity = np.ones((self.nx, self.ny, self.nz)) * {porosity}")
            code_lines.append(f"self.permeability = np.ones((self.nx, self.ny, self.nz)) * {permeability}  # mD")
            
            # Depth calculation if provided
            if "depth" in grid_params:
                depth = grid_params["depth"]
                code_lines.append("")
                code_lines.append(f"# Reservoir depth")
                code_lines.append(f"self.depth = {depth}  # m")
                code_lines.append(f"self.pressure_init = self.depth * 9.81 * 1000 + 101325  # Pa")
            
            return "\n        ".join(code_lines)
        
        @tool
        def generate_well_setup(well_params: Dict[str, Any], 
                              grid_dims: Tuple[int, int, int]) -> str:
            """Generate well setup code.
            
            Args:
                well_params: Well parameters
                grid_dims: Grid dimensions (nx, ny, nz)
                
            Returns:
                Generated well setup code
            """
            code_lines = []
            nx, ny, nz = grid_dims
            
            # Determine well pattern
            pattern = well_params.get("pattern", "single")
            n_injection = well_params.get("injection_wells", 1)
            n_production = well_params.get("production_wells", 1)
            
            code_lines.append("# Well locations")
            
            if pattern == "five-spot":
                # Five-spot pattern
                code_lines.append("# Five-spot pattern")
                code_lines.append("self.injection_wells = [")
                code_lines.append(f"    ('INJ1', (0, 0)),")
                code_lines.append(f"    ('INJ2', ({nx-1}, 0)),")
                code_lines.append(f"    ('INJ3', (0, {ny-1})),")
                code_lines.append(f"    ('INJ4', ({nx-1}, {ny-1}))")
                code_lines.append("]")
                code_lines.append(f"self.production_wells = [('PROD', ({nx//2}, {ny//2}))]")
                
            elif pattern == "line-drive":
                # Line drive pattern
                code_lines.append("# Line drive pattern")
                code_lines.append("self.injection_wells = [")
                for i in range(n_injection):
                    y_pos = int((i + 1) * ny / (n_injection + 1))
                    code_lines.append(f"    ('INJ{i+1}', (0, {y_pos})),")
                code_lines.append("]")
                code_lines.append("self.production_wells = [")
                for i in range(n_production):
                    y_pos = int((i + 1) * ny / (n_production + 1))
                    code_lines.append(f"    ('PROD{i+1}', ({nx-1}, {y_pos})),")
                code_lines.append("]")
                
            else:
                # Single injector/producer
                code_lines.append("# Single well pair")
                code_lines.append(f"self.injection_wells = [('INJ', (0, {ny//2}))]")
                code_lines.append(f"self.production_wells = [('PROD', ({nx-1}, {ny//2}))]")
            
            code_lines.append("")
            code_lines.append("# Add wells")
            code_lines.append("for name, (i, j) in self.injection_wells:")
            code_lines.append(f"    self.add_well(name, welltype='injector', i=i, j=j, k_range=(0, {nz-1}))")
            code_lines.append("")
            code_lines.append("for name, (i, j) in self.production_wells:")
            code_lines.append(f"    self.add_well(name, welltype='producer', i=i, j=j, k_range=(0, {nz-1}))")
            
            return "\n        ".join(code_lines)
        
        @tool
        def generate_well_controls(well_params: Dict[str, Any], 
                                 sim_params: Dict[str, Any],
                                 physics_type: str) -> str:
            """Generate well control code.
            
            Args:
                well_params: Well parameters
                sim_params: Simulation parameters
                physics_type: Physics type
                
            Returns:
                Generated well control code
            """
            code_lines = []
            
            total_time = sim_params.get("total_time", 365)
            
            code_lines.append("# Well controls")
            code_lines.append(f"self.simulation_time = {total_time}  # days")
            code_lines.append("")
            
            # Injection controls
            code_lines.append("# Injection wells")
            injection_rate = well_params.get("injection_rate", 100.0)  # m3/day
            
            if physics_type == "compositional":
                fluid = well_params.get("injection_fluid", "CO2")
                code_lines.append("for well_name, _ in self.injection_wells:")
                code_lines.append(f"    self.set_well_control(well_name, 'rate', {injection_rate})")
                code_lines.append(f"    self.set_injection_composition(well_name, '{fluid}', 1.0)")
            else:
                code_lines.append("for well_name, _ in self.injection_wells:")
                code_lines.append(f"    self.set_well_control(well_name, 'rate', {injection_rate})")
            
            code_lines.append("")
            
            # Production controls
            code_lines.append("# Production wells")
            bhp = well_params.get("production_bhp", 100e5)  # Pa
            code_lines.append("for well_name, _ in self.production_wells:")
            code_lines.append(f"    self.set_well_control(well_name, 'bhp', {bhp})")
            
            return "\n        ".join(code_lines)
        
        return [
            load_template_code,
            generate_physics_setup,
            generate_reservoir_setup,
            generate_well_setup,
            generate_well_controls
        ]
    
    def get_prompt_template(self) -> PromptTemplate:
        """Get the prompt template for code generation."""
        template = """You are an expert DARTS code generator.
        
Generate complete DARTS model code based on the template and parameters.

Selected Template: {template_info}
Extracted Parameters: {parameters}
Generated Code Sections:
- Physics Setup: {physics_code}
- Reservoir Setup: {reservoir_code}  
- Well Setup: {well_code}
- Well Controls: {control_code}

Template Code:
{template_code}

Instructions:
1. Fill the template with the generated code sections
2. Ensure all imports are correct
3. Add necessary helper methods
4. Make the code production-ready
5. Include comments explaining key sections
6. Generate both model.py and main.py
7. List any warnings or assumptions

{format_instructions}
"""
        
        return PromptTemplate(
            template=template,
            input_variables=[
                "template_info",
                "parameters",
                "physics_code",
                "reservoir_code",
                "well_code",
                "control_code",
                "template_code"
            ],
            partial_variables={
                "format_instructions": self.output_parser.get_format_instructions()
            }
        )
    
    async def process(self, input_data: AgentInput) -> AgentOutput:
        """Process parameters and generate DARTS code.
        
        Args:
            input_data: Input containing extracted parameters
            
        Returns:
            AgentOutput with generated code
        """
        if not self.validate_input(input_data):
            return AgentOutput(
                result={"error": "Invalid input"},
                reasoning="Empty or invalid input provided",
                confidence=0.0
            )
        
        # Get context from previous agents
        template_info = input_data.context.get("template_selection", {})
        parameters = input_data.context.get("parameter_extraction", {})
        template_name = template_info.get("template_name", "2ph_do")
        
        # Use tools to generate code sections
        tools = self.get_tools()
        load_tool = tools[0]
        physics_tool = tools[1]
        reservoir_tool = tools[2]
        well_tool = tools[3]
        control_tool = tools[4]
        
        # Load template
        template_code = load_tool.invoke({"template_name": template_name})
        
        # Generate code sections
        physics_code = physics_tool.invoke({
            "physics_type": TEMPLATES[template_name].physics_type,
            "parameters": parameters.get("custom_parameters", {})
        })
        
        reservoir_code = reservoir_tool.invoke({
            "grid_params": parameters.get("grid_parameters", {}),
            "rock_props": parameters.get("rock_properties", {})
        })
        
        grid_params = parameters.get("grid_parameters", {})
        grid_dims = (
            grid_params.get("nx", 50),
            grid_params.get("ny", 50),
            grid_params.get("nz", 10)
        )
        
        well_code = well_tool.invoke({
            "well_params": parameters.get("well_parameters", {}),
            "grid_dims": grid_dims
        })
        
        control_code = control_tool.invoke({
            "well_params": parameters.get("well_parameters", {}),
            "sim_params": parameters.get("simulation_parameters", {}),
            "physics_type": TEMPLATES[template_name].physics_type
        })
        
        # Prepare prompt
        prompt = self.get_prompt_template()
        formatted_prompt = prompt.format(
            template_info=json.dumps(template_info, indent=2),
            parameters=json.dumps(parameters, indent=2),
            physics_code=physics_code,
            reservoir_code=reservoir_code,
            well_code=well_code,
            control_code=control_code,
            template_code=json.dumps(template_code, indent=2)
        )
        
        # Get LLM response
        messages = [
            SystemMessage(content="You are an expert DARTS code generator."),
            HumanMessage(content=formatted_prompt)
        ]
        
        response = await self.llm.ainvoke(messages)
        
        # Parse response
        try:
            generation = self.output_parser.parse(response.content)
            
            # Update context
            updated_context = input_data.context.copy()
            updated_context["code_generation"] = generation.model_dump()
            
            return AgentOutput(
                result=generation.model_dump(),
                reasoning=generation.reasoning,
                confidence=0.9,
                next_agent="Validator"
            )
        except Exception as e:
            # Generate basic code as fallback
            model_code = template_code["model_code"].replace(
                "# Physics setup will be filled by code generator",
                physics_code
            ).replace(
                "# Reservoir setup will be filled by code generator",
                reservoir_code
            ).replace(
                "# Well setup will be filled by code generator",
                well_code
            ).replace(
                "# Well controls will be filled by code generator",
                control_code
            )
            
            fallback_generation = {
                "model_code": model_code,
                "main_code": template_code["main_code"],
                "imports": ["numpy", "darts"],
                "well_configuration": parameters.get("well_parameters", {}),
                "modifications": ["Filled template with generated sections"],
                "warnings": [f"Fallback generation: {str(e)}"],
                "reasoning": "Generated code using template filling approach"
            }
            
            updated_context = input_data.context.copy()
            updated_context["code_generation"] = fallback_generation
            
            return AgentOutput(
                result=fallback_generation,
                reasoning=fallback_generation["reasoning"],
                confidence=0.7,
                next_agent="Validator"
            )