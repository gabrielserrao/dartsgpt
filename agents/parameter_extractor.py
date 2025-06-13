"""Parameter Extractor Agent for DARTSGPT.

This agent extracts numerical parameters from natural language and:
- Identifies parameter values with units
- Converts units to DARTS standard (SI)
- Applies template-specific defaults
- Validates parameter ranges
"""

from typing import Dict, Any, List, Optional, Tuple, Union
import re
import json
from dataclasses import dataclass

from langchain.prompts import PromptTemplate
from langchain.tools import BaseTool, tool
from langchain.schema import HumanMessage, SystemMessage
from langchain_core.output_parsers import JsonOutputParser
from pydantic import BaseModel, Field

from ..agents.base_agent import BaseAgent, AgentInput, AgentOutput
from ..knowledge.templates.template_database import TEMPLATES
from ..knowledge.physics_features import UNIT_CONVERSIONS


@dataclass
class ExtractedParameter:
    """Extracted parameter with value and unit."""
    name: str
    value: Union[float, int, str]
    unit: Optional[str] = None
    original_text: str = ""
    confidence: float = 1.0


class ParameterExtraction(BaseModel):
    """Parameter extraction result."""
    grid_parameters: Dict[str, Any] = Field(description="Grid-related parameters")
    rock_properties: Dict[str, Any] = Field(description="Rock properties (porosity, permeability)")
    fluid_properties: Dict[str, Any] = Field(description="Fluid properties")
    well_parameters: Dict[str, Any] = Field(description="Well configuration parameters")
    simulation_parameters: Dict[str, Any] = Field(description="Simulation control parameters")
    custom_parameters: Dict[str, Any] = Field(description="Template-specific parameters")
    units_used: Dict[str, str] = Field(description="Units for each parameter")
    reasoning: str = Field(description="Explanation of extraction and conversions")


class ParameterExtractorAgent(BaseAgent):
    """Agent that extracts and converts parameters from user input."""
    
    def __init__(self, llm=None):
        """Initialize the Parameter Extractor Agent."""
        super().__init__(llm)
        self.output_parser = JsonOutputParser(pydantic_object=ParameterExtraction)
        
        # Common parameter patterns
        self.patterns = {
            # Grid dimensions
            "grid_3d": re.compile(r'(\d+)\s*[xX×]\s*(\d+)\s*[xX×]\s*(\d+)'),
            "grid_2d": re.compile(r'(\d+)\s*[xX×]\s*(\d+)(?!\s*[xX×])'),
            
            # Cell dimensions
            "cell_dims": re.compile(r'(\d+(?:\.\d+)?)\s*m\s*[xX×]\s*(\d+(?:\.\d+)?)\s*m\s*[xX×]\s*(\d+(?:\.\d+)?)\s*m'),
            
            # Depth
            "depth": re.compile(r'(\d+(?:\.\d+)?)\s*(?:m|meters?|ft|feet)\s*depth', re.I),
            
            # Porosity
            "porosity_percent": re.compile(r'(\d+(?:\.\d+)?)\s*%\s*porosity', re.I),
            "porosity_fraction": re.compile(r'porosity\s*[:=]?\s*(\d+(?:\.\d+)?)', re.I),
            
            # Permeability
            "permeability": re.compile(r'(\d+(?:\.\d+)?)\s*(?:md|mD|millidarcy|darcy)', re.I),
            
            # Temperature
            "temperature": re.compile(r'(\d+(?:\.\d+)?)\s*°?(?:C|F|celsius|fahrenheit|K|kelvin)', re.I),
            
            # Pressure
            "pressure": re.compile(r'(\d+(?:\.\d+)?)\s*(?:psi|bar|pa|pascal|atm)', re.I),
            
            # Rate
            "rate": re.compile(r'(\d+(?:\.\d+)?)\s*(?:m3/day|bbl/day|kg/s|stb/d)', re.I),
            
            # Time
            "time": re.compile(r'(\d+(?:\.\d+)?)\s*(?:days?|years?|months?)', re.I),
            
            # Well count
            "well_count": re.compile(r'(\d+)\s+(?:injection|production|injector|producer|well)s?', re.I)
        }
    
    def get_tools(self) -> List[BaseTool]:
        """Get tools for parameter extraction."""
        
        @tool
        def extract_numerical_values(text: str) -> List[ExtractedParameter]:
            """Extract numerical values with units from text.
            
            Args:
                text: Input text to analyze
                
            Returns:
                List of extracted parameters
            """
            extracted = []
            
            # Grid dimensions
            match = self.patterns["grid_3d"].search(text)
            if match:
                extracted.extend([
                    ExtractedParameter("nx", int(match.group(1)), None, match.group(0)),
                    ExtractedParameter("ny", int(match.group(2)), None, match.group(0)),
                    ExtractedParameter("nz", int(match.group(3)), None, match.group(0))
                ])
            else:
                match = self.patterns["grid_2d"].search(text)
                if match:
                    extracted.extend([
                        ExtractedParameter("nx", int(match.group(1)), None, match.group(0)),
                        ExtractedParameter("ny", int(match.group(2)), None, match.group(0)),
                        ExtractedParameter("nz", 1, None, "2D grid")
                    ])
            
            # Cell dimensions
            match = self.patterns["cell_dims"].search(text)
            if match:
                extracted.extend([
                    ExtractedParameter("dx", float(match.group(1)), "m", match.group(0)),
                    ExtractedParameter("dy", float(match.group(2)), "m", match.group(0)),
                    ExtractedParameter("dz", float(match.group(3)), "m", match.group(0))
                ])
            
            # Depth
            for match in re.finditer(self.patterns["depth"], text):
                unit = "ft" if "ft" in match.group(0).lower() or "feet" in match.group(0).lower() else "m"
                extracted.append(
                    ExtractedParameter("depth", float(match.group(1)), unit, match.group(0))
                )
            
            # Porosity
            for match in re.finditer(self.patterns["porosity_percent"], text):
                extracted.append(
                    ExtractedParameter("porosity", float(match.group(1)) / 100, "fraction", match.group(0))
                )
            for match in re.finditer(self.patterns["porosity_fraction"], text):
                value = float(match.group(1))
                if value > 1:  # Likely percentage
                    value = value / 100
                extracted.append(
                    ExtractedParameter("porosity", value, "fraction", match.group(0))
                )
            
            # Permeability
            for match in re.finditer(self.patterns["permeability"], text):
                unit = "md" if "d" in match.group(0).lower() else "darcy"
                extracted.append(
                    ExtractedParameter("permeability", float(match.group(1)), unit, match.group(0))
                )
            
            # Temperature
            for match in re.finditer(self.patterns["temperature"], text):
                text_lower = match.group(0).lower()
                if "f" in text_lower or "fahrenheit" in text_lower:
                    unit = "F"
                elif "k" in text_lower or "kelvin" in text_lower:
                    unit = "K"
                else:
                    unit = "C"
                extracted.append(
                    ExtractedParameter("temperature", float(match.group(1)), unit, match.group(0))
                )
            
            # Pressure
            for match in re.finditer(self.patterns["pressure"], text):
                text_lower = match.group(0).lower()
                if "psi" in text_lower:
                    unit = "psi"
                elif "bar" in text_lower:
                    unit = "bar"
                elif "atm" in text_lower:
                    unit = "atm"
                else:
                    unit = "pa"
                extracted.append(
                    ExtractedParameter("pressure", float(match.group(1)), unit, match.group(0))
                )
            
            # Rate
            for match in re.finditer(self.patterns["rate"], text):
                text_lower = match.group(0).lower()
                if "m3/day" in text_lower:
                    unit = "m3/day"
                elif "bbl/day" in text_lower or "stb/d" in text_lower:
                    unit = "bbl/day"
                elif "kg/s" in text_lower:
                    unit = "kg/s"
                else:
                    unit = "m3/day"
                extracted.append(
                    ExtractedParameter("rate", float(match.group(1)), unit, match.group(0))
                )
            
            # Time
            for match in re.finditer(self.patterns["time"], text):
                text_lower = match.group(0).lower()
                if "year" in text_lower:
                    unit = "years"
                elif "month" in text_lower:
                    unit = "months"
                else:
                    unit = "days"
                extracted.append(
                    ExtractedParameter("time", float(match.group(1)), unit, match.group(0))
                )
            
            # Well count
            for match in re.finditer(self.patterns["well_count"], text):
                well_type = "injection" if "inj" in match.group(0).lower() else "production"
                extracted.append(
                    ExtractedParameter(f"{well_type}_wells", int(match.group(1)), "count", match.group(0))
                )
            
            return extracted
        
        @tool
        def convert_to_si_units(value: float, from_unit: str, parameter_type: str) -> Tuple[float, str]:
            """Convert value to SI units used by DARTS.
            
            Args:
                value: Numerical value
                from_unit: Original unit
                parameter_type: Type of parameter (length, pressure, etc.)
                
            Returns:
                Tuple of (converted_value, si_unit)
            """
            from_unit = from_unit.lower()
            
            # Length conversions
            if from_unit in ["ft", "feet"]:
                return value * UNIT_CONVERSIONS["ft_to_m"], "m"
            
            # Pressure conversions
            if from_unit == "psi":
                return value * UNIT_CONVERSIONS["psi_to_pa"], "Pa"
            elif from_unit == "bar":
                return value * UNIT_CONVERSIONS["bar_to_pa"], "Pa"
            elif from_unit == "atm":
                return value * UNIT_CONVERSIONS["atm_to_pa"], "Pa"
            
            # Temperature conversions
            if from_unit == "f":
                return UNIT_CONVERSIONS["f_to_c"](value), "C"
            elif from_unit == "k":
                return value - 273.15, "C"
            
            # Volume rate conversions
            if from_unit in ["bbl/day", "stb/d"]:
                return value * UNIT_CONVERSIONS["stb_d_to_m3_s"], "m3/s"
            elif from_unit == "m3/day":
                return value / 86400, "m3/s"
            
            # Permeability conversions - keep in mD
            if from_unit in ["md", "millidarcy"]:
                return value, "mD"
            elif from_unit == "darcy":
                return value * 1000, "mD"  # Convert darcy to millidarcy
            
            # Time conversions
            if from_unit == "years":
                return value * 365, "days"
            elif from_unit == "months":
                return value * 30, "days"
            
            # No conversion needed
            return value, from_unit
        
        @tool
        def get_template_defaults(template_name: str) -> Dict[str, Any]:
            """Get default parameters for a template.
            
            Args:
                template_name: Name of the template
                
            Returns:
                Dictionary of default parameters
            """
            if template_name not in TEMPLATES:
                return {}
            
            template = TEMPLATES[template_name]
            
            # Common defaults
            defaults = {
                "grid": {
                    "nx": 50,
                    "ny": 50,
                    "nz": 10,
                    "dx": 10.0,  # m
                    "dy": 10.0,  # m
                    "dz": 2.0    # m
                },
                "rock": {
                    "porosity": 0.2,
                    "permeability": 100,  # mD
                    "compressibility": 1e-9   # 1/Pa
                },
                "simulation": {
                    "time_step": 1.0,      # days
                    "total_time": 365.0,   # days
                    "output_freq": 30.0    # days
                }
            }
            
            # Template-specific defaults
            if template.physics_type == "compositional":
                defaults["initial_composition"] = {"CO2": 0.1, "C1": 0.9}
                defaults["injection_composition"] = {"CO2": 1.0}
            elif template.physics_type == "geothermal":
                defaults["temperature"] = {
                    "initial": 150.0,      # C
                    "injection": 50.0      # C
                }
                defaults["thermal"] = {
                    "rock_heat_capacity": 1000.0,     # J/kg/K
                    "rock_thermal_cond": 2.0,         # W/m/K
                    "fluid_heat_capacity": 4200.0     # J/kg/K
                }
            elif template.physics_type == "dead_oil":
                defaults["initial_saturation"] = {
                    "water": 0.2,
                    "oil": 0.8
                }
            
            # Override with template-specific parameters
            if template.parameters:
                for category, params in template.parameters.items():
                    if category not in defaults:
                        defaults[category] = {}
                    defaults[category].update(params)
            
            return defaults
        
        @tool
        def validate_parameter_ranges(parameters: Dict[str, Any]) -> Dict[str, List[str]]:
            """Validate parameters are within reasonable ranges.
            
            Args:
                parameters: Dictionary of parameters to validate
                
            Returns:
                Dictionary of parameter names to validation warnings
            """
            warnings = {}
            
            # Grid validation
            if "nx" in parameters:
                if parameters["nx"] < 1 or parameters["nx"] > 1000:
                    warnings.setdefault("nx", []).append(f"Grid dimension nx={parameters['nx']} may be too {'small' if parameters['nx'] < 1 else 'large'}")
            
            # Porosity validation
            if "porosity" in parameters:
                if parameters["porosity"] < 0 or parameters["porosity"] > 1:
                    warnings.setdefault("porosity", []).append(f"Porosity must be between 0 and 1, got {parameters['porosity']}")
                elif parameters["porosity"] < 0.01 or parameters["porosity"] > 0.5:
                    warnings.setdefault("porosity", []).append(f"Porosity {parameters['porosity']} is unusual")
            
            # Permeability validation (in mD)
            if "permeability" in parameters:
                perm_md = parameters["permeability"]
                if perm_md < 0.01 or perm_md > 10000:
                    warnings.setdefault("permeability", []).append(f"Permeability {perm_md} mD is unusual")
            
            # Temperature validation (in C)
            if "temperature" in parameters:
                if parameters["temperature"] < -50 or parameters["temperature"] > 400:
                    warnings.setdefault("temperature", []).append(f"Temperature {parameters['temperature']}°C may be unrealistic")
            
            # Pressure validation (in Pa)
            if "pressure" in parameters:
                pressure_bar = parameters["pressure"] / UNIT_CONVERSIONS["bar_to_pa"]
                if pressure_bar < 1 or pressure_bar > 1000:
                    warnings.setdefault("pressure", []).append(f"Pressure {pressure_bar} bar may be unrealistic")
            
            return warnings
        
        return [
            extract_numerical_values,
            convert_to_si_units,
            get_template_defaults,
            validate_parameter_ranges
        ]
    
    def get_prompt_template(self) -> PromptTemplate:
        """Get the prompt template for parameter extraction."""
        template = """You are an expert DARTS parameter extractor.
        
Extract numerical parameters from the user request and template requirements.

User Request: {user_request}

Intent Classification: {intent_classification}

Selected Template: {template_info}

Extracted Values: {extracted_values}

Template Defaults: {template_defaults}

Instructions:
1. Use extracted values when available
2. Convert all units to DARTS standard (SI units)
3. Apply template defaults for missing parameters
4. Organize parameters by category:
   - grid_parameters: nx, ny, nz, dx, dy, dz
   - rock_properties: porosity, permeability, compressibility
   - fluid_properties: density, viscosity, etc.
   - well_parameters: locations, rates, controls
   - simulation_parameters: time_step, total_time, output
   - custom_parameters: template-specific parameters

5. Document all unit conversions performed
6. Explain any assumptions or defaults applied

{format_instructions}
"""
        
        return PromptTemplate(
            template=template,
            input_variables=[
                "user_request",
                "intent_classification",
                "template_info",
                "extracted_values",
                "template_defaults"
            ],
            partial_variables={
                "format_instructions": self.output_parser.get_format_instructions()
            }
        )
    
    async def process(self, input_data: AgentInput) -> AgentOutput:
        """Process the template selection and extract parameters.
        
        Args:
            input_data: Input containing template selection
            
        Returns:
            AgentOutput with extracted parameters
        """
        if not self.validate_input(input_data):
            return AgentOutput(
                result={"error": "Invalid input"},
                reasoning="Empty or invalid input provided",
                confidence=0.0
            )
        
        # Get context from previous agents
        template_info = input_data.context.get("template_selection", {})
        intent_info = input_data.context.get("intent_classification", {})
        template_name = template_info.get("template_name", "2ph_do")
        
        # Use tools to extract parameters
        tools = self.get_tools()
        extract_tool = tools[0]
        convert_tool = tools[1]
        defaults_tool = tools[2]
        validate_tool = tools[3]
        
        # Extract numerical values
        extracted_params = extract_tool.invoke(input_data.message)
        
        # Get template defaults
        template_defaults = defaults_tool.invoke({"template_name": template_name})
        
        # Convert extracted parameters to SI units
        converted_params = {}
        units_used = {}
        
        for param in extracted_params:
            if param.unit:
                converted_value, si_unit = convert_tool.invoke({
                    "value": param.value,
                    "from_unit": param.unit,
                    "parameter_type": param.name
                })
                converted_params[param.name] = converted_value
                units_used[param.name] = si_unit
            else:
                converted_params[param.name] = param.value
                units_used[param.name] = "dimensionless"
        
        # Prepare prompt
        prompt = self.get_prompt_template()
        formatted_prompt = prompt.format(
            user_request=input_data.message,
            intent_classification=json.dumps(intent_info, indent=2),
            template_info=json.dumps(template_info, indent=2),
            extracted_values=json.dumps(
                [{"name": p.name, "value": p.value, "unit": p.unit, "text": p.original_text} 
                 for p in extracted_params],
                indent=2
            ),
            template_defaults=json.dumps(template_defaults, indent=2)
        )
        
        # Get LLM response
        messages = [
            SystemMessage(content="You are an expert DARTS parameter extractor."),
            HumanMessage(content=formatted_prompt)
        ]
        
        response = await self.llm.ainvoke(messages)
        
        # Parse response
        try:
            extraction = self.output_parser.parse(response.content)
            
            # Validate parameters
            all_params = {}
            for category in ["grid_parameters", "rock_properties", "fluid_properties",
                           "well_parameters", "simulation_parameters", "custom_parameters"]:
                all_params.update(extraction.model_dump().get(category, {}))
            
            warnings = validate_tool.invoke({"parameters": all_params})
            
            # Add warnings to reasoning if any
            if warnings:
                extraction.reasoning += f"\n\nValidation warnings: {json.dumps(warnings, indent=2)}"
            
            # Update context
            updated_context = input_data.context.copy()
            updated_context["parameter_extraction"] = extraction.model_dump()
            
            return AgentOutput(
                result=extraction.model_dump(),
                reasoning=extraction.reasoning,
                confidence=0.9 if not warnings else 0.7,
                next_agent="CodeGenerator"
            )
        except Exception as e:
            # Fallback to basic parameter set
            fallback_params = {
                "grid_parameters": converted_params if converted_params else {"nx": 50, "ny": 50, "nz": 10},
                "rock_properties": {"porosity": 0.2, "permeability": 100},  # permeability in mD
                "fluid_properties": {},
                "well_parameters": {},
                "simulation_parameters": {"total_time": 365, "time_step": 1},
                "custom_parameters": {},
                "units_used": units_used,
                "reasoning": f"Fallback parameter extraction: {str(e)}"
            }
            
            updated_context = input_data.context.copy()
            updated_context["parameter_extraction"] = fallback_params
            
            return AgentOutput(
                result=fallback_params,
                reasoning=fallback_params["reasoning"],
                confidence=0.5,
                next_agent="CodeGenerator"
            )