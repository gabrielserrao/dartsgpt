"""Intent Classifier Agent for DARTSGPT.

This agent analyzes natural language prompts to determine:
- Physics type (compositional, dead oil, black oil, geothermal, etc.)
- Key features required
- Model complexity
- Special requirements
"""

from typing import Dict, Any, List, Optional
import re
import json

from langchain.prompts import PromptTemplate
from langchain.tools import BaseTool, tool
from langchain.schema import HumanMessage, SystemMessage
from langchain_core.output_parsers import JsonOutputParser
from pydantic import BaseModel, Field

from ..agents.base_agent import BaseAgent, AgentInput, AgentOutput
from ..knowledge.physics_features import PHYSICS_KEYWORDS, FEATURE_KEYWORDS


class IntentClassification(BaseModel):
    """Classification result from intent analysis."""
    physics_type: str = Field(description="Identified physics type")
    features: List[str] = Field(description="Required features")
    complexity: str = Field(description="Model complexity: simple, moderate, or complex")
    grid_info: Optional[Dict[str, Any]] = Field(description="Grid specifications if provided")
    well_info: Optional[Dict[str, Any]] = Field(description="Well configuration if provided")
    parameters: Dict[str, Any] = Field(description="Extracted numerical parameters")
    confidence: float = Field(description="Confidence score 0-1")
    reasoning: str = Field(description="Explanation of the classification")


class IntentClassifierAgent(BaseAgent):
    """Agent that classifies user intent and extracts requirements."""
    
    def __init__(self, llm=None):
        """Initialize the Intent Classifier Agent."""
        super().__init__(llm)
        self.output_parser = JsonOutputParser(pydantic_object=IntentClassification)
    
    def get_tools(self) -> List[BaseTool]:
        """Get tools for intent classification."""
        
        @tool
        def analyze_physics_keywords(text: str) -> Dict[str, List[str]]:
            """Analyze text for physics-related keywords.
            
            Args:
                text: Input text to analyze
                
            Returns:
                Dictionary mapping physics types to found keywords
            """
            text_lower = text.lower()
            found_physics = {}
            
            for physics_type, keywords in PHYSICS_KEYWORDS.items():
                found = [kw for kw in keywords if kw in text_lower]
                if found:
                    found_physics[physics_type] = found
            
            return found_physics
        
        @tool
        def analyze_feature_keywords(text: str) -> List[str]:
            """Analyze text for feature-related keywords.
            
            Args:
                text: Input text to analyze
                
            Returns:
                List of identified features
            """
            text_lower = text.lower()
            found_features = []
            
            for feature, keywords in FEATURE_KEYWORDS.items():
                if any(kw in text_lower for kw in keywords):
                    found_features.append(feature)
            
            return found_features
        
        @tool
        def extract_grid_specifications(text: str) -> Optional[Dict[str, Any]]:
            """Extract grid specifications from text.
            
            Args:
                text: Input text to analyze
                
            Returns:
                Dictionary with grid specifications or None
            """
            grid_info = {}
            
            # Look for grid dimensions (e.g., "50x50x10", "100 x 100 x 20")
            grid_pattern = r'(\d+)\s*[xX×]\s*(\d+)\s*[xX×]\s*(\d+)'
            match = re.search(grid_pattern, text)
            if match:
                grid_info['nx'] = int(match.group(1))
                grid_info['ny'] = int(match.group(2))
                grid_info['nz'] = int(match.group(3))
            
            # Look for cell dimensions
            cell_pattern = r'(\d+\.?\d*)\s*m\s*[xX×]\s*(\d+\.?\d*)\s*m\s*[xX×]\s*(\d+\.?\d*)\s*m'
            match = re.search(cell_pattern, text)
            if match:
                grid_info['dx'] = float(match.group(1))
                grid_info['dy'] = float(match.group(2))
                grid_info['dz'] = float(match.group(3))
            
            # Look for 2D grids
            grid_2d_pattern = r'(\d+)\s*[xX×]\s*(\d+)(?!\s*[xX×])'
            if 'nx' not in grid_info:
                match = re.search(grid_2d_pattern, text)
                if match:
                    grid_info['nx'] = int(match.group(1))
                    grid_info['ny'] = int(match.group(2))
                    grid_info['nz'] = 1
            
            return grid_info if grid_info else None
        
        @tool
        def extract_well_configuration(text: str) -> Optional[Dict[str, Any]]:
            """Extract well configuration from text.
            
            Args:
                text: Input text to analyze
                
            Returns:
                Dictionary with well configuration or None
            """
            well_info = {}
            text_lower = text.lower()
            
            # Well types
            if 'injection' in text_lower:
                well_info['type'] = 'injection'
                if 'co2' in text_lower:
                    well_info['fluid'] = 'CO2'
                elif 'water' in text_lower:
                    well_info['fluid'] = 'water'
                elif 'gas' in text_lower:
                    well_info['fluid'] = 'gas'
            elif 'production' in text_lower:
                well_info['type'] = 'production'
            
            # Well patterns
            if 'five-spot' in text_lower or '5-spot' in text_lower:
                well_info['pattern'] = 'five-spot'
            elif 'line drive' in text_lower:
                well_info['pattern'] = 'line-drive'
            elif 'doublet' in text_lower:
                well_info['pattern'] = 'doublet'
            
            # Number of wells
            well_count_pattern = r'(\d+)\s+(?:injection|production|well)'
            matches = re.findall(well_count_pattern, text_lower)
            if matches:
                well_info['count'] = int(matches[0])
            
            return well_info if well_info else None
        
        return [
            analyze_physics_keywords,
            analyze_feature_keywords,
            extract_grid_specifications,
            extract_well_configuration
        ]
    
    def get_prompt_template(self) -> PromptTemplate:
        """Get the prompt template for intent classification."""
        template = """You are an expert DARTS reservoir simulation intent classifier.
        
Analyze the user's request and classify their intent for creating a DARTS model.

User Request: {user_request}

Context:
{context}

Physics Keywords Found: {physics_keywords}
Feature Keywords Found: {feature_keywords}
Grid Specifications: {grid_specs}
Well Configuration: {well_config}

Classify the intent with the following structure:
1. Determine the physics type (must be one of: compositional, dead_oil, black_oil, geothermal, poroelastic, specialized)
2. Identify required features
3. Assess complexity (simple, moderate, or complex)
4. Extract any numerical parameters

Consider:
- CO2 injection/storage → compositional physics
- Water flooding → dead_oil physics
- Gas injection with oil → black_oil or compositional
- Temperature/heat → geothermal or thermal variants
- Stress/mechanics → poroelastic
- Chemical/foam → specialized

{format_instructions}
"""
        
        return PromptTemplate(
            template=template,
            input_variables=[
                "user_request",
                "context", 
                "physics_keywords",
                "feature_keywords",
                "grid_specs",
                "well_config"
            ],
            partial_variables={
                "format_instructions": self.output_parser.get_format_instructions()
            }
        )
    
    async def process(self, input_data: AgentInput) -> AgentOutput:
        """Process the user request and classify intent.
        
        Args:
            input_data: Input containing user request
            
        Returns:
            AgentOutput with classification results
        """
        if not self.validate_input(input_data):
            return AgentOutput(
                result={"error": "Invalid input"},
                reasoning="Empty or invalid input provided",
                confidence=0.0
            )
        
        # Use tools to analyze the request
        tools = self.get_tools()
        physics_keywords = tools[0].invoke(input_data.message)
        feature_keywords = tools[1].invoke(input_data.message)
        grid_specs = tools[2].invoke(input_data.message)
        well_config = tools[3].invoke(input_data.message)
        
        # Prepare prompt
        prompt = self.get_prompt_template()
        formatted_prompt = prompt.format(
            user_request=input_data.message,
            context=self.format_context(input_data.context),
            physics_keywords=json.dumps(physics_keywords, indent=2),
            feature_keywords=json.dumps(feature_keywords, indent=2),
            grid_specs=json.dumps(grid_specs, indent=2),
            well_config=json.dumps(well_config, indent=2)
        )
        
        # Get LLM response
        messages = [
            SystemMessage(content="You are an expert DARTS intent classifier."),
            HumanMessage(content=formatted_prompt)
        ]
        
        response = await self.llm.ainvoke(messages)
        
        # Parse response
        try:
            classification = self.output_parser.parse(response.content)
            
            # Add grid and well info if found
            if grid_specs:
                classification.grid_info = grid_specs
            if well_config:
                classification.well_info = well_config
            
            return AgentOutput(
                result=classification.model_dump(),
                reasoning=classification.reasoning,
                confidence=classification.confidence,
                next_agent="TemplateSelector"
            )
        except Exception as e:
            return AgentOutput(
                result={"error": str(e)},
                reasoning=f"Failed to parse classification: {str(e)}",
                confidence=0.0
            )