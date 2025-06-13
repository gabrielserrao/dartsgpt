"""Template Selector Agent for DARTSGPT.

This agent selects the most appropriate DARTS template based on:
- Physics type identified by Intent Classifier
- Required features
- Model complexity
- Historical success rates
"""

from typing import Dict, Any, List, Optional, Tuple
import json

from langchain.prompts import PromptTemplate
from langchain.tools import BaseTool, tool
from langchain.schema import HumanMessage, SystemMessage
from langchain_core.output_parsers import JsonOutputParser
from pydantic import BaseModel, Field

from ..agents.base_agent import BaseAgent, AgentInput, AgentOutput
from ..knowledge.templates.template_database import TEMPLATES, TemplateInfo
from ..knowledge.rag_system import RAGSystem


class TemplateSelection(BaseModel):
    """Template selection result."""
    template_name: str = Field(description="Selected template name")
    template_path: str = Field(description="Path to template files")
    physics_match: float = Field(description="Physics compatibility score 0-1")
    feature_match: float = Field(description="Feature compatibility score 0-1")
    overall_score: float = Field(description="Overall template score 0-1")
    reasoning: str = Field(description="Explanation of the selection")
    alternatives: List[Dict[str, Any]] = Field(description="Alternative templates considered")


class TemplateSelectorAgent(BaseAgent):
    """Agent that selects the best template for user requirements."""
    
    def __init__(self, llm=None, rag_system: Optional[RAGSystem] = None):
        """Initialize the Template Selector Agent.
        
        Args:
            llm: Language model to use
            rag_system: RAG system for similarity search
        """
        super().__init__(llm)
        self.rag_system = rag_system
        self.output_parser = JsonOutputParser(pydantic_object=TemplateSelection)
        self.templates = TEMPLATES
    
    def get_tools(self) -> List[BaseTool]:
        """Get tools for template selection."""
        
        @tool
        def score_physics_compatibility(template_physics: str, required_physics: str) -> float:
            """Score physics compatibility between template and requirements.
            
            Args:
                template_physics: Physics type of the template
                required_physics: Required physics type
                
            Returns:
                Compatibility score 0-1
            """
            # Exact match
            if template_physics == required_physics:
                return 1.0
            
            # Compatible physics
            compatibility_map = {
                "compositional": ["dead_oil"],  # Can simplify to dead oil
                "dead_oil": [],
                "black_oil": ["dead_oil"],  # Can simplify
                "geothermal": ["dead_oil"],  # With thermal
                "poroelastic": [],
                "specialized": []
            }
            
            if required_physics in compatibility_map.get(template_physics, []):
                return 0.7
            
            return 0.0
        
        @tool
        def score_feature_compatibility(template_features: List[str], 
                                      required_features: List[str]) -> float:
            """Score feature compatibility between template and requirements.
            
            Args:
                template_features: Features supported by template
                required_features: Required features
                
            Returns:
                Compatibility score 0-1
            """
            if not required_features:
                return 1.0
            
            if not template_features:
                return 0.0
            
            matched = sum(1 for feat in required_features if feat in template_features)
            return matched / len(required_features)
        
        @tool
        def score_complexity_match(template_complexity: str, 
                                 required_complexity: str) -> float:
            """Score complexity match between template and requirements.
            
            Args:
                template_complexity: Template complexity level
                required_complexity: Required complexity level
                
            Returns:
                Match score 0-1
            """
            levels = ["simple", "moderate", "complex"]
            
            if template_complexity == required_complexity:
                return 1.0
            
            try:
                template_idx = levels.index(template_complexity)
                required_idx = levels.index(required_complexity)
                distance = abs(template_idx - required_idx)
                return max(0, 1 - distance * 0.3)
            except ValueError:
                return 0.5
        
        @tool
        def get_template_metadata(template_name: str) -> Dict[str, Any]:
            """Get metadata for a specific template.
            
            Args:
                template_name: Name of the template
                
            Returns:
                Template metadata dictionary
            """
            if template_name not in TEMPLATES:
                return {"error": f"Template {template_name} not found"}
            
            template = TEMPLATES[template_name]
            return {
                "name": template.name,
                "path": template.path,
                "physics_type": template.physics_type,
                "description": template.description,
                "grid_type": template.grid_type,
                "features": template.features,
                "complexity": template.complexity,
                "use_cases": template.use_cases,
                "parameters": template.parameters
            }
        
        @tool
        def find_similar_examples(physics_type: str, features: List[str]) -> List[Dict[str, Any]]:
            """Find similar examples using RAG search.
            
            Args:
                physics_type: Physics type to search for
                features: Features to match
                
            Returns:
                List of similar examples with metadata
            """
            if not self.rag_system:
                return []
            
            # Construct search query
            query = f"physics:{physics_type}"
            if features:
                query += " " + " ".join(features)
            
            # Search for similar examples
            results = self.rag_system.search(query, k=5)
            
            return [
                {
                    "content": r.page_content,
                    "metadata": r.metadata,
                    "score": getattr(r, 'score', 0.0)
                }
                for r in results
            ]
        
        return [
            score_physics_compatibility,
            score_feature_compatibility,
            score_complexity_match,
            get_template_metadata,
            find_similar_examples
        ]
    
    def get_prompt_template(self) -> PromptTemplate:
        """Get the prompt template for template selection."""
        template = """You are an expert DARTS template selector.
        
Select the best template based on the classified intent and requirements.

Intent Classification:
{intent_classification}

Available Templates:
{available_templates}

Template Scores:
{template_scores}

Similar Examples:
{similar_examples}

Selection Criteria:
1. Physics compatibility is most important (weight: 0.5)
2. Feature support is critical (weight: 0.3)
3. Complexity match matters (weight: 0.2)
4. Prefer templates marked as "golden" (⭐) when scores are close

Calculate overall scores and select the best template. Also identify 2-3 alternatives.

{format_instructions}
"""
        
        return PromptTemplate(
            template=template,
            input_variables=[
                "intent_classification",
                "available_templates",
                "template_scores",
                "similar_examples"
            ],
            partial_variables={
                "format_instructions": self.output_parser.get_format_instructions()
            }
        )
    
    def _score_templates(self, physics_type: str, features: List[str], 
                        complexity: str) -> List[Tuple[str, Dict[str, float]]]:
        """Score all templates based on requirements.
        
        Args:
            physics_type: Required physics type
            features: Required features
            complexity: Required complexity
            
        Returns:
            List of (template_name, scores) tuples sorted by overall score
        """
        tools = self.get_tools()
        physics_scorer = tools[0]
        feature_scorer = tools[1]
        complexity_scorer = tools[2]
        
        template_scores = []
        
        for name, template in self.templates.items():
            # Calculate individual scores
            physics_score = physics_scorer.invoke({
                "template_physics": template.physics_type,
                "required_physics": physics_type
            })
            
            feature_score = feature_scorer.invoke({
                "template_features": template.features,
                "required_features": features
            })
            
            complexity_score = complexity_scorer.invoke({
                "template_complexity": template.complexity,
                "required_complexity": complexity
            })
            
            # Calculate weighted overall score
            overall_score = (
                physics_score * 0.5 +
                feature_score * 0.3 +
                complexity_score * 0.2
            )
            
            # Bonus for golden templates
            if template.is_golden:
                overall_score *= 1.1
            
            scores = {
                "physics": physics_score,
                "features": feature_score,
                "complexity": complexity_score,
                "overall": min(overall_score, 1.0)
            }
            
            template_scores.append((name, scores))
        
        # Sort by overall score
        template_scores.sort(key=lambda x: x[1]["overall"], reverse=True)
        
        return template_scores
    
    async def process(self, input_data: AgentInput) -> AgentOutput:
        """Process the intent classification and select a template.
        
        Args:
            input_data: Input containing intent classification
            
        Returns:
            AgentOutput with template selection
        """
        if not self.validate_input(input_data):
            return AgentOutput(
                result={"error": "Invalid input"},
                reasoning="Empty or invalid input provided",
                confidence=0.0
            )
        
        # Extract intent classification from context
        intent_data = input_data.context.get("intent_classification", {})
        physics_type = intent_data.get("physics_type", "dead_oil")
        features = intent_data.get("features", [])
        complexity = intent_data.get("complexity", "moderate")
        
        # Score all templates
        template_scores = self._score_templates(physics_type, features, complexity)
        
        # Get template metadata
        tools = self.get_tools()
        metadata_tool = tools[3]
        similar_tool = tools[4]
        
        # Get metadata for top templates
        available_templates = {}
        for name, _ in template_scores[:5]:
            available_templates[name] = metadata_tool.invoke({"template_name": name})
        
        # Find similar examples
        similar_examples = similar_tool.invoke({
            "physics_type": physics_type,
            "features": features
        })
        
        # Format scores for prompt
        scores_formatted = json.dumps(
            {name: scores for name, scores in template_scores[:5]},
            indent=2
        )
        
        # Prepare prompt
        prompt = self.get_prompt_template()
        formatted_prompt = prompt.format(
            intent_classification=json.dumps(intent_data, indent=2),
            available_templates=json.dumps(available_templates, indent=2),
            template_scores=scores_formatted,
            similar_examples=json.dumps(similar_examples, indent=2)
        )
        
        # Get LLM response
        messages = [
            SystemMessage(content="You are an expert DARTS template selector."),
            HumanMessage(content=formatted_prompt)
        ]
        
        response = await self.llm.ainvoke(messages)
        
        # Parse response
        try:
            selection = self.output_parser.parse(response.content)
            
            # Add template info to context
            updated_context = input_data.context.copy()
            updated_context["template_selection"] = selection.model_dump()
            
            return AgentOutput(
                result=selection.model_dump(),
                reasoning=selection.reasoning,
                confidence=selection.overall_score,
                next_agent="ParameterExtractor"
            )
        except Exception as e:
            # Fallback to best scoring template
            best_template = template_scores[0][0]
            best_scores = template_scores[0][1]
            
            fallback_selection = {
                "template_name": best_template,
                "template_path": self.templates[best_template].path,
                "physics_match": best_scores["physics"],
                "feature_match": best_scores["features"],
                "overall_score": best_scores["overall"],
                "reasoning": f"Fallback selection based on scoring: {str(e)}",
                "alternatives": [
                    {"name": name, "score": scores["overall"]}
                    for name, scores in template_scores[1:4]
                ]
            }
            
            updated_context = input_data.context.copy()
            updated_context["template_selection"] = fallback_selection
            
            return AgentOutput(
                result=fallback_selection,
                reasoning=fallback_selection["reasoning"],
                confidence=fallback_selection["overall_score"],
                next_agent="ParameterExtractor"
            )