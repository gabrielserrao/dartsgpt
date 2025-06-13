#!/usr/bin/env python
"""Test the full DARTSGPT agent pipeline with minimal LangChain usage."""

import sys
from pathlib import Path
import asyncio
from typing import Dict, Any

# Add dartsgpt to path
sys.path.insert(0, str(Path(__file__).parent))

from dartsgpt.config.settings import Settings
from dartsgpt.knowledge.templates.template_database import TEMPLATES
from dartsgpt.knowledge.physics_features import PHYSICS_KEYWORDS, UNIT_CONVERSIONS
from dartsgpt.knowledge.example_prompts import EXAMPLE_PROMPTS

# Import agents
from dartsgpt.agents.intent_classifier import IntentClassifierAgent
from dartsgpt.agents.template_selector import TemplateSelectorAgent
from dartsgpt.agents.parameter_extractor import ParameterExtractorAgent
from dartsgpt.agents.code_generator import CodeGeneratorAgent
from dartsgpt.agents.validator import ValidatorAgent


class SimpleOrchestrator:
    """Simple orchestrator for testing without full LangGraph."""
    
    def __init__(self):
        self.settings = Settings()
        self.intent_classifier = IntentClassifierAgent()
        self.template_selector = TemplateSelectorAgent()
        self.parameter_extractor = ParameterExtractorAgent()
        self.code_generator = CodeGeneratorAgent()
        self.validator = ValidatorAgent()
    
    async def run(self, prompt: str) -> Dict[str, Any]:
        """Run the agent pipeline."""
        print(f"\n{'='*60}")
        print(f"PROMPT: {prompt}")
        print(f"{'='*60}")
        
        state = {"prompt": prompt}
        
        # Step 1: Intent Classification
        print("\n1. INTENT CLASSIFICATION")
        try:
            state = await self.intent_classifier.classify(state)
            print(f"   Physics type: {state['intent']['physics_type']}")
            print(f"   Features: {state['intent']['features']}")
            print(f"   Complexity: {state['intent']['complexity']}")
        except Exception as e:
            print(f"   ❌ Error: {e}")
            # Fallback to simple classification
            state['intent'] = self._simple_intent_classification(prompt)
            print(f"   Using fallback classification")
        
        # Step 2: Template Selection
        print("\n2. TEMPLATE SELECTION")
        try:
            state = await self.template_selector.select(state)
            print(f"   Selected: {state['selected_template']}")
            template = TEMPLATES[state['selected_template']]
            print(f"   Description: {template.description}")
        except Exception as e:
            print(f"   ❌ Error: {e}")
            # Fallback to simple selection
            state['selected_template'] = self._simple_template_selection(state['intent'])
            print(f"   Using fallback selection: {state['selected_template']}")
        
        # Step 3: Parameter Extraction
        print("\n3. PARAMETER EXTRACTION")
        try:
            state = await self.parameter_extractor.extract(state)
            print(f"   Parameters extracted: {len(state.get('parameters', {}))} groups")
            if 'parameters' in state:
                for key, value in state['parameters'].items():
                    if value:
                        print(f"   - {key}: {value}")
        except Exception as e:
            print(f"   ❌ Error: {e}")
            # Fallback to simple extraction
            state['parameters'] = self._simple_parameter_extraction(prompt)
            print(f"   Using fallback extraction")
        
        # Step 4: Code Generation
        print("\n4. CODE GENERATION")
        try:
            state = await self.code_generator.generate(state)
            print(f"   Generated code: ✓")
            print(f"   Model code: {len(state.get('model_code', ''))} characters")
            print(f"   Main code: {len(state.get('main_code', ''))} characters")
        except Exception as e:
            print(f"   ❌ Error: {e}")
            # Fallback to simple generation
            self._simple_code_generation(state)
            print(f"   Using fallback generation")
        
        # Step 5: Validation (optional)
        print("\n5. VALIDATION")
        try:
            state = await self.validator.validate(state)
            if state.get('validation_passed'):
                print(f"   ✅ Validation passed")
            else:
                print(f"   ⚠️  Validation issues found")
                for issue in state.get('validation_issues', []):
                    print(f"      - {issue}")
        except Exception as e:
            print(f"   ❌ Validation error: {e}")
            print(f"   Skipping validation")
        
        return state
    
    def _simple_intent_classification(self, prompt: str) -> dict:
        """Fallback intent classification."""
        prompt_lower = prompt.lower()
        
        physics_type = "dead_oil"
        for ptype, keywords in PHYSICS_KEYWORDS.items():
            if any(kw in prompt_lower for kw in keywords):
                physics_type = ptype
                break
        
        features = []
        if "injection" in prompt_lower:
            features.append("injection")
        if "production" in prompt_lower:
            features.append("production")
        
        return {
            "physics_type": physics_type,
            "features": features,
            "complexity": "moderate"
        }
    
    def _simple_template_selection(self, intent: dict) -> str:
        """Fallback template selection."""
        physics_type = intent.get('physics_type', 'dead_oil')
        candidates = [
            (name, template) for name, template in TEMPLATES.items()
            if template.physics_type == physics_type
        ]
        
        if not candidates:
            return "2ph_do"
        
        # Prefer golden templates
        golden = [(n, t) for n, t in candidates if t.is_golden]
        if golden:
            return golden[0][0]
        
        return candidates[0][0]
    
    def _simple_parameter_extraction(self, prompt: str) -> dict:
        """Fallback parameter extraction."""
        import re
        
        params = {
            "grid_parameters": {},
            "rock_properties": {},
            "well_parameters": {}
        }
        
        # Grid dimensions
        grid_pattern = r'(\d+)\s*[xX×]\s*(\d+)\s*[xX×]\s*(\d+)'
        grid_match = re.search(grid_pattern, prompt)
        if grid_match:
            params["grid_parameters"]["nx"] = int(grid_match.group(1))
            params["grid_parameters"]["ny"] = int(grid_match.group(2))
            params["grid_parameters"]["nz"] = int(grid_match.group(3))
        
        return params
    
    def _simple_code_generation(self, state: dict) -> None:
        """Fallback code generation."""
        template_name = state.get('selected_template', '2ph_do')
        template = TEMPLATES[template_name]
        
        model_code = f'''"""Generated DARTS model using template: {template_name}"""

import numpy as np
from darts.models.cicd_model import CICDModel


class Model(CICDModel):
    """Generated {template.description}"""
    
    def __init__(self):
        super().__init__()
        self.physics_type = "{template.physics_type}"
'''
        
        main_code = '''"""Main script to run the DARTS model."""

from model import Model


def main():
    model = Model()
    model.init()
    model.run()
    model.save_results()


if __name__ == "__main__":
    main()
'''
        
        state['model_code'] = model_code
        state['main_code'] = main_code


async def test_agent_pipeline():
    """Test the full agent pipeline."""
    orchestrator = SimpleOrchestrator()
    
    # Test prompts
    test_prompts = [
        "Create a CO2 injection model for carbon storage in a 100x100x20 reservoir",
        "Build a waterflood simulation with five-spot pattern",
        "Model a geothermal reservoir at 150°C"
    ]
    
    for prompt in test_prompts:
        try:
            result = await orchestrator.run(prompt)
            
            # Save output
            output_dir = Path("agent_output")
            output_dir.mkdir(exist_ok=True)
            
            if 'model_code' in result:
                (output_dir / "model.py").write_text(result['model_code'])
            if 'main_code' in result:
                (output_dir / "main.py").write_text(result['main_code'])
            
            print(f"\n✅ Files saved to: {output_dir}/")
            
        except Exception as e:
            print(f"\n❌ Pipeline error: {e}")
            import traceback
            traceback.print_exc()


def main():
    """Run the test."""
    print("DARTSGPT Full Agent Pipeline Test")
    print("=" * 60)
    
    # Check if API key is set
    settings = Settings()
    if not settings.openai_api_key:
        print("\n⚠️  WARNING: No OpenAI API key found!")
        print("The agents will use fallback methods instead of LLM calls.")
        print("To use LLMs, set OPENAI_API_KEY in your .env file.")
    else:
        print(f"\n✓ OpenAI API key found")
        print(f"✓ Using model: {settings.openai_model}")
    
    print("\nStarting agent pipeline tests...\n")
    
    # Run async tests
    import nest_asyncio
    nest_asyncio.apply()
    
    asyncio.run(test_agent_pipeline())
    
    print("\n" + "=" * 60)
    print("Test complete!")


if __name__ == "__main__":
    main()