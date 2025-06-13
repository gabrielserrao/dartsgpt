#!/usr/bin/env python
"""Test script for DARTSGPT system components."""

import asyncio
import sys
from pathlib import Path

# Add dartsgpt to path
sys.path.insert(0, str(Path(__file__).parent))

from dartsgpt.config import get_settings
from dartsgpt.agents.intent_classifier import IntentClassifierAgent
from dartsgpt.agents.template_selector import TemplateSelectorAgent
from dartsgpt.agents.parameter_extractor import ParameterExtractorAgent
from dartsgpt.agents.code_generator import CodeGeneratorAgent
from dartsgpt.agents.validator import ValidatorAgent
from dartsgpt.agents.base_agent import AgentInput
from dartsgpt.knowledge.templates.template_database import TEMPLATES


async def test_basic_flow():
    """Test basic agent flow with a simple prompt."""
    print("=== Testing Basic Agent Flow ===\n")
    
    prompt = "Create a CO2 injection model for carbon storage"
    print(f"Test prompt: {prompt}\n")
    
    # Test settings
    try:
        settings = get_settings()
        print(f"✅ Settings loaded - Model: {settings.openai_model}")
    except Exception as e:
        print(f"❌ Settings error: {e}")
        return
    
    # Test agents
    context = {}
    
    # 1. Intent Classification
    print("\n1. Testing Intent Classifier...")
    try:
        intent_agent = IntentClassifierAgent()
        intent_input = AgentInput(message=prompt, context=context)
        intent_result = await intent_agent.process(intent_input)
        context['intent_classification'] = intent_result.result
        
        print(f"   ✅ Physics type: {intent_result.result.get('physics_type')}")
        print(f"   ✅ Features: {intent_result.result.get('features')}")
        print(f"   ✅ Confidence: {intent_result.confidence:.2f}")
    except Exception as e:
        print(f"   ❌ Intent classification failed: {e}")
        return
    
    # 2. Template Selection
    print("\n2. Testing Template Selector...")
    try:
        template_agent = TemplateSelectorAgent()
        template_input = AgentInput(message=prompt, context=context)
        template_result = await template_agent.process(template_input)
        context['template_selection'] = template_result.result
        
        print(f"   ✅ Selected template: {template_result.result.get('template_name')}")
        print(f"   ✅ Score: {template_result.result.get('overall_score', 0):.2f}")
    except Exception as e:
        print(f"   ❌ Template selection failed: {e}")
        return
    
    # 3. Parameter Extraction
    print("\n3. Testing Parameter Extractor...")
    try:
        param_agent = ParameterExtractorAgent()
        param_input = AgentInput(message=prompt, context=context)
        param_result = await param_agent.process(param_input)
        context['parameter_extraction'] = param_result.result
        
        grid = param_result.result.get('grid_parameters', {})
        print(f"   ✅ Grid extracted: {grid.get('nx', 'N/A')}x{grid.get('ny', 'N/A')}x{grid.get('nz', 'N/A')}")
    except Exception as e:
        print(f"   ❌ Parameter extraction failed: {e}")
        return
    
    # 4. Code Generation
    print("\n4. Testing Code Generator...")
    try:
        code_agent = CodeGeneratorAgent()
        code_input = AgentInput(message=prompt, context=context)
        code_result = await code_agent.process(code_input)
        context['code_generation'] = code_result.result
        
        model_code = code_result.result.get('model_code', '')
        print(f"   ✅ Model code generated: {len(model_code)} characters")
        print(f"   ✅ Main code generated: {len(code_result.result.get('main_code', ''))} characters")
    except Exception as e:
        print(f"   ❌ Code generation failed: {e}")
        return
    
    # 5. Validation
    print("\n5. Testing Validator...")
    try:
        validator_agent = ValidatorAgent()
        valid_input = AgentInput(message=prompt, context=context)
        valid_result = await validator_agent.process(valid_input)
        
        print(f"   ✅ Validation complete: {'Valid' if valid_result.result.get('is_valid') else 'Invalid'}")
        issues = valid_result.result.get('issues', [])
        if issues:
            print(f"   ⚠️  Issues found: {len(issues)}")
            for issue in issues[:3]:
                print(f"      - {issue.get('severity', 'unknown')}: {issue.get('message', 'No message')}")
    except Exception as e:
        print(f"   ❌ Validation failed: {e}")
    
    print("\n✅ Basic flow test complete!")


async def test_different_physics():
    """Test different physics types."""
    print("\n\n=== Testing Different Physics Types ===\n")
    
    test_prompts = {
        "dead_oil": "Create a waterflood model",
        "geothermal": "Build a geothermal reservoir model",
        "compositional": "Model CO2 injection with methane",
        "poroelastic": "Create a poroelastic model for reservoir compaction"
    }
    
    intent_agent = IntentClassifierAgent()
    
    for expected_physics, prompt in test_prompts.items():
        print(f"\nTesting: {prompt}")
        try:
            result = await intent_agent.process(AgentInput(message=prompt, context={}))
            detected_physics = result.result.get('physics_type')
            if detected_physics == expected_physics:
                print(f"   ✅ Correctly identified as {detected_physics}")
            else:
                print(f"   ⚠️  Identified as {detected_physics} (expected {expected_physics})")
        except Exception as e:
            print(f"   ❌ Error: {e}")


def test_templates():
    """Test template database."""
    print("\n\n=== Testing Template Database ===\n")
    
    print(f"Total templates available: {len(TEMPLATES)}\n")
    
    # Check golden templates
    golden_templates = [name for name, t in TEMPLATES.items() if t.is_golden]
    print(f"Golden templates (⭐): {len(golden_templates)}")
    for name in golden_templates[:5]:
        print(f"  - {name}: {TEMPLATES[name].description}")
    
    # Check physics coverage
    physics_types = set(t.physics_type for t in TEMPLATES.values())
    print(f"\nPhysics types covered: {len(physics_types)}")
    for ptype in physics_types:
        count = sum(1 for t in TEMPLATES.values() if t.physics_type == ptype)
        print(f"  - {ptype}: {count} templates")


async def test_complex_prompt():
    """Test a complex prompt."""
    print("\n\n=== Testing Complex Prompt ===\n")
    
    prompt = """
    Model a five-spot waterflood pattern in a 2D reservoir (50x50 cells) with 
    three layers of different permeabilities (100, 500, 50 mD). Include 25% 
    porosity and simulate for 5 years.
    """
    
    print(f"Complex prompt: {prompt.strip()}\n")
    
    # Just test intent and parameter extraction
    intent_agent = IntentClassifierAgent()
    param_agent = ParameterExtractorAgent()
    
    # Intent
    intent_result = await intent_agent.process(AgentInput(message=prompt, context={}))
    print("Intent Classification:")
    print(f"  Physics: {intent_result.result.get('physics_type')}")
    print(f"  Grid info: {intent_result.result.get('grid_info')}")
    print(f"  Well info: {intent_result.result.get('well_info')}")
    
    # Parameters
    context = {'intent_classification': intent_result.result}
    param_result = await param_agent.process(AgentInput(message=prompt, context=context))
    print("\nParameter Extraction:")
    print(f"  Grid: {param_result.result.get('grid_parameters')}")
    print(f"  Rock: {param_result.result.get('rock_properties')}")


async def main():
    """Run all tests."""
    print("DARTSGPT System Test Suite")
    print("=" * 50)
    
    # Run tests
    await test_basic_flow()
    await test_different_physics()
    test_templates()
    await test_complex_prompt()
    
    print("\n" + "=" * 50)
    print("All tests complete!")


if __name__ == "__main__":
    asyncio.run(main())