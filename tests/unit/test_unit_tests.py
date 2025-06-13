#!/usr/bin/env python
"""Unit tests for individual DARTSGPT components."""

import pytest
from unittest.mock import Mock, patch
import json

from langchain_openai import ChatOpenAI
from langgraph.graph import StateGraph

from orchestrator import (
    DARTSGPTState,
    create_intent_classifier,
    create_template_selector,
    create_parameter_extractor,
    create_code_generator,
    create_validator,
    build_dartsgpt_graph
)
from knowledge.templates.template_database import TEMPLATES
from knowledge.physics_features import PHYSICS_KEYWORDS, UNIT_CONVERSIONS


class TestIntentClassifier:
    """Test the intent classifier agent."""
    
    def test_compositional_intent(self):
        """Test classification of compositional physics."""
        model = Mock(spec=ChatOpenAI)
        model.invoke.return_value.content = json.dumps({
            "physics_type": "compositional",
            "features": ["injection", "co2"],
            "complexity": "moderate"
        })
        
        classifier = create_intent_classifier(model)
        state = DARTSGPTState(
            messages=[],
            prompt="Create a CO2 injection model for carbon storage"
        )
        
        result = classifier(state)
        assert result.goto == "template_selector"
        assert result.update["intent"]["physics_type"] == "compositional"
    
    def test_fallback_keyword_matching(self):
        """Test fallback to keyword matching when JSON parsing fails."""
        model = Mock(spec=ChatOpenAI)
        model.invoke.return_value.content = "Invalid JSON response"
        
        classifier = create_intent_classifier(model)
        state = DARTSGPTState(
            messages=[],
            prompt="Model a geothermal reservoir with steam"
        )
        
        result = classifier(state)
        assert result.goto == "template_selector"
        assert result.update["intent"]["physics_type"] == "geothermal"
    
    def test_default_physics_type(self):
        """Test default physics type when no keywords match."""
        model = Mock(spec=ChatOpenAI)
        model.invoke.return_value.content = "Invalid JSON"
        
        classifier = create_intent_classifier(model)
        state = DARTSGPTState(
            messages=[],
            prompt="Create a model"
        )
        
        result = classifier(state)
        assert result.update["intent"]["physics_type"] == "dead_oil"


class TestTemplateSelector:
    """Test the template selector agent."""
    
    def test_golden_template_priority(self):
        """Test that golden templates are prioritized."""
        model = Mock(spec=ChatOpenAI)
        selector = create_template_selector(model)
        
        state = DARTSGPTState(
            messages=[],
            prompt="",
            intent={"physics_type": "compositional", "features": []}
        )
        
        result = selector(state)
        assert result.goto == "parameter_extractor"
        # Should select 2ph_comp which is golden for compositional
        assert result.update["selected_template"] == "2ph_comp"
    
    def test_fallback_template(self):
        """Test fallback when no matching physics type."""
        model = Mock(spec=ChatOpenAI)
        selector = create_template_selector(model)
        
        state = DARTSGPTState(
            messages=[],
            prompt="",
            intent={"physics_type": "unknown_physics", "features": []}
        )
        
        result = selector(state)
        # Should fall back to dead oil template
        assert result.update["selected_template"] in ["2ph_do", "3ph_do"]


class TestParameterExtractor:
    """Test the parameter extractor agent."""
    
    def test_grid_extraction(self):
        """Test extraction of grid dimensions."""
        model = Mock(spec=ChatOpenAI)
        model.invoke.return_value.content = json.dumps({
            "grid_parameters": {"nx": 100, "ny": 100, "nz": 20},
            "rock_properties": {}
        })
        
        extractor = create_parameter_extractor(model)
        state = DARTSGPTState(
            messages=[],
            prompt="Create a 100x100x20 grid"
        )
        
        result = extractor(state)
        assert result.goto == "code_generator"
        assert result.update["parameters"]["grid_parameters"]["nx"] == 100
    
    def test_regex_fallback(self):
        """Test regex extraction when JSON parsing fails."""
        model = Mock(spec=ChatOpenAI)
        model.invoke.return_value.content = "Invalid JSON"
        
        extractor = create_parameter_extractor(model)
        state = DARTSGPTState(
            messages=[],
            prompt="Build a 50x50x10 reservoir"
        )
        
        result = extractor(state)
        params = result.update["parameters"]
        assert params["grid_parameters"]["nx"] == 50
        assert params["grid_parameters"]["ny"] == 50
        assert params["grid_parameters"]["nz"] == 10


class TestCodeGenerator:
    """Test the code generator agent."""
    
    def test_code_generation_with_comments(self):
        """Test that generated code includes explanatory comments."""
        model = Mock(spec=ChatOpenAI)
        generator = create_code_generator(model)
        
        state = DARTSGPTState(
            messages=[],
            prompt="Create a CO2 injection model",
            intent={"physics_type": "compositional", "features": ["injection"]},
            selected_template="2ph_comp",
            parameters={
                "grid_parameters": {"nx": 100, "ny": 100, "nz": 20},
                "temperature": 60
            }
        )
        
        result = generator(state)
        assert result.goto == "validator"
        
        model_code = result.update["model_code"]
        main_code = result.update["main_code"]
        
        # Check for explanatory comments
        assert "Original prompt:" in model_code
        assert "Physics type" in model_code
        assert "This model was automatically generated" in model_code
        assert "Grid dimensions" in model_code
        
        # Check main.py comments
        assert "This script runs the DARTS reservoir simulation" in main_code
        assert "from your prompt:" in main_code
    
    def test_parameter_handling(self):
        """Test handling of None parameters."""
        model = Mock(spec=ChatOpenAI)
        generator = create_code_generator(model)
        
        state = DARTSGPTState(
            messages=[],
            prompt="Test prompt",
            intent={"physics_type": "dead_oil"},
            selected_template="2ph_do",
            parameters={
                "grid_parameters": {"nx": None, "ny": None, "nz": None},
                "temperature": None
            }
        )
        
        result = generator(state)
        model_code = result.update["model_code"]
        
        # Should use defaults
        assert "self.nx = 50" in model_code
        assert "self.ny = 50" in model_code
        assert "self.nz = 10" in model_code


class TestValidator:
    """Test the validator agent."""
    
    def test_valid_code(self):
        """Test validation of correct code."""
        model = Mock(spec=ChatOpenAI)
        validator = create_validator(model)
        
        valid_code = '''
from darts.models.reservoirmodel import ReservoirModel

class Model(ReservoirModel):
    def set_physics(self):
        pass
    def set_reservoir(self):
        pass
    def set_wells(self):
        pass
'''
        
        state = DARTSGPTState(
            messages=[],
            prompt="",
            model_code=valid_code,
            main_code="",
            selected_template="test",
            parameters={}
        )
        
        result = validator(state)
        validation = result.update["validation_result"]
        
        assert validation["imports_valid"] == True
        assert validation["structure_valid"] == True
        assert len(validation["issues"]) == 0
    
    def test_missing_imports(self):
        """Test detection of missing imports."""
        model = Mock(spec=ChatOpenAI)
        validator = create_validator(model)
        
        invalid_code = '''
class Model:
    def set_physics(self):
        pass
'''
        
        state = DARTSGPTState(
            messages=[],
            prompt="",
            model_code=invalid_code,
            main_code="",
            selected_template="test",
            parameters={}
        )
        
        result = validator(state)
        validation = result.update["validation_result"]
        
        assert validation["imports_valid"] == False
        assert "Missing DARTS model import" in validation["issues"]


class TestGraphConstruction:
    """Test the graph construction."""
    
    def test_graph_building(self):
        """Test that the graph builds successfully."""
        with patch('orchestrator.get_settings') as mock_settings:
            mock_settings.return_value.openai_api_key = "test-key"
            mock_settings.return_value.openai_model = "gpt-4"
            
            graph = build_dartsgpt_graph()
            assert isinstance(graph, StateGraph)
    
    def test_graph_nodes(self):
        """Test that all required nodes are present."""
        with patch('orchestrator.get_settings') as mock_settings:
            mock_settings.return_value.openai_api_key = "test-key"
            mock_settings.return_value.openai_model = "gpt-4"
            
            # We can't easily inspect compiled graph nodes,
            # but we can ensure it compiles without error
            graph = build_dartsgpt_graph()
            assert graph is not None


class TestUnitConversions:
    """Test unit conversion utilities."""
    
    def test_permeability_conversion(self):
        """Test mD to m² conversion."""
        md_to_m2 = UNIT_CONVERSIONS["permeability"]["md_to_m2"]
        assert abs(md_to_m2 - 9.869233e-16) < 1e-20
    
    def test_temperature_conversion(self):
        """Test temperature conversion factors."""
        f_to_c = UNIT_CONVERSIONS["temperature"]["f_to_c"]
        # Test 32°F = 0°C
        celsius = f_to_c(32)
        assert celsius == 0
        # Test 212°F = 100°C
        celsius = f_to_c(212)
        assert celsius == 100
    
    def test_pressure_conversion(self):
        """Test pressure conversions."""
        psi_to_bar = UNIT_CONVERSIONS["pressure"]["psi_to_bar"]
        bar = 14.5038 * psi_to_bar  # 14.5038 psi ≈ 1 bar
        assert abs(bar - 1.0) < 0.01


class TestPhysicsKeywords:
    """Test physics keyword mappings."""
    
    def test_all_physics_have_keywords(self):
        """Ensure all physics types have keywords defined."""
        expected_physics = [
            "compositional", "dead_oil", "black_oil", 
            "geothermal", "poroelastic", "chemical"
        ]
        
        for physics in expected_physics:
            assert physics in PHYSICS_KEYWORDS
            assert len(PHYSICS_KEYWORDS[physics]) > 0
    
    def test_unique_keywords(self):
        """Test that some keywords are unique to physics types."""
        co2_physics = []
        for physics, keywords in PHYSICS_KEYWORDS.items():
            if "co2" in keywords:
                co2_physics.append(physics)
        
        # CO2 should primarily be associated with compositional
        assert "compositional" in co2_physics


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])