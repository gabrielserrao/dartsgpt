"""Test cases for agent evaluation."""

from typing import Dict, List, Optional
from pydantic import BaseModel, Field


class TestCase(BaseModel):
    """A single test case for agent evaluation."""
    
    id: str = Field(description="Unique identifier for the test case")
    name: str = Field(description="Human-readable name")
    prompt: str = Field(description="Natural language prompt")
    expected_physics: str = Field(description="Expected physics type")
    expected_template: str = Field(description="Expected template selection")
    expected_features: List[str] = Field(default_factory=list, description="Expected features")
    expected_parameters: Dict[str, any] = Field(default_factory=dict, description="Expected parameters")
    ground_truth: Optional[str] = Field(default=None, description="Ground truth code snippet")
    keywords: List[str] = Field(default_factory=list, description="Keywords that should be recognized")
    complexity: str = Field(default="basic", description="Expected complexity level")


class TestSuite(BaseModel):
    """Collection of test cases."""
    
    name: str
    description: str
    test_cases: List[TestCase]
    
    @classmethod
    def load_default_suite(cls) -> "TestSuite":
        """Load the default test suite."""
        return cls(
            name="DARTSGPT Default Test Suite",
            description="Comprehensive test cases for all physics types",
            test_cases=[
                # Compositional Flow Tests
                TestCase(
                    id="comp_001",
                    name="Basic CO2 Injection",
                    prompt="Create a CO2 injection model for carbon storage in a saline aquifer",
                    expected_physics="compositional",
                    expected_template="2ph_comp",
                    expected_features=[],
                    expected_parameters={
                        "nx": 50, "ny": 50, "nz": 10,
                        "well_type": "injection",
                        "fluid": "CO2"
                    },
                    keywords=["CO2", "injection", "carbon storage", "saline aquifer"],
                    complexity="basic"
                ),
                TestCase(
                    id="comp_002",
                    name="CO2 with Salt Precipitation",
                    prompt="Model CO2 injection with salt precipitation in a deep saline formation",
                    expected_physics="compositional",
                    expected_template="2ph_comp_solid",
                    expected_features=["solid", "precipitation"],
                    expected_parameters={
                        "solid_type": "salt",
                        "precipitation": True
                    },
                    keywords=["CO2", "salt", "precipitation", "solid"],
                    complexity="intermediate"
                ),
                
                # Dead Oil Tests
                TestCase(
                    id="do_001",
                    name="Water Flooding",
                    prompt="Create a simple water flooding model for oil recovery",
                    expected_physics="dead_oil",
                    expected_template="2ph_do",
                    expected_features=[],
                    expected_parameters={
                        "injection_fluid": "water",
                        "recovery_type": "secondary"
                    },
                    keywords=["water", "flooding", "oil", "recovery"],
                    complexity="basic"
                ),
                TestCase(
                    id="do_002",
                    name="Steam Injection",
                    prompt="Model steam injection for heavy oil recovery with thermal effects",
                    expected_physics="dead_oil",
                    expected_template="2ph_do_thermal",
                    expected_features=["thermal"],
                    expected_parameters={
                        "injection_fluid": "steam",
                        "thermal": True
                    },
                    keywords=["steam", "injection", "thermal", "heavy oil"],
                    complexity="intermediate"
                ),
                
                # Geothermal Tests
                TestCase(
                    id="geo_001",
                    name="Basic Geothermal",
                    prompt="Create a geothermal reservoir model with production and injection wells",
                    expected_physics="geothermal",
                    expected_template="2ph_geothermal",
                    expected_features=["thermal", "phase_change"],
                    expected_parameters={
                        "temperature_range": "high",
                        "phase_behavior": "two_phase"
                    },
                    keywords=["geothermal", "reservoir", "production", "injection"],
                    complexity="intermediate"
                ),
                TestCase(
                    id="geo_002",
                    name="Coaxial Well",
                    prompt="Model a closed-loop coaxial well heat exchanger for geothermal energy",
                    expected_physics="geothermal",
                    expected_template="CoaxWell",
                    expected_features=["thermal", "wellbore", "heat_exchange"],
                    expected_parameters={
                        "well_type": "coaxial",
                        "closed_loop": True
                    },
                    keywords=["coaxial", "well", "heat exchanger", "closed-loop"],
                    complexity="advanced"
                ),
                
                # Black Oil Tests
                TestCase(
                    id="bo_001",
                    name="Three-Phase Black Oil",
                    prompt="Create a three-phase black oil model for a conventional oil reservoir",
                    expected_physics="black_oil",
                    expected_template="3ph_bo",
                    expected_features=["pvt", "three_phase"],
                    expected_parameters={
                        "phases": ["oil", "gas", "water"],
                        "pvt_model": "black_oil"
                    },
                    keywords=["black oil", "three-phase", "conventional", "reservoir"],
                    complexity="intermediate"
                ),
                
                # Advanced Tests
                TestCase(
                    id="adv_001",
                    name="Chemical EOR",
                    prompt="Model chemical enhanced oil recovery with reactions and adsorption",
                    expected_physics="compositional",
                    expected_template="Chem_benchmark_new",
                    expected_features=["chemical", "reactions"],
                    expected_parameters={
                        "chemical_type": "surfactant",
                        "reactions": True,
                        "adsorption": True
                    },
                    keywords=["chemical", "EOR", "reactions", "adsorption"],
                    complexity="advanced"
                ),
                TestCase(
                    id="adv_002",
                    name="Poroelastic Analysis",
                    prompt="Create a poroelastic model to analyze stress changes during injection",
                    expected_physics="single_phase",
                    expected_template="1ph_1comp_poroelastic_analytics",
                    expected_features=["mechanical"],
                    expected_parameters={
                        "coupling": "poroelastic",
                        "stress_analysis": True
                    },
                    keywords=["poroelastic", "stress", "mechanical", "coupling"],
                    complexity="advanced"
                ),
            ]
        )
    
    def get_test_by_physics(self, physics_type: str) -> List[TestCase]:
        """Get all test cases for a specific physics type."""
        return [tc for tc in self.test_cases if tc.expected_physics == physics_type]
    
    def get_test_by_complexity(self, complexity: str) -> List[TestCase]:
        """Get all test cases of a specific complexity."""
        return [tc for tc in self.test_cases if tc.complexity == complexity]