"""Template database for DARTS models."""

from pathlib import Path
from typing import Dict, List, Optional

from pydantic import BaseModel


class TemplateInfo(BaseModel):
    """Information about a DARTS template."""

    name: str
    path: str
    physics_type: str
    description: str
    features: List[str] = []
    complexity: str = "basic"
    example_use_cases: List[str] = []
    required_parameters: List[str] = []
    optional_parameters: List[str] = []
    grid_type: str = "structured"
    use_cases: List[str] = []
    is_golden: bool = False
    parameters: Optional[Dict] = None


class TemplateDatabase:
    """Database of DARTS model templates."""

    def __init__(self):
        self.templates = self._initialize_templates()

    def _initialize_templates(self) -> Dict[str, TemplateInfo]:
        """Initialize the template database with known templates."""
        return {
            # Compositional Physics Templates
            "2ph_comp": TemplateInfo(
                name="2ph_comp",
                path="2ph_comp/model.py",
                physics_type="compositional",
                description="Basic two-phase compositional flow model",
                features=[],
                complexity="basic",
                example_use_cases=[
                    "CO2 injection in saline aquifer",
                    "Natural gas storage",
                    "Miscible gas injection",
                ],
                required_parameters=["nx", "ny", "nz", "dx", "dy", "dz", "poro", "perm"],
                optional_parameters=["rock_comp", "initial_pressure", "temperature"],
                is_golden=True,
                use_cases=["CO2 storage", "Gas injection", "EOR"],
            ),
            "2ph_comp_solid": TemplateInfo(
                name="2ph_comp_solid",
                path="2ph_comp_solid/model.py",
                physics_type="compositional",
                description="Two-phase compositional flow with solid precipitation",
                features=["solid", "precipitation"],
                complexity="intermediate",
                example_use_cases=[
                    "CO2 storage with salt precipitation",
                    "Scale formation modeling",
                ],
                required_parameters=["nx", "ny", "nz", "dx", "dy", "dz", "poro", "perm", "solid_density"],
                optional_parameters=["solid_solubility", "kinetic_rate"],
            ),
            "3ph_comp_w": TemplateInfo(
                name="3ph_comp_w",
                path="3ph_comp_w/model.py",
                physics_type="compositional",
                description="Three-phase compositional flow with water",
                features=["three_phase", "water"],
                complexity="intermediate",
                example_use_cases=[
                    "WAG injection",
                    "Three-phase compositional reservoir",
                ],
                required_parameters=["nx", "ny", "nz", "dx", "dy", "dz", "poro", "perm"],
                optional_parameters=["initial_water_sat", "water_properties"],
            ),
            
            # Dead Oil Templates
            "2ph_do": TemplateInfo(
                name="2ph_do",
                path="2ph_do/model.py",
                physics_type="dead_oil",
                description="Basic two-phase dead oil model",
                features=[],
                complexity="basic",
                example_use_cases=[
                    "Water flooding",
                    "Simple oil reservoir simulation",
                ],
                required_parameters=["nx", "ny", "nz", "dx", "dy", "dz", "poro", "perm"],
                optional_parameters=["oil_api", "initial_oil_sat"],
            ),
            "2ph_do_thermal": TemplateInfo(
                name="2ph_do_thermal",
                path="2ph_do_thermal/model.py",
                physics_type="dead_oil",
                description="Two-phase dead oil with thermal effects",
                features=["thermal"],
                complexity="intermediate",
                example_use_cases=[
                    "Steam injection",
                    "Thermal recovery",
                ],
                required_parameters=["nx", "ny", "nz", "dx", "dy", "dz", "poro", "perm", "rock_heat_capacity"],
                optional_parameters=["thermal_conductivity", "initial_temperature"],
            ),
            "3ph_do": TemplateInfo(
                name="3ph_do",
                path="3ph_do/model.py",
                physics_type="dead_oil",
                description="Three-phase dead oil model",
                features=["three_phase"],
                complexity="intermediate",
                example_use_cases=[
                    "Gas cap reservoirs",
                    "Solution gas drive",
                ],
                required_parameters=["nx", "ny", "nz", "dx", "dy", "dz", "poro", "perm"],
                optional_parameters=["initial_gas_sat", "bubble_point"],
            ),
            
            # Black Oil Templates
            "3ph_bo": TemplateInfo(
                name="3ph_bo",
                path="3ph_bo/model.py",
                physics_type="black_oil",
                description="Three-phase black oil model",
                features=["pvt", "three_phase"],
                complexity="intermediate",
                example_use_cases=[
                    "Conventional oil reservoir",
                    "Gas injection",
                    "Pressure depletion",
                ],
                required_parameters=["nx", "ny", "nz", "dx", "dy", "dz", "poro", "perm"],
                optional_parameters=["pvt_table", "initial_solution_gas"],
            ),
            
            # Geothermal Templates
            "2ph_geothermal": TemplateInfo(
                name="2ph_geothermal",
                path="2ph_geothermal/model.py",
                physics_type="geothermal",
                description="Two-phase water/steam geothermal model",
                features=["thermal", "phase_change"],
                complexity="intermediate",
                example_use_cases=[
                    "Geothermal reservoir",
                    "Steam production",
                    "Hot water injection",
                ],
                required_parameters=["nx", "ny", "nz", "dx", "dy", "dz", "poro", "perm", "rock_heat_capacity"],
                optional_parameters=["thermal_conductivity", "initial_temperature", "initial_pressure"],
            ),
            "2ph_geothermal_mass_flux": TemplateInfo(
                name="2ph_geothermal_mass_flux",
                path="2ph_geothermal_mass_flux/model.py",
                physics_type="geothermal",
                description="Geothermal model with mass flux boundary conditions",
                features=["thermal", "mass_flux", "boundary_conditions"],
                complexity="advanced",
                example_use_cases=[
                    "Natural convection",
                    "Geothermal systems with recharge",
                ],
                required_parameters=["nx", "ny", "nz", "dx", "dy", "dz", "poro", "perm", "rock_heat_capacity"],
                optional_parameters=["boundary_flux", "boundary_temperature"],
            ),
            "CoaxWell": TemplateInfo(
                name="CoaxWell",
                path="CoaxWell/model.py",
                physics_type="geothermal",
                description="Coaxial well heat exchanger model",
                features=["thermal", "wellbore", "heat_exchange"],
                complexity="advanced",
                example_use_cases=[
                    "Closed-loop geothermal",
                    "Coaxial heat exchanger",
                ],
                required_parameters=["well_length", "inner_radius", "outer_radius", "flow_rate"],
                optional_parameters=["pipe_conductivity", "fluid_properties"],
            ),
            
            # Advanced Templates
            "cpg_sloping_fault": TemplateInfo(
                name="cpg_sloping_fault",
                path="cpg_sloping_fault/model_cpg.py",
                physics_type="flexible",
                description="Corner-point grid with fault - supports multiple physics",
                features=["cpg", "fault", "flexible_physics"],
                complexity="advanced",
                example_use_cases=[
                    "Fault seal analysis",
                    "Complex geology",
                    "Multi-physics simulation",
                ],
                required_parameters=["grid_file", "physics_type"],
                optional_parameters=["fault_properties", "burden_layers"],
            ),
            
            # Specialized Templates
            "Chem_benchmark_new": TemplateInfo(
                name="Chem_benchmark_new",
                path="Chem_benchmark_new/model.py",
                physics_type="compositional",
                description="Chemical reaction benchmark model",
                features=["chemical", "reactions"],
                complexity="advanced",
                example_use_cases=[
                    "CO2 mineralization",
                    "Chemical EOR",
                ],
                required_parameters=["nx", "ny", "nz", "dx", "dy", "dz", "poro", "perm"],
                optional_parameters=["reaction_rates", "mineral_composition"],
            ),
            "CO2_foam_CCS": TemplateInfo(
                name="CO2_foam_CCS",
                path="CO2_foam_CCS/model.py",
                physics_type="compositional",
                description="CO2 foam for carbon capture and storage",
                features=["foam", "ccs", "mobility_control"],
                complexity="advanced",
                example_use_cases=[
                    "CO2 storage with foam",
                    "Mobility control in CCS",
                ],
                required_parameters=["nx", "ny", "nz", "dx", "dy", "dz", "poro", "perm"],
                optional_parameters=["foam_parameters", "surfactant_concentration"],
            ),
        }

    def get_template(self, name: str) -> Optional[TemplateInfo]:
        """Get a template by name."""
        return self.templates.get(name)

    def get_templates_by_physics(self, physics_type: str) -> List[TemplateInfo]:
        """Get all templates for a physics type."""
        return [t for t in self.templates.values() if t.physics_type == physics_type]

    def get_templates_by_feature(self, feature: str) -> List[TemplateInfo]:
        """Get all templates with a specific feature."""
        return [t for t in self.templates.values() if feature in t.features]

    def get_templates_by_complexity(self, complexity: str) -> List[TemplateInfo]:
        """Get all templates of a specific complexity."""
        return [t for t in self.templates.values() if t.complexity == complexity]

    def find_best_template(
        self,
        physics_type: str,
        features: List[str] = [],
        complexity: str = "basic"
    ) -> Optional[TemplateInfo]:
        """Find the best matching template."""
        candidates = self.get_templates_by_physics(physics_type)
        
        if not candidates:
            # Try flexible templates
            candidates = [t for t in self.templates.values() if t.physics_type == "flexible"]
        
        # Filter by features
        if features:
            scored_candidates = []
            for template in candidates:
                score = sum(1 for f in features if f in template.features)
                scored_candidates.append((score, template))
            
            # Sort by score (descending)
            scored_candidates.sort(key=lambda x: x[0], reverse=True)
            
            # Return the best match if it has at least one required feature
            if scored_candidates and scored_candidates[0][0] > 0:
                return scored_candidates[0][1]
        
        # Filter by complexity
        complexity_candidates = [t for t in candidates if t.complexity == complexity]
        if complexity_candidates:
            return complexity_candidates[0]
        
        # Return simplest available template
        for comp in ["basic", "intermediate", "advanced"]:
            comp_candidates = [t for t in candidates if t.complexity == comp]
            if comp_candidates:
                return comp_candidates[0]
        
        return None

# Create global instance
template_db = TemplateDatabase()
TEMPLATES = template_db.templates
