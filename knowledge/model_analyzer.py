"""Analyze DARTS models and extract metadata."""

import ast
import json
from pathlib import Path
from typing import Dict, List, Optional, Set

from pydantic import BaseModel


class ModelMetadata(BaseModel):
    """Metadata for a DARTS model."""

    name: str
    path: str
    physics_type: Optional[str] = None
    grid_type: Optional[str] = None
    features: List[str] = []
    imports: List[str] = []
    well_types: List[str] = []
    parameters: Dict[str, any] = {}
    description: Optional[str] = None
    complexity: str = "basic"  # basic, intermediate, advanced


class DARTSModelAnalyzer:
    """Analyze DARTS models to extract metadata."""

    def __init__(self, models_path: Path):
        self.models_path = models_path
        self.physics_patterns = {
            "Compositional": "compositional",
            "DeadOil": "dead_oil",
            "BlackOil": "black_oil",
            "Geothermal": "geothermal",
            "SinglePhase": "single_phase",
        }
        self.grid_patterns = {
            "StructReservoir": "structured",
            "UnstructReservoir": "unstructured",
            "CPG_Reservoir": "cpg",
        }

    def analyze_model(self, model_path: Path) -> Optional[ModelMetadata]:
        """Analyze a single model file."""
        try:
            with open(model_path, "r") as f:
                content = f.read()

            # Parse the AST
            tree = ast.parse(content)

            metadata = ModelMetadata(
                name=model_path.parent.name,
                path=str(model_path),
            )

            # Extract imports
            metadata.imports = self._extract_imports(tree)

            # Determine physics type
            metadata.physics_type = self._detect_physics_type(content, metadata.imports)

            # Determine grid type
            metadata.grid_type = self._detect_grid_type(content, metadata.imports)

            # Extract features
            metadata.features = self._extract_features(content, metadata.imports)

            # Extract well types
            metadata.well_types = self._extract_well_types(content)

            # Determine complexity
            metadata.complexity = self._assess_complexity(metadata)

            return metadata

        except Exception as e:
            print(f"Error analyzing {model_path}: {e}")
            return None

    def _extract_imports(self, tree: ast.AST) -> List[str]:
        """Extract all imports from the AST."""
        imports = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    imports.append(alias.name)
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ""
                for alias in node.names:
                    imports.append(f"{module}.{alias.name}")
        return imports

    def _detect_physics_type(self, content: str, imports: List[str]) -> Optional[str]:
        """Detect the physics type from imports and content."""
        # Check imports
        for pattern, physics_type in self.physics_patterns.items():
            if any(pattern.lower() in imp.lower() for imp in imports):
                return physics_type

        # Check content
        for pattern, physics_type in self.physics_patterns.items():
            if pattern in content:
                return physics_type

        # Check for specific physics modules
        if "comp_physics" in str(imports):
            return "compositional"
        elif "dead_oil" in str(imports):
            return "dead_oil"
        elif "black_oil" in str(imports):
            return "black_oil"
        elif "geothermal" in str(imports):
            return "geothermal"

        return None

    def _detect_grid_type(self, content: str, imports: List[str]) -> Optional[str]:
        """Detect the grid type from imports and content."""
        for pattern, grid_type in self.grid_patterns.items():
            if pattern in content or any(pattern in imp for imp in imports):
                return grid_type
        return "structured"  # Default

    def _extract_features(self, content: str, imports: List[str]) -> List[str]:
        """Extract features from the model."""
        features = []

        # Thermal
        if "thermal" in content.lower() or any("thermal" in imp.lower() for imp in imports):
            features.append("thermal")

        # Mechanical
        if "mech" in content.lower() or any("mech" in imp.lower() for imp in imports):
            features.append("mechanical")

        # Chemical
        if "chemical" in content.lower() or "phreeqc" in content.lower():
            features.append("chemical")

        # Solid
        if "solid" in content.lower():
            features.append("solid")

        # Foam
        if "foam" in content.lower():
            features.append("foam")

        # Adjoint
        if "adjoint" in content.lower():
            features.append("adjoint")

        # MPFA
        if "mpfa" in content.lower():
            features.append("mpfa")

        return features

    def _extract_well_types(self, content: str) -> List[str]:
        """Extract well types used in the model."""
        well_types = []

        if "add_well(" in content:
            well_types.append("vertical")

        if "injector" in content.lower():
            well_types.append("injection")

        if "producer" in content.lower():
            well_types.append("production")

        if "rate" in content.lower():
            well_types.append("rate_control")

        if "bhp" in content.lower():
            well_types.append("bhp_control")

        return well_types

    def _assess_complexity(self, metadata: ModelMetadata) -> str:
        """Assess the complexity of the model."""
        feature_count = len(metadata.features)

        if feature_count >= 3 or "adjoint" in metadata.features:
            return "advanced"
        elif feature_count >= 1 or metadata.grid_type != "structured":
            return "intermediate"
        else:
            return "basic"

    def analyze_all_models(self) -> List[ModelMetadata]:
        """Analyze all models in the directory."""
        models = []

        for model_dir in self.models_path.iterdir():
            if model_dir.is_dir() and not model_dir.name.startswith("_"):
                model_file = model_dir / "model.py"
                if model_file.exists():
                    metadata = self.analyze_model(model_file)
                    if metadata:
                        models.append(metadata)

        return models

    def save_metadata(self, output_path: Path):
        """Save metadata to JSON file."""
        models = self.analyze_all_models()
        
        output_data = {
            "models": [model.model_dump() for model in models],
            "summary": {
                "total_models": len(models),
                "physics_types": list(set(m.physics_type for m in models if m.physics_type)),
                "grid_types": list(set(m.grid_type for m in models if m.grid_type)),
                "features": list(set(f for m in models for f in m.features)),
            }
        }

        with open(output_path, "w") as f:
            json.dump(output_data, f, indent=2)

        return output_data