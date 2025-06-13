"""Test script for knowledge base initialization."""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))

from dartsgpt.knowledge.rag_system import DARTSKnowledgeBase
from dartsgpt.knowledge.model_analyzer import DARTSModelAnalyzer
from dartsgpt.config import settings


def test_model_analyzer():
    """Test the model analyzer."""
    print("Testing Model Analyzer...")
    print(f"Models path: {settings.darts_models_path}")
    
    analyzer = DARTSModelAnalyzer(settings.darts_models_path)
    
    # Test analyzing a single model
    test_model = settings.darts_models_path / "2ph_comp" / "model.py"
    if test_model.exists():
        metadata = analyzer.analyze_model(test_model)
        if metadata:
            print(f"\nAnalyzed {metadata.name}:")
            print(f"  Physics: {metadata.physics_type}")
            print(f"  Grid: {metadata.grid_type}")
            print(f"  Features: {metadata.features}")
            print(f"  Complexity: {metadata.complexity}")
    
    # Analyze all models
    print("\nAnalyzing all models...")
    all_models = analyzer.analyze_all_models()
    print(f"Found {len(all_models)} models")
    
    # Group by physics type
    physics_types = {}
    for model in all_models:
        if model.physics_type:
            if model.physics_type not in physics_types:
                physics_types[model.physics_type] = []
            physics_types[model.physics_type].append(model.name)
    
    print("\nModels by physics type:")
    for physics, models in physics_types.items():
        print(f"  {physics}: {models}")


def test_knowledge_base():
    """Test the knowledge base initialization."""
    print("\n\nTesting Knowledge Base...")
    
    kb = DARTSKnowledgeBase()
    
    # Initialize (this will build the vector store)
    kb.initialize(force_rebuild=True)
    
    # Test searching
    print("\nTesting search...")
    
    # Search for CO2 injection models
    results = kb.search("CO2 injection carbon storage")
    print(f"\nSearch for 'CO2 injection': {len(results)} results")
    for i, doc in enumerate(results[:3]):
        print(f"  {i+1}. {doc.metadata.get('model_name', 'Unknown')} - {doc.metadata.get('type', 'Unknown')}")
    
    # Search for geothermal models
    results = kb.search_models_by_physics("geothermal")
    print(f"\nGeothermal models: {len(results)} results")
    for i, doc in enumerate(results[:3]):
        print(f"  {i+1}. {doc.metadata.get('model_name', 'Unknown')}")
    
    # Search for templates
    results = kb.search_templates("water injection")
    print(f"\nTemplate search: {len(results)} results")
    for i, doc in enumerate(results[:3]):
        print(f"  {i+1}. {doc.metadata.get('template_name', 'Unknown')}")


if __name__ == "__main__":
    print("DARTSGPT Knowledge Base Test")
    print("=" * 50)
    
    # Check if paths exist
    if not settings.darts_models_path.exists():
        print(f"ERROR: DARTS models path not found: {settings.darts_models_path}")
        print("Please update the path in your .env file")
        sys.exit(1)
    
    test_model_analyzer()
    test_knowledge_base()
    
    print("\n\nTest completed!")