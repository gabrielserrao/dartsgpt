# DARTSGPT: AI-Powered DARTS Model Generation System

<div align="center">
  
  [![Python Version](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
  [![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
  [![Documentation](https://img.shields.io/badge/docs-comprehensive-brightgreen.svg)](docs/)
  [![Architecture](https://img.shields.io/badge/architecture-multi--agent-purple.svg)](docs/DARTSGPT_Architecture_Diagrams.md)
  [![Status](https://img.shields.io/badge/status-production%20ready-success.svg)](docs/PROJECT_STATUS.md)
</div>

## 🚀 Overview

**DARTSGPT** is an innovative AI-powered system that automatically generates [DARTS (Delft Advanced Research Terra Simulator)](https://darts.citg.tudelft.nl/) models from natural language descriptions. Using state-of-the-art Large Language Models and a sophisticated multi-agent architecture built with LangChain and LangGraph, DARTSGPT transforms complex reservoir engineering requirements into ready-to-run simulation code with detailed explanatory comments.

### 🎯 Project Goal

The ultimate goal of DARTSGPT is to democratize reservoir simulation by allowing users to create complex DARTS models using natural language descriptions instead of manual coding. The system understands the DARTS codebase, executes simulations, and generates complete Python scripts with all necessary imports and configurations.

### 🌟 Key Features

- **Natural Language Interface**: Describe your reservoir model in plain English
- **Multi-Agent Orchestration**: Specialized agents work together using LangGraph
- **Comprehensive Physics Support**: All DARTS physics types (compositional, dead_oil, black_oil, geothermal, poroelastic, chemical)
- **Smart Code Generation**: Includes detailed comments explaining rationale based on your prompts
- **Parameter Extraction**: Extracts numerical parameters with automatic unit conversion
- **Rich CLI**: Beautiful terminal interface with progress indicators and formatted output
- **Extensive Templates**: Pre-validated model templates for various simulation scenarios
- **Evaluation Framework**: Comprehensive testing with vector store and ground truth evaluation

## 📋 Table of Contents

- [Features](#-features)
- [Architecture](#-architecture)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Usage Examples](#-usage-examples)
- [Physics Types](#-physics-types)
- [Templates & Examples](#-templates--examples)
- [Agent System](#-agent-system)
- [Development](#-development)
- [Testing](#-testing)
- [Contributing](#-contributing)
- [Project Status](#-project-status)
- [License](#-license)

## 🌟 Features

### Supported Physics Types

| Physics Type | Description | Templates | Use Cases |
|-------------|-------------|-----------|-----------|
| **Compositional** | Multi-component, multi-phase flow | `2ph_comp`, `2ph_comp_solid`, `3ph_comp_w` | CO2 storage, gas injection, EOR |
| **Dead Oil** | Simplified two/three-phase flow | `2ph_do`, `2ph_do_thermal`, `3ph_do` | Water flooding, steam injection |
| **Black Oil** | Industry-standard PVT model | `3ph_bo` | Conventional reservoirs |
| **Geothermal** | Water/steam systems | `2ph_geothermal`, `CoaxWell` | Geothermal energy extraction |
| **Poroelastic** | Coupled flow and mechanics | `1ph_1comp_poroelastic` | Stress analysis, subsidence |
| **Chemical** | Reactive transport | `Chem_benchmark`, `phreeqc_dissolution` | Chemical EOR, reactive flow |
| **Specialized** | Advanced physics | `CO2_foam_CCS`, `SPE11b` | CCS, foam flooding, benchmarks |

### Natural Language Understanding

DARTSGPT understands various ways to describe reservoir models:

```text
Simple: "Create a CO2 injection model"
Detailed: "Model CO2 injection into a 50x50x10 saline aquifer at 1000m depth with 20% porosity"
Technical: "Build a compositional model with WAG injection using 3-phase relative permeability"
Complex: "Create a geothermal doublet system with 180°C reservoir temperature and 150 kg/s production rate"
```

## 🏗️ Architecture

DARTSGPT uses a sophisticated multi-agent architecture powered by LangGraph:

```
┌─────────────────────────────────────────────────────────┐
│                    User Input                           │
│        "Create a geothermal reservoir model..."         │
└─────────────────────────┬───────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│                 🎯 Supervisor Agent                     │
│         Orchestrates the entire workflow                │
└─────────────────────────┬───────────────────────────────┘
                          │
        ┌─────────────────┴─────────────────┐
        │                                   │
        ▼                                   ▼
┌──────────────────┐              ┌──────────────────┐
│ 🔍 Intent        │              │ 📋 Template      │
│    Classifier    │              │    Selector      │
└──────────────────┘              └──────────────────┘
        │                                   │
        ▼                                   ▼
┌──────────────────┐              ┌──────────────────┐
│ 🔢 Parameter     │              │ 💻 Code          │
│    Extractor     │              │    Generator     │
└──────────────────┘              └──────────────────┘
                          │
                          ▼
                ┌──────────────────┐
                │ ✅ Validator     │
                └──────────────────┘
                          │
                          ▼
                ┌──────────────────┐
                │  Generated Files │
                │ • model.py       │
                │ • main.py        │
                │ • metadata.json  │
                └──────────────────┘
```

### Agent Responsibilities

1. **Supervisor Agent**: Central coordinator using LangGraph's StateGraph
2. **Intent Classifier**: Analyzes prompts to identify physics type and features
3. **Template Selector**: Matches intent to available templates with golden priority
4. **Parameter Extractor**: Extracts values with unit conversion (mD→m², °F→°C)
5. **Code Generator**: Creates complete models with explanatory comments
6. **Validator**: Ensures code correctness and structure

## 💻 Installation

### Prerequisites

- Python 3.11 or higher
- Git
- OpenAI API key (for LLM access)

### Step 1: Clone the Repository

```bash
git clone https://github.com/yourusername/dartsgpt.git
cd dartsgpt
```

### Step 2: Install Dependencies

We recommend using `uv` for fast package management:

```bash
# Install dependencies with uv
uv sync

# Or use pip with requirements.txt
pip install -r requirements.txt
```

### Step 3: Configure Environment

Create a `.env` file in the project root:

```env
# Required
OPENAI_API_KEY=your-api-key-here
OPENAI_MODEL=gpt-4o-mini

# Optional
LANGCHAIN_TRACING_V2=false
LOG_LEVEL=INFO
```

## 🚀 Quick Start

### Streamlit Web Interface

DARTSGPT includes a beautiful web interface built with Streamlit that provides an interactive way to generate DARTS models.

#### Starting the Web App

```bash
uv run streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`

#### Features

- **🎨 Beautiful UI**: Modern, responsive design with gradient backgrounds and smooth animations
- **🎯 Real-time Agent Visualization**: Watch as each agent processes your request in real-time
- **📊 Interactive Dashboard**: Control panel with model selection, creativity settings, and statistics
- **📚 Template Browser**: Browse all available templates organized by physics type
- **📜 Generation History**: Track all your generated models with easy access to code
- **💾 Export Options**: Download generated files or preview them directly in the browser
- **📖 Built-in Documentation**: Comprehensive guide on how to use DARTSGPT effectively

#### Screenshots

The Streamlit interface provides an intuitive way to interact with DARTSGPT:

1. **Main Interface** - Clean design with gradient background, sidebar controls, and input area
2. **Templates Browser** - Browse all available DARTS templates organized by physics type
3. **Documentation Tab** - Built-in guide explaining how to use DARTSGPT effectively

The interface features:
- Sidebar with AI model selection (GPT-4/GPT-3.5)
- Creativity level slider for controlling output variation
- Example prompts dropdown for quick testing
- Real-time statistics showing total generations and available templates
- About section with project information

### Command Line Interface

```bash
# Generate a DARTS model
uv run python main.py generate "Create a CO2 injection model for carbon storage"

# Generate with custom output directory
uv run python main.py generate "Build a waterflood model" -o my_model

# Preview without saving (dry run)
uv run python main.py generate "Model a geothermal reservoir" --dry-run

# Verbose output showing agent decisions
uv run python main.py generate "Create a black oil model" --verbose

# Generate AND execute the model (requires Docker on macOS)
uv run python main.py generate "Create a simple dead oil model" --execute

# List available templates
uv run python main.py list-templates

# Show example prompts
uv run python main.py show-examples --physics geothermal

# Check system configuration
uv run python main.py check-system

# Run interactive demo
uv run python main.py demo --interactive
```

### Python API

```python
from orchestrator import run_dartsgpt_orchestrator

# Generate a model programmatically
result = run_dartsgpt_orchestrator(
    prompt="Create a geothermal reservoir model at 200°C with injection and production wells",
    output_dir="geothermal_model"
)

# Check results
if result['success']:
    print(f"Model generated successfully!")
    print(f"Template used: {result['template_used']}")
    print(f"Files saved to: geothermal_model/")

# Generate and execute model
result = run_dartsgpt_orchestrator(
    prompt="Create a simple waterflood model",
    output_dir="waterflood_model",
    execute=True  # This will run the model after generation
)

if result['success'] and result.get('execution'):
    print(f"Execution completed: {result['execution']['success']}")
    print(f"Simulation time: {result['execution']['execution_time']:.2f}s")
```

## 🚀 Model Execution

DARTSGPT can automatically execute the generated models, providing immediate validation and results.

### Execution Features

- **Automatic execution** after generation with `--execute` flag
- **Real-time output streaming** during simulation
- **Resource monitoring** (memory, CPU usage)
- **Execution report** with convergence status and metrics

### Platform Compatibility

⚠️ **Important Note**: DARTS (open-darts) only runs on Linux x86_64. 
- **Linux**: Full functionality including model execution
- **macOS/Windows**: Model generation and validation only
- **Solution**: Use a Linux machine or VM for model execution

### Using Model Execution (Linux Only)

#### Command Line
```bash
# Generate and execute a model
uv run python main.py generate "Create a CO2 injection model" --execute

# Execute an existing model
uv run python main.py execute output/model_folder/
```

#### Streamlit Interface
Simply check the "Execute generated model" checkbox before clicking "Generate DARTS Model".

#### Execution Output
The execution agent provides:
- Simulation progress and convergence status
- Time step information
- Warning and error messages
- Final simulation metrics
- Location of result files

### Docker Setup for Cross-Platform Development

Since DARTS requires Linux libraries, macOS users must use Docker:

```bash
# Build and run with Docker Compose
docker-compose up -d

# Access Streamlit UI
open http://localhost:8502

# Run CLI commands in Docker
docker-compose run dartsgpt-cli uv run python main.py generate "Your prompt" --execute
```

## 📚 Usage Examples

### Example 1: CO2 Storage

```python
prompt = """
Create a CO2 injection model for carbon capture and storage:
- Saline aquifer at 1000m depth
- 3D grid: 100x100x20 cells, each 50m x 50m x 5m
- Porosity: 0.2, Permeability: 100 mD
- One injection well at center, rate: 1000 m3/day
- Include caprock with low permeability
- Simulate for 30 years
"""

result = await generator.generate(prompt)
```

### Example 2: Geothermal Energy

```python
prompt = """
Model a geothermal doublet system:
- High-temperature reservoir at 3000m depth
- Initial temperature: 180°C
- Production rate: 150 kg/s
- Reinjection temperature: 70°C
- Include thermal breakthrough analysis
"""

result = await generator.generate(prompt)
```

### Example 3: Enhanced Oil Recovery

```python
prompt = """
Design a polymer flooding EOR model:
- Oil reservoir with 30% initial oil saturation
- Polymer injection for mobility control
- Five-spot well pattern
- Include adsorption and degradation
"""

result = await generator.generate(prompt)
```

## 📦 DARTS Code Examples

### Compositional Flow Model

Here's an example of a generated compositional model for CO2 injection:

```python
from darts.reservoirs.struct_reservoir import StructReservoir
from darts.models.cicd_model import CICDModel
from darts.physics.super.physics import Compositional
from darts.physics.super.property_container import PropertyContainer
from darts.physics.properties.flash import ConstantK
from darts.physics.properties.basic import ConstFunc, PhaseRelPerm
from darts.physics.properties.density import DensityBasic
import numpy as np

class Model(CICDModel):
    def __init__(self):
        super().__init__()
        self.timer.node["initialization"].start()
        
        self.set_reservoir()
        self.set_physics()
        self.set_sim_params(first_ts=0.001, mult_ts=2, max_ts=1, 
                           runtime=1000, tol_newton=1e-2, tol_linear=1e-3,
                           it_newton=10, it_linear=50)
        
        self.timer.node["initialization"].stop()

    def set_reservoir(self):
        # Create 1D reservoir for testing
        nx = 1000
        self.reservoir = StructReservoir(self.timer, nx=nx, ny=1, nz=1, 
                                       dx=1, dy=10, dz=10,
                                       permx=100, permy=100, permz=10, 
                                       poro=0.3, depth=1000)

    def set_physics(self):
        # Define components and phases
        components = ['CO2', 'C1', 'H2O']
        phases = ['gas', 'oil']
        Mw = [44.01, 16.04, 18.015]  # Molecular weights
        
        property_container = PropertyContainer(phases_name=phases, 
                                             components_name=components,
                                             Mw=Mw, min_z=1e-9, temperature=1.)
        
        # Define property correlations
        property_container.flash_ev = ConstantK(len(components), [4, 2, 1e-1], 1e-8)
        property_container.density_ev = dict([
            ('gas', DensityBasic(compr=1e-3, dens0=200)),
            ('oil', DensityBasic(compr=1e-5, dens0=600))
        ])
        property_container.viscosity_ev = dict([
            ('gas', ConstFunc(0.05)),
            ('oil', ConstFunc(0.5))
        ])
        property_container.rel_perm_ev = dict([
            ('gas', PhaseRelPerm("gas")),
            ('oil', PhaseRelPerm("oil"))
        ])
        
        # Create physics object
        self.physics = Compositional(components, phases, self.timer,
                                   state_spec=Compositional.StateSpecification.P,
                                   n_points=200, min_p=1, max_p=300)
        self.physics.add_property_region(property_container)
```

### Geothermal Model with Temperature

Here's an example of a geothermal model with temperature tracking:

```python
from darts.physics.properties.enthalpy import EnthalpyBasic

class GeothermalModel(CICDModel):
    def set_physics(self):
        # Define water component for geothermal
        components = ['w']
        phases = ['wat', 'gas']
        
        property_container = PropertyContainer(phases_name=phases, 
                                             components_name=components)
        
        # Define density and viscosity for water and steam
        property_container.density_ev = dict([
            ('wat', DensityBasic(compr=1e-5, dens0=1014)),
            ('gas', DensityBasic(compr=5e-3, dens0=50))
        ])
        property_container.viscosity_ev = dict([
            ('wat', ConstFunc(0.3)),
            ('gas', ConstFunc(0.03))
        ])
        
        # Add enthalpy for temperature tracking
        property_container.enthalpy_ev = dict([
            ('wat', EnthalpyBasic(hcap=4.18)),  # Water heat capacity
            ('gas', EnthalpyBasic(hcap=0.035))  # Steam heat capacity
        ])
        property_container.rock_energy_ev = EnthalpyBasic(hcap=1.0)
        
        # Enable thermal mode
        thermal = True
        state_spec = Compositional.StateSpecification.PT
        self.physics = Compositional(components, phases, self.timer,
                                   state_spec=state_spec, n_points=400,
                                   min_p=0, max_p=1000, 
                                   min_t=293.15, max_t=473.15)  # 20-200°C
        self.physics.add_property_region(property_container)

    def set_initial_conditions(self):
        # Set pressure and temperature
        input_distribution = {
            self.physics.vars[0]: 200.,    # Pressure [bar]
            self.physics.vars[1]: 350.,     # Temperature [K] = 77°C
        }
        return self.physics.set_initial_conditions_from_array(
            mesh=self.reservoir.mesh,
            input_distribution=input_distribution
        )
```

### CO2 Storage with Advanced EoS

Here's an example using advanced equations of state for CO2 storage:

```python
from dartsflash.libflash import NegativeFlash, CubicEoS, AQEoS, FlashParams
from dartsflash.components import CompData
from darts.physics.properties.eos_properties import EoSDensity, EoSEnthalpy
from darts.physics.properties.density import Garcia2001
from darts.physics.properties.viscosity import Fenghour1998, Islam2012

class CCSModel(DartsModel):
    def set_physics(self, zero=1e-13, n_points=400, temperature=None):
        # Components for CO2-water system
        components = ["H2O", "CO2"]
        phases = ["Aq", "V"]  # Aqueous and Vapor
        
        comp_data = CompData(components, setprops=True)
        
        # Set up equations of state
        pr = CubicEoS(comp_data, CubicEoS.PR)  # Peng-Robinson for gas
        aq = AQEoS(comp_data, AQEoS.Ziabakhsh2012)  # Aqueous EoS
        
        flash_params = FlashParams(comp_data)
        flash_params.add_eos("PR", pr)
        flash_params.add_eos("AQ", aq)
        flash_params.eos_order = ["AQ", "PR"]
        
        # Property container with advanced correlations
        property_container = PropertyContainer(phases_name=phases, 
                                             components_name=components,
                                             Mw=comp_data.Mw)
        
        # Use negative flash for phase equilibrium
        property_container.flash_ev = NegativeFlash(flash_params, ["AQ", "PR"])
        
        # Advanced density correlations
        property_container.density_ev = dict([
            ('V', EoSDensity(pr, comp_data.Mw)),
            ('Aq', Garcia2001(components))
        ])
        
        # Advanced viscosity correlations
        property_container.viscosity_ev = dict([
            ('V', Fenghour1998()),
            ('Aq', Islam2012(components))
        ])
        
        # Output properties for monitoring
        property_container.output_props = {
            "satA": lambda: property_container.sat[0],
            "satV": lambda: property_container.sat[1],
            "xCO2": lambda: property_container.x[0, 1],  # CO2 in water
            "yH2O": lambda: property_container.x[1, 0]   # H2O in gas
        }
        
        # Create physics with thermal option
        thermal = temperature is None
        state_spec = Compositional.StateSpecification.PT if thermal else Compositional.StateSpecification.P
        
        self.physics = Compositional(components, phases, self.timer, n_points,
                                   min_p=1, max_p=400, min_z=zero/10, max_z=1-zero/10,
                                   min_t=273.15, max_t=373.15, state_spec=state_spec)
        self.physics.add_property_region(property_container)
```

### Well Control Examples

```python
from darts.engines import well_control_iface

def set_well_controls(self):
    for i, w in enumerate(self.reservoir.wells):
        if 'I' in w.name:  # Injection well
            # BHP control with composition
            self.physics.set_well_controls(
                wctrl=w.control, 
                control_type=well_control_iface.BHP,
                is_inj=True, 
                target=140.,  # Target pressure [bar]
                inj_composition=[0.99, 0.01],  # CO2, H2O
                inj_temp=350.  # Temperature [K] for thermal
            )
        else:  # Production well
            # Rate control
            self.physics.set_well_controls(
                wctrl=w.control,
                control_type=well_control_iface.MOLAR_RATE,
                is_inj=False, 
                target=5.,  # Production rate [mol/s]
                phase_name='wat'
            )
```

These examples show the key components of DARTS models:
- **Reservoir Definition**: Grid structure, properties (porosity, permeability)
- **Physics Setup**: Components, phases, property correlations
- **Initial Conditions**: Pressure, temperature, composition distributions
- **Well Controls**: Injection/production with various control types

## 🔬 Physics Types

### Compositional Flow
Best for: CO2 storage, gas injection, miscible displacement

```python
# Example prompt
"Model CO2 injection into a depleted gas reservoir with CH4 and C2H6"
```

### Dead Oil
Best for: Water flooding, immiscible displacement, thermal recovery

```python
# Example prompt
"Create a water flooding model for heavy oil recovery"
```

### Black Oil
Best for: Conventional oil reservoirs with dissolved gas

```python
# Example prompt
"Model a three-phase reservoir with solution gas drive"
```

### Geothermal
Best for: Geothermal energy, heat extraction, thermal systems

```python
# Example prompt
"Design a binary cycle geothermal plant with reinjection"
```

### Poroelastic
Best for: Coupled flow-mechanics, subsidence, stress analysis

```python
# Example prompt
"Create a poroelastic model for studying reservoir compaction"
```

### Chemical/EOR
Best for: Surfactant flooding, polymer injection, reactive transport

```python
# Example prompt
"Build a chemical EOR model with ASP flooding"
```

## 📂 Templates & Examples

### Template System

DARTSGPT includes 28+ pre-validated templates across all physics types:

- **Golden Templates** (⭐): Verified and recommended for production use
- **Basic Templates**: Simple models for learning and prototyping
- **Advanced Templates**: Complex physics and specialized applications

### Example Directory Structure

```
examples/
├── compositional/     # Multi-component flow examples
├── dead_oil/         # Two-phase flow examples
├── geothermal/       # Thermal reservoir examples
├── ccs/              # CO2 storage examples
├── poroelastic/      # Geomechanics examples
├── chemical/         # Reactive transport examples
├── example_catalog.md    # Comprehensive example documentation
└── browse_examples.py    # Interactive example browser
```

### Using Examples

1. **Browse interactively**: `python examples/browse_examples.py`
2. **View catalog**: Open `examples/example_catalog.md`
3. **Explore specific physics**: Check individual README files in subdirectories

## 🤖 Agent System

### Intent Classification
The Intent Classifier uses advanced NLP to understand:
- Physics type requirements
- Grid specifications
- Well configurations
- Special features needed

### Template Selection
Intelligent matching based on:
- Physics compatibility
- Feature requirements
- Complexity level
- Best practices

### Parameter Extraction
Handles:
- Unit conversion (ft to m, psi to Pa, etc.)
- Default values for missing parameters
- Range validation
- Consistency checks

### Code Generation
Produces:
- Complete `model.py` with proper inheritance
- Configured `main.py` with run parameters
- Proper imports and dependencies
- Documentation comments explaining choices

### Enhanced Code Generation Features

The latest version includes:

1. **Header docstrings** showing:
   - Original user prompt
   - Model configuration details
   - Template selection rationale
   
2. **Method-level comments** explaining:
   - Why specific physics were chosen
   - Parameter interpretation
   - Configuration choices
   
3. **Inline annotations** providing:
   - Unit conversions with original values
   - Suggestions for post-processing
   - Context for each decision

## 📚 Knowledge Base & Template Management

### How DARTSGPT Learns About DARTS

DARTSGPT uses an advanced RAG (Retrieval-Augmented Generation) system powered by vector embeddings to understand DARTS code patterns and best practices.

#### Vector Store Architecture

The system uses **Chroma** as a vector database with **OpenAI embeddings** to create a searchable knowledge base:

```
Knowledge Base Components:
├── Vector Store (Chroma DB)
│   ├── DARTS Model Embeddings
│   ├── Template Embeddings
│   ├── Documentation Embeddings
│   └── Test Suite Embeddings
├── Model Analyzer
│   └── Extracts metadata from DARTS models
├── Template Database
│   └── Pre-defined model templates
└── RAG System
    └── Semantic search interface
```

#### Data Sources Ingested

1. **DARTS Models** (`knowledge/darts-code/open-darts-main/models/`):
   - Complete model.py files from all example models
   - Main.py runner scripts
   - Extracted metadata: physics type, grid type, features, complexity

2. **Template Database** (`knowledge/templates/`):
   - 28+ pre-validated templates across all physics types
   - Template metadata including features, parameters, use cases
   - Golden templates marked for best practices

3. **Documentation**:
   - README files from DARTS repository
   - Physics module documentation
   - Tutorial content

4. **Test Suite** (`Run Test Suite 2.py`):
   - Test patterns and validation logic
   - Performance benchmarks
   - Model initialization patterns

#### Knowledge Base Initialization

```python
# How the knowledge base is built (from rag_system.py)
class DARTSKnowledgeBase:
    def __init__(self):
        self.embeddings = OpenAIEmbeddings()
        self.vector_store = Chroma(
            persist_directory="./embeddings/darts_knowledge",
            embedding_function=self.embeddings
        )
    
    def initialize(self):
        # 1. Analyze and index all DARTS models
        model_docs = self._index_models()  # Extracts AST, physics, features
        
        # 2. Index template information
        template_docs = self._index_templates()
        
        # 3. Index documentation
        doc_docs = self._index_documentation()
        
        # 4. Create vector embeddings
        self.vector_store = Chroma.from_documents(
            documents=all_docs,
            embedding=self.embeddings
        )
```

#### Model Analysis Process

The system uses AST (Abstract Syntax Tree) analysis to understand DARTS models:

```python
# From model_analyzer.py
class DARTSModelAnalyzer:
    def analyze_model(self, model_path):
        # Parse Python AST
        tree = ast.parse(content)
        
        # Extract:
        # - Physics type from imports and class usage
        # - Grid type (structured, unstructured, CPG)
        # - Features (thermal, chemical, mechanical)
        # - Well configurations
        # - Complexity assessment
        
        return ModelMetadata(
            physics_type="compositional",
            features=["thermal", "reactions"],
            grid_type="structured",
            complexity="advanced"
        )
```

#### How Agents Query the Knowledge Base

Agents use semantic search to find relevant examples and patterns:

```python
# Template Selector Agent using RAG
def find_similar_examples(physics_type, features):
    # Construct semantic query
    query = f"physics:{physics_type} {' '.join(features)}"
    
    # Search vector store
    results = self.rag_system.search(query, k=5)
    
    # Results include:
    # - Similar model code
    # - Relevance scores
    # - Metadata (model name, physics type, features)
    return results
```

### Template System

Templates are stored in a structured database with comprehensive metadata:

```python
# From template_database.py
class TemplateInfo(BaseModel):
    name: str                    # Template identifier
    path: str                    # Path to template files
    physics_type: str            # Physics classification
    description: str             # Human-readable description
    features: List[str]          # Supported features
    complexity: str              # simple/moderate/complex
    is_golden: bool              # Best practice template
    required_parameters: List    # Must-have parameters
    optional_parameters: List    # Optional parameters
    example_use_cases: List      # Example applications
```

### Building and Updating the Knowledge Base

#### Initial Build
```bash
# Build vector store from scratch
python scripts/test_knowledge_base.py

# This will:
# 1. Scan all DARTS models
# 2. Extract metadata
# 3. Create embeddings
# 4. Store in Chroma DB
```

#### Adding New Knowledge

To add new DARTS models or templates to the knowledge base:

1. **Add Model Files**: Place new models in the DARTS models directory
2. **Update Templates**: Add new templates to `knowledge/templates/`
3. **Rebuild Index**: Run with `force_rebuild=True`

```python
kb = DARTSKnowledgeBase()
kb.initialize(force_rebuild=True)  # Rebuilds from scratch
```

### Search Capabilities

The knowledge base supports various search types:

```python
# Search by query
results = kb.search("CO2 injection in saline aquifer")

# Search by physics type
models = kb.search_models_by_physics("compositional")

# Search for templates
templates = kb.search_templates("water injection")

# Get parameter examples
examples = kb.get_parameter_examples("2ph_comp")
```

Each search returns:
- Relevant code snippets
- Metadata (physics type, features, complexity)
- Similarity scores
- Source file information

## 🤖 Multi-Agent System Architecture

DARTSGPT employs a sophisticated multi-agent architecture where each agent has specialized responsibilities:

### Agent Overview

The system uses LangGraph for orchestration and consists of five specialized agents that work together in a pipeline:

```
User Prompt → Intent Classifier → Template Selector → Parameter Extractor → Code Generator → Validator → Output
```

### Detailed Agent Descriptions

#### 1. Intent Classifier Agent

**Purpose**: Analyzes natural language prompts to understand user requirements

**Key Responsibilities**:
- Identifies physics type (compositional, dead oil, black oil, geothermal, etc.)
- Extracts key features (thermal, reactions, tracers, etc.)
- Determines model complexity (simple, moderate, complex)
- Extracts grid specifications if provided
- Identifies well configurations

**Tools Used**:
- `analyze_physics_keywords`: Scans for physics-related keywords
- `analyze_feature_keywords`: Identifies required features
- `extract_grid_specifications`: Parses grid dimensions (e.g., "50x50x10")
- `extract_well_configuration`: Identifies well patterns and types

**How It Works**:
```python
# Example processing
Input: "Create a CO2 injection model in a 100x100x20 grid"
Output: {
    "physics_type": "compositional",
    "features": ["co2_injection", "two_phase"],
    "complexity": "moderate",
    "grid_info": {"nx": 100, "ny": 100, "nz": 20},
    "confidence": 0.95
}
```

#### 2. Template Selector Agent

**Purpose**: Selects the most appropriate DARTS template based on classified intent

**Key Responsibilities**:
- Scores all available templates against requirements
- Considers physics compatibility (weight: 0.5)
- Evaluates feature support (weight: 0.3)
- Assesses complexity match (weight: 0.2)
- Prioritizes "golden" templates (best practices)
- Uses RAG system to find similar examples

**Tools Used**:
- `score_physics_compatibility`: Calculates physics match score
- `score_feature_compatibility`: Evaluates feature coverage
- `score_complexity_match`: Assesses complexity alignment
- `find_similar_examples`: Queries vector store for similar models

**Selection Process**:
```python
# Scoring algorithm
overall_score = (
    physics_score * 0.5 +
    feature_score * 0.3 +
    complexity_score * 0.2
)
if template.is_golden:
    overall_score *= 1.1  # 10% bonus for golden templates
```

#### 3. Parameter Extractor Agent

**Purpose**: Extracts numerical parameters from natural language and converts units

**Key Responsibilities**:
- Identifies numerical values with units
- Converts units to DARTS standard (SI units)
- Applies template-specific defaults
- Validates parameter ranges
- Handles ambiguous values intelligently

**Tools Used**:
- `extract_numerical_values`: Regex-based parameter extraction
- `convert_units`: Unit conversion to SI
- `apply_defaults`: Fills missing parameters
- `validate_parameters`: Range and consistency checks

**Unit Conversions**:
```python
# Common conversions
- Permeability: mD → m² (multiply by 1e-15)
- Temperature: °F → K (fahrenheit_to_kelvin)
- Pressure: psi → Pa (multiply by 6894.76)
- Depth: ft → m (multiply by 0.3048)
- Rate: bbl/day → m³/s (multiply by 1.84e-6)
```

**Example Processing**:
```python
Input: "100 mD permeability at 1000 feet depth"
Output: {
    "rock_properties": {
        "permeability": 1e-13,  # m²
        "depth": 304.8          # m
    },
    "units_used": {
        "permeability": "m²",
        "depth": "m"
    }
}
```

#### 4. Code Generator Agent

**Purpose**: Generates complete DARTS model code from templates and parameters

**Key Responsibilities**:
- Loads selected template code
- Fills template with extracted parameters
- Generates physics-specific setup
- Creates well configurations
- Produces both model.py and main.py
- Adds explanatory comments

**Tools Used**:
- `load_template_code`: Retrieves template files
- `generate_physics_setup`: Creates physics initialization
- `generate_reservoir_setup`: Configures grid and properties
- `generate_well_setup`: Creates well locations and controls
- `add_documentation`: Adds helpful comments

**Code Generation Features**:
- Header docstrings with original prompt
- Method-level comments explaining choices
- Unit conversion annotations
- Parameter source tracking
- Best practice suggestions

#### 5. Validator Agent

**Purpose**: Ensures generated code is syntactically correct and physics-consistent

**Key Responsibilities**:
- Checks Python syntax validity
- Verifies all required methods exist
- Validates physics consistency
- Ensures import completeness
- Checks parameter ranges
- Provides fix suggestions

**Tools Used**:
- `check_python_syntax`: AST-based syntax validation
- `check_required_methods`: Ensures DARTS interface compliance
- `check_physics_consistency`: Validates physics setup
- `check_imports`: Verifies all imports present
- `suggest_fixes`: Provides correction recommendations

**Validation Checks**:
```python
# Required methods for DARTS models
required_methods = [
    "set_physics",
    "set_reservoir",
    "set_wells",
    "set_initial_conditions",
    "set_well_controls"
]

# Physics-specific validations
- Compositional: Flash calculations, EOS parameters
- Dead Oil: Relative permeability, PVT data
- Geothermal: Temperature ranges, enthalpy
- Black Oil: Three-phase consistency
```

### Agent Communication

Agents communicate through a shared state managed by LangGraph:

```python
class DARTSGPTState(MessagesState):
    current_agent: str
    prompt: str
    intent: Dict[str, Any]
    selected_template: str
    parameters: Dict[str, Any]
    model_code: str
    main_code: str
    validation_result: Dict[str, Any]
    final_output: Dict[str, Any]
```

Each agent:
1. Receives the current state
2. Performs its specialized task
3. Updates the state with results
4. Passes control to the next agent

### Error Handling and Fallbacks

Each agent includes robust error handling:
- **Intent Classifier**: Falls back to keyword matching if LLM parsing fails
- **Template Selector**: Selects highest-scoring template if LLM fails
- **Parameter Extractor**: Uses template defaults for missing parameters
- **Code Generator**: Provides basic template if generation fails
- **Validator**: Attempts auto-fixes for common issues

## 🔧 Development

### Project Structure

```
dartsgpt/
├── orchestrator.py        # Multi-agent orchestrator with LangGraph
├── main.py               # CLI interface with rich output
├── config.py             # Configuration settings
├── knowledge/            # Knowledge base and templates
│   ├── templates/
│   │   └── template_database.py
│   ├── physics_features.py
│   └── example_prompts.py
├── examples/             # Example models for each physics type
├── tests/                # Comprehensive test suite
├── docs/                 # Documentation
│   ├── presentation/     # LaTeX presentation files
│   ├── DARTSGPT_Architecture_Diagrams.md
│   ├── DARTSGPT_Implementation_Plan.md
│   └── PROJECT_STATUS.md
└── pyproject.toml       # Project configuration
```

### Running Tests

```bash
# Run comprehensive test suite
uv run python -m pytest tests/

# Test specific physics type
uv run python -m pytest tests/test_physics.py::test_compositional

# Test with coverage
uv run python -m pytest --cov=dartsgpt tests/
```

## 🧪 Testing

### Test Coverage

- **Unit Tests**: All agents and tools individually tested
- **Integration Tests**: Multi-agent workflows validated
- **Edge Cases**: Complex scenarios and error conditions
- **Performance Tests**: Generation time and scalability

### Evaluation Framework

DARTSGPT includes dual evaluation methods:

1. **Vector Store Evaluation**: Uses RAG to find similar examples
2. **Ground Truth Evaluation**: Compares against known correct answers

```python
from dartsgpt.evaluation import TestSuite, AgentEvaluator

# Load and run test suite
test_suite = TestSuite.load_default_suite()
evaluator = AgentEvaluator()
results = evaluator.evaluate_test_suite(agent, test_suite, method="both")

print(f"Success rate: {results['summary']['success_rate']:.2%}")
print(f"Physics accuracy: {results['summary']['physics_accuracy']:.2%}")
```

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Areas for Contribution

- Additional physics templates
- Improved parameter extraction
- Enhanced validation rules
- Documentation improvements
- Test cases
- UI enhancements

### Development Setup

```bash
# Install dev dependencies
uv add --dev pytest black ruff mypy

# Run formatter
uv run black .

# Run linter
uv run ruff check .

# Run type checker
uv run mypy .
```

## 📊 Project Status

### Current State: Production Ready ✅

DARTSGPT is now a fully functional multi-agent system that generates DARTS reservoir simulation models from natural language descriptions.

### Completed Features
- ✅ Multi-Agent Orchestration with LangGraph
- ✅ Intent Classification with physics identification
- ✅ Template Selection with golden template priority
- ✅ Parameter Extraction with unit conversions
- ✅ Enhanced Code Generation with explanatory comments
- ✅ Comprehensive Validation
- ✅ Rich CLI Interface
- ✅ 28+ Templates across all physics types
- ✅ Extensive Documentation
- ✅ Comprehensive Test Suite

### System Capabilities
- **Template Coverage**: All physics types supported with 28+ templates
- **Physics Support**: All DARTS physics modules covered
- **Generation Features**: Complete model generation with explanatory comments
- **Multi-Agent System**: Specialized agents for each task

### Known Limitations
1. Physics misclassification for ambiguous prompts
2. Limited to basic well patterns
3. Primarily Cartesian grids
4. Complex parameter distributions not fully supported

### Future Enhancements
1. Improved physics classification with structured output
2. Support for more complex well patterns
3. Enhanced parameter extraction for heterogeneous properties
4. Additional specialized templates

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- TU Delft DARTS Team for the simulation framework
- OpenAI for language models
- LangChain community for the agent framework

## 📚 Additional Resources

- [DARTS Official Documentation](https://darts.citg.tudelft.nl/)
- [Architecture Diagrams](docs/DARTSGPT_Architecture_Diagrams.md)
- [Implementation Plan](docs/DARTSGPT_Implementation_Plan.md)
- [Model Categorization](docs/DARTS_Model_Categorization.md)
- [Project Status](docs/PROJECT_STATUS.md)
- [Test Results](docs/TEST_RESULTS.md)

## 🏗️ Project Structure

This document describes the organization of the DARTSGPT codebase.

### Directory Structure

```
dartsgpt/
├── app.py                     # Streamlit web application
├── main.py                    # CLI entry point
├── orchestrator.py            # Multi-agent orchestrator
├── cli.py                     # CLI implementation
├── config.py                  # Configuration management
├── setup.py                   # Package setup
│
├── agents/                    # Agent implementations
│   ├── __init__.py
│   ├── intent_classifier.py   # Intent classification agent
│   ├── template_selector.py   # Template selection agent
│   ├── parameter_extractor.py # Parameter extraction agent
│   ├── code_generator.py      # Code generation agent
│   └── validator.py           # Code validation agent
│
├── knowledge/                 # Knowledge base and DARTS code
│   ├── __init__.py
│   ├── darts-code/           # DARTS source code and models
│   │   ├── open-darts-main/
│   │   └── darts-models-development/
│   ├── templates/            # Code generation templates
│   ├── physics_features.py   # Physics type definitions
│   ├── example_prompts.py    # Example prompts
│   ├── rag_system.py        # RAG implementation
│   └── model_analyzer.py     # DARTS model analysis
│
├── tests/                    # Test suite
│   ├── __init__.py
│   ├── run_all_tests.py     # Test runner
│   ├── unit/                # Unit tests
│   │   ├── __init__.py
│   │   ├── test_cli.py
│   │   ├── test_system.py
│   │   └── test_unit_tests.py
│   └── integration/         # Integration tests
│       ├── __init__.py
│       ├── test_comprehensive.py
│       ├── test_edge_cases.py
│       ├── test_full_agent.py
│       ├── test_integration.py
│       └── test_knowledge_base.py
│
├── demos/                   # Demo scripts
│   ├── __init__.py
│   ├── demo_dartsgpt.py    # Main demo
│   ├── demo_simple.py      # Simple demo
│   ├── run_example.py      # Example runner
│   └── run_cli.py          # CLI runner
│
├── scripts/                 # Utility scripts
│   └── utils/              # Utility functions
│
├── docs/                   # Documentation
│   ├── presentation/       # LaTeX presentation
│   │   └── dartsgpt_presentation.tex
│   └── architecture/       # Architecture diagrams
│
├── examples/               # Example outputs
│   └── generated_models/   # Generated DARTS models
│
├── embeddings/            # Vector store data
│   └── chroma/           # Chroma DB files
│
├── output/               # Generated model outputs
├── logs/                 # Application logs
│
├── .env                  # Environment variables (not in git)
├── .env.example          # Example environment file
├── pyproject.toml        # Project dependencies
├── uv.lock              # Locked dependencies
├── Dockerfile           # Docker image definition
├── docker-compose.yml   # Docker compose configuration
├── README.md            # Main documentation
└── LICENSE              # MIT License
```

### Key Files

#### Core Application
- `app.py` - Streamlit web interface for interactive model generation
- `main.py` - Command-line interface entry point
- `orchestrator.py` - LangGraph multi-agent orchestration system
- `config.py` - Centralized configuration management

#### Agent System
- `agents/` - Five specialized agents for the generation pipeline:
  - Intent Classifier - Understands user requirements
  - Template Selector - Chooses appropriate DARTS template
  - Parameter Extractor - Extracts and converts parameters
  - Code Generator - Generates DARTS model code
  - Validator - Validates generated code

#### Knowledge Base
- `knowledge/darts-code/` - Complete DARTS source code and examples
- `knowledge/templates/` - Pre-validated DARTS model templates
- `knowledge/rag_system.py` - Chroma vector store for semantic search
- `knowledge/model_analyzer.py` - AST-based model analysis

#### Testing
- `tests/unit/` - Component-level tests
- `tests/integration/` - End-to-end system tests
- `tests/run_all_tests.py` - Test suite runner

#### Documentation
- `README.md` - Project overview and usage
- `docs/presentation/` - LaTeX presentation materials

### Environment Setup

1. Copy `.env.example` to `.env`
2. Add your OpenAI API key
3. Install dependencies: `uv sync`
4. Run tests: `uv run python tests/run_all_tests.py`
5. Start app: `uv run streamlit run app.py`

### Development Workflow

1. **Adding Templates**: Place new templates in `knowledge/templates/`
2. **Running Tests**: Use `uv run python tests/run_all_tests.py`
3. **Debugging**: Set `DEBUG_MODE=true` in `.env`
4. **Docker**: Build with `docker-compose up -d`

### Output Locations

- **Generated Models**: `output/` directory
- **Vector Store**: `embeddings/` directory
- **Logs**: `logs/` directory
- **Test Results**: `tests/test_output/`

## 🐳 Docker Deployment

This guide explains how to run DARTSGPT using Docker.

### Prerequisites

- Docker 20.10+ installed
- Docker Compose 2.0+ installed
- At least 4GB of RAM available
- OpenAI API key

### Quick Start

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd dartsgpt
   ```

2. **Set up environment variables**:
   ```bash
   cp .env.example .env
   # Edit .env and add your OPENAI_API_KEY
   ```

3. **Build and run with Docker Compose**:
   ```bash
   docker-compose up -d
   ```

4. **Access the application**:
   - Open http://localhost:8502 in your browser
   - The Streamlit UI should be available

### Docker Commands

#### Build the image
```bash
docker build -t dartsgpt .
```

#### Run the Streamlit app
```bash
docker run -p 8502:8502 \
  -e OPENAI_API_KEY=your-key-here \
  -v $(pwd)/output:/app/output \
  -v $(pwd)/embeddings:/app/embeddings \
  dartsgpt
```

#### Run CLI commands
```bash
# Generate a model
docker run --rm \
  -e OPENAI_API_KEY=your-key-here \
  -v $(pwd)/output:/app/output \
  dartsgpt \
  uv run python main.py generate "Create a CO2 injection model"

# List templates
docker run --rm \
  -e OPENAI_API_KEY=your-key-here \
  dartsgpt \
  uv run python main.py list-templates

# Check system
docker run --rm \
  -e OPENAI_API_KEY=your-key-here \
  dartsgpt \
  uv run python main.py check-system
```

### Docker Compose Usage

#### Start services
```bash
# Start the Streamlit app
docker-compose up -d

# Start with rebuild
docker-compose up -d --build

# View logs
docker-compose logs -f dartsgpt
```

#### Run CLI commands via docker-compose
```bash
# Run a one-off CLI command
docker-compose run --rm dartsgpt-cli uv run python main.py generate "Create a geothermal model"

# Or use the CLI profile
docker-compose --profile cli run --rm dartsgpt-cli uv run python main.py list-templates
```

#### Stop services
```bash
docker-compose down
```

#### Clean up
```bash
# Remove containers and networks
docker-compose down

# Remove everything including volumes
docker-compose down -v
```

### Volume Mounts

The Docker setup uses several volumes:

- `./embeddings:/app/embeddings` - Vector store persistence
- `./output:/app/output` - Generated models
- `./logs:/app/logs` - Application logs

### Environment Variables

Key environment variables (set in `.env` file):

- `OPENAI_API_KEY` - Your OpenAI API key (required)
- `OPENAI_MODEL` - Model to use (default: gpt-4o-mini)
- `LOG_LEVEL` - Logging level (default: INFO)
- `DEBUG_MODE` - Enable debug mode (default: false)

### Resource Requirements

- **Minimum**: 2GB RAM, 1 CPU core
- **Recommended**: 4GB RAM, 2 CPU cores
- **Disk space**: ~2GB for Docker image + space for embeddings/outputs

### Troubleshooting

#### Container won't start
Check logs:
```bash
docker-compose logs dartsgpt
```

#### Permission issues
The container runs as non-root user (UID 1000). Ensure output directories are writable:
```bash
chmod -R 777 embeddings output logs
```

#### Out of memory
Increase Docker memory allocation or adjust resource limits in `docker-compose.yml`

#### API key not working
Ensure the API key is properly set:
```bash
docker-compose exec dartsgpt env | grep OPENAI
```

### Production Deployment

For production use:

1. Use environment-specific `.env` files
2. Set up proper volume backups
3. Configure reverse proxy (nginx/traefik)
4. Enable HTTPS
5. Set resource limits appropriately
6. Use Docker secrets for API keys

Example production compose override:
```yaml
# docker-compose.prod.yml
version: '3.8'

services:
  dartsgpt:
    restart: always
    deploy:
      resources:
        limits:
          cpus: '4'
          memory: 8G
    environment:
      - LOG_LEVEL=WARNING
      - DEBUG_MODE=false
```

Run with:
```bash
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d
```

---

<div align="center">
  Made with ❤️ by the open-DARTS team
</div>