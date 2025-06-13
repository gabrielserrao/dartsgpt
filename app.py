"""
DARTSGPT Streamlit Application
=============================

An interactive web interface for the DARTSGPT multi-agent system.
Shows real-time agent orchestration and code generation.
"""

import streamlit as st
import time
from pathlib import Path
import json
from datetime import datetime
import sys
import os
from typing import Dict, Any, List
import subprocess

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from orchestrator import build_dartsgpt_graph, DARTSGPTState
from langchain_core.messages import HumanMessage
from knowledge.templates.template_database import TEMPLATES
from knowledge.example_prompts import EXAMPLE_PROMPTS
from config import get_settings

# Configure Streamlit page
st.set_page_config(
    page_title="DARTSGPT - AI-Powered DARTS Model Generator",
    page_icon="🛢️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for beautiful UI
st.markdown("""
<style>
    /* Main container styling */
    .stApp {
        background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
    }
    
    /* Headers */
    h1 {
        color: #2c3e50;
        font-family: 'Arial Black', sans-serif;
        text-align: center;
        padding: 20px;
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }
    
    h2 {
        color: #34495e;
        border-bottom: 2px solid #667eea;
        padding-bottom: 10px;
        margin-top: 30px;
    }
    
    h3 {
        color: #7f8c8d;
    }
    
    /* Agent status cards */
    .agent-card {
        background: white;
        border-radius: 10px;
        padding: 15px;
        margin: 10px 0;
        box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        border-left: 4px solid #667eea;
        transition: all 0.3s ease;
    }
    
    .agent-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 5px 15px rgba(0,0,0,0.2);
    }
    
    .agent-active {
        border-left-color: #27ae60;
        background: #e8f8f5;
    }
    
    .agent-completed {
        border-left-color: #3498db;
        background: #ebf5fb;
    }
    
    .agent-error {
        border-left-color: #e74c3c;
        background: #fadbd8;
    }
    
    /* Metrics styling */
    .metric-card {
        background: white;
        border-radius: 15px;
        padding: 20px;
        text-align: center;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
    
    .metric-value {
        font-size: 2em;
        font-weight: bold;
        color: #667eea;
    }
    
    .metric-label {
        color: #7f8c8d;
        font-size: 0.9em;
        margin-top: 5px;
    }
    
    /* Code blocks */
    .stCode {
        background: #2c3e50 !important;
        border-radius: 10px;
        padding: 20px !important;
    }
    
    /* Buttons */
    .stButton > button {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        padding: 0.5rem 2rem;
        font-weight: bold;
        border-radius: 25px;
        transition: all 0.3s ease;
    }
    
    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 5px 15px rgba(102, 126, 234, 0.4);
    }
    
    /* Sidebar */
    .css-1d391kg {
        background: #2c3e50;
    }
    
    .sidebar .sidebar-content {
        background: #34495e;
    }
    
    /* Progress bars */
    .stProgress > div > div > div > div {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
    }
    
    /* Success/Error messages */
    .success-message {
        background: #d4edda;
        border: 1px solid #c3e6cb;
        color: #155724;
        padding: 15px;
        border-radius: 5px;
        margin: 10px 0;
    }
    
    .error-message {
        background: #f8d7da;
        border: 1px solid #f5c6cb;
        color: #721c24;
        padding: 15px;
        border-radius: 5px;
        margin: 10px 0;
    }
    
    /* Loading animation */
    @keyframes pulse {
        0% { opacity: 1; }
        50% { opacity: 0.5; }
        100% { opacity: 1; }
    }
    
    .loading {
        animation: pulse 1.5s ease-in-out infinite;
    }
</style>
""", unsafe_allow_html=True)

# Initialize session state
if 'agents_status' not in st.session_state:
    st.session_state.agents_status = {
        'supervisor': {'status': 'waiting', 'message': 'Ready to orchestrate'},
        'intent_classifier': {'status': 'waiting', 'message': 'Ready to analyze intent'},
        'template_selector': {'status': 'waiting', 'message': 'Ready to select template'},
        'parameter_extractor': {'status': 'waiting', 'message': 'Ready to extract parameters'},
        'code_generator': {'status': 'waiting', 'message': 'Ready to generate code'},
        'validator': {'status': 'waiting', 'message': 'Ready to validate'},
        'executor': {'status': 'waiting', 'message': 'Ready to execute'}
    }

if 'generation_history' not in st.session_state:
    st.session_state.generation_history = []

if 'current_output' not in st.session_state:
    st.session_state.current_output = None

# Sidebar
with st.sidebar:
    st.markdown("## 🛢️ DARTSGPT Control Panel")
    st.markdown("---")
    
    # Model selection
    model_choice = st.selectbox(
        "🤖 AI Model",
        ["gpt-4", "gpt-3.5-turbo"],
        help="Choose the AI model for code generation"
    )
    
    # Temperature slider
    temperature = st.slider(
        "🌡️ Creativity Level",
        min_value=0.0,
        max_value=1.0,
        value=0.7,
        step=0.1,
        help="Higher values make output more creative"
    )
    
    st.markdown("---")
    
    # Example prompts
    st.markdown("### 📝 Example Prompts")
    examples = [
        "CO2 injection for carbon storage",
        "Geothermal reservoir simulation",
        "Water flooding in oil reservoir",
        "Steam injection for heavy oil",
        "Polymer flooding for EOR"
    ]
    
    selected_example = st.selectbox(
        "Choose an example:",
        [""] + examples
    )
    
    st.markdown("---")
    
    # Statistics
    st.markdown("### 📊 Statistics")
    st.metric("Total Generations", len(st.session_state.generation_history))
    st.metric("Available Templates", len(TEMPLATES))
    st.metric("Physics Types", 6)
    
    # About section
    with st.expander("ℹ️ About DARTSGPT"):
        st.markdown("""
        **DARTSGPT** is an AI-powered system that generates 
        DARTS reservoir simulation models from natural language 
        descriptions.
        
        It uses a multi-agent architecture with specialized 
        agents for:
        - Intent understanding
        - Template selection
        - Parameter extraction
        - Code generation
        - Validation
        
        Created with ❤️ for reservoir engineers.
        """)

# Main content
st.title("🛢️ DARTSGPT - AI-Powered DARTS Model Generator")
st.markdown("### Transform natural language into reservoir simulation models")

# Create tabs
tab1, tab2, tab3, tab4 = st.tabs(["🚀 Generate Model", "📜 History", "📚 Templates", "📖 Documentation"])

with tab1:
    # Input section
    col1, col2 = st.columns([3, 1])
    
    with col1:
        prompt = st.text_area(
            "Describe your reservoir simulation model:",
            value=selected_example if selected_example else "",
            height=120,
            placeholder="Example: Create a CO2 injection model for carbon storage in a 100x100x20 reservoir with 25% porosity and 500 mD permeability",
            help="Describe the physics, geometry, properties, and objectives of your simulation"
        )
    
    with col2:
        st.markdown("### Quick Options")
        include_comments = st.checkbox("Include detailed comments", value=True)
        dry_run = st.checkbox("Preview only (don't save)", value=False)
        verbose = st.checkbox("Show detailed logs", value=True)
        execute_model = st.checkbox("Execute generated model", value=False, help="Run the model after generation (requires Docker)")
    
    # Function to render agent cards
    def render_agent_card(placeholder, agent_name: str, status: dict):
        """Render an agent status card."""
        status_class = {
            'waiting': '',
            'active': 'agent-active loading',
            'completed': 'agent-completed',
            'error': 'agent-error'
        }.get(status['status'], '')
        
        emoji = {
            'supervisor': '🎯',
            'intent_classifier': '🔍',
            'template_selector': '📋',
            'parameter_extractor': '🔢',
            'code_generator': '💻',
            'validator': '✅',
            'executor': '🚀'
        }.get(agent_name, '🤖')
        
        agent_display = agent_name.replace('_', ' ').title()
        
        placeholder.markdown(f"""
        <div class="agent-card {status_class}">
            <h4>{emoji} {agent_display}</h4>
            <p>{status['message']}</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Generate button
    if st.button("🎯 Generate DARTS Model", type="primary", use_container_width=True):
        if prompt:
            # Reset agent status
            for agent in st.session_state.agents_status:
                st.session_state.agents_status[agent] = {'status': 'waiting', 'message': 'Waiting...'}
            
            # Create containers for real-time updates
            status_container = st.container()
            progress_bar = st.progress(0)
            log_container = st.container()
            
            with status_container:
                st.markdown("## 🤖 Agent Orchestration")
                agent_cols = st.columns(3)
                
                # Create placeholders for each agent
                agent_placeholders = {}
                agents = list(st.session_state.agents_status.keys())
                
                for i, agent in enumerate(agents):
                    with agent_cols[i % 3]:
                        placeholder = st.empty()
                        agent_placeholders[agent] = placeholder
                        render_agent_card(placeholder, agent, st.session_state.agents_status[agent])
            
            # Log container
            with log_container:
                log_placeholder = st.empty()
                logs = []
                
                def add_log(message: str, level: str = "info"):
                    timestamp = datetime.now().strftime("%H:%M:%S")
                    emoji = {"info": "ℹ️", "success": "✅", "error": "❌", "warning": "⚠️"}.get(level, "📝")
                    logs.append(f"[{timestamp}] {emoji} {message}")
                    if verbose:
                        log_text = "\n".join(logs[-10:])  # Show last 10 logs
                        log_placeholder.code(log_text, language="text")
            
            try:
                # Initialize the graph
                add_log("Initializing DARTSGPT multi-agent system...", "info")
                progress_bar.progress(0.1)
                
                settings = get_settings()
                from langchain_openai import ChatOpenAI
                model = ChatOpenAI(
                    model=model_choice,
                    api_key=settings.openai_api_key,
                    temperature=temperature
                )
                
                graph = build_dartsgpt_graph(model)
                add_log(f"Graph initialized with {model_choice} model", "success")
                
                # Create initial state
                initial_state = {
                    "messages": [HumanMessage(content=prompt)],
                    "prompt": prompt,
                    "execute_model": execute_model and not dry_run,
                    "output_dir": "output"
                }
                
                # Simulate agent workflow with updates
                agents_sequence = [
                    ("supervisor", "Starting orchestration...", 0.2),
                    ("intent_classifier", "Analyzing your prompt...", 0.3),
                    ("template_selector", "Selecting best template...", 0.4),
                    ("parameter_extractor", "Extracting parameters...", 0.5),
                    ("code_generator", "Generating DARTS code...", 0.7),
                    ("validator", "Validating output...", 0.8)
                ]
                
                # Add executor if requested
                if execute_model and not dry_run:
                    agents_sequence.append(("executor", "Executing DARTS model...", 0.95))
                
                # Run through agents with visual updates
                for agent_name, message, progress in agents_sequence:
                    # Update agent status
                    st.session_state.agents_status[agent_name] = {
                        'status': 'active',
                        'message': message
                    }
                    render_agent_card(agent_placeholders[agent_name], agent_name, st.session_state.agents_status[agent_name])
                    
                    add_log(f"{agent_name}: {message}", "info")
                    progress_bar.progress(progress)
                    time.sleep(0.5)  # Visual delay
                    
                    # Mark as completed
                    st.session_state.agents_status[agent_name] = {
                        'status': 'completed',
                        'message': 'Completed successfully'
                    }
                    render_agent_card(agent_placeholders[agent_name], agent_name, st.session_state.agents_status[agent_name])
                
                # Run the actual graph
                add_log("Executing multi-agent workflow...", "info")
                result = graph.invoke(initial_state)
                
                # Extract final output
                final_output = result.get('final_output', {})
                
                if final_output.get('success'):
                    progress_bar.progress(1.0)
                    add_log("Model generation completed successfully!", "success")
                    
                    # Store in session state
                    st.session_state.current_output = final_output
                    st.session_state.generation_history.append({
                        'timestamp': datetime.now(),
                        'prompt': prompt,
                        'output': final_output
                    })
                    
                    # Display results
                    st.markdown("---")
                    st.success("✅ Model generated successfully!")
                    
                    # Show metrics
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.markdown("""
                        <div class="metric-card">
                            <div class="metric-value">✅</div>
                            <div class="metric-label">Status</div>
                        </div>
                        """, unsafe_allow_html=True)
                    with col2:
                        st.markdown(f"""
                        <div class="metric-card">
                            <div class="metric-value">{final_output.get('template_used', 'N/A')}</div>
                            <div class="metric-label">Template Used</div>
                        </div>
                        """, unsafe_allow_html=True)
                    with col3:
                        params = final_output.get('parameters', {})
                        param_count = sum(len(v) for v in params.values() if isinstance(v, dict))
                        st.markdown(f"""
                        <div class="metric-card">
                            <div class="metric-value">{param_count}</div>
                            <div class="metric-label">Parameters Extracted</div>
                        </div>
                        """, unsafe_allow_html=True)
                    with col4:
                        code_lines = len(final_output.get('model_code', '').splitlines())
                        st.markdown(f"""
                        <div class="metric-card">
                            <div class="metric-value">{code_lines}</div>
                            <div class="metric-label">Lines of Code</div>
                        </div>
                        """, unsafe_allow_html=True)
                    
                    # Display generated code
                    st.markdown("### 📄 Generated Files")
                    
                    # Add execution tab if model was executed
                    if execute_model and final_output.get('execution'):
                        tab_model, tab_main, tab_params, tab_exec = st.tabs(["model.py", "main.py", "parameters.json", "execution.log"])
                    else:
                        tab_model, tab_main, tab_params = st.tabs(["model.py", "main.py", "parameters.json"])
                    
                    with tab_model:
                        st.code(final_output.get('model_code', ''), language='python')
                        if not dry_run:
                            st.download_button(
                                "📥 Download model.py",
                                final_output.get('model_code', ''),
                                "model.py",
                                "text/x-python"
                            )
                    
                    with tab_main:
                        st.code(final_output.get('main_code', ''), language='python')
                        if not dry_run:
                            st.download_button(
                                "📥 Download main.py",
                                final_output.get('main_code', ''),
                                "main.py",
                                "text/x-python"
                            )
                    
                    with tab_params:
                        params_json = json.dumps(final_output.get('parameters', {}), indent=2)
                        st.code(params_json, language='json')
                        if not dry_run:
                            st.download_button(
                                "📥 Download parameters.json",
                                params_json,
                                "parameters.json",
                                "application/json"
                            )
                    
                    # Display execution results if available
                    if execute_model and final_output.get('execution'):
                        with tab_exec:
                            exec_result = final_output['execution']
                            
                            # Show execution status
                            if exec_result.get('success'):
                                st.success(f"✅ Model executed successfully in {exec_result.get('execution_time', 0):.2f} seconds")
                            else:
                                st.error("❌ Execution failed")
                            
                            # Show execution output
                            if exec_result.get('output'):
                                st.markdown("**Output:**")
                                st.code(exec_result['output'], language='text')
                            
                            # Show execution errors if any
                            if exec_result.get('error'):
                                st.markdown("**Errors:**")
                                st.code(exec_result['error'], language='text')
                            
                            # Show summary metrics
                            if exec_result.get('summary'):
                                summary = exec_result['summary']
                                st.markdown("**Execution Summary:**")
                                sum_col1, sum_col2 = st.columns(2)
                                with sum_col1:
                                    st.metric("Simulation Completed", "✅ Yes" if summary.get('simulation_completed') else "❌ No")
                                    st.metric("Convergence", "✅ Yes" if summary.get('convergence_achieved') else "❌ No")
                                with sum_col2:
                                    if summary.get('final_time'):
                                        st.metric("Final Time", f"{summary['final_time']} days")
                                    if summary.get('warnings'):
                                        st.metric("Warnings", len(summary['warnings']))
                    
                    # Save files if not dry run
                    if not dry_run:
                        output_dir = Path("streamlit_output") / datetime.now().strftime("%Y%m%d_%H%M%S")
                        output_dir.mkdir(parents=True, exist_ok=True)
                        
                        (output_dir / "model.py").write_text(final_output.get('model_code', ''))
                        (output_dir / "main.py").write_text(final_output.get('main_code', ''))
                        (output_dir / "metadata.json").write_text(json.dumps({
                            "prompt": prompt,
                            "template_used": final_output.get('template_used'),
                            "parameters": final_output.get('parameters'),
                            "timestamp": datetime.now().isoformat()
                        }, indent=2))
                        
                        st.info(f"📁 Files saved to: {output_dir}")
                
                else:
                    st.error("❌ Model generation failed. Please check the logs.")
                    add_log("Generation failed", "error")
                    
            except Exception as e:
                st.error(f"❌ Error: {str(e)}")
                add_log(f"Error: {str(e)}", "error")
                # Update agent status to show error
                for agent in st.session_state.agents_status:
                    if st.session_state.agents_status[agent]['status'] == 'active':
                        st.session_state.agents_status[agent] = {
                            'status': 'error',
                            'message': 'Error occurred'
                        }
                        render_agent_card(agent_placeholders[agent], agent, st.session_state.agents_status[agent])
        else:
            st.warning("⚠️ Please enter a prompt to generate a model.")

with tab2:
    st.markdown("## 📜 Generation History")
    
    if st.session_state.generation_history:
        for i, item in enumerate(reversed(st.session_state.generation_history)):
            with st.expander(f"Generation {len(st.session_state.generation_history) - i}: {item['timestamp'].strftime('%Y-%m-%d %H:%M:%S')}"):
                st.markdown(f"**Prompt:** {item['prompt']}")
                st.markdown(f"**Template:** {item['output'].get('template_used', 'N/A')}")
                
                col1, col2 = st.columns(2)
                with col1:
                    if st.button(f"View Code", key=f"view_{i}"):
                        st.code(item['output'].get('model_code', ''), language='python')
                with col2:
                    if st.button(f"Reuse Prompt", key=f"reuse_{i}"):
                        st.session_state.reuse_prompt = item['prompt']
                        st.rerun()
    else:
        st.info("No generation history yet. Generate your first model!")

with tab3:
    st.markdown("## 📚 Available Templates")
    
    # Group templates by physics type
    physics_types = {}
    for name, template in TEMPLATES.items():
        if template.physics_type not in physics_types:
            physics_types[template.physics_type] = []
        physics_types[template.physics_type].append((name, template))
    
    # Display templates by physics type
    for physics_type, templates in physics_types.items():
        st.markdown(f"### {physics_type.replace('_', ' ').title()}")
        
        cols = st.columns(3)
        for i, (name, template) in enumerate(templates):
            with cols[i % 3]:
                with st.container():
                    st.markdown(f"""
                    <div class="agent-card">
                        <h4>{name}</h4>
                        <p>{template.description}</p>
                        <p><strong>Physics Type:</strong> {template.physics_type}</p>
                        <p><strong>Complexity:</strong> {template.complexity}</p>
                        {'<p>⭐ Golden Template</p>' if template.is_golden else ''}
                    </div>
                    """, unsafe_allow_html=True)

with tab4:
    st.markdown("## 📖 Documentation")
    
    st.markdown("""
    ### How DARTSGPT Works
    
    DARTSGPT uses a sophisticated multi-agent system to transform your natural language 
    descriptions into working DARTS reservoir simulation models:
    
    1. **🎯 Intent Classifier**: Analyzes your prompt to understand the physics type and requirements
    2. **📋 Template Selector**: Chooses the most appropriate template from our curated database
    3. **🔢 Parameter Extractor**: Extracts numerical values and units from your description
    4. **💻 Code Generator**: Creates complete, commented Python code for your model
    5. **✅ Validator**: Ensures the generated code is syntactically correct and complete
    
    ### Supported Physics Types
    
    - **Compositional**: CO2 injection, gas injection, EOS modeling
    - **Dead Oil**: Water flooding, immiscible displacement
    - **Black Oil**: PVT relations, solution gas, volatile oil
    - **Geothermal**: Temperature modeling, heat transfer
    - **Poroelastic**: Geomechanics, stress-strain coupling
    - **Chemical**: Surfactant/polymer flooding, EOR chemicals
    
    ### Tips for Best Results
    
    1. Be specific about your requirements
    2. Include numerical values with units
    3. Mention the physics type if known
    4. Describe well configurations
    5. Specify simulation objectives
    
    ### Example Prompts
    
    - "Create a CO2 injection model for carbon storage in a 100x100x20 reservoir with 25% porosity"
    - "Build a geothermal model at 250°C with production at 150 kg/s"
    - "Simulate water flooding in a 50x50x10 grid with horizontal wells"
    """)


# Handle prompt reuse
if 'reuse_prompt' in st.session_state:
    st.session_state.prompt = st.session_state.reuse_prompt
    del st.session_state.reuse_prompt
    st.rerun()

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #7f8c8d; padding: 20px;">
    <p>DARTSGPT - Democratizing Reservoir Simulation with AI</p>
    <p>Made with ❤️ using Streamlit and LangChain</p>
</div>
""", unsafe_allow_html=True)