"""Command Line Interface for DARTSGPT."""

import asyncio
import json
import sys
from pathlib import Path
from typing import Optional
import click

from .agents.intent_classifier import IntentClassifierAgent
from .agents.template_selector import TemplateSelectorAgent
from .agents.parameter_extractor import ParameterExtractorAgent
from .agents.code_generator import CodeGeneratorAgent
from .agents.validator import ValidatorAgent
from .agents.base_agent import AgentInput
from .knowledge.rag_system import RAGSystem
from .config import get_settings


@click.group()
def cli():
    """DARTSGPT - AI-powered DARTS model generation."""
    pass


@cli.command()
@click.argument('prompt')
@click.option('--output-dir', '-o', default='output', help='Output directory for generated files')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--no-validation', is_flag=True, help='Skip validation step')
async def generate(prompt: str, output_dir: str, verbose: bool, no_validation: bool):
    """Generate DARTS model from natural language prompt."""
    settings = get_settings()
    
    if verbose:
        click.echo(f"Using model: {settings.openai_model}")
        click.echo(f"Processing prompt: {prompt}\n")
    
    # Initialize agents
    intent_agent = IntentClassifierAgent()
    template_agent = TemplateSelectorAgent()
    param_agent = ParameterExtractorAgent()
    code_agent = CodeGeneratorAgent()
    validator_agent = ValidatorAgent()
    
    # Process through agent pipeline
    context = {}
    
    # 1. Intent Classification
    click.echo("🔍 Analyzing intent...")
    intent_input = AgentInput(message=prompt, context=context)
    intent_result = await intent_agent.process(intent_input)
    context['intent_classification'] = intent_result.result
    
    if verbose:
        click.echo(f"  Physics type: {intent_result.result.get('physics_type')}")
        click.echo(f"  Features: {intent_result.result.get('features')}")
        click.echo(f"  Confidence: {intent_result.confidence:.2f}\n")
    
    # 2. Template Selection
    click.echo("📋 Selecting template...")
    template_input = AgentInput(message=prompt, context=context)
    template_result = await template_agent.process(template_input)
    context['template_selection'] = template_result.result
    
    if verbose:
        click.echo(f"  Selected: {template_result.result.get('template_name')}")
        click.echo(f"  Score: {template_result.result.get('overall_score'):.2f}\n")
    
    # 3. Parameter Extraction
    click.echo("🔢 Extracting parameters...")
    param_input = AgentInput(message=prompt, context=context)
    param_result = await param_agent.process(param_input)
    context['parameter_extraction'] = param_result.result
    
    if verbose:
        grid_params = param_result.result.get('grid_parameters', {})
        click.echo(f"  Grid: {grid_params.get('nx', 'N/A')}x{grid_params.get('ny', 'N/A')}x{grid_params.get('nz', 'N/A')}")
        rock_props = param_result.result.get('rock_properties', {})
        click.echo(f"  Porosity: {rock_props.get('porosity', 'N/A')}")
        click.echo(f"  Permeability: {rock_props.get('permeability', 'N/A')} m²\n")
    
    # 4. Code Generation
    click.echo("💻 Generating code...")
    code_input = AgentInput(message=prompt, context=context)
    code_result = await code_agent.process(code_input)
    context['code_generation'] = code_result.result
    
    # 5. Validation (optional)
    if not no_validation:
        click.echo("✅ Validating code...")
        valid_input = AgentInput(message=prompt, context=context)
        valid_result = await validator_agent.process(valid_input)
        context['validation'] = valid_result.result
        
        if verbose:
            click.echo(f"  Valid: {valid_result.result.get('is_valid')}")
            issues = valid_result.result.get('issues', [])
            if issues:
                click.echo(f"  Issues: {len(issues)} found")
                for issue in issues[:3]:  # Show first 3 issues
                    click.echo(f"    - {issue['severity']}: {issue['message']}")
    
    # Save output
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Save model.py
    model_code = code_result.result.get('model_code', '')
    model_file = output_path / 'model.py'
    model_file.write_text(model_code)
    click.echo(f"\n✨ Generated model.py saved to: {model_file}")
    
    # Save main.py
    main_code = code_result.result.get('main_code', '')
    main_file = output_path / 'main.py'
    main_file.write_text(main_code)
    click.echo(f"✨ Generated main.py saved to: {main_file}")
    
    # Save metadata
    metadata = {
        'prompt': prompt,
        'intent_classification': context.get('intent_classification'),
        'template_selection': context.get('template_selection'),
        'parameter_extraction': context.get('parameter_extraction'),
        'validation': context.get('validation')
    }
    metadata_file = output_path / 'metadata.json'
    metadata_file.write_text(json.dumps(metadata, indent=2))
    
    if verbose:
        click.echo(f"\n📊 Metadata saved to: {metadata_file}")
    
    click.echo("\n✅ Generation complete!")


@cli.command()
def list_templates():
    """List available DARTS templates."""
    from .knowledge.templates.template_database import TEMPLATES
    
    click.echo("Available DARTS Templates:\n")
    
    for name, template in TEMPLATES.items():
        star = "⭐" if template.is_golden else "  "
        click.echo(f"{star} {name:20} - {template.description}")
        click.echo(f"    Physics: {template.physics_type:15} Complexity: {template.complexity}")
        if template.use_cases:
            click.echo(f"    Use cases: {', '.join(template.use_cases[:2])}")
        click.echo()


@cli.command()
@click.argument('physics_type', required=False)
def show_examples(physics_type: Optional[str]):
    """Show example prompts for physics types."""
    from .knowledge.example_prompts import EXAMPLE_PROMPTS
    
    if physics_type:
        if physics_type in EXAMPLE_PROMPTS:
            examples = EXAMPLE_PROMPTS[physics_type]
            click.echo(f"\nExamples for {physics_type.upper()} physics:\n")
            
            click.echo("Simple prompts:")
            for prompt in examples.get('simple', []):
                click.echo(f"  • {prompt}")
            
            click.echo("\nDetailed prompts:")
            for prompt in examples.get('detailed', []):
                click.echo(f"  • {prompt}")
        else:
            click.echo(f"Unknown physics type: {physics_type}")
            click.echo("Available types: " + ", ".join(EXAMPLE_PROMPTS.keys()))
    else:
        click.echo("Available physics types:")
        for ptype in EXAMPLE_PROMPTS.keys():
            click.echo(f"  • {ptype}")
        click.echo("\nUse 'dartsgpt show-examples <physics_type>' for specific examples")


@cli.command()
@click.argument('prompt')
async def test_intent(prompt: str):
    """Test intent classification for a prompt."""
    settings = get_settings()
    agent = IntentClassifierAgent()
    
    click.echo(f"Testing intent classification for: {prompt}\n")
    
    input_data = AgentInput(message=prompt, context={})
    result = await agent.process(input_data)
    
    click.echo("Results:")
    click.echo(f"  Physics type: {result.result.get('physics_type')}")
    click.echo(f"  Features: {result.result.get('features')}")
    click.echo(f"  Complexity: {result.result.get('complexity')}")
    click.echo(f"  Confidence: {result.confidence:.2f}")
    
    if result.result.get('grid_info'):
        click.echo(f"  Grid info: {result.result['grid_info']}")
    if result.result.get('well_info'):
        click.echo(f"  Well info: {result.result['well_info']}")


@cli.command()
@click.option('--test-db', is_flag=True, help='Test database connection')
@click.option('--test-agents', is_flag=True, help='Test agent initialization')
@click.option('--test-templates', is_flag=True, help='Test template loading')
async def check_system(test_db: bool, test_agents: bool, test_templates: bool):
    """Check system configuration and dependencies."""
    from .config import get_settings
    
    click.echo("DARTSGPT System Check\n")
    
    # Check environment
    settings = get_settings()
    click.echo("✅ Environment loaded")
    click.echo(f"  Model: {settings.openai_model}")
    click.echo(f"  API Key: {'*' * 20}...{settings.openai_api_key[-4:]}")
    
    # Check templates
    if test_templates:
        from .knowledge.templates.template_database import TEMPLATES
        click.echo(f"\n✅ Templates loaded: {len(TEMPLATES)} available")
    
    # Check agents
    if test_agents:
        try:
            intent = IntentClassifierAgent()
            template = TemplateSelectorAgent()
            param = ParameterExtractorAgent()
            code = CodeGeneratorAgent()
            valid = ValidatorAgent()
            click.echo("\n✅ All agents initialized successfully")
        except Exception as e:
            click.echo(f"\n❌ Agent initialization failed: {e}")
    
    # Check database
    if test_db:
        try:
            rag = RAGSystem()
            click.echo("\n✅ RAG system initialized")
        except Exception as e:
            click.echo(f"\n❌ RAG system failed: {e}")


def main():
    """Main entry point."""
    cli()


if __name__ == '__main__':
    # Handle async commands
    if len(sys.argv) > 1 and sys.argv[1] in ['generate', 'test-intent', 'check-system']:
        # Get the command function
        ctx = cli.make_context('dartsgpt', sys.argv[1:])
        cmd = ctx.command.callback
        
        # Run async command
        if asyncio.iscoroutinefunction(cmd):
            asyncio.run(cmd(**ctx.params))
        else:
            cmd(**ctx.params)
    else:
        cli()