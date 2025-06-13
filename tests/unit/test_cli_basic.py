#!/usr/bin/env python
"""Basic CLI test for DARTSGPT functionality."""

import click
from pathlib import Path
import json

from knowledge.templates.template_database import TEMPLATES
from knowledge.physics_features import PHYSICS_KEYWORDS
from knowledge.example_prompts import EXAMPLE_PROMPTS
from config import get_settings


@click.group()
def cli():
    """DARTSGPT - AI-powered DARTS model generation (Basic Test)."""
    pass


@cli.command()
def list_templates():
    """List available DARTS templates."""
    click.echo("\nAvailable DARTS Templates:")
    click.echo("=" * 60)
    
    for name, template in TEMPLATES.items():
        star = "⭐" if template.is_golden else "  "
        click.echo(f"{star} {name:20} - {template.description}")
        click.echo(f"    Physics: {template.physics_type:15} Complexity: {template.complexity}")
    
    click.echo(f"\nTotal: {len(TEMPLATES)} templates")


@cli.command()
@click.option('--physics', default='compositional', help='Physics type')
def show_examples(physics):
    """Show example prompts for a physics type."""
    if physics not in EXAMPLE_PROMPTS:
        click.echo(f"No examples found for physics type: {physics}")
        click.echo(f"Available types: {', '.join(EXAMPLE_PROMPTS.keys())}")
        return
    
    examples = EXAMPLE_PROMPTS[physics]
    click.echo(f"\nExample prompts for {physics.upper()} physics:")
    click.echo("=" * 60)
    
    if 'simple' in examples:
        click.echo("\nSimple prompts:")
        for prompt in examples['simple']:
            click.echo(f"  • {prompt}")
    
    if 'detailed' in examples:
        click.echo("\nDetailed prompts:")
        for prompt in examples['detailed']:
            click.echo(f"  • {prompt}")


@cli.command()
def check_config():
    """Check configuration and API settings."""
    settings = get_settings()
    
    click.echo("\nDARTSGPT Configuration:")
    click.echo("=" * 60)
    click.echo(f"OpenAI Model: {settings.openai_model}")
    
    api_key_set = bool(settings.openai_api_key)
    if api_key_set:
        click.echo("OpenAI API Key: ✓ Set")
    else:
        click.echo("OpenAI API Key: ✗ Not set")
        click.echo("\nTo use LLM features, set OPENAI_API_KEY in your .env file")
    
    click.echo(f"\nVector Store Path: {settings.vector_store_path}")
    click.echo(f"Log Level: {settings.log_level}")


@cli.command()
@click.argument('prompt')
def test_intent(prompt):
    """Test intent classification on a prompt (simple rule-based)."""
    prompt_lower = prompt.lower()
    
    # Detect physics type
    physics_type = "dead_oil"  # default
    for ptype, keywords in PHYSICS_KEYWORDS.items():
        if any(kw in prompt_lower for kw in keywords):
            physics_type = ptype
            break
    
    # Detect features
    features = []
    if "injection" in prompt_lower:
        features.append("injection")
    if "production" in prompt_lower:
        features.append("production")
    if "thermal" in prompt_lower or "temperature" in prompt_lower:
        features.append("thermal")
    
    click.echo(f"\nIntent Analysis for: {prompt}")
    click.echo("=" * 60)
    click.echo(f"Detected physics type: {physics_type}")
    click.echo(f"Detected features: {features}")


if __name__ == "__main__":
    cli()