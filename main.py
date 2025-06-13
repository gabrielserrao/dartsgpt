#!/usr/bin/env python
"""Main CLI script for DARTSGPT with multi-agent orchestration."""

import click
import asyncio
from pathlib import Path
import json
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table
from rich.panel import Panel
from rich.syntax import Syntax

from orchestrator import run_dartsgpt_orchestrator, build_dartsgpt_graph
from config import get_settings
from knowledge.templates.template_database import TEMPLATES
from knowledge.example_prompts import EXAMPLE_PROMPTS

console = Console()


@click.group()
def cli():
    """🎯 DARTSGPT - AI-powered DARTS model generation with multi-agent orchestration."""
    pass


@cli.command()
@click.argument('prompt')
@click.option('--output-dir', '-o', default='output', help='Output directory for generated files')
@click.option('--verbose', '-v', is_flag=True, help='Verbose output')
@click.option('--dry-run', is_flag=True, help='Show what would be generated without saving files')
@click.option('--execute', '-e', is_flag=True, help='Execute the generated model after creation')
def generate(prompt: str, output_dir: str, verbose: bool, dry_run: bool, execute: bool):
    """Generate DARTS model from natural language prompt using multi-agent orchestration."""
    
    # Check configuration
    settings = get_settings()
    if not settings.openai_api_key:
        console.print("[red]❌ Error: OpenAI API key not found![/red]")
        console.print("Please set OPENAI_API_KEY in your .env file")
        return
    
    console.print(f"[blue]🤖 Using model: {settings.openai_model}[/blue]\n")
    
    # Run the orchestrator with progress indicator
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        task = progress.add_task("Running multi-agent orchestration...", total=None)
        
        try:
            result = run_dartsgpt_orchestrator(
                prompt, 
                output_dir if not dry_run else None,
                execute=execute and not dry_run
            )
            progress.update(task, completed=True)
        except Exception as e:
            progress.stop()
            console.print(f"[red]❌ Error: {str(e)}[/red]")
            if verbose:
                import traceback
                console.print(traceback.format_exc())
            return
    
    # Display results
    if result.get('success'):
        console.print("\n[green]✅ Model generation completed successfully![/green]\n")
        
        # Show generated code preview
        if verbose or dry_run:
            console.print(Panel(
                Syntax(result['model_code'][:500] + "...", "python", theme="monokai"),
                title="Generated Model Code (preview)",
                border_style="green"
            ))
        
        # Show metadata
        table = Table(title="Generation Metadata")
        table.add_column("Property", style="cyan")
        table.add_column("Value", style="yellow")
        
        table.add_row("Template Used", result['template_used'])
        table.add_row("Physics Type", str(result['parameters'].get('physics_type', 'N/A')))
        table.add_row("Grid Size", str(result['parameters'].get('grid_parameters', {})))
        table.add_row("Validation", "✅ Passed" if result['validation'].get('issues', []) == [] else "⚠️  Has issues")
        
        # Add execution info if available
        if execute and result.get('execution'):
            exec_result = result['execution']
            exec_status = "✅ Success" if exec_result.get('success') else "❌ Failed"
            exec_time = f"{exec_result.get('execution_time', 0):.2f}s"
            table.add_row("Execution", f"{exec_status} ({exec_time})")
            
            if exec_result.get('summary', {}).get('simulation_completed'):
                table.add_row("Simulation", "✅ Completed")
        
        console.print(table)
        
        if not dry_run:
            console.print(f"\n[green]📁 Files saved to: {output_dir}/[/green]")
    else:
        console.print("[red]❌ Model generation failed![/red]")


@cli.command()
def list_templates():
    """List all available DARTS templates."""
    
    table = Table(title="Available DARTS Templates", show_lines=True)
    table.add_column("Template", style="cyan", width=20)
    table.add_column("Physics Type", style="yellow", width=15)
    table.add_column("Description", style="white", width=40)
    table.add_column("Golden", style="green", width=6, justify="center")
    
    for name, template in TEMPLATES.items():
        table.add_row(
            name,
            template.physics_type,
            template.description,
            "⭐" if template.is_golden else ""
        )
    
    console.print(table)
    console.print(f"\n[blue]Total templates: {len(TEMPLATES)}[/blue]")


@cli.command()
@click.option('--physics', default='compositional', help='Physics type to show examples for')
def show_examples(physics: str):
    """Show example prompts for a physics type."""
    
    if physics not in EXAMPLE_PROMPTS:
        console.print(f"[red]No examples found for physics type: {physics}[/red]")
        console.print(f"Available types: {', '.join(EXAMPLE_PROMPTS.keys())}")
        return
    
    examples = EXAMPLE_PROMPTS[physics]
    
    console.print(Panel(f"Example Prompts for [bold]{physics.upper()}[/bold] Physics", style="blue"))
    
    if 'simple' in examples:
        console.print("\n[yellow]Simple prompts:[/yellow]")
        for i, prompt in enumerate(examples['simple'], 1):
            console.print(f"  {i}. {prompt}")
    
    if 'detailed' in examples:
        console.print("\n[yellow]Detailed prompts:[/yellow]")
        for i, prompt in enumerate(examples['detailed'], 1):
            console.print(f"  {i}. {prompt}")


@cli.command()
def check_system():
    """Check system configuration and agent status."""
    
    settings = get_settings()
    
    console.print(Panel("DARTSGPT System Status", style="blue"))
    
    # Configuration status
    config_table = Table(title="Configuration", show_header=False)
    config_table.add_column("Setting", style="cyan")
    config_table.add_column("Value", style="yellow")
    
    config_table.add_row("OpenAI Model", settings.openai_model)
    config_table.add_row("API Key", "✅ Set" if settings.openai_api_key else "❌ Not set")
    config_table.add_row("Vector Store Path", str(settings.vector_store_path))
    
    console.print(config_table)
    
    # Agent status
    console.print("\n[blue]Multi-Agent System Components:[/blue]")
    agents = [
        ("🎯 Supervisor", "Orchestrates the entire workflow"),
        ("🔍 Intent Classifier", "Identifies physics type and features"),
        ("📋 Template Selector", "Chooses optimal DARTS template"),
        ("🔢 Parameter Extractor", "Extracts numerical parameters"),
        ("💻 Code Generator", "Generates DARTS model code"),
        ("✅ Validator", "Validates generated code")
    ]
    
    for agent, desc in agents:
        console.print(f"  {agent}: {desc}")
    
    # Test graph building
    try:
        graph = build_dartsgpt_graph()
        console.print("\n[green]✅ Multi-agent graph built successfully![/green]")
    except Exception as e:
        console.print(f"\n[red]❌ Error building graph: {str(e)}[/red]")


@cli.command()
@click.option('--interactive', '-i', is_flag=True, help='Interactive mode')
def demo(interactive: bool):
    """Run a demonstration of DARTSGPT capabilities."""
    
    demo_prompts = [
        "Create a CO2 injection model for carbon storage in a 100x100x20 reservoir with 25% porosity",
        "Build a waterflood model with five-spot pattern in a 50x50x10 grid",
        "Model a geothermal reservoir at 150°C with injection and production wells"
    ]
    
    if interactive:
        console.print("[blue]Interactive Demo Mode[/blue]")
        console.print("Select a demo prompt or enter your own:\n")
        
        for i, prompt in enumerate(demo_prompts, 1):
            console.print(f"  {i}. {prompt}")
        
        console.print(f"  {len(demo_prompts) + 1}. Enter custom prompt")
        
        choice = click.prompt("\nSelect option", type=int)
        
        if 1 <= choice <= len(demo_prompts):
            prompt = demo_prompts[choice - 1]
        else:
            prompt = click.prompt("Enter your prompt")
    else:
        # Run all demos
        console.print("[blue]Running all demo prompts...[/blue]\n")
        
        for i, prompt in enumerate(demo_prompts, 1):
            console.print(f"\n[yellow]Demo {i}/{len(demo_prompts)}[/yellow]")
            output_dir = f"demo_output/demo_{i}"
            
            try:
                result = run_dartsgpt_orchestrator(prompt, output_dir, execute=False)
                if result.get('success'):
                    console.print(f"[green]✅ Demo {i} completed![/green]")
            except Exception as e:
                console.print(f"[red]❌ Demo {i} failed: {str(e)}[/red]")
        
        console.print("\n[green]All demos completed! Check demo_output/ directory.[/green]")
        return
    
    # Run single demo
    output_dir = "demo_output/interactive"
    result = run_dartsgpt_orchestrator(prompt, output_dir, execute=False)
    
    if result.get('success'):
        console.print(f"\n[green]✅ Demo completed! Files saved to {output_dir}/[/green]")


if __name__ == "__main__":
    cli()