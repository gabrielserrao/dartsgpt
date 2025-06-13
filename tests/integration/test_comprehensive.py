#!/usr/bin/env python
"""Comprehensive tests for DARTSGPT multi-agent system."""

import json
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

from orchestrator import run_dartsgpt_orchestrator
from knowledge.templates.template_database import TEMPLATES, TemplateDatabase
from knowledge.physics_features import PHYSICS_KEYWORDS
from config import get_settings

console = Console()


def test_all_physics_types():
    """Test generation for all physics types."""
    console.print(Panel("Testing All Physics Types", style="bold blue"))
    
    test_cases = {
        "compositional": [
            "Create a CO2 injection model for carbon storage with 5 components",
            "Build a compositional gas injection model in a 100x100x10 reservoir",
            "Model a multi-component hydrocarbon system with 25% porosity"
        ],
        "dead_oil": [
            "Create a waterflood model with five-spot pattern",
            "Build a simple oil reservoir with 20% porosity and 100 mD permeability",
            "Model an oil field with 50x50x5 grid and water injection"
        ],
        "black_oil": [
            "Create a three-phase black oil model with gas cap",
            "Build a black oil reservoir with PVT properties",
            "Model a black oil field with solution gas and water drive"
        ],
        "geothermal": [
            "Model a geothermal reservoir at 200°C with injection wells",
            "Create a two-phase geothermal system with steam",
            "Build an enhanced geothermal system with cold water injection"
        ],
        "poroelastic": [
            "Create a poroelastic model for geomechanics studies",
            "Model reservoir compaction with poroelastic effects",
            "Build a coupled flow-geomechanics model"
        ],
        "chemical": [
            "Create a chemical EOR model with surfactant flooding",
            "Model polymer injection for enhanced oil recovery",
            "Build an alkaline-surfactant-polymer flooding model"
        ]
    }
    
    results = []
    
    for physics_type, prompts in test_cases.items():
        console.print(f"\n[yellow]Testing {physics_type.upper()} physics:[/yellow]")
        
        for i, prompt in enumerate(prompts):
            try:
                output_dir = f"test_output/{physics_type}_{i+1}"
                result = run_dartsgpt_orchestrator(prompt, output_dir)
                
                success = result.get('success', False)
                template = result.get('template_used', 'N/A')
                
                results.append({
                    'physics': physics_type,
                    'prompt': prompt[:50] + '...',
                    'success': success,
                    'template': template
                })
                
                console.print(f"  ✓ Test {i+1}: {'✅ Passed' if success else '❌ Failed'}")
                
            except Exception as e:
                console.print(f"  ✗ Test {i+1}: ❌ Error - {str(e)}")
                results.append({
                    'physics': physics_type,
                    'prompt': prompt[:50] + '...',
                    'success': False,
                    'template': 'ERROR'
                })
    
    # Display results summary
    display_test_results(results)
    return results


def test_parameter_extraction():
    """Test parameter extraction capabilities."""
    console.print(Panel("Testing Parameter Extraction", style="bold blue"))
    
    test_prompts = [
        {
            "prompt": "Create a 100x100x20 reservoir with 25% porosity and 500 mD permeability",
            "expected": {
                "grid": {"nx": 100, "ny": 100, "nz": 20},
                "rock": {"porosity": 0.25, "permeability": 500}
            }
        },
        {
            "prompt": "Model a reservoir at 150°C with injection rate of 1000 m3/day",
            "expected": {
                "temperature": 150,
                "well": {"injection_rate": 1000}
            }
        },
        {
            "prompt": "Build a model with pressure from 200 to 300 bar",
            "expected": {
                "pressure": {"min": 200, "max": 300}
            }
        }
    ]
    
    for i, test in enumerate(test_prompts):
        console.print(f"\n[yellow]Test {i+1}:[/yellow] {test['prompt']}")
        
        try:
            result = run_dartsgpt_orchestrator(test['prompt'], f"test_output/param_test_{i+1}")
            params = result.get('parameters', {})
            
            console.print(f"  Extracted: {json.dumps(params, indent=2)}")
            console.print(f"  Expected: {json.dumps(test['expected'], indent=2)}")
            
        except Exception as e:
            console.print(f"  ❌ Error: {str(e)}")


def test_template_coverage():
    """Test that all templates can be selected and used."""
    console.print(Panel("Testing Template Coverage", style="bold blue"))
    
    template_db = TemplateDatabase()
    
    table = Table(title="Template Coverage Test")
    table.add_column("Template", style="cyan")
    table.add_column("Physics", style="yellow")
    table.add_column("Tested", style="green")
    table.add_column("Result", style="magenta")
    
    for name, template in TEMPLATES.items():
        # Create a prompt that should trigger this template
        physics = template.physics_type
        prompt = f"Create a {physics.replace('_', ' ')} model using advanced features"
        
        try:
            result = run_dartsgpt_orchestrator(prompt, f"test_output/template_{name}")
            selected = result.get('template_used', '')
            success = selected == name or result.get('success', False)
            
            table.add_row(
                name,
                physics,
                "✓",
                "✅ Matched" if selected == name else f"➡️ {selected}"
            )
            
        except Exception as e:
            table.add_row(name, physics, "✗", "❌ Error")
    
    console.print(table)


def test_edge_cases():
    """Test edge cases and error handling."""
    console.print(Panel("Testing Edge Cases", style="bold blue"))
    
    edge_cases = [
        "Create a model",  # Minimal prompt
        "Build a reservoir with negative porosity of -10%",  # Invalid parameter
        "Model a system at 1000000°C",  # Extreme value
        "Create a model with 1x1x1 grid",  # Tiny grid
        "Build a reservoir with 1000x1000x1000 grid",  # Huge grid
        "",  # Empty prompt
        "Generate Python code for data analysis",  # Off-topic
        "Model with πορosity of 0.2",  # Unicode
    ]
    
    for i, prompt in enumerate(edge_cases):
        console.print(f"\n[yellow]Edge case {i+1}:[/yellow] '{prompt}'")
        
        try:
            result = run_dartsgpt_orchestrator(prompt, f"test_output/edge_{i+1}")
            success = result.get('success', False)
            validation = result.get('validation', {})
            
            console.print(f"  Result: {'✅ Handled' if success else '⚠️  Issues found'}")
            if validation.get('issues'):
                for issue in validation['issues']:
                    console.print(f"    - {issue}")
                    
        except Exception as e:
            console.print(f"  ❌ Exception: {type(e).__name__}: {str(e)}")


def test_complex_scenarios():
    """Test complex, real-world scenarios."""
    console.print(Panel("Testing Complex Scenarios", style="bold blue"))
    
    scenarios = [
        {
            "name": "CO2 Storage Project",
            "prompt": """Create a CO2 injection model for carbon storage project with:
            - 200x200x50 grid representing 10km x 10km x 100m
            - Heterogeneous porosity from 0.15 to 0.25
            - Permeability ranging from 10 to 1000 mD
            - Multiple injection wells at 2000 m depth
            - Monitoring wells for pressure observation
            - Initial pressure of 200 bar at datum
            - Temperature of 60°C""",
        },
        {
            "name": "Geothermal Power Plant",
            "prompt": """Model a geothermal reservoir for power generation:
            - Two-phase water/steam system at 250°C
            - Fractured reservoir with dual porosity
            - Production wells targeting steam zones
            - Reinjection wells for sustainability
            - Natural recharge from boundaries
            - 30-year production forecast""",
        },
        {
            "name": "Enhanced Oil Recovery",
            "prompt": """Build a chemical EOR model for mature field:
            - Black oil reservoir with 15% current recovery
            - ASP (alkaline-surfactant-polymer) flooding
            - Five-spot pattern with 40-acre spacing
            - Heterogeneous layers with flow barriers
            - Economic optimization of chemical slug size
            - History matching with 20 years production data""",
        }
    ]
    
    for scenario in scenarios:
        console.print(f"\n[bold green]{scenario['name']}:[/bold green]")
        console.print(f"[dim]{scenario['prompt'][:100]}...[/dim]")
        
        try:
            output_dir = f"test_output/scenario_{scenario['name'].replace(' ', '_')}"
            result = run_dartsgpt_orchestrator(scenario['prompt'], output_dir)
            
            if result.get('success'):
                console.print("  ✅ Successfully generated")
                console.print(f"  Template: {result.get('template_used')}")
                console.print(f"  Parameters extracted: {len(result.get('parameters', {}))}")
            else:
                console.print("  ⚠️  Generation completed with issues")
                
        except Exception as e:
            console.print(f"  ❌ Error: {str(e)}")


def display_test_results(results):
    """Display test results in a summary table."""
    table = Table(title="Test Results Summary")
    table.add_column("Physics Type", style="cyan")
    table.add_column("Tests Run", style="yellow")
    table.add_column("Passed", style="green")
    table.add_column("Failed", style="red")
    table.add_column("Success Rate", style="magenta")
    
    # Group by physics type
    physics_results = {}
    for result in results:
        physics = result['physics']
        if physics not in physics_results:
            physics_results[physics] = {'total': 0, 'passed': 0}
        
        physics_results[physics]['total'] += 1
        if result['success']:
            physics_results[physics]['passed'] += 1
    
    # Display results
    total_tests = 0
    total_passed = 0
    
    for physics, stats in physics_results.items():
        total = stats['total']
        passed = stats['passed']
        failed = total - passed
        rate = (passed / total * 100) if total > 0 else 0
        
        table.add_row(
            physics.title(),
            str(total),
            str(passed),
            str(failed),
            f"{rate:.1f}%"
        )
        
        total_tests += total
        total_passed += passed
    
    # Add total row
    table.add_row(
        "[bold]TOTAL[/bold]",
        f"[bold]{total_tests}[/bold]",
        f"[bold]{total_passed}[/bold]",
        f"[bold]{total_tests - total_passed}[/bold]",
        f"[bold]{(total_passed/total_tests*100) if total_tests > 0 else 0:.1f}%[/bold]"
    )
    
    console.print("\n")
    console.print(table)


def run_all_tests():
    """Run all test suites."""
    console.print(Panel(
        "[bold]DARTSGPT Comprehensive Test Suite[/bold]\n" +
        "Testing multi-agent orchestration system",
        style="bold blue"
    ))
    
    # Check configuration
    settings = get_settings()
    if not settings.openai_api_key:
        console.print("[red]⚠️  Warning: No OpenAI API key found![/red]")
        console.print("Tests will use fallback methods.\n")
    
    # Create test output directory
    Path("test_output").mkdir(exist_ok=True)
    
    # Run test suites
    test_suites = [
        ("Physics Types", test_all_physics_types),
        ("Parameter Extraction", test_parameter_extraction),
        ("Template Coverage", test_template_coverage),
        ("Edge Cases", test_edge_cases),
        ("Complex Scenarios", test_complex_scenarios)
    ]
    
    for name, test_func in test_suites:
        console.print(f"\n[bold cyan]Running: {name}[/bold cyan]")
        console.print("=" * 60)
        
        try:
            test_func()
        except Exception as e:
            console.print(f"[red]Test suite failed: {str(e)}[/red]")
        
        console.print("\n")
    
    console.print(Panel(
        "[bold green]Test Suite Completed![/bold green]\n" +
        "Check test_output/ directory for generated files",
        style="bold green"
    ))


if __name__ == "__main__":
    run_all_tests()