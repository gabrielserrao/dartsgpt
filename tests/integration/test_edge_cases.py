#!/usr/bin/env python
"""Edge case tests for DARTSGPT system."""

import json
from pathlib import Path
from rich.console import Console

from orchestrator import run_dartsgpt_orchestrator

console = Console()


def test_minimal_prompts():
    """Test with very minimal prompts."""
    console.print("\n[bold blue]Testing Minimal Prompts[/bold blue]")
    
    minimal_prompts = [
        "model",
        "reservoir",
        "oil",
        "water",
        "gas",
        "CO2",
        "temperature",
        "injection",
        "production",
        "simulate"
    ]
    
    results = []
    for prompt in minimal_prompts:
        console.print(f"\n[yellow]Testing: '{prompt}'[/yellow]")
        try:
            result = run_dartsgpt_orchestrator(prompt, f"test_output/minimal_{prompt}")
            results.append({
                'prompt': prompt,
                'success': result.get('success', False),
                'template': result.get('template_used', 'N/A')
            })
            console.print(f"  Result: {'✅' if result.get('success') else '❌'}")
        except Exception as e:
            console.print(f"  Error: {str(e)}")
            results.append({'prompt': prompt, 'success': False, 'template': 'ERROR'})
    
    return results


def test_conflicting_physics():
    """Test prompts with conflicting physics requirements."""
    console.print("\n[bold blue]Testing Conflicting Physics[/bold blue]")
    
    conflicting_prompts = [
        "Create a dead oil model with compositional CO2 injection",
        "Build a black oil reservoir with geothermal temperature effects",
        "Model a compositional system using dead oil physics",
        "Design a geothermal reservoir with black oil PVT",
        "Create a chemical EOR model with poroelastic effects"
    ]
    
    for prompt in conflicting_prompts:
        console.print(f"\n[yellow]Testing: {prompt[:60]}...[/yellow]")
        try:
            result = run_dartsgpt_orchestrator(prompt, f"test_output/conflict_{len(prompt)}")
            physics = result.get('intent', {}).get('physics_type', 'unknown')
            template = result.get('template_used', 'none')
            console.print(f"  Physics selected: {physics}")
            console.print(f"  Template used: {template}")
        except Exception as e:
            console.print(f"  Handled gracefully: {type(e).__name__}")


def test_extreme_parameters():
    """Test with extreme parameter values."""
    console.print("\n[bold blue]Testing Extreme Parameters[/bold blue]")
    
    extreme_prompts = [
        "Create a model with 0% porosity",
        "Build a reservoir at 1000°C temperature",
        "Model with permeability of 0.0001 mD",
        "Design a grid with 1000x1000x1000 cells",
        "Create a well with -100 m3/day injection rate",
        "Model a reservoir at -50°C",
        "Build a system with 200% porosity",
        "Create a model with 1,000,000 bar pressure"
    ]
    
    for prompt in extreme_prompts:
        console.print(f"\n[yellow]Testing: {prompt}[/yellow]")
        try:
            result = run_dartsgpt_orchestrator(prompt, f"test_output/extreme_{len(prompt)}")
            params = result.get('parameters', {})
            console.print(f"  Parameters extracted: {json.dumps(params, indent=2)}")
        except Exception as e:
            console.print(f"  Error handled: {str(e)}")


def test_unit_variations():
    """Test various unit specifications."""
    console.print("\n[bold blue]Testing Unit Variations[/bold blue]")
    
    unit_prompts = [
        "Create a model with 100 millidarcy permeability",
        "Build a reservoir with 500 md perm",
        "Model at 150 degrees celsius",
        "Design with 200 F temperature",
        "Create a well at 5000 feet depth",
        "Model with 2000 psi pressure",
        "Build with 10 km x 10 km area",
        "Create cells of 100m × 100m × 10m"
    ]
    
    for prompt in unit_prompts:
        console.print(f"\n[yellow]Testing: {prompt}[/yellow]")
        try:
            result = run_dartsgpt_orchestrator(prompt, f"test_output/units_{len(prompt)}")
            params = result.get('parameters', {})
            # Check if units were converted
            if 'rock_properties' in params:
                console.print(f"  Rock properties: {params['rock_properties']}")
            if 'temperature' in params:
                console.print(f"  Temperature: {params['temperature']}")
        except Exception as e:
            console.print(f"  Error: {str(e)}")


def test_multilingual_keywords():
    """Test with non-English characters and mixed languages."""
    console.print("\n[bold blue]Testing Multilingual Keywords[/bold blue]")
    
    multilingual_prompts = [
        "Create a reservoir with φ = 0.2",  # Greek phi for porosity
        "Model with μ = 1 cp viscosity",    # Greek mu
        "Build a system with Δp = 100 bar", # Greek delta
        "Create modelo with inyección de CO2",  # Spanish
        "Modèle de réservoir géothermique",     # French
        "Model with 渗透率 = 100 mD"              # Chinese for permeability
    ]
    
    for prompt in multilingual_prompts:
        console.print(f"\n[yellow]Testing: {prompt}[/yellow]")
        try:
            result = run_dartsgpt_orchestrator(prompt, f"test_output/multi_{len(prompt)}")
            console.print(f"  Success: {'✅' if result.get('success') else '❌'}")
        except Exception as e:
            console.print(f"  Handled: {type(e).__name__}")


def test_complex_well_patterns():
    """Test various well pattern descriptions."""
    console.print("\n[bold blue]Testing Complex Well Patterns[/bold blue]")
    
    well_patterns = [
        "Create a five-spot pattern with 4 injectors and 1 producer",
        "Build a line drive with 3 injectors and 3 producers",
        "Model with inverted nine-spot pattern",
        "Design peripheral water injection with 8 wells",
        "Create a horizontal well at 2000m depth",
        "Model with multilateral wells",
        "Build a smart well with ICVs",
        "Create deviated wells at 45 degree angle"
    ]
    
    for prompt in well_patterns:
        console.print(f"\n[yellow]Testing: {prompt}[/yellow]")
        try:
            result = run_dartsgpt_orchestrator(prompt, f"test_output/wells_{len(prompt)}")
            # Check generated well configuration
            if result.get('success'):
                model_code = result.get('model_code', '')
                well_count = model_code.count('add_well')
                console.print(f"  Wells added: {well_count}")
        except Exception as e:
            console.print(f"  Error: {str(e)}")


def test_time_dependent_scenarios():
    """Test time-dependent and dynamic scenarios."""
    console.print("\n[bold blue]Testing Time-Dependent Scenarios[/bold blue]")
    
    time_scenarios = [
        "Model 10 years of production followed by 5 years of injection",
        "Create a WAG cycle with 3 months water, 3 months gas",
        "Build a model with declining production rate over time",
        "Design injection that starts after 2 years",
        "Model with seasonal temperature variations",
        "Create a system with time-varying BHP",
        "Build a model simulating 50 years"
    ]
    
    for prompt in time_scenarios:
        console.print(f"\n[yellow]Testing: {prompt[:60]}...[/yellow]")
        try:
            result = run_dartsgpt_orchestrator(prompt, f"test_output/time_{len(prompt)}")
            main_code = result.get('main_code', '')
            if 'final_time' in main_code:
                console.print("  ✅ Simulation time parameters found")
            else:
                console.print("  ⚠️  Default simulation time used")
        except Exception as e:
            console.print(f"  Error: {str(e)}")


def test_invalid_combinations():
    """Test physically invalid or impossible combinations."""
    console.print("\n[bold blue]Testing Invalid Combinations[/bold blue]")
    
    invalid_prompts = [
        "Create a model with injection rate > production rate in closed system",
        "Build a reservoir with porosity > 1",
        "Model with negative absolute permeability",
        "Design a system with temperature below absolute zero",
        "Create incompressible fluid with gas",
        "Model water injection into gas reservoir at 100% gas saturation",
        "Build a reservoir above ground surface"
    ]
    
    for prompt in invalid_prompts:
        console.print(f"\n[yellow]Testing: {prompt[:60]}...[/yellow]")
        try:
            result = run_dartsgpt_orchestrator(prompt, f"test_output/invalid_{len(prompt)}")
            validation = result.get('validation', {})
            if validation.get('issues'):
                console.print(f"  ⚠️  Validation issues: {validation['issues']}")
            else:
                console.print("  ✅ Generated (may need manual validation)")
        except Exception as e:
            console.print(f"  Error caught: {type(e).__name__}")


def run_all_edge_tests():
    """Run all edge case tests."""
    console.print("[bold cyan]DARTSGPT Edge Case Test Suite[/bold cyan]")
    console.print("=" * 60)
    
    Path("test_output").mkdir(exist_ok=True)
    
    test_functions = [
        ("Minimal Prompts", test_minimal_prompts),
        ("Conflicting Physics", test_conflicting_physics),
        ("Extreme Parameters", test_extreme_parameters),
        ("Unit Variations", test_unit_variations),
        ("Multilingual Keywords", test_multilingual_keywords),
        ("Complex Well Patterns", test_complex_well_patterns),
        ("Time-Dependent Scenarios", test_time_dependent_scenarios),
        ("Invalid Combinations", test_invalid_combinations)
    ]
    
    summary = {}
    
    for name, test_func in test_functions:
        console.print(f"\n[bold green]Running: {name}[/bold green]")
        console.print("-" * 40)
        try:
            test_func()
            summary[name] = "✅ Completed"
        except Exception as e:
            summary[name] = f"❌ Failed: {str(e)}"
            console.print(f"[red]Test suite failed: {str(e)}[/red]")
    
    # Print summary
    console.print("\n[bold cyan]Test Summary[/bold cyan]")
    console.print("=" * 60)
    for test_name, status in summary.items():
        console.print(f"{test_name}: {status}")
    
    console.print("\n[bold green]Edge case testing completed![/bold green]")


if __name__ == "__main__":
    run_all_edge_tests()