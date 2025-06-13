#!/usr/bin/env python
"""Run all DARTSGPT test suites."""

import sys
import subprocess
from pathlib import Path
from rich.console import Console
from rich.panel import Panel
from datetime import datetime

console = Console()


def run_test_file(test_file: str, description: str) -> dict:
    """Run a single test file and capture results."""
    console.print(f"\n[bold blue]Running: {description}[/bold blue]")
    console.print(f"File: {test_file}")
    console.print("-" * 60)
    
    start_time = datetime.now()
    
    try:
        # Run the test file
        result = subprocess.run(
            [sys.executable, test_file],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        # Display output
        if result.stdout:
            console.print(result.stdout)
        
        if result.stderr:
            console.print("[red]Errors:[/red]")
            console.print(result.stderr)
        
        return {
            "file": test_file,
            "description": description,
            "success": result.returncode == 0,
            "duration": duration,
            "return_code": result.returncode
        }
        
    except subprocess.TimeoutExpired:
        console.print("[red]Test timed out after 5 minutes[/red]")
        return {
            "file": test_file,
            "description": description,
            "success": False,
            "duration": 300,
            "error": "Timeout"
        }
    except Exception as e:
        console.print(f"[red]Error running test: {str(e)}[/red]")
        return {
            "file": test_file,
            "description": description,
            "success": False,
            "duration": 0,
            "error": str(e)
        }


def main():
    """Run all test suites."""
    console.print(Panel(
        "[bold cyan]DARTSGPT Complete Test Suite Runner[/bold cyan]\n" +
        f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        style="bold blue"
    ))
    
    # Ensure test output directory exists
    Path("test_output").mkdir(exist_ok=True)
    
    # Define test suites to run
    test_suites = [
        ("integration/test_comprehensive.py", "Comprehensive Physics Tests"),
        ("integration/test_edge_cases.py", "Edge Case Tests"),
        ("unit/test_unit_tests.py", "Unit Tests for Components"),
        ("integration/test_integration.py", "Integration Tests"),
        ("unit/test_cli.py", "CLI Tests"),
        ("unit/test_system.py", "System Tests"),
        ("integration/test_full_agent.py", "Full Agent Tests"),
    ]
    
    results = []
    total_duration = 0
    
    # Run each test suite
    for test_file, description in test_suites:
        if Path(test_file).exists():
            result = run_test_file(test_file, description)
            results.append(result)
            total_duration += result["duration"]
        else:
            console.print(f"[yellow]Warning: {test_file} not found, skipping[/yellow]")
            results.append({
                "file": test_file,
                "description": description,
                "success": False,
                "duration": 0,
                "error": "File not found"
            })
    
    # Display summary
    console.print("\n" + "=" * 80)
    console.print(Panel("[bold]Test Suite Summary[/bold]", style="cyan"))
    
    passed = sum(1 for r in results if r["success"])
    total = len(results)
    
    # Summary table
    from rich.table import Table
    
    table = Table(title="Test Results")
    table.add_column("Test Suite", style="cyan")
    table.add_column("Status", style="green")
    table.add_column("Duration", style="yellow")
    table.add_column("Notes", style="magenta")
    
    for result in results:
        status = "✅ PASS" if result["success"] else "❌ FAIL"
        duration = f"{result['duration']:.2f}s"
        notes = result.get("error", "OK" if result["success"] else f"Exit code: {result.get('return_code', 'N/A')}")
        
        table.add_row(
            result["description"],
            status,
            duration,
            notes
        )
    
    console.print(table)
    
    # Final summary
    success_rate = (passed / total * 100) if total > 0 else 0
    
    console.print(f"\n[bold]Total Test Suites: {total}[/bold]")
    console.print(f"[bold green]Passed: {passed}[/bold green]")
    console.print(f"[bold red]Failed: {total - passed}[/bold red]")
    console.print(f"[bold]Success Rate: {success_rate:.1f}%[/bold]")
    console.print(f"[bold]Total Duration: {total_duration:.2f} seconds[/bold]")
    
    # Save summary to file
    summary_file = Path("test_output/test_summary.txt")
    with open(summary_file, "w") as f:
        f.write(f"DARTSGPT Test Summary\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"=" * 50 + "\n\n")
        
        for result in results:
            f.write(f"{result['description']}: {'PASS' if result['success'] else 'FAIL'}\n")
            f.write(f"  Duration: {result['duration']:.2f}s\n")
            if not result['success']:
                f.write(f"  Error: {result.get('error', 'Unknown')}\n")
            f.write("\n")
        
        f.write(f"\nSummary:\n")
        f.write(f"Total: {total}, Passed: {passed}, Failed: {total - passed}\n")
        f.write(f"Success Rate: {success_rate:.1f}%\n")
        f.write(f"Total Duration: {total_duration:.2f}s\n")
    
    console.print(f"\n[dim]Test summary saved to: {summary_file}[/dim]")
    
    # Exit with appropriate code
    exit_code = 0 if passed == total else 1
    console.print(f"\n[{'green' if exit_code == 0 else 'red'}]Exiting with code: {exit_code}[/]")
    sys.exit(exit_code)


if __name__ == "__main__":
    main()