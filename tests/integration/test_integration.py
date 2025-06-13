#!/usr/bin/env python
"""Integration tests for DARTSGPT system."""

import json
import asyncio
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich.progress import track

from orchestrator import run_dartsgpt_orchestrator
from main import cli
from click.testing import CliRunner

console = Console()


class IntegrationTestSuite:
    """Integration tests for the complete system."""
    
    def __init__(self):
        self.results = []
        self.runner = CliRunner()
    
    def test_cli_commands(self):
        """Test all CLI commands."""
        console.print("\n[bold blue]Testing CLI Commands[/bold blue]")
        
        commands = [
            (["check-system"], "System check"),
            (["list-templates"], "List templates"),
            (["show-examples"], "Show examples"),
            (["show-examples", "--physics", "geothermal"], "Show geothermal examples"),
        ]
        
        for cmd, description in commands:
            result = self.runner.invoke(cli, cmd)
            status = "✅" if result.exit_code == 0 else "❌"
            console.print(f"{description}: {status}")
            
            self.results.append({
                "test": f"CLI: {' '.join(cmd)}",
                "passed": result.exit_code == 0,
                "output": result.output[:100] if result.output else "No output"
            })
    
    def test_generation_workflow(self):
        """Test complete generation workflow."""
        console.print("\n[bold blue]Testing Generation Workflow[/bold blue]")
        
        test_cases = [
            {
                "name": "Simple CO2 injection",
                "prompt": "Create a CO2 injection model",
                "expected_physics": "compositional",
                "check_for": ["CO2", "injection", "compositional"]
            },
            {
                "name": "Waterflood with parameters",
                "prompt": "Build a waterflood model in 100x100x10 reservoir with 25% porosity",
                "expected_physics": "dead_oil",
                "check_for": ["water", "100", "0.25"]
            },
            {
                "name": "Geothermal with temperature",
                "prompt": "Model a geothermal system at 250°C",
                "expected_physics": "geothermal",
                "check_for": ["geothermal", "250", "temperature"]
            }
        ]
        
        for test in test_cases:
            console.print(f"\n[yellow]Testing: {test['name']}[/yellow]")
            
            # Test via CLI
            result = self.runner.invoke(cli, [
                "generate", test["prompt"], 
                "-o", f"test_output/integration_{test['name'].replace(' ', '_')}",
                "--dry-run"
            ])
            
            output = result.output
            success = result.exit_code == 0
            
            # Check for expected content
            checks_passed = []
            for keyword in test["check_for"]:
                if keyword.lower() in output.lower():
                    checks_passed.append(keyword)
            
            console.print(f"  Exit code: {result.exit_code}")
            console.print(f"  Keywords found: {checks_passed}")
            
            self.results.append({
                "test": test["name"],
                "passed": success and len(checks_passed) > 0,
                "details": {
                    "exit_code": result.exit_code,
                    "keywords_found": checks_passed,
                    "expected_keywords": test["check_for"]
                }
            })
    
    def test_error_handling(self):
        """Test error handling and edge cases."""
        console.print("\n[bold blue]Testing Error Handling[/bold blue]")
        
        error_cases = [
            {
                "name": "Empty prompt",
                "args": ["generate", ""],
                "should_fail": True
            },
            {
                "name": "Invalid output directory",
                "args": ["generate", "test", "-o", "/invalid/path/that/does/not/exist"],
                "should_fail": False  # Should create directory
            },
            {
                "name": "Very long prompt",
                "args": ["generate", "x" * 1000],
                "should_fail": False
            }
        ]
        
        for test in error_cases:
            console.print(f"\n[yellow]Testing: {test['name']}[/yellow]")
            
            result = self.runner.invoke(cli, test["args"])
            
            if test["should_fail"]:
                passed = result.exit_code != 0
                console.print(f"  Expected failure: {'✅' if passed else '❌'}")
            else:
                passed = result.exit_code == 0
                console.print(f"  Expected success: {'✅' if passed else '❌'}")
            
            self.results.append({
                "test": f"Error handling: {test['name']}",
                "passed": passed,
                "exit_code": result.exit_code
            })
    
    def test_file_generation(self):
        """Test actual file generation."""
        console.print("\n[bold blue]Testing File Generation[/bold blue]")
        
        # Generate actual files
        test_prompt = "Create a simple dead oil reservoir model"
        output_dir = Path("test_output/integration_files")
        
        result = self.runner.invoke(cli, [
            "generate", test_prompt,
            "-o", str(output_dir)
        ])
        
        if result.exit_code == 0:
            # Check generated files
            expected_files = ["model.py", "main.py", "metadata.json"]
            files_exist = []
            
            for filename in expected_files:
                filepath = output_dir / filename
                if filepath.exists():
                    files_exist.append(filename)
                    console.print(f"  ✅ {filename} created")
                    
                    # Check file content
                    if filename == "model.py":
                        content = filepath.read_text()
                        has_comments = "Original prompt:" in content
                        has_imports = "from darts.models" in content
                        console.print(f"    - Has comments: {'✅' if has_comments else '❌'}")
                        console.print(f"    - Has imports: {'✅' if has_imports else '❌'}")
                else:
                    console.print(f"  ❌ {filename} missing")
            
            self.results.append({
                "test": "File generation",
                "passed": len(files_exist) == len(expected_files),
                "files_created": files_exist
            })
        else:
            console.print(f"  ❌ Generation failed: {result.exit_code}")
            self.results.append({
                "test": "File generation",
                "passed": False,
                "error": "Generation failed"
            })
    
    def test_dry_run_mode(self):
        """Test dry run mode doesn't create files."""
        console.print("\n[bold blue]Testing Dry Run Mode[/bold blue]")
        
        output_dir = Path("test_output/dry_run_test")
        
        # Ensure directory doesn't exist
        if output_dir.exists():
            import shutil
            shutil.rmtree(output_dir)
        
        result = self.runner.invoke(cli, [
            "generate", "Test model",
            "-o", str(output_dir),
            "--dry-run"
        ])
        
        # Check that no files were created
        files_created = output_dir.exists() and any(output_dir.iterdir())
        
        self.results.append({
            "test": "Dry run mode",
            "passed": result.exit_code == 0 and not files_created,
            "files_created": files_created
        })
        
        console.print(f"  Exit code: {result.exit_code}")
        console.print(f"  Files created: {'❌ Yes' if files_created else '✅ No'}")
    
    def test_verbose_mode(self):
        """Test verbose output mode."""
        console.print("\n[bold blue]Testing Verbose Mode[/bold blue]")
        
        result = self.runner.invoke(cli, [
            "generate", "Create a model",
            "--verbose",
            "--dry-run"
        ])
        
        output = result.output
        
        # Check for agent outputs in verbose mode
        verbose_indicators = [
            "Intent classified:",
            "Template selected:",
            "Parameters extracted:",
            "Code generated:",
            "Validation"
        ]
        
        found_indicators = [ind for ind in verbose_indicators if ind in output]
        
        self.results.append({
            "test": "Verbose mode",
            "passed": len(found_indicators) >= 3,
            "indicators_found": found_indicators
        })
        
        console.print(f"  Verbose indicators found: {len(found_indicators)}/{len(verbose_indicators)}")
    
    def test_demo_command(self):
        """Test demo command."""
        console.print("\n[bold blue]Testing Demo Command[/bold blue]")
        
        # Test non-interactive demo (should run quickly)
        result = self.runner.invoke(cli, ["demo"], input="n\n")
        
        self.results.append({
            "test": "Demo command",
            "passed": result.exit_code == 0,
            "output_length": len(result.output) if result.output else 0
        })
        
        console.print(f"  Exit code: {result.exit_code}")
    
    def run_all_tests(self):
        """Run all integration tests."""
        console.print("[bold cyan]DARTSGPT Integration Test Suite[/bold cyan]")
        console.print("=" * 60)
        
        # Create test output directory
        Path("test_output").mkdir(exist_ok=True)
        
        # Run test suites
        test_methods = [
            self.test_cli_commands,
            self.test_generation_workflow,
            self.test_error_handling,
            self.test_file_generation,
            self.test_dry_run_mode,
            self.test_verbose_mode,
            self.test_demo_command
        ]
        
        for test_method in track(test_methods, description="Running tests..."):
            try:
                test_method()
            except Exception as e:
                console.print(f"[red]Test failed: {str(e)}[/red]")
                self.results.append({
                    "test": test_method.__name__,
                    "passed": False,
                    "error": str(e)
                })
        
        # Display results
        self.display_results()
    
    def display_results(self):
        """Display test results in a table."""
        console.print("\n[bold cyan]Test Results Summary[/bold cyan]")
        
        table = Table(title="Integration Test Results")
        table.add_column("Test", style="cyan")
        table.add_column("Status", style="green")
        table.add_column("Details", style="yellow")
        
        passed = 0
        total = len(self.results)
        
        for result in self.results:
            status = "✅ PASS" if result["passed"] else "❌ FAIL"
            if result["passed"]:
                passed += 1
            
            details = ""
            if "error" in result:
                details = f"Error: {result['error'][:50]}"
            elif "details" in result:
                details = str(result["details"])[:50]
            elif "output" in result:
                details = result["output"][:50]
            
            table.add_row(result["test"], status, details)
        
        console.print(table)
        
        # Summary
        success_rate = (passed / total * 100) if total > 0 else 0
        console.print(f"\n[bold]Total: {total} | Passed: {passed} | Failed: {total - passed}[/bold]")
        console.print(f"[bold]Success Rate: {success_rate:.1f}%[/bold]")
        
        # Save results to file
        results_file = Path("test_output/integration_results.json")
        with open(results_file, "w") as f:
            json.dump(self.results, f, indent=2)
        console.print(f"\nDetailed results saved to: {results_file}")


def main():
    """Run integration tests."""
    suite = IntegrationTestSuite()
    suite.run_all_tests()


if __name__ == "__main__":
    main()