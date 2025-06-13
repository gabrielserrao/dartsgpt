"""Execution agent for running generated DARTS models."""

import subprocess
import sys
import os
from pathlib import Path
from typing import Dict, Any, Optional, Tuple
import json
import time
from datetime import datetime

from .base_agent import BaseAgent, AgentInput, AgentOutput


class ExecutorAgent(BaseAgent):
    """Agent responsible for executing generated DARTS models."""
    
    def __init__(self):
        super().__init__(
            name="Executor",
            description="Executes generated DARTS models and collects results"
        )
        
    def process(self, input_data: AgentInput) -> AgentOutput:
        """Execute the generated DARTS model."""
        
        # Get the output directory where model files are saved
        output_dir = input_data.get("output_dir", "output")
        model_path = Path(output_dir)
        
        # Check if required files exist
        model_file = model_path / "model.py"
        main_file = model_path / "main.py"
        
        if not model_file.exists() or not main_file.exists():
            return AgentOutput(
                success=False,
                data={
                    "error": "Model files not found",
                    "details": f"Expected files at {model_path}"
                },
                metadata={"agent": self.name}
            )
        
        # Execute the model
        execution_result = self._execute_model(model_path, main_file)
        
        return AgentOutput(
            success=execution_result["success"],
            data=execution_result,
            metadata={
                "agent": self.name,
                "execution_time": execution_result.get("execution_time", 0),
                "output_dir": str(model_path)
            }
        )
    
    def _execute_model(self, model_path: Path, main_file: Path) -> Dict[str, Any]:
        """Execute the DARTS model using subprocess."""
        
        start_time = time.time()
        result = {
            "success": False,
            "output": "",
            "error": "",
            "return_code": None,
            "execution_time": 0,
            "timestamp": datetime.now().isoformat()
        }
        
        try:
            # Prepare environment
            env = os.environ.copy()
            env["PYTHONPATH"] = str(model_path) + ":" + env.get("PYTHONPATH", "")
            
            # Run the model
            print(f"🚀 Executing DARTS model at {main_file}...")
            
            # Use uv run to execute with proper dependencies
            cmd = ["uv", "run", "python", str(main_file)]
            
            # Execute with real-time output streaming
            process = subprocess.Popen(
                cmd,
                cwd=str(model_path),
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                universal_newlines=True
            )
            
            # Collect output
            stdout_lines = []
            stderr_lines = []
            
            # Stream output in real-time
            while True:
                # Check if process is still running
                poll = process.poll()
                
                # Read available output
                if process.stdout:
                    line = process.stdout.readline()
                    if line:
                        stdout_lines.append(line)
                        print(f"  {line.rstrip()}")
                
                # Check if process has finished
                if poll is not None:
                    # Read any remaining output
                    if process.stdout:
                        remaining_stdout = process.stdout.read()
                        if remaining_stdout:
                            stdout_lines.append(remaining_stdout)
                            for line in remaining_stdout.splitlines():
                                print(f"  {line}")
                    
                    if process.stderr:
                        stderr_output = process.stderr.read()
                        if stderr_output:
                            stderr_lines.append(stderr_output)
                    
                    break
            
            # Get return code
            return_code = process.returncode
            
            # Prepare result
            result["output"] = "".join(stdout_lines)
            result["error"] = "".join(stderr_lines)
            result["return_code"] = return_code
            result["success"] = (return_code == 0)
            result["execution_time"] = time.time() - start_time
            
            # Check for results directory
            results_dir = model_path / "results"
            if results_dir.exists():
                result["results_path"] = str(results_dir)
                result["results_files"] = [f.name for f in results_dir.iterdir()]
            
            # Parse execution summary if available
            result["summary"] = self._parse_execution_summary(result["output"])
            
        except subprocess.TimeoutExpired:
            result["error"] = "Execution timed out after 5 minutes"
            result["success"] = False
            
        except Exception as e:
            result["error"] = f"Execution failed: {str(e)}"
            result["success"] = False
            
        return result
    
    def _parse_execution_summary(self, output: str) -> Dict[str, Any]:
        """Parse execution output to extract summary information."""
        
        summary = {
            "simulation_completed": False,
            "convergence_achieved": False,
            "time_steps": 0,
            "final_time": None,
            "warnings": [],
            "errors": []
        }
        
        lines = output.splitlines()
        for line in lines:
            line_lower = line.lower()
            
            # Check completion
            if "simulation completed" in line_lower or "results saved" in line_lower:
                summary["simulation_completed"] = True
            
            # Check convergence
            if "converged" in line_lower:
                summary["convergence_achieved"] = True
            
            # Extract simulation time
            if "final time" in line_lower or "simulation time" in line_lower:
                try:
                    # Try to extract number from line
                    import re
                    numbers = re.findall(r'\d+\.?\d*', line)
                    if numbers:
                        summary["final_time"] = float(numbers[-1])
                except:
                    pass
            
            # Collect warnings and errors
            if "warning" in line_lower:
                summary["warnings"].append(line.strip())
            elif "error" in line_lower:
                summary["errors"].append(line.strip())
        
        return summary
    
    def validate_input(self, input_data: AgentInput) -> Tuple[bool, Optional[str]]:
        """Validate that model files exist before execution."""
        
        output_dir = input_data.get("output_dir", "output")
        model_path = Path(output_dir)
        
        if not model_path.exists():
            return False, f"Output directory does not exist: {model_path}"
        
        model_file = model_path / "model.py"
        main_file = model_path / "main.py"
        
        if not model_file.exists():
            return False, f"Model file not found: {model_file}"
            
        if not main_file.exists():
            return False, f"Main file not found: {main_file}"
        
        return True, None