"""Utilities for executing DARTS models."""

import os
import sys
import subprocess
import threading
import queue
from pathlib import Path
from typing import Dict, Any, Optional, Callable
import time
import psutil
import json


class DARTSExecutor:
    """Utility class for executing DARTS models with monitoring."""
    
    def __init__(self, timeout: int = 300, memory_limit_gb: float = 4.0):
        """
        Initialize DARTS executor.
        
        Args:
            timeout: Maximum execution time in seconds (default: 5 minutes)
            memory_limit_gb: Maximum memory usage in GB (default: 4 GB)
        """
        self.timeout = timeout
        self.memory_limit_gb = memory_limit_gb
        self.process = None
        self.monitor_thread = None
        
    def execute(
        self, 
        model_path: Path, 
        callback: Optional[Callable[[str], None]] = None
    ) -> Dict[str, Any]:
        """
        Execute a DARTS model with monitoring.
        
        Args:
            model_path: Path to the model directory
            callback: Optional callback for output streaming
            
        Returns:
            Execution result dictionary
        """
        main_file = model_path / "main.py"
        if not main_file.exists():
            return {
                "success": False,
                "error": f"Main file not found: {main_file}",
                "execution_time": 0
            }
        
        # Prepare execution
        start_time = time.time()
        output_lines = []
        error_lines = []
        
        # Setup environment
        env = os.environ.copy()
        env["PYTHONPATH"] = str(model_path) + ":" + env.get("PYTHONPATH", "")
        
        # Command to run
        cmd = [sys.executable, str(main_file)]
        
        try:
            # Start process
            self.process = subprocess.Popen(
                cmd,
                cwd=str(model_path),
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1
            )
            
            # Start monitoring thread
            monitor_queue = queue.Queue()
            self.monitor_thread = threading.Thread(
                target=self._monitor_process,
                args=(self.process.pid, monitor_queue)
            )
            self.monitor_thread.start()
            
            # Read output with timeout
            import select
            
            while self.process.poll() is None:
                # Check timeout
                if time.time() - start_time > self.timeout:
                    self.process.terminate()
                    raise TimeoutError(f"Execution exceeded {self.timeout} seconds")
                
                # Check memory from monitor
                try:
                    memory_gb = monitor_queue.get_nowait()
                    if memory_gb > self.memory_limit_gb:
                        self.process.terminate()
                        raise MemoryError(f"Memory usage exceeded {self.memory_limit_gb} GB")
                except queue.Empty:
                    pass
                
                # Read available output
                if sys.platform != 'win32':
                    # Unix-like systems
                    readable, _, _ = select.select([self.process.stdout, self.process.stderr], [], [], 0.1)
                    
                    if self.process.stdout in readable:
                        line = self.process.stdout.readline()
                        if line:
                            output_lines.append(line)
                            if callback:
                                callback(line.rstrip())
                    
                    if self.process.stderr in readable:
                        line = self.process.stderr.readline()
                        if line:
                            error_lines.append(line)
                else:
                    # Windows - simpler approach
                    line = self.process.stdout.readline()
                    if line:
                        output_lines.append(line)
                        if callback:
                            callback(line.rstrip())
            
            # Get remaining output
            remaining_out, remaining_err = self.process.communicate()
            if remaining_out:
                output_lines.append(remaining_out)
            if remaining_err:
                error_lines.append(remaining_err)
            
            # Stop monitoring
            if self.monitor_thread:
                self.monitor_thread.join(timeout=1)
            
            execution_time = time.time() - start_time
            
            return {
                "success": self.process.returncode == 0,
                "return_code": self.process.returncode,
                "output": "".join(output_lines),
                "error": "".join(error_lines),
                "execution_time": execution_time,
                "memory_peak_gb": self._get_peak_memory(monitor_queue)
            }
            
        except Exception as e:
            return {
                "success": False,
                "error": str(e),
                "output": "".join(output_lines),
                "execution_time": time.time() - start_time
            }
        finally:
            # Cleanup
            if self.process and self.process.poll() is None:
                self.process.terminate()
                self.process.wait()
    
    def _monitor_process(self, pid: int, memory_queue: queue.Queue):
        """Monitor process memory usage."""
        try:
            process = psutil.Process(pid)
            while process.is_running():
                memory_gb = process.memory_info().rss / (1024 ** 3)
                memory_queue.put(memory_gb)
                time.sleep(0.5)
        except psutil.NoSuchProcess:
            pass
    
    def _get_peak_memory(self, memory_queue: queue.Queue) -> float:
        """Get peak memory usage from monitoring queue."""
        peak = 0.0
        while not memory_queue.empty():
            try:
                memory = memory_queue.get_nowait()
                peak = max(peak, memory)
            except queue.Empty:
                break
        return peak


def check_darts_installation() -> bool:
    """Check if DARTS is properly installed."""
    try:
        import darts
        return True
    except ImportError:
        return False


def parse_darts_output(output: str) -> Dict[str, Any]:
    """Parse DARTS simulation output for key metrics."""
    
    metrics = {
        "converged": False,
        "timesteps": 0,
        "iterations": 0,
        "final_time": None,
        "cpu_time": None,
        "warnings": [],
        "performance": {}
    }
    
    lines = output.splitlines()
    for i, line in enumerate(lines):
        line_lower = line.lower()
        
        # Check convergence
        if "converged" in line_lower and "not" not in line_lower:
            metrics["converged"] = True
        
        # Extract iterations
        if "iterations" in line_lower:
            try:
                import re
                numbers = re.findall(r'\d+', line)
                if numbers:
                    metrics["iterations"] = int(numbers[0])
            except:
                pass
        
        # Extract timesteps
        if "time step" in line_lower or "timestep" in line_lower:
            try:
                metrics["timesteps"] += 1
            except:
                pass
        
        # Extract final time
        if "final time" in line_lower:
            try:
                import re
                numbers = re.findall(r'\d+\.?\d*', line)
                if numbers:
                    metrics["final_time"] = float(numbers[-1])
            except:
                pass
        
        # Extract CPU time
        if "cpu time" in line_lower or "elapsed time" in line_lower:
            try:
                import re
                numbers = re.findall(r'\d+\.?\d*', line)
                if numbers:
                    metrics["cpu_time"] = float(numbers[0])
            except:
                pass
        
        # Collect warnings
        if "warning" in line_lower:
            metrics["warnings"].append(line.strip())
        
        # Performance metrics
        if "performance" in line_lower or "statistics" in line_lower:
            # Look for key performance indicators in next few lines
            for j in range(i, min(i + 10, len(lines))):
                perf_line = lines[j]
                if "linear solver" in perf_line.lower():
                    metrics["performance"]["linear_solver"] = perf_line.strip()
                elif "newton" in perf_line.lower():
                    metrics["performance"]["newton"] = perf_line.strip()
    
    return metrics


def format_execution_report(result: Dict[str, Any]) -> str:
    """Format execution result into a readable report."""
    
    report = []
    report.append("=" * 60)
    report.append("DARTS Model Execution Report")
    report.append("=" * 60)
    
    # Basic info
    report.append(f"Status: {'✅ SUCCESS' if result.get('success') else '❌ FAILED'}")
    report.append(f"Execution Time: {result.get('execution_time', 0):.2f} seconds")
    
    if result.get('memory_peak_gb'):
        report.append(f"Peak Memory: {result['memory_peak_gb']:.2f} GB")
    
    # Parse output for metrics
    if result.get('output'):
        metrics = parse_darts_output(result['output'])
        
        report.append("\nSimulation Metrics:")
        report.append(f"  - Converged: {'Yes' if metrics['converged'] else 'No'}")
        report.append(f"  - Time Steps: {metrics['timesteps']}")
        report.append(f"  - Iterations: {metrics['iterations']}")
        
        if metrics['final_time']:
            report.append(f"  - Final Time: {metrics['final_time']} days")
        
        if metrics['cpu_time']:
            report.append(f"  - CPU Time: {metrics['cpu_time']:.2f} seconds")
        
        if metrics['warnings']:
            report.append(f"\nWarnings ({len(metrics['warnings'])}):")
            for warning in metrics['warnings'][:5]:  # Show first 5 warnings
                report.append(f"  - {warning}")
            if len(metrics['warnings']) > 5:
                report.append(f"  ... and {len(metrics['warnings']) - 5} more")
    
    # Errors
    if result.get('error'):
        report.append("\nError Details:")
        for line in result['error'].splitlines()[:10]:  # First 10 lines
            report.append(f"  {line}")
    
    # Results location
    if result.get('results_path'):
        report.append(f"\nResults saved to: {result['results_path']}")
        if result.get('results_files'):
            report.append(f"Generated files: {', '.join(result['results_files'][:5])}")
    
    report.append("=" * 60)
    
    return "\n".join(report)