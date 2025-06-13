"""Validation Agent for DARTSGPT.

This agent validates the generated DARTS code by:
- Checking syntax validity
- Verifying required methods exist
- Validating physics consistency
- Checking parameter ranges
- Ensuring import completeness
"""

from typing import Dict, Any, List, Optional, Tuple
import ast
import json
import re

from langchain.prompts import PromptTemplate
from langchain.tools import BaseTool, tool
from langchain.schema import HumanMessage, SystemMessage
from langchain_core.output_parsers import JsonOutputParser
from pydantic import BaseModel, Field

from ..agents.base_agent import BaseAgent, AgentInput, AgentOutput
from ..knowledge.templates.template_database import TEMPLATES


class ValidationIssue(BaseModel):
    """A validation issue found in the code."""
    severity: str = Field(description="Severity: error, warning, or info")
    category: str = Field(description="Category: syntax, logic, physics, or style")
    line: Optional[int] = Field(description="Line number if applicable")
    message: str = Field(description="Description of the issue")
    suggestion: Optional[str] = Field(description="Suggested fix")


class ValidationResult(BaseModel):
    """Validation result for generated code."""
    is_valid: bool = Field(description="Whether code is valid and runnable")
    syntax_valid: bool = Field(description="Whether code has valid Python syntax")
    physics_valid: bool = Field(description="Whether physics setup is consistent")
    imports_complete: bool = Field(description="Whether all imports are present")
    methods_complete: bool = Field(description="Whether all required methods exist")
    issues: List[ValidationIssue] = Field(description="List of validation issues")
    fixed_code: Optional[Dict[str, str]] = Field(description="Fixed code if auto-fixable")
    reasoning: str = Field(description="Explanation of validation")


class ValidatorAgent(BaseAgent):
    """Agent that validates generated DARTS code."""
    
    def __init__(self, llm=None):
        """Initialize the Validator Agent."""
        super().__init__(llm)
        self.output_parser = JsonOutputParser(pydantic_object=ValidationResult)
        
        # Required methods for DARTS models
        self.required_methods = [
            "set_physics",
            "set_reservoir", 
            "set_wells",
            "set_initial_conditions",
            "set_well_controls"
        ]
        
        # Common DARTS imports
        self.common_imports = {
            "numpy": "import numpy as np",
            "darts.models": "from darts.models.cicd_model import CICDModel",
            "darts.physics.compositional": "from darts.physics.compositional import Compositional",
            "darts.physics.dead_oil": "from darts.physics.dead_oil import DeadOil",
            "darts.physics.geothermal": "from darts.physics.geothermal import Geothermal",
            "darts.physics.black_oil": "from darts.physics.black_oil import BlackOil"
        }
    
    def get_tools(self) -> List[BaseTool]:
        """Get tools for code validation."""
        
        @tool
        def check_python_syntax(code: str) -> Tuple[bool, List[str]]:
            """Check if code has valid Python syntax.
            
            Args:
                code: Python code to validate
                
            Returns:
                Tuple of (is_valid, error_messages)
            """
            try:
                ast.parse(code)
                return True, []
            except SyntaxError as e:
                errors = [f"Syntax error at line {e.lineno}: {e.msg}"]
                return False, errors
            except Exception as e:
                errors = [f"Parse error: {str(e)}"]
                return False, errors
        
        @tool
        def check_required_methods(code: str, required: List[str]) -> Tuple[bool, List[str]]:
            """Check if all required methods are present.
            
            Args:
                code: Python code to check
                required: List of required method names
                
            Returns:
                Tuple of (all_present, missing_methods)
            """
            try:
                tree = ast.parse(code)
                defined_methods = set()
                
                for node in ast.walk(tree):
                    if isinstance(node, ast.ClassDef):
                        for item in node.body:
                            if isinstance(item, ast.FunctionDef):
                                defined_methods.add(item.name)
                
                missing = [m for m in required if m not in defined_methods]
                return len(missing) == 0, missing
            except:
                return False, required
        
        @tool
        def check_physics_consistency(code: str, physics_type: str) -> List[ValidationIssue]:
            """Check physics setup consistency.
            
            Args:
                code: Model code
                physics_type: Expected physics type
                
            Returns:
                List of validation issues
            """
            issues = []
            
            # Check physics imports
            physics_imports = {
                "compositional": ["Compositional", "ConstantK", "Flash"],
                "dead_oil": ["DeadOil", "WaterOil"],
                "geothermal": ["Geothermal", "IAPWS", "Enthalpy"],
                "black_oil": ["BlackOil", "PVT", "ThreePhase"]
            }
            
            expected_imports = physics_imports.get(physics_type, [])
            for imp in expected_imports:
                if imp not in code and imp.lower() not in code.lower():
                    issues.append(ValidationIssue(
                        severity="warning",
                        category="physics",
                        line=None,
                        message=f"Expected {imp} for {physics_type} physics",
                        suggestion=f"Add import for {imp}"
                    ))
            
            # Check phase definitions
            if physics_type == "compositional" and "phases" not in code:
                issues.append(ValidationIssue(
                    severity="error",
                    category="physics",
                    line=None,
                    message="Compositional model missing phase definitions",
                    suggestion="Add self.phases = ['gas', 'oil']"
                ))
            
            # Check thermal properties for geothermal
            if physics_type == "geothermal":
                if "thermal" not in code and "temperature" not in code:
                    issues.append(ValidationIssue(
                        severity="error",
                        category="physics",
                        line=None,
                        message="Geothermal model missing thermal properties",
                        suggestion="Add thermal conductivity and heat capacity"
                    ))
            
            return issues
        
        @tool
        def check_parameter_ranges(code: str) -> List[ValidationIssue]:
            """Check if parameters are within reasonable ranges.
            
            Args:
                code: Model code
                
            Returns:
                List of validation issues
            """
            issues = []
            
            # Extract numerical values with regex
            patterns = {
                "porosity": (r'porosity\s*=\s*([0-9.]+)', 0.0, 1.0),
                "nx": (r'nx\s*=\s*(\d+)', 1, 1000),
                "ny": (r'ny\s*=\s*(\d+)', 1, 1000),
                "nz": (r'nz\s*=\s*(\d+)', 1, 100),
                "time_step": (r'time_step\s*=\s*([0-9.]+)', 0.001, 100),
            }
            
            for param, (pattern, min_val, max_val) in patterns.items():
                matches = re.findall(pattern, code)
                for match in matches:
                    value = float(match)
                    if value < min_val or value > max_val:
                        issues.append(ValidationIssue(
                            severity="warning",
                            category="logic",
                            line=None,
                            message=f"{param} value {value} outside range [{min_val}, {max_val}]",
                            suggestion=f"Check if {param}={value} is intended"
                        ))
            
            return issues
        
        @tool
        def check_imports(code: str, physics_type: str) -> Tuple[bool, List[str]]:
            """Check if all necessary imports are present.
            
            Args:
                code: Model code
                physics_type: Physics type
                
            Returns:
                Tuple of (complete, missing_imports)
            """
            missing = []
            
            # Always need numpy
            if "import numpy" not in code and "from numpy" not in code:
                missing.append("numpy")
            
            # Need base model
            if "CICDModel" not in code:
                missing.append("CICDModel")
            
            # Physics-specific imports
            physics_map = {
                "compositional": "compositional",
                "dead_oil": "dead_oil", 
                "geothermal": "geothermal",
                "black_oil": "black_oil"
            }
            
            physics_module = physics_map.get(physics_type)
            if physics_module and f"darts.physics.{physics_module}" not in code:
                missing.append(f"darts.physics.{physics_module}")
            
            return len(missing) == 0, missing
        
        @tool
        def suggest_fixes(issues: List[ValidationIssue]) -> Dict[str, str]:
            """Suggest fixes for common issues.
            
            Args:
                issues: List of validation issues
                
            Returns:
                Dictionary of fix suggestions
            """
            fixes = {}
            
            for issue in issues:
                if issue.category == "syntax" and "import" in issue.message.lower():
                    fixes["imports"] = "Add missing imports at the top of the file"
                elif issue.category == "physics" and "phase" in issue.message.lower():
                    fixes["phases"] = "Define phases in set_physics method"
                elif issue.category == "logic" and "method" in issue.message.lower():
                    fixes["methods"] = "Implement missing required methods"
            
            return fixes
        
        return [
            check_python_syntax,
            check_required_methods,
            check_physics_consistency,
            check_parameter_ranges,
            check_imports,
            suggest_fixes
        ]
    
    def get_prompt_template(self) -> PromptTemplate:
        """Get the prompt template for validation."""
        template = """You are an expert DARTS code validator.
        
Validate the generated DARTS model code for correctness and completeness.

Generated Code:
{generated_code}

Template Info: {template_info}

Validation Results:
- Syntax Check: {syntax_result}
- Required Methods: {methods_result}
- Physics Issues: {physics_issues}
- Parameter Issues: {param_issues}
- Import Check: {import_result}

Instructions:
1. Determine if the code is valid and runnable
2. Identify all issues by severity (error, warning, info)
3. Suggest fixes for any issues found
4. If code has minor fixable issues, provide fixed version
5. Explain the validation results clearly

Consider:
- Syntax must be valid Python
- All required DARTS methods must be present
- Physics setup must match the selected physics type
- Parameters should be within reasonable ranges
- All necessary imports must be included

{format_instructions}
"""
        
        return PromptTemplate(
            template=template,
            input_variables=[
                "generated_code",
                "template_info",
                "syntax_result",
                "methods_result", 
                "physics_issues",
                "param_issues",
                "import_result"
            ],
            partial_variables={
                "format_instructions": self.output_parser.get_format_instructions()
            }
        )
    
    async def process(self, input_data: AgentInput) -> AgentOutput:
        """Validate the generated DARTS code.
        
        Args:
            input_data: Input containing generated code
            
        Returns:
            AgentOutput with validation results
        """
        if not self.validate_input(input_data):
            return AgentOutput(
                result={"error": "Invalid input"},
                reasoning="Empty or invalid input provided",
                confidence=0.0
            )
        
        # Get generated code from context
        code_gen = input_data.context.get("code_generation", {})
        model_code = code_gen.get("model_code", "")
        template_info = input_data.context.get("template_selection", {})
        template_name = template_info.get("template_name", "2ph_do")
        physics_type = TEMPLATES.get(template_name, TEMPLATES["2ph_do"]).physics_type
        
        # Use tools to validate
        tools = self.get_tools()
        syntax_tool = tools[0]
        methods_tool = tools[1]
        physics_tool = tools[2]
        params_tool = tools[3]
        imports_tool = tools[4]
        fixes_tool = tools[5]
        
        # Run validation checks
        syntax_valid, syntax_errors = syntax_tool.invoke({"code": model_code})
        methods_valid, missing_methods = methods_tool.invoke({
            "code": model_code,
            "required": self.required_methods
        })
        physics_issues = physics_tool.invoke({
            "code": model_code,
            "physics_type": physics_type
        })
        param_issues = params_tool.invoke({"code": model_code})
        imports_valid, missing_imports = imports_tool.invoke({
            "code": model_code,
            "physics_type": physics_type
        })
        
        # Collect all issues
        all_issues = []
        
        # Syntax issues
        for error in syntax_errors:
            all_issues.append(ValidationIssue(
                severity="error",
                category="syntax",
                line=None,
                message=error,
                suggestion="Fix syntax error"
            ))
        
        # Method issues
        for method in missing_methods:
            all_issues.append(ValidationIssue(
                severity="error",
                category="logic",
                line=None,
                message=f"Missing required method: {method}",
                suggestion=f"Implement {method} method"
            ))
        
        # Physics issues
        all_issues.extend([issue.model_dump() for issue in physics_issues])
        
        # Parameter issues
        all_issues.extend([issue.model_dump() for issue in param_issues])
        
        # Import issues
        for imp in missing_imports:
            all_issues.append(ValidationIssue(
                severity="error",
                category="syntax",
                line=1,
                message=f"Missing import: {imp}",
                suggestion=f"Add import for {imp}"
            ))
        
        # Get fix suggestions
        fixes = fixes_tool.invoke({"issues": all_issues})
        
        # Prepare prompt
        prompt = self.get_prompt_template()
        formatted_prompt = prompt.format(
            generated_code=model_code,
            template_info=json.dumps(template_info, indent=2),
            syntax_result=json.dumps({"valid": syntax_valid, "errors": syntax_errors}),
            methods_result=json.dumps({"valid": methods_valid, "missing": missing_methods}),
            physics_issues=json.dumps([i.model_dump() for i in physics_issues], indent=2),
            param_issues=json.dumps([i.model_dump() for i in param_issues], indent=2),
            import_result=json.dumps({"valid": imports_valid, "missing": missing_imports})
        )
        
        # Get LLM response
        messages = [
            SystemMessage(content="You are an expert DARTS code validator."),
            HumanMessage(content=formatted_prompt)
        ]
        
        response = await self.llm.ainvoke(messages)
        
        # Parse response
        try:
            validation = self.output_parser.parse(response.content)
            
            # Update context with validation results
            updated_context = input_data.context.copy()
            updated_context["validation"] = validation.model_dump()
            
            # Determine confidence based on validation
            if validation.is_valid:
                confidence = 0.95
            elif len([i for i in validation.issues if i["severity"] == "error"]) == 0:
                confidence = 0.8
            else:
                confidence = 0.5
            
            return AgentOutput(
                result=validation.model_dump(),
                reasoning=validation.reasoning,
                confidence=confidence,
                next_agent=None  # End of pipeline
            )
        except Exception as e:
            # Fallback validation result
            error_count = len(syntax_errors) + len(missing_methods) + len(missing_imports)
            is_valid = error_count == 0
            
            fallback_validation = {
                "is_valid": is_valid,
                "syntax_valid": syntax_valid,
                "physics_valid": len(physics_issues) == 0,
                "imports_complete": imports_valid,
                "methods_complete": methods_valid,
                "issues": [i.model_dump() if hasattr(i, 'model_dump') else i for i in all_issues],
                "fixed_code": None,
                "reasoning": f"Validation {'passed' if is_valid else 'failed'} with {error_count} errors"
            }
            
            updated_context = input_data.context.copy()
            updated_context["validation"] = fallback_validation
            
            return AgentOutput(
                result=fallback_validation,
                reasoning=fallback_validation["reasoning"],
                confidence=0.9 if is_valid else 0.3,
                next_agent=None
            )