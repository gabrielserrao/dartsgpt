"""Agent evaluation framework."""

import json
from datetime import datetime
from typing import Dict, List, Optional, Tuple
from pydantic import BaseModel, Field
from langchain_core.messages import HumanMessage, AIMessage
from langchain_openai import ChatOpenAI

from dartsgpt.knowledge.rag_system import DARTSKnowledgeBase
from dartsgpt.evaluation.test_cases import TestCase, TestSuite


class EvaluationResult(BaseModel):
    """Result of evaluating an agent response."""
    
    test_case_id: str
    success: bool
    physics_match: bool = False
    template_match: bool = False
    features_match: float = 0.0  # Percentage of features matched
    parameters_match: float = 0.0  # Percentage of parameters matched
    keyword_coverage: float = 0.0  # Percentage of keywords found
    explanation: str = ""
    agent_response: Dict = Field(default_factory=dict)
    evaluation_time: datetime = Field(default_factory=datetime.now)
    

class AgentEvaluator:
    """Evaluates agent responses against test cases."""
    
    def __init__(self, knowledge_base: Optional[DARTSKnowledgeBase] = None):
        self.knowledge_base = knowledge_base
        self.llm = ChatOpenAI(temperature=0, model="gpt-4")
        self.evaluator_llm = ChatOpenAI(temperature=0, model="gpt-4")
        
    def evaluate_with_vector_store(self, agent_response: Dict, test_case: TestCase) -> EvaluationResult:
        """Evaluate using only vector store search (no ground truth in prompt)."""
        # Search for similar examples
        similar_docs = []
        if self.knowledge_base:
            similar_docs = self.knowledge_base.search(
                f"{test_case.expected_physics} {' '.join(test_case.expected_features)}",
                k=3
            )
        
        # Build evaluation prompt with similar examples
        eval_prompt = f"""
Evaluate if the agent's response correctly addresses the user prompt.

User Prompt: {test_case.prompt}

Agent Response:
- Physics Type: {agent_response.get('physics_type', 'Not specified')}
- Template: {agent_response.get('template', 'Not specified')}
- Features: {agent_response.get('features', [])}
- Parameters: {json.dumps(agent_response.get('parameters', {}), indent=2)}

Similar Examples from Knowledge Base:
{self._format_similar_docs(similar_docs)}

Evaluate based on:
1. Does the physics type make sense for the prompt?
2. Is the template appropriate?
3. Are the features correctly identified?
4. Are the parameters reasonable?

Provide a JSON response with:
{{
    "physics_correct": true/false,
    "template_correct": true/false,
    "features_score": 0-100,
    "parameters_score": 0-100,
    "explanation": "Brief explanation"
}}
"""
        
        eval_response = self.evaluator_llm.invoke([HumanMessage(content=eval_prompt)])
        eval_data = json.loads(eval_response.content)
        
        return self._create_evaluation_result(agent_response, test_case, eval_data, method="vector_store")
    
    def evaluate_with_ground_truth(self, agent_response: Dict, test_case: TestCase) -> EvaluationResult:
        """Evaluate with ground truth provided in the prompt."""
        eval_prompt = f"""
Evaluate if the agent's response correctly addresses the user prompt by comparing to the ground truth.

User Prompt: {test_case.prompt}

Ground Truth:
- Expected Physics Type: {test_case.expected_physics}
- Expected Template: {test_case.expected_template}
- Expected Features: {test_case.expected_features}
- Expected Parameters: {json.dumps(test_case.expected_parameters, indent=2)}
- Keywords: {test_case.keywords}

Agent Response:
- Physics Type: {agent_response.get('physics_type', 'Not specified')}
- Template: {agent_response.get('template', 'Not specified')}
- Features: {agent_response.get('features', [])}
- Parameters: {json.dumps(agent_response.get('parameters', {}), indent=2)}

Evaluate the match and provide a JSON response with:
{{
    "physics_correct": true/false,
    "template_correct": true/false,
    "features_score": 0-100 (percentage of expected features found),
    "parameters_score": 0-100 (percentage of expected parameters correctly identified),
    "keyword_coverage": 0-100 (percentage of keywords recognized),
    "explanation": "Brief explanation of evaluation"
}}
"""
        
        eval_response = self.evaluator_llm.invoke([HumanMessage(content=eval_prompt)])
        eval_data = json.loads(eval_response.content)
        
        return self._create_evaluation_result(agent_response, test_case, eval_data, method="ground_truth")
    
    def compare_evaluation_methods(self, agent_response: Dict, test_case: TestCase) -> Tuple[EvaluationResult, EvaluationResult]:
        """Compare both evaluation methods."""
        vector_eval = self.evaluate_with_vector_store(agent_response, test_case)
        ground_truth_eval = self.evaluate_with_ground_truth(agent_response, test_case)
        
        return vector_eval, ground_truth_eval
    
    def evaluate_test_suite(self, agent, test_suite: TestSuite, method: str = "both") -> Dict:
        """Evaluate an agent against an entire test suite."""
        results = {
            "suite_name": test_suite.name,
            "total_tests": len(test_suite.test_cases),
            "results": [],
            "summary": {}
        }
        
        for test_case in test_suite.test_cases:
            # Get agent response (this would be replaced with actual agent call)
            agent_response = self._mock_agent_response(test_case)
            
            if method == "both":
                vector_eval, ground_truth_eval = self.compare_evaluation_methods(agent_response, test_case)
                results["results"].append({
                    "test_case": test_case.dict(),
                    "vector_eval": vector_eval.dict(),
                    "ground_truth_eval": ground_truth_eval.dict()
                })
            elif method == "vector_store":
                eval_result = self.evaluate_with_vector_store(agent_response, test_case)
                results["results"].append({
                    "test_case": test_case.dict(),
                    "evaluation": eval_result.dict()
                })
            else:  # ground_truth
                eval_result = self.evaluate_with_ground_truth(agent_response, test_case)
                results["results"].append({
                    "test_case": test_case.dict(),
                    "evaluation": eval_result.dict()
                })
        
        # Calculate summary statistics
        results["summary"] = self._calculate_summary(results["results"])
        
        return results
    
    def _format_similar_docs(self, docs: List) -> str:
        """Format similar documents for the prompt."""
        if not docs:
            return "No similar examples found."
        
        formatted = []
        for i, doc in enumerate(docs[:3]):
            formatted.append(f"""
Example {i+1}:
- Model: {doc.metadata.get('model_name', 'Unknown')}
- Physics: {doc.metadata.get('physics_type', 'Unknown')}
- Features: {doc.metadata.get('features', 'None')}
""")
        return "\n".join(formatted)
    
    def _create_evaluation_result(self, agent_response: Dict, test_case: TestCase, eval_data: Dict, method: str) -> EvaluationResult:
        """Create an evaluation result from the evaluation data."""
        return EvaluationResult(
            test_case_id=test_case.id,
            success=eval_data.get("physics_correct", False) and eval_data.get("template_correct", False),
            physics_match=eval_data.get("physics_correct", False),
            template_match=eval_data.get("template_correct", False),
            features_match=eval_data.get("features_score", 0) / 100.0,
            parameters_match=eval_data.get("parameters_score", 0) / 100.0,
            keyword_coverage=eval_data.get("keyword_coverage", 0) / 100.0,
            explanation=f"[{method}] " + eval_data.get("explanation", ""),
            agent_response=agent_response
        )
    
    def _mock_agent_response(self, test_case: TestCase) -> Dict:
        """Mock agent response for testing (replace with actual agent call)."""
        # This simulates an agent response with some variations
        import random
        
        # Simulate some errors
        if random.random() < 0.2:  # 20% error rate
            return {
                "physics_type": random.choice(["compositional", "dead_oil", "geothermal"]),
                "template": "wrong_template",
                "features": random.sample(["thermal", "mechanical", "chemical"], k=1),
                "parameters": {"nx": 100, "ny": 100}
            }
        
        # Return mostly correct response
        return {
            "physics_type": test_case.expected_physics,
            "template": test_case.expected_template,
            "features": test_case.expected_features,
            "parameters": test_case.expected_parameters
        }
    
    def _calculate_summary(self, results: List[Dict]) -> Dict:
        """Calculate summary statistics from evaluation results."""
        total = len(results)
        if total == 0:
            return {}
        
        # Extract evaluation results
        evaluations = []
        for result in results:
            if "evaluation" in result:
                evaluations.append(result["evaluation"])
            else:
                # Handle comparison results
                if "vector_eval" in result:
                    evaluations.append(result["vector_eval"])
                if "ground_truth_eval" in result:
                    evaluations.append(result["ground_truth_eval"])
        
        if not evaluations:
            return {}
        
        return {
            "total_evaluations": len(evaluations),
            "success_rate": sum(1 for e in evaluations if e["success"]) / len(evaluations),
            "physics_accuracy": sum(1 for e in evaluations if e["physics_match"]) / len(evaluations),
            "template_accuracy": sum(1 for e in evaluations if e["template_match"]) / len(evaluations),
            "avg_features_match": sum(e["features_match"] for e in evaluations) / len(evaluations),
            "avg_parameters_match": sum(e["parameters_match"] for e in evaluations) / len(evaluations),
            "avg_keyword_coverage": sum(e.get("keyword_coverage", 0) for e in evaluations) / len(evaluations),
        }