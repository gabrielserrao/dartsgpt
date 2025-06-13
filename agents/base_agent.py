"""Base agent class for DARTSGPT agents."""

from abc import ABC, abstractmethod
from typing import Dict, Any, List, Optional
from pydantic import BaseModel

from langchain.schema import BaseMessage
from langchain.tools import BaseTool
from langchain.prompts import PromptTemplate
from langchain_openai import ChatOpenAI
from langchain.agents import AgentExecutor


class AgentInput(BaseModel):
    """Input for agent processing."""
    message: str
    context: Dict[str, Any] = {}
    history: List[BaseMessage] = []


class AgentOutput(BaseModel):
    """Output from agent processing."""
    result: Dict[str, Any]
    reasoning: str
    confidence: float
    next_agent: Optional[str] = None


class BaseAgent(ABC):
    """Base class for all DARTSGPT agents."""
    
    def __init__(self, llm: Optional[ChatOpenAI] = None):
        """Initialize the agent.
        
        Args:
            llm: Language model to use. Defaults to GPT-4.
        """
        self.llm = llm or ChatOpenAI(
            model="gpt-4-turbo-preview",
            temperature=0,
            max_tokens=2000
        )
        self._tools = None
        self._executor = None
    
    @property
    def name(self) -> str:
        """Agent name."""
        return self.__class__.__name__
    
    @abstractmethod
    def get_tools(self) -> List[BaseTool]:
        """Get tools available to this agent."""
        pass
    
    @abstractmethod
    def get_prompt_template(self) -> PromptTemplate:
        """Get the prompt template for this agent."""
        pass
    
    @abstractmethod
    async def process(self, input_data: AgentInput) -> AgentOutput:
        """Process the input and return output.
        
        Args:
            input_data: Input data for processing
            
        Returns:
            AgentOutput with results
        """
        pass
    
    def validate_input(self, input_data: AgentInput) -> bool:
        """Validate input data.
        
        Args:
            input_data: Input to validate
            
        Returns:
            True if valid, False otherwise
        """
        return bool(input_data.message.strip())
    
    def format_context(self, context: Dict[str, Any]) -> str:
        """Format context for inclusion in prompts.
        
        Args:
            context: Context dictionary
            
        Returns:
            Formatted context string
        """
        if not context:
            return "No additional context provided."
        
        formatted = []
        for key, value in context.items():
            formatted.append(f"{key}: {value}")
        
        return "\n".join(formatted)