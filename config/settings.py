"""Configuration settings for DARTSGPT."""

import os
from pathlib import Path
from typing import Optional
from dotenv import load_dotenv

# Load environment variables
load_dotenv()


class Settings:
    """Application settings."""
    
    def __init__(self):
        # OpenAI Configuration
        self.openai_api_key = os.getenv("OPENAI_API_KEY", "")
        self.openai_model = os.getenv("OPENAI_MODEL", "gpt-4o-mini")
        
        # Vector Store Configuration
        self.vector_store_type = os.getenv("VECTOR_STORE_TYPE", "chroma")
        self.vector_store_path = Path(os.getenv("VECTOR_STORE_PATH", "./embeddings"))
        
        # LangChain Configuration
        self.langchain_api_key = os.getenv("LANGCHAIN_API_KEY")
        self.langchain_tracing_v2 = os.getenv("LANGCHAIN_TRACING_V2", "true").lower() == "true"
        self.langchain_project = os.getenv("LANGCHAIN_PROJECT", "dartsgpt")
        
        # DARTS Configuration
        self.darts_models_path = Path(os.getenv("DARTS_MODELS_PATH", "./knowledge/darts-code/open-darts-main/models"))
        self.darts_development_path = Path(os.getenv("DARTS_DEVELOPMENT_PATH", "./knowledge/darts-code/darts-models-development"))
        self.darts_templates_path = Path(os.getenv("DARTS_TEMPLATES_PATH", "./knowledge/templates"))
        
        # Application Configuration
        self.log_level = os.getenv("LOG_LEVEL", "INFO")
        self.debug_mode = os.getenv("DEBUG_MODE", "false").lower() == "true"
        self.max_generation_attempts = int(os.getenv("MAX_GENERATION_ATTEMPTS", "3"))
        self.validation_timeout = int(os.getenv("VALIDATION_TIMEOUT", "30"))
        
        # Streamlit Configuration
        self.streamlit_server_port = int(os.getenv("STREAMLIT_SERVER_PORT", "8501"))
        self.streamlit_server_address = os.getenv("STREAMLIT_SERVER_ADDRESS", "localhost")
        
        # Ensure paths are absolute
        self.vector_store_path = self.vector_store_path.resolve()
        self.darts_models_path = self.darts_models_path.resolve()
        self.darts_development_path = self.darts_development_path.resolve()
        self.darts_templates_path = self.darts_templates_path.resolve()
        
        # Create directories if they don't exist
        self.vector_store_path.mkdir(parents=True, exist_ok=True)
        self.darts_templates_path.mkdir(parents=True, exist_ok=True)


# Global settings instance
settings = Settings()