"""Enhanced RAG system for DARTS knowledge."""

import json
from pathlib import Path
from typing import Dict, List, Optional

from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.document_loaders import DirectoryLoader, PythonLoader, TextLoader
from langchain_community.vectorstores import Chroma
from langchain_core.documents import Document
from langchain_openai import OpenAIEmbeddings

from config import settings
from knowledge.model_analyzer import DARTSModelAnalyzer
from knowledge.templates.template_database import TemplateDatabase


class DARTSKnowledgeBase:
    """Enhanced knowledge base for DARTS models and documentation."""

    def __init__(self):
        self.embeddings = OpenAIEmbeddings(openai_api_key=settings.openai_api_key)
        self.vector_store_path = settings.vector_store_path / "darts_knowledge"
        self.vector_store = None
        self.model_analyzer = DARTSModelAnalyzer(settings.darts_models_path)
        self.template_db = TemplateDatabase()
        self.text_splitter = RecursiveCharacterTextSplitter(
            chunk_size=2000,
            chunk_overlap=200,
            separators=["\n\n", "\n", " ", ""],
        )

    def initialize(self, force_rebuild: bool = False):
        """Initialize or load the knowledge base."""
        if self.vector_store_path.exists() and not force_rebuild:
            print("Loading existing vector store...")
            self.vector_store = Chroma(
                persist_directory=str(self.vector_store_path),
                embedding_function=self.embeddings,
            )
        else:
            print("Building new vector store...")
            self._build_vector_store()

    def _build_vector_store(self):
        """Build the vector store from scratch."""
        documents = []

        # 1. Index DARTS models
        print("Indexing DARTS models...")
        model_docs = self._index_models()
        documents.extend(model_docs)

        # 2. Index templates
        print("Indexing templates...")
        template_docs = self._index_templates()
        documents.extend(template_docs)

        # 3. Index documentation
        print("Indexing documentation...")
        doc_docs = self._index_documentation()
        documents.extend(doc_docs)

        # 4. Index test suite
        print("Indexing test suite...")
        test_docs = self._index_test_suite()
        documents.extend(test_docs)

        # Create vector store
        self.vector_store = Chroma.from_documents(
            documents=documents,
            embedding=self.embeddings,
            persist_directory=str(self.vector_store_path),
        )
        self.vector_store.persist()
        print(f"Vector store created with {len(documents)} documents")

    def _index_models(self) -> List[Document]:
        """Index all DARTS models."""
        documents = []
        
        # Analyze all models
        model_metadata = self.model_analyzer.analyze_all_models()
        
        for metadata in model_metadata:
            # Index the model file
            model_path = Path(metadata.path)
            if model_path.exists():
                with open(model_path, "r") as f:
                    content = f.read()
                
                # Create document with metadata
                doc = Document(
                    page_content=content,
                    metadata={
                        "source": metadata.path,
                        "type": "model",
                        "model_name": metadata.name,
                        "physics_type": metadata.physics_type,
                        "grid_type": metadata.grid_type,
                        "features": ",".join(metadata.features),
                        "complexity": metadata.complexity,
                    }
                )
                documents.append(doc)
                
                # Also index the main.py if it exists
                main_path = model_path.parent / "main.py"
                if main_path.exists():
                    with open(main_path, "r") as f:
                        main_content = f.read()
                    
                    main_doc = Document(
                        page_content=main_content,
                        metadata={
                            "source": str(main_path),
                            "type": "main",
                            "model_name": metadata.name,
                            "physics_type": metadata.physics_type,
                        }
                    )
                    documents.append(main_doc)
        
        return documents

    def _index_templates(self) -> List[Document]:
        """Index template information."""
        documents = []
        
        for template_name, template_info in self.template_db.templates.items():
            # Create a comprehensive template document
            content = f"""
Template: {template_info.name}
Physics Type: {template_info.physics_type}
Description: {template_info.description}
Complexity: {template_info.complexity}
Features: {', '.join(template_info.features)}

Use Cases:
{chr(10).join(f'- {uc}' for uc in template_info.example_use_cases)}

Required Parameters:
{chr(10).join(f'- {param}' for param in template_info.required_parameters)}

Optional Parameters:
{chr(10).join(f'- {param}' for param in template_info.optional_parameters)}
"""
            
            doc = Document(
                page_content=content,
                metadata={
                    "source": f"template_db/{template_name}",
                    "type": "template_info",
                    "template_name": template_name,
                    "physics_type": template_info.physics_type,
                    "complexity": template_info.complexity,
                    "features": ",".join(template_info.features),
                }
            )
            documents.append(doc)
        
        return documents

    def _index_documentation(self) -> List[Document]:
        """Index DARTS documentation."""
        documents = []
        
        # Index README files
        readme_patterns = ["**/README.md", "**/README.rst", "**/readme.md"]
        for pattern in readme_patterns:
            for readme_path in settings.darts_models_path.parent.rglob(pattern):
                try:
                    with open(readme_path, "r") as f:
                        content = f.read()
                    
                    doc = Document(
                        page_content=content,
                        metadata={
                            "source": str(readme_path),
                            "type": "documentation",
                            "doc_type": "readme",
                        }
                    )
                    documents.append(doc)
                except Exception as e:
                    print(f"Error reading {readme_path}: {e}")
        
        return documents

    def _index_test_suite(self) -> List[Document]:
        """Index the test suite for patterns."""
        documents = []
        
        # Index Run Test Suite 2.py
        test_suite_path = settings.darts_models_path.parent.parent / "Run Test Suite 2.py"
        if test_suite_path.exists():
            with open(test_suite_path, "r") as f:
                content = f.read()
            
            # Split into chunks
            chunks = self.text_splitter.split_text(content)
            
            for i, chunk in enumerate(chunks):
                doc = Document(
                    page_content=chunk,
                    metadata={
                        "source": str(test_suite_path),
                        "type": "test_suite",
                        "chunk": i,
                    }
                )
                documents.append(doc)
        
        return documents

    def search(self, query: str, k: int = 5, filter_dict: Optional[Dict] = None) -> List[Document]:
        """Search the knowledge base."""
        if not self.vector_store:
            raise ValueError("Knowledge base not initialized")
        
        if filter_dict:
            return self.vector_store.similarity_search(query, k=k, filter=filter_dict)
        else:
            return self.vector_store.similarity_search(query, k=k)

    def search_models_by_physics(self, physics_type: str, k: int = 10) -> List[Document]:
        """Search for models by physics type."""
        return self.search(
            f"physics type {physics_type} model",
            k=k,
            filter_dict={"physics_type": physics_type}
        )

    def search_templates(self, query: str, k: int = 5) -> List[Document]:
        """Search for templates."""
        return self.search(
            query,
            k=k,
            filter_dict={"type": "template_info"}
        )

    def get_model_examples(self, physics_type: str, features: List[str] = []) -> List[Document]:
        """Get example models for a physics type and features."""
        query = f"{physics_type} model with features: {' '.join(features)}"
        return self.search(query, k=5, filter_dict={"type": "model"})

    def get_parameter_examples(self, template_name: str) -> List[Document]:
        """Get parameter examples for a template."""
        return self.search(
            f"parameters for {template_name}",
            k=5,
            filter_dict={"model_name": template_name}
        )