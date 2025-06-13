"""Configuration module for DARTSGPT."""

from .settings import settings

def get_settings():
    """Get settings instance."""
    return settings

__all__ = ["settings", "get_settings"]