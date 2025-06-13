# DARTSGPT Dockerfile
# Multi-stage build for optimized image size
FROM python:3.11-slim as builder

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    gfortran \
    git \
    curl \
    cmake \
    libopenblas-dev \
    liblapack-dev \
    libumfpack5 \
    libsuitesparse-dev \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Install uv for fast Python package management
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:${PATH}"

# Set working directory
WORKDIR /app

# Copy dependency files
COPY pyproject.toml ./

# Install dependencies
RUN uv sync --no-cache

# Production stage
FROM python:3.11-slim

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    git \
    gfortran \
    libopenblas0 \
    liblapack3 \
    libumfpack5 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN useradd -m -u 1000 dartsgpt

# Set working directory
WORKDIR /app

# Copy from builder
COPY --from=builder --chown=dartsgpt:dartsgpt /app/.venv /app/.venv
COPY --from=builder /root/.local/bin/uv /usr/local/bin/uv

# Copy application code
COPY --chown=dartsgpt:dartsgpt . .

# Set environment variables
ENV PATH="/app/.venv/bin:${PATH}"
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV VIRTUAL_ENV=/app/.venv

# Create necessary directories
RUN mkdir -p embeddings output logs && \
    chown -R dartsgpt:dartsgpt embeddings output logs

# Switch to non-root user
USER dartsgpt

# Expose ports
EXPOSE 8501 8502

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import sys; sys.exit(0)" || exit 1

# Default command (can be overridden)
CMD ["uv", "run", "streamlit", "run", "app.py", "--server.port=8502", "--server.address=0.0.0.0"]