#!/bin/bash
# Helper script for running DARTSGPT in Docker

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "Docker is not running. Starting Docker Desktop..."
    open -a Docker
    echo "Waiting for Docker to start..."
    while ! docker info > /dev/null 2>&1; do
        sleep 1
    done
    print_success "Docker is now running"
fi

# Build the image if needed
if [[ "$1" == "build" ]] || [[ "$1" == "--build" ]]; then
    print_info "Building Docker image..."
    docker compose build
    print_success "Docker image built successfully"
    shift
fi

# Run the command
if [[ "$1" == "generate" ]] || [[ "$1" == "validate" ]] || [[ "$1" == "execute" ]]; then
    # CLI command
    print_info "Running DARTSGPT CLI command..."
    docker compose run --rm dartsgpt-cli uv run python main.py "$@"
elif [[ "$1" == "ui" ]] || [[ "$1" == "streamlit" ]]; then
    # Streamlit UI
    print_info "Starting DARTSGPT UI..."
    docker compose up dartsgpt
elif [[ "$1" == "shell" ]]; then
    # Interactive shell
    print_info "Starting interactive shell..."
    docker compose run --rm dartsgpt-cli /bin/bash
elif [[ "$1" == "test" ]]; then
    # Run tests
    print_info "Running tests..."
    docker compose run --rm dartsgpt-cli uv run python -m pytest tests/
else
    # Show help
    echo "DARTSGPT Docker Helper"
    echo ""
    echo "Usage: ./docker-run.sh [command] [options]"
    echo ""
    echo "Commands:"
    echo "  build                    Build the Docker image"
    echo "  generate [prompt] [opts] Generate a DARTS model"
    echo "  validate [path]          Validate a DARTS model"
    echo "  execute [path]           Execute a DARTS model (requires Linux with DARTS installed)"
    echo "  ui/streamlit            Start the Streamlit UI"
    echo "  shell                   Start an interactive shell"
    echo "  test                    Run the test suite"
    echo ""
    echo "Options for generate:"
    echo "  --validate             Validate the model after generation"
    echo ""
    echo "Note: Model execution (--execute) requires DARTS to be installed,"
    echo "which only works on Linux x86_64. Use a Linux machine for execution."
    echo ""
    echo "Examples:"
    echo "  ./docker-run.sh build"
    echo "  ./docker-run.sh generate \"Create a water injection model\""
    echo "  ./docker-run.sh ui"
fi