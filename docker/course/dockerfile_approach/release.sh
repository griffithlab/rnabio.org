#!/bin/bash

# Release script for RNA-seq Toolkit Docker image

# Default Docker registry (change as needed)
REGISTRY="your-registry.com"  # Change to your Docker registry
IMAGE_NAME="rnaseq-toolkit"

# Function to show usage
show_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -r, --registry REGISTRY    Docker registry to push to (default: $REGISTRY)"
    echo "  -n, --name NAME            Image name (default: $IMAGE_NAME)"
    echo "  -t, --tag-only             Only create tags, don't push to registry"
    echo "  -h, --help                 Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                                    # Build and push to default registry"
    echo "  $0 --registry docker.io/username     # Push to Docker Hub"
    echo "  $0 --tag-only                        # Only create local tags"
    echo "  $0 --registry ghcr.io/username       # Push to GitHub Container Registry"
}

# Parse command line arguments
TAG_ONLY=false
while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--registry)
            REGISTRY="$2"
            shift 2
            ;;
        -n|--name)
            IMAGE_NAME="$2"
            shift 2
            ;;
        -t|--tag-only)
            TAG_ONLY=true
            shift
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Get version from VERSION file
if [ -f "VERSION" ]; then
    VERSION=$(cat VERSION | tr -d '\n')
else
    echo "ERROR: VERSION file not found. Run './version.sh current' to check or set a version."
    exit 1
fi

echo "=== RNA-seq Toolkit Release Script ==="
echo "Version: $VERSION"
echo "Registry: $REGISTRY"
echo "Image name: $IMAGE_NAME"
echo "Tag only: $TAG_ONLY"
echo ""

# Check if we're in a git repository and get commit info
if git rev-parse --git-dir > /dev/null 2>&1; then
    VCS_REF=$(git rev-parse --short HEAD)
    echo "VCS ref: $VCS_REF"
    
    # Check if there are uncommitted changes
    if ! git diff-index --quiet HEAD --; then
        echo "WARNING: There are uncommitted changes in the repository."
        read -p "Continue anyway? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Aborted."
            exit 1
        fi
    fi
else
    echo "WARNING: Not in a git repository"
fi

echo ""

# Build the image
echo "Building Docker image..."
./build.sh

if [ $? -ne 0 ]; then
    echo "ERROR: Build failed"
    exit 1
fi

# Tag for registry
if [ "$TAG_ONLY" = false ]; then
    echo "Tagging for registry: $REGISTRY/$IMAGE_NAME"
    
    # Tag with version
    docker tag rnaseq_toolkit:$VERSION $REGISTRY/$IMAGE_NAME:$VERSION
    docker tag rnaseq_toolkit:$VERSION $REGISTRY/$IMAGE_NAME:latest
    
    echo "Tagged as:"
    echo "  $REGISTRY/$IMAGE_NAME:$VERSION"
    echo "  $REGISTRY/$IMAGE_NAME:latest"
    echo ""
    
    # Push to registry
    echo "Pushing to registry..."
    docker push $REGISTRY/$IMAGE_NAME:$VERSION
    docker push $REGISTRY/$IMAGE_NAME:latest
    
    if [ $? -eq 0 ]; then
        echo ""
        echo "=== Release Complete ==="
        echo "Successfully pushed version $VERSION to $REGISTRY/$IMAGE_NAME"
        echo ""
        echo "You can now pull the image with:"
        echo "  docker pull $REGISTRY/$IMAGE_NAME:$VERSION"
        echo "  docker pull $REGISTRY/$IMAGE_NAME:latest"
        
        # Create GitHub release notes template if in git repo
        if git rev-parse --git-dir > /dev/null 2>&1; then
            echo ""
            echo "=== Release Notes Template ==="
            echo "Copy this for GitHub release:"
            echo ""
            echo "## RNA-seq Toolkit v$VERSION"
            echo ""
            echo "### Docker Images"
            echo "- \`$REGISTRY/$IMAGE_NAME:$VERSION\`"
            echo "- \`$REGISTRY/$IMAGE_NAME:latest\`"
            echo ""
            echo "### Usage"
            echo "\`\`\`bash"
            echo "docker pull $REGISTRY/$IMAGE_NAME:$VERSION"
            echo "docker run -it $REGISTRY/$IMAGE_NAME:$VERSION"
            echo "\`\`\`"
            echo ""
            echo "### Changes"
            echo "See [CHANGELOG.md](CHANGELOG.md) for detailed changes."
        fi
    else
        echo "ERROR: Failed to push to registry"
        exit 1
    fi
else
    echo "Tag-only mode: Created local tags but did not push to registry"
    echo "Local tags created:"
    echo "  rnaseq_toolkit:$VERSION"
    echo "  rnaseq_toolkit:latest"
fi
