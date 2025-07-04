#!/bin/bash

# Build the Docker image for x86_64 architecture with version tagging

# Read version from VERSION file
if [ -f "VERSION" ]; then
    VERSION=$(cat VERSION | tr -d '\n')
else
    echo "VERSION file not found. Using 'latest' as default."
    VERSION="latest"
fi

# Get build metadata
BUILD_DATE=$(date -u +'%Y-%m-%dT%H:%M:%SZ')
VCS_REF=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "Building RNA-seq Toolkit Docker image version: $VERSION"
echo "Build date: $BUILD_DATE"
echo "VCS ref: $VCS_REF"

# Build with both version tag and latest tag
docker build --platform linux/amd64 \
    --build-arg VERSION="$VERSION" \
    --build-arg BUILD_DATE="$BUILD_DATE" \
    --build-arg VCS_REF="$VCS_REF" \
    -t rnaseq_toolkit:$VERSION \
    -t rnaseq_toolkit:latest \
    .

echo "Build completed successfully!"
echo "Tagged as: rnaseq_toolkit:$VERSION and rnaseq_toolkit:latest"
