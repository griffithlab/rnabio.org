#!/bin/bash

# Build the Docker image for x86_64 architecture
docker build --platform linux/amd64 -t rnaseq_toolkit:latest .
