#!/bin/bash

# Version management script for RNA-seq Toolkit Docker image

VERSION_FILE="VERSION"
CHANGELOG_FILE="CHANGELOG.md"

# Function to get current version
get_current_version() {
    if [ -f "$VERSION_FILE" ]; then
        cat "$VERSION_FILE" | tr -d '\n'
    else
        echo "0.0.0"
    fi
}

# Function to increment version based on type
increment_version() {
    local version=$1
    local increment_type=$2
    
    # Split version into components
    IFS='.' read -r major minor patch <<< "$version"
    
    case $increment_type in
        "major")
            major=$((major + 1))
            minor=0
            patch=0
            ;;
        "minor")
            minor=$((minor + 1))
            patch=0
            ;;
        "patch")
            patch=$((patch + 1))
            ;;
        *)
            echo "Invalid increment type. Use: major, minor, or patch"
            exit 1
            ;;
    esac
    
    echo "$major.$minor.$patch"
}

# Function to update version file
update_version() {
    local new_version=$1
    echo "$new_version" > "$VERSION_FILE"
    echo "Version updated to: $new_version"
}

# Function to update changelog
update_changelog() {
    local new_version=$1
    local date=$(date '+%Y-%m-%d')
    
    # Create backup of changelog
    cp "$CHANGELOG_FILE" "$CHANGELOG_FILE.bak"
    
    # Add new version entry
    {
        echo "# Changelog"
        echo ""
        echo "All notable changes to the RNA-seq Bioinformatics Toolkit Docker image will be documented in this file."
        echo ""
        echo "The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),"
        echo "and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)."
        echo ""
        echo "## [Unreleased]"
        echo ""
        echo "### Added"
        echo "- "
        echo ""
        echo "### Changed"
        echo "- "
        echo ""
        echo "### Fixed"
        echo "- "
        echo ""
        echo "## [$new_version] - $date"
        echo ""
        echo "### Added"
        echo "- Version bump to $new_version"
        echo ""
        # Add rest of changelog, skipping the header
        tail -n +8 "$CHANGELOG_FILE.bak"
    } > "$CHANGELOG_FILE"
    
    rm "$CHANGELOG_FILE.bak"
    echo "Changelog updated with version $new_version"
}

# Function to create git tag
create_git_tag() {
    local version=$1
    if git rev-parse --git-dir > /dev/null 2>&1; then
        git tag -a "v$version" -m "Release version $version"
        echo "Git tag v$version created"
    else
        echo "Not in a git repository, skipping tag creation"
    fi
}

# Main script logic
case "$1" in
    "current")
        echo "Current version: $(get_current_version)"
        ;;
    "major"|"minor"|"patch")
        current_version=$(get_current_version)
        new_version=$(increment_version "$current_version" "$1")
        update_version "$new_version"
        update_changelog "$new_version"
        create_git_tag "$new_version"
        ;;
    "set")
        if [ -z "$2" ]; then
            echo "Usage: $0 set <version>"
            exit 1
        fi
        update_version "$2"
        update_changelog "$2"
        create_git_tag "$2"
        ;;
    *)
        echo "Usage: $0 {current|major|minor|patch|set <version>}"
        echo ""
        echo "Commands:"
        echo "  current     - Show current version"
        echo "  major       - Increment major version (1.0.0 -> 2.0.0)"
        echo "  minor       - Increment minor version (1.0.0 -> 1.1.0)"
        echo "  patch       - Increment patch version (1.0.0 -> 1.0.1)"
        echo "  set <ver>   - Set specific version"
        echo ""
        echo "Current version: $(get_current_version)"
        exit 1
        ;;
esac
