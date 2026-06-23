# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Jekyll-based static site ([rnabio.org](https://www.rnabio.org/)) serving as an RNA-seq bioinformatics workshop course. It uses the Alembic Jekyll theme and is hosted on GitHub Pages. Course content is licensed CC-BY-SA 4.0; theme code is MIT.

The live site deploys from the `master` branch. Active development happens on the `dev` branch.

## Audience and Pedagogical Goals

The course targets **advanced trainees** (PhD candidates and postdocs) with strong experimental biology backgrounds but limited computational experience. The aim is not just to run working examples but to build genuine computational intuition and scientific rigor — producing bioinformaticians capable of applying these skills independently to large, multi-dimensional datasets.

Content is organized around **transcriptomics as the exemplar domain**, but the underlying goal is transferable general bioinformatics skills. Hands-on exercises use R and Unix tools in pre-configured cloud environments so students can experiment without fighting setup issues. These exercises are paired with lectures and class discussions that carry the conceptual weight: explaining the genomics and bioinformatics principles that transform a working pipeline into genuine understanding.

## Local Development

**Standard (Ruby/Bundler):**
```bash
bundle install
bundle exec jekyll serve --watch
```
Site runs at http://localhost:4000. Note: changes to `_config.yml` require a server restart.

**Docker (preferred for reproducibility):**
```bash
docker run -p 4000:4000 -v ~/git/rnabio.org/:/opt/git/rnabio.org -it griffithlab/rnabiodev:0.0.4
# Inside container:
cd /opt/git/rnabio.org && bundle exec jekyll serve --watch --host 0.0.0.0
```

## Content Architecture

All course content lives in `_posts/` as Markdown files with a strict naming convention:

```
YYYY-MM-DD-ModuleName-LessonName.md
```

The date prefix controls ordering and module grouping. Module numbers are encoded in the date:
- `0001-01-XX` = Module 1 (Inputs)
- `0002-01-XX` = Module 2 (Alignment)
- `0003-03-XX` = Module 3 (Expression/DESeq2)
- etc.

When renaming/reordering posts, update the date prefix accordingly. The git status shows pending renames like `0003-03-03 → 0003-03-04` — this pattern is normal when inserting new lessons.

Each post uses Jekyll front matter:
```yaml
---
feature_image: "path/to/image"
title: "Lesson Title"
categories:
  - Module Name
tags:
  - tag1
---
```

## Site Structure

- `_posts/` — All course lessons (Markdown)
- `_layouts/` — Page templates (`post.html` for lessons, `page.html` for static pages)
- `_includes/` — Reusable partials (navigation, post headers/footers, etc.)
- `_sass/` — SCSS styles (Alembic theme base + customizations)
- `assets/` — Static files organized by module (`module_0/` through `module_8/`), plus `lectures/`, `scripts/`, `ResourceFiles/`
- `docker/` — Dockerfiles for the dev site environment (versioned)
- `_config.yml` — Jekyll configuration, navigation structure, pagination settings

## Key Configuration

`_config.yml` controls site-wide settings. Navigation items are defined here under the `navigation_pages` key. Pagination is set to 2 posts per page. Jekyll plugins in use: jekyll-sitemap, jekyll-mentions, jekyll-paginate, jekyll-seo-tag, jekyll-redirect-from, jekyll-feed, jemoji.

## Docker Image Maintenance

The dev Docker image is maintained at `docker/site/` with versioned Dockerfiles (`Dockerfile.0.0.1` through `Dockerfile.0.0.4`). The current version is `0.0.4` (Ruby 3.x base). When updating dependencies, increment the version and build/push to `griffithlab/rnabiodev`.
