## Overview

This is the github repo for the RNA-Seq Workshop [website](http://rnabio.org). The website uses the static site generator [jekyll](https://jekyllrb.com/) and is based on the [Alembic theme](https://github.com/daviddarnes/alembic). Development occurs on the dev branch, the live site is located on the master branch.

## Installation

To install this site locally run the following commands:

1. Clone the repo and cd into it `$ git clone git@github.com:griffithlab/rnabio.org.git`
2. Install the bundler `$ gem install bundler`
3. Install gems `$ bundle install`
4. run jekyll and watch for changes `$ bundle exec jekyll serve --watch`

The site should now be running on localhost port 4000. Changes to files will show up interactively on localhost:4000

**Note:** The _config.yml file is only read during the initial serve, changing this file will require re-running step 4 for changes to appear.

## Docker Installation

To install the site locally with docker run the following commands:

1. Clone the repo `$ git clone git@github.com:griffithlab/rnabio.org.git`
2. Pull the docker image `$ docker pull griffithlab/rnabiodev:0.0.2`
3. Run the docker `$ docker run -p 4000:4000 -v ~/git/rnabio.org/:/opt/git/rnabio.org -it griffithlab/rnabiodev:0.0.2`

Make sure that the above command has correct path to cloned git repo (The first part specified with -v option). The site should now be running on localhost port 4000. Changes to files will show up interactively on localhost:4000.

**Note:** The Dockerfile for this image is maintained in [this repo](https://github.com/griffithlab/rnabio.org/blob/master/docker/site/Dockerfile) and the image itself is hosted on [dockerhub](https://hub.docker.com/r/griffithlab/rnabiodev).

**Note:** The _config.yml file is only read during the initial serve, changing this file will require re-running step 3 for changes to appear.

## Adding Course Content

Course content is located in the _posts directory, course pages must be named following the format: year-month-day-name.md for example `0000-00-00-name.md`. This naming is important for ordering course content, further the front matter tag in the markdown file should include date as well for the same reason. Additional useful tags are `categories` for locating a course in a specific module. An example of front matter tags in a markdown file is supplied below:
```
---
title: Introduction to RNA-Seq
categories:
    - Day 1
feature_image: "https://unsplash.it/1200/400?image=200"
date: 0001-01-01
---
```
### Adjusting The Navigation Header and Pages
The Navigation Header is defined in the `_config.yml` file. Navigation header pages are markdown files and are located in the root directory. For example the about page is located at /about.md. All markdown files should have front matter tags, see jekyll docs.

### _layouts and Page Display

In General all markdown files should have a `layout:` front matter tag. The exception to this are files in the _posts directory which default to a post layout format. Layouts are defined in the _layouts directory which consist of html files specifying content. Jekyll uses liquid templating and this syntax is commonly seen throughout these html pages. A command such as `{% include post-list.html %}` will look in the _includes directory and add in the html content for post-list.html.

### Hints and Tips

### License

The code base for this website originates from the [alembic jekyll theme](https://github.com/daviddarnes/alembic) by David Darnes and is licensed under an [MIT license](https://github.com/griffithlab/pmbio.org/blob/master/LICENSE). Course content, is licensed under [CC-BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) unless otherwise stated.
