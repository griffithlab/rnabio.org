################################################################################
##################### Set Inital Image to work from ############################

FROM ruby:3.1.2

################################################################################
##################### Add Container Labels #####################################

LABEL "Description"="ruby environment for development of rnabio.org site"

################################################################################
##################### Install Necessary gems ############################

WORKDIR /tmp
RUN git clone https://github.com/griffithlab/rnabio.org.git
WORKDIR /tmp/rnabio.org
RUN gem install bundler
RUN bundle install
WORKDIR /usr/local/bin
RUN echo "bundle exec jekyll serve --host 0.0.0.0 --watch" > start_jekyll.sh
RUN chmod +x /usr/local/bin/start_jekyll.sh

################################################################################
##################### Set Entrypoint ############################
ENTRYPOINT ["/usr/local/bin/start_jekyll.sh"]

