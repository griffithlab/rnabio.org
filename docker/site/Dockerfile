################################################################################
##################### Set Inital Image to work from ############################

# work from supported ruby version
FROM ruby:3.1.2

################################################################################
##################### Install necessary gems ############################

WORKDIR /tmp
RUN git clone https://github.com/griffithlab/rnabio.org.git
WORKDIR /tmp/rnabio.org
RUN gem install bundler
RUN bundle install

################################################################################
##################### Add Container Labels #####################################
LABEL "Description"="ruby environment for development of rnabio.org site"


