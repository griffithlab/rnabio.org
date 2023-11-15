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
#RUN wget https://raw.githubusercontent.com/griffithlab/rnabio.org/master/Gemfile
#RUN wget https://raw.githubusercontent.com/griffithlab/rnabio.org/master/Gemfile.lock
WORKDIR /tmp/rnabio.org
RUN gem install bundler &&\
    bundle install &&\
    rm -rf /tmp/rnabio.org
WORKDIR /usr/local/bin
RUN echo "#!/bin/bash\nbundle exec jekyll serve --host 0.0.0.0 --watch" > start_jekyll.sh &&\
    chmod +x start_jekyll.sh &&\
    mkdir -p /opt/git/rnabio.org
WORKDIR /opt/git/rnabio.org

################################################################################
##################### Set Entrypoint ############################
ENTRYPOINT ["/usr/local/bin/start_jekyll.sh"]

