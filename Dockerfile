FROM python:3.8-slim as base

# define shell
SHELL ["/bin/bash", "-c", "-l"]

# switch to root in order to install dependencies
USER root
ENV DEBIAN_FRONTEND noninteractive

# retrieve build args
ARG GITLAB_ACCESS_TOKEN
ARG WRAPPER_DEPS_INSTALL
ARG WRAPPER_DEPS_EXTRA

# install linux-dependencies
RUN apt update
RUN apt install -y git

# add code base
WORKDIR /home/reaxpro/reaxpro_wrappers
COPY osp osp 
COPY tests tests
COPY examples examples
COPY setup.* README.md ./
RUN chmod -R 0777 .

# # install wrappers and their python-dependencies
RUN pip install --upgrade pip
RUN pip install osp-core $WRAPPER_DEPS_INSTALL
RUN pip install .[develop] $WRAPPER_DEPS_EXTRA

################################## target: dev ##################################
from base as dev

# run pytests
RUN pytest tests/test_energy_landscape_refinement.py

# go /app dir
WORKDIR /app
USER root

################################## target: production ##################################

from base as production

# add reaxpro user
RUN useradd -ms /bin/bash reaxpro
USER reaxpro
