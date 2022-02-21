#
# Dockerfile for quickstart containerised scipy development
#
#
# Usage:
# ------
# If you make changes to the Dockerfile you can rebuild the container:
#
# docker build .
# docker push <build-hash> scipy/scipy-dev
#
# To use the container use the following command. It assumes that you are in
# the root folder of the scipy git repository, making it available as
# /home/scipy in the container. Whatever changes you make to that directory
# are visible in the host and container.
# The docker image is retrieved from the scipy dockerhub repository
#
# docker run -it --rm -v $PWD/:/home/scipy scipy/scipy-dev /bin/bash
#
# The available Python executables are python3.8 (along with
# pip3.8.
# The image does not come bundled with cython/numpy/pytest/pybind11, those need to
# be installed with pip before trying to build/run scipy.
#
# pip3 install cython numpy pytest pybind11 pythran matplotlib sphinx asv
#
# You can find the detailed guide in:
# Development environment quickstart guide (Docker)
# http://scipy.github.io/devdocs/dev/contributor/quickstart_docker.html#quickstart-docker


# ubuntu focal has python 3.8 as default
FROM ubuntu:focal

RUN apt-get update && apt-get install -y \
	python3-pip \
	build-essential \
	vim \
	libatlas-base-dev \
	gfortran \
	libgfortran4 \
	liblapack-dev \
	curl \
	libgmp-dev \
	libmpfr-dev \
	libsuitesparse-dev \
	libmpc-dev \
	git


# setup pips and pip3.8
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3.8 get-pip.py && rm get-pip.py
