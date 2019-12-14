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
# The available Python executables are python3.6 and python3.7 (along with
# pip3.6 and pip3.7.
# The image does not come bundled with cython/numpy/pytest/pybind11, those need to
# be installed with pip before trying to build/run scipy.
#
# pip3 install cython numpy pytest pybind11
#


FROM ubuntu:bionic

RUN apt-get update

RUN apt-get install -y python3.7 python3.7-dev python3-pip git build-essential vim libatlas-base-dev gfortran liblapack-dev curl libgmp-dev libmpfr-dev libsuitesparse-dev libmpc-dev

# setup pips, pip3.6 and pip3.7
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3.7 get-pip.py && python3.6 get-pip.py && rm get-pip.py
