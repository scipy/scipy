# 
# Dockerfile for Scipy development 
# 
# Usage: 
# -------
# 
# To make a local build of the container, from the 'Docker-dev' directory:
# docker build  --rm -f "Dockerfile" -t <build-tag> "."    
# 
# To use the container use the following command. It assumes that you are in
# the root folder of the scipy git repository, making it available as
# /home/scipy in the container. Whatever changes you make to that directory
# are visible in the host and container.
# The docker image is retrieved from the scipy dockerhub repository
#
# docker run --rm -it -v $(pwd):/home/scipy scipy:<image-tag>
# 
# Once the container is running, you'll need to activate the conda environment
# conda activate scipydev
# 
# And run the tests with
# python runtests.py
# 
# Ubuntu 20.04 (focal)
# https://hub.docker.com/_/ubuntu/?tab=tags&name=focal
# OS/ARCH: linux/amd64
ARG ROOT_CONTAINER=ubuntu:focal-20210325
ARG BASE_CONTAINER=${ROOT_CONTAINER}

# hadolint ignore=DL3006 as base container is tagged
FROM ${BASE_CONTAINER}

# Change default shell - this avoids issues with Conda later
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# -----------------------------------------------------------------------------
# ---- Miniforge installer ----
# Default values can be overridden at build time
# (ARGS are in lower case to distinguish them from ENV)
# Check https://github.com/conda-forge/miniforge/releases
# Conda version
ARG conda_version="4.9.2"
# Miniforge installer patch version
ARG miniforge_patch_number="5"
# Miniforge installer architecture
ARG miniforge_arch="x86_64"
# Python implementation to use 
# can be either Miniforge3 to use Python or Miniforge-pypy3 to use PyPy
ARG miniforge_python="Miniforge3"

# Miniforge archive to install
ARG miniforge_version="${conda_version}-${miniforge_patch_number}"
# Miniforge installer
ARG miniforge_installer="${miniforge_python}-${miniforge_version}-Linux-${miniforge_arch}.sh"
# Miniforge checksum
ARG miniforge_checksum="49dddb3998550e40adc904dae55b0a2aeeb0bd9fc4306869cc4a600ec4b8b47c"

# -----------------------------------------------------------------------------
# ---- Python version to install ----
# Currently Python 3.8
ARG PYTHON_VERSION=default

# ---- Configure environment ----
ENV CONDA_DIR=/opt/conda \
    SHELL=/bin/bash 

ENV CONDA_VERSION="${conda_version}" \
    MINIFORGE_VERSION="${miniforge_version}" \
    CONDA_ENV=scipydev \
    # ---- Gitpod user ----
    GP_USER=gitpod \
    GP_GROUP=gitpod \
    UID=33333

ENV HOME="/home/${GP_USER}" \
    PROJECT_DIR="${GP_HOME}/scipy"\
    PATH=${CONDA_DIR}/bin:$PATH 

# -----------------------------------------------------------------------------
# ---- OS dependencies ----
ENV DEBIAN_FRONTEND noninteractive

# hadolint ignore=DL3008
RUN apt-get update && \ 
    apt-get install -yq --no-install-recommends \
    build-essential \
    ca-certificates \
    git \
    libatlas-base-dev \
    software-properties-common \
    sudo \
    vim \
    wget && \
    apt-get clean && \
    rm -rf /tmp/* && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/cache/apt/* 

# -----------------------------------------------------------------------------
# ---- Copy needed files ----
# Copy multiple scripts - fix directory permissions and 
# basic workspace configurations
COPY ./workspace_config /usr/local/bin/workspace_config
COPY ./fix_permissions /usr/local/bin/fix_permissions

RUN chmod a+rx /usr/local/bin/workspace_config && \
    chmod a+rx /usr/local/bin/fix_permissions 

# -----------------------------------------------------------------------------
# ---- Create gitpod user and group: needed for permissions later ----
# '-l': see https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#user
RUN echo "auth requisite pam_deny.so" >> /etc/pam.d/su && \
    sed -i.bak -e 's/^%admin/#%admin/' /etc/sudoers && \
    sed -i.bak -e 's/^%sudo/#%sudo/' /etc/sudoers && \
    # create /home/gitpod too for later use
    useradd -l -u ${UID} -G sudo -md ${HOME} -s /bin/bash -p ${GP_GROUP} ${GP_USER} && \
    # passwordless sudo for users in the 'sudo' group
    sed -i.bkp -e 's/%sudo\s\+ALL=(ALL\(:ALL\)\?)\s\+ALL/%sudo ALL=NOPASSWD:ALL/g' /etc/sudoers && \
    mkdir -p ${CONDA_DIR} && \
    # transfer conda path ownership and ensure it is user writable
    chown -R ${GP_USER}:${GP_GROUP} ${CONDA_DIR} && \
    fix_permissions ${HOME} && \
    fix_permissions ${CONDA_DIR} && \
    workspace_config

# Setup work directory for backward-compatibility
RUN mkdir ${PROJECT_DIR} && \
    fix_permissions ${HOME}

USER ${GP_USER}

WORKDIR /tmp

# -----------------------------------------------------------------------------
# ---- Installing conda  ----
RUN wget --quiet "https://github.com/conda-forge/miniforge/releases/download/${miniforge_version}/${miniforge_installer}" && \
    echo "${miniforge_checksum} *${miniforge_installer}" | sha256sum --check && \
    /bin/bash "${miniforge_installer}" -f -b -p $CONDA_DIR && \
    rm "${miniforge_installer}" && \
    # Conda configuration see https://conda.io/projects/conda/en/latest/configuration.html
    echo "conda ${CONDA_VERSION}" >> $CONDA_DIR/conda-meta/pinned && \
    conda config --system --set auto_update_conda false && \
    conda config --system --set show_channel_urls true && \
    # This allows to change the Python version installed by passing the arg `PYTHON_VERSION` at build time
    # then version is added to conda-meta/pinned 
    if [ ! $PYTHON_VERSION = 'default' ]; then conda install --yes python=$PYTHON_VERSION; fi && \
    conda list python | grep '^python ' | tr -s ' ' | cut -d '.' -f 1,2 | sed 's/$/.*/' >> $CONDA_DIR/conda-meta/pinned && \
    conda install --quiet --yes \
    "conda=${CONDA_VERSION}" \
    'pip' && \
    conda update --all --quiet --yes && \
    conda clean --all -f -y && \
    fix_permissions ${CONDA_DIR}


# -----------------------------------------------------------------------------
# ---- Create conda environment ----
# Install SciPy dependencies - since using miniforge no need to add 
# conda-forge channel
# Clean up temporary files 
COPY ./environment.yml /tmp/environment.yml

RUN conda env create -f /tmp/environment.yml && \
    fix_permissions ${CONDA_DIR} &&\
    conda clean --all -f -y 

# -----------------------------------------------------------------------------
# Create a Scipy directory inside /home/ 
# this will be used to mount the local volume
# RUN mkdir ${PROJECT_DIR}
WORKDIR ${PROJECT_DIR}
USER ${GP_USER}