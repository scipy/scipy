# 
# Dockerfile for Scipy development 
# Meson build system - https://scipy.github.io/devdocs/dev/contributor/meson.html#full-details-and-explanation
# https://scipy.github.io/devdocs/dev/contributor/conda_guide.html#conda-guide
# 
# Usage: 
# -------
# 
# To make a local build of the container, from the root directory:
# docker build  --rm -f "./tools/docker_dev/meson.Dockerfile" -t <build-tag> "."  
# 
# To use the container use the following command. It assumes that you are in
# the root folder of the scipy git repository, making it available as
# /home/scipy in the container. Whatever changes you make to that directory
# are visible in the host and container.
# The docker image is retrieved from the scipy dockerhub repository
#
# docker run --rm -it -v $(pwd):/home/scipy scipy:<image-tag>
# 
# By default the container will activate the conda environment scipy-dev
# which contains all the dependencies needed for SciPy development
# 
# To build Scipy run: python dev.py --build-only -j2 
# For the all-in-one (configure,build,test SciPy and docs) use: python dev.py
#
# To run the tests use: python dev.py -n 
# 
# This image is based on: Ubuntu 20.04 (focal)
# https://hub.docker.com/_/ubuntu/?tab=tags&name=focal
# OS/ARCH: linux/amd64
ARG ROOT_CONTAINER=gitpod/workspace-base:latest
ARG BASE_CONTAINER=${ROOT_CONTAINER}

# hadolint ignore=DL3006
FROM ${BASE_CONTAINER}

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
    SHELL=/bin/bash  \
    GP_USER=gitpod \
    GP_GROUP=gitpod \
    GP_UID=33333 

ENV CONDA_VERSION="${conda_version}" \
    MINIFORGE_VERSION="${miniforge_version}" \
    CONDA_ENV=scipy-dev \
    PATH=${CONDA_DIR}/bin:$PATH 

# -----------------------------------------------------------------------------
# ---- OS dependencies ----
ENV DEBIAN_FRONTEND noninteractive

USER root

# Change default shell - this avoids issues with Conda later
SHELL ["/bin/bash", "--login", "-o", "pipefail", "-c"]

# hadolint ignore=DL3008
RUN apt-get update && \ 
    apt-get install -yq --no-install-recommends \
    ca-certificates \
    ccache \
    dirmngr \
    gnupg \
    gpg-agent \
    libatlas-base-dev \
    vim \
    wget && \
    # this needs to be done after installing dirmngr
    apt-key adv --keyserver keyserver.ubuntu.com --recv-key C99B11DEB97541F0 && \ 
    apt-add-repository https://cli.github.com/packages && \ 
    apt-get install -yq --no-install-recommends \
    gh && \ 
    locale-gen en_US.UTF-8 && \
    apt-get clean && \
    rm -rf /var/cache/apt/* &&\
    rm -rf /var/lib/apt/lists/* &&\
    rm -rf /tmp/*

# -----------------------------------------------------------------------------
# ---- Copy needed files ----
# Copy multiple scripts - fix directory permissions and 
# basic workspace configurations
COPY ./tools/docker_dev/workspace_config /usr/local/bin/workspace_config
COPY ./tools/docker_dev/fix_permissions /usr/local/bin/fix_permissions

RUN chmod a+rx /usr/local/bin/workspace_config && \
    chmod a+rx /usr/local/bin/fix_permissions 

# -----------------------------------------------------------------------------
RUN mkdir -p "${CONDA_DIR}" && \
    # transfer conda path ownership and ensure it is user writable
    chown -R ${GP_USER}:${GP_GROUP} ${CONDA_DIR} && \
    fix_permissions ${CONDA_DIR} && \
    workspace_config

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
COPY environment.yml /tmp/environment.yml

RUN conda env create -f /tmp/environment.yml && \
    conda activate ${CONDA_ENV} && \
    # needed for docs rendering later on
    python -m pip install --no-cache-dir sphinx-autobuild && \
    conda install ccache -y && \
    # need to use sudo to remove tmp files
    sudo rm -rf /tmp/* && \
    conda clean --all -f -y && \
    # for good measure after installing things
    fix_permissions ${CONDA_DIR} 

# -----------------------------------------------------------------------------
USER ${GP_USER}

WORKDIR $HOME/scipy
