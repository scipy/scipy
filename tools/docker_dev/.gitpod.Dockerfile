# Using the Scipy-dev Docker image as a base
# This way, we ensure we have all the needed compilers and dependencies
# while reducing the build time
ARG BASE_CONTAINER=scipy/scipy-dev:latest
FROM ${BASE_CONTAINER} as build

# -----------------------------------------------------------------------------
USER root

# -----------------------------------------------------------------------------
# ---- ENV variables ----
# ---- Directories needed ----
ENV HOME=/home/scipy \
    CONDA_DIR=/opt/conda \ 
    WORKSPACE=/workspace/scipy/ 

# -----------------------------------------------------------------------------
# Change default shell - this avoids issues with Conda later - note we do need
# login bash here as we are building SciPy inside
# Fix DL4006
SHELL ["/bin/bash","--login", "-o", "pipefail", "-c"]

RUN mkdir -p ${WORKSPACE} 

# -----------------------------------------------------------------------------
# ---- Build Scipy here ----
# Ensure the following happens in the workspace
WORKDIR ${WORKSPACE}

RUN git clone https://github.com/scipy/scipy.git  --depth 1 --single-branch . && \
    conda activate scipydev && \
    python setup.py build_ext --inplace 


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ******************************************************************************
# Ubuntu 20.04 (focal)
# https://hub.docker.com/_/ubuntu/?tab=tags&name=focal
# OS/ARCH: linux/amd64

FROM scipy/scipy-dev:latest as runtime

ARG DEBIAN_FRONTEND=noninteractive

USER root

# -----------------------------------------------------------------------------
# ---- OS dependencies ----
# hadolint ignore=DL3008
RUN apt-get update && \ 
    apt-get install -yq --no-install-recommends \
    bash-completion \
    dirmngr \
    gpg-agent \
    gnupg \
    jq \ 
    less \
    locales \
    lsof \
    man-db \
    ssl-cert \
    time && \
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
# ---- ENV variables ----
# by default we are using gitpod as user and group
ENV LANG=en_US.UTF-8 \
    # ---- Gitpod user ----
    GP_USER=gitpod \
    GP_GROUP=gitpod \
    UID=33333
# ---- Directories needed ----
ENV HOME=/home/gitpod \
    CONDA_DIR=/opt/conda \ 
    WORKSPACE=/workspace/scipy \
    CONDA_ENV=scipydev

ENV PATH=${CONDA_DIR}/bin:$PATH 

# -----------------------------------------------------------------------------
# ---- Copy needed files ----
# Copy multiple scripts - fix directory permissions and 
# basic workspace configurations
COPY ./tools/docker_dev/fix_permissions /usr/local/bin/fix_permissions
COPY ./tools/docker_dev/workspace_config /usr/local/bin/workspace_config

WORKDIR ${WORKSPACE}

RUN chmod a+rx /usr/local/bin/fix_permissions && \
    chmod a+rx /usr/local/bin/workspace_config && \
    workspace_config && \
    chown -R ${GP_USER}:${GP_GROUP} ${WORKSPACE} && \ 
    fix_permissions ${WORKSPACE}  

USER ${GP_USER}

# install Sphinx autobuild - needed for rst preview
RUN python -m pip install --no-cache-dir sphinx-autobuild && \
    rm -rf /tmp/* && \
    ${CONDA_DIR}/bin/conda clean -afy && \
    fix_permissions ${CONDA_DIR}

# Copy build directory from build stage - doing this to avoid issues
# with how gitpod clones repos later
COPY --chown=${GP_USER}:${GP_GROUP} --from=build ${WORKSPACE}/build/ ${HOME}/build/
