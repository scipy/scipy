# Using the Scipy-dev Docker image as a base
# This way, we ensure we have all the needed compilers and dependencies
# while reducing the build time
ARG BASE_CONTAINER=scipy/scipy-dev:latest
FROM ${BASE_CONTAINER}

# -----------------------------------------------------------------------------
USER root

# -----------------------------------------------------------------------------
# ---- ENV variables ----
# ---- Directories needed ----
ENV WORKSPACE=/workspace/scipy/ \
    CONDA_ENV=scipy-dev

# -----------------------------------------------------------------------------
# Change default shell - this avoids issues with Conda later - note we do need
# login bash here as we are building SciPy inside
# Fix DL4006
SHELL ["/bin/bash","--login", "-o", "pipefail", "-c"]

# -----------------------------------------------------------------------------
# ---- Build Scipy here ----
# Ensure the following happens in the workspace
USER gitpod

WORKDIR ${WORKSPACE}

RUN git clone https://github.com/scipy/scipy.git  --depth 1 --single-branch . && \
    conda activate ${CONDA_ENV} && \
    python setup.py build_ext --inplace && \
    ccache -s
