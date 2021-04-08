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
    BUILD_DIR=/opt/conda/envs/scipy-dev/lib/python3.8/site-packages

# -----------------------------------------------------------------------------
# Change default shell - this avoids issues with Conda later - note we do need
# login bash here as we are building SciPy inside
# Fix DL4006
SHELL ["/bin/bash","--login", "-o", "pipefail", "-c"]

# -----------------------------------------------------------------------------
# ---- Build Scipy here ----
# Ensure the following happens in the workspace
USER gitpod

RUN git clone https://github.com/scipy/scipy.git  --depth 1 --single-branch . && \
    conda activate scipy-dev && \
    python setup.py build --build-base=${BUILD_DIR}/scipy && \
    python setup.py install --install-lib ${BUILD_DIR}/scipy


