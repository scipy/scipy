ARG BASE_CONTAINER=trallard/scipy-dev:latest
FROM ${BASE_CONTAINER} as clone

COPY --chown=gitpod . /tmp/scipy_repo
RUN git clone --depth 1 file:////tmp/scipy_repo /tmp/scipy

# Using the Scipy-dev Docker image as a base
# This way, we ensure we have all the needed compilers and dependencies
# while reducing the build time
ARG BASE_CONTAINER=trallard/scipy-dev:latest
FROM ${BASE_CONTAINER} as build

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
# Install scipy dev dependencies
COPY --from=clone --chown=gitpod /tmp/scipy ${WORKSPACE}

WORKDIR ${WORKSPACE}

# Build scipy to populate the cache used by ccache
RUN conda activate ${CONDA_ENV} && \
    python setup.py build_ext --inplace && \
    ccache -s

# gitpod will load the repository into /workspace/scipy. We remove the
# directoy from the image to prevent conflicts
RUN sudo rm -rf /workspace/scipy

USER gitpod
