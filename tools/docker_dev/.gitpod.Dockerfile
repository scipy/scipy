# Using the Scipy-dev Docker image as a base
# This way, we ensure we have all the needed compilers and dependencies
# while reducing the build time
ARG BASE_CONTAINER=scipy/scipy-dev:latest
FROM ${BASE_CONTAINER}
# -----------------------------------------------------------------------------

USER root

# -----------------------------------------------------------------------------
# ---- OS dependencies ----
# hadolint ignore=DL3008
RUN apt-get update && \ 
    apt-get install -yq --no-install-recommends \
    bash-completion \
    dirmngr \
    gpg-agent \
    git \
    jq \
    less \
    locales \
    lsof \
    man-db \
    sudo \
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
    UID=3333
# ---- Directories needed ----
ENV HOME=/home/gitpod \
    CONDA_DIR=/opt/conda \ 
    WORKSPACE=/workspace/scipy/ 

# -----------------------------------------------------------------------------
# Change default shell - this avoids issues with Conda later
# Fix DL4006
SHELL ["/bin/bash","--login", "-o", "pipefail", "-c"]

# -----------------------------------------------------------------------------
# ---- Copy needed files ----
# Copy multiple scripts - fix directory permissions and 
# basic workspace configurations
COPY ./tools/docker_dev/fix_permissions /usr/local/bin/fix_permissions
COPY ./tools/docker_dev/workspace_config /usr/local/bin/workspace_config

RUN chmod a+rx /usr/local/bin/fix_permissions && \
    chmod a+rx /usr/local/bin/workspace_config

# -----------------------------------------------------------------------------
# ---- Create gitpod user and group: needed for permissions later ----
# '-l': see https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#user
RUN echo "auth requisite pam_deny.so" >> /etc/pam.d/su && \
    sed -i.bak -e 's/^%admin/#%admin/' /etc/sudoers && \
    sed -i.bak -e 's/^%sudo/#%sudo/' /etc/sudoers && \
    useradd -l -u ${UID} -G sudo -md "${HOME}" -s /bin/bash -p ${GP_GROUP} ${GP_USER} && \
    # passwordless sudo for users in the 'sudo' group
    sed -i.bkp -e 's/%sudo\s\+ALL=(ALL\(:ALL\)\?)\s\+ALL/%sudo ALL=NOPASSWD:ALL/g' /etc/sudoers && \
    # transfer ownership to gitpod user 
    chown -R ${GP_USER}:${GP_GROUP} ${CONDA_DIR} && \ 
    fix_permissions ${CONDA_DIR} && \
    mkdir -p ${WORKSPACE} && \
    chown -R ${GP_USER}:${GP_GROUP} ${WORKSPACE} && \ 
    fix_permissions ${WORKSPACE}  

# workspace configs here
RUN workspace_config 

# Always favour the least-privileged user, in this case gitpod
USER ${GP_USER}
# -----------------------------------------------------------------------------
# ---- Build Scipy here ----
# Ensure the following happens in the workspace
WORKDIR ${WORKSPACE}

RUN git clone https://github.com/scipy/scipy.git  --depth 1 --single-branch . && \
    conda activate scipydev && \
    python setup.py build_ext --inplace && \
    conda develop . && \ 
    python -m pip install --no-cache-dir sphinx-autobuild && \
    rm -rf /tmp/* && \
    # for good measure as we installed things
    fix_permissions ${WORKSPACE} && \
    fix_permissions ${CONDA_DIR}

# Ensure we are in the correct directory
WORKDIR ${WORKSPACE}
