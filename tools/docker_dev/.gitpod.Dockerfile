# Using the Scipy-dev Docker image as a base
# This way, we ensure we have all the needed compilers and dependencies
# while reducing the build time
ARG BASE_CONTAINER=scipy/scipy-dev:latest
FROM ${BASE_CONTAINER}

# -----------------------------------------------------------------------------
# ---- OS dependencies ----
# hadolint ignore=DL3008
RUN apt-get update && \ 
    apt-get install -yq --no-install-recommends \
    bash-completion \
    dirmngr \
    gpg-agent \
    git \
    htop \
    jq \
    less \
    locales \
    lsof \
    man-db \
    sudo \
    ssl-cert \
    time \
    zip \
    unzip && \
    # this needs to be done after installing dirmngr
    apt-key adv --keyserver keyserver.ubuntu.com --recv-key C99B11DEB97541F0 && \ 
    apt-add-repository https://cli.github.com/packages && \ 
    apt-get install -yq --no-install-recommends \
    gh \ 
    && locale-gen en_US.UTF-8 \
    && apt-get clean \
    && rm -rf /var/cache/apt/* \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/*

# -----------------------------------------------------------------------------
# ---- Gitpod user ----
# by default we are using gitpod as user and group
ENV LANG=en_US.UTF-8 \
    HOME=/home/gitpod \
    GP_USER=gitpod \
    GP_GROUP=gitpod \
    CONDA_DIR=/opt/conda \ 
    WORKSPACE=/workspace/scipy/

# Change default shell - this avoids issues with Conda later
SHELL ["/bin/bash","--login", "-o", "pipefail", "-c"]

# Copy multiple scripts - fix directory permissions and 
# basic workspace configurations
COPY ./tools/docker_dev/fix_permissions /usr/local/bin/fix_permissions
COPY ./tools/docker_dev/workspace_config /usr/local/bin/workspace_config

RUN chmod a+rx /usr/local/bin/fix_permissions && \
    chmod a+rx /usr/local/bin/workspace_config

# '-l': see https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#user
RUN useradd -l -u 33333 -G sudo -md "${HOME}" -s /bin/bash -p ${GP_GROUP} ${GP_USER} && \
    # passwordless sudo for users in the 'sudo' group
    sed -i.bkp -e 's/%sudo\s\+ALL=(ALL\(:ALL\)\?)\s\+ALL/%sudo ALL=NOPASSWD:ALL/g' /etc/sudoers && \
    # ensuring we can use conda develop
    chown -R ${GP_USER}:${GP_GROUP} ${CONDA_DIR} 

# Using root to fix permissions - make sure we change to 'gitpod' always
USER root

RUN mkdir -p ${WORKSPACE} && \ 
    echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile && \ 
    conda init bash && \
    # Making the directories human writable - this is needed to use `conda develop`
    # $HOME is by default user writable
    fix_permissions ${CONDA_DIR}/envs/scipydev/ && \
    fix_permissions $CONDA_DIR && \
    workspace_config 


# Ensure the following happens in the workspace
WORKDIR ${WORKSPACE}
RUN git clone https://github.com/scipy/scipy.git  --depth 1 --single-branch . && \
    conda activate scipydev && \
    python setup.py build_ext --inplace && \
    conda develop . && \ 
    python -m pip install --upgrade pip --no-cache-dir && \ 
    python -m pip install sphinx-autobuild --no-cache-dir

# Always favour the least-privileged user, in this case gitpod
USER gitpod

# Ensure we are in the correct directory
WORKDIR ${WORKSPACE}
