ARG BASE_CONTAINER=trallard/scipy:20210127-feature-Docker-b7f5aac1ce8009af05e9011152180cc3112cf13f
FROM ${BASE_CONTAINER}

# -----------------------------------------------------------------------------
# ---- OS dependencies ----
# hadolint ignore=DL3008
RUN apt-get update && \ 
    apt-get install -yq --no-install-recommends \
    zip \
    unzip \
    bash-completion \
    build-essential \
    htop \
    jq \
    less \
    locales \
    man-db \
    sudo \
    time \
    multitail \
    lsof \
    ssl-cert \
    && locale-gen en_US.UTF-8 \
    && apt-get clean \
    && rm -rf /var/cache/apt/* \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/*

RUN add-apt-repository -y ppa:git-core/ppa \
    && apt-get install -yq git \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# ---- Gitpod users ----
ENV LANG=en_US.UTF-8 \
    HOME=/home/scipy

# '-l': see https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#user
RUN useradd -l -u 33333 -G sudo -md "${HOME}" -s /bin/bash -p gitpod gitpod \
    # passwordless sudo for users in the 'sudo' group
    && sed -i.bkp -e 's/%sudo\s\+ALL=(ALL\(:ALL\)\?)\s\+ALL/%sudo ALL=NOPASSWD:ALL/g' /etc/sudoers
WORKDIR $HOME
# custom Bash prompt
RUN { echo && echo "PS1='\[\e]0;\u \w\a\]\[\033[01;32m\]\u\[\033[00m\] \[\033[01;34m\]\w\[\033[00m\] \\\$ '" ; } >> .bashrc

### Gitpod user (2) ###
USER gitpod
# use sudo so that user does not get sudo usage info on (the first) login
RUN sudo echo "Running 'sudo' for Gitpod: success" && \
    # create .bashrc.d folder and source it in the bashrc
    mkdir /home/scipy/.bashrc.d && \
    (echo; echo "for i in \$(ls \$HOME/.bashrc.d/*); do source \$i; done"; echo) >> /home/scipy/.bashrc


USER gitpod
