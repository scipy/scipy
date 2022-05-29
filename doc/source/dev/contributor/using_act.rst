.. _using-act:

========================================
`act` for running GitHub Actions locally
========================================

``act`` is a tool offered by Nektos which provides a handy way to run GitHub
Actions locally using Docker. ``act`` provides a quick way to validate your
changes on the CI locally, without committing/pushing your changes to the
workflows to trigger and validate the same. It leads to fast feedback and its
compatibility as a local task runner, to validate all our CI jobs makes it a
handy tool.

``act`` can be set up locally with Homebrew, Chocolatey or even a simple BASH
script. To set it up using the BASH script, just push the following command on
your terminal::

  curl https://raw.githubusercontent.com/nektos/act/master/install.sh | sudo bash

Using Homebrew you can set it up via: ``brew install act``.

The next step is to define the custom image that we can use to run our actions
locally. ``act`` provides a micro, medium and larger Docker image for Ubuntu
GitHub runner. ``act`` does not support Windows and macOS images yet.

While running ``act`` for the first time, we can define the image that we would
like to utilize for our local CI runs. The configuration is saved inside the
``~/.actrc`` file.

In a GitHub repository, while running ``act`` for the first time, it will find
the ``./.github/workflows`` and all the workflows present. To checkout the jobs
listed as part of the GitHub Actions CI, push the following command::

  act -l

It will list all the jobs and you can pick up the particular jobs you wish to
run. If you are looking to run a particular job, push in the following
command::

  act -j <JOB_NAME>

To run the job in dry run, push in the following command::

  act -n

To run the job with verbose logging, push in the following command::

  act -v

To reuse the containers in `act` to maintain state, push in the following command::

  act -j <JOB_NAME> --bind --reuse

It is recommended to comment out GitHub specific events like
``github.repository`` or ``github.event.head_commit.message``. If you are using
environment variables, in your action, it is recommended to have a
``my.secrets`` file and supply these environment variables to the ``act`` by
pushing the following command::

  act --secret-file my.secrets

If the environment variables are supplied via `.env` file, use the following
command::

  act --env-file my.env

