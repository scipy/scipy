.. _git-start:

====================================
Getting started with Git development
====================================

This section and the next describe in detail how to set up git for working
with the SciPy source code.  If you have git already set up, skip to
:ref:`development-workflow`.

Basic Git setup
###############

* :ref:`git-intro`.
* Introduce yourself to Git::

      git config --global user.email you@yourdomain.example.com
      git config --global user.name "Your Name Comes Here"

.. _forking:

Making your own copy (fork) of SciPy
####################################

You need to do this only once.  The instructions here are very similar
to the instructions at http://help.github.com/forking/ - please see that
page for more detail.  We're repeating some of it here just to give the
specifics for the SciPy_ project, and to suggest some default names.

Set up and configure a github_ account
======================================

If you don't have a github_ account, go to the github_ page, and make one.

You then need to configure your account to allow write access - see the
``Generating SSH keys`` help on `github help`_.

Create your own forked copy of SciPy_
=========================================

#. Log into your github_ account.
#. Go to the SciPy_ github home at `SciPy github`_.
#. Click on the *fork* button:

   .. image:: forking_button.png

   After a short pause, you should find yourself at the home page for
   your own forked copy of SciPy_.

.. include:: git_links.inc


.. _set-up-fork:

Set up your fork
################

First you follow the instructions for :ref:`forking`.

Overview
========

::

   git clone https://github.com/your-user-name/scipy.git
   cd scipy
   git remote add upstream https://github.com/scipy/scipy.git

In detail
=========

Clone your fork
---------------

#. Clone your fork to the local computer with ``git clone
   https://github.com/your-user-name/scipy.git``
#. Investigate.  Change directory to your new repo: ``cd scipy``. Then
   ``git branch -a`` to show you all branches.  You'll get something
   like::

      * main
      remotes/origin/main

   This tells you that you are currently on the ``main`` branch, and
   that you also have a ``remote`` connection to ``origin/main``.
   What remote repository is ``remote/origin``? Try ``git remote -v`` to
   see the URLs for the remote.  They will point to your github_ fork.

   Now you want to connect to the upstream `SciPy github`_ repository, so
   you can merge in changes from trunk.

.. _linking-to-upstream:

Linking your repository to the upstream repo
--------------------------------------------

::

   cd scipy
   git remote add upstream https://github.com/scipy/scipy.git

``upstream`` here is just the arbitrary name we're using to refer to the
main SciPy_ repository at `SciPy github`_.

Just for your own satisfaction, show yourself that you now have a new
'remote', with ``git remote -v show``, giving you something like::

   upstream	https://github.com/scipy/scipy.git (fetch)
   upstream	https://github.com/scipy/scipy.git (push)
   origin	https://github.com/your-user-name/scipy.git (fetch)
   origin	https://github.com/your-user-name/scipy.git (push)

To keep in sync with changes in SciPy, you want to set up your repository
so it pulls from ``upstream`` by default.  This can be done with::

   git config branch.main.remote upstream
   git config branch.main.merge refs/heads/main

Your config file should now look something like (from
``$ cat .git/config``)::

   [core]
           repositoryformatversion = 0
           filemode = true
           bare = false
           logallrefupdates = true
           ignorecase = true
           precomposeunicode = false
   [remote "origin"]
           url = https://github.com/your-user-name/scipy.git
           fetch = +refs/heads/*:refs/remotes/origin/*
   [remote "upstream"]
           url = https://github.com/scipy/scipy.git
           fetch = +refs/heads/*:refs/remotes/upstream/*
   [branch "main"]
           remote = upstream
           merge = refs/heads/main

.. include:: git_links.inc
