:orphan:

.. highlight:: console

.. _quickstart-gitpod:


=======================================================
Development environment quickstart guide (gitpod)
=======================================================

This quickstart guide covers:

* using GitPod for your SciPy development environment
* creating a personal fork of the SciPy repository on GitHub
* a quick tour of Gitpod and VSCode
* working on the SciPy documentation in SciPy

Gitpod
-------

`Gitpod`_  is an open source platform for automated and ready-to-code development environments. It enables developers to describe their dev environment as code and start instant and fresh development environments for each new task directly from your browser. This reduces the need to install local development environments and deal with incompatible dependencies.

To be able to use Gitpod you will need to login with your GitHub account, so if you do not have an account yet you will need to create one first. 

Head over to the `Gitpod`_ website and click on the **Continue with GitHub**, you will need to authenticate to your GitHub account in the next window.

Forking the SciPy repository
-----------------------------

The best way to work on Gitpod as a contributor is to first make a fork of the SciPy repository and work on it.

#. Browse to the `SciPy repository on GitHub`_ and `create your own fork`_.

#. Browse to your fork. Your fork will have a URL like
    https://github.com/andyfaff/scipy, except with your GitHub username
    in place of "andyfaff".

Starting Gitpod
----------------
After you have authenticated to Gitpod through GitHub you should see a **Gitpod** button next to the "Code" button.

#. Click the "Gitpod" button. This will open a new tab on your web browser and start building your development environment. Please note this can take a few minutes.

#. Once the build is complete, you will be directed to your workspace, which includes VSCode and all the dependencies you need to work on SciPy. The first time you start your  workspace you will notice that there might be some actions running, this will ensure that you have a development version of SciPy installed and that the docs are prebuilt for you.

#. Once the build is complete you can test the build by entering::

        python runtests.py -v

    ``runtests.py`` is another script in the SciPy root directory. It runs a
    suite of tests that make sure SciPy is working as it should, and ``-v``
    activates the ``–verbose`` option to show all the test output.

Quick workspace tour
---------------------

Your workspace will look similar to the image below:

..image:: _static/gitpod/gitpod-workspace.png 
    :alt: Gitpod workspace screenshot
    :figclass: align-center

Some important parts of the editor have been marked in the image:
1. Your current Python interpreter - by default the ``scipy-dev`` is marked as ``base`` on the status bar but it should be displayed as ``scipy-dev`` on your terminal. You do not need activate the conda environment as this will always be activated for you.
2. Your current branch is always displayed in the status bar.
3. GitHub Pull Requests extension - you can use this to work with Pull Requests from your workspace.
4. Marketplace extensions - we have added some basic extensions to the SciPy Gitpod but you can also install other extensions or syntax highlighting themes for your user and these will be preserved
5. Your workspace directory - by default is is ``/workspace/scipy`` **do not change this** as this is the only directory that is preserved in Gitpod.

Some tools and extensions installed in your workspace are:

*  `GitHub CLI <https://cli.github.com/>`_
*  vim
*  `VSCode rst extension <https://marketplace.visualstudio.com/items?itemName=lextudio.restructuredtext>`_
*  `VSCode Liveserver extension <https://marketplace.visualstudio.com/items?itemName=ritwickdey.LiveServer>`_
*  `VSCode Gitlens extension <https://marketplace.visualstudio.com/items?itemName=eamodio.gitlens>`_
*  `VSCode autodocstrings extension <https://marketplace.visualstudio.com/items?itemName=njpwerner.autodocstring>`_ set to "numpy" by default

Development workflow
-----------------------
:ref:`development-workflow` goes in detail into the SciPy development workflow.

.. note:: You do not need to configure your git username and email as this should be done for you as you authenticated through GitHub. You can check the git configuration with the command ``git config --list`` in your terminal

Also, since you started your workspace from your own SciPy fork your workspace will also have ``upstream`` and ``origin`` added as remotes. You can verify this by typing ``git remote`` on your terminal or by clicking on the "branch name" on the status bar (see image below).

..image:: _static/gitpod/scipy-gitpod-branches.png 
    :alt: Gitpod workspace with branches menu screenshot
    :figclass: align-center

From here on you can follow the workflow described in :ref:`development-workflow`. 

Rendering the documentation
----------------------------


.. _Gitpod: https://www.gitpod.io/
.. _Scipy repository on GitHub: https://github.com/scipy/scipy
.. _create your own fork: https://help.github.com/en/articles/fork-a-repo