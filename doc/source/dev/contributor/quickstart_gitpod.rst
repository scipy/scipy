:orphan:

.. _quickstart-gitpod:

=======================================================
Development environment quickstart guide (Gitpod)
=======================================================

This quick start guide covers:

*  using GitPod for your SciPy development environment
*  creating a personal fork of the SciPy repository on GitHub
*  a quick tour of Gitpod and VSCode
*  working on the SciPy documentation in Gitpod

Gitpod
-------

`Gitpod`_  is an open-source platform for automated and ready-to-code development environments. It enables developers to describe their dev environment as code and start instant and fresh development environments for each new task directly from your browser. This reduces the need to install local development environments and deal with incompatible dependencies.

Gitpod GitHub integration
--------------------------

To be able to use Gitpod, you will need to have the Gitpod app installed on your GitHub account, so if
you do not have an account yet, you will need to create one first.

Head over to the `Gitpod`_ website and click on the **Continue with GitHub** button. You will be redirected to the GitHub authentication page.
You will then be asked to install the `Gitpod GitHub app <https://github.com/marketplace/gitpod-io>`_.

Make sure to select **All repositories** access option to avoid issues with permissions later on. Click on the green **Install** button

.. image:: ../../_static/gitpod/installing-gitpod-io.png
    :alt: Gitpod application repository access and installation screenshot - All repositories checked (read and write access)

This will install the necessary hooks for the integration.

Forking the SciPy repository
-----------------------------

The best way to work on SciPy as a contributor is by making a fork of the repository first.

#. Browse to the `SciPy repository on GitHub`_ and `create your own fork`_.

#. Browse to your fork. Your fork will have a URL like https://github.com/andyfaff/scipy, except with your GitHub username in place of "andyfaff".

Starting Gitpod
----------------
Once you have authenticated to Gitpod through GitHub, you can install the `Gitpod browser extension <https://www.gitpod.io/docs/browser-extension>`_  which will add a **Gitpod** button next to the **Code** button in the repository:

.. image:: ../../_static/gitpod/scipy-github.png
    :alt: SciPy repository with Gitpod button screenshot

#. If you install the extension - you can click the **Gitpod** button to start a new workspace.
#. Alternatively, if you do not want to install the browser extension you can visit https://gitpod.io/#https://github.com/USERNAME/scipy replacing ``USERNAME`` with your GitHub username.

#. In both cases, this will open a new tab on your web browser and start building your development environment. Please note this can take a few minutes.

#. Once the build is complete, you will be directed to your workspace, including VSCode and all the dependencies you need to work on SciPy. The first time you start your workspace, you will notice that there might be some actions running. This will ensure that you have a development version of SciPy installed and that the docs are being pre-built for you.

#. Once the build is complete, you can test the build by entering::

        python dev.py

This command make use of Meson, for more details on how this affects your workflow migrating from ``distutils`` check the :ref:`meson` section in the Developer documentation.

   .. note::

      If you are a regular VSCode user you can install the Gitpod Remote extension on your local VSCode editor and sync your local editor to the remote Gitpod server. To do this click <kbd> F1 </kbd> and search for Gitpod open on VSCode.

Quick workspace tour
---------------------
Gitpod uses VSCode as the editor. If you have not used this editor before, you can check the Getting started `VSCode docs`_ to familiarise yourself with it.

Your workspace will look similar to the image below:

.. image:: ../../_static/gitpod/gitpod-workspace.png
    :alt: Gitpod workspace showing the VSCode editor - some parts are highlighted: 1- Python interpreter, 2- Git branch, 3- GitHub Pull Requests extension on the side bar 4- VsCode Marketplace 5- Current directory

.. note::
    By default VSCode initialises with a light theme, you can change to a dark theme by with the keyboard shortcut :kbd:`Cmd-K Cmd-T` in Mac or :kbd:`Ctrl-K Ctrl-T` in Linux and Windows.

We have marked some important sections in the editor:

#. Your current Python interpreter - by default, the ``scipy-dev`` is marked as ``base`` on the status bar, but it should be displayed as ``scipy-dev`` on your terminal. You do not need to activate the conda environment as this will always be activated for you.
#. Your current branch is always displayed in the status bar. You can also use this button to change or create branches.
#. GitHub Pull Requests extension - you can use this to work with Pull Requests from your workspace.
#. Marketplace extensions - we have added some essential extensions to the SciPy Gitpod. Still, you can also install other extensions or syntax highlighting themes for your user, and these will be preserved for you.
#. Your workspace directory - by default is ``/workspace/scipy`` **do not change this** as this is the only directory preserved in Gitpod.

We have also pre-installed a few tools and VSCode extensions to help with the development experience:

*  `GitHub CLI <https://cli.github.com/>`_
*  `VSCode rst extension <https://marketplace.visualstudio.com/items?itemName=lextudio.restructuredtext>`_
*  `VSCode Live server extension <https://marketplace.visualstudio.com/items?itemName=ritwickdey.LiveServer>`_
*  `VSCode Gitlens extension <https://marketplace.visualstudio.com/items?itemName=eamodio.gitlens>`_
*  `VSCode autodocstrings extension <https://marketplace.visualstudio.com/items?itemName=njpwerner.autodocstring>`_
*  `VSCode Git Graph extension <https://marketplace.visualstudio.com/items?itemName=mhutchie.git-graph>`_

Development workflow
-----------------------
The  :ref:`development-workflow` section of this documentation contains information regarding the SciPy development workflow. Make sure to check this before working on your contributions.

When using Gitpod, note these main differences with the setup described in :ref:`development-workflow`.

#. You **do not** need to configure your git username, and email as this should be done for you as you authenticated through GitHub. You can check the git configuration with the command ``git config --list`` in your terminal.
#. As you started your workspace from your own SciPy fork, you will by default have both "upstream" and "origin" added as remotes. You can verify this by typing ``git remote`` on your terminal or by clicking on the **branch name** on the status bar (see image below).

.. image:: ../../_static/gitpod/scipy-gitpod-branches.png
    :alt: Gitpod VSCode editor - git branches dropdown expanded

Rendering the SciPy documentation
----------------------------------
You can find the detailed documentation on how rendering the documentation with Sphinx works in the :ref:`rendering-documentation` section.

The documentation is pre-built during your workspace initialization. So once this task is completed, you have two main options to render the documentation in Gitpod.

Option 1: Using Liveserve
***************************

#. View the documentation in ``scipy/doc/build/html``. You can start with ``index.html`` and browse, or you can jump straight to the file you're interested in.
#. To see the rendered version of a page, you can right-click on the ``.html`` file and click on **Open with Live Serve**. Alternatively, you can open the file in the editor and click on the **Go live** button on the status bar.

    .. image:: ../../_static/gitpod/vscode-statusbar.png
        :alt: Gitpod VSCode editor - status bar zoom with "Go Live" tab highligthed by a teal rectangle

#. A simple browser will open to the right-hand side of the editor. We recommend closing it and click on the **Open in browser** button in the pop-up.
#. To stop the server click on the **Port: 5500** button on the status bar.

Option 2: Using the rst extension
***********************************

A quick and easy way to see live changes in a ``.rst`` file as you work on it uses the rst extension with docutils.

.. note::
    This will generate a simple live preview of the document without the ``html`` theme, and some backlinks might not be added correctly. But it is an easy and lightweight way to get instant feedback on your work.

#. Open any of the source documentation files located in ``doc/source`` in the editor.
#. Open VSCode Command Palette with :kbd:`Cmd-Shift-P` in Mac or :kbd:`Ctrl-Shift-P` in Linux and Windows. Start typing "restructured" and choose either "Open preview" or "Open preview to the Side".

    .. image:: ../../_static/gitpod/vscode-rst.png
        :alt: Gitpod VSCode editor - command palette expanded showing autocomplete options for "restruct"

#. As you work on the document, you will see a live rendering of it on the editor.

    .. image:: ../../_static/gitpod/rst-rendering.png
        :alt: Gitpod workspace - Quickstart Docker rst file opened on the left and rendered version on the right group of the editor.

If you want to see the final output with the ``html`` theme you will need to rebuild the docs with ``make html`` and use Live Serve as described in option 1.

FAQ's
-----

#. How long is my Gitpod workspace kept for?
    Your stopped workspace will be kept for 14 days and deleted afterwards if you do not use them.

#. Can I come back to a previous workspace?
    Yes, let's say you stepped away for a while and you want to carry on working on your SciPy contributions. You need to visit https://gitpod.io/workspaces and click on the workspace you want to spin up again. All your changes will be there as you last left them.

#. Can I install additional VSCode extensions?
    Absolutely! Any extensions you installed will be installed in your own workspace and preserved.

#. I registered on Gitpod but I still cannot see a **Gitpod** button in my repositories
    Head to https://gitpod.io/integrations and make sure you are logged in. Hover over GitHub and click on the three buttons that appear on the right. Click on edit permissions and make sure you have ``user:email``, ``read:user``, and ``public_repo`` checked.
    Click on **Update Permissions** and confirm the changes in the GitHub application page.

    .. image:: ../../_static/gitpod/gitpod-edit-permissions-gh.png
        :alt: Gitpod dashboard integrations section - edit GitHub permissions dropdown expanded

#. How long does my workspace stay active if I'm not using it?
    If you keep your workspace open in a browser tab but don't interact with it, it will shut down after 30 minutes. If you close the browser tab, it will shut down after 3 minutes.

.. _Gitpod: https://www.gitpod.io/
.. _SciPy repository on GitHub: https://github.com/scipy/scipy
.. _create your own fork: https://help.github.com/en/articles/fork-a-repo
.. _VSCode docs: https://code.visualstudio.com/docs/getstarted/tips-and-tricks


.. |br| raw:: html

    <br>
