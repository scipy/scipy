.. _following-latest:

These are the instructions if you just want to follow the latest
*SciPy* source, but you don't need to do any development for now.
If you do want to contribute a patch (excellent!) or do more extensive
SciPy development, see :ref:`development-workflow`.

The steps are:

* :ref:`git-intro`
* get local copy of the git repository from Github_
* update local copy from time to time

Get the local copy of the code
==============================

From the command line::

   git clone git://github.com/scipy/scipy.git

You now have a copy of the code tree in the new ``scipy`` directory.
If this doesn't work you can try the alternative read-only url::

   git clone https://github.com/scipy/scipy.git

Updating the code
=================

From time to time you may want to pull down the latest code.  Do this with::

   cd scipy
   git fetch
   git merge --ff-only

The tree in ``scipy`` will now have the latest changes from the initial
repository.

.. _Github: https://github.com/scipy
