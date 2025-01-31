.. _making-a-release:

Making a SciPy release
======================

At the highest level, this is what the release manager does to release a new
SciPy version:

#. Propose a release schedule in the SciPy forum at https://discuss.scientific-python.org/.
#. Create the maintenance branch for the release.
#. Tag the release.
#. Build all release artifacts (sources, installers, docs).
#. Upload the release artifacts.
#. Announce the release.
#. Port relevant changes to release notes and build scripts to main.

In this guide we attempt to describe in detail how to perform each of the above
steps.  In addition to those steps, which have to be performed by the release
manager, here are descriptions of release-related activities and conventions of
interest:

- :ref:`backporting`
- :ref:`labels-and-milestones`
- :ref:`version-numbering`
- :ref:`supported-py-numpy-versions`
- :ref:`deprecations`


Proposing a release schedule
----------------------------
A typical release cycle looks like:

- Create the maintenance branch
- Release a beta version
- Release a "release candidate" (RC)
- If needed, release one or more new RCs
- Release the final version once there are no issues with the last release
  candidate

There's usually at least one week between each of the above steps.  Experience
shows that a cycle takes between 4 and 8 weeks for a new minor version.
Bug-fix versions don't need a beta or RC, and can be done much quicker.

Ideally the final release is identical to the last RC, however there may be
minor difference - it's up to the release manager to judge the risk of that.
Typically, if compiled code or complex pure Python code changes then a new RC
is needed, while a simple bug-fix that's backported from main doesn't require
a new RC.

To propose a schedule, send a list with estimated dates for branching and
beta/rc/final releases to the SciPy forum at
https://discuss.scientific-python.org/. In the same message, ask everyone to check
if there are important issues/PRs that need to be included and aren't tagged
with the Milestone for the release or the "backport-candidate" label.


Creating the maintenance branch
-------------------------------
Before branching, ensure that the release notes are updated as far as possible.
Include the output of ``tools/gh_lists.py`` and ``tools/authors.py`` in the
release notes.

Maintenance branches are named ``maintenance/<major>.<minor>.x`` (e.g. 0.19.x).
To create one, simply push a branch with the correct name to the scipy repo.
Immediately after, push a commit where you increment the version number on the
main branch and add release notes for that new version.  Send an email to
the SciPy forum at https://discuss.scientific-python.org/ to let people know
that you've done this.


Updating the version switcher
------------------------------
The version switcher dropdown needs to be updated with the new release
information on the ``main`` branch only.

- ``doc/source/_static/version_switcher.json``: add the new release,
  the new development version, and transfer ``"preferred": true`` from the old release
  to the new release.


Updating upper bounds of dependencies
-------------------------------------
In main we do not set upper bounds, because we want to test new releases or
development versions of dependencies there. In a maintenance branch however,
the goal is to be able to create releases that stay working for years. Hence
correct upper bounds must be set. The following places must be updated after
creating a maintenance branch:

- ``pyproject.toml``: all build-time dependencies, as well as  supported Python
                      and NumPy versions
- ``scipy/__init__.py``: for NumPy version check

Each file has comments describing how to set the correct upper bounds.


Tagging a release
-----------------
First ensure that you have set up GPG correctly.  See
https://github.com/scipy/scipy/issues/4919 for a discussion of signing release
tags, and https://keyring.debian.org/creating-key.html for instructions on
creating a GPG key if you do not have one. Note that on some platforms
it may be more suitable to use ``gpg2`` instead of ``gpg`` so that
passwords may be stored by ``gpg-agent`` as discussed in
https://github.com/scipy/scipy/issues/10189.
When preparing a release remotely, it may be necessary to set
``pinentry-mode loopback`` in ``~/.gnupg/gpg-agent.conf`` because
use of ``gpg2`` will otherwise proceed via an inaccessible graphical
password prompt.

To make your key more readily identifiable as you, consider sending your key
to public key servers, with a command such as::

    gpg --send-keys <yourkeyid>

Check that all relevant commits are in the branch.  In particular, check issues
and PRs under the Milestone for the release
(https://github.com/scipy/scipy/milestones), PRs labeled "backport-candidate",
and that the release notes are up-to-date and included in the html docs.

Then update the ``version`` in ``pyproject.toml`` and commit the change with a
message like ``REL: set version to <version-number>``.  Don't push this commit
to the SciPy repo yet.

Finally tag the release locally with ``git tag -s <v1.x.y>`` (the ``-s`` ensures
the tag is signed). If ``gpg2`` is preferred, then
``git config --global gpg.program gpg2`` may be appropriate. Continue with
building release artifacts (next section).
Only push the release commit to the scipy repo once you have built the
sdists and docs successfully.  Then continue with building wheels.  Only push
the release tag to the repo once all wheels have been built successfully on
TravisCI and Appveyor (if it fails, you have to move the tag otherwise - which
is bad practice).  Finally, after pushing the tag, also push a second
commit which increments the version number and appends ``.dev0`` for
``version:`` and sets ``ISRELEASED`` to False again. This also applies with new
release candidates, and for removing the ``rc`` affix when switching from
release candidate to release proper.


Building release artifacts
--------------------------
Here is a complete list of artifacts created for a release:

- sdist (``scipy-x.y.y.tar.gz``, for both PyPI and GitHub Releases)
- Binary wheels for Windows, Linux and macOS
- Documentation (html)
- A ``README.txt`` file
- A ``Changelog`` file

An ``sdist`` is generated by running ``python -m build --sdist`` (note: we
still need to move this into a CI job!), and the Changelog and README are built
by running ``python dev.py notes`` (with tags, see ``python dev.py notes
--help``) in the repo root, and end up in ``REPO_ROOT/release/``.  Do this
after you've created the signed tag locally.  If this completes without issues,
push the release commit (not the tag, see section above) to the scipy repo.

To build wheels, push a commit containing the text ``[wheel build]`` to the branch
used for the current release. This triggers ``cibuildwheel`` builds for all needed
Python versions and platforms. The appropriate version pins for NumPy and other
dependencies should have been updated in ``pyproject.toml`` just after branching.
If the wheel builds reveal issues that need to be fixed
with backports on the maintenance branch, you may remove the local tags (for example
``git tag -d v1.2.0rc1``) and restart with tagging above on the new candidate commit.

The ``cibuildwheel`` infrastructure runs the tests from the built wheels and if they pass,
uploads the wheels to https://anaconda.org/multibuild-wheels-staging/scipy.

From there you can download them for uploading to PyPI. This can be
done in an automated fashion using ``tools/download-wheels.py``::

  $ python tools/download-wheels.py 1.5.0rc1 -w REPO_ROOT/release/installers

After this, we want to regenerate the README file, in order to have the MD5 and
SHA256 checksums of the just downloaded wheels in it.  Run ``python dev.py
notes`` again.


Uploading release artifacts
---------------------------
For a release there are currently five places on the web to upload things to:

- PyPI (sdist, wheels)
- GitHub Releases (sdist, release notes, Changelog)
- scipy.org (an announcement of the release)
- docs.scipy.org (html docs)

**PyPI:**

Upload first the wheels and then the sdist::

  twine upload REPO_ROOT/release/installers/*.whl
  twine upload REPO_ROOT/release/installers/scipy-1.x.y.tar.gz

**Github Releases:**

Use GUI on https://github.com/scipy/scipy/releases to create release and
upload all release artifacts. At this stage, it is appropriate to push
the tag and associate the new release (candidate) with this tag in the GUI.
For example, ``git push upstream v1.2.0rc1``, where ``upstream`` represents
``scipy/scipy``. It is useful to check a previous
release to determine exactly which artifacts should be included in the GUI
upload process. Also, note that the release notes are not automatically populated
into the release description on GitHub, and some manual reformatting to markdown
can be quite helpful to match the formatting of previous releases on the site.
We generally do not include Issue and Pull Request lists in these GUI
descriptions.

**scipy.org:**

Sources for the site are in https://github.com/scipy/scipy.org.
Update the News section in ``content/en/news.md`` via a PR.  This is only for
proper releases, not release candidates.

**docs.scipy.org:**

First build the scipy docs, by running ``make dist`` in ``scipy/doc/``.  Verify
that they look OK, then upload them to the doc server with
``make upload USERNAME=rgommers RELEASE=0.19.0``.  Note that SSH access to the
doc server is needed; ask @pv (server admin), @tylerjereddy or @rgommers (can
upload) if you don't have that.

The sources for the website itself are maintained in
https://github.com/scipy/docs.scipy.org/.  Add the new SciPy version in the
table of releases in ``index.rst``.  Push that commit, then do ``make upload
USERNAME=yourusername``. This is only for proper releases,
not release candidates.

.. note::

  The version switcher on the deployed docs will not consider the new release
  as the stable release until the changes to ``version_switcher.json`` have
  been merged into the ``main`` branch
  (see `gh-22305 <https://github.com/scipy/scipy/pull/22305>`__). Hence forward
  porting the release notes and version switcher changes to ``main`` should
  ideally be done at the same time as deploying the docs for the new release.

Wrapping up
-----------
Send a message announcing the release to
https://discuss.scientific-python.org/c/announcements/.

For beta and rc versions, ask people to test (run the scipy tests
and test against their own code) and report issues on Github or Discourse.

After the final release is done, port relevant changes to release notes, build
scripts, author name mapping in ``tools/authors.py`` and any other changes that
were only made on the maintenance branch to main.

Enable interactive examples by editing the runtime configuration file,
``try_examples.json``, in the root folder of the uploaded documentation on the
release server. One must remove the regular expression pattern ``".*"`` from
the ``ignore_patterns`` list.

.. code-block:: console

  $ ssh your-username@docs.scipy.org
  $ cd /srv/docs_scipy_org/doc/scipy-1.13.1
  $ vim try_examples.json  # edit the ignore list to remove: ".*"
