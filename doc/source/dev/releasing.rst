Making a Scipy release
======================

At the highest level, this is what the release manager does to release a new
Scipy version:

#. Propose a release schedule on the scipy-dev mailing list.
#. Create the maintenance branch for the release.
#. Tag the release.
#. Build all release artifacts (sources, installers, docs).
#. Upload the release artifacts.
#. Announce the release.
#. Port relevant changes to release notes and build scripts to master.

In this guide we attempt to describe in detail how to perform each of the above
steps.  In addition to those steps, which have to be performed by the release
manager, here are descriptions of release-related activities and conventions of
interest:

- `Backporting`_
- `Labels and Milestones`_
- `Version numbering`_
- `Supported Python and Numpy versions`_
- `Deprecations`_


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
is needed, while a simple bug-fix that's backported from master doesn't require
a new RC.

To propose a schedule, send a list with estimated dates for branching and
beta/rc/final releases to scipy-dev. In the same email, ask everyone to check
if there are important issues/PRs that need to be included and aren't tagged
with the Milestone for the release or the "backport-candidate" label.


Creating the maintenance branch
-------------------------------
Before branching, ensure that the release notes are updated as far as possible.
Include the output of ``tools/gh_lists.py`` and ``tools/authors.py`` in the
release notes.

Maintenance branches are named ``maintenance/0.<minor-version>.x``.  To create
one, simply push a branch with the correct name to the scipy repo.  Immediately
after, push a commit where you increment the version number on the master
branch and add release notes for that new version.
Send an email to scipy-dev to let people know that you've done this.


Tagging a release
-----------------
First ensure that you have set up GPG correctly.  See
https://github.com/scipy/scipy/issues/4919 for a discussion of signing release
tags, and http://keyring.debian.org/creating-key.html for instructions on
creating a GPG key if you do not have one.

To make your key more readily identifiable as you, consider sending your key
to public keyservers, with a command such as::

    gpg --send-keys <yourkeyid>

Check that all relevant commits are in the branch.  In particular, check issues
and PRs under the Milestone for the release
(https://github.com/scipy/scipy/milestones), PRs labeled "backport-candidate",
and that the release notes are up-to-date and included in the html docs.

Then edit ``setup.py`` to get the correct version number (set
``ISRELEASED = True``) and commit it with a message like ``REL: set version to
<version-number>``.  Don't push this commit to the Scipy repo yet.

Finally tag the release locally with ``git tag -s <v0.x.y>`` (the ``-s`` ensures
the tag is signed).  Continue with building release artifacts (next section).
Only push the release commit and tag to the scipy repo once you have built the
docs and Windows installers successfully.  After that push, also push a second
commit which increment the version number and sets ``ISRELEASED`` to False
again.


Building release artifacts
--------------------------
Here is a complete list of artifacts created for a release:

- source archives (``.tar.gz``, ``.zip`` and ``.tar.xz``)
- Binary wheels, for OS X only at the moment
- Documentation (``html``, ``pdf``)
- A ``README`` file
- A ``Changelog`` file

All of these except the OS X wheel are built by running ``paver release`` in
the repo root.  Do this after you've created the signed tag.  If this completes
without issues, push the release tag to the scipy repo.  This is needed because
the OS X wheel build automatically builds the last tag.

To build wheels for OS X, push a commit to the master branch of
https://github.com/MacPython/scipy-wheels.  This triggers builds for all needed
Python versions on TravisCI.  Check in the ``.travis.yml`` config file what
version of Python and Numpy are used for the builds (it needs to be the lowest
supported Numpy version for each Python version).  See the README file in the
scipy-wheels repo for more details.

The TravisCI builds run the tests from the built wheels and if they pass upload
the wheels to http://wheels.scipy.org/.  From there you can copy them for
uploading to PyPI (see next section).


Uploading release artifacts
---------------------------
For a release there are currently five places on the web to upload things to:

- PyPI (tarballs, OS X wheels)
- Github releases (tarballs, release notes, Changelog)
- scipy.org (an announcement of the release)
- docs.scipy.org (html/pdf docs)

**PyPI:**

``twine upload -s <tarballs or wheels to upload>``

**Github Releases:**

Use GUI on https://github.com/scipy/scipy/releases to create release and
upload all release artifacts.

**SourceForge:**

The main download sites are PyPI and Github Releases.  Older releases are stored
on SourceForce ( http://sourceforge.net/projects/scipy/files/scipy).  That
download site has a "Latest" folder which redirects users to PyPI/GitHub, so
it's not needed to upload anything to SourceForge for new releases.

**scipy.org:**

Sources for the site are in https://github.com/scipy/scipy.org.
Update the News section in ``www/index.rst`` and then do
``make upload USERNAME=yourusername``.

**docs.scipy.org:**

First build the scipy docs, by running ``make dist`` in ``scipy/doc/``.  Verify
that they look OK, then upload them to the doc server with
``make upload USERNAME=rgommers RELEASE=0.17.0``.  Note that SSH access to the
doc server is needed; ask @pv (server admin) or @rgommers (can upload) if you
don't have that.

The sources for the website itself are maintained in
https://github.com/scipy/docs.scipy.org/.  Add the new Scipy version in the
table of releases in ``index.rst``.  Push that commit, then do ``make upload
USERNAME=yourusername``.


Wrapping up
-----------
Send an email announcing the release to the following mailing lists:

- scipy-dev
- scipy-user
- numpy-discussion
- python-announce (not for beta/rc releases)

For beta and rc versions, ask people in the email to test (run the scipy tests
and test against their own code) and report issues on Github or scipy-dev.

After the final release is done, port relevant changes to release notes, build
scripts, author name mapping in ``tools/authors.py`` and any other changes that
were only made on the maintenance branch to master.
