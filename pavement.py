"""
This paver file is intended to help with the release process, and build sdist,
documentation, release notes, and generate checksums for them.

For details on the release process, see
http://scipy.github.io/devdocs/dev/core-dev/index.html#making-a-scipy-release

Building changelog + notes
==========================

Assumes you have git and the binaries/tarballs in installers/::

    paver write_release_and_log

This automatically puts the checksum into NOTES.txt and writes the Changelog,
which can be uploaded to Github Releases.

"""

import os
import sys
import subprocess
import re
import shutil
import warnings
from hashlib import md5
from hashlib import sha256

try:
    from paver.tasks import VERSION as _PVER
    if not _PVER >= '1.0':
        raise RuntimeError("paver version >= 1.0 required (was %s)" % _PVER)
except ImportError as e:
    raise RuntimeError("paver version >= 1.0 required") from e

import paver
import paver.doctools
import paver.path
from paver.easy import options, Bunch, task, needs, dry, sh, call_task, cmdopts

sys.path.insert(0, os.path.dirname(__file__))
try:
    setup_py = __import__("setup")
    FULLVERSION = setup_py.VERSION
    # This is duplicated from setup.py
    if os.path.exists('.git'):
        GIT_REVISION = setup_py.git_version()
    else:
        GIT_REVISION = "Unknown"

    if not setup_py.ISRELEASED:
        if GIT_REVISION == "Unknown":
            FULLVERSION += '.dev0+Unknown'
        else:
            FULLVERSION += '.dev0+' + GIT_REVISION[:7]
finally:
    sys.path.pop(0)

try:
    # Ensure sensible file permissions
    os.umask(0o022)
except AttributeError:
    # No umask on non-posix
    pass


#-----------------------------------
# Things to be changed for a release
#-----------------------------------

# Source of the release notes
RELEASE = 'doc/release/1.7.0-notes.rst'

# Start/end of the log (from git)
LOG_START = 'v1.6.0'
LOG_END = 'master'


#-------------------------------------------------------
# Hardcoded build/install dirs, virtualenv options, etc.
#-------------------------------------------------------

# Default Python version
PYVER="3.6"

# Paver options object, holds all default dirs
options(bootstrap=Bunch(bootstrap_dir="bootstrap"),
        virtualenv=Bunch(packages_to_install=["sphinx==1.8.5", "numpydoc"],
                         no_site_packages=False),
        sphinx=Bunch(builddir="build", sourcedir="source", docroot='doc'),
        superpack=Bunch(builddir="build-superpack",
                        bindir=os.path.join("build-superpack","binaries")),
        installers=Bunch(releasedir="release",
                         installersdir=os.path.join("release", "installers")),
        doc=Bunch(doc_root="doc",
            sdir=os.path.join("doc", "source"),
            bdir=os.path.join("doc", "build"),
            bdir_latex=os.path.join("doc", "build", "latex"),
            destdir_pdf=os.path.join("build_doc", "pdf")),
        html=Bunch(builddir=os.path.join("build", "html")),
        dmg=Bunch(python_version=PYVER),
        bdist_wininst_simple=Bunch(python_version=PYVER),)


#--------------------
# Documentation tasks
#--------------------

@task
def html(options):
    """Build SciPy documentation and put it into build/docs"""
    # Don't use paver html target because of scipy bootstrapping problems
    subprocess.check_call(["make", "html"], cwd="doc")
    builtdocs = paver.path.path("doc") / options.sphinx.builddir / "html"
    options.html.builddir.rmtree()
    builtdocs.copytree(options.html.builddir)

@task
def latex():
    """Build SciPy documentation in latex format."""
    subprocess.check_call(["make", "latex"], cwd="doc")

@task
@needs('latex')
def pdf():
    bdir_latex = options.doc.bdir_latex
    destdir_pdf = options.doc.destdir_pdf

    def build_pdf():
        subprocess.check_call(["make", "all-pdf"], cwd=str(bdir_latex))
    dry("Build pdf doc", build_pdf)

    if os.path.exists(destdir_pdf):
        shutil.rmtree(destdir_pdf)
    os.makedirs(destdir_pdf)

    ref = os.path.join(bdir_latex, "scipy-ref.pdf")
    shutil.copy(ref, os.path.join(destdir_pdf, "reference.pdf"))

def tarball_name(type='gztar'):
    root = 'scipy-%s' % FULLVERSION
    if type == 'gztar':
        return root + '.tar.gz'
    elif type == 'xztar':
        return root + '.tar.xz'
    elif type == 'tar':
        return root + '.tar'
    elif type == 'zip':
        return root + '.zip'
    raise ValueError("Unknown type %s" % type)

@task
def sdist():
    # First, clean the repo and update submodules (for up-to-date doc html theme
    # and Sphinx extensions)
    sh('git clean -xdf')
    sh('git submodule init')
    sh('git submodule update')

    # Fix file permissions
    sh('chmod -R a+rX *')

    # To be sure to bypass paver when building sdist... paver + scipy.distutils
    # do not play well together.
    # Cython is run over all Cython files in setup.py, so generated C files
    # will be included.
    sh('python setup.py sdist --formats=gztar,zip')
    sh('python setup.py sdist --formats=tar')
    if os.path.exists(os.path.join('dist', tarball_name("xztar"))):
        os.unlink(os.path.join('dist', tarball_name("xztar")))
    sh('xz %s' % os.path.join('dist', tarball_name("tar")), ignore_error=True)

    # Copy the sdists into installers dir
    if not os.path.exists(options.installers.installersdir):
        os.makedirs(options.installers.installersdir)

    if not os.path.exists(os.path.join('dist', tarball_name("xztar"))):
        warnings.warn("Could not create tar.xz! Do you have xz installed?")
    else:
        t = 'xztar'
        source = os.path.join('dist', tarball_name(t))
        target = os.path.join(options.installers.installersdir, tarball_name(t))
        shutil.copy(source, target)

    for t in ['gztar', 'zip']:
        source = os.path.join('dist', tarball_name(t))
        target = os.path.join(options.installers.installersdir, tarball_name(t))
        shutil.copy(source, target)

@task
def release(options):
    """sdists, release notes and changelog. Docs and wheels are built in
    separate steps (see doc/source/dev/releasing.rst).
    """
    # Source tarballs
    sdist()

    # README (gpg signed) and Changelog
    write_release_and_log()


#----------------------------
# Release notes and Changelog
#----------------------------

def compute_md5(idirs):
    released = paver.path.path(idirs).listdir()
    checksums = []
    for fn in sorted(released):
        with open(fn, 'rb') as f:
            m = md5(f.read())
        checksums.append('%s  %s' % (m.hexdigest(), os.path.basename(fn)))

    return checksums

def compute_sha256(idirs):
    # better checksum so gpg signed README.txt containing the sums can be used
    # to verify the binaries instead of signing all binaries
    released = paver.path.path(idirs).listdir()
    checksums = []
    for fn in sorted(released):
        with open(fn, 'rb') as f:
            m = sha256(f.read())
        checksums.append('%s  %s' % (m.hexdigest(), os.path.basename(fn)))

    return checksums

def write_release_task(options, filename='NOTES.txt'):
    idirs = options.installers.installersdir
    source = paver.path.path(RELEASE)
    target = paver.path.path(filename)
    if target.exists():
        target.remove()

    tmp_target = paver.path.path(filename + '.tmp')
    source.copy(tmp_target)

    with open(str(tmp_target), 'a') as ftarget:
        ftarget.writelines("""
Checksums
=========

MD5
~~~

""")
        ftarget.writelines(['%s\n' % c for c in compute_md5(idirs)])
        ftarget.writelines("""
SHA256
~~~~~~

""")
        ftarget.writelines(['%s\n' % c for c in compute_sha256(idirs)])

    # Sign release; on some platforms gpg2 may actually
    # be named gpg
    cmd = ['gpg2', '--clearsign', '--armor']
    if hasattr(options, 'gpg_key'):
        cmd += ['--default-key', options.gpg_key]
    cmd += ['--output', str(target), str(tmp_target)]
    subprocess.check_call(cmd)
    print("signed %s" % (target,))
    tmp_target.remove()


def write_log_task(filename='Changelog'):
    st = subprocess.Popen(
            ['git', 'log',  '%s..%s' % (LOG_START, LOG_END)],
            stdout=subprocess.PIPE)

    out = st.communicate()[0].decode()
    with open(filename, 'w') as a:
        a.writelines(out)

@task
@cmdopts([('gpg_key=', 'g', 'GPG key to use for signing')])
def write_release(options):
    write_release_task(options)

@task
def write_log():
    write_log_task()

@task
@cmdopts([('gpg_key=', 'g', 'GPG key to use for signing')])
def write_release_and_log(options):
    write_release_task(options, os.path.join(options.installers.releasedir, 'README'))
    write_log_task(os.path.join(options.installers.releasedir, 'Changelog'))
