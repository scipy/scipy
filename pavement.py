"""
This paver file is intended to help with the release process as much as
possible. It relies on virtualenv to generate 'bootstrap' environments as
independent from the user system as possible (e.g. to make sure the sphinx doc
is built against the built scipy, not an installed one).

The release is assumed to be done on OS X. See release.sh for a script that
employs the Paver tasks defined in this file, and builds everything required
for a release at once.


Building a Windows installer from Wine
======================================

The Python version the installer is for can be specified with the ``-p 2.6``
switch (this works for all installer tasks). To build a simple (no SSE
instructions) installer, do::

    paver bdist_wininst_simple

This assumes that blas/lapack are in c:\local\lib inside drive_c.  You will
have to make sure your Wine python locations (WINE_PYS) are configured
correctly.

The superpack requires all the Atlas libraries for every arch to be installed
(see SITECFG), and can then be built as follows::

    paver bdist_superpack


Building an installer for OS X
==============================

For a simple installer, which is just an mpkg inside a dmg, do::

  paver simple_dmg

For a more fancy installer which includes documentation and looks better, do::

  paver pdf  # needs to be done only once
  paver dmg


Building changelog + notes
==========================

Assumes you have git and the binaries/tarballs in installers/::

    paver write_release_and_log

This automatically put the checksum into NOTES.txt, and write the Changelog
which can be uploaded to Github Releases (and maybe sourceforge for historical
reasons, see gh-4939).


TODO
====
    - make it more easily customizable (through command line args)
    - missing targets: install & test, sdist test, debian packaging
    - fix bdist_mpkg: we build the same source twice -> how to make sure we use
      the same underlying python for egg install in venv and for bdist_mpkg
"""

import os
import sys
import subprocess
import re
import shutil
import warnings
from hashlib import md5
from hashlib import sha256

import distutils

try:
    from paver.tasks import VERSION as _PVER
    if not _PVER >= '1.0':
        raise RuntimeError("paver version >= 1.0 required (was %s)" % _PVER)
except ImportError, e:
    raise RuntimeError("paver version >= 1.0 required")

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
RELEASE = 'doc/release/1.2.0-notes.rst'

# Start/end of the log (from git)
LOG_START = 'v1.1.0'
LOG_END = 'master'


#-------------------------------------------------------
# Hardcoded build/install dirs, virtualenv options, etc.
#-------------------------------------------------------

# Default python version
PYVER="2.7"

# Paver options object, holds all default dirs
options(bootstrap=Bunch(bootstrap_dir="bootstrap"),
        virtualenv=Bunch(packages_to_install=["sphinx==1.1.3", "numpydoc"],
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

# Where we can find BLAS/LAPACK/ATLAS on Windows/Wine
SITECFG = {"sse3" : {'BLAS': 'None', 'LAPACK': 'None',
                     'ATLAS': r'C:\local\lib\atlas\sse3'},
           "sse2" : {'BLAS': 'None', 'LAPACK': 'None',
                     'ATLAS': r'C:\local\lib\atlas\sse2'},
           "nosse" : {'ATLAS': 'None', 'BLAS': r'C:\local\lib\atlas\nosse',
                      'LAPACK': r'C:\local\lib\atlas\nosse'}}

# Wine config for win32 builds
if sys.platform == "win32":
    WINE_PY26 = [r"C:\Python26\python.exe"]
    WINE_PY27 = [r"C:\Python27\python.exe"]
    WINE_PY32 = [r"C:\Python32\python.exe"]
    WINE_PY33 = [r"C:\Python33\python.exe"]
    WINE_PY34 = [r"C:\Python34\python.exe"]
    WINDOWS_ENV = os.environ
    MAKENSIS = ["makensis"]
elif sys.platform == "darwin":
    WINE_PY26 = ["wine", os.environ['HOME'] + "/.wine/drive_c/Python26/python.exe"]
    WINE_PY27 = ["wine", os.environ['HOME'] + "/.wine/drive_c/Python27/python.exe"]
    WINE_PY32 = ["wine", os.environ['HOME'] + "/.wine/drive_c/Python32/python.exe"]
    WINE_PY33 = ["wine", os.environ['HOME'] + "/.wine/drive_c/Python33/python.exe"]
    WINE_PY34 = ["wine", os.environ['HOME'] + "/.wine/drive_c/Python34/python.exe"]
    WINDOWS_ENV = os.environ
    WINDOWS_ENV["DYLD_FALLBACK_LIBRARY_PATH"] = "/usr/X11/lib:/usr/lib"
    MAKENSIS = ["wine", "makensis"]
else:
    WINE_PY26 = [os.environ['HOME'] + "/.wine/drive_c/Python26/python.exe"]
    WINE_PY27 = [os.environ['HOME'] + "/.wine/drive_c/Python27/python.exe"]
    WINE_PY32 = [os.environ['HOME'] + "/.wine/drive_c/Python32/python.exe"]
    WINE_PY33 = [os.environ['HOME'] + "/.wine/drive_c/Python33/python.exe"]
    WINE_PY34 = [os.environ['HOME'] + "/.wine/drive_c/Python34/python.exe"]
    WINDOWS_ENV = os.environ
    MAKENSIS = ["wine", "makensis"]
WINE_PYS = {'3.4':WINE_PY34, '3.3':WINE_PY33, '3.2':WINE_PY32,
            '2.7':WINE_PY27, '2.6':WINE_PY26}

# Framework Python locations on OS X
MPKG_PYTHON = {
        "2.6": "/Library/Frameworks/Python.framework/Versions/2.6/bin/python",
        "2.7": "/Library/Frameworks/Python.framework/Versions/2.7/bin/python",
        "3.2": "/Library/Frameworks/Python.framework/Versions/3.2/bin/python3",
        "3.3": "/Library/Frameworks/Python.framework/Versions/3.3/bin/python3",
        "3.4": "/Library/Frameworks/Python.framework/Versions/3.4/bin/python3"
        }
# Full path to the *static* gfortran runtime
LIBGFORTRAN_A_PATH = "/usr/local/lib/libgfortran.a"


#--------------------------------------
# Utility functions and bootstrap stuff
#--------------------------------------

def parse_numpy_version(pyexec):
    if isinstance(pyexec, str):
        cmd = [pyexec, "-c", "'import numpy; print(numpy.version.version)'"]
    else:
        # sequence for pyexec
        cmd = pyexec + ["-c", "'import numpy; print(numpy.version.version)'"]

    # Execute in shell because launching python from python does not work
    # (hangs)
    p = subprocess.Popen(" ".join(cmd), stdout=subprocess.PIPE, shell=True)
    out = p.communicate()[0]
    if p.returncode:
        raise RuntimeError("Command %s failed" % " ".join(cmd))

    a = re.compile("^([0-9]+)\.([0-9]+)\.([0-9]+)")
    if a:
        return tuple([int(i) for i in a.match(out).groups()[:3]])
    else:
        raise ValueError("Could not parse version (%s)" % out)

@task
def bootstrap():
    """create virtualenv in ./install"""
    try:
        import virtualenv
    except ImportError:
        raise RuntimeError("virtualenv is needed for bootstrap")

    bdir = options.bootstrap_dir
    if not os.path.exists(bdir):
        os.makedirs(bdir)
    bscript = "boostrap.py"

    options.virtualenv.script_name = os.path.join(options.bootstrap_dir,
                                                  bscript)
    options.bootstrap.no_site_packages = False
    call_task('paver.virtual.bootstrap')
    sh('cd %s; %s %s' % (bdir, sys.executable, bscript))

@task
def clean():
    """Remove build, dist, egg-info garbage."""
    d = ['build', 'dist', 'scipy.egg-info']
    for i in d:
        if os.path.exists(i):
            shutil.rmtree(i)

    bdir = os.path.join('doc', options.sphinx.builddir)
    if os.path.exists(bdir):
        shutil.rmtree(bdir)

@task
def clean_bootstrap():
    bdir = os.path.join(options.bootstrap.bootstrap_dir)
    if os.path.exists(bdir):
        shutil.rmtree(bdir)

@task
@needs('clean', 'clean_bootstrap')
def nuke():
    """Remove everything: build dir, installers, bootstrap dirs, etc..."""
    for d in [options.superpack.builddir, options.installers.releasedir]:
        if os.path.exists(d):
            shutil.rmtree(d)

#--------------------
# Documentation tasks
#--------------------

@task
def html(options):
    """Build scipy documentation and put it into build/docs"""
    # Don't use paver html target because of scipy bootstrapping problems
    subprocess.check_call(["make", "html"], cwd="doc")
    builtdocs = paver.path.path("doc") / options.sphinx.builddir / "html"
    options.html.builddir.rmtree()
    builtdocs.copytree(options.html.builddir)

@task
def latex():
    """Build scipy documentation in latex format."""
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
    # First clean the repo and update submodules (for up-to-date doc html theme
    # and Sphinx extensions)
    sh('git clean -xdf')
    sh('git submodule init')
    sh('git submodule update')

    # Fix file permissions
    sh('chmod a+rX -R *')

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
    """sdists, release notes and changelog.  Docs and wheels are built in
    separate steps (see doc/source/dev/releasing.rst).
    """
    # Source tarballs
    sdist()

    # README (gpg signed) and Changelog
    write_release_and_log()


#---------------------------------------
# Windows installers (Wine-based builds)
#---------------------------------------

def superpack_name(pyver, numver):
    """Return the filename of the superpack installer."""
    return 'scipy-%s-win32-superpack-python%s.exe' % (numver, pyver)

def internal_wininst_name(arch):
    """Return the name of the wininst as it will be inside the superpack (i.e.
    with the arch encoded."""
    ext = '.exe'
    return "scipy-%s-%s%s" % (FULLVERSION, arch, ext)

def wininst_name(pyver):
    """Return the name of the installer built by wininst command."""
    ext = '.exe'
    return "scipy-%s.win32-py%s%s" % (FULLVERSION, pyver, ext)

def bdist_wininst_arch(pyver, arch):
    """Arch specific wininst build."""
    if os.path.exists("build"):
        shutil.rmtree("build")
    _bdist_wininst(pyver, SITECFG[arch])

def prepare_nsis_script(pyver, numver):
    if not os.path.exists(options.superpack.builddir):
        os.makedirs(options.superpack.builddir)

    tpl = os.path.join('tools/win32/build_scripts/nsis_scripts', 'scipy-superinstaller.nsi.in')
    source = open(tpl, 'r')
    target = open(os.path.join(options.superpack.builddir, 'scipy-superinstaller.nsi'), 'w')

    installer_name = superpack_name(pyver, numver)
    cnt = "".join(source.readlines())
    cnt = cnt.replace('@SCIPY_INSTALLER_NAME@', installer_name)
    for arch in ['nosse', 'sse2', 'sse3']:
        cnt = cnt.replace('@%s_BINARY@' % arch.upper(),
                          internal_wininst_name(arch))
    target.write(cnt)

@task
def bdist_wininst_nosse(options):
    """Build the nosse wininst installer."""
    bdist_wininst_arch(options.python_version, 'nosse')

@task
def bdist_wininst_sse2(options):
    """Build the sse2 wininst installer."""
    bdist_wininst_arch(options.python_version, 'sse2')

@task
def bdist_wininst_sse3(options):
    """Build the sse3 wininst installer."""
    bdist_wininst_arch(options.python_version, 'sse3')

@task
@cmdopts([("python-version=", "p", "python version")])
def bdist_superpack(options):
    """Build all arch specific wininst installers."""
    pyver = options.python_version
    def copy_bdist(arch):
        # Copy the wininst in dist into the release directory
        source = os.path.join('dist', wininst_name(pyver))
        target = os.path.join(options.superpack.bindir, internal_wininst_name(arch))
        if os.path.exists(target):
            os.remove(target)
        if not os.path.exists(os.path.dirname(target)):
            os.makedirs(os.path.dirname(target))

        try:
            os.rename(source, target)
        except OSError:
            # May be due to dev version having 'Unknown' in name, if git isn't
            # found.  This can be the case when compiling under Wine.
            ix = source.find('.dev0+') + 6
            source = source[:ix] + 'Unknown' + source[ix+7:]
            os.rename(source, target)

    bdist_wininst_arch(pyver, 'nosse')
    copy_bdist("nosse")
    bdist_wininst_arch(pyver, 'sse2')
    copy_bdist("sse2")
    bdist_wininst_arch(pyver, 'sse3')
    copy_bdist("sse3")

    prepare_nsis_script(pyver, FULLVERSION)
    subprocess.check_call(MAKENSIS + ['scipy-superinstaller.nsi'],
                          cwd=options.superpack.builddir)

    # Copy the superpack into installers dir
    if not os.path.exists(options.installers.installersdir):
        os.makedirs(options.installers.installersdir)

    source = os.path.join(options.superpack.builddir,
                          superpack_name(pyver, FULLVERSION))
    target = os.path.join(options.installers.installersdir,
                          superpack_name(pyver, FULLVERSION))
    shutil.copy(source, target)

@task
@cmdopts([('python_version=', 'p', 'Python version to build the installer against')])
def bdist_wininst_simple():
    """Simple wininst-based installer."""
    call_task("clean")
    env = os.environ.copy()
    for k, v in SITECFG['nosse'].items():
        env[k] = v
    _bdist_wininst(options.bdist_wininst_simple.python_version, env)

def _bdist_wininst(pyver, cfg_env=None):
    cmd = WINE_PYS[pyver] + ['setup.py', 'build', '-c', 'mingw32', 'bdist_wininst']
    if cfg_env:
        for k, v in WINDOWS_ENV.items():
            cfg_env[k] = v
    else:
        cfg_env = WINDOWS_ENV
    try:
        subprocess.check_call(cmd, env=cfg_env)
    except subprocess.CalledProcessError:
        # Too many open files to compile in one go, so re-run.
        print('RESTART WINDOWS BUILD.  See gh-2709.')
        subprocess.check_call(cmd, env=cfg_env)


#--------------------
# Mac OS X installers
#--------------------

def dmg_name(fullversion, pyver, osxver=None):
    """Return name for dmg installer.

    Notes
    -----
    Python 2.7 has two binaries, one for 10.3 (ppc, i386) and one for 10.6
    (i386, x86_64). All other Python versions at python.org at the moment
    have binaries for 10.3 only. The "macosx%s" part of the dmg name should
    correspond to the python.org naming scheme.
    """
    # assume that for the py2.7/osx10.6 build the deployment target is set
    # (should be done in the release script).
    if not osxver:
        osxver = os.environ.get('MACOSX_DEPLOYMENT_TARGET', '10.3')
    return "scipy-%s-py%s-python.org-macosx%s.dmg" % (fullversion, pyver,
                                                      osxver)

def macosx_version():
    if not sys.platform == 'darwin':
        raise ValueError("Not darwin ??")
    st = subprocess.Popen(["sw_vers"], stdout=subprocess.PIPE)
    out = st.stdout.readlines()
    ver = re.compile("ProductVersion:\s+([0-9]+)\.([0-9]+)\.([0-9]+)")
    for i in out:
        m = ver.match(i)
        if m:
            return m.groups()

def mpkg_name(pyver):
    maj, min = macosx_version()[:2]
    return "scipy-%s-py%s-macosx%s.%s.mpkg" % (FULLVERSION, pyver, maj, min)

def prepare_static_gfortran_runtime(d):
    if not os.path.exists(d):
        os.makedirs(d)
    shutil.copy(LIBGFORTRAN_A_PATH, d)

@task
@cmdopts([('python_version=', 'p', 'Python version to build the installer against')])
def bdist_mpkg():
    call_task("clean")

    try:
        pyver = options.bdist_mpkg.python_version
    except AttributeError:
        pyver = PYVER

    _build_mpkg(pyver)

def _build_mpkg(pyver):
    numver = parse_numpy_version(MPKG_PYTHON[pyver])
    numverstr = ".".join(["%i" % i for i in numver])
    if pyver < "3.3":
        if not numver == (1, 8, 2):
            raise ValueError("Scipy 0.19.x should be built against numpy "
                             "1.8.2, (detected %s) for Python >= 3.4" % numverstr)

    prepare_static_gfortran_runtime("build")
    # account for differences between Python 2.7.1 versions from python.org
    if os.environ.get('MACOSX_DEPLOYMENT_TARGET', None) == "10.6":
        ldflags = "-undefined dynamic_lookup -bundle -arch i386 -arch x86_64 -Wl,-search_paths_first"
    else:
        ldflags = "-undefined dynamic_lookup -bundle -arch i386 -arch ppc -Wl,-search_paths_first"
    ldflags += " -L%s" % os.path.join(os.path.dirname(__file__), "build")

    sh("LDFLAGS='%s' %s setup.py bdist_mpkg" % (ldflags, MPKG_PYTHON[pyver]))


@task
@cmdopts([("python_version=", "p", "python version")])
def dmg():
    try:
        pyver = options.dmg.python_version
    except Exception:
        pyver = PYVER
    idirs = options.installers.installersdir

    # Check if doc exists. If not, say so and quit.
    docpath = options.doc.destdir_pdf
    ref = os.path.join(docpath, "reference.pdf")
    if not os.path.exists(ref):
        warnings.warn("Docs need to be built first! \n%s not found." % docpath)

    dmg_n = dmg_name(FULLVERSION, pyver)
    dmg = paver.path.path('scipy-macosx-installer') / dmg_n
    if dmg.exists():
        dmg.remove()

    call_task("clean")
    _build_mpkg(pyver)

    macosx_installer_dir = "tools/scipy-macosx-installer"
    dmg = os.path.join(macosx_installer_dir, dmg_name(FULLVERSION, pyver))
    if os.path.exists(dmg):
        os.remove(dmg)

    # Clean the image source
    content = os.path.join(macosx_installer_dir, 'content')
    if os.path.exists(content):
        shutil.rmtree(content)
    os.makedirs(content)

    # Copy mpkg into image source
    mpkg_source = os.path.join("dist", mpkg_name(pyver))
    mpkg_target = os.path.join(content, "scipy-%s-py%s.mpkg" % (FULLVERSION, pyver))
    shutil.copytree(mpkg_source, mpkg_target)

    # Copy docs into image source
    pdf_docs = os.path.join(content, "Documentation")
    if os.path.exists(pdf_docs):
        shutil.rmtree(pdf_docs)
    os.makedirs(pdf_docs)
    shutil.copy(ref, os.path.join(pdf_docs, "reference.pdf"))

    # Build the dmg
    cmd = ["./new-create-dmg", "--pkgname", os.path.basename(mpkg_target),
        "--volname", "scipy", os.path.basename(dmg), "./content"]
    st = subprocess.check_call(cmd, cwd=macosx_installer_dir)

    source = dmg
    target = os.path.join(idirs, os.path.basename(dmg))
    if not os.path.exists(os.path.dirname(target)):
        os.makedirs(os.path.dirname(target))
    shutil.copy(source, target)


@task
@cmdopts([('python_version=', 'p', 'Python version to build the installer against')])
def simple_dmg():
    try:
        pyver = options.simple_dmg.python_version
    except AttributeError:
        pyver = PYVER

    src_dir = "dmg-source"

    # Clean the source dir
    if os.path.exists(src_dir):
        shutil.rmtree(src_dir)
    os.makedirs(src_dir)

    # Build the mpkg
    clean()
    _build_mpkg(pyver)

    # Build the dmg
    shutil.copytree(os.path.join("dist", mpkg_name(pyver)),
                    os.path.join(src_dir, mpkg_name(pyver)))
    _create_dmg(pyver, src_dir, "Scipy Universal %s" % FULLVERSION)

def _create_dmg(pyver, src_dir, volname=None):
    # Build the dmg
    image_name = dmg_name(FULLVERSION, pyver)
    image = paver.path.path(image_name)
    image.remove()
    cmd = ["hdiutil", "create", image_name, "-srcdir", src_dir]
    if volname:
        cmd.extend(["-volname", "'%s'" % volname])
    sh(" ".join(cmd))


#----------------------------
# Release notes and Changelog
#----------------------------

def compute_md5(idirs):
    released = paver.path.path(idirs).listdir()
    checksums = []
    for f in sorted(released):
        m = md5(open(f, 'r').read())
        checksums.append('%s  %s' % (m.hexdigest(), os.path.basename(f)))

    return checksums

def compute_sha256(idirs):
    # better checksum so gpg signed README.txt containing the sums can be used
    # to verify the binaries instead of signing all binaries
    released = paver.path.path(idirs).listdir()
    checksums = []
    for f in sorted(released):
        m = sha256(open(f, 'r').read())
        checksums.append('%s  %s' % (m.hexdigest(), os.path.basename(f)))

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

    # Sign release
    cmd = ['gpg', '--clearsign', '--armor']
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

    out = st.communicate()[0]
    a = open(filename, 'w')
    a.writelines(out)
    a.close()

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
