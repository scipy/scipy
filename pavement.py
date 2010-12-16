"""
This paver file is intented to help with the release process as much as
possible. It relies on virtualenv to generate 'bootstrap' environments as
independent from the user system as possible (e.g. to make sure the sphinx doc
is built against the built scipy, not an installed one).

Building a simple (no-superpack) windows installer from wine
============================================================

It assumes that blas/lapack are in c:\local\lib inside drive_c. Build python
2.5 and python 2.6 installers.

    paver bdist_wininst_simple

You will have to configure your wine python locations (WINE_PYS).

The superpack requires all the atlas libraries for every arch to be installed
(see SITECFG), and can then be built as follows::

    paver bdist_superpack -p 2.6

Building changelog + notes
==========================

Assumes you have git and the binaries/tarballs in installers/::

    paver write_release_and_log

This automatically put the checksum into NOTES.txt, and write the Changelog
which can be uploaded to sourceforge.

TODO
====
    - the script is messy, lots of global variables
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
try:
    from hash import md5
except ImportError:
    import md5

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
    FULLVERSION = setup_py.FULLVERSION
finally:
    sys.path.pop(0)

# Default python version
PYVER="2.6"
DMG_DIR = "dmg-source"


SSE3_CFG = {'BLAS': 'None', 'LAPACK': 'None', 'ATLAS': r'C:\local\lib\yop\sse3'}
SSE2_CFG = {'BLAS': 'None', 'LAPACK': 'None', 'ATLAS': r'C:\local\lib\yop\sse2'}
NOSSE_CFG = {'ATLAS': 'None', 'BLAS': r'C:\local\lib\yop\nosse', 'LAPACK': r'C:\local\lib\yop\nosse'}
SITECFG = {"sse2" : SSE2_CFG, "sse3" : SSE3_CFG, "nosse" : NOSSE_CFG}

# Wine config for win32 builds
if sys.platform == "win32":
    WINE_PY25 = [r"C:\Python25\python.exe"]
    WINE_PY26 = [r"C:\Python26\python26.exe"]
    WINE_PY27 = [r"C:\Python27\python27.exe"]
    WINE_PY31 = [r"C:\Python31\python.exe"]
    WINDOWS_ENV = os.environ
    MAKENSIS = ["makensis"]
elif sys.platform == "darwin":
    WINE_PY25 = ["wine", os.environ['HOME'] + "/.wine/drive_c/Python25/python.exe"]
    WINE_PY26 = ["wine", os.environ['HOME'] + "/.wine/drive_c/Python26/python.exe"]
    WINE_PY27 = ["wine", os.environ['HOME'] + "/.wine/drive_c/Python27/python.exe"]
    WINE_PY31 = ["wine", os.environ['HOME'] + "/.wine/drive_c/Python31/python.exe"]
    WINDOWS_ENV = os.environ
    WINDOWS_ENV["DYLD_FALLBACK_LIBRARY_PATH"] = "/usr/X11/lib:/usr/lib"
    MAKENSIS = ["wine", "makensis"]
else:
    WINE_PY25 = [os.environ['HOME'] + "/.wine/drive_c/Python25/python.exe"]
    WINE_PY26 = [os.environ['HOME'] + "/.wine/drive_c/Python26/python.exe"]
    WINE_PY27 = [os.environ['HOME'] + "/.wine/drive_c/Python27/python.exe"]
    WINE_PY31 = [os.environ['HOME'] + "/.wine/drive_c/Python31/python.exe"],
    WINDOWS_ENV = os.environ
    MAKENSIS = ["wine", "makensis"]
WINE_PYS = {'3.1':WINE_PY31, '2.7':WINE_PY27, '2.6':WINE_PY26, '2.5':WINE_PY25}
SUPERPACK_BUILD = 'build-superpack'
SUPERPACK_BINDIR = os.path.join(SUPERPACK_BUILD, 'binaries')

# XXX: fix this in a sane way
MPKG_PYTHON = {"2.5": "/Library/Frameworks/Python.framework/Versions/2.5/bin/python",
        "2.6": "/Library/Frameworks/Python.framework/Versions/2.6/bin/python",
        "2.7": "/Library/Frameworks/Python.framework/Versions/2.7/bin/python",
        "3.1": "/Library/Frameworks/Python.framework/Versions/3.1/bin/python3"}
# Full path to the *static* gfortran runtime
LIBGFORTRAN_A_PATH = "/usr/local/lib/libgfortran.a"

# Where to put built documentation (where it will picked up for copy to
# binaries)
PDF_DESTDIR = paver.path.path('build') / 'pdf'
HTML_DESTDIR = paver.path.path('build') / 'html'
DOC_ROOT = paver.path.path("doc")
DOC_SRC = DOC_ROOT / "source"
DOC_BLD = DOC_ROOT / "build"
DOC_BLD_LATEX = DOC_BLD / "latex"

# Source of the release notes
RELEASE = 'doc/release/0.9.0-notes.rst'

# Start/end of the log (from git)
LOG_START = 'svn/tags/0.8.0'
LOG_END = 'svn/0.9.x'

# Virtualenv bootstrap stuff
BOOTSTRAP_DIR = "bootstrap"
BOOTSTRAP_PYEXEC = "%s/bin/python" % BOOTSTRAP_DIR
BOOTSTRAP_SCRIPT = "%s/bootstrap.py" % BOOTSTRAP_DIR

# Where to put the final installers, as put on sourceforge
RELEASE_DIR = 'release'
INSTALLERS_DIR = os.path.join(RELEASE_DIR, 'installers')


options(sphinx=Bunch(builddir="build", sourcedir="source", docroot='doc'),
        virtualenv=Bunch(script_name=BOOTSTRAP_SCRIPT,
        packages_to_install=["sphinx==1.0.4"]),
        wininst=Bunch(pyver=PYVER))

def parse_numpy_version(pyexec):
    if isinstance(pyexec, str):
        cmd = [pyexec, "-c", "'import numpy; print numpy.version.version'"]
    else:
        # sequence for pyexec
        cmd = pyexec + ["-c", "'import numpy; print numpy.version.version'"]

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

# Bootstrap stuff
@task
def bootstrap():
    """create virtualenv in ./install"""
    install = paver.path.path(BOOTSTRAP_DIR)
    if not install.exists():
        install.mkdir()
    call_task('paver.virtual.bootstrap')
    sh('cd %s; %s bootstrap.py' % (BOOTSTRAP_DIR, sys.executable))

@task
def clean():
    """Remove build, dist, egg-info garbage."""
    d = ['build', 'dist', 'scipy.egg-info']
    for i in d:
        paver.path.path(i).rmtree()

    (paver.path.path('doc') / options.sphinx.builddir).rmtree()

@task
def clean_bootstrap():
    paver.path.path('bootstrap').rmtree()

@task
@needs('clean', 'clean_bootstrap')
def nuke():
    """Remove everything: build dir, installers, bootstrap dirs, etc..."""
    d = [SUPERPACK_BUILD, INSTALLERS_DIR, DMG_DIR]
    for i in d:
        paver.path.path(i).rmtree()



#------------
# Doc tasks
#------------
@task
def html(options):
    """Build scipy documentation and put it into build/docs"""
    # Don't use paver html target because of scipy bootstrapping problems
    subprocess.check_call(["make", "html"], cwd="doc")
    builtdocs = paver.path.path("doc") / options.sphinx.builddir / "html"
    HTML_DESTDIR.rmtree()
    builtdocs.copytree(HTML_DESTDIR)

@task
def latex():
    """Build scipy documentation in latex format."""
    subprocess.check_call(["make", "latex"], cwd="doc")

@task
@needs('latex')
def pdf():
    def build_pdf():
        subprocess.check_call(["make", "all-pdf"], cwd=str(DOC_BLD_LATEX))
    dry("Build pdf doc", build_pdf)

    PDF_DESTDIR.rmtree()
    PDF_DESTDIR.makedirs()

    user = DOC_BLD_LATEX / "scipy-user.pdf"
    user.copy(PDF_DESTDIR / "userguide.pdf")
    ref =  DOC_BLD_LATEX / "scipy-ref.pdf"
    ref.copy(PDF_DESTDIR / "reference.pdf")

def tarball_name(type='gztar'):
    root = 'scipy-%s' % FULLVERSION
    if type == 'gztar':
        return root + '.tar.gz'
    elif type == 'zip':
        return root + '.zip'
    raise ValueError("Unknown type %s" % type)

@task
def sdist():
    # To be sure to bypass paver when building sdist... paver + scipy.distutils
    # do not play well together.
    sh('python setup.py sdist --formats=gztar,zip')

    # Copy the superpack into installers dir
    if not os.path.exists(INSTALLERS_DIR):
        os.makedirs(INSTALLERS_DIR)

    for t in ['gztar', 'zip']:
        source = os.path.join('dist', tarball_name(t))
        target = os.path.join(INSTALLERS_DIR, tarball_name(t))
        shutil.copy(source, target)

#------------------
# Wine-based builds
#------------------

def internal_wininst_name(arch, ismsi=False):
    """Return the name of the wininst as it will be inside the superpack (i.e.
    with the arch encoded."""
    if ismsi:
        ext = '.msi'
    else:
        ext = '.exe'
    return "scipy-%s-%s%s" % (FULLVERSION, arch, ext)

def wininst_name(pyver, ismsi=False):
    """Return the name of the installer built by wininst command."""
    # Yeah, the name logic is harcoded in distutils. We have to reproduce it
    # here
    if ismsi:
        ext = '.msi'
    else:
        ext = '.exe'
    name = "scipy-%s.win32-py%s%s" % (FULLVERSION, pyver, ext)
    return name

def bdist_wininst_arch(pyver, arch):
    """Arch specific wininst build."""
    if os.path.exists("build"):
        shutil.rmtree("build")

    _bdist_wininst(pyver, SITECFG[arch])

def superpack_name(pyver, numver):
    """Return the filename of the superpack installer."""
    return 'scipy-%s-win32-superpack-python%s.exe' % (numver, pyver)

def prepare_nsis_script(pyver, numver):
    if not os.path.exists(SUPERPACK_BUILD):
        os.makedirs(SUPERPACK_BUILD)

    tpl = os.path.join('tools/win32/build_scripts/nsis_scripts', 'scipy-superinstaller.nsi.in')
    source = open(tpl, 'r')
    target = open(os.path.join(SUPERPACK_BUILD, 'scipy-superinstaller.nsi'), 'w')

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
    bdist_wininst_arch(options.wininst.pyver, 'nosse')

@task
def bdist_wininst_sse2(options):
    """Build the sse2 wininst installer."""
    bdist_wininst_arch(options.wininst.pyver, 'sse2')

@task
def bdist_wininst_sse3(options):
    """Build the sse3 wininst installer."""
    bdist_wininst_arch(options.wininst.pyver, 'sse3')

@task
@cmdopts([("python-version=", "p", "python version")])
def bdist_superpack(options):
    """Build all arch specific wininst installers."""
    pyver = options.python_version
    def copy_bdist(arch):
        # Copy the wininst in dist into the release directory
        if int(pyver[0]) >= 3:
            source = os.path.join('build', 'py3k', 'dist', wininst_name(pyver))
        else:
            source = os.path.join('dist', wininst_name(pyver))
        target = os.path.join(SUPERPACK_BINDIR, internal_wininst_name(arch))
        if os.path.exists(target):
            os.remove(target)
        if not os.path.exists(os.path.dirname(target)):
            os.makedirs(os.path.dirname(target))
        os.rename(source, target)

    bdist_wininst_arch(pyver, 'nosse')
    copy_bdist("nosse")
    bdist_wininst_arch(pyver, 'sse2')
    copy_bdist("sse2")
    bdist_wininst_arch(pyver, 'sse3')
    copy_bdist("sse3")

    prepare_nsis_script(options.wininst.pyver, FULLVERSION)
    subprocess.check_call(MAKENSIS + ['scipy-superinstaller.nsi'],
                          cwd=SUPERPACK_BUILD)

    # Copy the superpack into installers dir
    if not os.path.exists(INSTALLERS_DIR):
        os.makedirs(INSTALLERS_DIR)

    source = os.path.join(SUPERPACK_BUILD, superpack_name(pyver, FULLVERSION))
    target = os.path.join(INSTALLERS_DIR, superpack_name(pyver, FULLVERSION))
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
    subprocess.check_call(cmd, env=cfg_env)


#-------------------
# Mac OS X installer
#-------------------
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
    return "scipy-%s-py%s-macosx%s.%s.mpkg" % \
            (FULLVERSION, pyver, maj, min)

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
    return "numpy-%s-py%s-python.org-macosx%s.dmg" % (fullversion, pyver,
                                                      osxver)

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
    if not numver == (1, 5, 1):
        raise ValueError("Scipy 0.9.x should be built against numpy 1.5.1, (detected %s)" % numverstr)

    prepare_static_gfortran_runtime("build")
    # account for differences between Python 2.7.1 versions from python.org
    if os.environ.get('MACOSX_DEPLOYMENT_TARGET', None) == "10.6":
        ldflags = "-undefined dynamic_lookup -bundle -arch i386 -arch x86_64 -Wl,-search_paths_first"
    else:
        ldflags = "-undefined dynamic_lookup -bundle -arch i386 -arch ppc -Wl,-search_paths_first"
    ldflags += " -L%s" % os.path.join(os.path.dirname(__file__), "build")

    if pyver == "2.5":
        sh("CC=gcc-4.0 LDFLAGS='%s' %s setupegg.py bdist_mpkg" % (ldflags, MPKG_PYTHON[pyver]))
    sh("LDFLAGS='%s' %s setupegg.py bdist_mpkg" % (ldflags, MPKG_PYTHON[pyver]))


# TODO: fix this!
@task
#@needs("pdf")
@cmdopts([("python-version=", "p", "python version")])
def dmg():
    try:
        pyver = options.dmg.python_version
    except:
        pyver = PYVER

    dmg_n = dmg_name(FULLVERSION, pyver)
    dmg = paver.path.path('scipy-macosx-installer') / dmg_n
    if dmg.exists():
        dmg.remove()

    call_task("clean")
    _build_mpkg(pyver)

    # FIXME: copy/adapt tools from numpy
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

    user = os.path.join(options.doc.destdir_pdf, "userguide.pdf")
    shutil.copy(user, os.path.join(pdf_docs, "userguide.pdf"))
    ref = os.path.join(options.doc.destdir_pdf, "reference.pdf")
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

    src_dir = DMG_DIR

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


#------------------------
# NOTES/Changelog stuff
#------------------------

def compute_md5():
    released = paver.path.path(INSTALLERS_DIR).listdir()
    checksums = []
    for f in released:
        if not f.endswith('DS_Store'):
            m = md5.md5(open(f, 'r').read())
            checksums.append('%s  %s' % (m.hexdigest(), f))

    return checksums

def write_release_task(filename='NOTES.txt'):
    source = paver.path.path(RELEASE)
    target = paver.path.path(filename)
    if target.exists():
        target.remove()
    source.copy(target)
    ftarget = open(str(target), 'a')
    ftarget.writelines("""
Checksums
=========

""")
    ftarget.writelines(['%s\n' % c for c in compute_md5()])


def write_log_task(filename='Changelog'):
    st = subprocess.Popen(
            ['git', 'svn', 'log',  '%s..%s' % (LOG_START, LOG_END)],
            stdout=subprocess.PIPE)

    out = st.communicate()[0]
    a = open(filename, 'w')
    a.writelines(out)
    a.close()

@task
def write_release():
    write_release_task()

@task
def write_log():
    write_log_task()

@task
def write_release_and_log():
    write_release_task(os.path.join(RELEASE_DIR, 'NOTES.txt'))
    write_log_task(os.path.join(RELEASE_DIR, 'Changelog'))
