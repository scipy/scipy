import os
from os.path import join as pjoin, normpath, exists as pexists, dirname
import subprocess
from shutil import rmtree, move as shmove
import re
from zipfile import ZipFile

from lib import get_svn_version, get_scipy_version

BUILD_MSI = False
SRC_ROOT = normpath(pjoin(os.getcwd(), os.pardir, os.pardir, os.pardir))
BUILD_ROOT = os.getcwd()

PYVER = '2.5'
ARCH = 'nosse'

PYEXECS = {"2.5": "C:\python25\python.exe",
        "2.4": "C:\python24\python24.exe",
        "2.3": "C:\python23\python23.exe"}

_SSE3_CFG = r"""[atlas]
library_dirs = C:\local\lib\yop\sse3"""
_SSE2_CFG = r"""[atlas]
library_dirs = C:\local\lib\yop\sse2"""
_NOSSE_CFG = r"""[atlas]
library_dirs = fakedirectorywhichhopefullydoesnotexist
[DEFAULT]
library_dirs = C:\local\lib\yop\nosse"""

SITECFG = {"sse2": _SSE2_CFG, "sse3": _SSE3_CFG, "nosse": _NOSSE_CFG}

options(
    clean=Bunch(
        src_dir = SRC_ROOT,
        pyver = PYVER
    ),
    clean_bootstrap=Bunch(
        src_dir = SRC_ROOT,
        pyver = PYVER
    ),
    build_sdist=Bunch(
        src_dir = SRC_ROOT
    ),
    build_binary=Bunch(
        pyver = PYVER,
        arch = ARCH,
        src_root = SRC_ROOT
    ),
    bootstrap=Bunch(
        pyver = PYVER,
        src_root = SRC_ROOT
    ),
    bootstrap_arch=Bunch(
        pyver = PYVER,
        arch = ARCH
    ),
    bootstrap_nsis=Bunch(
        pyver = PYVER,
        src_root = SRC_ROOT
    )
)

# Clean everything, including bootstrap source tree
@task
def clean():
    raw_clean(options.src_dir, options.pyver)

# Clean the bootstrap source tree for a clean build from scratch
@task
def clean_bootstrap():
    raw_clean_bootstrap(options.pyver)

@task
def build_sdist():
    raw_build_sdist(options.src_dir)

@task
@needs('build_sdist')
def bootstrap():
    raw_bootstrap(options.pyver, options.src_dir)

@task
def bootstrap_arch():
    pyver = options.pyver
    arch = options.arch
    set_bootstrap_sources(arch, pyver)

@task
def bootstrap_nsis():
    pyver = options.pyver
    bdir = bootstrap_dir(options.pyver)
    prepare_nsis_script(bdir, pyver, get_scipy_version(options.src_root))

@task
def build_binary():
    pyver = options.pyver
    arch = options.arch
    raw_build_arch(pyver, arch, options.src_root)

@task
@needs('bootstrap')
@needs('clean')
def build_nsis():
    scipy_verstr = get_scipy_version(options.src_root)
    bdir = bootstrap_dir(options.pyver)

    prepare_nsis_script(bdir, options.pyver, scipy_verstr)
    for arch in ['nosse', 'sse2', 'sse3']:
        raw_clean_bootstrap(options.pyver)
        set_bootstrap_sources(arch, options.pyver)
        raw_build_arch(options.pyver, arch, options.src_root)

    raw_build_nsis(options.pyver)

# Helpers
def set_bootstrap_sources(arch, pyver):
    bdir = bootstrap_dir(pyver)
    write_site_cfg(arch, cwd=bdir)

def get_sdist_tarball(src_root):
    """Return the name of the installer built by sdist command."""
    # Yeah, the name logic is harcoded in distutils. We have to reproduce it
    # here
    name = "scipy-%s.zip" % get_scipy_version(src_root)
    return name

def prepare_scipy_sources(src_root, bootstrap):
    zid = ZipFile(pjoin(src_root, 'dist', get_sdist_tarball(src_root)))
    root = 'scipy-%s' % get_scipy_version(src_root)

    # From the sdist-built tarball, extract all files into bootstrap directory,
    # but removing the scipy-VERSION head path
    for name in zid.namelist():
        cnt = zid.read(name)
        if name.startswith(root):
            # XXX: even on windows, the path sep in zip is '/' ?
            name = name.split('/', 1)[1]
        newname = pjoin(bootstrap, name)

        if not pexists(dirname(newname)):
            os.makedirs(dirname(newname))
        fid = open(newname, 'wb')
        fid.write(cnt)

def prepare_nsis_script(bdir, pyver, numver):
    tpl = pjoin('nsis_scripts', 'scipy-superinstaller.nsi.in')
    source = open(tpl, 'r')
    target = open(pjoin(bdir, 'scipy-superinstaller.nsi'), 'w')

    installer_name = 'scipy-%s-win32-superpack-python%s.exe' % (numver, pyver)
    cnt = "".join(source.readlines())
    cnt = cnt.replace('@SCIPY_INSTALLER_NAME@', installer_name)
    for arch in ['nosse', 'sse2', 'sse3']:
        cnt = cnt.replace('@%s_BINARY@' % arch.upper(),
                          get_binary_name(arch, numver))

    target.write(cnt)

def bootstrap_dir(pyver):
    return pjoin(BUILD_ROOT, "bootstrap-%s" % pyver)

def get_python_exec(ver):
    """Return the executable of python for the given version."""
    # XXX Check that the file actually exists
    try:
        return PYEXECS[ver]
    except KeyError:
        raise ValueError("Version %s not supported/recognized" % ver)

def write_site_cfg(arch, cwd=None):
    if not cwd:
        cwd = os.getcwd()

    scfg = pjoin(cwd, "site.cfg")
    if pexists(scfg):
        os.remove(scfg)
    f = open(scfg, 'w')
    f.writelines(SITECFG[arch])
    f.close()

def move_binary(arch, pyver, cwd, scipy_verstr):
    if not pexists(pjoin(cwd, "binaries")):
        os.makedirs(pjoin(cwd, "binaries"))

    shmove(pjoin(cwd, 'dist', get_windist_exec(pyver, scipy_verstr)),
           pjoin(cwd, 'binaries', get_binary_name(arch, scipy_verstr)))

def get_binary_name(arch, scipy_verstr):
    if BUILD_MSI:
        ext = '.msi'
    else:
        ext = '.exe'
    return "scipy-%s-%s%s" % (scipy_verstr, arch, ext)

def get_windist_exec(pyver, scipy_verstr):
    """Return the name of the installer built by wininst command."""
    # Yeah, the name logic is harcoded in distutils. We have to reproduce it
    # here
    if BUILD_MSI:
        ext = '.msi'
    else:
        ext = '.exe'
    name = "scipy-%s.win32-py%s%s" % (scipy_verstr, pyver, ext)
    return name

def raw_clean(src_dir, pyver):
    # Clean sdist
    sdir = pjoin(src_dir, "dist")
    if pexists(sdir):
        rmtree(sdir)
    mani = pjoin(src_dir, "MANIFEST")
    if pexists(mani):
        os.remove(mani)

    # Clean bootstrap directory
    bdir = bootstrap_dir(pyver)
    if pexists(bdir):
        rmtree(bdir)

def raw_clean_bootstrap(pyver):
    bdir = bootstrap_dir(pyver)
    for d in ["build", "dist"]:
        if pexists(pjoin(bdir, d)):
            rmtree(pjoin(bdir, d))

    if pexists(pjoin(bdir, "site.cfg")):
        os.remove(pjoin(bdir, "site.cfg"))

def raw_build_sdist(cwd):
    cmd = ["python", "setup.py", "sdist", "--format=zip"]

    build_log = "sdist.log"
    f = open(build_log, 'w')
    try:
        try:
            st = subprocess.call(cmd, #shell = True,
                            stderr = subprocess.STDOUT, stdout = f,
                            cwd=cwd)
            if st:
                raise RuntimeError("The cmd failed with status %d" % st)
        finally:
            f.close()
    except (subprocess.CalledProcessError, RuntimeError), e:
        print e
        msg = """
There was an error while executing the following command:

    %s

Error was : %s

Look at the log (%s).""" % (cmd, str(e), build_log)
        raise Exception(msg)

def raw_bootstrap(pyver, src_dir):
    bdir = bootstrap_dir(pyver)
    prepare_scipy_sources(src_dir, bdir)

def raw_build_arch(pyver, arch, src_root):
    scipy_verstr = get_scipy_version(src_root)
    bdir = bootstrap_dir(pyver)

    print "Building scipy (version %s) binary for python %s, arch is %s" % \
          (scipy_verstr, get_python_exec(pyver), arch)

    if BUILD_MSI:
        cmd = [get_python_exec(pyver), "setup.py", "build", "-c", "mingw32",
               "bdist_msi"]
    else:
        cmd = [get_python_exec(pyver), "setup.py", "build", "-c", "mingw32",
               "bdist_wininst"]
    build_log = "build-%s-%s.log" % (arch, pyver)
    f = open(build_log, 'w')

    try:
        try:
            st = subprocess.call(cmd, #shell = True,
                            stderr = subprocess.STDOUT, stdout = f,
                            cwd=bdir)
            if st:
                raise RuntimeError("The cmd failed with status %d" % st)
        finally:
            f.close()
    except (subprocess.CalledProcessError, RuntimeError), e:
        print e
        msg = """
There was an error while executing the following command:

    %s

Error was : %s

Look at the build log (%s).""" % (cmd, str(e), build_log)
        raise Exception(msg)

    move_binary(arch, pyver, bdir, scipy_verstr)

def raw_build_nsis(pyver):
    bdir = bootstrap_dir(options.pyver)
    st = subprocess.call(['makensis', 'scipy-superinstaller.nsi'],
                         cwd=bdir)
    if st:
        raise RuntimeError("Error while executing makensis command")
