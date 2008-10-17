import os
from os.path import join as pjoin, normpath, exists as pexists, dirname
import subprocess
from shutil import rmtree
import re
from zipfile import ZipFile

SRC_ROOT = normpath(pjoin(os.getcwd(), os.pardir, os.pardir, os.pardir))
BUILD_ROOT = os.getcwd()

PYVER = 2.5

options(
    clean=Bunch(
        src_dir = SRC_ROOT,
        pyver = PYVER
    ),
    build_sdist=Bunch(
        src_dir = SRC_ROOT
    )
)

@task
def clean():
    # Clean sdist
    sdir = pjoin(options.src_dir, "dist")
    if pexists(sdir):
        rmtree(sdir)
    mani = pjoin(options.src_dir, "MANIFEST")
    if pexists(mani):
        os.remove(mani)

    # Clean bootstrap directory
    bdir = bootstrap_dir(options.pyver)
    if pexists(bdir):
        rmtree(bdir)

@task
def build_sdist():
    cmd = ["python", "setup.py", "sdist", "--format=zip"]
    st = subprocess.call(cmd, cwd=options.src_dir)

@task
@needs('build_sdist')
def bootstrap():
    prepare_scipy_sources(options.src_dir, bootstrap_dir(options.pyver))

# Helpers
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

def bootstrap_dir(pyver):
    return pjoin(BUILD_ROOT, "bootstrap-%s" % pyver)

def get_scipy_version(src_root):
    version_file = pjoin(src_root, "scipy", "version.py")
    if not pexists(version_file):
        raise IOError("file %s not found" % version_file)

    fid = open(version_file, "r")
    vregex = re.compile("version\s*=\s*'(\d+)\.(\d+)\.(\d+)'")
    isrelregex = re.compile("release\s*=\s*True")
    isdevregex = re.compile("release\s*=\s*False")
    isdev = None
    version = None
    for line in fid.readlines():
        m = vregex.match(line)
        if m:
            version = [int(i) for i in m.groups()]
        if isrelregex.match(line):
            if isdev is None:
                isdev = False
            else:
                raise RuntimeError("isdev already set ?")
        if isdevregex.match(line):
            if isdev is None:
                isdev = True
            else:
                raise RuntimeError("isdev already set ?")

    verstr = ".".join([str(i) for i in version])
    if isdev:
        verstr += ".dev"
        verstr += get_svn_version(src_root)
    return verstr

def get_svn_version(chdir):
    out = subprocess.Popen(['svn', 'info'],
                           stdout = subprocess.PIPE, cwd = chdir).communicate()[0]
    r = re.compile('Revision: ([0-9]+)')
    svnver = None
    for line in out.split('\n'):
        m = r.match(line)
        if m:
            svnver = m.group(1)

    if not svnver:
        raise ValueError("Error while parsing svn version ?")

    return svnver
