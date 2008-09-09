import os
import shutil
import subprocess
from os.path import join as pjoin, split as psplit, dirname, exists as pexists
import re
from zipfile import ZipFile

def get_sdist_tarball(src_root):
    """Return the name of the installer built by sdist command."""
    # Yeah, the name logic is harcoded in distutils. We have to reproduce it
    # here
    name = "scipy-%s.zip" % get_scipy_version(src_root)
    return name

def build_sdist(chdir):
    cwd = os.getcwd()
    try:
        os.chdir(chdir)
        cmd = ["python", "setup.py", "sdist", "--format=zip"]
        subprocess.call(cmd)
    except Exception, e:
        raise RuntimeError("Error while executing cmd (%s)" % e)
    finally:
        os.chdir(cwd)

def prepare_scipy_sources(src_root, bootstrap = 'bootstrap'):
    zid = ZipFile(pjoin(src_root, 'dist', get_sdist_tarball(src_root)))
    root = 'scipy-%s' % get_scipy_version(src_root)

    # From the sdist-built tarball, extract all files into bootstrap directory,
    # but removing the numpy-VERSION head path
    for name in zid.namelist():
        cnt = zid.read(name)
        if name.startswith(root):
            # XXX: even on windows, the path sep in zip is '/' ?
            name = name.split('/', 1)[1]
        newname = pjoin(bootstrap, name)

        if not os.path.exists(dirname(newname)):
            os.makedirs(dirname(newname))
        fid = open(newname, 'wb')
        fid.write(cnt)

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

def prepare_nsis_script(bootstrap, pyver, numver):
    tpl = os.path.join('nsis_scripts', 'scipy-superinstaller.nsi.in')
    source = open(tpl, 'r')
    target = open(pjoin(bootstrap, 'scipy-superinstaller.nsi'), 'w')

    installer_name = 'scipy-%s-win32-superpack-python%s.exe' % (numver, pyver)
    cnt = "".join(source.readlines())
    cnt = cnt.replace('@SCIPY_INSTALLER_NAME@', installer_name)
    for arch in ['nosse', 'sse2', 'sse3']:
        cnt = cnt.replace('@%s_BINARY@' % arch.upper(),
                          get_binary_name(arch))

    target.write(cnt)

def get_binary_name(arch):
    return "scipy-%s-%s.exe" % (get_scipy_version(ROOT), arch)

def prepare_bootstrap(src_root, pyver):
    bootstrap = "bootstrap-%s" % pyver
    if os.path.exists(bootstrap):
        shutil.rmtree(bootstrap)
    os.makedirs(bootstrap)

    build_sdist(src_root)
    prepare_scipy_sources(src_root, bootstrap)

    shutil.copy('build.py', bootstrap)
    prepare_nsis_script(bootstrap, pyver, get_numpy_version())

if __name__ == '__main__':
    ROOT = os.path.join("..", "..", "..")
    prepare_bootstrap(ROOT, "2.5")
