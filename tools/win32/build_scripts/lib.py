from os.path import join, exists
import re
import subprocess

VREGEX = re.compile("version\s*=\s*'(\d+)\.(\d+)\.(\d+)([a-zA-Z0-9]*)'")
ISRELREGEX = re.compile("release\s*=\s*True")
ISDEVREGEX = re.compile("release\s*=\s*False")
RCEX = re.compile("rc(\d)")

def get_svn_version(chdir):
    out = subprocess.Popen(['svn', 'info'],
                           stdout = subprocess.PIPE,
                           cwd = chdir).communicate()[0]
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
    version_file = join(src_root, "scipy", "version.py")
    if not exists(version_file):
        raise IOError("file %s not found" % version_file)

    fid = open(version_file, "r")
    version, rc, isdev = parse_verstring(fid.readlines())

    verstr = ".".join([str(i) for i in version])
    if rc > 0:
        verstr += "rc%d" % rc

    if isdev:
        verstr += ".dev"
        verstr += get_svn_version(src_root)
    return verstr

def parse_verstring(lines):
    isdev = None
    version = None
    for line in lines:
        m = VREGEX.match(line)
        if m:
            gr = m.groups()
            version = [int(i) for i in gr[:3]]
            m = RCEX.match(gr[3])
            if m:
                rc = int(m.group(1))
            else:
                rc = 0
            break

    if not version:
        raise ValueError("Error parsing %s" % "".join(lines))

    for line in lines:
        if ISRELREGEX.match(line):
            if isdev is None:
                isdev = False
            else:
                raise RuntimeError("isdev already set ?")
        if ISDEVREGEX.match(line):
            if isdev is None:
                isdev = True
            else:
                raise RuntimeError("isdev already set ?")

    return version, rc, isdev
