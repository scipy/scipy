import os
import sys
import glob
import shutil
import textwrap
import platform
import hashlib

from tempfile import mkstemp, gettempdir
from urllib.request import urlopen, Request
from urllib.error import HTTPError
import zipfile
import tarfile

OPENBLAS_V = 'v0.3.9'
OPENBLAS_LONG = 'v0.3.9'
BASE_LOC = ''
ANACONDA = 'https://anaconda.org/multibuild-wheels-staging/openblas-libs'
ARCHITECTURES = ['', 'windows', 'darwin', 'aarch64', 'x86', 'ppc64le', 's390x']
sha256_vals = {
'openblas64_-v0.3.9-macosx_10_9_x86_64-gf_1becaaa.tar.gz':
'53f606a7da75d390287f1c51b2af7866b8fe7553a26d2474f827daf0e5c8a886',
'openblas64_-v0.3.9-manylinux1_x86_64.tar.gz':
'6fe5b1e2a4baa16833724bcc94a80b22e9c99fc1b9a2ddbce4f1f82a8002d906',
'openblas64_-v0.3.9-win_amd64-gcc_7_1_0.zip':
'15d24a66c5b22cc7b3120e831658f491c7a063804c33813235044a6f8b56686d',
'openblas-v0.3.9-macosx_10_9_x86_64-gf_1becaaa.tar.gz': 
'8221397b9cfb8cb22f3efb7f228ef901e13f9fd89c7d7d0cb7b8a79b0610bf33',
'openblas-v0.3.9-manylinux1_i686.tar.gz': 
'31abf8eccb697a320a998ce0f59045edc964602f815d78690c5a23839819261c',
'openblas-v0.3.9-manylinux1_x86_64.tar.gz':
'd9c39acbafae9b1daef19c2738ec938109a59e9322f93eb9a3c50869d220deff',
'openblas-v0.3.9-win32-gcc_7_1_0.zip':
'69a7dc265e8a8e45b358637d11cb1710ce88c4456634c7ce37d429b1d9bc9aaa',
'openblas-v0.3.9-win_amd64-gcc_7_1_0.zip': 
'0cea06f4a2afebaa6255854f73f237802fc6b58eaeb1a8b1c22d87cc399e0d48',
'openblas-v0.3.9-manylinux2014_aarch64.tar.gz':
'10d5ef5e9e19af5c199b59a17f43763e0c85ecf13cbc8f2d91e076f7847cdb5e'
}

IS_32BIT = sys.maxsize < 2**32
def get_arch():
    if platform.system() == 'Windows':
        ret = 'windows'
    elif platform.system() == 'Darwin':
        ret = 'darwin'
    else:
        ret = platform.uname().machine
        # What do 32 bit machines report?
        # If they are a docker, they report x86_64 or i686
        if 'x86' in ret or ret == 'i686':
            ret = 'x86'
    assert ret in ARCHITECTURES
    return ret

def get_ilp64():
    if os.environ.get("NPY_USE_BLAS_ILP64", "0") == "0":
        return None
    if IS_32BIT:
        raise RuntimeError("NPY_USE_BLAS_ILP64 set on 32-bit arch")
    return "64_"

def download_openblas(target, arch, ilp64):
    fnsuffix = {None: "", "64_": "64_"}[ilp64]
    filename = ''
    if arch in ('aarch64', 'ppc64le', 's390x'):
        suffix = f'manylinux2014_{arch}.tar.gz'
        filename = f'{ANACONDA}/{OPENBLAS_LONG}/download/openblas{fnsuffix}-{OPENBLAS_LONG}-{suffix}'
        typ = 'tar.gz'
        typ = 'tar.gz'
    elif arch == 'darwin':
        suffix = 'macosx_10_9_x86_64-gf_1becaaa.tar.gz'
        filename = f'{ANACONDA}/{OPENBLAS_LONG}/download/openblas{fnsuffix}-{OPENBLAS_LONG}-{suffix}'
        typ = 'tar.gz'
    elif arch == 'windows':
        if IS_32BIT:
            suffix = 'win32-gcc_7_1_0.zip'
        else:
            suffix = 'win_amd64-gcc_7_1_0.zip'
        filename = f'{ANACONDA}/{OPENBLAS_LONG}/download/openblas{fnsuffix}-{OPENBLAS_LONG}-{suffix}'
        typ = 'zip'
    elif 'x86' in arch:
        if IS_32BIT:
            suffix = 'manylinux1_i686.tar.gz'
        else:
            suffix = 'manylinux1_x86_64.tar.gz'
        filename = f'{ANACONDA}/{OPENBLAS_LONG}/download/openblas{fnsuffix}-{OPENBLAS_LONG}-{suffix}'
        typ = 'tar.gz'
    if not filename:
        return None
    print("Downloading:", filename, file=sys.stderr)
    try:
        with open(target, 'wb') as fid:
            # anaconda.org download location guards against
            # scraping so trick it with a fake browser header
            # see: https://medium.com/@speedforcerun/python-crawler-http-error-403-forbidden-1623ae9ba0f
            headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.3'}
            req = Request(url=filename, headers=headers)
            fid.write(urlopen(req).read())
        with open(target, 'rb') as binary_to_check:
            data = binary_to_check.read()
            sha256_returned = hashlib.sha256(data).hexdigest()
            sha256_expected = sha256_vals[os.path.basename(filename)]
            if sha256_returned != sha256_expected:
                raise ValueError('sha256 hash mismatch for downloaded OpenBLAS')

    except HTTPError as e:
        print(f'Could not download "{filename}"')
        print(f'Error message: {e}')
        return None
    return typ

def setup_openblas(arch=get_arch(), ilp64=get_ilp64()):
    '''
    Download and setup an openblas library for building. If successful,
    the configuration script will find it automatically.

    Returns
    -------
    msg : str
        path to extracted files on success, otherwise indicates what went wrong
        To determine success, do ``os.path.exists(msg)``
    '''
    _, tmp = mkstemp()
    if not arch:
        raise ValueError('unknown architecture')
    typ = download_openblas(tmp, arch, ilp64)
    if not typ:
        return ''
    if arch == 'windows':
        if not typ == 'zip':
            return 'expecting to download zipfile on windows, not %s' % str(typ)
        return unpack_windows_zip(tmp)
    else:
        if not typ == 'tar.gz':
            return 'expecting to download tar.gz, not %s' % str(typ)
        return unpack_targz(tmp)

def unpack_windows_zip(fname):
    with zipfile.ZipFile(fname, 'r') as zf:
        # Get the openblas.a file, but not openblas.dll.a nor openblas.dev.a
        lib = [x for x in zf.namelist() if OPENBLAS_LONG in x and
                  x.endswith('a') and not x.endswith('dll.a') and
                  not x.endswith('dev.a')]
        if not lib:
            return 'could not find libopenblas_%s*.a ' \
                    'in downloaded zipfile' % OPENBLAS_LONG
        target = os.path.join(gettempdir(), 'openblas.a')
        with open(target, 'wb') as fid:
            fid.write(zf.read(lib[0]))
    return target

def unpack_targz(fname):
    target = os.path.join(gettempdir(), 'openblas')
    if not os.path.exists(target):
        os.mkdir(target)
    with tarfile.open(fname, 'r') as zf:
        # Strip common prefix from paths when unpacking
        prefix = os.path.commonpath(zf.getnames())
        extract_tarfile_to(zf, target, prefix)
        return target

def extract_tarfile_to(tarfileobj, target_path, archive_path):
    """Extract TarFile contents under archive_path/ to target_path/"""

    target_path = os.path.abspath(target_path)

    def get_members():
        for member in tarfileobj.getmembers():
            if archive_path:
                norm_path = os.path.normpath(member.name)
                if norm_path.startswith(archive_path + os.path.sep):
                    member.name = norm_path[len(archive_path)+1:]
                else:
                    continue

            dst_path = os.path.abspath(os.path.join(target_path, member.name))
            if os.path.commonpath([target_path, dst_path]) != target_path:
                # Path not under target_path, probably contains ../
                continue

            yield member

    tarfileobj.extractall(target_path, members=get_members())

def make_init(dirname):
    '''
    Create a _distributor_init.py file for OpenBlas
    '''
    with open(os.path.join(dirname, '_distributor_init.py'), 'wt') as fid:
        fid.write(textwrap.dedent("""
            '''
            Helper to preload windows dlls to prevent dll not found errors.
            Once a DLL is preloaded, its namespace is made available to any
            subsequent DLL. This file originated in the numpy-wheels repo,
            and is created as part of the scripts that build the wheel.
            '''
            import os
            from ctypes import WinDLL
            import glob
            if os.name == 'nt':
                # convention for storing / loading the DLL from
                # numpy/.libs/, if present
                try:
                    basedir = os.path.dirname(__file__)
                except:
                    pass
                else:
                    libs_dir = os.path.abspath(os.path.join(basedir, '.libs'))
                    DLL_filenames = []
                    if os.path.isdir(libs_dir):
                        for filename in glob.glob(os.path.join(libs_dir,
                                                             '*openblas*dll')):
                            # NOTE: would it change behavior to load ALL
                            # DLLs at this path vs. the name restriction?
                            WinDLL(os.path.abspath(filename))
                            DLL_filenames.append(filename)
                if len(DLL_filenames) > 1:
                    import warnings
                    warnings.warn("loaded more than 1 DLL from .libs:\\n%s" %
                              "\\n".join(DLL_filenames),
                              stacklevel=1)
    """))

def test_setup(arches):
    '''
    Make sure all the downloadable files exist and can be opened
    '''
    def items():
        for arch in arches:
            yield arch, None
            if arch in ('x86', 'darwin', 'windows'):
                yield arch, '64_'

    for arch, ilp64 in items():
        if arch == '':
            continue

        target = None
        try:
            try:
                target = setup_openblas(arch, ilp64)
            except Exception:
                print(f'Could not setup {arch}')
                raise
            if not target:
                raise RuntimeError(f'Could not setup {arch}')
            print(target)
            if arch == 'windows':
                if not target.endswith('.a'):
                    raise RuntimeError("Not .a extracted!")
            else:
                files = glob.glob(os.path.join(target, "lib", "*.a"))
                if not files:
                    raise RuntimeError("No lib/*.a unpacked!")
        finally:
            if target is not None:
                if os.path.isfile(target):
                    os.unlink(target)
                else:
                    shutil.rmtree(target)

def test_version(expected_version, ilp64=get_ilp64()):
    """
    Assert that expected OpenBLAS version is
    actually available via SciPy
    """
    import scipy
    import scipy.linalg
    import ctypes

    dll = ctypes.CDLL(scipy.linalg.cython_blas.__file__)
    if ilp64 == "64_":
        get_config = dll.openblas_get_config64_
    else:
        get_config = dll.openblas_get_config
    get_config.restype = ctypes.c_char_p
    res = get_config()
    print('OpenBLAS get_config returned', str(res))
    check_str = b'OpenBLAS %s' % expected_version[0].encode()
    assert check_str in res

    if 'dev' not in expected_version[0]:
        assert b'dev' not in res

    if ilp64:
        assert b"USE64BITINT" in res
    else:
        assert b"USE64BITINT" not in res


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Download and expand an OpenBLAS archive for this '
                    'architecture')
    parser.add_argument('--test', nargs='*', default=None,
        help='Test different architectures. "all", or any of %s' % ARCHITECTURES)
    parser.add_argument('--check_version', nargs=1, default=None,
        help='Check provided OpenBLAS version string against available OpenBLAS')
    args = parser.parse_args()
    if args.check_version is not None:
        test_version(args.check_version)
    elif args.test is None:
        print(setup_openblas())
    else:
        if len(args.test) == 0 or 'all' in args.test:
            test_setup(ARCHITECTURES)
        else:
            test_setup(args.test)
