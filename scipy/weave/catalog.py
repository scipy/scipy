""" Track relationships between compiled extension functions & code fragments

    catalog keeps track of which compiled(or even standard) functions are
    related to which code fragments.  It also stores these relationships
    to disk so they are remembered between Python sessions.  When

        a = 1
        compiler.inline('printf("printed from C: %d",a);',['a'] )

    is called, inline() first looks to see if it has seen the code
    'printf("printed from C");' before.  If not, it calls

        catalog.get_functions('printf("printed from C: %d", a);')

    which returns a list of all the function objects that have been compiled
    for the code fragment.  Multiple functions can occur because the code
    could be compiled for different types for 'a' (although not likely in
    this case). The catalog first looks in its cache and quickly returns
    a list of the functions if possible.  If the cache lookup fails, it then
    looks through possibly multiple catalog files on disk and fills its
    cache with all the functions that match the code fragment.

    In case where the code fragment hasn't been compiled, inline() compiles
    the code and then adds it to the catalog:

        function = <code to compile function>
        catalog.add_function('printf("printed from C: %d", a);',function)

    add_function() adds function to the front of the cache.  function,
    along with the path information to its module, are also stored in a
    persistent catalog for future use by python sessions.
"""
from __future__ import absolute_import, print_function

import os
import sys
import stat
import pickle
import socket
import tempfile

try:
    import dbhash
    import shelve
    dumb = 0
except ImportError:
    from . import _dumb_shelve as shelve
    dumb = 1

# For testing...
# import scipy.io.dumb_shelve as shelve
# dumb = 1

# import shelve
# dumb = 0


def getmodule(object):
    """ Discover the name of the module where object was defined.

        This is an augmented version of inspect.getmodule that can discover
        the parent module for extension functions.
    """
    import inspect
    value = inspect.getmodule(object)
    if value is None:
        # walk trough all modules looking for function
        for name,mod in sys.modules.items():
            # try except used because of some comparison failures
            # in wxPoint code.  Need to review this
            try:
                if mod and object in mod.__dict__.values():
                    value = mod
                    # if it is a built-in module, keep looking to see
                    # if a non-builtin also has it.  Otherwise quit and
                    # consider the module found. (ain't perfect, but will
                    # have to do for now).
                    if str(mod) not in '(built-in)':
                        break

            except (TypeError, KeyError, ImportError):
                pass
    return value


def expr_to_filename(expr):
    """ Convert an arbitrary expr string to a valid file name.

        The name is based on the SHA-256 check sum for the string and
        Something that was a little more human readable would be
        nice, but the computer doesn't seem to care.
    """
    from hashlib import sha256
    base = 'sc_'
    # 32 chars is enough for unique filenames; too long names don't work for
    # MSVC (see gh-3216).  Don't use md5, gives a FIPS warning.
    return base + sha256(expr).hexdigest()[:32]


def unique_file(d,expr):
    """ Generate a unqiue file name based on expr in directory d

        This is meant for use with building extension modules, so
        a file name is considered unique if none of the following
        extension '.cpp','.o','.so','module.so','.py', or '.pyd'
        exists in directory d.  The fully qualified path to the
        new name is returned.  You'll need to append your own
        extension to it before creating files.
    """
    files = os.listdir(d)
    # base = 'scipy_compile'
    base = expr_to_filename(expr)
    for i in xrange(1000000):
        fname = base + repr(i)
        if not (fname+'.cpp' in files or
                fname+'.o' in files or
                fname+'.so' in files or
                fname+'module.so' in files or
                fname+'.py' in files or
                fname+'.pyd' in files):
            break
    return os.path.join(d,fname)


def is_writable(dir):
    """Determine whether a given directory is writable in a portable manner.

    Parameters
    ----------
    dir : str
        A string represeting a path to a directory on the filesystem.

    Returns
    -------
    res : bool
        True or False.
    """
    if not os.path.isdir(dir):
        return False

    # Do NOT use a hardcoded name here due to the danger from race conditions
    # on NFS when multiple processes are accessing the same base directory in
    # parallel.  We use both hostname and process id for the prefix in an
    # attempt to ensure that there can really be no name collisions (tempfile
    # appends 6 random chars to this prefix).
    prefix = 'dummy_%s_%s_' % (socket.gethostname(),os.getpid())
    try:
        tmp = tempfile.TemporaryFile(prefix=prefix,dir=dir)
    except OSError:
        return False
    # The underlying file is destroyed upon closing the file object (under
    # *nix, it was unlinked at creation time)
    tmp.close()
    return True


def whoami():
    """return a string identifying the user."""
    return os.environ.get("USER") or os.environ.get("USERNAME") or "unknown"


def _create_dirs(path):
    """ create provided path, ignore errors """
    try:
        os.makedirs(path, mode=0o700)
    except OSError:
        pass


def default_dir_posix(tmp_dir=None):
    """
    Create or find default catalog store for posix systems

    purpose of 'tmp_dir' is to enable way how to test this function easily
    """
    path_candidates = []
    python_name = "python%d%d_compiled" % tuple(sys.version_info[:2])

    if tmp_dir:
        home_dir = tmp_dir
    else:
        home_dir = os.path.expanduser('~')
    tmp_dir = tmp_dir or tempfile.gettempdir()

    xdg_cache = (os.environ.get("XDG_CACHE_HOME", None) or
                 os.path.join(home_dir, '.cache'))
    xdg_temp_dir = os.path.join(xdg_cache, 'scipy', python_name)
    path_candidates.append(xdg_temp_dir)

    home_temp_dir_name = '.' + python_name
    home_temp_dir = os.path.join(home_dir, home_temp_dir_name)
    path_candidates.append(home_temp_dir)

    temp_dir_name = repr(os.getuid()) + '_' + python_name
    temp_dir_path = os.path.join(tmp_dir, temp_dir_name)
    path_candidates.append(temp_dir_path)

    for path in path_candidates:
        _create_dirs(path)
        if check_dir(path):
            return path

    # since we got here, both dirs are not useful
    tmp_dir_path = find_valid_temp_dir(temp_dir_name, tmp_dir)
    if not tmp_dir_path:
        tmp_dir_path = create_temp_dir(temp_dir_name, tmp_dir=tmp_dir)
    return tmp_dir_path


def default_dir_win(tmp_dir=None):
    """
    Create or find default catalog store for Windows systems

    purpose of 'tmp_dir' is to enable way how to test this function easily
    """
    def create_win_temp_dir(prefix, inner_dir=None, tmp_dir=None):
        """
        create temp dir starting with 'prefix' in 'tmp_dir' or
        'tempfile.gettempdir'; if 'inner_dir' is specified, it should be
        created inside
        """
        tmp_dir_path = find_valid_temp_dir(prefix, tmp_dir)
        if tmp_dir_path:
            if inner_dir:
                tmp_dir_path = os.path.join(tmp_dir_path, inner_dir)
                if not os.path.isdir(tmp_dir_path):
                    os.mkdir(tmp_dir_path, 0o700)
        else:
            tmp_dir_path = create_temp_dir(prefix, inner_dir, tmp_dir)
        return tmp_dir_path

    python_name = "python%d%d_compiled" % tuple(sys.version_info[:2])
    tmp_dir = tmp_dir or tempfile.gettempdir()

    temp_dir_name = "%s" % whoami()
    temp_root_dir = os.path.join(tmp_dir, temp_dir_name)
    temp_dir_path = os.path.join(temp_root_dir, python_name)
    _create_dirs(temp_dir_path)
    if check_dir(temp_dir_path) and check_dir(temp_root_dir):
        return temp_dir_path
    else:
        if check_dir(temp_root_dir):
            return create_win_temp_dir(python_name, tmp_dir=temp_root_dir)
        else:
            return create_win_temp_dir(temp_dir_name, python_name, tmp_dir)


def default_dir():
    """ Return a default location to store compiled files and catalogs.

        XX is the Python version number in all paths listed below
        On windows, the default location is the temporary directory
        returned by gettempdir()/pythonXX.

        On Unix, ~/.pythonXX_compiled is the default location.  If it doesn't
        exist, it is created.  The directory is marked rwx------.

        If for some reason it isn't possible to build a default directory
        in the user's home, /tmp/<uid>_pythonXX_compiled is used.  If it
        doesn't exist, it is created.  The directory is marked rwx------
        to try and keep people from being able to sneak a bad module
        in on you. If the directory already exists in /tmp/ and is not
        secure, new one is created.
    """
    # Use a cached value for fast return if possible
    if hasattr(default_dir, "cached_path") and \
       check_dir(default_dir.cached_path):
        return default_dir.cached_path

    if sys.platform == 'win32':
        path = default_dir_win()
    else:
        path = default_dir_posix()

    # Cache the default dir path so that this function returns quickly after
    # being called once (nothing in it should change after the first call)
    default_dir.cached_path = path

    return path


def check_dir(im_dir):
    """
    Check if dir is safe; if it is, return True.
    These checks make sense only on posix:
     * directory has correct owner
     * directory has correct permissions (0700)
     * directory is not a symlink
    """
    def check_is_dir():
        return os.path.isdir(im_dir)

    def check_permissions():
        """ If on posix, permissions should be 0700. """
        writable = is_writable(im_dir)
        if sys.platform != 'win32':
            try:
                im_dir_stat = os.stat(im_dir)
            except OSError:
                return False
            writable &= stat.S_IMODE(im_dir_stat.st_mode) == 0o0700
        return writable

    def check_ownership():
        """ Intermediate dir owner should be same as owner of process. """
        if sys.platform != 'win32':
            try:
                im_dir_stat = os.stat(im_dir)
            except OSError:
                return False
            proc_uid = os.getuid()
            return proc_uid == im_dir_stat.st_uid
        return True

    def check_is_symlink():
        """ Check if intermediate dir is symlink. """
        try:
            return not os.path.islink(im_dir)
        except OSError:
            return False

    checks = [check_is_dir, check_permissions,
              check_ownership, check_is_symlink]

    for check in checks:
        if not check():
            return False

    return True


def create_temp_dir(prefix, inner_dir=None, tmp_dir=None):
    """
    Create intermediate dirs <tmp>/<prefix+random suffix>/<inner_dir>/

    argument 'tmp_dir' is used in unit tests
    """
    if not tmp_dir:
        tmp_dir_path = tempfile.mkdtemp(prefix=prefix)
    else:
        tmp_dir_path = tempfile.mkdtemp(prefix=prefix, dir=tmp_dir)
    if inner_dir:
        tmp_dir_path = os.path.join(tmp_dir_path, inner_dir)
        os.mkdir(tmp_dir_path, 0o700)
    return tmp_dir_path


def intermediate_dir_prefix():
    """ Prefix of root intermediate dir (<tmp>/<root_im_dir>). """
    return "%s-%s-" % ("scipy", whoami())


def find_temp_dir(prefix, tmp_dir=None):
    """ Find temp dirs in 'tmp_dir' starting with 'prefix'"""
    matches = []
    tmp_dir = tmp_dir or tempfile.gettempdir()
    for tmp_file in os.listdir(tmp_dir):
        if tmp_file.startswith(prefix):
            matches.append(os.path.join(tmp_dir, tmp_file))
    return matches


def find_valid_temp_dir(prefix, tmp_dir=None):
    """
    Try to look for existing temp dirs.
    If there is one suitable found, return it, otherwise return None.
    """
    matches = find_temp_dir(prefix, tmp_dir)
    for match in matches:
        if check_dir(match):
            # as soon as we find correct dir, we can stop searching
            return match


def py_intermediate_dir():
    """
    Name of intermediate dir for current python interpreter:
    <temp dir>/<name>/pythonXY_intermediate/
    """
    name = "python%d%d_intermediate" % tuple(sys.version_info[:2])
    return name


def create_intermediate_dir(tmp_dir=None):
    py_im_dir = py_intermediate_dir()
    return create_temp_dir(intermediate_dir_prefix(), py_im_dir, tmp_dir)


def intermediate_dir(tmp_dir=None):
    """
    Temporary directory for storing .cpp and .o files during builds.

    First, try to find the dir and if it exists, verify it is safe.
    Otherwise, create it.
    """
    im_dir = find_valid_temp_dir(intermediate_dir_prefix(), tmp_dir)
    py_im_dir = py_intermediate_dir()
    if im_dir is None:
        py_im_dir = py_intermediate_dir()
        im_dir = create_intermediate_dir(tmp_dir)
    else:
        im_dir = os.path.join(im_dir, py_im_dir)
        if not os.path.isdir(im_dir):
            os.mkdir(im_dir, 0o700)
    return im_dir


def default_temp_dir():
    path = os.path.join(default_dir(),'temp')
    if not os.path.exists(path):
        os.makedirs(path, mode=0o700)
    if not is_writable(path):
        print('warning: default directory is not write accessible.')
        print('default:', path)
    return path


def os_dependent_catalog_name():
    """ Generate catalog name dependent on OS and Python version being used.

        This allows multiple platforms to have catalog files in the
        same directory without stepping on each other.  For now, it
        bases the name of the value returned by sys.platform and the
        version of python being run.  If this isn't enough to descriminate
        on some platforms, we can try to add other info.  It has
        occurred to me that if we get fancy enough to optimize for different
        architectures, then chip type might be added to the catalog name also.
    """
    version = '%d%d' % sys.version_info[:2]
    return sys.platform+version+'compiled_catalog'


def catalog_path(module_path):
    """ Return the full path name for the catalog file in the given directory.

        module_path can either be a file name or a path name.  If it is a
        file name, the catalog file name in its parent directory is returned.
        If it is a directory, the catalog file in that directory is returned.

        If module_path doesn't exist, None is returned.  Note though, that the
        catalog file does *not* have to exist, only its parent.  '~', shell
        variables, and relative ('.' and '..') paths are all acceptable.

        catalog file names are os dependent (based on sys.platform), so this
        should support multiple platforms sharing the same disk space
        (NFS mounts). See os_dependent_catalog_name() for more info.
    """
    module_path = os.path.expanduser(module_path)
    module_path = os.path.expandvars(module_path)
    module_path = os.path.abspath(module_path)
    if not os.path.exists(module_path):
        catalog_file = None
    elif not os.path.isdir(module_path):
        module_path,dummy = os.path.split(module_path)
        catalog_file = os.path.join(module_path,os_dependent_catalog_name())
    else:
        catalog_file = os.path.join(module_path,os_dependent_catalog_name())
    return catalog_file


def get_catalog(module_path,mode='r'):
    """ Return a function catalog (shelve object) from the path module_path

        If module_path is a directory, the function catalog returned is
        from that directory.  If module_path is an actual module_name,
        then the function catalog returned is from its parent directory.
        mode uses the standard 'c' = create, 'n' = new, 'r' = read,
        'w' = write file open modes available for anydbm databases.

        Well... it should be.  Stuck with dumbdbm for now and the modes
        almost don't matter.  We do some checking for 'r' mode, but that
        is about it.

        See catalog_path() for more information on module_path.
    """
    if mode not in ['c','r','w','n']:
        msg = " mode must be 'c', 'n', 'r', or 'w'.  See anydbm for more info"
        raise ValueError(msg)
    catalog_file = catalog_path(module_path)
    if (catalog_file is not None) \
           and ((dumb and os.path.exists(catalog_file+'.dat'))
                or os.path.exists(catalog_file)):
        sh = shelve.open(catalog_file,mode)
    else:
        if mode == 'r':
            sh = None
        else:
            sh = shelve.open(catalog_file,mode)
    return sh


class catalog(object):
    """ Stores information about compiled functions both in cache and on disk.

        catalog stores (code, list_of_function) pairs so that all the functions
        that have been compiled for code are available for calling (usually in
        inline or blitz).

        catalog keeps a dictionary of previously accessed code values cached
        for quick access.  It also handles the looking up of functions compiled
        in previously called Python sessions on disk in function catalogs.
        catalog searches the directories in the PYTHONCOMPILED environment
        variable in order loading functions that correspond to the given code
        fragment.  A default directory is also searched for catalog functions.
        On unix, the default directory is usually '~/.pythonxx_compiled' where
        xx is the version of Python used. On windows, it is the directory
        returned by temfile.gettempdir().  Functions closer to the front are of
        the variable list are guaranteed to be closer to the front of the
        function list so that they will be called first.  See
        get_cataloged_functions() for more info on how the search order is
        traversed.

        Catalog also handles storing information about compiled functions to
        a catalog.  When writing this information, the first writable catalog
        file in PYTHONCOMPILED path is used.  If a writable catalog is not
        found, it is written to the catalog in the default directory.  This
        directory should always be writable.
    """
    def __init__(self,user_path_list=None):
        """ Create a catalog for storing/searching for compiled functions.

            user_path_list contains directories that should be searched
            first for function catalogs.  They will come before the path
            entries in the PYTHONCOMPILED environment varilable.
        """
        if isinstance(user_path_list, str):
            self.user_path_list = [user_path_list]
        elif user_path_list:
            self.user_path_list = user_path_list
        else:
            self.user_path_list = []
        self.cache = {}
        self.module_dir = None
        self.paths_added = 0
        # unconditionally append the default dir for auto-generated compiled
        # extension modules, so that pickle.load()s don't fail.
        sys.path.append(default_dir())

    def set_module_directory(self,module_dir):
        """ Set the path that will replace 'MODULE' in catalog searches.

            You should call clear_module_directory() when your finished
            working with it.
        """
        self.module_dir = module_dir

    def get_module_directory(self):
        """ Return the path used to replace the 'MODULE' in searches.
        """
        return self.module_dir

    def clear_module_directory(self):
        """ Reset 'MODULE' path to None so that it is ignored in searches.
        """
        self.module_dir = None

    def get_environ_path(self):
        """ Return list of paths from 'PYTHONCOMPILED' environment variable.

            On Unix the path in PYTHONCOMPILED is a ':' separated list of
            directories.  On Windows, a ';' separated list is used.
        """
        paths = []
        if 'PYTHONCOMPILED' in os.environ:
            path_string = os.environ['PYTHONCOMPILED']
            paths = path_string.split(os.path.pathsep)
        return paths

    def build_search_order(self):
        """ Returns a list of paths that are searched for catalogs.

            Values specified in the catalog constructor are searched first,
            then values found in the PYTHONCOMPILED environment variable.
            The directory returned by default_dir() is always returned at
            the end of the list.

            There is a 'magic' path name called 'MODULE' that is replaced
            by the directory defined by set_module_directory().  If the
            module directory hasn't been set, 'MODULE' is ignored.
        """

        paths = self.user_path_list + self.get_environ_path()
        search_order = []
        for path in paths:
            if path == 'MODULE':
                if self.module_dir:
                    search_order.append(self.module_dir)
            else:
                search_order.append(path)
        search_order.append(default_dir())
        return search_order

    def get_catalog_files(self):
        """ Returns catalog file list in correct search order.

            Some of the catalog files may not currently exists.
            However, all will be valid locations for a catalog
            to be created (if you have write permission).
        """
        files = map(catalog_path,self.build_search_order())
        files = filter(lambda x: x is not None,files)
        return files

    def get_existing_files(self):
        """ Returns all existing catalog file list in correct search order.
        """
        files = self.get_catalog_files()
        # open every stinking file to check if it exists.
        # This is because anydbm doesn't provide a consistent naming
        # convention across platforms for its files
        existing_files = []
        for file in files:
            cat = get_catalog(os.path.dirname(file),'r')
            if cat is not None:
                existing_files.append(file)
                cat.close()
        # This is the non-portable (and much faster) old code
        # existing_files = filter(os.path.exists,files)
        return existing_files

    def get_writable_file(self,existing_only=0):
        """ Return the name of the first writable catalog file.

            Its parent directory must also be writable.  This is so that
            compiled modules can be written to the same directory.
        """
        # note: both file and its parent directory must be writeable
        if existing_only:
            files = self.get_existing_files()
        else:
            files = self.get_catalog_files()
        # filter for (file exists and is writable) OR directory is writable

        def file_test(x):
            from os import access, F_OK, W_OK
            return (access(x,F_OK) and access(x,W_OK) or
                    access(os.path.dirname(x),W_OK))
        writable = filter(file_test,files)
        if writable:
            file = writable[0]
        else:
            file = None
        return file

    def get_writable_dir(self):
        """ Return the parent directory of first writable catalog file.

            The returned directory has write access.
        """
        return os.path.dirname(self.get_writable_file())

    def unique_module_name(self,code,module_dir=None):
        """ Return full path to unique file name that in writable location.

            The directory for the file is the first writable directory in
            the catalog search path.  The unique file name is derived from
            the code fragment.  If, module_dir is specified, it is used
            to replace 'MODULE' in the search path.
        """
        if module_dir is not None:
            self.set_module_directory(module_dir)
        try:
            d = self.get_writable_dir()
        finally:
            if module_dir is not None:
                self.clear_module_directory()
        return unique_file(d, code)

    def path_key(self,code):
        """ Return key for path information for functions associated with code.
        """
        return '__path__' + code

    def configure_path(self,cat,code):
        """ Add the python path for the given code to the sys.path

            unconfigure_path() should be called as soon as possible after
            imports associated with code are finished so that sys.path
            is restored to normal.
        """
        try:
            paths = cat[self.path_key(code)]
            self.paths_added = len(paths)
            sys.path = paths + sys.path
        except:
            self.paths_added = 0

    def unconfigure_path(self):
        """ Restores sys.path to normal after calls to configure_path()

            Remove the previously added paths from sys.path
        """
        sys.path = sys.path[self.paths_added:]
        self.paths_added = 0

    def get_cataloged_functions(self,code):
        """ Load all functions associated with code from catalog search path.

            Sometimes there can be trouble loading a function listed in a
            catalog file because the actual module that holds the function
            has been moved or deleted.  When this happens, that catalog file
            is "repaired", meaning the entire entry for this function is
            removed from the file.  This only affects the catalog file that
            has problems -- not the others in the search path.

            The "repair" behavior may not be needed, but I'll keep it for now.
        """
        mode = 'r'
        cat = None
        function_list = []
        for path in self.build_search_order():
            cat = get_catalog(path,mode)
            if cat is not None and code in cat:
                # set up the python path so that modules for this
                # function can be loaded.
                self.configure_path(cat,code)
                try:
                    function_list += cat[code]
                except:  # SystemError and ImportError so far seen
                    # problems loading a function from the catalog.  Try to
                    # repair the cause.
                    cat.close()
                    self.repair_catalog(path,code)
                self.unconfigure_path()
            if cat is not None:
                # ensure that the catalog is properly closed
                cat.close()
        return function_list

    def repair_catalog(self,catalog_path,code):
        """ Remove entry for code from catalog_path

            Occasionally catalog entries could get corrupted. An example
            would be when a module that had functions in the catalog was
            deleted or moved on the disk.  The best current repair method is
            just to trash the entire catalog entry for this piece of code.
            This may loose function entries that are valid, but thats life.

            catalog_path must be writable for repair.  If it isn't, the
            function exists with a warning.
        """
        writable_cat = None
        if (catalog_path is not None) and (not os.path.exists(catalog_path)):
            return
        try:
            writable_cat = get_catalog(catalog_path,'w')
        except:
            print('warning: unable to repair catalog entry\n %s\n in\n %s' %
                  (code,catalog_path))
            # shelve doesn't guarantee flushing, so it's safest to explicitly
            # close the catalog
            writable_cat.close()
            return
        if code in writable_cat:
            print('repairing catalog by removing key')
            del writable_cat[code]

        # it is possible that the path key doesn't exist (if the function
        # registered was a built-in function), so we have to check if the path
        # exists before arbitrarily deleting it.
        path_key = self.path_key(code)
        if path_key in writable_cat:
            del writable_cat[path_key]
        writable_cat.close()

    def get_functions_fast(self,code):
        """ Return list of functions for code from the cache.

            Return an empty list if the code entry is not found.
        """
        return self.cache.get(code,[])

    def get_functions(self,code,module_dir=None):
        """ Return the list of functions associated with this code fragment.

            The cache is first searched for the function.  If an entry
            in the cache is not found, then catalog files on disk are
            searched for the entry.  This is slooooow, but only happens
            once per code object.  All the functions found in catalog files
            on a cache miss are loaded into the cache to speed up future calls.
            The search order is as follows:

                1. user specified path (from catalog initialization)
                2. directories from the PYTHONCOMPILED environment variable
                3. The temporary directory on your platform.

            The path specified by module_dir will replace the 'MODULE'
            place holder in the catalog search path. See build_search_order()
            for more info on the search path.
        """
        # Fast!! try cache first.
        if code in self.cache:
            return self.cache[code]

        # 2. Slow!! read previously compiled functions from disk.
        try:
            self.set_module_directory(module_dir)
            function_list = self.get_cataloged_functions(code)
            # put function_list in cache to save future lookups.
            if function_list:
                self.cache[code] = function_list
            # return function_list, empty or otherwise.
        finally:
            self.clear_module_directory()
        return function_list

    def add_function(self,code,function,module_dir=None):
        """ Adds a function to the catalog.

            The function is added to the cache as well as the first
            writable file catalog found in the search path.  If no
            code entry exists in the cache, the on disk catalogs
            are loaded into the cache and function is added to the
            beginning of the function list.

            The path specified by module_dir will replace the 'MODULE'
            place holder in the catalog search path. See build_search_order()
            for more info on the search path.
        """

        # 1. put it in the cache.
        if code in self.cache:
            if function not in self.cache[code]:
                self.cache[code].insert(0,function)
            else:
                # if it is in the cache, then it is also
                # been persisted
                return
        else:
            # Load functions and put this one up front
            self.cache[code] = self.get_functions(code)
            self.fast_cache(code,function)
        # 2. Store the function entry to disk.
        try:
            self.set_module_directory(module_dir)
            self.add_function_persistent(code,function)
        finally:
            self.clear_module_directory()

    def add_function_persistent(self,code,function):
        """ Store the code->function relationship to disk.

            Two pieces of information are needed for loading functions
            from disk -- the function pickle (which conveniently stores
            the module name, etc.) and the path to its module's directory.
            The latter is needed so that the function can be loaded no
            matter what the user's Python path is.
        """
        # add function to data in first writable catalog
        mode = 'c'  # create if doesn't exist, otherwise, use existing
        cat_dir = self.get_writable_dir()
        cat = get_catalog(cat_dir,mode)
        if cat is None:
            cat_dir = default_dir()
            cat = get_catalog(cat_dir,mode)
        if cat is None:
            cat_dir = default_dir()
            cat_file = catalog_path(cat_dir)
            print('problems with default catalog -- removing')
            import glob
            files = glob.glob(cat_file+'*')
            for f in files:
                os.remove(f)
            cat = get_catalog(cat_dir,mode)
        if cat is None:
            raise ValueError('Failed to access a catalog for storing functions')
        # Prabhu was getting some corrupt catalog errors.  I'll put a try/except
        # to protect against this, but should really try and track down the issue.
        function_list = [function]
        try:
            function_list = function_list + cat.get(code,[])
        except pickle.UnpicklingError:
            pass
        cat[code] = function_list
        # now add needed path information for loading function
        module = getmodule(function)
        try:
            # built in modules don't have the __file__ extension, so this
            # will fail.  Just pass in this case since path additions aren't
            # needed for built-in modules.
            mod_path,f = os.path.split(os.path.abspath(module.__file__))
            pkey = self.path_key(code)
            cat[pkey] = [mod_path] + cat.get(pkey,[])
        except:
            pass
        cat.close()

    def fast_cache(self,code,function):
        """ Move function to the front of the cache entry for code

            If future calls to the function have the same type signature,
            this will speed up access significantly because the first
            function call is correct.

            Note:  The cache added to the inline_tools module is significantly
                   faster than always calling get_functions, so this isn't
                   as necessary as it used to be.  Still, it's probably worth
                   doing.
        """
        try:
            if self.cache[code][0] == function:
                return
        except:  # KeyError, IndexError
            pass
        try:
            self.cache[code].remove(function)
        except ValueError:
            pass
        # put new function at the beginning of the list to search.
        self.cache[code].insert(0,function)
