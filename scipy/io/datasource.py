"""A generic interface for importing data.  Data sources can originate from
URLs or local paths and can be in compressed or uncompressed form.

"""

# TODO: Make DataSource and Repository the public interface.
#       _Cache will be used internally.  Add methods in DataSource to expose
#       some of the functionality that exists only in _Cache currently.

__docformat__ = "restructuredtext en"

import os
import gzip
import bz2
from urlparse import urlparse
from urllib2 import urlopen, URLError
from tempfile import mkstemp

import warnings

# datasource has been used for a while in the NIPY project for analyzing
# large fmri imaging files hosted over a network.  Data would be fetched
# via URLs, cached locally and analyzed. Under these conditions the code
# worked well, however it needs to be documented, tested and reviewed
# before being fully exposed to SciPy.  We hope to do this before the
# 0.7 release.
_api_warning = "The datasource API will be changing frequently before \
the 0.7 release as the code is ported from the NIPY project to SciPy. \
Some of the current public interface may become private during the port! \
Use this module minimally, if at all, until it is stabilized."

warnings.warn(_api_warning)

# TODO: Don't check file extensions, just try and open zip files?
#       Follow the, it's easier to get forgiveness than permission?
#       And because this shouldn't succeed:
#           In [134]: datasource._iszip('/my/fake/dir/foobar.gz')
#           Out[134]: True

# TODO: .zip support, .tar support?
_file_openers = {".gz":gzip.open, ".bz2":bz2.BZ2File, None:file}

def _iszip(filename):
    """Test if the given file is a zip file.

    Currently only looks at the file extension, not very robust.
    """

    fname, ext = os.path.splitext(filename)
    return ext in _file_openers.keys()

def _unzip(filename):
    """Unzip the given file and return the path object to the new file."""

    # This function is not used in datasource.  Appears it was created
    # so users could unzip a file without having to import the corresponding
    # compression module.  Should this be part of datasource?
    if not _iszip(filename):
        raise ValueError("file %s is not zipped"%filename)
    unzip_name, zipext = _splitzipext(filename)
    opener = _file_openers[zipext]
    outfile = file(unzip_name, 'w')
    outfile.write(opener(filename).read())
    outfile.close()
    return unzip_name

def _iswritemode(mode):
    """Test if the given mode will open a file for writing."""

    _writemodes = ("w", "+", "a")
    for c in mode:
        if c in _writemodes: return True
    return False

def _splitzipext(filename):
    """Split the filename into a path object and a zip extension.

    If the filename does not have a zip extension the zip_ext in the
    return will be None.

    Parameters:

        filename : {string}
            Filename to split.

    Returns:

        base, zip_ext : {tuple}
            Tuple containing a path object to the file and the zip extension.

    """

    if _iszip(filename):
        return os.path.splitext(filename)
    else:
        return filename, None

def _isurl(pathstr):
    """Test whether a given string can be parsed as a URL.

    A pathstr with a valid network scheme (http, ftp, ...) will return true.
    A pathstr with a 'file' scheme will return false.

    """

    scheme, netloc, upath, uparams, uquery, ufrag = urlparse(pathstr)
    return bool(scheme and netloc)

def _ensuredirs(directory):
    """Ensure that the given directory path exists.  If not, create it."""

    if not os.path.exists(directory):
        os.makedirs(directory)

class _Cache (object):
    """A local file cache for URL datasources.

    The path of the cache can be specified on intialization.  The default
    path is ~/.scipy/cache


    """

    def __init__(self, cachepath=None):
        if cachepath is not None:
            self.path = cachepath
        elif os.name == 'posix':
            self.path = os.path.join(os.environ["HOME"], ".scipy", "cache")
        elif os.name == 'nt':
            self.path = os.path.join(os.environ["HOMEPATH"], ".scipy", "cache")
        if not os.path.exists(self.path):
            _ensuredirs(self.path)

    def tempfile(self, suffix='', prefix=''):
        """Create and return a temporary file in the cache.

        Parameters:
            suffix : {''}, optional

            prefix : {''}, optional

        Returns:
            tmpfile : {string}
                String containing the full path to the temporary file.

        Examples

            >>> mycache = datasource._Cache()
            >>> mytmpfile = mycache.tempfile()
            >>> mytmpfile
            '/home/guido/.scipy/cache/GUPhDv'

        """

        _tmp, fname = mkstemp(suffix, prefix, self.path)
        return fname

    def filepath(self, uri):
        """Return a path object to the uri in the cache.

        Parameters:
            uri : {string}
                Filename to use in the returned path object.

        Returns:
            path : {string}
                Complete path for the given uri.

        Examples

            >>> mycache = datasource._Cache()
            >>> mycache.filepath('xyzcoords.txt')
            '/home/guido/.scipy/cache/xyzcoords.txt'

        """

        scheme, netloc, upath, uparams, uquery, ufrag = urlparse(uri)
        return os.path.join(self.path, netloc, upath.strip('/'))

    def cache(self, uri):
        """Copy the file at uri into the cache.

        Parameters:
            uri : {string}
                path or url of source file to cache.

        Returns:
            None

        """

        if self.iscached(uri):
            return

        upath = self.filepath(uri)
        _ensuredirs(os.path.dirname(upath))

        print 'cache - source:', uri
        print 'cache - destination:', upath

        if _isurl(uri):
            try:
                openedurl = urlopen(uri)
                file(upath, 'w').write(openedurl.read())
            except URLError:
                raise URLError("URL not found: " + str(uri))
        else:
            try:
                fp = file(uri, 'r')
                file(upath, 'w').write(fp.read())
            except IOError:
                raise IOError("File not founcd: " + str(uri))

    def clear(self):
        """Delete all files in the cache."""

        # TODO: This deletes all files in the cache directory, regardless
        #       of if this instance created them.  Too destructive and
        #       unexpected behavior.
        #for _file in self.path.files():
        #    os.remove(file)
        raise NotImplementedError

    def iscached(self, uri):
        """ Check if a file exists in the cache.

        Returns
            boolean

        """

        upath = self.filepath(uri)
        return os.path.exists(upath)

    def retrieve(self, uri):
        """Retrieve a file from the cache.
        If not already there, create the file and add it to the cache.

        Returns
            open file object

        """

        self.cache(uri)
        return file(self.filepath(uri))


class DataSource (object):
    """A generic data source class.

    Data sets could be from a file, a URL, or a cached file.  They may also
    be compressed or uncompressed.

    TODO: Improve DataSource docstring

    """

    def __init__(self, cachepath=os.curdir):
        self._cache = _Cache(cachepath)

    def tempfile(self, suffix='', prefix=''):
        """Create a temporary file in the DataSource cache.

        Parameters:
            suffix : {''}, optional

            prefix : {''}, optional

        Returns:
            tmpfile : {string}
                String containing the full path to the temporary file.

        Examples

            >>> datasrc = datasource.DataSource()
            >>> tmpfile = datasrc.tempfile()
            >>> tmpfile
            '/home/guido/src/scipy-trunk/scipy/io/PZTuKo'

        """
        return self._cache.tempfile(suffix, prefix)

    def _possible_names(self, filename):
        """Return a tuple containing compressed filenames."""
        names = [filename]
        if not _iszip(filename):
            for zipext in _file_openers.keys():
                names.append(filename+zipext)
        return tuple(names)

    def cache(self, pathstr):
        """Cache the file specified by pathstr.

        Creates a copy of file pathstr in the datasource cache.

        """

        self._cache.cache(pathstr)

    def clear(self):
        # TODO: Implement a clear interface for deleting tempfiles.
        #       There's a problem with the way this is handled in the _Cache,
        #       All files in the cache directory will be deleted.  In the
        #       default instance, that's the os.curdir.  I doubt this is what
        #       people would want.  The instance should only delete files that
        #       it created!
        raise NotImplementedError

    def filename(self, pathstr):
        """Searches for pathstr file and returns full path if found.

        If pathstr is an URL, filename will cache a local copy and return
        the path to the cached file.
        If pathstr is a local file, filename will return a path to that local
        file.
        BUG:  This should be modified so the behavior is identical for both
              types of files!

        The search will include possible compressed versions of the files.
        BUG:  Will return the first occurence found, regardless of which
              version is newer.

        """

        found = None
        for name in self._possible_names(pathstr):
            try:
                if _isurl(name):
                    self.cache(name)
                    found = self._cache.filepath(name)
                else:
                    raise Exception
            except:
                if os.path.exists(name):
                    found = name
            if found:
                break
        if found is None:
            raise IOError("%s not found"%pathstr)
        return found

    def exists(self, pathstr):
        """Test if pathstr exists in the cache or the current directory.

        If pathstr is an URL, it will be fetched and cached.

        """

        # Is this method doing to much?  At very least may want an option to
        # not fetch and cache URLs.

        try:
            _datafile = self.filename(pathstr)
            return True
        except IOError:
            return False

    def open(self, pathstr, mode='r'):
        """Open pathstr and return file object.

        If pathstr is an URL, it will be fetched and cached.

        """

        # Is this method doing to much?  Should it be fetching and caching?

        if _isurl(pathstr) and _iswritemode(mode):
            raise ValueError("URLs are not writeable")
        found = self.filename(pathstr)
        _fname, ext = _splitzipext(found)
        if ext == 'bz2':
            mode.replace("+", "")
        return _file_openers[ext](found, mode=mode)


class Repository (DataSource):
    """Multiple DataSource's that share one base url.

    TODO: Improve Repository docstring.

    """

    def __init__(self, baseurl, cachepath=None):
        DataSource.__init__(self, cachepath=cachepath)
        self._baseurl = baseurl

    def _fullpath(self, pathstr):
        return os.path.join(self._baseurl, pathstr)

    def filename(self, pathstr):
        return DataSource.filename(self, self._fullpath(pathstr))

    def exists(self, pathstr):
        return DataSource.exists(self, self._fullpath(pathstr))

    def open(self, pathstr, mode='r'):
        return DataSource.open(self, self._fullpath(pathstr), mode)
