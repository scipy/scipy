"""Utilities for importing (possibly compressed) data sets from an URL
(or file) and possibly caching them.

"""

__docformat__ = "restructuredtext en"

import os
import gzip
import bz2
from urlparse import urlparse
from urllib2 import urlopen
from tempfile import mkstemp

# TODO: replace with newer tuple-based path module
from scipy.io.path import path

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

# TODO: .zip support
zipexts = (".gz",".bz2")
file_openers = {".gz":gzip.open, ".bz2":bz2.BZ2File, None:file}

def iszip(filename):
    """Test if the given file is a zip file.

    *Parameters*:

        filename : {string}
            Filename to test.

    *Returns*:

        bool
            Results of test.

    """

    _tmp, ext = path(filename).splitext()
    return ext in zipexts

def unzip(filename):
    """Unzip the given file and return the path object to the new file.

    *Parameters*:

        filename : string
            Filename to unzip.

    *Returns*:

        path_obj
            Path object of the unzipped file.

    """

    if not iszip(filename):
        raise ValueError("file %s is not zipped"%filename)
    unzip_name, zipext = splitzipext(filename)
    opener = file_openers[zipext]
    outfile = file(unzip_name, 'w')
    outfile.write(opener(filename).read())
    outfile.close()
    return unzip_name

def iswritemode(mode):
    """Test if the given mode will open a file for writing.

    *Parameters*:

        mode : {string}
            The mode to be checked.

    *Returns*:

        bool
            Result of test.

    """

    _writemodes = ("w", "+", "a")
    for c in mode:
        if c in _writemodes: return True
    return False

def splitzipext(filename):
    """Split the filename into a path object and a zip extension.

    If the filename does not have a zip extension the zip_ext in the
    return will be None.
    
    *Parameters*:

        filename : {string}
            Filename to split.

    *Returns*:

        base, zip_ext : {tuple}
            Tuple containing a path object to the file and the zip extension.

    """

    if iszip(filename):
        return path(filename).splitext()
    else:
        return filename, None

def isurl(pathstr):
    """Test whether a given string can be parsed as a URL.

    *Parameters*:
    
        pathstr : {string}
            The string to be checked.

    *Returns*:
        bool
            Results of test.

    """

    scheme, netloc, _tmp, _tmp, _tmp, _tmp = urlparse(pathstr)
    return bool(scheme and netloc)

def ensuredirs(directory):
    """Ensure that the given directory path exists.  If not, create it.

    *Parameters*:
        directory : {path object}
        
    *Returns*:
        None

    """

    if not isinstance(directory, path):
        directory = path(directory)
    if not directory.exists():
        directory.makedirs()

class Cache (object):
    """A local file cache for URL datasources.

    The path of the cache can be specified on intialization.  The default
    path is ~/.scipy/cache


    """

    def __init__(self, cachepath=None):
        if cachepath is not None:
            self.path = path(cachepath)
        elif os.name == 'posix':
            self.path = path(os.environ["HOME"]).joinpath(".scipy","cache")
        elif os.name == 'nt':
            self.path = path(os.environ["HOMEPATH"]).joinpath(".scipy","cache")
        if not self.path.exists():
            ensuredirs(self.path)

    def tempfile(self, suffix='', prefix=''):
        """Create and return a temporary file in the cache.

        *Parameters*:
            suffix : {''}, optional

            prefix : {''}, optional

        *Returns*:
            tmpfile : {string}
                String containing the full path to the temporary file.

        *Examples*

            >>> mycache = datasource.Cache()
            >>> mytmpfile = mycache.tempfile()
            >>> mytmpfile
            '/home/guido/.scipy/cache/GUPhDv'

        """

        _tmp, fname = mkstemp(suffix, prefix, self.path)
        return fname

    def filepath(self, uri):
        """Return a path object to the uri in the cache.

        *Parameters*:
            uri : {string}
                Filename to use in the returned path object.

        *Returns*:
            path_obj
                Path object for the given uri.
                
        *Examples*

            >>> mycache = datasource.Cache()
            >>> mycache.filepath('xyzcoords.txt')
            path('/home/guido/.scipy/cache/xyzcoords.txt')

        """
        # TODO: Change to non-public?

        #       It appears the Cache is designed to work with URLs only!
        (_tmp, netloc, upath, _tmp, _tmp, _tmp) = urlparse(uri)
        return self.path.joinpath(netloc, upath.strip('/'))

    def filename(self, uri):
        """Return the complete path + filename within the cache.

        *Parameters*:
            uri : {string}
                Filename to usein the returned path.
                
        *Returns*:
            filename

        *Examples*

            >>> mycache = datasource.Cache()
            >>> mycache.filename('xyzcoords.txt')
            '/home/guido/.scipy/cache/yzcoords.txt'

        """
        # TODO: Change to non-public?
                
        return str(self.filepath(uri))

    def cache(self, uri):
        """Copy a file into the cache.

        :Returns: ``None``
        """
        if self.iscached(uri):
            return

        print 'cache uri:', uri
        
        upath = self.filepath(uri)

        print 'upath:', upath
        
        ensuredirs(upath.dirname())
        try:
            print 'uri:', uri
            openedurl = urlopen(uri)
        except:
            raise IOError("url not found: "+str(uri))
        file(upath, 'w').write(openedurl.read())

    def clear(self):
        """ Delete all files in the cache.

        :Returns: ``None``
        """
        for _file in self.path.files():
            _file.rm()

    def iscached(self, uri):
        """ Check if a file exists in the cache.

        :Returns: ``bool``
        """
        return self.filepath(uri).exists()

    def retrieve(self, uri):
        """
        Retrieve a file from the cache.
        If not already there, create the file and
        add it to the cache.

        :Returns: ``file``
        """
        self.cache(uri)
        return file(self.filename(uri))


class DataSource (object):
    """A generic data source class.

    Data sets could be from a file, a URL, or a cached file.  They may also
    be compressed or uncompressed.

    TODO: Improve DataSource docstring

    """

    def __init__(self, cachepath=os.curdir):
        self._cache = Cache(cachepath)

    def tempfile(self, suffix='', prefix=''):
        """Return an temporary file name in the cache."""
        return self._cache.tempfile(suffix, prefix)

    def _possible_names(self, filename):
        names = [filename]
        if not iszip(filename):
            for zipext in zipexts:
                names.append(filename+zipext)
        return tuple(names)

    def cache(self, pathstr):
        if isurl(pathstr):
            self._cache.cache(pathstr)

    def filename(self, pathstr):
        found = None
        for name in self._possible_names(pathstr):
            try:
                if isurl(name):
                    self.cache(name)
                    found = self._cache.filename(name)
                else:
                    raise Exception
            except:
                if path(name).exists():
                    found = name
            if found:
                break
        if found is None:
            raise IOError("%s not found"%pathstr)
        return found

    def exists(self, pathstr):
        try:
            _ = self.filename(pathstr)
            return True
        except IOError:
            return False

    def open(self, pathstr, mode='r'):
        if isurl(pathstr) and iswritemode(mode):
            raise ValueError("URLs are not writeable")
        found = self.filename(pathstr)
        _, ext = splitzipext(found)
        if ext == 'bz2':
            mode.replace("+", "")
        return file_openers[ext](found, mode=mode)

    def _fullpath(self, pathstr):
        return pathstr


class Repository (DataSource):
    """Multiple DataSource's that share one base url.

    TODO: Improve Repository docstring.

    """

    #"""DataSource with an implied root."""

    def __init__(self, baseurl, cachepath=None):
        DataSource.__init__(self, cachepath=cachepath)
        self._baseurl = baseurl

    def _fullpath(self, pathstr):
        return path(self._baseurl).joinpath(pathstr)

    def filename(self, pathstr):
        return DataSource.filename(self, str(self._fullpath(pathstr)))

    def exists(self, pathstr):
        return DataSource.exists(self, self._fullpath(pathstr))

    def open(self, pathstr, mode='r'):
        return DataSource.open(self, self._fullpath(pathstr), mode)
