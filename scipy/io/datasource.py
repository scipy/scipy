"""A file interface for handling local and remote data files.
The goal of datasource is to abstract some of the file system operations when 
dealing with data files so the researcher doesn't have to know all the
low-level details.  Through datasource, a researcher can obtain and use a 
file with one function call, regardless of location of the file.

DataSource files can originate locally or remotely:

- local files : '/home/guido/src/local/data.txt'
- URLs (http, ftp, ...) : 'http://www.scipy.org/not/real/data.txt'

DataSource files can also be compressed or uncompressed.  Currently only gzip
and bz2 are supported.

Example:

    >>> # Create a DataSource and use '/home/guido/tmpdata/' for local storage.
    >>> ds = datasource.DataSource('/home/guido/tmpdata/')
    >>>
    >>> # Open a remote, gzipped file.
    >>> # DataSource downloads the file, stores it locally in:
    >>> #     '/home/guido/tmpdata/www.scipy.org/not/real/data.txt.gz'
    >>> # opens the file with the gzip module and returns a file-like object.
    >>>
    >>> fp = ds.open('http://www.scipy.org/not/real/data.txt.gz')
    >>> fp.read()    # Use the file
    >>> fp.close()
    >>> del ds, fp

"""

__docformat__ = "restructuredtext en"

import bz2
import gzip
import os
import tempfile
from shutil import rmtree
from urllib2 import urlopen, URLError
from urlparse import urlparse

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

# TODO: .zip support, .tar support?
_file_openers = {".gz":gzip.open, ".bz2":bz2.BZ2File, None:file}


class DataSource (object):
    """A generic data source (file, http, ftp, ...).

    DataSource could be from a local file or remote file/URL.  The file may
    also be compressed or uncompressed.

    Ex URL DataSources:
        Initialize DataSource with a local directory.  Default is os.curdir

        >>> ds = DataSource('/home/guido')
        >>> ds.open('http://fake.xyz.web/site/xyz.txt')

        Opened file exists here:  /home/guido/site/xyz.txt

    Ex using DataSource for temporary files:
        Initialize DataSource with 'None' for local directory.

        >>> ds = DataSource(None)
        >>> ds.open('/home/guido/foobar.txt')

        Opened file exists in tempdir like: /tmp/tmpUnhcvM/foobar.txt

    """

    def __init__(self, destpath=os.curdir):
        if destpath:
            self._destpath = destpath
            self._istmpdest = False
        else:
            self._destpath = tempfile.mkdtemp()
            self._istmpdest = True

    def __del__(self):
        if self._istmpdest:
            rmtree(self._destpath)

    def _iszip(self, filename):
        """Test if the filename is a zip file by looking at the file extension.
        """
        fname, ext = os.path.splitext(filename)
        return ext in _file_openers.keys()

    def _iswritemode(self, mode):
        """Test if the given mode will open a file for writing."""

        # Currently only used to test the bz2 files.  Not thoroughly tested!
        _writemodes = ("w", "+")
        for c in mode:
            if c in _writemodes:
                return True
        return False

    def _splitzipext(self, filename):
        """Split zip extension from filename and return filename.

        Returns:
            base, zip_ext : {tuple}

        """

        if self._iszip(filename):
            return os.path.splitext(filename)
        else:
            return filename, None

    def _possible_names(self, filename):
        """Return a tuple containing compressed filename variations."""
        names = [filename]
        if not self._iszip(filename):
            for zipext in _file_openers.keys():
                if zipext:
                    names.append(filename+zipext)
        return names

    def _isurl(self, path):
        """Test if path is a net location.  Tests the scheme and netloc."""
        scheme, netloc, upath, uparams, uquery, ufrag = urlparse(path)
        return bool(scheme and netloc)

    def _cache(self, path):
        """Cache the file specified by path.

        Creates a copy of the file in the datasource cache.

        """

        upath = self.abspath(path)

        # ensure directory exists
        if not os.path.exists(os.path.dirname(upath)):
            os.makedirs(os.path.dirname(upath))

        # TODO: Doesn't handle compressed files!
        if self._isurl(path):
            try:
                openedurl = urlopen(path)
                file(upath, 'w').write(openedurl.read())
            except URLError:
                raise URLError("URL not found: ", path)
        else:
            try:
                # TODO: Why not just copy the file with shutils.copyfile?
                fp = file(path, 'r')
                file(upath, 'w').write(fp.read())
            except IOError:
                raise IOError("File not found: ", path)
        return upath

    def _findfile(self, path):
        """Searches for path and returns full path if found.

        If path is an URL, _findfile will cache a local copy and return
        the path to the cached file.
        If path is a local file, _findfile will return a path to that local
        file.

        The search will include possible compressed versions of the file and
        return the first occurence found.

        """

        # Build list of possible local file paths
        filelist = self._possible_names(self.abspath(path))
        if self._isurl(path):
            # Add list of possible remote urls
            filelist = filelist + self._possible_names(path)

        for name in filelist:
            if self.exists(name):
                if self._isurl(name):
                    name = self._cache(name)
                return name
        return None

    def abspath(self, path):
        """Return an absolute path in the DataSource destination directory.

        """

        # handle case where path includes self._destpath
        splitpath = path.split(self._destpath, 2)
        if len(splitpath) > 1:
            path = splitpath[1]
        scheme, netloc, upath, uparams, uquery, ufrag = urlparse(path)
        return os.path.join(self._destpath, netloc, upath.strip(os.sep))

    def exists(self, path):
        """Test if path exists.

        Tests for local files, locally cached URLs and remote URLs.

        """

        upath = self.abspath(path)
        if os.path.exists(upath):
            return True
        elif self._isurl(path):
            try:
                netfile = urlopen(path)
                # just validate existence, nothing more.
                del(netfile)
                return True
            except URLError:
                return False
        else:
            return False

    def open(self, path, mode='r'):
        """Open path and return file object.

        If path is an URL, it will be downloaded, stored in the DataSource
        directory and opened.

        TODO: Currently only opening for reading has been tested.  There is no
              support for opening a file for writing which doesn't exist yet
              (creating a file).

        """

        if self._isurl(path) and self._iswritemode(mode):
            raise ValueError("URLs are not writeable")

        # NOTE: _findfile will fail on a new file opened for writing.
        found = self._findfile(path)
        if found:
            _fname, ext = self._splitzipext(found)
            if ext == 'bz2':
                mode.replace("+", "")
            return _file_openers[ext](found, mode=mode)
        else:
            raise IOError("%s not found." % path)


class Repository (DataSource):
    """A data repository where multiple DataSource's share one base URL.

    Use a Repository when you will be working with multiple files from one
    base URL or directory.  Initialize the Respository with the base URL,
    then refer to each file only by it's filename.

    >>> repos = Repository('/home/user/data/dir/')
    >>> fp = repos.open('data01.txt')
    >>> fp.analyze()
    >>> fp.close()

    Similarly you could use a URL for a repository:
    >>> repos = Repository('http://www.xyz.edu/data')

    """

    def __init__(self, baseurl, destpath=os.curdir):
        DataSource.__init__(self, destpath=destpath)
        self._baseurl = baseurl

    def _fullpath(self, path):
        '''Return complete path for path.  Prepends baseurl if necessary.'''
        #print 'Repository._fullpath:', path
        #print '          ._baseurl: ', self._baseurl
        splitpath = path.split(self._baseurl, 2)
        if len(splitpath) == 1:
            result = os.path.join(self._baseurl, path)
        else:
            result = path    # path contains baseurl already
        return result

    def _findfile(self, path):
        #print 'Repository._findfile:', path
        return DataSource._findfile(self, self._fullpath(path))

    def abspath(self, path):
        return DataSource.abspath(self, self._fullpath(path))

    def exists(self, path):
        #print 'Respository.exists:', path
        return DataSource.exists(self, self._fullpath(path))

    def open(self, path, mode='r'):
        #print 'Repository.open:', path
        return DataSource.open(self, self._fullpath(path), mode)

    def listdir(self):
        '''List files in the source Repository.'''
        if self._isurl(self._baseurl):
            raise NotImplementedError
        else:
            return os.listdir(self._baseurl)
