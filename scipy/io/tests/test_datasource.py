
import bz2
import gzip
import os
import sys
import struct
from tempfile import mkdtemp, mkstemp
from shutil import rmtree
from urlparse import urlparse

from numpy.testing import *

# HACK: import the src datasource
#       Until datasource is robust enough to include in scipy.io package
sys.path.insert(0, os.path.abspath('..'))
import datasource
del sys.path[0]

#set_package_path()
#from scipy.io import datasource
#restore_path()

# Can rebind urlopen for testing, so we don't open real net connection.
#def urlopen(url, data=None):
#    print 'test_datasource urlopen(', url, ')'
#datasource.urlopen = urlopen

#http_baseurl = 'http://nifti.nimh.nih.gov/nifti-1/data/'
#http_filename = 'minimal.nii.gz'
#http_abspath = os.path.join(http_baseurl, http_filename)

# Temporarily use one of our files so we don't abuse someone elses server.
http_path = 'https://cirl.berkeley.edu/twiki/pub/BIC/ImagingDocuments/'
http_file = 'dork.pdf'

http_fakepath = 'http://fake.abc.web/site/'
http_fakefile = 'fake.txt'

magic_line = 'three is the magic number'


# Utility functions used by many TestCases
def valid_textfile(filedir):
    # Generate and return a valid temporary file.
    fd, path = mkstemp(suffix='.txt', dir=filedir, text=True)
    os.close(fd)
    return path

def invalid_textfile(filedir):
    # Generate and return an invalid filename.
    fd, path = mkstemp(suffix='.txt', dir=filedir)
    os.close(fd)
    os.remove(path)
    return path

def valid_httpurl():
    return http_path+http_file

def invalid_httpurl():
    return http_fakepath+http_fakefile

def valid_baseurl():
    return http_path

def invalid_baseurl():
    return http_fakepath

def valid_httpfile():
    return http_file

def invalid_httpfile():
    return http_fakefile

class TestDataSourceOpen(NumpyTestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.ds = datasource.DataSource(self.tmpdir)

    def tearDown(self):
        rmtree(self.tmpdir)
        del self.ds

    def test_ValidHTTP(self):
        assert self.ds.open(valid_httpurl())
    
    def test_InvalidHTTP(self):
        self.assertRaises(IOError, self.ds.open, invalid_httpurl())

    def test_ValidFile(self):
        local_file = valid_textfile(self.tmpdir)
        #print '\nDataSourceOpen test_ValidFile:', local_file
        assert self.ds.open(local_file)

    def test_InvalidFile(self):
        invalid_file = invalid_textfile(self.tmpdir)
        #print '\nDataSourceOpen test_InvalidFile:', invalid_file
        self.assertRaises(IOError, self.ds.open, invalid_file)

    def test_ValidGzipFile(self):
        # Test datasource's internal file_opener for Gzip files.
        filepath = os.path.join(self.tmpdir, 'foobar.txt.gz')
        #print '\nDataSourceOpen test_ValidGzipFile:', filepath
        fp = gzip.open(filepath, 'w')
        fp.write(magic_line)
        fp.close()
        fp = self.ds.open(filepath)
        result = fp.readline()
        fp.close()
        self.assertEqual(magic_line, result)

    def test_ValidBz2File(self):
        # Test datasource's internal file_opener for BZip2 files.
        filepath = os.path.join(self.tmpdir, 'foobar.txt.bz2')
        #print '\nDataSourceOpen test_ValidBZ2File:', filepath
        fp = bz2.BZ2File(filepath, 'w')
        fp.write(magic_line)
        fp.close()
        fp = self.ds.open(filepath)
        result = fp.readline()
        fp.close()
        self.assertEqual(magic_line, result)


class TestDataSourceExists(NumpyTestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.ds = datasource.DataSource(self.tmpdir)

    def tearDown(self):
        rmtree(self.tmpdir)
        del self.ds

    def test_ValidHTTP(self):
        #print 'DataSourceExists test_ValidHTTP'
        assert self.ds.exists(valid_httpurl())
    
    def test_InvalidHTTP(self):
        #print 'DataSourceExists test_InvalidHTTP'
        self.assertEqual(self.ds.exists(invalid_httpurl()), False)

    def test_ValidFile(self):
        tmpfile = valid_textfile(self.tmpdir)
        #print 'DataSourceExists test_ValidFile:', tmpfile
        assert self.ds.exists(tmpfile)

    def test_InvalidFile(self):
        tmpfile = invalid_textfile(self.tmpdir)
        #print 'DataSourceExists test_InvalidFile:', tmpfile
        self.assertEqual(self.ds.exists(tmpfile), False)


class TestDataSourceAbspath(NumpyTestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.ds = datasource.DataSource(self.tmpdir)

    def tearDown(self):
        rmtree(self.tmpdir)
        del self.ds

    def test_ValidHTTP(self):
        scheme, netloc, upath, pms, qry, frg = urlparse(valid_httpurl())
        local_path = os.path.join(self.tmpdir, netloc, upath.strip(os.sep))
        self.assertEqual(local_path, self.ds.abspath(valid_httpurl()))

    def test_ValidFile(self):
        tmpfile = valid_textfile(self.tmpdir)
        tmpfilename = os.path.split(tmpfile)[-1]
        # Test with filename only
        self.assertEqual(tmpfile, self.ds.abspath(os.path.split(tmpfile)[-1]))
        # Test filename with complete path
        self.assertEqual(tmpfile, self.ds.abspath(tmpfile))

    def test_InvalidHTTP(self):
        scheme, netloc, upath, pms, qry, frg = urlparse(invalid_httpurl())
        invalidhttp = os.path.join(self.tmpdir, netloc, upath.strip(os.sep))
        self.assertNotEqual(invalidhttp, self.ds.abspath(valid_httpurl()))

    def test_InvalidFile(self):
        fd, invalidfile = mkstemp(suffix='.txt')
        os.close(fd)
        tmpfile = valid_textfile(self.tmpdir)
        tmpfilename = os.path.split(tmpfile)[-1]
        # Test with filename only
        self.assertNotEqual(invalidfile, self.ds.abspath(tmpfilename))
        # Test filename with complete path
        self.assertNotEqual(invalidfile, self.ds.abspath(tmpfile))


class TestRespositoryAbspath(NumpyTestCase):
    def setUp(self):
        self.repos = datasource.Repository(valid_baseurl(), None)

    def tearDown(self):
        del self.repos

    def test_ValidHTTP(self):
        scheme, netloc, upath, pms, qry, frg = urlparse(valid_httpurl())
        local_path = os.path.join(self.repos._destpath, netloc, \
                                  upath.strip(os.sep))
        filepath = self.repos.abspath(valid_httpfile())
        self.assertEqual(local_path, filepath)


if __name__ == "__main__":
    NumpyTest().run()

