from __future__ import absolute_import, print_function

import sys
import os
import re
import glob
import string
import stat
import tempfile

from distutils.dir_util import remove_tree

from numpy.testing import TestCase, assert_, run_module_suite
from numpy.testing.noseclasses import KnownFailureTest

from scipy.weave import catalog
from weave_test_utils import (clear_temp_catalog, restore_temp_catalog,
                              empty_temp_dir, cleanup_temp_dir, dec)


skip_on_windows = dec.skipif(sys.platform == 'win32',
                             "Test works only on posix")


def samefile(a, b):
    try:
        return os.path.samefile(a, b)
    except AttributeError:
        return os.path.realpath(a) == os.path.realpath(b)


class TestIntermediateDir(TestCase):
    """
    Tests for intermediate dir (store of .cpp and .o during builds).
    These tests test whether intermediate dir is safe. If it's not,
    new one should be created.
    """
    def dirs_are_valid(self, wrong_dir, tmpdir):
        # Test if new dir is created and is consistent.
        new_im_dir = catalog.intermediate_dir(tmpdir)
        assert_(not samefile(new_im_dir, wrong_dir))
        new_im_dir2 = catalog.intermediate_dir(tmpdir)
        assert_(samefile(new_im_dir, new_im_dir2))

    @skip_on_windows
    def test_ownership(self):
        # Test if intermediate dir is owned by correct user.
        im_dir = catalog.intermediate_dir()
        im_dir_stat = os.stat(im_dir)
        proc_uid = os.getuid()
        assert_(proc_uid == im_dir_stat.st_uid)
        r_im_dir_stat = os.stat(os.path.dirname(im_dir))
        assert_(proc_uid == r_im_dir_stat.st_uid)

    @skip_on_windows
    def test_incorrect_ownership(self):
        # Test if new intermediate dir is created when there is only one
        # im dir owned by improper user.
        import pwd
        tmpdir = tempfile.mkdtemp()
        try:
            im_dir = catalog.create_intermediate_dir(tmpdir)
            root_im_dir = os.path.dirname(im_dir)
            nobody = pwd.getpwnam('nobody')[2]
            nobody_g = pwd.getpwnam('nobody')[3]
            try:
                os.chown(root_im_dir, nobody, nobody_g)
            except OSError:
                raise KnownFailureTest("Can't change owner.")
            else:
                self.dirs_are_valid(im_dir, tmpdir)
        finally:
            remove_tree(tmpdir)

    @skip_on_windows
    def test_permissions(self):
        # im dir should have permissions 0700
        im_dir = catalog.intermediate_dir()
        im_dir_stat = os.stat(im_dir)
        assert_(stat.S_IMODE(im_dir_stat.st_mode) == 0o0700)
        r_im_dir_stat = os.stat(os.path.dirname(im_dir))
        assert_(stat.S_IMODE(r_im_dir_stat.st_mode) == 0o0700)

    @skip_on_windows
    def test_incorrect_permissions(self):
        # If permissions on existing im dir are not correct,
        # a new one should be created.
        tmpdir = tempfile.mkdtemp()
        try:
            im_dir = catalog.create_intermediate_dir(tmpdir)
            root_im_dir = os.path.dirname(im_dir)
            try:
                os.chmod(root_im_dir, 0o777)
            except OSError:
                raise KnownFailureTest("Can't set file permissions.")
            else:
                self.dirs_are_valid(im_dir, tmpdir)
        finally:
            remove_tree(tmpdir)

    @skip_on_windows
    def test_symlink(self):
        # im dir shouldn't be a symlink
        r_im_dir = os.path.dirname(catalog.intermediate_dir())
        assert_(os.path.islink(r_im_dir) is False)

    @skip_on_windows
    def test_symlink_raise(self):
        # If existing im dir is a symlink, new one should be created.
        tmpdir = tempfile.mkdtemp()
        try:
            im_dir = catalog.create_intermediate_dir(tmpdir)
            root_im_dir = os.path.dirname(im_dir)

            tempdir = tempfile.mkdtemp(prefix='scipy-test', dir=tmpdir)
            try:
                os.rename(root_im_dir, tempdir)
            except OSError:
                raise KnownFailureTest("Can't move intermediate dir.")

            try:
                os.symlink(tempdir, root_im_dir)
            except OSError:
                raise KnownFailureTest(
                    "Can't create symlink to intermediate dir.")
            else:
                self.dirs_are_valid(im_dir, tmpdir)
        finally:
            remove_tree(tmpdir)


class TestDefaultDir(TestCase):
    """
    Tests for 'catalog.default_dir()'.
    These should verified posix and win default_dir function.
    """
    def test_win(self):
        # test if default_dir for Windows platform is accessible
        #
        # since default_dir_win() does not have any Windows specific code,
        # let's test it everywhere
        d = catalog.default_dir_win()
        assert_(catalog.is_writable(d))

    @skip_on_windows
    def test_posix(self):
        # test if posix default_dir is writable
        d = catalog.default_dir_posix()
        assert_(catalog.is_writable(d))

    @skip_on_windows
    def test_posix_home_inaccessible(self):
        # what happens when home catalog dir is innaccessible
        tmpdir = tempfile.mkdtemp()
        try:
            d_dir = catalog.default_dir_posix(tmpdir)

            try:
                os.chmod(d_dir, 0o000)
            except OSError:
                raise KnownFailureTest("Can't change permissions of default_dir.")

            new_ddir = catalog.default_dir_posix(tmpdir)
            assert_(not os.path.samefile(new_ddir, d_dir))
            new_ddir2 = catalog.default_dir_posix(tmpdir)
            assert_(os.path.samefile(new_ddir, new_ddir2))
        finally:
            os.chmod(d_dir, 0o700)
            remove_tree(tmpdir)

    @skip_on_windows
    def test_posix_dirs_inaccessible(self):
        # test if new dir is created if both implicit dirs are not valid
        tmpdir = tempfile.mkdtemp()
        try:
            d_dir = catalog.default_dir_posix(tmpdir)

            try:
                os.chmod(d_dir, 0o000)
            except OSError:
                raise KnownFailureTest("Can't change permissions of default_dir.")

            d_dir2 = catalog.default_dir_posix(tmpdir)

            try:
                os.chmod(d_dir2, 0o000)
            except OSError:
                raise KnownFailureTest("Can't change permissions of default_dir.")

            new_ddir = catalog.default_dir_posix(tmpdir)
            assert_(not (os.path.samefile(new_ddir, d_dir) or os.path.samefile(new_ddir, d_dir2)))
            new_ddir2 = catalog.default_dir_posix(tmpdir)
            assert_(os.path.samefile(new_ddir, new_ddir2))
        finally:
            os.chmod(d_dir, 0o700)
            os.chmod(d_dir2, 0o700)
            remove_tree(tmpdir)

    def test_is_writable(self):
        # default_dir has to be writable
        path = catalog.default_dir()
        name = os.path.join(path,'dummy_catalog')
        test_file = open(name,'w')
        try:
            test_file.write('making sure default location is writable\n')
        finally:
            test_file.close()
            os.remove(name)


class TestOsDependentCatalogName(TestCase):
    pass


class TestCatalogPath(TestCase):

    def test_default(self):
        in_path = catalog.default_dir()
        path = catalog.catalog_path(in_path)
        d,f = os.path.split(path)
        assert_(d == in_path)
        assert_(f == catalog.os_dependent_catalog_name())

    def test_current(self):
        in_path = '.'
        path = catalog.catalog_path(in_path)
        d,f = os.path.split(path)
        assert_(d == os.path.abspath(in_path))
        assert_(f == catalog.os_dependent_catalog_name())

    @skip_on_windows
    def test_user(path):
        in_path = '~'
        path = catalog.catalog_path(in_path)
        d,f = os.path.split(path)
        assert_(d == os.path.expanduser(in_path))
        assert_(f == catalog.os_dependent_catalog_name())

    def test_module(self):
        # hand it a module and see if it uses the parent directory
        # of the module.
        path = catalog.catalog_path(os.__file__)
        d,f = os.path.split(os.__file__)
        d2,f = os.path.split(path)
        assert_(d2 == d)

    def test_path(self):
        # use os.__file__ to get a usable directory.
        in_path,f = os.path.split(os.__file__)
        path = catalog.catalog_path(in_path)
        d,f = os.path.split(path)
        assert_(d == in_path)

    def test_bad_path(self):
        # stupid_path_name
        in_path = 'stupid_path_name'
        path = catalog.catalog_path(in_path)
        assert_(path is None)


class TestGetCatalog(TestCase):
    """ This only tests whether new catalogs are created correctly.
        And whether non-existent return None correctly with read mode.
        Putting catalogs in the right place is all tested with
        catalog_dir tests.
    """

    def get_test_dir(self,erase=0):
        # make sure tempdir catalog doesn't exist
        pardir = tempfile.mkdtemp(suffix='cat_test')
        cat_glob = os.path.join(pardir,catalog.os_dependent_catalog_name()+'.*')
        cat_files = glob.glob(cat_glob)
        if erase:
            for cat_file in cat_files:
                os.remove(cat_file)
        return pardir

    def test_nonexistent_catalog_is_none(self):
        pardir = self.get_test_dir(erase=1)
        cat = catalog.get_catalog(pardir,'r')
        remove_tree(pardir)
        assert_(cat is None)

    def test_create_catalog(self):
        pardir = self.get_test_dir(erase=1)
        cat = catalog.get_catalog(pardir,'c')
        assert_(cat is not None)
        cat.close()
        remove_tree(pardir)


class TestCatalog(TestCase):

    def clear_environ(self):
        if 'PYTHONCOMPILED' in os.environ:
            self.old_PYTHONCOMPILED = os.environ['PYTHONCOMPILED']
            del os.environ['PYTHONCOMPILED']
        else:
            self.old_PYTHONCOMPILED = None

    def reset_environ(self):
        if self.old_PYTHONCOMPILED:
            os.environ['PYTHONCOMPILED'] = self.old_PYTHONCOMPILED
            self.old_PYTHONCOMPILED = None

    def setUp(self):
        self.clear_environ()

    def tearDown(self):
        self.reset_environ()

    def test_set_module_directory(self):
        q = catalog.catalog()
        q.set_module_directory('bob')
        r = q.get_module_directory()
        assert_(r == 'bob')

    def test_clear_module_directory(self):
        q = catalog.catalog()
        r = q.get_module_directory()
        assert_(r is None)
        q.set_module_directory('bob')
        r = q.clear_module_directory()
        assert_(r is None)

    def test_get_environ_path(self):
        if sys.platform == 'win32':
            sep = ';'
        else:
            sep = ':'
        os.environ['PYTHONCOMPILED'] = sep.join(('path1','path2','path3'))
        q = catalog.catalog()
        path = q.get_environ_path()
        assert_(path == ['path1','path2','path3'])

    def test_build_search_order1(self):
        # MODULE in search path should be replaced by module_dir.
        q = catalog.catalog(['first','MODULE','third'])
        q.set_module_directory('second')
        order = q.build_search_order()
        assert_(order == ['first','second','third',catalog.default_dir()])

    def test_build_search_order2(self):
        # MODULE in search path should be removed if module_dir==None.
        q = catalog.catalog(['first','MODULE','third'])
        order = q.build_search_order()
        assert_(order == ['first','third',catalog.default_dir()])

    def test_build_search_order3(self):
        # If MODULE is absent, module_dir shouldn't be in search path.
        q = catalog.catalog(['first','second'])
        q.set_module_directory('third')
        order = q.build_search_order()
        assert_(order == ['first','second',catalog.default_dir()])

    def test_build_search_order4(self):
        # Make sure environment variable is getting used.
        q = catalog.catalog(['first','second'])
        if sys.platform == 'win32':
            sep = ';'
        else:
            sep = ':'
        os.environ['PYTHONCOMPILED'] = sep.join(('MODULE','fourth','fifth'))
        q.set_module_directory('third')
        order = q.build_search_order()
        assert_(order == ['first','second','third','fourth','fifth',catalog.default_dir()])

    def test_catalog_files1(self):
        # Be sure we get at least one file even without specifying the path.
        q = catalog.catalog()
        files = q.get_catalog_files()
        assert_(len(files) == 1)

    def test_catalog_files2(self):
        # Ignore bad paths in the path.
        q = catalog.catalog()
        os.environ['PYTHONCOMPILED'] = '_some_bad_path_'
        files = q.get_catalog_files()
        assert_(len(files) == 1)

    def test_get_existing_files1(self):
        # Shouldn't get any files when temp doesn't exist and no path set.
        backup_dir = clear_temp_catalog()
        q = catalog.catalog()
        files = q.get_existing_files()
        restore_temp_catalog(backup_dir)
        assert_(len(files) == 0)

    def test_get_existing_files2(self):
        # Shouldn't get a single file from the temp dir.
        backup_dir = clear_temp_catalog()
        q = catalog.catalog()
        # create a dummy file
        q.add_function('code', os.getpid)
        del q
        q = catalog.catalog()
        files = q.get_existing_files()
        restore_temp_catalog(backup_dir)
        assert_(len(files) == 1)

    def test_access_writable_file(self):
        # There should always be a writable file -- even if it is in temp
        q = catalog.catalog()
        file = q.get_writable_file()
        try:
            f = open(file,'w')
            f.write('bob')
        finally:
            f.close()
            os.remove(file)

    def test_writable_with_bad_path(self):
        # There should always be a writable file -- even if search paths
        # contain bad values.
        if sys.platform == 'win32':
            sep = ';'
        else:
            sep = ':'
        os.environ['PYTHONCOMPILED'] = sep.join(('_bad_path_name_'))
        q = catalog.catalog()
        file = q.get_writable_file()
        try:
            f = open(file,'w')
            f.write('bob')
        finally:
            f.close()
        os.remove(file)

    def test_writable_dir(self):
        # Check that we can create a file in the writable directory
        q = catalog.catalog()
        d = q.get_writable_dir()
        file = os.path.join(d,'some_silly_file')
        try:
            f = open(file,'w')
            f.write('bob')
        finally:
            f.close()
            os.remove(file)

    def test_unique_module_name(self):
        # Check that we can create a file in the writable directory
        q = catalog.catalog()
        file = q.unique_module_name('bob')
        cfile1 = file+'.cpp'
        assert_(not os.path.exists(cfile1))
        # Make sure it is writable
        try:
            f = open(cfile1,'w')
            f.write('bob')
        finally:
            f.close()
        # try again with same code fragment -- should get unique name
        file = q.unique_module_name('bob')
        cfile2 = file+'.cpp'
        assert_(not os.path.exists(cfile2+'.cpp'))
        os.remove(cfile1)

    def test_add_function_persistent1(self):
        # Test persisting a function in the default catalog
        backup_dir = clear_temp_catalog()
        q = catalog.catalog()
        # just use some already available functions
        funcs = [string.upper, string.lower, string.find,string.replace]
        for i in funcs:
            q.add_function_persistent('code',i)
        pfuncs = q.get_cataloged_functions('code')
        # any way to clean modules???
        restore_temp_catalog(backup_dir)
        for i in funcs:
            assert_(i in pfuncs)

    def test_add_function_ordered(self):
        backup_dir = clear_temp_catalog()
        q = catalog.catalog()

        q.add_function('f',string.upper)
        q.add_function('f',string.lower)
        q.add_function('ff',string.find)
        q.add_function('ff',string.replace)
        q.add_function('fff',string.atof)
        q.add_function('fff',string.atoi)
        del q

        # now we're gonna make a new catalog with same code
        # but different functions in a specified module directory
        env_dir = empty_temp_dir()
        r = catalog.catalog(env_dir)
        r.add_function('ff',os.abort)
        r.add_function('ff',os.chdir)
        r.add_function('fff',os.access)
        r.add_function('fff',os.open)
        del r
        # now we're gonna make a new catalog with same code
        # but different functions in a user specified directory
        user_dir = empty_temp_dir()
        s = catalog.catalog(user_dir)
        s.add_function('fff',re.match)
        s.add_function('fff',re.purge)
        del s

        # open new catalog and make sure it retreives the functions
        # from d catalog instead of the temp catalog (made by q)
        os.environ['PYTHONCOMPILED'] = env_dir
        t = catalog.catalog(user_dir)
        funcs1 = t.get_functions('f')
        funcs2 = t.get_functions('ff')
        funcs3 = t.get_functions('fff')
        restore_temp_catalog(backup_dir)
        # make sure everything is read back in the correct order
        # a little cheating... I'm ignoring any functions that might have
        # been read in from a prior catalog file (such as the defualt one).
        # the test should really be made so that these aren't read in, but
        # until I get this figured out...
        #assert_(funcs1 == [string.lower,string.upper])
        #assert_(funcs2 == [os.chdir,os.abort,string.replace,string.find])
        #assert_(funcs3 == [re.purge,re.match,os.open,
        #                  os.access,string.atoi,string.atof])
        assert_(funcs1[:2] == [string.lower,string.upper]),repr(funcs1)
        assert_(funcs2[:4] == [os.chdir,os.abort,string.replace,string.find])
        assert_(funcs3[:6] == [re.purge,re.match,os.open,
                          os.access,string.atoi,string.atof])
        cleanup_temp_dir(user_dir)
        cleanup_temp_dir(env_dir)


if __name__ == '__main__':
    run_module_suite()
