from __future__ import absolute_import, print_function

import os

def remove_whitespace(in_str):
    out = in_str.replace(" ","")
    out = out.replace("\t","")
    out = out.replace("\n","")
    return out

###################################################
# mainly used by catalog tests
###################################################

from scipy.weave import catalog

import glob

def temp_catalog_files(prefix=''):
    # might need to add some more platform specific catalog file
    # suffixes to remove.  The .pag was recently added for SunOS
    d = catalog.default_dir()
    f = catalog.os_dependent_catalog_name()
    return glob.glob(os.path.join(d,prefix+f+'*'))

import tempfile

def clear_temp_catalog():
    """ Remove any catalog from the temp dir
    """
    global backup_dir
    backup_dir =tempfile.mktemp()
    os.mkdir(backup_dir)
    for file in temp_catalog_files():
        move_file(file,backup_dir)
        #d,f = os.path.split(file)
        #backup = os.path.join(backup_dir,f)
        #os.rename(file,backup)

def restore_temp_catalog():
    """ Remove any catalog from the temp dir
    """
    global backup_dir
    cat_dir = catalog.default_dir()
    for file in os.listdir(backup_dir):
        file = os.path.join(backup_dir,file)
        d,f = os.path.split(file)
        dst_file = os.path.join(cat_dir, f)
        if os.path.exists(dst_file):
            os.remove(dst_file)
        #os.rename(file,dst_file)
        move_file(file,dst_file)
    os.rmdir(backup_dir)
    backup_dir = None

def empty_temp_dir():
    """ Create a sub directory in the temp directory for use in tests
    """
    import tempfile
    d = catalog.default_dir()
    for i in range(10000):
        new_d = os.path.join(d,tempfile.gettempprefix()[1:-1]+repr(i))
        if not os.path.exists(new_d):
            os.mkdir(new_d)
            break
    return new_d

def cleanup_temp_dir(d):
    """ Remove a directory created by empty_temp_dir
        should probably catch errors
    """
    files = map(lambda x,d=d: os.path.join(d,x),os.listdir(d))
    for i in files:
        try:
            if os.path.isdir(i):
                cleanup_temp_dir(i)
            else:
                os.remove(i)
        except OSError:
            pass # failed to remove file for whatever reason
                 # (maybe it is a DLL Python is currently using)
    try:
        os.rmdir(d)
    except OSError:
        pass


# from distutils -- old versions had bug, so copying here to make sure
# a working version is available.
from distutils.errors import DistutilsFileError
import distutils.file_util
def move_file (src, dst,
               verbose=0,
               dry_run=0):

    """Move a file 'src' to 'dst'.  If 'dst' is a directory, the file will
    be moved into it with the same name; otherwise, 'src' is just renamed
    to 'dst'.  Return the new full name of the file.

    Handles cross-device moves on Unix using 'copy_file()'.  What about
    other systems???
    """
    from os.path import exists, isfile, isdir, basename, dirname
    import errno

    if verbose:
        print("moving %s -> %s" % (src, dst))

    if dry_run:
        return dst

    if not isfile(src):
        raise DistutilsFileError("can't move '%s': not a regular file" % src)

    if isdir(dst):
        dst = os.path.join(dst, basename(src))
    elif exists(dst):
        raise DistutilsFileError("can't move '%s': destination '%s' already "
                                    "exists" % (src, dst))

    if not isdir(dirname(dst)):
        raise DistutilsFileError("can't move '%s': destination '%s' not a "
                                    "valid path" % (src, dst))

    copy_it = 0
    try:
        os.rename(src, dst)
    except os.error as xxx_todo_changeme1:
        (num, msg) = xxx_todo_changeme1.args
        if num == errno.EXDEV:
            copy_it = 1
        else:
            raise DistutilsFileError("couldn't move '%s' to '%s': %s" %
                                        (src, dst, msg))

    if copy_it:
        distutils.file_util.copy_file(src, dst)
        try:
            os.unlink(src)
        except os.error as xxx_todo_changeme:
            (num, msg) = xxx_todo_changeme.args
            try:
                os.unlink(dst)
            except os.error:
                pass
            raise DistutilsFileError("couldn't move '%s' to '%s' by copy/delete: "
                            "delete '%s' failed: %s" % (src, dst, src, msg))

    return dst
