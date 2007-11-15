import distutils
from distutils.core import setup, Extension
import distutils.dep_util
import os

def fortran_extension(module_name, c_files, fortran_files, library_dirs,
                      libraries):
    fcompiler = f90_compiler()
    library_name = module_name + '_fortran'
    fcompiler.build_library(library_name,fortran_files)
    libraries.append(library_name)
    ext = Extension(module_name, c_files,
                             library_dirs=library_dirs_list,
                             libraries = libraries_list)
    return ext

class f90_compiler:
    def __init__(self):
        self.compiler_name = 'g77'
    def to_object(self,dirty_files):
        files = " ".join(dirty_files)
        cmd = self.compiler_name + ' -c ' + files
        print cmd
        failure = os.system(cmd)
        if failure:
            raise ValueError, 'failure during compile'
    def object_to_library(self,library_name,object_files):
        objects = " ".join(object_files)
        cmd = 'ar -cr lib%s.a %s' % (library_name,objects)
        print cmd
        os.system(cmd)
        cmd = 'ranlib lib%s.a' % library_name
        print cmd
        os.system(cmd)
    def build_library(self,library_name,source_list):

        object_list = map(lambda x: x[:-1] +'o',source_list)
        file_pairs = zip(source_list,object_list)
        dirty_files = []
        for source,object in file_pairs:
            if distutils.dep_util.newer(source,object):
                dirty_files.append(source)
        if dirty_files != []:
            self.to_object(dirty_files)
        self.object_to_library(library_name,object_list)

if __name__ == "__main__":
    import setup # this file
    d,f = os.path.split(setup.__file__)
    print d,f
    files = os.listdir(os.path.abspath(d))
    source_files = filter(lambda x: x[-1:] == 'f',files)
    source_files = map(lambda x: os.path.abspath(x),source_files)
    compiler = f90_compiler()
    compiler.build_library('common',source_files)
