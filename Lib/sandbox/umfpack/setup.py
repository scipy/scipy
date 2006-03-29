# 05.12.2005, c
# last change: 27.03.2006
def setup_package():
    from distutils.core import setup, Extension
    from glob import glob
    import os
    import os.path as op
    import numpy

    moduleName = 'umfpack'

    umfpackDir = op.normpath( '/home/share/software/packages/UMFPACK/UMFPACK/' )
    umfpackInclude = op.join( umfpackDir, 'Include' )
    umfpackLib = op.join( umfpackDir, 'Lib' )
    amdLib = op.join( umfpackDir, '../AMD/Lib' )

    scipyInclude = numpy.get_numpy_include()

    _umfpack \
      = Extension( '__umfpack',
                   sources = [op.join( 'umfpack', ii ) for ii in ['umfpack.i']],
                   swig_opts = ['-I' + umfpackInclude],
                   include_dirs = [umfpackInclude, scipyInclude],
                   libraries = ['cblas'],
                   extra_objects = [op.join( umfpackLib, 'libumfpack.a' ),
                                    op.join( amdLib, 'libamd.a' )] )

    setup( name = moduleName,
           description = 'Python bindings for UMFPACK v4.4',
           version = '4.4.0',
           author = 'Robert Cimrman',
           author_email = 'cimrman3 (at) ntc (dot) zcu (dot) cz',
           scripts = ['test_umfpack.py'],
           ext_package  = 'umfpack',
           ext_modules = [_umfpack],
           package_dir  = {'umfpack' : 'umfpack'},
           packages = ['umfpack'] )

if __name__ == '__main__':
    setup_package()
