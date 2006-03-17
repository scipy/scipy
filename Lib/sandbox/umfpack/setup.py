# 05.12.2005, c
def setup_package():
    from distutils.core import setup, Extension
    from glob import glob
    import os
    import os.path as op

    moduleName = 'umfpack'

    umfpackDir = op.normpath( '/home/share/software/packages/UMFPACK/UMFPACK/' )
    umfpackInclude = op.join( umfpackDir, 'Include' )
    umfpackLib = op.join( umfpackDir, 'Lib' )
    amdLib = op.join( umfpackDir, '../AMD/Lib' )

    scipyDir = op.normpath( '/home/share/software/' )
    scipyInclude = op.join( scipyDir, 'include' )

    ##
    # The code below does not work - I don't know how to tell distutils
    # to include the generated shadow class file _umfpack,py.
##     _umfpack \
##       = Extension( '__umfpack',
##                    sources = ['umfpack.i'],
##                    include_dirs = [umfpackInclude, scipyInclude],
##                    swig_opts = ['-noruntime',  '-python', '-I'
##                                 + umfpackInclude],
##                    libraries = ['swigpy'],
##                    extra_objects = [op.join( umfpackLib, 'libumfpack.a' )] )

    ##
    # So this is the work-around... (*)
    swig_cmd = 'swig -noruntime -python -I%s %s.i' % (umfpackInclude,
                                                      moduleName)
    print 'running SWIG:', swig_cmd
    os.system( swig_cmd )

    _umfpack\
      = Extension( '__umfpack',
                   sources = ['umfpack_wrap.c'],
                   include_dirs = [umfpackInclude, scipyInclude],
                   libraries = ['swigpy', 'cblas'],
                   extra_objects = [op.join( umfpackLib, 'libumfpack.a' ),
                                    op.join( amdLib, 'libamd.a' )] )


    setup( name = moduleName,
           description = 'Python bindings for UMFPACK v4.4',
           version = '4.4.0',
           author = 'Robert Cimrman',
           author_email = 'cimrman3 (at) ntc (dot) zcu (dot) cz',
           py_modules = ['umfpack', '_umfpack'], # _umfpack (*)
           ext_modules = [_umfpack] )

if __name__ == '__main__':
    setup_package()
