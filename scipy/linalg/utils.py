from os.path import dirname, realpath

__all__ = ['get_linalg_dir']

def get_linalg_dir():
    """
    Return the directory containing scipy.linalg.
    
    This function should be used to get the proper include directory
    to include the proper cython pxd files when using the Cython
    wrappers for BLAS and LAPACK.
    
    Notes
    -----
    When using ``distutils``, for example, in ``setup.py``.
    ::
    
        from scipy.linalg import get_linalg_dir
        ...
        Extension('extension_name', ...
                  include_dirs=[get_linalg_dir()])
        ...
    
    """
    # Used to get the proper include directory to include f2pyptr.h
    # and the cython pxd files when using the Cython wrappers
    # for BLAS and LAPACK.
    return dirname(realpath(__file__))
