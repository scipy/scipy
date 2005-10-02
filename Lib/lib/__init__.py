
from info_lib import __doc__
__all__ = ['lapack','blas']


from scipy.base import ppimport
lapack = ppimport('lapack')
blas = ppimport('blas')
