""" Linear algebra routines.

 solve
 inv
 lu
 cholesky
 qr
 schur
 rsf2csf
 norm
 eig
 eigvals
 svd
 pinv
 pinv2
 det
 lstsq
 expm
 cosm
 sinm
 tanm
 coshm
 sinhm
 tanhm
 funm 
 
"""

_modules = ['fblas', 'flapack', 'cblas', 'clapack']
_namespaces = ['linear_algebra']

__all__ = []
import scipy
scipy.modules2all(__all__, _modules, globals())
scipy.names2all(__all__, _namespaces, globals())
