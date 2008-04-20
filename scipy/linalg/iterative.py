__all__ = ['bicg','bicgstab','cg','cgs','gmres','qmr']

# Deprecated on January 26, 2008

from scipy.sparse.linalg import isolve
from numpy import deprecate

for name in __all__:
    oldfn = getattr(isolve, name)
    oldname='scipy.linalg.' + name
    newname='scipy.sparse.linalg.' + name
    newfn = deprecate(oldfn, oldname=oldname, newname=newname)
    exec(name + ' = newfn')
