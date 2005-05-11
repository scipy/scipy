import scipy_base.numerix as _nx

if _nx.which[0] == "numeric":
    from _cephes import *
else:
    from na_cephes import *

