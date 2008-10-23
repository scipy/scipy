import _nnls
from numpy import asarray_chkfinite, zeros, double

def nnls(A,b):
    """ 
          Solve  || Ax - b ||_2 -> min  with  x>=0 

          Inputs:
                    A     --   matrix as above
                    b     --   vector as above

          Outputs:
                    x     --   solution vector
                    rnorm --  residual || Ax-b ||_2


    wrapper around NNLS.F code below nnls/ directory

    """

    A,b = map(asarray_chkfinite, (A,b))

    if len(A.shape)!=2:
        raise ValueError, "expected matrix"
    if len(b.shape)!=1:
        raise ValueError, "expected vector"
        
    m,n = A.shape

    if m != b.shape[0]:
        raise ValueError, "incompatible dimensions"

    w   = zeros((n,), dtype=double)
    zz  = zeros((m,), dtype=double)
    index=zeros((n,), dtype=int)

    x,rnorm,mode = _nnls.nnls(A,m,n,b,w,zz,index)
    if mode != 1: raise RuntimeError, "too many iterations"

    return x, rnorm

