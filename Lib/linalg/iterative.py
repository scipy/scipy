
# Iterative methods using reverse-communication raw material
#   These methods solve
#   Ax = b  for x

#   where A must have A.matvec(x,*args) defined
#    or be a numeric array
#    or be a callable object matvec(x,*args)

__all__ = ['bicg']
import _iterative

_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}

def type1(self,x):
    return dot(self.obj, x)

def type2(self, x):
    return self.obj(x,*self.args)

class get_matvec:
    def __init__(self, obj, *args):
        self.obj = obj
        self.args = args
        if isinstance(obj, ArrayType):
            self.callfunc = new.instancemethod(type1,self,get_matvec)
            return
        try:
            meth = obj.matvec
        except AttributeError:
            meth = obj
        if not callable(meth):
            raise ValueError, "Object must be an array or a callable object"\
                  "or have a callable matvec attribute."

        self.obj = meth
        self.callfunc = new.instancemethod(type2, self, get_matvec)
        
    def __call__(self, x):
        return self.callfunc(self, x)
        


def bicg(A,b,tol=1e-5,maxiter=None):
    b = asarray(b)
    n = len(b)
    typ = b.typecode()
    if maxiter is None:
        maxiter = n*5

    matvec, psolve, matvecT, psolveT = (None,)*4
    ltr = _type_conf[typ]
    revcom = _iterative.__dict__[ltr+'bicgrevcom']
    stoptest = _iterative.__dict__[ltr+'stoptest2']

    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = zeros(6*n,typ)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    while 1:
        x, work, iter, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter, resid, info, ndx1, ndx2, ijob)
        slice1 = slice(ndx1, ndx1+n)
        slice2 = slice(ndx2, ndx2+n)
        if (ijob == -1):
            break
        elif (ijob == 1):
            if matvec is None:
                matvec = get_matvec(A)
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
        elif (ijob == 2):
            if matvectrans is None:
                matvectrans = get_matvectrans(A)
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvectrans(work[slice1])
        elif (ijob == 3):
            if psolve is None:
                psolve = get_psolve(A)
            work[slice1] = psolve(work[slice2])
        elif (ijob == 4):
            if psolvetrans is None:
                psolvetrans = get_psolvetrans(A)
            work[slice1] = psolvetrans(work[slice2])
        elif (ijob == 5):
            if matvec is None:
                matvec = get_matvec(A)
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 6):
            if ftflag:
                info = -1
                ftflag = False
            bnrm2, resid, info = stoptest(work[slice1], b, bnrm2, tol)
        ijob = 2

    return x, info
