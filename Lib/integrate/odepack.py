import _odepack

def odeint(func, y0, t, args=(), Dfun=None, col_deriv=0, ml=None, mu=None,
           full_output=0, rtol=None, atol=None, tcrit=None, h0=0.0, hmax=0.0,
           hmin=0.0, ixpr=0, mxstep=0, mxhnil=0, mxordn=12, mxords=5):

    kwdargs = {}
    return _odepack.odeint(func, y0, t, args, Dfun, col_deriv, ml, mu, full_output,
                           rtol, atol, tcrit, h0, hmax, hmin, ixpr, mxstep, mxhnil,
                           mxordn, mxords)



