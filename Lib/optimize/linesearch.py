import minpack2
import scipy_base
import sys
if sys.version[:3] < "2.3":
    True = 1
    False = 0

import __builtin__
pymin = __builtin__.min
    

def line_search(f, myfprime, xk, pk, gfk, old_fval, old_old_fval,
                args=(), c1=1e-4, c2=0.9, amax=50):

    fc = 0
    gc = 0
    phi0 = old_fval
    derphi0 = scipy_base.dot(gfk,pk)
    alpha1 = pymin(1.0,1.01*2*(phi0-old_old_fval)/derphi0)
    
    if isinstance(myfprime,type(())):
        eps = myfprime[1]
        fprime = myfprime[0]
        newargs = (f,eps) + args
        gradient = False
    else:
        fprime = myfprime
        newargs = args
        gradient = True

    xtol = 1e-14
    amin = 1e-8
    isave = scipy_base.zeros((2,),'i')
    dsave = scipy_base.zeros((13,),'d')
    task = 'START'
    fval = old_fval
    gval = gfk
    
    while 1:
        stp,fval,derphi,task = minpack2.dcsrch(alpha1, phi0, derphi0, c1, c2,
                                               xtol, task, amin, amax,isave,dsave)
        
        if task[:2] == 'FG':
            alpha1 = stp
            fval = f(xk+stp*pk,*args)
            fc += 1
            gval = fprime(xk+stp*pk,*newargs)
            if gradient: gc += 1
            else: fc += len(xk) + 1            
            phi0 = fval
            derphi0 = scipy_base.dot(gval,pk)
        else:
            break                

    if task[:5] == 'ERROR' or task[1:4] == 'WARN':
        stp = None  # failed
    return stp, fc, gc, fval, old_fval, gval
