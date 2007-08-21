from base import *
from log import Log

##
# 19.04.2006, c
# 26.04.2006
def convTest( conf, it, of, of0, ofgNorm = None ):
    """
    -1 ... continue
     0 ... small OF -> stop
     1 ... iMax reached -> stop
     2 ... small OFG -> stop
     3 ... small relative decrase of OF
     """

    status = -1
    print 'opt: iter: %d, of: %e (||ofg||: %e)' % (it, of, ofgNorm)
#    print (of0 - of), (conf.epsRD * of0)

    if (abs( of ) < conf.epsOF):
        status = 0
    elif ofgNorm and (ofgNorm < conf.epsOFG):
        status = 2
    elif (it > 0) and (abs(of0 - of) < (conf.epsRD * abs( of0 ))):
        status = 3
        
    if (status == -1) and (it >= conf.iMax):
        status = 1

    return status

##
# 19.04.2006, from scipy.optimize
# 21.04.2006
def wrapFunction(function, args):
    ncalls = [0]
    times = []
    def function_wrapper(x):
        ncalls[0] += 1
        tt = time.clock()
        out = function(x, *args)
        times.append( time.clock() - tt )
        return out
    return ncalls, times, function_wrapper

##
# 19.04.2006, from scipy.optimize
def vecNorm(x, ord=2):
    if ord == nm.Inf:
        return nm.amax(abs(x))
    elif ord == -nm.Inf:
        return nm.amin(abs(x))
    else:
        return nm.sum(abs(x)**ord)**(1.0/ord)

##
# 20.04.2006, c
def checkGradient( xit, aofg, fn_of, delta, check ):

    dofg = nm.zeros_like( aofg )
    xd = xit.copy()
    for ii in xrange( xit.shape[0] ):
        xd[ii] = xit[ii] + delta
        ofp = fn_of( xd )

        xd[ii] = xit[ii] - delta
        ofm = fn_of( xd )

        xd[ii] = xit[ii]

        dofg[ii] = 0.5 * (ofp - ofm) / delta

        print '**********', ii, aofg[ii], dofg[ii]

    diff = abs( aofg - dofg )
    aux = nm.concatenate( (aofg[:,nm.NewAxis], dofg[:,nm.NewAxis],
                           diff[:,nm.NewAxis]), 1 )
    print aux
    print vecNorm( diff, nm.Inf )
##     aofg.tofile( 'aofg.txt', ' ' )
##     dofg.tofile( 'dofg.txt', ' ' )
##     diff.tofile( 'diff.txt', ' ' )
    if check == 2:
        pylab.plot( aofg )
        pylab.plot( dofg )
        pylab.legend( ('analytical', 'finite difference') )
        pylab.show()
    print 'gradient checking done'
    
##
# 19.04.2006, c
# 20.04.2006
# 21.04.2006
# 26.04.2006
def fmin_sd( conf, x0, fn_of, fn_ofg, args = () ):

    nc_of, tt_of, fn_of = wrapFunction( fn_of, args )
    nc_ofg, tt_ofg, fn_ofg = wrapFunction( fn_ofg, args )

    timeStats = {'of' : tt_of, 'ofg': tt_ofg, 'check' : []}

    if conf.log:
        log = Log.fromConf( conf, (['of'], ['ofgNorm'], ['alpha']) )

    ofg = None

    it = 0
    xit = x0.copy()
    while 1:

        of = fn_of( xit )

        if it == 0:
            of0 = ofit0 = ofPrev = of
            ofPrevPrev = of + 5000.0

        if ofg is None:
            ofg = fn_ofg( xit )

        if conf.check:
            tt = time.clock()
            checkGradient( xit, ofg, fn_of, conf.delta, conf.check )
            timeStats['check'].append( time.clock() - tt )

        ofgNorm = vecNorm( ofg, conf.norm )

        status = convTest( conf, it, of, ofit0, ofgNorm )
        if status >= 0:
            break

        # These values are modified by the line search, even if it fails
        ofPrev_bak = ofPrev
        ofPrevPrev_bak = ofPrevPrev

        if conf.ls:
            alpha, fc, gc, ofPrev, ofPrevPrev, ofg1 = \
                   linesearch.line_search(fn_of,fn_ofg,xit,\
                                          -ofg,ofg,ofPrev,ofPrevPrev,c2=0.4)
            if alpha is None:  # line search failed -- use different one.
                alpha, fc, gc, ofPrev, ofPrevPrev, ofg1 = \
                       sopt.line_search(fn_of,fn_ofg,xit,\
                                        -ofg,ofg,ofPrev_bak,ofPrevPrev_bak)
                if alpha is None or alpha == 0:
                    # This line search also failed to find a better solution.
                    status = 3
                    break
        else:
            alpha = 1.0
            ofg1 = None
        
        if conf.log:
            log( of, ofgNorm, alpha )

        xit = xit - alpha * ofg
        if ofg1 is None:
            ofg = None
        else:
            ofg = ofg1.copy()
        
        for key, val in timeStats.iteritems():
            if len( val ):
                print '%10s: %7.2f [s]' % (key, val[-1])

        ofit0 = of

        it = it + 1

    print 'status:               %d' % status
    print 'initial value:        %.8e' % of0
    print 'current value:        %.8e' % of
    print 'iterations:           %d' % it
    print 'function evaluations: %d in %.2f [s]' \
          % (nc_of[0], nm.sum( timeStats['of'] ) )
    print 'gradient evaluations: %d in %.2f [s]' \
          % (nc_ofg[0], nm.sum( timeStats['ofg'] ) )

    if conf.log:
        log( of, ofgNorm, alpha, finished = True )
        return xit, log
    else:
        return xit
