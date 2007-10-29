##
# 16.02.2005, c
import numpy as nm
import scipy.optimize as sopt
import scipy.optimize.linesearch as linesearch
import pylab

import glob, re, time, sys

##
# 22.09.2005, c
# 24.10.2005
if sys.version[:5] < '2.4.0':
    def sorted( sequence ):
        tmp = copy( sequence )
        tmp.sort()
        return tmp


# Some usefull definitions.

##
# 02.01.2005
class Struct( object ):
    # 03.10.2005, c
    # 26.10.2005
    def __init__( self, **kwargs ):
        if kwargs:
            self.__dict__.update( kwargs )

    # 08.03.2005
    def __str__( self ):
        ss = "%s\n" % self.__class__
        for key, val in self.__dict__.iteritems():
            if (issubclass( self.__dict__[key].__class__, Struct )):
                ss += "  %s:\n    %s\n" % (key, self.__dict__[key].__class__)
            else:
                aux = "\n" + str( val )
                aux = aux.replace( "\n", "\n    " );
                ss += "  %s:\n%s\n" % (key, aux[1:])
        return( ss.rstrip() )

    # 08.03.2005, c
    def strAll( self ):
        ss = "%s\n" % self.__class__
        for key, val in self.__dict__.iteritems():
            if (issubclass( self.__dict__[key].__class__, Struct )):
                ss += "  %s:\n" % key
                aux = "\n" + self.__dict__[key].strAll()
                aux = aux.replace( "\n", "\n    " );
                ss += aux[1:] + "\n"
            else:
                aux = "\n" + str( val )
                aux = aux.replace( "\n", "\n    " );
                ss += "  %s:\n%s\n" % (key, aux[1:])
        return( ss.rstrip() )
