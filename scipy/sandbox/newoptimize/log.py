from base import *

##
# 20.04.2006, c
class Log( Struct ):

    ##
    # 20.04.2006, c
    # 21.04.2006
    def fromConf( conf, dataNames ):
        obj = Log( dataNames = dataNames, seqDataNames = [], igs = [],
                   data = {}, nCalls = 0 )

        for ig, names in enumerate( obj.dataNames ):
            for name in names:
                obj.data[name] = []
                obj.igs.append( ig )
                obj.seqDataNames.append( name )
        obj.nArg = len( obj.igs )
        obj.nGr = len( obj.dataNames )

        try:
            obj.isPlot = conf.isPlot
        except:
            obj.isPlot = True

        return obj
    fromConf = staticmethod( fromConf )

    ##
    # 20.04.2006, c
    # 21.04.2006
    # 26.04.2006
    def __call__( self, *args, **kwargs ):

        finished = False
        if kwargs:
            if kwargs.has_key( 'finished' ):
                finished = kwargs['finished']

        ls = len( args ), self.nArg
        if ls[0] != ls[1]:
            raise IndexError, '%d == %d' % ls

        for ii, name in enumerate( self.seqDataNames ):
            self.data[name].append( args[ii] )

        if self.isPlot:
            pylab.ion()
            if self.nCalls == 0:
                self.fig = pylab.figure()
                self.ax = []
                for ig in range( self.nGr ):
                    self.ax.append( pylab.subplot( 100 * self.nGr + 11 + ig ) )


            for ig in range( self.nGr ):
                self.ax[ig].clear()
            for ii, name in enumerate( self.seqDataNames ):
                try:
                    self.ax[self.igs[ii]].plot( self.data[name] )
                except:
                    pass
            for ig in range( self.nGr ):
                self.ax[ig].legend( self.dataNames[ig] )
            pylab.draw()
            pylab.ioff()

        self.nCalls += 1

        if finished and self.isPlot:
            pylab.show()
            return
