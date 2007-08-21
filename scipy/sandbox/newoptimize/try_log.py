from base import Struct
from log import Log

conf = Struct( isPlot = True )
log = Log.fromConf( conf, (['data1', 'data2'], ['data3']) )

log( 1, 2, 3 )
log( 1, 2, 3 )
log( 3, 2, 1 )
log( 0, 4, 1, finished = True )
