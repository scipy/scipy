# $Id$
#  -----------------------------------------------------------------
#  LLNL-specific file
#  -----------------------------------------------------------------

from posix import system
try:
   system ( "/usr/apps/tracker/bin/tracker -s -n PYGIST -v %s" % __version__ )
except:
   pass
