#
# Title: Using chaco.wxplt from Python prompt under Linux.
# Author: Pearu Peterson <pearu@cens.ioc.ee>
# Created: October 2003
#
# Usage: from Python prompt import this module as
#
#  execfile('/path/to/chaco_wxplt.py')
#
# The obvious way `from chaco_wxplt import *` will hang
# when importing wxPython.wx for reasons that are mystery for me.
#
# For example:
#   >>> import gui_thread
#   >>> execfile(gui_thread.__path__[0]+'/chaco_wxplt.py')
#   >>> plot([1,2])
#

import sys
import inspect
import threading
from scipy_base import migrate, ParallelExec

# Create a wxPython thread. Before check for no earlier wxPython
# imports, otherwise they would be in a wrong thread:
assert not sys.modules.has_key('wxPython.wx'),\
       'wxPython is already imported (probably from a wrong thread)'
pexec=ParallelExec()

# Wait until bg_app appears in the current namespace:
pexec('from gui_thread.wximporter import *;bg_app=wxBackgroundApp()',wait=1)

# This will run in background forever... until bg_app.shutdown()
# is called by atexit module:
pexec('bg_app.run()')

# Chaco must be imported from wxPython thread because it already
# executes some wxPyhton functions.
bg_app('import chaco.wxplt')

# Now create wrappers to chaco functions so that they will
# be executed in the wxPython thread:
for _n,_v in inspect.getmembers(chaco.wxplt,inspect.isroutine):
    if sys.version[:3]>'2.2' and _v.__module__[:5]!='chaco' \
       or _n in ['bin_search1','bin_search2','col','hstack',
                 'row','v','remove_bad_vals','split_name_for',
                 'process_format','amax','amin','user_name_for',
                 'is1D']:
        # Avoid calling non chaco functions (that are in chaco
        # namespace) from wxPython thread.
        # What chaco.wxplt really needs is a __all__ list.
        continue
    exec '%s = migrate(_v,bg_app)' % (_n)
del _n,_v

# Warning: You must use only the above defined wrappers!
# Using anything else from chaco or wxPython will cause
# `Xlib: unexpected async reply` crashes.
# So, remove chaco from this namespace to avoid any misuses:
del chaco

# BUT if you really need to use other chaco or wxPython
# function then you can call them through bg_app. For example:
#   some_result = bg_app('chaco.some_function(..)')
# Using bg_app will ensure that the command will be exexuted
# in wxPython thread.
#
# You might need (read: I am not sure; if there will be a crash then
# you'll certainly need) to use bg_app also in the following cases:
#
#   >>> f=figure()
#   >>> f_title = bg_app('f.GetTitle()')
#
# for instance.
