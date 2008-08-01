""" 
    regression test:

    This script runs a simple regression test on the functionality of
    the interpolation module.  Currently, when run, it times each
    unit test in interpolate1d.py and stores those times in a dict
    of dicts; outer keys are time test was performed, and inner
    keys are names of tests run.

"""

# hack to test on Field's computer
import sys
sys.path.append('c:/home/python/Interpolate1d')

import shelve, time
from test_interpolate1d import Test

# name of log file to which all data is stored.
filename = 'regression_test.dbm'

log_total = shelve.open(filename)
current_time = str(time.localtime()[0:5]) # specified up to the minute

# run all tests in interpolate1d's test class
test_list = [name for name in dir(Test) if name.find('test_') == 0]
log_now = {}

# record time taken for each test
for test_name in test_list:
    t1 = time.clock()
    eval('Test.%s' % test_name)
    t2 = time.clock()
    log_now[test_name] = t2-t1

log_total[current_time] = log_now
log_total.close()
