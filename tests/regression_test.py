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
sys.path.append('c:/home/python/Interpolate1d/tests')

import shelve, time
from test_interpolate1d import Test as Test1
from test_interpolate2d import Test as Test2
from test_interpolateNd import Test as TestN

# name of log file to which all data is stored.
filename = 'regression_test.dbm'

log_total = shelve.open(filename)
current_time = str(time.localtime()[0:5]) # specified up to the minute
current_dict = {} # holds results for each dimensionality

# run all tests in interpolate1d's test class
if True:
    test_list = [name for name in dir(Test1) if name.find('test_') == 0]
    Test1_dict = {}

    # record time taken for each test
    for test_name in test_list:
        t1 = time.clock()
        eval('Test1.%s' % test_name)
        t2 = time.clock()
        Test1_dict[test_name] = t2-t1

    current_dict['Test1'] = Test1_dict

# 2-dimensional
if True:
    test_list = [name for name in dir(Test2) if name.find('test_') == 0]
    Test2_dict = {}

    # record time taken for each test
    for test_name in test_list:
        t1 = time.clock()
        eval('Test2.%s' % test_name)
        t2 = time.clock()
        Test2_dict[test_name] = t2-t1

    current_dict['Test2'] = Test2_dict

# N-dmensional
if True:
    test_list = [name for name in dir(TestN) if name.find('test_') == 0]
    TestN_dict = {}

    # record time taken for each test
    for test_name in test_list:
        t1 = time.clock()
        eval('TestN.%s' % test_name)
        t2 = time.clock()
        TestN_dict[test_name] = t2-t1

    current_dict['TestN'] = TestN_dict

log_total[current_time] = current_dict
log_total.close()
