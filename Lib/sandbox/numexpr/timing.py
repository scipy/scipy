import timeit

array_size = 1e6
iterations = 10

def compare_times(setup, expr):
    print "Expression:", expr

    numpy_timer = timeit.Timer(expr, setup)
    numpy_time = numpy_timer.timeit(number=iterations)
    print 'numpy:', numpy_time / iterations

    try:
        weave_timer = timeit.Timer('blitz("result=%s")' % expr, setup)
        weave_time = weave_timer.timeit(number=iterations)
        print "Weave:", weave_time/iterations
    
        print "Speed-up of weave over numpy:", numpy_time/weave_time
    except:
        print "Skipping weave timing"

    numexpr_timer = timeit.Timer('evaluate("%s")' % expr, setup)
    numexpr_time = numexpr_timer.timeit(number=iterations)
    print "numexpr:", numexpr_time/iterations

    print "Speed-up of numexpr over numpy:", numpy_time/numexpr_time

setup1 = """\
from numpy import arange
try: from scipy.weave import blitz
except: pass
from numexpr import evaluate
result = arange(%f)
b = arange(%f)
c = arange(%f)
d = arange(%f)
e = arange(%f)
""" % ((array_size,)*5)
expr1 = 'b*c+d*e'

setup2 = """\
from numpy import arange
try: from scipy.weave import blitz
except: pass
from numexpr import evaluate
a = arange(%f)
b = arange(%f)
result = arange(%f)
""" % ((array_size,)*3)
expr2 = '2*a+3*b'



def compare():
    compare_times(setup1, expr1)
    print
    compare_times(setup2, expr2)

if __name__ == '__main__':
    compare()
