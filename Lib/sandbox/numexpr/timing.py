import timeit, numpy

array_size = 1e6
iterations = 10

def compare_times(setup, expr):
    print "Expression:", expr
    namespace = {}
    exec setup in namespace    
    if not numpy.alltrue(eval(expr, namespace) == eval('evaluate("%s")' % expr, namespace)):
        print eval(expr, namespace)[:10], 'versus', eval('evaluate("%s")' % expr, namespace)[:10]
        raise RuntimeError("numexpr returned incorrect value")
        
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


setup3 = """\
from numpy import arange, sin, cos, sinh
try: from scipy.weave import blitz
except: pass
from numexpr import evaluate
a = arange(%f)
b = arange(%f)
result = arange(%f)
""" % ((array_size,)*3)
expr3 = '2*a + (cos(3)+5)*sinh(cos(b))'


setup4 = """\
from numpy import arange, sin, cos, sinh, arctan2
try: from scipy.weave import blitz
except: pass
from numexpr import evaluate
a = arange(%f)
b = arange(%f)
result = arange(%f)
""" % ((array_size,)*3)
expr4 = '2*a + arctan2(a, b)'


setup5 = """\
from numpy import arange, sin, cos, sinh, arctan2, where
try: from scipy.weave import blitz
except: pass
from numexpr import evaluate
a = arange(%f)
b = arange(%f)
result = arange(%f)
""" % ((array_size,)*3)
expr5 = 'where(0.1*a > arctan2(a, b), 2*a, arctan2(a,b))'

expr6 = 'where(a, 2, b)'

expr7 = 'where(a-10, a, 2)'

expr8 = 'where(a%2, b+5, 2)'

expr9 = 'where(a%2, 2, b+5)'

expr10 = 'a**2 + b**0.5'

def compare(check_only=False):
    #~ compare_times(setup1, expr1)
    #~ print
    #~ compare_times(setup2, expr2)
    #~ print
    #~ compare_times(setup3, expr3)
    #~ print
    #~ compare_times(setup4, expr4)
    #~ print
    #~ compare_times(setup5, expr6)
    #~ print
    #~ compare_times(setup5, expr7)
    #~ print
    #~ compare_times(setup5, expr8)
    #~ print
    #~ compare_times(setup5, expr9)
    print
    compare_times(setup5, expr10)

if __name__ == '__main__':
    compare()
