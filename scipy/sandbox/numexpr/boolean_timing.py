import sys
import timeit
import numpy

array_size = 1000*1000
iterations = 10

numpy_ttime = 0
numpy_sttime = 0
numpy_nttime = 0
numexpr_ttime = 0
numexpr_sttime = 0
numexpr_nttime = 0

def compare_times(expr, nexpr):
    global numpy_ttime
    global numpy_sttime
    global numpy_nttime
    global numexpr_ttime
    global numexpr_sttime
    global numexpr_nttime

    print "******************* Expression:", expr

    setup_contiguous = setupNP_contiguous
    setup_strided = setupNP_strided
    setup_unaligned = setupNP_unaligned

    numpy_timer = timeit.Timer(expr, setup_contiguous)
    numpy_time = round(numpy_timer.timeit(number=iterations), 4)
    numpy_ttime += numpy_time
    print 'numpy:', numpy_time / iterations

    numpy_timer = timeit.Timer(expr, setup_strided)
    numpy_stime = round(numpy_timer.timeit(number=iterations), 4)
    numpy_sttime += numpy_stime
    print 'numpy strided:', numpy_stime / iterations

    numpy_timer = timeit.Timer(expr, setup_unaligned)
    numpy_ntime = round(numpy_timer.timeit(number=iterations), 4)
    numpy_nttime += numpy_ntime
    print 'numpy unaligned:', numpy_ntime / iterations

    evalexpr = 'evaluate("%s", optimization="aggressive")' % expr
    numexpr_timer = timeit.Timer(evalexpr, setup_contiguous)
    numexpr_time = round(numexpr_timer.timeit(number=iterations), 4)
    numexpr_ttime += numexpr_time
    print "numexpr:", numexpr_time/iterations,
    print "Speed-up of numexpr over numpy:", round(numpy_time/numexpr_time, 4)

    evalexpr = 'evaluate("%s", optimization="aggressive")' % expr
    numexpr_timer = timeit.Timer(evalexpr, setup_strided)
    numexpr_stime = round(numexpr_timer.timeit(number=iterations), 4)
    numexpr_sttime += numexpr_stime
    print "numexpr strided:", numexpr_stime/iterations,
    print "Speed-up of numexpr strided over numpy:", \
          round(numpy_stime/numexpr_stime, 4)

    evalexpr = 'evaluate("%s", optimization="aggressive")' % expr
    numexpr_timer = timeit.Timer(evalexpr, setup_unaligned)
    numexpr_ntime = round(numexpr_timer.timeit(number=iterations), 4)
    numexpr_nttime += numexpr_ntime
    print "numexpr unaligned:", numexpr_ntime/iterations,
    print "Speed-up of numexpr unaligned over numpy:", \
          round(numpy_ntime/numexpr_ntime, 4)



setupNP = """\
from numpy import arange, where, arctan2, sqrt
from numpy import rec as records
from numexpr import evaluate

# Initialize a recarray of 16 MB in size
r=records.array(None, formats='a%s,i4,f8', shape=%s)
c1 = r.field('f0')%s
i2 = r.field('f1')%s
f3 = r.field('f2')%s
c1[:] = "a"
i2[:] = arange(%s)/1000
f3[:] = i2/2.
"""

setupNP_contiguous = setupNP % (4, array_size,
                                ".copy()", ".copy()", ".copy()",
                                array_size)
setupNP_strided = setupNP % (4, array_size, "", "", "", array_size)
setupNP_unaligned = setupNP % (1, array_size, "", "", "", array_size)


expressions = []
expressions.append('i2 > 0')
expressions.append('i2 < 0')
expressions.append('i2 < f3')
expressions.append('i2-10 < f3')
expressions.append('i2*f3+f3*f3 > i2')
expressions.append('0.1*i2 > arctan2(i2, f3)')
expressions.append('i2%2 > 3')
expressions.append('i2%10 < 4')
expressions.append('i2**2 + (f3+1)**-2.5 < 3')
expressions.append('(f3+1)**50 > i2')
expressions.append('sqrt(i2**2 + f3**2) > 1')
expressions.append('(i2>2) | ((f3**2>3) & ~(i2*f3<2))')

def compare(expression=False):
    if expression:
        compare_times(expression, 1)
        sys.exit(0)
    nexpr = 0
    for expr in expressions:
        nexpr += 1
        compare_times(expr, nexpr)
    print

if __name__ == '__main__':

    print 'Python version:    %s' % sys.version
    print "NumPy version:     %s" % numpy.__version__

    if len(sys.argv) > 1:
        expression = sys.argv[1]
        print "expression-->", expression
        compare(expression)
    else:
        compare()

    print "*************** TOTALS **************************"
    print "numpy total:", numpy_ttime/iterations
    print "numpy strided total:", numpy_sttime/iterations
    print "numpy unaligned total:", numpy_nttime/iterations
    print "numexpr total:", numexpr_ttime/iterations
    print "Speed-up of numexpr over numpy:", \
          round(numpy_ttime/numexpr_ttime, 3)
    print "numexpr strided total:", numexpr_sttime/iterations
    print "Speed-up of numexpr strided over numpy:", \
          round(numpy_sttime/numexpr_sttime, 3)
    print "numexpr unaligned total:", numexpr_nttime/iterations
    print "Speed-up of numexpr unaligned over numpy:", \
          round(numpy_nttime/numexpr_nttime, 3)
