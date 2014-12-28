from __future__ import division, print_function, absolute_import

from math import sqrt

from numpy.testing import *


# Import testing parameters
from scipy.optimize._tstutils import (methods, mstrings, functions,
     fstrings, description)


class BenchZeros(TestCase):
    def bench_run(self):
        a = .5
        b = sqrt(3)
        repeat = 2000

        print(description)

        print('TESTING SPEED\n')
        print('times in seconds for %d iterations \n' % repeat)
        for i in range(len(functions)):
            print('function %s\n' % fstrings[i])
            func = functions[i]
            for j in range(len(methods)):
                meth = methods[j]
                try:
                    t = measure("meth(func,a,b)",repeat)
                except:
                    print('%s : failed' % mstrings[j])
                else:
                    print('%s : %5.3f' % (mstrings[j],t))
            print('\n\n')

if __name__ == '__main__':
    Tester().bench()
