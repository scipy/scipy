
from scipy_test.testing import *

set_package_path()
from optimize import cobyla as co
restore_path()

class test_cobyla(ScipyTestCase):
    def check_simple(self, level=1):

        function = lambda x: x[0]**2 + abs(x[1])**3
	con1 = lambda x: x[0]**2 + x[1]**2 - 25
	con2 = lambda x: -con1(x)

        x = co.fmin_cobyla(function, [4.95,0.66], [con1, con2], rhobeg=1, 
			   rhoend=1e-5, iprint=0, maxfun=100)
        print 'Result:',x,'(exact result = 4.955356249106168, 0.666666666666666)'

if __name__ == "__main__":
    ScipyTest().run()
