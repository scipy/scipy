# -*- coding: iso-8859-1 -*-
import unittest
from sparray import sparray

class test_sparray(unittest.TestCase):
    def testinit(self):
        a = sparray((5,5))
        a = sparray((5,5), {(3,3):13}, (0,0))
        a = sparray((5,5), {(123,-48):13}, (120,-50))
        self.assertEqual(a.length, 25)

        #tableaux 1D
        a=sparray(5)
        a=sparray(5, {13:10}, 9)
        self.assertEqual(a.is1D,True)

    def testsetitem(self):
        a = sparray((5,5))
        a[5] = 2
        a[3,3] = 7
        self.assertEqual(a.comp((3,3)), 18)
        self.assertEqual(a.decomp(17), (3,2))
        self.assertEqual(a[3,3], a[18])
        a[:,3] = (1,2,3,4,5)
        self.assertEqual(a[3,3], 4)
        
    def testgetitem(self):
        a = sparray((5,5))
        a[:,3] = (1,2,3,4,5)
        self.assertEqual(a[:,3], map(float, (1,2,3,4,5)))
        a[3,3] = 7
        self.assertEqual(a[3,3], 7)
        self.assertEqual(a[3,:], map(float, (0,0,0,7,0)))
    
 


if __name__=='__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(test_sparray))
    unittest.TextTestRunner(verbosity=2).run(suite)

