from __future__ import division, print_function, absolute_import

__all__ = ["SPool"]

class SPool(object):
    '''
    Serial pool which wraps an objective function into a pool object.
    Created for the unification of the code between serial and parallel cases.
    
    '''

    def map(self, func, task):
        return map(func, task)
        
    def poolsize(self):
        return 1
