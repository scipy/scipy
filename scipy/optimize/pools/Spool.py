from __future__ import division, print_function, absolute_import

__all__ = ["Spool"]

class Spool(object):
    '''
    Serial pool which wraps an objective function into a pool object.
    Created for unification of the code between serial and parallel cases.
    
    '''

    def map(self, func, task):
        return func(task)
        
    def poolsize(self):
        return 1
