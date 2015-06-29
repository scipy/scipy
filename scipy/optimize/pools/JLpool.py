from __future__ import division, print_function, absolute_import

__all__ = ["JLpool"]

# conditional import of joblib
try:
    from joblib import Parallel, delayed, cpu_count
except ImportError:
    Parallel = None

class JLpool(object):
    '''
    A joblib based pool which evaluates objective function using multiprocessing
    of the joblib.Parallel. Can be initialized with the same parameters as the
    Parallel itself. The pool only support the :func:`map` method.
    
    Example 1: joblib pool

    >>> from JLpool import JLpool
    >>> from scipy.optimize import rosen
    >>> import differential_evolution as de
    >>> bounds = [(0,2), (0, 2), (0,2)]
    >>> #if n_jobs is not specified -- sprawns cpu_count() processes by default
    >>> pool = JLpool()
    >>> 
    >>> #run DE 
    >>> result = de(rosen, bounds, pool=pool, aggressive=True, disp=True)
    '''
    
    def __init__(self, *args, **kwargs):
        if Parallel is None:
            raise ImportError("joblib")
        
        self.args = args
        self.kwargs = kwargs
        
    
        self.nworkers = None
        if 'n_jobs' in kwargs:
            n_jobs = kwargs.pop('n_jobs')
            self.nworkers = n_jobs
        else:
            self.nworkers = cpu_count()
            

    def map(self, func, tasks):
        dfunc = delayed(func)
        return Parallel(n_jobs=self.nworkers, *(self.args), **(self.kwargs))(dfunc(a) for a in tasks)
        
    def poolsize(self):
        return self.nworkers
