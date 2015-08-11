from __future__ import division, print_function, absolute_import
from multiprocessing import Pool, cpu_count

__all__ = ["PPool"]

class PPool(object):
    '''
    Parallel pool which wraps an objective function into a pool object for
    parallel execution. This pool is based on multiprocessing.Pool.map 
    augumented with `poolsize` method to tell the user the number of active 
    workers.
    
    Parameters
    ----------
    n_jobs : int (Default)
        Sets the number of workers. In order to keep the responsivenes of the
        system, sometimes, this number could be set less then cpu_count().
        The default is cpu_count()

    '''

    def __init__(self, *args, **kwargs):
        
        self.args = args
        self.kwargs = kwargs
        
        self.nworkers = None
        if 'n_jobs' in kwargs:
            n_jobs = kwargs.pop('n_jobs')
            self.nworkers = n_jobs
        else:
            self.nworkers = cpu_count()
        
        self.pool = Pool(processes=self.nworkers)   
    
    def map(self, func, task):
        return self.pool.map(func, task)
        
    def poolsize(self):
        return self.nworkers
