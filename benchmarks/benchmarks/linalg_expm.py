import numpy as np

from .common import Benchmark, safe_import

with safe_import():
    from scipy.linalg import expm



class Expm_2d(Benchmark):

    def setup(self):
        self.small_2d_matrices = np.array([
            [[1, 4], [1, 1]],
            [[1, 3], [1, -1]],
            [[1, 3], [4, 5]],
            [[1, 3], [5, 3]],
            [[4, 5], [-3, -4]]
            ], order='F')
        
        self.large_2d_matrices = np.array([
            [[-494.08845191, 0], [12566.3706, - 12566.3706]],
            [[200, 1000], [3, 4]],
            [[90, 100], [500, 5]],
            [[80, 30.5], [20.17, 200.5]]
            ], order='F')
        
    
    def time_2d_small(self):
        expm(self.small_2d_matrices)
    
    def time_2d_large(self):
        expm(self.large_2d_matrices)
    

class Expm_nd(Benchmark):
    params = [[500, 1000],
              [20, 50]
              ]
    param_names = ['n', 'num_matrices']

    def setup(self, n, num_matrices):
        self.nd_matrices = np.random.rand(num_matrices, n, n)
    
    def time_nd_matrices(self, n, num_matrices):
        expm(self.nd_matrices)




    


