import sys
from numpy.testing import *
from scipy.spatial import cKDTree, KDTree
import numpy as np

class TestBuild(TestCase):

    def bench_build(self):
        print
        print '        Constructing kd-tree'
        print '====================================='
        print ' dim | # points |  KDTree  | cKDTree '

        for (m, n, repeat) in [(3,10000,3), (8,10000,3), (16,10000,3)]:
            print '%4s | %7s ' % (m, n),
            sys.stdout.flush()
            
            data = np.concatenate((np.random.randn(n//2,m),
                                   np.random.randn(n-n//2,m)+np.ones(m)))
            
            print '| %6.3fs ' % (measure('T1 = KDTree(data)', repeat) / repeat),
            sys.stdout.flush()
            print '| %6.3fs' % (measure('T2 = cKDTree(data)', repeat) / repeat),
            sys.stdout.flush()
            print ''

class TestQuery(TestCase):

    def bench_query(self):
        print
        print '                       Querying kd-tree'
        print '==============================================================='
        print ' dim | # points | # queries |  KDTree  | cKDTree | flat cKDTree'

        for (m, n, r, repeat) in [(3,10000,1000,3),
                                  (8,10000,1000,3),
                                  (16,10000,1000,3)]:
            print '%4s | %8s | %8s ' % (m, n, r),
            sys.stdout.flush()
            
            data = np.concatenate((np.random.randn(n//2,m),
                                   np.random.randn(n-n//2,m)+np.ones(m)))
            queries = np.concatenate((np.random.randn(r//2,m),
                                      np.random.randn(r-r//2,m)+np.ones(m)))

            T1 = KDTree(data)
            T2 = cKDTree(data)
            T3 = cKDTree(data,leafsize=n)
            print '| %6.3fs ' % (measure('T1.query(queries)', 1) / 1),
            sys.stdout.flush()
            print '| %6.3fs' % (measure('T2.query(queries)', repeat) / repeat),
            sys.stdout.flush()
            print '| %6.3fs' % (measure('T3.query(queries)', repeat) / repeat),
            sys.stdout.flush()
            print ''

class TestQueryBallPoint(TestCase):
    def bench_query_ball_point(self):
        print
        print '                   Query ball point kd-tree'
        print '==============================================================='
        print ' dim | # points | # queries | probe radius |  KDTree  | cKDTree | flat cKDTree'

        for (m, n, r, repeat) in [(3,10000,1000,3)]:#,
#                                  (8,10000,1000,3),
#                                  (16,10000,1000,3)]:
            for probe_radius in (0.2, 0.5):
                print '%4s | %8s | %9s | %11.1f ' % (m, n, r, probe_radius),
                sys.stdout.flush()
            
                data = np.concatenate((np.random.randn(n//2,m),
                                       np.random.randn(n-n//2,m)+np.ones(m)))
                queries = np.concatenate((np.random.randn(r//2,m),
                                          np.random.randn(r-r//2,m)+np.ones(m)))
                
                T1 = KDTree(data)
                T2 = cKDTree(data)
                T3 = cKDTree(data,leafsize=n)
                print '| %6.3fs ' % (measure('T1.query_ball_point(queries, probe_radius)', 1) / 1),
                sys.stdout.flush()
                print '| %6.3fs' % (measure('T2.query_ball_point(queries, probe_radius)', repeat) / repeat),
                sys.stdout.flush()
                print '| %6.3fs' % (measure('T3.query_ball_point(queries, probe_radius)', repeat) / repeat),
                sys.stdout.flush()
                print ''

if __name__ == "__main__":
    run_module_suite()
