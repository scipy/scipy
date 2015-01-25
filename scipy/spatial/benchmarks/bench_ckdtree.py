from __future__ import division, print_function, absolute_import

import sys
from numpy.testing import *
from scipy.spatial import cKDTree, KDTree
import numpy as np


class TestBuild(TestCase):

    def bench_build(self):
        print()
        print('        Constructing kd-tree')
        print('=====================================')
        print(' dim | # points |  KDTree  | cKDTree ')

        for (m, n, repeat) in [(3,10000,3), (8,10000,3), (16,10000,3)]:
            print('%4s | %7s ' % (m, n), end=' ')
            sys.stdout.flush()

            data = np.concatenate((np.random.randn(n//2,m),
                                   np.random.randn(n-n//2,m)+np.ones(m)))

            print('| %6.3fs ' % (measure('T1 = KDTree(data)', repeat) / repeat), end=' ')
            sys.stdout.flush()
            print('| %6.3fs' % (measure('T2 = cKDTree(data)', repeat) / repeat), end=' ')
            sys.stdout.flush()
            print('')


class TestQuery(TestCase):

    def bench_query(self):
        print()
        print('                       Querying kd-tree')
        print('===============================================================')
        print(' dim | # points | # queries |  KDTree  | cKDTree | flat cKDTree')

        for (m, n, r, repeat) in [(3,10000,1000,3),
                                  (8,10000,1000,3),
                                  (16,10000,1000,3)]:
            print('%4s | %8s | %8s ' % (m, n, r), end=' ')
            sys.stdout.flush()

            data = np.concatenate((np.random.randn(n//2,m),
                                   np.random.randn(n-n//2,m)+np.ones(m)))
            queries = np.concatenate((np.random.randn(r//2,m),
                                      np.random.randn(r-r//2,m)+np.ones(m)))

            T1 = KDTree(data)
            T2 = cKDTree(data)
            T3 = cKDTree(data,leafsize=n)
            print('| %6.3fs ' % (measure('T1.query(queries)', 1) / 1), end=' ')
            sys.stdout.flush()
            print('| %6.3fs' % (measure('T2.query(queries)', repeat) / repeat), end=' ')
            sys.stdout.flush()
            print('| %6.3fs' % (measure('T3.query(queries)', repeat) / repeat), end=' ')
            sys.stdout.flush()
            print('')


class TestQueryBallPoint(TestCase):
    def bench_query_ball_point(self):
        print()
        print('                   Query ball point kd-tree')
        print('===============================================================')
        print(' dim | # points | # queries | probe radius |  KDTree  | cKDTree | flat cKDTree')

        for (m, n, r, repeat) in [(3,10000,1000,3)]:
            for probe_radius in (0.2, 0.5):
                print('%4s | %8s | %9s | %11.1f ' % (m, n, r, probe_radius), end=' ')
                sys.stdout.flush()

                data = np.concatenate((np.random.randn(n//2,m),
                                       np.random.randn(n-n//2,m)+np.ones(m)))
                queries = np.concatenate((np.random.randn(r//2,m),
                                          np.random.randn(r-r//2,m)+np.ones(m)))

                T1 = KDTree(data)
                T2 = cKDTree(data)
                T3 = cKDTree(data,leafsize=n)
                print('| %6.3fs ' % (measure('T1.query_ball_point(queries, probe_radius)', 1) / 1), end=' ')
                sys.stdout.flush()
                print('| %6.3fs' % (measure('T2.query_ball_point(queries, probe_radius)', repeat) / repeat), end=' ')
                sys.stdout.flush()
                print('| %6.3fs' % (measure('T3.query_ball_point(queries, probe_radius)', repeat) / repeat), end=' ')
                sys.stdout.flush()
                print('')


class TestQueryPairs(TestCase):
    def bench_query_pairs(self):
        print()
        print('                     Query pairs kd-tree')
        print('==================================================================')
        print(' dim | # points | probe radius |  KDTree  | cKDTree | flat cKDTree')

        for (m, n, repeat) in [(3,1000,30),
                               (8,1000,30),
                               (16,1000,30)]:
            for probe_radius in (0.2, 0.5):
                print('%4s | %8s | %11.1f ' % (m, n, probe_radius), end=' ')
                sys.stdout.flush()

                data = np.concatenate((np.random.randn(n//2,m),
                                       np.random.randn(n-n//2,m)+np.ones(m)))

                T1 = KDTree(data)
                T2 = cKDTree(data)
                T3 = cKDTree(data,leafsize=n)
                print('| %6.3fs ' % (measure('T1.query_pairs(probe_radius)', 1) / 1), end=' ')
                sys.stdout.flush()
                print('| %6.3fs' % (measure('T2.query_pairs(probe_radius)', repeat) / repeat), end=' ')
                sys.stdout.flush()
                print('| %6.3fs' % (measure('T3.query_pairs(probe_radius)', repeat) / repeat), end=' ')
                sys.stdout.flush()
                print('')


class TestSparseDistanceMatrix(TestCase):
    def bench_sparse_distance_matrix(self):
        print()
        print('                   Sparse distance matrix kd-tree')
        print('====================================================================')
        print(' dim | # points T1 | # points T2 | probe radius |  KDTree  | cKDTree')

        for (m, n1, n2, repeat) in [(3,1000,1000,30),
                                    (8,1000,1000,30),
                                    (16,1000,1000,30)]:

            data1 = np.concatenate((np.random.randn(n1//2,m),
                                    np.random.randn(n1-n1//2,m)+np.ones(m)))
            data2 = np.concatenate((np.random.randn(n2//2,m),
                                    np.random.randn(n2-n2//2,m)+np.ones(m)))

            T1 = KDTree(data1)
            T2 = KDTree(data2)
            cT1 = cKDTree(data1)
            cT2 = cKDTree(data2)

            for probe_radius in (0.2, 0.5):
                print('%4s | %11s | %11s | %11.1f ' % (m, n1, n2, probe_radius), end=' ')
                sys.stdout.flush()

                print('| %6.3fs ' % (measure('T1.sparse_distance_matrix(T2, probe_radius)', 1) / 1), end=' ')
                sys.stdout.flush()
                print('| %6.3fs ' % (measure('cT1.sparse_distance_matrix(cT2, probe_radius)', repeat) / repeat), end=' ')
                sys.stdout.flush()
                print('')


class TestCountNeighbors(TestCase):
    def bench_count_neighbors(self):
        print()
        print('                     Count neighbors kd-tree')
        print('====================================================================')
        print(' dim | # points T1 | # points T2 | probe radius |  KDTree  | cKDTree')

        for (m, n1, n2, repeat) in [(3,1000,1000,30),
                                    (8,1000,1000,30),
                                    (16,1000,1000,30)]:

            data1 = np.concatenate((np.random.randn(n1//2,m),
                                    np.random.randn(n1-n1//2,m)+np.ones(m)))
            data2 = np.concatenate((np.random.randn(n2//2,m),
                                    np.random.randn(n2-n2//2,m)+np.ones(m)))

            T1 = KDTree(data1)
            T2 = KDTree(data2)
            cT1 = cKDTree(data1)
            cT2 = cKDTree(data2)

            for probe_radius in (0.2, 0.5):
                print('%4s | %11s | %11s | %11.1f ' % (m, n1, n2, probe_radius), end=' ')
                sys.stdout.flush()

                print('| %6.3fs ' % (measure('T1.count_neighbors(T2, probe_radius)', 1) / 1), end=' ')
                sys.stdout.flush()
                print('| %6.3fs ' % (measure('cT1.count_neighbors(cT2, probe_radius)', repeat) / repeat), end=' ')
                sys.stdout.flush()
                print('')

if __name__ == "__main__":
    Tester().bench()
