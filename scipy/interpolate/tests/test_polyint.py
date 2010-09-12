
from numpy.testing import assert_almost_equal, assert_array_equal, \
        TestCase, run_module_suite
from scipy.interpolate import KroghInterpolator, krogh_interpolate, \
        BarycentricInterpolator, barycentric_interpolate, \
        PiecewisePolynomial, piecewise_polynomial_interpolate, \
        approximate_taylor_polynomial
import scipy
import numpy as np
from scipy.interpolate import splrep, splev

class CheckKrogh(TestCase):
    def setUp(self):
        self.true_poly = scipy.poly1d([-2,3,1,5,-4])
        self.test_xs = np.linspace(-1,1,100)
        self.xs = np.linspace(-1,1,5)
        self.ys = self.true_poly(self.xs)

    def test_lagrange(self):
        P = KroghInterpolator(self.xs,self.ys)
        assert_almost_equal(self.true_poly(self.test_xs),P(self.test_xs))

    def test_scalar(self):
        P = KroghInterpolator(self.xs,self.ys)
        assert_almost_equal(self.true_poly(7),P(7))
        assert_almost_equal(self.true_poly(np.array(7)), P(np.array(7)))

    def test_derivatives(self):
        P = KroghInterpolator(self.xs,self.ys)
        D = P.derivatives(self.test_xs)
        for i in xrange(D.shape[0]):
            assert_almost_equal(self.true_poly.deriv(i)(self.test_xs),
                                D[i])

    def test_low_derivatives(self):
        P = KroghInterpolator(self.xs,self.ys)
        D = P.derivatives(self.test_xs,len(self.xs)+2)
        for i in xrange(D.shape[0]):
            assert_almost_equal(self.true_poly.deriv(i)(self.test_xs),
                                D[i])

    def test_derivative(self):
        P = KroghInterpolator(self.xs,self.ys)
        m = 10
        r = P.derivatives(self.test_xs,m)
        for i in xrange(m):
            assert_almost_equal(P.derivative(self.test_xs,i),r[i])

    def test_high_derivative(self):
        P = KroghInterpolator(self.xs,self.ys)
        for i in xrange(len(self.xs),2*len(self.xs)):
            assert_almost_equal(P.derivative(self.test_xs,i),
                                np.zeros(len(self.test_xs)))

    def test_hermite(self):
        xs = [0,0,0,1,1,1,2]
        ys = [self.true_poly(0),
              self.true_poly.deriv(1)(0),
              self.true_poly.deriv(2)(0),
              self.true_poly(1),
              self.true_poly.deriv(1)(1),
              self.true_poly.deriv(2)(1),
              self.true_poly(2)]
        P = KroghInterpolator(self.xs,self.ys)
        assert_almost_equal(self.true_poly(self.test_xs),P(self.test_xs))

    def test_vector(self):
        xs = [0, 1, 2]
        ys = np.array([[0,1],[1,0],[2,1]])
        P = KroghInterpolator(xs,ys)
        Pi = [KroghInterpolator(xs,ys[:,i]) for i in xrange(ys.shape[1])]
        test_xs = np.linspace(-1,3,100)
        assert_almost_equal(P(test_xs),
                np.rollaxis(np.asarray([p(test_xs) for p in Pi]),-1))
        assert_almost_equal(P.derivatives(test_xs),
                np.transpose(np.asarray([p.derivatives(test_xs) for p in Pi]),
                    (1,2,0)))

    def test_empty(self):
        P = KroghInterpolator(self.xs,self.ys)
        assert_array_equal(P([]), [])

    def test_shapes_scalarvalue(self):
        P = KroghInterpolator(self.xs,self.ys)
        assert_array_equal(np.shape(P(0)), ())
        assert_array_equal(np.shape(P(np.array(0))), ())
        assert_array_equal(np.shape(P([0])), (1,))
        assert_array_equal(np.shape(P([0,1])), (2,))

    def test_shapes_scalarvalue_derivative(self):
        P = KroghInterpolator(self.xs,self.ys)
        n = P.n
        assert_array_equal(np.shape(P.derivatives(0)), (n,))
        assert_array_equal(np.shape(P.derivatives(np.array(0))), (n,))
        assert_array_equal(np.shape(P.derivatives([0])), (n,1))
        assert_array_equal(np.shape(P.derivatives([0,1])), (n,2))

    def test_shapes_vectorvalue(self):
        P = KroghInterpolator(self.xs,np.outer(self.ys,np.arange(3)))
        assert_array_equal(np.shape(P(0)), (3,))
        assert_array_equal(np.shape(P([0])), (1,3))
        assert_array_equal(np.shape(P([0,1])), (2,3))

    def test_shapes_1d_vectorvalue(self):
        P = KroghInterpolator(self.xs,np.outer(self.ys,[1]))
        assert_array_equal(np.shape(P(0)), (1,))
        assert_array_equal(np.shape(P([0])), (1,1))
        assert_array_equal(np.shape(P([0,1])), (2,1))

    def test_shapes_vectorvalue_derivative(self):
        P = KroghInterpolator(self.xs,np.outer(self.ys,np.arange(3)))
        n = P.n
        assert_array_equal(np.shape(P.derivatives(0)), (n,3))
        assert_array_equal(np.shape(P.derivatives([0])), (n,1,3))
        assert_array_equal(np.shape(P.derivatives([0,1])), (n,2,3))

    def test_wrapper(self):
        P = KroghInterpolator(self.xs,self.ys)
        assert_almost_equal(P(self.test_xs),krogh_interpolate(self.xs,self.ys,self.test_xs))
        assert_almost_equal(P.derivative(self.test_xs,2),krogh_interpolate(self.xs,self.ys,self.test_xs,der=2))
        assert_almost_equal(P.derivatives(self.test_xs,2),krogh_interpolate(self.xs,self.ys,self.test_xs,der=[0,1]))

class CheckTaylor(TestCase):
    def test_exponential(self):
        degree = 5
        p = approximate_taylor_polynomial(np.exp, 0, degree, 1, 15)
        for i in xrange(degree+1):
            assert_almost_equal(p(0),1)
            p = p.deriv()
        assert_almost_equal(p(0),0)

class CheckBarycentric(TestCase):
    def setUp(self):
        self.true_poly = scipy.poly1d([-2,3,1,5,-4])
        self.test_xs = np.linspace(-1,1,100)
        self.xs = np.linspace(-1,1,5)
        self.ys = self.true_poly(self.xs)

    def test_lagrange(self):
        P = BarycentricInterpolator(self.xs,self.ys)
        assert_almost_equal(self.true_poly(self.test_xs),P(self.test_xs))

    def test_scalar(self):
        P = BarycentricInterpolator(self.xs,self.ys)
        assert_almost_equal(self.true_poly(7),P(7))
        assert_almost_equal(self.true_poly(np.array(7)),P(np.array(7)))

    def test_delayed(self):
        P = BarycentricInterpolator(self.xs)
        P.set_yi(self.ys)
        assert_almost_equal(self.true_poly(self.test_xs),P(self.test_xs))

    def test_append(self):
        P = BarycentricInterpolator(self.xs[:3],self.ys[:3])
        P.add_xi(self.xs[3:],self.ys[3:])
        assert_almost_equal(self.true_poly(self.test_xs),P(self.test_xs))

    def test_vector(self):
        xs = [0, 1, 2]
        ys = np.array([[0,1],[1,0],[2,1]])
        P = BarycentricInterpolator(xs,ys)
        Pi = [BarycentricInterpolator(xs,ys[:,i]) for i in xrange(ys.shape[1])]
        test_xs = np.linspace(-1,3,100)
        assert_almost_equal(P(test_xs),
                np.rollaxis(np.asarray([p(test_xs) for p in Pi]),-1))

    def test_shapes_scalarvalue(self):
        P = BarycentricInterpolator(self.xs,self.ys)
        assert_array_equal(np.shape(P(0)), ())
        assert_array_equal(np.shape(P(np.array(0))), ())
        assert_array_equal(np.shape(P([0])), (1,))
        assert_array_equal(np.shape(P([0,1])), (2,))

    def test_shapes_vectorvalue(self):
        P = BarycentricInterpolator(self.xs,np.outer(self.ys,np.arange(3)))
        assert_array_equal(np.shape(P(0)), (3,))
        assert_array_equal(np.shape(P([0])), (1,3))
        assert_array_equal(np.shape(P([0,1])), (2,3))

    def test_shapes_1d_vectorvalue(self):
        P = BarycentricInterpolator(self.xs,np.outer(self.ys,[1]))
        assert_array_equal(np.shape(P(0)), (1,))
        assert_array_equal(np.shape(P([0])), (1,1))
        assert_array_equal(np.shape(P([0,1])), (2,1))

    def test_wrapper(self):
        P = BarycentricInterpolator(self.xs,self.ys)
        assert_almost_equal(P(self.test_xs),barycentric_interpolate(self.xs,self.ys,self.test_xs))

class CheckPiecewise(TestCase):
    def setUp(self):
        self.tck = splrep([0,1,2,3,4,5], [0,10,-1,3,7,2], s=0)
        self.test_xs = np.linspace(-1,6,100)
        self.spline_ys = splev(self.test_xs, self.tck)
        self.spline_yps = splev(self.test_xs, self.tck, der=1)
        self.xi = np.unique(self.tck[0])
        self.yi = [[splev(x, self.tck, der=j) for j in xrange(3)] for x in self.xi]

    def test_construction(self):
        P = PiecewisePolynomial(self.xi, self.yi, 3)
        assert_almost_equal(P(self.test_xs), self.spline_ys)

    def test_scalar(self):
        P = PiecewisePolynomial(self.xi,self.yi,3)
        assert_almost_equal(P(self.test_xs[0]),self.spline_ys[0])
        assert_almost_equal(P.derivative(self.test_xs[0],1),self.spline_yps[0])
        assert_almost_equal(P(np.array(self.test_xs[0])),self.spline_ys[0])
        assert_almost_equal(P.derivative(np.array(self.test_xs[0]),1),
                            self.spline_yps[0])

    def test_derivative(self):
        P = PiecewisePolynomial(self.xi,self.yi,3)
        assert_almost_equal(P.derivative(self.test_xs,1),self.spline_yps)

    def test_derivatives(self):
        P = PiecewisePolynomial(self.xi,self.yi,3)
        m = 4
        r = P.derivatives(self.test_xs,m)
        #print r.shape, r
        for i in xrange(m):
            assert_almost_equal(P.derivative(self.test_xs,i),r[i])

    def test_vector(self):
        xs = [0, 1, 2]
        ys = [[[0,1]],[[1,0],[-1,-1]],[[2,1]]]
        P = PiecewisePolynomial(xs,ys)
        Pi = [PiecewisePolynomial(xs,[[yd[i] for yd in y] for y in ys])
            for i in xrange(len(ys[0][0]))]
        test_xs = np.linspace(-1,3,100)
        assert_almost_equal(P(test_xs),
                np.rollaxis(np.asarray([p(test_xs) for p in Pi]),-1))
        assert_almost_equal(P.derivative(test_xs,1),
                np.transpose(np.asarray([p.derivative(test_xs,1) for p in Pi]),
                    (1,0)))

    def test_incremental(self):
        P = PiecewisePolynomial([self.xi[0]], [self.yi[0]], 3)
        for i in xrange(1,len(self.xi)):
            P.append(self.xi[i],self.yi[i],3)
        assert_almost_equal(P(self.test_xs),self.spline_ys)

    def test_shapes_scalarvalue(self):
        P = PiecewisePolynomial(self.xi,self.yi,4)
        assert_array_equal(np.shape(P(0)), ())
        assert_array_equal(np.shape(P(np.array(0))), ())
        assert_array_equal(np.shape(P([0])), (1,))
        assert_array_equal(np.shape(P([0,1])), (2,))

    def test_shapes_scalarvalue_derivative(self):
        P = PiecewisePolynomial(self.xi,self.yi,4)
        n = 4
        assert_array_equal(np.shape(P.derivative(0,1)), ())
        assert_array_equal(np.shape(P.derivative(np.array(0),1)), ())
        assert_array_equal(np.shape(P.derivative([0],1)), (1,))
        assert_array_equal(np.shape(P.derivative([0,1],1)), (2,))

    def test_shapes_vectorvalue(self):
        yi = np.multiply.outer(np.asarray(self.yi),np.arange(3))
        P = PiecewisePolynomial(self.xi,yi,4)
        assert_array_equal(np.shape(P(0)), (3,))
        assert_array_equal(np.shape(P([0])), (1,3))
        assert_array_equal(np.shape(P([0,1])), (2,3))

    def test_shapes_vectorvalue_1d(self):
        yi = np.multiply.outer(np.asarray(self.yi),np.arange(1))
        P = PiecewisePolynomial(self.xi,yi,4)
        assert_array_equal(np.shape(P(0)), (1,))
        assert_array_equal(np.shape(P([0])), (1,1))
        assert_array_equal(np.shape(P([0,1])), (2,1))

    def test_shapes_vectorvalue_derivative(self):
        P = PiecewisePolynomial(self.xi,np.multiply.outer(self.yi,np.arange(3)),4)
        n = 4
        assert_array_equal(np.shape(P.derivative(0,1)), (3,))
        assert_array_equal(np.shape(P.derivative([0],1)), (1,3))
        assert_array_equal(np.shape(P.derivative([0,1],1)), (2,3))

    def test_wrapper(self):
        P = PiecewisePolynomial(self.xi,self.yi)
        assert_almost_equal(P(self.test_xs),piecewise_polynomial_interpolate(self.xi,self.yi,self.test_xs))
        assert_almost_equal(P.derivative(self.test_xs,2),piecewise_polynomial_interpolate(self.xi,self.yi,self.test_xs,der=2))
        assert_almost_equal(P.derivatives(self.test_xs,2),piecewise_polynomial_interpolate(self.xi,self.yi,self.test_xs,der=[0,1]))


if __name__=='__main__':
    run_module_suite()
