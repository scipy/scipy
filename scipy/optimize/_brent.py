from __future__ import division

from math import sqrt

from numpy import finfo

_sqrt_epsilon = sqrt(finfo(float).eps)
# about 0.3819660112501051
_golden = 0.5 * (3.0 - sqrt(5.0))


def find_minimum(f, a, b, rtol=_sqrt_epsilon, atol=_sqrt_epsilon, maxiter=500):
    """Seeks a local minimum of a function f in a closed interval [a, b] via
    Brent's method.

    Given a function f with a minimum in the interval the a <= b,
    seeks a local minima using a combination of golden section search and
    successive parabolic interpolation.

    Let ```tol = rtol * abs(x0) + atol```, where ```x0``` is the best guess
    found so far. It converges if evaluating a next guess would imply
    evaluating ```f``` at a point that is closer than ```tol``` to a previously
    evaluated one or if the number of iterations reaches ```maxiter```.

    Args:
        f (object): Objective function to be minimized.
        a, b (float): endpoints of the interval a <= b.
        rtol (float): relative tolerance. Defaults to 1.4901161193847656e-08.
        atol (float): absolute tolerance. Defaults to 1.4901161193847656e-08.
        maxiter (int): maximum number of iterations


    Returns:
        float: best guess for the minimum of f.
        float: value of f evaluated at the best guess.
        int: number of iterations performed.

    References:
        - http://people.sc.fsu.edu/~jburkardt/c_src/brent/brent.c
        - Numerical Recipes 3rd Edition: The Art of Scientific Computing
        - https://en.wikipedia.org/wiki/Brent%27s_method
    """
    # a, b: interval within the minimum should lie
    #       no function evaluation will be requested
    #       outside that range.
    # x0: least function value found so far (or the most recent one in
    #                                            case of a tie)
    # x1: second least function value
    # x2: previous value of x1
    # (x0, x1, x2): Memory triple, updated at the end of each interation.
    # u : point at which the function was evaluated most recently.
    # m : midpoint between the current interval (a, b).
    # d : step size and direction.
    # e : memorizes the step size (and direction) taken two iterations ago
    #      and it is used to (definitively) fall-back to golden-section steps
    #      when its value is too small (indicating that the polynomial fitting
    #      is not helping to speedup the convergence.)
    #
    #
    # References: Numerical Recipes: The Art of Scientific Computing
    # http://people.sc.fsu.edu/~jburkardt/c_src/brent/brent.c

    assert a <= b
    x0 = a + _golden * (b - a)
    x1 = x0
    x2 = x1
    niters = -1
    d = 0.0
    e = 0.0
    fx0 = f(x0)
    fx1 = fx0
    fx2 = fx1

    for niters in range(maxiter):

        m = 0.5 * (a + b)
        tol = rtol * abs(x0) + atol
        tol2 = 2.0 * tol

        # Check the stopping criterion.
        if abs(x0 - m) <= tol2 - 0.5 * (b - a):
            break

        r = 0.0
        q = r
        p = q

        # "To be acceptable, the parabolic step must (i) fall within the
        # bounding interval (a, b), and (ii) imply a movement from the best
        # current value x0 that is less than half the movement of the step
        # before last."
        #   - Numerical Recipes 3rd Edition: The Art of Scientific Computing.

        if tol < abs(e):
            # Compute the polynomial of the least degree (Lagrange polynomial)
            # that goes through (x0, fx0), (x1, fx1), (x2, fx2).
            r = (x0 - x1) * (fx0 - fx2)
            q = (x0 - x2) * (fx0 - fx1)
            p = (x0 - x2) * q - (x0 - x1) * r
            q = 2.0 * (q - r)
            if 0.0 < q:
                p = - p
            q = abs(q)
            r = e
            e = d

        if abs(p) < abs(0.5*q*r) and q*(a-x0) < p and p < q * (b-x0):
            # Take the polynomial interpolation step.
            d = p / q
            u = x0 + d

            # Function must not be evaluated too close to a or b.
            if (u - a) < tol2 or (b - u) < tol2:
                if x0 < m:
                    d = tol
                else:
                    d = - tol
        else:
            # Take the golden-section step.
            if x0 < m:
                e = b - x0
            else:
                e = a - x0
            d = _golden * e

        # Function must not be evaluated too close to x0.
        if tol <= abs(d):
            u = x0 + d
        elif 0.0 < d:
            u = x0 + tol
        else:
            u = x0 - tol

        # Notice that we have u \in [a+tol, x0-tol] or
        #                     u \in [x0+tol, b-tol],
        # (if one ignores rounding errors.)
        fu = f(u)

        # Housekeeping.

        # Is the most recently evaluated point better (or equal) than the
        # best so far?
        if fu <= fx0:

            # Decrease interval size.
            if u < x0:
                b = x0
            else:
                a = x0

            # Shift: drop the previous third best point out and
            # include the newest point (found to be the best so far).
            x2 = x1
            fx2 = fx1
            x1 = x0
            fx1 = fx0
            x0 = u
            fx0 = fu

        else:
            # Decrease interval size.
            if u < x0:
                a = u
            else:
                b = u

            # Is the most recently evaluated point at better (or equal)
            # than the second best one?
            if fu <= fx1 or x1 == x0:
                # Insert u between (rank-wise) x0 and x1 in the triple
                # (x0, x1, x2).
                x2 = x1
                fx2 = fx1
                x1 = u
                fx1 = fu
            elif fu <= fx2 or x2 == x0 or x2 == x1:
                # Insert u in the last position of the triple (x0, x1, x2).
                x2 = u
                fx2 = fu

    return x0, fx0, niters+1
