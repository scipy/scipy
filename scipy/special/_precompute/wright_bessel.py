"""Precompute coefficients of several series expansions
of Wright's generalized Bessel function Phi(a, b, x).

See https://dlmf.nist.gov/10.46.E1 with rho=a, beta=b, z=x.
"""
from argparse import ArgumentParser, RawTextHelpFormatter

try:
    import sympy
    from sympy import EulerGamma, Rational, S, Sum, \
        factorial, gamma, gammasimp, pi, polygamma, symbols, zeta
except ImportError:
    pass


def series_small_a():
    """Tylor series expansion of Phi(a, b, x) in a=0 up to order 5.
    """
    order = 5
    a, b, x, k = symbols("a b x k")
    A = []  # terms with a
    X = []  # terms with x
    B = []  # terms with b (polygammas)
    # Phi(a, b, x) = exp(x)/gamma(b) * sum(A[i] * X[i] * B[i])
    expression = Sum(x**k/factorial(k)/gamma(a*k+b), (k, 0, S.Infinity))
    expression = gamma(b)/sympy.exp(x) * expression

    # nth term of taylor series in a=0: a^n/n! * (d^n Phi(a, b, x)/da^n at a=0)
    for n in range(0, order+1):
        term = expression.diff(a, n).subs(a, 0).simplify().doit()
        # set the whole bracket involving polygammas to 1
        x_part = (term.subs(polygamma(0, b), 1)
                  .replace(polygamma, lambda *args: 0))
        # sign convetion: x part always positive
        x_part *= (-1)**n

        A.append(a**n/factorial(n))
        X.append(x_part)
        B.append((term/x_part).simplify())

    s = "Tylor series expansion of Phi(a, b, x) in a=0 up to order 5.\n"
    s += "Phi(a, b, x) = exp(x)/gamma(b) * sum(A[i] * X[i] * B[i], i=0..5)\n"
    for name, c in zip(['A', 'X', 'B'], [A, X, B]):
        for i in range(len(c)):
            s += f"\n{name}[{i}] = " + str(c[i])
    return s


# expansion of digamma
def dg_series(z, n):
    """Symbolic expansion of digamma(z) in z=0 to order n.

    See https://dlmf.nist.gov/5.7.E4 and with https://dlmf.nist.gov/5.5.E2
    """
    k = symbols("k")
    return -1/z - EulerGamma + \
        sympy.summation((-1)**k * zeta(k) * z**(k-1), (k, 2, n+1))


def pg_series(k, z, n):
    """Symbolic expansion of polygamma(k, z) in z=0 to order n."""
    return sympy.diff(dg_series(z, n+k), z, k)


def series_small_a_small_b():
    """Tylor series expansion of Phi(a, b, x) in a=0 and b=0 up to order 5.

    Be aware of cancellation of poles in b=0 of digamma(b)/Gamma(b) and
    polygamma functions.

    digamma(b)/Gamma(b) = −1 − 2*M_EG*b + O(b^2)
    digamma(b)^2/Gamma(b) = 1/b + 3*M_EG + b*(−5/12*PI^2+7/2*M_EG^2) + O(b^2)
    polygamma(1, b)/Gamma(b) = 1/b + M_EG + b*(1/12*PI^2 + 1/2*M_EG^2) + O(b^2)
    and so on.
    """
    order = 5
    a, b, x, k = symbols("a b x k")
    M_PI, M_EG, M_Z3 = symbols("M_PI M_EG M_Z3")
    c_subs = {pi: M_PI, EulerGamma: M_EG, zeta(3): M_Z3}
    A = []  # terms with a
    X = []  # terms with x
    B = []  # terms with b (polygammas expanded)
    C = []  # terms that generate B
    # Phi(a, b, x) = exp(x) * sum(A[i] * X[i] * B[i])
    # B[0] = 1
    # B[k] = sum(C[k] * b**k/k!, k=0..)
    # Note: C[k] can be obtained from a series expansion of 1/gamma(b).
    expression = gamma(b)/sympy.exp(x) * \
        Sum(x**k/factorial(k)/gamma(a*k+b), (k, 0, S.Infinity))

    # nth term of taylor series in a=0: a^n/n! * (d^n Phi(a, b, x)/da^n at a=0)
    for n in range(0, order+1):
        term = expression.diff(a, n).subs(a, 0).simplify().doit()
        # set the whole bracket involving polygammas to 1
        x_part = (term.subs(polygamma(0, b), 1)
                  .replace(polygamma, lambda *args: 0))
        # sign convetion: x part always positive
        x_part *= (-1)**n
        # expansion of polygamma part with 1/gamma(b)
        pg_part = term/x_part/gamma(b)
        if n >= 1:
            # Note: highest term is digamma^n
            pg_part = pg_part.replace(polygamma,
                                      lambda k, x: pg_series(k, x, order+1+n))
            pg_part = (pg_part.series(b, 0, n=order+1-n)
                       .removeO()
                       .subs(polygamma(2, 1), -2*zeta(3))
                       .simplify()
                       )

        A.append(a**n/factorial(n))
        X.append(x_part)
        B.append(pg_part)

    # Calculate C and put in the k!
    C = sympy.Poly(B[1].subs(c_subs), b).coeffs()
    C.reverse()
    for i in range(len(C)):
        C[i] = (C[i] * factorial(i)).simplify()

    s = "Tylor series expansion of Phi(a, b, x) in a=0 and b=0 up to order 5."
    s += "\nPhi(a, b, x) = exp(x) * sum(A[i] * X[i] * B[i], i=0..5)\n"
    s += "B[0] = 1\n"
    s += "B[i] = sum(C[k+i-1] * b**k/k!, k=0..)\n"
    s += "\nM_PI = pi"
    s += "\nM_EG = EulerGamma"
    s += "\nM_Z3 = zeta(3)"
    for name, c in zip(['A', 'X', 'C'], [A, X, C]):
        for i in range(len(c)):
            s += f"\n{name}[{i}] = "
            s += str(c[i])

    # Does B have the assumed structure?
    s += "\n\nTest if B[i] does have the assumed structure."
    s += "\nC[i] are derived from B[1] allone."
    s += "\nTest B[2] == C[1] + b*C[2] + b^2/2*C[3] + b^3/6*C[4] + .."
    test = sum([b**k/factorial(k) * C[k+1] for k in range(order-1)])
    test = (test - B[2].subs(c_subs)).simplify()
    s += f"\ntest successful = {test==S(0)}"
    s += "\nTest B[3] == C[2] + b*C[3] + b^2/2*C[4] + .."
    test = sum([b**k/factorial(k) * C[k+2] for k in range(order-2)])
    test = (test - B[3].subs(c_subs)).simplify()
    s += f"\ntest successful = {test==S(0)}"
    return s


def asymptotic_series():
    """Asymptotic expansion for large x

    Phi(a, b, x) ~ Z^(1/2-b) * exp((1+a)/a * Z) * sum_k (-1)^k * C_k / Z^k
    Z = (a*x)^(1/(1+a))

    Wright (1935) lists the coefficients C_0 and C_1 (he calls them a_0 and
    a_1). With slightly different notation, Paris (2017) lists coefficients
    c_k up to order k=3.
    Paris (2017) uses ZP = (1+a)/a * Z  (ZP = Z of Paris) and
    C_k = C_0 * (-a/(1+a))^k * c_k
    """
    order = 6

    class g(sympy.Function):
        """Helper function g according to Wright (1935)

        g(n, rho, v) = (1 + (rho+2)/3 * v + (rho+2)*(rho+3)/(2*3) * v^2 + ...)

        Note: Wright (1935) uses square root of above definition.
        """
        nargs = 3

        @classmethod
        def eval(cls, n, rho, v):
            if not n >= 0:
                raise ValueError("must have n >= 0")
            elif n == 0:
                return 1
            else:
                return g(n-1, rho, v) \
                    + gammasimp(gamma(rho+2+n)/gamma(rho+2)) \
                    / gammasimp(gamma(3+n)/gamma(3))*v**n

    class coef_C(sympy.Function):
        """Calculate coefficients C_m for integer m.

        C_m is the coefficient of v^(2*m) in the Taylor expansion in v=0 of
        Gamma(m+1/2)/(2*pi) * (2/(rho+1))^(m+1/2) * (1-v)^(-b)
            * g(rho, v)^(-m-1/2)
        """
        nargs = 3

        @classmethod
        def eval(cls, m, rho, beta):
            if not m >= 0:
                raise ValueError("must have m >= 0")

            v = symbols("v")
            expression = (1-v)**(-beta) * g(2*m, rho, v)**(-m-Rational(1, 2))
            res = expression.diff(v, 2*m).subs(v, 0) / factorial(2*m)
            res = res * (gamma(m + Rational(1, 2)) / (2*pi)
                         * (2/(rho+1))**(m + Rational(1, 2)))
            return res

    # in order to have nice ordering/sorting of expressions, we set a = xa.
    xa, b, xap1 = symbols("xa b xap1")
    C0 = coef_C(0, xa, b)
    # a1 = a(1, rho, beta)
    s = "Asymptotic expansion for large x\n"
    s += "Phi(a, b, x) = Z**(1/2-b) * exp((1+a)/a * Z) \n"
    s += "               * sum((-1)**k * C[k]/Z**k, k=0..6)\n\n"
    s += "Z      = pow(a * x, 1/(1+a))\n"
    s += "A[k]   = pow(a, k)\n"
    s += "B[k]   = pow(b, k)\n"
    s += "Ap1[k] = pow(1+a, k)\n\n"
    s += "C[0] = 1./sqrt(2. * M_PI * Ap1[1])\n"
    for i in range(1, order+1):
        expr = (coef_C(i, xa, b) / (C0/(1+xa)**i)).simplify()
        factor = [x.denominator() for x in sympy.Poly(expr).coeffs()]
        factor = sympy.lcm(factor)
        expr = (expr * factor).simplify().collect(b, sympy.factor)
        expr = expr.xreplace({xa+1: xap1})
        s += f"C[{i}] = C[0] / ({factor} * Ap1[{i}])\n"
        s += f"C[{i}] *= {str(expr)}\n\n"
    import re
    re_a = re.compile(r'xa\*\*(\d+)')
    s = re_a.sub(r'A[\1]', s)
    re_b = re.compile(r'b\*\*(\d+)')
    s = re_b.sub(r'B[\1]', s)
    s = s.replace('xap1', 'Ap1[1]')
    s = s.replace('xa', 'a')
    # max integer = 2^31-1 = 2,147,483,647. Solution: Put a point after 10
    # or more digits.
    re_digits = re.compile(r'(\d{10,})')
    s = re_digits.sub(r'\1.', s)
    return s


def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument('action', type=int, choices=[1, 2, 3],
                        help='chose what expansion to precompute\n'
                             '1 : Series for small a\n'
                             '2 : Series for small a and small b\n'
                             '3 : Asymptotic series for large x\n'
                             '    This may take some time.')
    args = parser.parse_args()

    switch = {1: lambda: print(series_small_a()),
              2: lambda: print(series_small_a_small_b()),
              3: lambda: print(asymptotic_series()),
              }
    switch.get(args.action, lambda: print("Invalid input."))()


if __name__ == '__main__':
    main()
