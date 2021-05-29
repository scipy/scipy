import numpy as np


try:
    import mpmath
except ImportError:
    pass


def mp_hyp2f1(a, b, c, z):
    """Return mpmath hyp2f1 calculated on same branch as scipy hyp2f1.

    For most values of a,b,c mpmath returns the x - 0j branch of hyp2f1 on the
    branch cut x=(1,inf) whereas scipy's hyp2f1 calculates the x + 0j branch.
    Thus, to generate the right comparison values on the branch cut, we
    evaluate mpmath.hyp2f1 at x + 1e-15*j.

    The exception to this occurs when c-a=-m in which case both mpmath and
    scipy calculate the x + 0j branch on the branch cut. When this happens
    mpmath.hyp2f1 will be evaluated at the original z point.
    """
    on_branch_cut = z.real > 1.0 and abs(z.imag) < 1.0e-15
    cond1 = abs(c - a - round(c - a)) < 1.0e-15 and round(c - a) <= 0
    cond2 = abs(c - b - round(c - b)) < 1.0e-15 and round(c - b) <= 0
    # Make sure imaginary part is *exactly* zero
    if on_branch_cut:
        z = z.real + 0.0j
    if on_branch_cut and not (cond1 or cond2):
        z_mpmath = z.real + 1.0e-15j
    else:
        z_mpmath = z
    return complex(mpmath.hyp2f1(a, b, c, z_mpmath))




def _get_rand(minv, maxv, only_int=False, exclude_int=False):
    """Gets a random value in the range (minv, maxv).
    When only_int=True, this function only returns integers in (minv, maxv).
    When exclude_int=True, this function only returns *non*-integers in
    (minv, maxv).
    """
    if only_int:
        randv = np.random.randint(minv, maxv + 1)
    else:
        keep_going = True
        count = 0
        while keep_going and count < 100:
            randv = minv + (maxv - minv) * np.random.rand()
            # Get a new randv if we are excluding integers and this value is
            # too close to an int
            keep_going = (
                exclude_int and
                np.abs(np.round(randv) - randv) < 1.0e-2
            )
            count += 1

    return randv


def _get_rnd_tuple(*args):
    """Returns a tuple containing randomly generated values.

    The returned tuple will contain as many random values as there are
    arguments to this function. Each argument to _get_rnd_tuple should be a
    3-tuple of the form:

    (flag, minv, maxv)

    flag: i -> this value will be a random integer
            f -> this value will be a floating point number BUT not an integer
            a -> this value will be any random number between minv, maxv
                    (integers included)
    """
    curr_params = []

    for p in args:
        flag = p[0]
        flag = flag.capitalize()

        minv = p[1]
        maxv = p[2]

        # I -> pick a random integer
        if flag == "I":
            randv = _get_rand(minv, maxv, only_int=True)
        # F -> pick a floating point number (exclude integers)
        elif flag == "F":
            randv = _get_rand(minv, maxv, exclude_int=True)
        else:
            randv = _get_rand(minv, maxv)

        curr_params.append(randv)

    return tuple(curr_params)


def _add_params(plist, *args):
    rnd_tuple = _get_rnd_tuple(*args)
    plist.append(rnd_tuple)
    return


# Build a list of a, b, c values that covers each code branch
# in specfun.f
def _build_abc_list():
    small = 1.0e-2
    plist = []

    # +a, +b, +c test case
    _add_params(
        plist, ("f", small, 3.0), ("f", small, 3.0), ("f", small, 3.0)
    )

    # -a, +b, +c test case
    _add_params(
        plist, ("f", -3.0, -small), ("f", small, 3.0), ("f", small, 3.0)
    )

    # -a, -b, +c test case
    _add_params(
        plist, ("f", -3.0, -small), ("f", -3.0, -small), ("f", small, 3.0)
    )

    # -a, -b, -c test case
    _add_params(
        plist, ("f", -3.0, -small), ("f", -3.0, -small), ("f", -3.0, -small)
    )

    # Re(c-a-b)>0 test case
    _add_params(
        plist, ("f", small, 2.0), ("f", small, 2.0), ("f", 4.0 + small, 6.0)
    )

    # Re(c-a-b)<-1 test case
    _add_params(
        plist, ("f", 2.0, 4.0), ("f", 2.0, 4.0), ("f", small, 3.0 - small)
    )

    # c-a-b=m>0
    a, b, m = _get_rnd_tuple(("f", small, 2.0), ("f", small, 2.0), ("i", 1, 6))
    c = a + b + m
    plist.append((a, b, c))

    # c-a-b=0
    a, b = _get_rnd_tuple(("f", small, 2.0), ("f", small, 2.0))
    c = a + b
    plist.append((a, b, c))

    # c-a-b=-m<0
    a, b, m = _get_rnd_tuple(("f", 3.0, 5.0), ("f", 3.0, 5.0), ("i", 1, 6))
    c = a + b - m
    plist.append((a, b, c))

    # b-a=m>0
    a, m, c = _get_rnd_tuple(("f", small, 2.0), ("i", 1, 5), ("f", small, 2.0))
    b = a + m
    plist.append((a, b, c))

    # b-a=0
    a, c = _get_rnd_tuple(("f", small, 3.0), ("f", small, 3.0))
    b = a
    plist.append((a, b, c))

    # b-a=-m<0
    a, m, c = _get_rnd_tuple(("f", 4.0, 6.0), ("i", 1, 4), ("f", small, 2.0))
    b = a - m
    plist.append((a, b, c))

    # c-a=m>0
    a, b, m = _get_rnd_tuple(("f", small, 3.0), ("f", small, 3.0), ("i", 1, 5))
    c = a + m
    plist.append((a, b, c))

    # c-a=0
    a, b = _get_rnd_tuple(("f", small, 3.0), ("f", small, 3.0))
    c = a
    plist.append((a, b, c))

    # c-a=-m<0
    a, b, m = _get_rnd_tuple(("f", 4.0, 6.0), ("f", small, 3.0), ("i", 1, 4))
    c = a - m
    plist.append((a, b, c))

    # a=-m
    m, b, c = _get_rnd_tuple(
        ("i", 1, 10), ("f", small, 3.0), ("f", small, 3.0)
    )
    a = -m
    plist.append((a, b, c))

    # c=-m, a=-n, |m|>|n|

    # b, m, n = _get_rnd_tuple(('f',small,3.0),('i',5,8),('i',2,4))
    # a = -n
    # c = -m
    # plist.append((a,b,c))

    return plist


# Build the test case data that FuncData needs
def _build_test_cases(rho, phi):
    abc_list = _build_abc_list()

    # The total number of test cases
    N = 2 * len(abc_list)
    dataset = np.zeros((N, 5), dtype=complex)

    count = 0
    for a, b, c in abc_list:
        z = rho * np.exp(1.0j * phi)
        right = mpmath_hyp2f1_wrap(a, b, c, z)
        dataset[count, :] = a, b, c, z, right
        count += 1

        # Make sure that scipy's hyp2f1 returns same result
        # when we swap a, b
        dataset[count, :] = b, a, c, z, right
        count += 1

    return dataset
