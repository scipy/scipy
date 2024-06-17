import collections
import inspect
import numbers
import numpy as np

from ._input_validation import _nonneg_int_or_fail

from ._special_ufuncs import (legendre_p, assoc_legendre_p,
                              multi_assoc_legendre_p, sph_legendre_p,
                              sph_harm_y)

from ._gufuncs import (legendre_p_all, assoc_legendre_p_all,
                       multi_assoc_legendre_p_all,
                       sph_legendre_p_all, sph_harm_y_all)

__all__ = [
    "assoc_legendre_p",
    "assoc_legendre_p_all",
    "legendre_p",
    "legendre_p_all",
    "multi_assoc_legendre_p",
    "multi_assoc_legendre_p_all",
    "sph_harm_y",
    "sph_harm_y_all",
    "sph_legendre_p",
    "sph_legendre_p_all",
]


class MultiUFunc:
    def __init__(self, ufuncs_map, doc = None, *, force_complex_output=False):
        # Gather leaf level ufuncs from ufuncs_map.
        ufuncs = []
        def traverse(obj):
            if isinstance(obj, collections.abc.Mapping):
                for value in obj.values():
                    traverse(value)
            elif (isinstance(obj, collections.abc.Iterable)
                and not isinstance(obj, (str, bytes))):
                for item in obj:
                    traverse(item)
            else:
                ufuncs.append(obj)
        traverse(ufuncs_map)

        # Perform input validation to ensure all ufuncs in ufuncs_map are
        # actually ufuncs and all take the same input types.
        seen_input_types = set()
        for ufunc in ufuncs:
            if not isinstance(ufunc, np.ufunc):
                raise ValueError("All leaf elements of ufuncs_map must have"
                                 f" type `numpy.ufunc`. Received {ufuncs_map}")
            seen_input_types.add(frozenset(x.split("->")[0] for x in ufunc.types))
        if len(seen_input_types) > 1:
            raise ValueError("All ufuncs in ufuncs_map must take the same"
                             " input types.")

        self._ufuncs = tuple(ufuncs)
        self._ufuncs_map = ufuncs_map
        self.resolve_out_shapes = None
        self.resolve_ufunc = None
        self.__force_complex_output = force_complex_output
        self.__doc = doc

    @property
    def force_complex_output(self):
        return self.__force_complex_output

    def register_resolve_ufunc(self, func):
        """Set `resolve_ufunc` method by decorating a function.

        The decorated function's first argument should be a JSON-like
        `ufuncs_map`, and additional arguments should be keywords which 
        are used to dispatch to the correct ufunc at the leaf level of
        `ufuncs_map`.
        """
        # Given func, we construct a wrapper which no longer takes the
        # `ufuncs_map` as input, but instead gets it from the
        # class. Use inspect to add an informative signature. Add a
        # a docstring if none exists.
        sig = inspect.signature(func)
        params = list(inspect.signature(func).parameters.values())[1:]
        new_sig = sig.replace(parameters=params)

        def resolve_ufunc(**kwargs):
            return func(self._ufuncs_map, **kwargs)

        resolve_ufunc.__signature__ = new_sig
        docstring = func.__doc__
        if docstring is not None:
            resolve_ufunc.__doc__ = docstring
        resolve_ufunc.__doc__ = \
            """Resolve to a ufunc based on keyword arguments."""
        self.resolve_ufunc = resolve_ufunc

    def register_resolve_out_shapes(self, func):
        """Set `resolve_out_shapes` method by decorating a function."""
        if func.__doc__ is None:
            func.__doc__ = \
                """Resolve to output shapes based on relevant inputs."""
        func.__name__ = "resolve_out_shapes"
        self.resolve_out_shapes = func

    @property
    def __doc__(self):
        return self.__doc

    def __call__(self, *args, **kwargs):
        ufunc = self.resolve_ufunc(**kwargs)
        if ((ufunc.nout == 0) or (self.resolve_out_shapes is None)):
            return ufunc(*args)

        ufunc_args = args[-ufunc.nin:] # array arguments to be passed to the ufunc

        ufunc_arg_shapes = tuple(np.shape(ufunc_arg) for ufunc_arg in ufunc_args)
        ufunc_out_shapes = self.resolve_out_shapes(*args[:-ufunc.nin],
            *ufunc_arg_shapes, ufunc.nout)

        ufunc_arg_dtypes = tuple((ufunc_arg.dtype if hasattr(ufunc_arg, 'dtype')
            else np.dtype(type(ufunc_arg))) for ufunc_arg in ufunc_args)
        if hasattr(ufunc, 'resolve_dtypes'):
            ufunc_dtypes = ufunc_arg_dtypes + ufunc.nout * (None,)
            ufunc_dtypes = ufunc.resolve_dtypes(ufunc_dtypes) 
            ufunc_out_dtypes = ufunc_dtypes[-ufunc.nout:]
        else:
            ufunc_out_dtype = np.result_type(*ufunc_arg_dtypes)
            if (not np.issubdtype(ufunc_out_dtype, np.inexact)):
                ufunc_out_dtype = np.float64

            ufunc_out_dtypes = ufunc.nout * (ufunc_out_dtype,)

        if self.force_complex_output:
            ufunc_out_dtypes = tuple(np.result_type(1j, ufunc_out_dtype)
                for ufunc_out_dtype in ufunc_out_dtypes)

        b = np.broadcast(*ufunc_args)
        ufunc_out_new_dims = tuple(len(ufunc_out_shape) - b.ndim
            for ufunc_out_shape in ufunc_out_shapes)

        out = tuple(np.empty(ufunc_out_shape, dtype = ufunc_out_dtype)
            for ufunc_out_shape, ufunc_out_dtype
            in zip(ufunc_out_shapes, ufunc_out_dtypes))

        ufunc_out = tuple(np.moveaxis(out[i],
            tuple(range(axis)), tuple(range(-axis, 0)))
            for i, axis in enumerate(ufunc_out_new_dims))
        ufunc(*ufunc_args, out = ufunc_out)

        if (len(out) == 1):
            out, = out

        return out


assoc_legendre_p = MultiUFunc(assoc_legendre_p,
    r"""assoc_legendre_p(n, m, z, *, norm=False, diff_n=0)

    Associated Legendre polynomial of the first kind.

    Parameters
    ----------
    n : array_like, int
        Order of the Legendre polynomial, must have ``n >= 0``.
    m : array_like, int
        Degree of the Legendre polynomial.
    z : array_like, float
        Input value.
    norm : Optional[bool]
        If True, compute the normalized associated Legendre polynomial.
        Default is False.
    diff_n : Optional[int]
        A non-negative integer. Compute and return all derivatives up
        to order ``diff_n``.

    Returns
    -------
    p : ndarray or tuple[ndarray]
        The assocated Legendre polynomial with ``diff_n`` derivatives.

    Notes
    -----
    With respect to the associated Legendre polynomial :math:`P_{m}^{n}(x)`,
    the normalised associated Legendre polynomials is defined as

    .. math::

        \sqrt{\frac{(2 n + 1) (n - m)!}{2 (n + m)!}} P_{n}^{m}(x)
    """
)


@assoc_legendre_p.register_resolve_ufunc
def _(ufuncs, norm=False, diff_n=0):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return ufuncs[norm][diff_n]


assoc_legendre_p_all = MultiUFunc(assoc_legendre_p_all,
    r"""assoc_legendre_p_all(n, m, z, *, norm=False, diff_n=0)

    Associated Legendre polynomial of the first kind.

    Parameters
    ----------
    n : int
        Order (max) of the Legendre polynomial, must have ``n >= 0``.
    m : int
        Degree (max) of the Legendre polynomial, must have ``m >= 0``.
    z : array_like, float
        Input value
    norm : Optional[bool]
        If True, compute the normalized associated Legendre polynomial.
        Default is False.
    diff_n : Optional[int]
        A non-negative integer. Compute and return all derivatives up
        to order ``diff_n``.

    Returns
    -------
    p : ndarray or tuple[ndarray]
        The assocated Legendre polynomial with ``diff_n`` derivatives,
        each having shape (n + 1, 2 * m + 1, ...).

    Notes
    -----
    With respect to the associated Legendre polynomial :math:`P_{m}^{n}(x)`,
    the normalised associated Legendre polynomial is defined as

    .. math::

        \sqrt{\frac{(2 n + 1) (n - m)!}{2 (n + m)!}} P_{n}^{m}(x)
    """
)

@assoc_legendre_p_all.register_resolve_ufunc
def _(ufuncs, norm=False, diff_n=0):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    return ufuncs[norm][diff_n]


@assoc_legendre_p_all.register_resolve_out_shapes
def _(n, m, z_shape, nout):
    if ((not np.isscalar(m)) or (abs(m) > n)):
        raise ValueError("m must be <= n.")

    if ((not np.isscalar(n)) or (n < 0)):
        raise ValueError("n must be a non-negative integer.")

    return nout * ((n + 1, 2 * abs(m) + 1) + z_shape,)


sph_legendre_p = MultiUFunc(sph_legendre_p,
    r"""sph_legendre_p(n, m, phi, *, diff_n=0)

    Spherical Legendre polynomial of the first kind.

    Parameters
    ----------
    n : array_like, int
        Order of the Legendre polynomial, must have ``n >= 0``.
    m : array_like, int
        Degree of the Legendre polynomial.
    phi : array_like, float
        Input value.
    diff_n : Optional[int]
        A non-negative integer. Compute and return all derivatives up
        to order ``diff_n``.

    Returns
    -------
    p : ndarray or tuple[ndarray]
        The spherical Legendre polynomial with ``diff_n`` derivatives.

    Notes
    -----
    With respect to the associated Legendre polynomial :math:`P_{m}^{n}(\cos \phi)`,
    the spherical Legendre polynomial is defined as

    .. math::

        \sqrt{\frac{(2 n + 1) (n - m)!}{4 \pi (n + m)!}} P_{n}^{m}(\cos \phi)

    This is the same as the spherical harmonic :math:`Y_{m}^{n}(\theta, \phi)`
    with :math:`\theta = 0`.
    """
)


@sph_legendre_p.register_resolve_ufunc
def _(ufuncs, diff_n = 0):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return ufuncs[diff_n]


sph_legendre_p_all = MultiUFunc(sph_legendre_p_all,
    r"""sph_legendre_p_all(n, m, phi, *, diff_n=0)

    Spherical Legendre polynomial of the first kind.

    Parameters
    ----------
    n : array_like, int
        Order (max) of the Legendre polynomial, must have ``n >= 0``.
    m : array_like, int
        Degree (max) of the Legendre polynomial, must have ``m >= 0``
    phi : array_like, float
        Input value.
    diff_n : Optional[int]
        A non-negative integer. Compute and return all derivatives up
        to order ``diff_n``.

    Returns
    -------
    p : ndarray or tuple[ndarray]
        The spherical Legendre polynomial with ``diff_n`` derivatives,
        each having shape (n + 1, 2 * m + 1, ...).

    Notes
    -----
    With respect to the associated Legendre polynomial :math:`P_{m}^{n}(\cos \phi)`,
    the spherical Legendre polynomial is defined as

    .. math::

        \sqrt{\frac{(2 n + 1) (n - m)!}{4 \pi (n + m)!}} P_{n}^{m}(\cos \phi)

    This is the same as the spherical harmonic :math:`Y_{m}^{n}(\theta, \phi)`
    with :math:`\theta = 0`.
    """
)


@sph_legendre_p_all.register_resolve_ufunc
def _(ufuncs, diff_n=0):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return ufuncs[diff_n]


@sph_legendre_p_all.register_resolve_out_shapes
def _(n, m, z_shape, nout):
    if ((not np.isscalar(n)) or (n < 0)):
        raise ValueError("n must be a non-negative integer.")

    return nout * ((n + 1, 2 * abs(m) + 1) + z_shape,)


multi_assoc_legendre_p = MultiUFunc(multi_assoc_legendre_p,
    """
    multi_assoc_legendre_p_all(n, m, type, z, *, norm=False, diff_n=0)

    Associated Legendre function of the first kind for complex arguments.
    """, force_complex_output=True
)


@multi_assoc_legendre_p.register_resolve_ufunc
def _(ufuncs, norm=False, diff_n=0):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return ufuncs[norm][diff_n]


multi_assoc_legendre_p_all = MultiUFunc(multi_assoc_legendre_p_all,
    """
    multi_assoc_legendre_p_all(n, m, type, z, *, norm=False, diff_n=0)

    Table of associated Legendre functions of the first kind for complex
    arguments.
    """, force_complex_output=True
)


@multi_assoc_legendre_p_all.register_resolve_ufunc
def _(ufuncs, norm=False, diff_n=0):
    if not ((isinstance(diff_n, int) or np.issubdtype(diff_n, np.integer))
            and diff_n >= 0):
        raise ValueError(
            f"diff_n must be a non-negative integer, received: {diff_n}."
        )
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return ufuncs[norm][diff_n]


@multi_assoc_legendre_p_all.register_resolve_out_shapes
def _(n, m, type_shape, z_shape, nout):
    if not isinstance(m, numbers.Integral) or (abs(m) > n):
        raise ValueError("m must be <= n.")
    if not isinstance(n, numbers.Integral) or (n < 0):
        raise ValueError("n must be a non-negative integer.")

    return nout * ((n + 1, 2 * abs(m) + 1) + np.broadcast_shapes(type_shape, z_shape),)


legendre_p = MultiUFunc(legendre_p,
    """
    legendre_p(n, z, *, diff_n=0)

    Legendre polynomials of the first kind.

    Compute Legendre polynomials of the first kind 
    :math:`P_{n}(z)`.

    Parameters
    ----------
    n : array_like, int
        Order of the Legendre polynomial, must have ``n >= 0``.
    z : array_like, float
        Input value.
    diff_n : Optional[int]
        A non-negative integer. Compute and return all derivatives up
        to order ``diff_n``.

    Returns
    -------
    p : ndarray or tuple[ndarray]
        Legendre polynomial with ``diff_n`` derivatives.

    See also special.legendre for polynomial class.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           https://people.sc.fsu.edu/~jburkardt/f77_src/special_functions/special_functions.html
    """
)


@legendre_p.register_resolve_ufunc
def _(ufuncs, diff_n=0):
    if not ((isinstance(diff_n, int) or np.issubdtype(diff_n, np.integer))
            and diff_n >= 0):
        raise ValueError(
            f"diff_n must be a non-negative integer, received: {diff_n}."
        )
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return ufuncs[diff_n]


legendre_p_all = MultiUFunc(legendre_p_all,
    """
    Legendre polynomials of the first kind.

    Compute sequence of Legendre functions of the first kind (polynomials),
    Pn(z) and derivatives for all degrees from 0 to n (inclusive).

    See also special.legendre for polynomial class.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           https://people.sc.fsu.edu/~jburkardt/f77_src/special_functions/special_functions.html
    """
)


@legendre_p_all.register_resolve_ufunc
def _(ufuncs, diff_n=0):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return ufuncs[diff_n]


@legendre_p_all.register_resolve_out_shapes
def _(n, z_shape, nout):
    n = _nonneg_int_or_fail(n, 'n', strict=False)

    return nout * ((n + 1,) + z_shape,)


sph_harm_y = MultiUFunc(sph_harm_y,
    r"""
    sph_harm_y(n, m, theta, phi)

    Compute spherical harmonics.

    The spherical harmonics are defined as

    .. math::

        Y^m_n(\theta,\phi) = \sqrt{\frac{2n+1}{4\pi} \frac{(n-m)!}{(n+m)!}}
          e^{i m \theta} P^m_n(\cos(\phi))

    where :math:`P_n^m` are the associated (spherical) Legendre functions;
    see `sph_legendre_p`.

    Parameters
    ----------
    n : array_like
        Degree of the harmonic (int); must have ``n >= 0``. This is
        often denoted by ``l`` (lower case L) in descriptions of
        spherical harmonics.
    m : array_like
        Order of the harmonic (int); must have ``|m| <= n``.
    theta : array_like
        Azimuthal (longitudinal) coordinate; must be in ``[0, 2*pi]``.
    phi : array_like
        Polar (colatitudinal) coordinate; must be in ``[0, pi]``.

    Returns
    -------
    y : complex ndarray
       The harmonic :math:`Y^m_n` sampled at ``theta`` and ``phi``.

    Notes
    -----
    There are different conventions for the meanings of the input
    arguments ``theta`` and ``phi``. In SciPy ``theta`` is the
    azimuthal angle and ``phi`` is the polar angle. It is common to
    see the opposite convention, that is, ``theta`` as the polar angle
    and ``phi`` as the azimuthal angle.

    Note that SciPy's spherical harmonics include the Condon-Shortley
    phase [2]_ because it is part of `lpmv`.

    With SciPy's conventions, the first several spherical harmonics
    are

    .. math::

        Y_0^0(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{1}{\pi}} \\
        Y_1^{-1}(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{3}{2\pi}}
                                    e^{-i\theta} \sin(\phi) \\
        Y_1^0(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{3}{\pi}}
                                 \cos(\phi) \\
        Y_1^1(\theta, \phi) &= -\frac{1}{2} \sqrt{\frac{3}{2\pi}}
                                 e^{i\theta} \sin(\phi).

    References
    ----------
    .. [1] Digital Library of Mathematical Functions, 14.30.
           https://dlmf.nist.gov/14.30
    .. [2] https://en.wikipedia.org/wiki/Spherical_harmonics#Condon.E2.80.93Shortley_phase
    """, force_complex_output=True
)


@sph_harm_y.register_resolve_ufunc
def _(ufuncs, diff_n=0):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return ufuncs[diff_n]


sph_harm_y_all = MultiUFunc(sph_harm_y_all, 
    r"""
    sph_harm_y_all(n, m, theta, phi, *, diff_n=0)

    Compute all spherical harmonics up to degree n and order m.

    The spherical harmonics are defined as

    .. math::

        Y^m_n(\theta,\phi) = e^{i m \theta} P^m_n(\cos(\phi))

    where :math:`P_n^m` are the associated (spherical) Legendre polynomials;
    see `sph_legendre_p`.

    Parameters
    ----------
    n : int
        Degree (max) of the harmonics; must have ``n >= 0``. This is
        often denoted by ``l`` (lower case L) in descriptions of
        spherical harmonics.
    m : int
        Order (max) of the harmonics; must have ``m >= 0``.
    theta : array_like
        Azimuthal (longitudinal) coordinate; must be in ``[0, 2*pi]``.
    phi : array_like
        Polar (colatitudinal) coordinate; must be in ``[0, pi]``.
    diff_n : Optional[int]
        A non-negative integer. Return all derivatives up to
        order ``diff_n``.

    Returns
    -------
    y : ndarray or tuple of ndarray
        The complex harmonic :math:`Y^m_n` sampled at ``theta`` and ``phi``
        with n and m as the leading dimensions. 

    Notes
    -----
    There are different conventions for the meanings of the input
    arguments ``theta`` and ``phi``. In SciPy ``theta`` is the
    azimuthal angle and ``phi`` is the polar angle. It is common to
    see the opposite convention, that is, ``theta`` as the polar angle
    and ``phi`` as the azimuthal angle.

    Note that SciPy's spherical harmonics include the Condon-Shortley
    phase [2]_ because it is part of `lpmv`.

    With SciPy's conventions, the first several spherical harmonics
    are

    .. math::

        Y_0^0(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{1}{\pi}} \\
        Y_1^{-1}(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{3}{2\pi}}
                                    e^{-i\theta} \sin(\phi) \\
        Y_1^0(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{3}{\pi}}
                                 \cos(\phi) \\
        Y_1^1(\theta, \phi) &= -\frac{1}{2} \sqrt{\frac{3}{2\pi}}
                                 e^{i\theta} \sin(\phi).

    References
    ----------
    .. [1] Digital Library of Mathematical Functions, 14.30.
           https://dlmf.nist.gov/14.30
    .. [2] https://en.wikipedia.org/wiki/Spherical_harmonics#Condon.E2.80.93Shortley_phase
    """, force_complex_output=True
)


@sph_harm_y_all.register_resolve_ufunc
def _(ufuncs, diff_n=0):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 2,"
            f" received: {diff_n}."
        )
    return ufuncs[diff_n]


@sph_harm_y_all.register_resolve_out_shapes
def _(n, m, theta_shape, phi_shape, nout):
    if ((not np.isscalar(n)) or (n < 0)):
        raise ValueError("n must be a non-negative integer.")

    return tuple(diff_ndims * (2,) + (n + 1, 2 * abs(m) + 1) +
        np.broadcast_shapes(theta_shape, phi_shape) for diff_ndims in range(nout))