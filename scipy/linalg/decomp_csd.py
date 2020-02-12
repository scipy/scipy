from numpy.linalg import LinAlgError

from .lapack import _compute_lwork, get_lapack_funcs

__all__ = ['cs_decomp']

def cs_decomp(X, p, q,
              comp_u1=True, comp_u2=True,
              comp_v1h=True, comp_v2h=True,
              sign_conv=''):
    """
        Computes cosine-sine (CS) decomposition of an M-by-M partitioned unitary
        matrix X, s.t.

                                         [  I  0  0 |  0  0  0 ]
                                         [  0  C  0 |  0 -S  0 ]
             [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**H
         X = [-----------] = [---------] [---------------------] [---------]   .
             [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
                                         [  0  S  0 |  0  C  0 ]
                                         [  0  0  I |  0  0  0 ]

        for the identity matrices see https://www.nag.com/numeric/fl/nagdoc_fl26/pdf/f08/f08rnf.pdf

        :param X:           square unitary to decompose
        :param p:           number of rows in upper left block
        :param q:           number of columns in the upper left block
        :param comp_u1:     compute unitary u1
        :param comp_u2:     compute unitary u2
        :param comp_v1h:    compute unitary v1^H
        :param comp_v2h:    compute unitary v2^H
        :param sign_conv:   'low'       minus sign is put in lower left block
                            otherwise:  minus sign is put in upper right block

        :returns u1:        (optional) upper left unitary block u1
        :returns u2:        (optional) lower right unitary block u2
        :returns v1h:       (optional) upper left unitary block v1^H
        :returns v2h:       (optional) lower right unitary block v2^H
        :returns theta:     angles in radians to construct
                            C = diag(cos(theta_1), ..., cos(theta_2))
                            S = diag(sin(theta_1), ..., sin(theta_2))

        Attribution to @anuesseller with his work in
        https://github.com/anuesseler/pycsd
    """

    m = X.shape[0]

    x11 = X[:p, :q]
    ldx11 = p

    x12 = X[:p, q:]
    ldx12 = p

    x21 = X[p:, :q]
    ldx21 = m - p

    x22 = X[p:, q:]
    ldx22 = m - p

    ldu1 = p
    ldu2 = m - p

    ldv1t = q
    ldv2t = m - q

    # TODO(balopat): maybe we can make these booleans instead?
    jobu1 = 'Y' if comp_u1 else ''
    jobu2 = 'Y' if comp_u2 else ''
    jobv1t = 'Y' if comp_v1h else ''
    jobv2t = 'Y' if comp_v2h else ''
    signs = 'O' if sign_conv == 'low' else ''

    kind = X.dtype.kind
    if kind == 'f':
        return _orthogonal_csd(X, jobu1, jobu2, jobv1t, jobv2t, ldu1, ldu2,
                               ldv1t,
                               ldv2t, ldx11, ldx12, ldx21, ldx22, m, p, q,
                               signs, x11,
                               x12, x21, x22)
    elif kind == 'c':
        return _unitary_csd(X, jobu1, jobu2, jobv1t, jobv2t, ldu1, ldu2, ldv1t,
                            ldv2t, ldx11, ldx12, ldx21, ldx22, m, p, q, signs,
                            x11,
                            x12, x21, x22)
    raise Exception("can't kind '{}' of datatype {}".format(kind, X.dtype))


def _orthogonal_csd(X, jobu1, jobu2, jobv1t, jobv2t, ldu1, ldu2, ldv1t, ldv2t,
                    ldx11, ldx12, ldx21, ldx22, m, p, q, signs, x11, x12, x21,
                    x22):
    csd, csd_lwork = get_lapack_funcs(["orcsd", "orcsd_lwork"], [X[0][0]])
    lwork = _compute_lwork(csd_lwork, m, p, q, x11, ldx11, x12, ldx12, x21,
                           ldx21,
                           x22, ldx22, ldu1, ldu2, ldv1t, ldv2t, jobu1=jobu1,
                           jobu2=jobu2, jobv1t=jobv1t, jobv2t=jobv2t, trans='',
                           signs=signs)

    theta, u1, u2, v1h, v2h, _, info = \
        csd(m, p, q,
            x11, ldx11,
            x12, ldx12,
            x21, ldx21,
            x22, ldx22,
            ldu1, ldu2,
            ldv1t, ldv2t,
            lwork[0],
            jobu1=jobu1,
            jobu2=jobu2,
            jobv1t=jobv1t,
            jobv2t=jobv2t,
            trans='',
            signs=signs)
    method_name = csd.typecode + "orcsd"
    if info < 0:
        raise ValueError('illegal value in argument %d of internal %s'
                         % (-info, method_name))
    if info > 0:
        raise LinAlgError("%s did not converge: %d"
                          % (info, method_name))
    return u1, u2, v1h, v2h, theta


def _unitary_csd(X, jobu1, jobu2, jobv1t, jobv2t, ldu1, ldu2, ldv1t, ldv2t,
                 ldx11, ldx12, ldx21, ldx22, m, p, q, signs, x11, x12, x21,
                 x22):
    csd, csd_lwork = get_lapack_funcs(["uncsd", "uncsd_lwork"], [X[0][0]])
    print("UNITARY {}".format(csd.typecode))

    lwork = _compute_lwork(csd_lwork, m, p, q, x11, ldx11, x12, ldx12, x21,
                           ldx21,
                           x22, ldx22, ldu1, ldu2, ldv1t, ldv2t, jobu1=jobu1,
                           jobu2=jobu2, jobv1t=jobv1t, jobv2t=jobv2t, trans='',
                           signs=signs)

    theta, u1, u2, v1h, v2h, info = \
        csd(m, p, q,
            x11, ldx11,
            x12, ldx12,
            x21, ldx21,
            x22, ldx22,
            ldu1, ldu2,
            ldv1t, ldv2t,
            # TODO(balopat) I guess this works?
            lwork[0], lwork[1],
            jobu1=jobu1,
            jobu2=jobu2,
            jobv1t=jobv1t,
            jobv2t=jobv2t,
            trans='',
            signs=signs)
    method_name = csd.typecode + "uncsd"

    if info < 0:
        raise ValueError('illegal value in argument %d of internal %s'
                         % (-info, method_name))
    if info > 0:
        raise LinAlgError("%s did not converge: %d"
                          % (info, method_name))

    return u1, u2, v1h, v2h, theta
