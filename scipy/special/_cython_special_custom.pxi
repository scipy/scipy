from . cimport _spherical_bessel


cpdef number_t spherical_jn(long n, number_t z, bint derivative=0) noexcept nogil:
    """See the documentation for scipy.special.spherical_jn"""
    if derivative:
        if number_t is double:
            return _spherical_bessel.spherical_jn_d_real(n, z)
        else:
            return _spherical_bessel.spherical_jn_d_complex(n, z)

    if number_t is double:
        return _spherical_bessel.spherical_jn_real(n, z)
    else:
        return _spherical_bessel.spherical_jn_complex(n, z)


cpdef number_t spherical_yn(long n, number_t z, bint derivative=0) noexcept nogil:
    """See the documentation for scipy.special.spherical_yn"""
    if derivative:
        if number_t is double:
            return _spherical_bessel.spherical_yn_d_real(n, z)
        else:
            return _spherical_bessel.spherical_yn_d_complex(n, z)

    if number_t is double:
        return _spherical_bessel.spherical_yn_real(n, z)
    else:
        return _spherical_bessel.spherical_yn_complex(n, z)


cpdef number_t spherical_in(long n, number_t z, bint derivative=0) noexcept nogil:
    """See the documentation for scipy.special.spherical_in"""
    if derivative:
        if number_t is double:
            return _spherical_bessel.spherical_in_d_real(n, z)
        else:
            return _spherical_bessel.spherical_in_d_complex(n, z)

    if number_t is double:
        return _spherical_bessel.spherical_in_real(n, z)
    else:
        return _spherical_bessel.spherical_in_complex(n, z)


cpdef number_t spherical_kn(long n, number_t z, bint derivative=0) noexcept nogil:
    """See the documentation for scipy.special.spherical_kn"""
    if derivative:
        if number_t is double:
            return _spherical_bessel.spherical_kn_d_real(n, z)
        else:
            return _spherical_bessel.spherical_kn_d_complex(n, z)

    if number_t is double:
        return _spherical_bessel.spherical_kn_real(n, z)
    else:
        return _spherical_bessel.spherical_kn_complex(n, z)
