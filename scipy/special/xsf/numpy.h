#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <array>
#include <cassert>
#include <cstring>
#include <iostream>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <numpy/arrayobject.h>
#include <numpy/npy_3kcompat.h>
#include <numpy/npy_math.h>
#include <numpy/ufuncobject.h>

#include "dual.h"
#include "error.h"
#include "third_party/kokkos/mdspan.hpp"

/* PyUFunc_getfperr gets bits for current floating point error (fpe) status codes so we
 * can check for floating point errors and make proper calls to set_error in ufunc loops.
 * Define a wrapper so it can be given C linkage within this C++ header. */
extern "C" int wrap_PyUFunc_getfperr() { return PyUFunc_getfperr(); }

namespace xsf {
namespace numpy {

    void set_error_check_fpe(const char *func_name) {
	int status = wrap_PyUFunc_getfperr();
	if (status & NPY_FPE_DIVIDEBYZERO) {
	    xsf::set_error(func_name, SF_ERROR_SINGULAR, "floating point division by zero");
	}
	if (status & NPY_FPE_OVERFLOW) {
	    xsf::set_error(func_name, SF_ERROR_UNDERFLOW, "floating point underflow");
	}
	if (status & NPY_FPE_UNDERFLOW) {
	    xsf::set_error(func_name, SF_ERROR_OVERFLOW, "floating point overflow");
	}
	if (status & NPY_FPE_INVALID) {
	    xsf::set_error(func_name, SF_ERROR_DOMAIN, "floating point invalid value");
	}
    }

    namespace detail {

        // This is std::accumulate, but that is not constexpr until C++20
        template <typename InputIt, typename T>
        constexpr T initializer_accumulate(InputIt first, InputIt last, T init) {
            for (InputIt it = first; it != last; ++it) {
                init = std::move(init) + *it;
            }

            return init;
        }

    } // namespace detail

    using cfloat = std::complex<float>;
    using cdouble = std::complex<double>;

    using float_1d = std::mdspan<float, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using float_2d = std::mdspan<float, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using float_3d = std::mdspan<float, std::dextents<ptrdiff_t, 3>, std::layout_stride>;
    using float_4d = std::mdspan<float, std::dextents<ptrdiff_t, 4>, std::layout_stride>;
    using double_1d = std::mdspan<double, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using double_2d = std::mdspan<double, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using double_3d = std::mdspan<double, std::dextents<ptrdiff_t, 3>, std::layout_stride>;
    using double_4d = std::mdspan<double, std::dextents<ptrdiff_t, 4>, std::layout_stride>;
    using cfloat_1d = std::mdspan<cfloat, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using cfloat_2d = std::mdspan<cfloat, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using cfloat_3d = std::mdspan<cfloat, std::dextents<ptrdiff_t, 3>, std::layout_stride>;
    using cfloat_4d = std::mdspan<cfloat, std::dextents<ptrdiff_t, 4>, std::layout_stride>;
    using cdouble_1d = std::mdspan<cdouble, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using cdouble_2d = std::mdspan<cdouble, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using cdouble_3d = std::mdspan<cdouble, std::dextents<ptrdiff_t, 3>, std::layout_stride>;
    using cdouble_4d = std::mdspan<cdouble, std::dextents<ptrdiff_t, 4>, std::layout_stride>;

    using autodiff0_float = dual<float, 0>;
    using autodiff0_double = dual<double, 0>;
    using autodiff0_cfloat = dual<cfloat, 0>;
    using autodiff0_cdouble = dual<cdouble, 0>;
    using autodiff0_float_1d = std::mdspan<autodiff0_float, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff0_double_1d = std::mdspan<autodiff0_double, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff0_cfloat_1d = std::mdspan<autodiff0_cfloat, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff0_cdouble_1d = std::mdspan<autodiff0_cdouble, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff0_float_2d = std::mdspan<autodiff0_float, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff0_double_2d = std::mdspan<autodiff0_double, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff0_cfloat_2d = std::mdspan<autodiff0_cfloat, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff0_cdouble_2d = std::mdspan<autodiff0_cdouble, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff1_float = dual<float, 1>;
    using autodiff1_double = dual<double, 1>;
    using autodiff1_cfloat = dual<cfloat, 1>;
    using autodiff1_cdouble = dual<cdouble, 1>;
    using autodiff1_float_1d = std::mdspan<autodiff1_float, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff1_double_1d = std::mdspan<autodiff1_double, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff1_cfloat_1d = std::mdspan<autodiff1_cfloat, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff1_cdouble_1d = std::mdspan<autodiff1_cdouble, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff1_float_2d = std::mdspan<autodiff1_float, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff1_double_2d = std::mdspan<autodiff1_double, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff1_cfloat_2d = std::mdspan<autodiff1_cfloat, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff1_cdouble_2d = std::mdspan<autodiff1_cdouble, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff2_float = dual<float, 2>;
    using autodiff2_double = dual<double, 2>;
    using autodiff2_cfloat = dual<cfloat, 2>;
    using autodiff2_cdouble = dual<cdouble, 2>;
    using autodiff2_float_1d = std::mdspan<autodiff2_float, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff2_double_1d = std::mdspan<autodiff2_double, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff2_cfloat_1d = std::mdspan<autodiff2_cfloat, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff2_cdouble_1d = std::mdspan<autodiff2_cdouble, std::dextents<ptrdiff_t, 1>, std::layout_stride>;
    using autodiff2_float_2d = std::mdspan<autodiff2_float, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff2_double_2d = std::mdspan<autodiff2_double, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff2_cfloat_2d = std::mdspan<autodiff2_cfloat, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff2_cdouble_2d = std::mdspan<autodiff2_cdouble, std::dextents<ptrdiff_t, 2>, std::layout_stride>;

    using autodiff00_float = dual<float, 0, 0>;
    using autodiff00_double = dual<double, 0, 0>;
    using autodiff00_cfloat = dual<cfloat, 0, 0>;
    using autodiff00_cdouble = dual<cdouble, 0, 0>;
    using autodiff00_cfloat_2d = std::mdspan<autodiff00_cfloat, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff00_cdouble_2d = std::mdspan<autodiff00_cdouble, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff11_float = dual<float, 1, 1>;
    using autodiff11_double = dual<double, 1, 1>;
    using autodiff11_cfloat = dual<cfloat, 1, 1>;
    using autodiff11_cdouble = dual<cdouble, 1, 1>;
    using autodiff11_cfloat_2d = std::mdspan<autodiff11_cfloat, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff11_cdouble_2d = std::mdspan<autodiff11_cdouble, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff22_float = dual<float, 2, 2>;
    using autodiff22_double = dual<double, 2, 2>;
    using autodiff22_cfloat = dual<cfloat, 2, 2>;
    using autodiff22_cdouble = dual<cdouble, 2, 2>;
    using autodiff22_cfloat_2d = std::mdspan<autodiff22_cfloat, std::dextents<ptrdiff_t, 2>, std::layout_stride>;
    using autodiff22_cdouble_2d = std::mdspan<autodiff22_cdouble, std::dextents<ptrdiff_t, 2>, std::layout_stride>;

    // The following are based off NumPy's dtype type codes and functions like PyUFunc_dd_d

    // 1 input, 1 output
    using f_f = float (*)(float);
    using d_d = double (*)(double);
    using F_F = cfloat (*)(cfloat);
    using D_D = cdouble (*)(cdouble);

    // autodiff, 1 input, 1 output
    using autodiff0_f_f1 = void (*)(autodiff0_float, autodiff0_float_1d);
    using autodiff0_d_d1 = void (*)(autodiff0_double, autodiff0_double_1d);
    using autodiff0_F_F1 = void (*)(autodiff0_cfloat, autodiff0_cfloat_1d);
    using autodiff0_D_D1 = void (*)(autodiff0_cdouble, autodiff0_cdouble_1d);
    using autodiff1_f_f1 = void (*)(autodiff1_float, autodiff1_float_1d);
    using autodiff1_d_d1 = void (*)(autodiff1_double, autodiff1_double_1d);
    using autodiff1_F_F1 = void (*)(autodiff1_cfloat, autodiff1_cfloat_1d);
    using autodiff1_D_D1 = void (*)(autodiff1_cdouble, autodiff1_cdouble_1d);
    using autodiff2_f_f1 = void (*)(autodiff2_float, autodiff2_float_1d);
    using autodiff2_d_d1 = void (*)(autodiff2_double, autodiff2_double_1d);
    using autodiff2_F_F1 = void (*)(autodiff2_cfloat, autodiff2_cfloat_1d);
    using autodiff2_D_D1 = void (*)(autodiff2_cdouble, autodiff2_cdouble_1d);

    using autodiff0_f_f2 = void (*)(autodiff0_float, autodiff0_float_2d);
    using autodiff0_d_d2 = void (*)(autodiff0_double, autodiff0_double_2d);
    using autodiff0_F_F2 = void (*)(autodiff0_cfloat, autodiff0_cfloat_2d);
    using autodiff0_D_D2 = void (*)(autodiff0_cdouble, autodiff0_cdouble_2d);
    using autodiff1_f_f2 = void (*)(autodiff1_float, autodiff1_float_2d);
    using autodiff1_d_d2 = void (*)(autodiff1_double, autodiff1_double_2d);
    using autodiff1_F_F2 = void (*)(autodiff1_cfloat, autodiff1_cfloat_2d);
    using autodiff1_D_D2 = void (*)(autodiff1_cdouble, autodiff1_cdouble_2d);
    using autodiff2_f_f2 = void (*)(autodiff2_float, autodiff2_float_2d);
    using autodiff2_d_d2 = void (*)(autodiff2_double, autodiff2_double_2d);
    using autodiff2_F_F2 = void (*)(autodiff2_cfloat, autodiff2_cfloat_2d);
    using autodiff2_D_D2 = void (*)(autodiff2_cdouble, autodiff2_cdouble_2d);

    // 1 input, 2 outputs
    using f_ff = void (*)(float, float &, float &);
    using d_dd = void (*)(double, double &, double &);
    using f_FF = void (*)(float, cfloat &, cfloat &);
    using d_DD = void (*)(double, cdouble &, cdouble &);
    using F_FF = void (*)(cfloat, cfloat &, cfloat &);
    using D_DD = void (*)(cdouble, cdouble &, cdouble &);

    // 1 input, 4 outputs
    using f_ffff = void (*)(float, float &, float &, float &, float &);
    using d_dddd = void (*)(double, double &, double &, double &, double &);
    using f_FFFF = void (*)(float, cfloat &, cfloat &, cfloat &, cfloat &);
    using d_DDDD = void (*)(double, cdouble &, cdouble &, cdouble &, cdouble &);
    using F_FFFF = void (*)(cfloat, cfloat &, cfloat &, cfloat &, cfloat &);
    using D_DDDD = void (*)(cdouble, cdouble &, cdouble &, cdouble &, cdouble &);

    // 2 inputs, 1 output
    using qf_f = float (*)(long long int, float);
    using qd_d = double (*)(long long int, double);
    using ff_f = float (*)(float, float);
    using dd_d = double (*)(double, double);
    using FF_F = cfloat (*)(cfloat, cfloat);
    using DD_D = cdouble (*)(cdouble, cdouble);
    using fF_F = cfloat (*)(float, cfloat);
    using dD_D = cdouble (*)(double, cdouble);
    using lf_f = float (*)(long int, float);
    using ld_d = double (*)(long int, double);
    using lF_F = cfloat (*)(long int, cfloat);
    using lD_D = cdouble (*)(long int, cdouble);
    using Dd_D = cdouble (*) (cdouble, double);
    using Ff_F = cfloat (*) (cfloat, float);

    // autodiff, 2 inputs, 1 output
    using autodiff0_if_f = autodiff0_float (*)(int, autodiff0_float);
    using autodiff0_id_d = autodiff0_double (*)(int, autodiff0_double);
    using autodiff1_if_f = autodiff1_float (*)(int, autodiff1_float);
    using autodiff1_id_d = autodiff1_double (*)(int, autodiff1_double);
    using autodiff2_if_f = autodiff2_float (*)(int, autodiff2_float);
    using autodiff2_id_d = autodiff2_double (*)(int, autodiff2_double);

    // 2 inputs, 2 outputs
    using qf_ff = void (*)(long long int, float, float &, float &);
    using qd_dd = void (*)(long long int, double, double &, double &);

    // 2 inputs, 3 outputs
    using qf_fff = void (*)(long long int, float, float &, float &, float &);
    using qd_ddd = void (*)(long long int, double, double &, double &, double &);

    // 2 inputs, 2 outputs
    using ff_ff = void (*)(float, float, float &, float &);
    using dd_dd = void (*)(double, double, double &, double &);
    using lf_ff = void (*)(long int, float, float &, float &);
    using ld_dd = void (*)(long int, double, double &, double &);

    // 2 inputs, 3 outputs
    using lf_fff = void (*)(long int, float, float &, float &, float &);
    using ld_ddd = void (*)(long int, double, double &, double &, double &);

    // 2 inputs, 4 outputs
    using ff_ffff = void (*)(float, float, float &, float &, float &, float &);
    using dd_dddd = void (*)(double, double, double &, double &, double &, double &);

    // 3 inputs, 1 output
    using fff_f = float (*)(float, float, float);
    using ddd_d = double (*)(double, double, double);
    using Flf_F = cfloat (*)(cfloat, long int, float);
    using Dld_D = cdouble (*)(cdouble, long int, double);

    // 3 inputs, 2 outputs
    using fff_ff = void (*)(float, float, float, float &, float &);
    using ddd_dd = void (*)(double, double, double, double &, double &);

    // 3 inputs, 1 output
    using qqf_f = float (*)(long long int, long long int, float);
    using qqd_d = double (*)(long long int, long long int, double);

    // autodiff, 3 inputs, 1 ouput
    using autodiff0_iif_f = autodiff0_float (*)(int, int, autodiff0_float);
    using autodiff0_iid_d = autodiff0_double (*)(int, int, autodiff0_double);
    using autodiff1_iif_f = autodiff1_float (*)(int, int, autodiff1_float);
    using autodiff1_iid_d = autodiff1_double (*)(int, int, autodiff1_double);
    using autodiff2_iif_f = autodiff2_float (*)(int, int, autodiff2_float);
    using autodiff2_iid_d = autodiff2_double (*)(int, int, autodiff2_double);

    // 3 inputs, 2 outputs
    using qqf_ff = void (*)(long long int, long long int, float, float &, float &);
    using qqd_dd = void (*)(long long int, long long int, double, double &, double &);

    // 3 inputs, 3 outputs
    using qqf_fff = void (*)(long long int, long long int, float, float &, float &, float &);
    using qqd_ddd = void (*)(long long int, long long int, double, double &, double &, double &);

    // 4 inputs, 1 output
    using qqqF_F = cfloat (*)(long long int, long long int, long long int, cfloat);
    using qqqD_D = cdouble (*)(long long int, long long int, long long int, cdouble);
    using qqff_F = cfloat (*)(long long int, long long int, float, float);
    using qqdd_D = cdouble (*)(long long int, long long int, double, double);
    using ffff_f = float (*)(float, float, float, float);
    using dddd_d = double (*)(double, double, double, double);
    using fffF_F = cfloat (*)(float, float, float, cfloat);
    using dddD_D = cdouble (*)(double, double, double, cdouble);
    using ffff_F = cfloat (*)(float, float, float, float);
    using dddd_D = cdouble (*)(double, double, double, double);

    // autodiff, 4 inputs, 1 output
    using autodiff00_iiff_F = autodiff00_cfloat (*)(int, int, autodiff00_float, autodiff00_float);
    using autodiff00_iidd_D = autodiff00_cdouble (*)(int, int, autodiff00_double, autodiff00_double);
    using autodiff11_iiff_F = autodiff11_cfloat (*)(int, int, autodiff11_float, autodiff11_float);
    using autodiff11_iidd_D = autodiff11_cdouble (*)(int, int, autodiff11_double, autodiff11_double);
    using autodiff22_iiff_F = autodiff22_cfloat (*)(int, int, autodiff22_float, autodiff22_float);
    using autodiff22_iidd_D = autodiff22_cdouble (*)(int, int, autodiff22_double, autodiff22_double);

    // 4 inputs, 2 outputs
    using qqqf_ff = void (*)(long long int, long long int, long long int, float, float &, float &);
    using ffff_ff = void (*)(float, float, float, float, float &, float &);
    using qqqd_dd = void (*)(long long int, long long int, long long int, double, double &, double &);
    using dddd_dd = void (*)(double, double, double, double, double &, double &);
    using qqqF_FF = void (*)(long long int, long long int, long long int, cfloat, cfloat &, cfloat &);
    using qqqD_DD = void (*)(long long int, long long int, long long int, cdouble, cdouble &, cdouble &);
    using qqff_FF = void (*)(long long int, long long int, float, float, cfloat &, cfloat &);
    using qqdd_DD = void (*)(long long int, long long int, double, double, cdouble &, cdouble &);

    // 4 inputs, 3 outputs
    using qqqf_fff = void (*)(long long int, long long int, long long int, float, float &, float &, float &);
    using qqqd_ddd = void (*)(long long int, long long int, long long int, double, double &, double &, double &);
    using qqqF_FFF = void (*)(long long int, long long int, long long int, cfloat, cfloat &, cfloat &, cfloat &);
    using qqqD_DDD = void (*)(long long int, long long int, long long int, cdouble, cdouble &, cdouble &, cdouble &);
    using qqff_FFF = void (*)(long long int, long long int, float, float, cfloat &, cfloat &, cfloat &);
    using qqdd_DDD = void (*)(long long int, long long int, double, double, cdouble &, cdouble &, cdouble &);

    // 5 inputs, 2 outputs
    using fffff_ff = void (*)(float, float, float, float, float, float &, float &);
    using ddddd_dd = void (*)(double, double, double, double, double, double &, double &);

#if (NPY_SIZEOF_LONGDOUBLE == NPY_SIZEOF_DOUBLE)
    using g_g = double (*)(double);
    using gg_g = double (*)(double);
#else
    using g_g = long double (*)(long double);
    using gg_g = long double (*)(long double);
#endif

    // 1 input, 1 output
    using f_f1 = void (*)(float, float_1d);
    using f_f2 = void (*)(float, float_2d);
    using d_d1 = void (*)(double, double_1d);
    using d_d2 = void (*)(double, double_2d);
    using F_F1 = void (*)(cfloat, cfloat_1d);
    using D_D1 = void (*)(cdouble, cdouble_1d);

    // 1 input, 2 outputs
    using f_f1f1 = void (*)(float, float_1d, float_1d);
    using f_f2f2 = void (*)(float, float_2d, float_2d);
    using d_d1d1 = void (*)(double, double_1d, double_1d);
    using d_d2d2 = void (*)(double, double_2d, double_2d);
    using F_F1F1 = void (*)(cfloat, cfloat_1d, cfloat_1d);
    using F_F2F2 = void (*)(cfloat, cfloat_2d, cfloat_2d);
    using D_D1D1 = void (*)(cdouble, cdouble_1d, cdouble_1d);
    using D_D2D2 = void (*)(cdouble, cdouble_2d, cdouble_2d);

    // 1 input, 3 outputs
    using f_f1f1f1 = void (*)(float, float_1d, float_1d, float_1d);
    using f_f2f2f2 = void (*)(float, float_2d, float_2d, float_2d);
    using d_d1d1d1 = void (*)(double, double_1d, double_1d, double_1d);
    using d_d2d2d2 = void (*)(double, double_2d, double_2d, double_2d);
    using F_F1F1F1 = void (*)(cfloat, cfloat_1d, cfloat_1d, cfloat_1d);
    using D_D1D1D1 = void (*)(cdouble, cdouble_1d, cdouble_1d, cdouble_1d);

    // 2 inputs, 1 output
    using ff_F2 = void (*)(float, float, cfloat_2d);
    using dd_D2 = void (*)(double, double, cdouble_2d);
    using qF_F2 = void (*)(long long int, cfloat, cfloat_2d);
    using qD_D2 = void (*)(long long int, cdouble, cdouble_2d);

    using autodiff00_ff_F2 = void (*)(autodiff00_float, autodiff00_float, autodiff00_cfloat_2d);
    using autodiff00_dd_D2 = void (*)(autodiff00_double, autodiff00_double, autodiff00_cdouble_2d);
    using autodiff11_ff_F2 = void (*)(autodiff11_float, autodiff11_float, autodiff11_cfloat_2d);
    using autodiff11_dd_D2 = void (*)(autodiff11_double, autodiff11_double, autodiff11_cdouble_2d);
    using autodiff22_ff_F2 = void (*)(autodiff22_float, autodiff22_float, autodiff22_cfloat_2d);
    using autodiff22_dd_D2 = void (*)(autodiff22_double, autodiff22_double, autodiff22_cdouble_2d);

    // 2 inputs, 2 outputs
    using qF_F2F2 = void (*)(long long int, cfloat, cfloat_2d, cfloat_2d);
    using qD_D2D2 = void (*)(long long int, cdouble, cdouble_2d, cdouble_2d);
    using ff_F2F3 = void (*)(float, float, cfloat_2d, cfloat_3d);
    using dd_D2D3 = void (*)(double, double, cdouble_2d, cdouble_3d);

    // 2 inputs, 3 outputs
    using qF_F2F2F2 = void (*)(long long int, cfloat, cfloat_2d, cfloat_2d, cfloat_2d);
    using qD_D2D2D2 = void (*)(long long int, cdouble, cdouble_2d, cdouble_2d, cdouble_2d);

    // 2 inputs, 4 outputs
    using ff_F2F3F4 = void (*)(float, float, cfloat_2d, cfloat_3d, cfloat_4d);
    using dd_D2D3D4 = void (*)(double, double, cdouble_2d, cdouble_3d, cdouble_4d);

    template <typename Func>
    struct signature_of {
        using type = typename signature_of<decltype(&Func::operator())>::type;
    };

    template <typename Res, typename... Args>
    struct signature_of<Res(Args...)> {
        using type = Res(Args...);
    };

    template <typename Res, typename... Args>
    struct signature_of<Res (*)(Args...)> {
        using type = Res(Args...);
    };

    template <typename T, typename Res, typename... Args>
    struct signature_of<Res (T::*)(Args...)> {
        using type = Res(Args...);
    };

    template <typename T, typename Res, typename... Args>
    struct signature_of<Res (T::*)(Args...) const> {
        using type = Res(Args...);
    };

    template <typename Func>
    using signature_of_t = typename signature_of<Func>::type;

    // Deduces the number of arguments of a callable F.
    template <typename Func>
    struct arity_of {
        static constexpr size_t value = arity_of<signature_of_t<Func>>::value;
    };

    template <typename Res, typename... Args>
    struct arity_of<Res(Args...)> {
        static constexpr size_t value = sizeof...(Args);
    };

    template <typename Func>
    constexpr size_t arity_of_v = arity_of<Func>::value;

    template <typename Func>
    struct has_return {
        static constexpr bool value = has_return<signature_of_t<Func>>::value;
    };

    template <typename Res, typename... Args>
    struct has_return<Res(Args...)> {
        static constexpr bool value = true;
    };

    template <typename... Args>
    struct has_return<void(Args...)> {
        static constexpr bool value = false;
    };

    template <typename Func>
    constexpr size_t has_return_v = has_return<Func>::value;

    template <typename T>
    struct rank_of {
        static constexpr size_t value = 0;
    };

    template <typename T, typename Extents, typename LayoutPolicy, typename AccessorPolicy>
    struct rank_of<std::mdspan<T, Extents, LayoutPolicy, AccessorPolicy>> {
        static constexpr size_t value = Extents::rank() + rank_of<T>::value;
    };

    template <typename T, size_t... Orders>
    struct rank_of<dual<T, Orders...>> {
        static constexpr size_t value = sizeof...(Orders);
    };

    template <typename T>
    inline constexpr size_t rank_of_v = rank_of<T>::value;

    // Maps a C++ type to a NumPy type
    template <typename T>
    struct npy_type;

    template <>
    struct npy_type<bool> {
        using type = npy_bool;
    };

    template <>
    struct npy_type<char> {
        using type = npy_byte;
    };

    template <>
    struct npy_type<short> {
        using type = npy_short;
    };

    template <>
    struct npy_type<int> {
        using type = npy_int;
    };

    template <>
    struct npy_type<long> {
        using type = npy_long;
    };

    template <>
    struct npy_type<long long> {
        using type = npy_longlong;
    };

    template <>
    struct npy_type<unsigned char> {
        using type = npy_ubyte;
    };

    template <>
    struct npy_type<unsigned short> {
        using type = npy_ushort;
    };

    template <>
    struct npy_type<unsigned int> {
        using type = npy_uint;
    };

    template <>
    struct npy_type<unsigned long> {
        using type = npy_ulong;
    };

    template <>
    struct npy_type<unsigned long long> {
        using type = npy_ulonglong;
    };

    template <>
    struct npy_type<float> {
        using type = npy_float;
    };

    template <>
    struct npy_type<double> {
        using type = npy_double;
    };

    template <>
    struct npy_type<long double> {
        using type = npy_longdouble;
    };

    template <>
    struct npy_type<std::complex<float>> {
        using type = npy_cfloat;
    };

    template <>
    struct npy_type<std::complex<double>> {
        using type = npy_cdouble;
    };

    template <>
    struct npy_type<std::complex<long double>> {
        using type = npy_clongdouble;
    };

    template <typename T>
    using npy_type_t = typename npy_type<T>::type;

    // Maps a C++ type to a NumPy type number
    template <typename T>
    struct npy_typenum {
        static constexpr NPY_TYPES value = npy_typenum<npy_type_t<T>>::value;
    };

    // We need to specialise for bool as npy_bool is defined as npy_ubyte
    template <>
    struct npy_typenum<bool> {
        static constexpr NPY_TYPES value = NPY_BOOL;
    };

    template <>
    struct npy_typenum<npy_byte> {
        static constexpr NPY_TYPES value = NPY_BYTE;
    };

    template <>
    struct npy_typenum<npy_short> {
        static constexpr NPY_TYPES value = NPY_SHORT;
    };

    template <>
    struct npy_typenum<npy_int> {
        static constexpr NPY_TYPES value = NPY_INT;
    };

    template <>
    struct npy_typenum<npy_long> {
        static constexpr NPY_TYPES value = NPY_LONG;
    };

    template <>
    struct npy_typenum<npy_longlong> {
        static constexpr NPY_TYPES value = NPY_LONGLONG;
    };

    template <>
    struct npy_typenum<npy_ubyte> {
        static constexpr NPY_TYPES value = NPY_UBYTE;
    };

    template <>
    struct npy_typenum<npy_ushort> {
        static constexpr NPY_TYPES value = NPY_USHORT;
    };

    template <>
    struct npy_typenum<npy_uint> {
        static constexpr NPY_TYPES value = NPY_UINT;
    };

    template <>
    struct npy_typenum<npy_ulong> {
        static constexpr NPY_TYPES value = NPY_ULONG;
    };

    template <>
    struct npy_typenum<npy_ulonglong> {
        static constexpr NPY_TYPES value = NPY_ULONGLONG;
    };

    template <>
    struct npy_typenum<npy_float> {
        static constexpr NPY_TYPES value = NPY_FLOAT;
    };

    template <>
    struct npy_typenum<npy_double> {
        static constexpr NPY_TYPES value = NPY_DOUBLE;
    };

// When NPY_SIZEOF_LONGDOUBLE == NPY_SIZEOF_DOUBLE, npy_longdouble is defined as npy_double
// See https://github.com/numpy/numpy/blob/main/numpy/_core/include/numpy/npy_common.h
#if (NPY_SIZEOF_LONGDOUBLE != NPY_SIZEOF_DOUBLE)
    template <>
    struct npy_typenum<npy_longdouble> {
        static constexpr NPY_TYPES value = NPY_LONGDOUBLE;
    };
#endif

    template <>
    struct npy_typenum<npy_cfloat> {
        static constexpr NPY_TYPES value = NPY_CFLOAT;
    };

    template <>
    struct npy_typenum<npy_cdouble> {
        static constexpr NPY_TYPES value = NPY_CDOUBLE;
    };

    template <>
    struct npy_typenum<npy_clongdouble> {
        static constexpr NPY_TYPES value = NPY_CLONGDOUBLE;
    };

    template <typename T>
    struct npy_typenum<T *> {
        static constexpr NPY_TYPES value = npy_typenum<T>::value;
    };

    template <typename T>
    struct npy_typenum<T &> {
        static constexpr NPY_TYPES value = npy_typenum<T>::value;
    };

    template <typename T, size_t N>
    struct npy_typenum<T[N]> {
        static constexpr NPY_TYPES value = npy_typenum<T>::value;
    };

    template <typename T, typename Extents, typename LayoutPolicy, typename AccessorPolicy>
    struct npy_typenum<std::mdspan<T, Extents, LayoutPolicy, AccessorPolicy>> {
        static constexpr NPY_TYPES value = npy_typenum<T>::value;
    };

    template <typename T, size_t... Orders>
    struct npy_typenum<dual<T, Orders...>> {
        static constexpr NPY_TYPES value = npy_typenum<T>::value;
    };

    template <typename T>
    inline constexpr NPY_TYPES npy_typenum_v = npy_typenum<T>::value;

    template <typename T>
    struct npy_traits {
        static T get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
            return *reinterpret_cast<T *>(src);
        }

        static void set(char *dst, const T &src) { *reinterpret_cast<npy_type_t<T> *>(dst) = src; }
    };

    template <typename T>
    struct npy_traits<std::complex<T>> {
        static std::complex<T> get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
            return *reinterpret_cast<std::complex<T> *>(src);
        }

        static void set(char *dst, const std::complex<T> &src) {
            *reinterpret_cast<npy_type_t<T> *>(dst) = std::real(src);
            *reinterpret_cast<npy_type_t<T> *>(dst + sizeof(T)) = std::imag(src);
        }
    };

    template <typename T>
    struct npy_traits<T *> {
        static T *get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
            static_assert(sizeof(T) == sizeof(npy_type_t<T>), "NumPy type has different size than argument type");

            return reinterpret_cast<T *>(src);
        }
    };

    template <typename T>
    struct npy_traits<T &> {
        static T &get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
            static_assert(sizeof(T) == sizeof(npy_type_t<T>), "NumPy type has different size than argument type");

            return *reinterpret_cast<T *>(src);
        }
    };

    template <typename T, size_t N>
    struct npy_traits<T (&)[N]> {
        using dst_type = T (&)[N];

        static dst_type get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
            static_assert(sizeof(T) == sizeof(npy_type_t<T>), "NumPy type has different size than argument type");

            return *reinterpret_cast<T(*)[N]>(src);
        }
    };

    template <typename T, size_t N>
    struct npy_traits<T (&)[N][N]> {
        using dst_type = T (&)[N][N];

        static dst_type get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
            static_assert(sizeof(T) == sizeof(npy_type_t<T>), "NumPy type has different size than argument type");

            return *reinterpret_cast<T(*)[N][N]>(src);
        }
    };

    template <typename T, typename Extents, typename AccessorPolicy>
    struct npy_traits<std::mdspan<T, Extents, std::layout_stride, AccessorPolicy>> {
        static std::mdspan<T, Extents, std::layout_stride, AccessorPolicy>
        get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
            //            static_assert(sizeof(T) == sizeof(npy_type_t<T>), "NumPy type has different size than argument
            //            type");

            std::array<ptrdiff_t, Extents::rank()> strides;
            for (npy_uintp i = 0; i < strides.size(); ++i) {
                strides[i] = steps[i] / sizeof(T);
            }

            std::array<ptrdiff_t, Extents::rank()> exts;
            for (npy_uintp i = 0; i < exts.size(); ++i) {
                exts[i] = dimensions[i];
            }

            return {reinterpret_cast<T *>(src), {exts, strides}};
        }
    };

    template <typename T, size_t... Orders>
    struct npy_traits<dual<T, Orders...>> {
        static dual<T, Orders...> get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
            return *reinterpret_cast<dual<T, Orders...> *>(src);
        }

        static void set(char *dst, const dual<T, Orders...> &src) {
            *reinterpret_cast<dual<T, Orders...> *>(dst) = src;
        }
    };

    using map_dims_type = void (*)(const npy_intp *, npy_intp *);

    struct base_ufunc_data {
        const char *name;
        map_dims_type map_dims;
        int flags;
    };

    template <typename Func>
    struct ufunc_data : base_ufunc_data {
        Func func;
    };

    template <
        typename Func, typename Signature = signature_of_t<Func>,
        typename Indices = std::make_index_sequence<arity_of_v<Signature>>>
    struct ufunc_traits;

    template <typename Func, typename Res, typename... Args, size_t... I>
    struct ufunc_traits<Func, Res(Args...), std::index_sequence<I...>> {
        static constexpr char types[sizeof...(Args) + 1] = {npy_typenum_v<Args>..., npy_typenum_v<Res>};

        static constexpr size_t ranks[sizeof...(Args) + 1] = {rank_of_v<Args>..., rank_of_v<Res>};

        static constexpr size_t ranks_scan[sizeof...(Args) + 2] = {
            detail::initializer_accumulate(ranks, ranks + I, 0)...,
            detail::initializer_accumulate(ranks, ranks + sizeof...(Args), 0),
            detail::initializer_accumulate(ranks, ranks + sizeof...(Args) + 1, 0)
        };

        static void loop(char **args, const npy_intp *dims, const npy_intp *steps, void *data) {
            std::array<npy_intp, ranks_scan[sizeof...(Args) + 1]> new_dims;

            map_dims_type map_dims = static_cast<ufunc_data<Func> *>(data)->map_dims;
            map_dims(dims + 1, new_dims.data());

            Func func = static_cast<ufunc_data<Func> *>(data)->func;
            for (npy_intp i = 0; i < dims[0]; ++i) {
                Res res = func(npy_traits<Args>::get(
                    args[I], new_dims.data() + ranks_scan[I], steps + ranks_scan[I] + sizeof...(Args) + 1
                )...);
                npy_traits<Res>::set(args[sizeof...(Args)], res); // assign to the output pointer

                for (npy_uintp j = 0; j <= sizeof...(Args); ++j) {
                    args[j] += steps[j];
                }
            }

            const char *name = static_cast<ufunc_data<Func> *>(data)->name;
	    set_error_check_fpe(name);
        }
    };

    template <typename Func, typename... Args, size_t... I>
    struct ufunc_traits<Func, void(Args...), std::index_sequence<I...>> {
        static constexpr char types[sizeof...(Args)] = {npy_typenum_v<Args>...};

        static constexpr size_t ranks[sizeof...(Args)] = {rank_of_v<Args>...};

        static constexpr size_t ranks_scan[sizeof...(Args) + 1] = {
            detail::initializer_accumulate(ranks, ranks + I, 0)...,
            detail::initializer_accumulate(ranks, ranks + sizeof...(Args), 0)
        };

        static void loop(char **args, const npy_intp *dims, const npy_intp *steps, void *data) {
            std::array<npy_intp, ranks_scan[sizeof...(Args)]> new_dims;

            map_dims_type map_dims = static_cast<ufunc_data<Func> *>(data)->map_dims;
            map_dims(dims + 1, new_dims.data());

            Func func = static_cast<ufunc_data<Func> *>(data)->func;
            for (npy_intp i = 0; i < dims[0]; ++i) {
                func(npy_traits<Args>::get(
                    args[I], new_dims.data() + ranks_scan[I], steps + ranks_scan[I] + sizeof...(Args)
                )...);

                for (npy_uintp j = 0; j < sizeof...(Args); ++j) {
                    args[j] += steps[j];
                }
            }

            const char *name = static_cast<ufunc_data<Func> *>(data)->name;
	    set_error_check_fpe(name);
        }
    };

    namespace detail {

        template <typename Func, typename Tr>
        decltype(auto) compose(Func func, Tr tr) {
            return tr(func);
        }

        template <typename Func, typename Tr0, typename Tr1, typename... Trs>
        decltype(auto) compose(Func func, Tr0 tr0, Tr1 tr1, Trs... trs) {
            return compose(tr0(func), tr1, trs...);
        }

    } // namespace detail

    template <typename... Trs>
    class compose {
        std::tuple<Trs...> m_trs;

      public:
        compose(Trs... trs) : m_trs(trs...) {}

        template <typename Func>
        decltype(auto) operator()(Func func) const {
            return std::apply([func](auto... trs) { return detail::compose(func, trs...); }, m_trs);
        }
    };

    struct ufunc_wraps {
        bool has_return;
        int nin_and_nout;
        PyUFuncGenericFunction func;
        void *data;
        void (*data_deleter)(void *);
        const char *types;

        template <typename Func>
        ufunc_wraps(Func func)
            : has_return(has_return_v<Func>), nin_and_nout(arity_of_v<Func> + has_return),
              func(ufunc_traits<Func>::loop), data(new ufunc_data<Func>{{nullptr}, func}),
              data_deleter([](void *ptr) { delete static_cast<ufunc_data<Func> *>(ptr); }),
              types(ufunc_traits<Func>::types) {}
    };

    class ufunc_overloads {
      public:
        using data_handle_type = void *;
        using data_deleter_type = void (*)(void *);

      private:
        int m_ntypes;
        bool m_has_return;
        int m_nin_and_nout;
        std::unique_ptr<PyUFuncGenericFunction[]> m_func;
        std::unique_ptr<data_handle_type[]> m_data;
        std::unique_ptr<data_deleter_type[]> m_data_deleters;
        std::unique_ptr<char[]> m_types;

      public:
        template <typename Func0, typename... Funcs>
        ufunc_overloads(Func0 func0, Funcs... funcs)
            : m_ntypes(sizeof...(Funcs) + 1), m_has_return(has_return_v<Func0>),
              m_nin_and_nout(arity_of_v<Func0> + m_has_return), m_func(new PyUFuncGenericFunction[m_ntypes]),
              m_data(new data_handle_type[m_ntypes]), m_data_deleters(new data_deleter_type[m_ntypes]),
              m_types(new char[m_ntypes * m_nin_and_nout]) {
            ufunc_wraps func[sizeof...(Funcs) + 1] = {func0, funcs...};
            for (auto it = std::begin(func); it != std::end(func); ++it) {
                if (it->nin_and_nout != m_nin_and_nout) {
                    PyErr_SetString(PyExc_RuntimeError, "all functions must have the same number of arguments");
                }
                if (it->has_return != m_has_return) {
                    PyErr_SetString(PyExc_RuntimeError, "all functions must be void if any function is");
                }

                size_t i = it - std::begin(func);
                m_func[i] = it->func;
                m_data[i] = it->data;
                m_data_deleters[i] = it->data_deleter;
                std::memcpy(m_types.get() + i * m_nin_and_nout, it->types, m_nin_and_nout);
            }
        }

        template <typename... Trs, typename... Funcs>
        ufunc_overloads(compose<Trs...> trs, Funcs... funcs) : ufunc_overloads(trs(funcs)...) {}

        ufunc_overloads(ufunc_overloads &&other) = default;

        ~ufunc_overloads() {
            if (m_data) {
                for (int i = 0; i < m_ntypes; ++i) {
                    data_deleter_type data_deleter = m_data_deleters[i];
                    data_deleter(m_data[i]);
                }
            }
        }

        int ntypes() const { return m_ntypes; }

        bool has_return() const { return m_has_return; }

        int nin_and_nout() const { return m_nin_and_nout; }

        PyUFuncGenericFunction *func() const { return m_func.get(); }

        data_handle_type *data() const { return m_data.get(); }

        char *types() const { return m_types.get(); }

        void set_name(const char *name) {
            for (int i = 0; i < m_ntypes; ++i) {
                static_cast<base_ufunc_data *>(m_data[i])->name = name;
            }
        }

        void set_map_dims(map_dims_type map_dims) {
            for (int i = 0; i < m_ntypes; ++i) {
                static_cast<base_ufunc_data *>(m_data[i])->map_dims = map_dims;
            }
        }
    };

    PyObject *ufunc(ufunc_overloads func, int nout, const char *name, const char *doc) {
        static std::vector<ufunc_overloads> ufuncs;

        if (PyErr_Occurred()) {
            return nullptr;
        }

        ufunc_overloads &ufunc = ufuncs.emplace_back(std::move(func));
        ufunc.set_name(name);
        ufunc.set_map_dims([](const npy_intp *dims, npy_intp *new_dims) {});

        return PyUFunc_FromFuncAndData(
            ufunc.func(), ufunc.data(), ufunc.types(), ufunc.ntypes(), ufunc.nin_and_nout() - nout, nout, PyUFunc_None,
            name, doc, 0
        );
    }

    PyObject *ufunc(ufunc_overloads overloads, const char *name, const char *doc) {
        int nout = overloads.has_return();

        return ufunc(std::move(overloads), nout, name, doc);
    }

    PyObject *gufunc(
        ufunc_overloads overloads, int nout, const char *name, const char *doc, const char *signature,
        map_dims_type map_dims
    ) {
        static std::vector<ufunc_overloads> ufuncs;

        if (PyErr_Occurred()) {
            return nullptr;
        }

        ufunc_overloads &ufunc = ufuncs.emplace_back(std::move(overloads));
        ufunc.set_name(name);
        ufunc.set_map_dims(map_dims);

        return PyUFunc_FromFuncAndDataAndSignature(
            ufunc.func(), ufunc.data(), ufunc.types(), ufunc.ntypes(), ufunc.nin_and_nout() - nout, nout, PyUFunc_None,
            name, doc, 0, signature
        );
    }

    PyObject *gufunc(
        ufunc_overloads overloads, const char *name, const char *doc, const char *signature, map_dims_type map_dims
    ) {
        int nout = overloads.has_return();

        return gufunc(std::move(overloads), nout, name, doc, signature, map_dims);
    }

    // rename to autodiff_var?
    template <typename T>
    struct autodiff_traits {
        static T to_var(T arg, size_t i) { return arg; }
    };

    template <typename T, size_t... Orders>
    struct autodiff_traits<dual<T, Orders...>> {
        static dual<T, Orders...> to_var(T arg, size_t i) { return dual_var<Orders...>(arg, i); }
    };

    template <
        typename Func, typename Signature = signature_of_t<Func>,
        typename Indices = std::make_index_sequence<arity_of_v<Signature>>>
    struct autodiff_wrapper;

    template <typename Func, typename Res, typename... Args, size_t... I>
    struct autodiff_wrapper<Func, Res(Args...), std::index_sequence<I...>> {
        Func func;

        Res operator()(remove_dual_t<Args>... args) {
            return func(autodiff_traits<Args>::to_var(args, I - i_scan[I])...);
        }

        static constexpr size_t is_autodiff[sizeof...(Args)] = {std::is_same_v<Args, remove_dual_t<Args>>...};

        static constexpr size_t i_scan[sizeof...(Args)] = {
            detail::initializer_accumulate(is_autodiff, is_autodiff + I, 0)...,
        };
    };

    template <typename Func>
    autodiff_wrapper(Func func) -> autodiff_wrapper<Func>;

    struct autodiff {
        template <typename Func>
        decltype(auto) operator()(Func f) {
            return autodiff_wrapper{f};
        }
    };

    template <typename Func, typename Signature = signature_of_t<Func>>
    struct use_long_long_int_wrapper;

    template <typename Func, typename Res, typename... Args>
    struct use_long_long_int_wrapper<Func, Res(Args...)> {
        Func func;

        Res operator()(
            std::conditional_t<std::is_integral_v<Args> && !std::is_same_v<Args, bool>, long long int, Args>... args
        ) {
            return func(args...);
        }
    };

    template <typename Func>
    use_long_long_int_wrapper(Func func) -> use_long_long_int_wrapper<Func>;

    struct use_long_long_int {
        template <typename Func>
        decltype(auto) operator()(Func f) {
            return use_long_long_int_wrapper{f};
        }
    };

} // namespace numpy
} // namespace xsf
