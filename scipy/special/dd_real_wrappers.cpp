/* These wrappers exist so that double-double extended precision arithmetic
 * now translated to in C++ in xsf/cephes/dd_real.h can be used in Cython.
 *  The original API of the C implementation which existed prior to gh-20390 has
 *  been replicated to avoid the need to modify downstream Cython files.
 */

#include "dd_real_wrappers.h"
#include <xsf/cephes/dd_real.h>

using xsf::cephes::detail::double_double;

double2 dd_create_d(double x) {
    return {x, 0.0};
}

double2 dd_create(double x, double y) {
    return {x, y};
}

double2 dd_add(const double2* a, const double2* b) {
    double_double dd_a(a->hi, a->lo);
    double_double dd_b(b->hi, b->lo);
    double_double result = dd_a + dd_b;
    return {result.hi, result.lo};
}

double2 dd_mul(const double2* a, const double2* b) {
    double_double dd_a(a->hi, a->lo);
    double_double dd_b(b->hi, b->lo);
    double_double result = dd_a * dd_b;
    return {result.hi, result.lo};
}

double2 dd_div(const double2* a, const double2* b) {
    double_double dd_a(a->hi, a->lo);
    double_double dd_b(b->hi, b->lo);
    double_double result = dd_a / dd_b;
    return {result.hi, result.lo};
}

double2 dd_exp(const double2* x) {
    double_double dd_x(x->hi, x->lo);
    double_double result = xsf::cephes::detail::exp(dd_x);
    return {result.hi, result.lo};
}

double2 dd_log(const double2* x) {
    double_double dd_x(x->hi, x->lo);
    double_double result = xsf::cephes::detail::log(dd_x);
    return {result.hi, result.lo};
}

double dd_to_double(const double2* a) {
    return a->hi;
}


