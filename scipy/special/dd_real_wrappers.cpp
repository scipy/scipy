extern "C" {
#include "dd_real_wrappers.h"
}

#include "special/cephes/dd_real.h"


extern "C" double2 dd_create_d(double x) {
    return {x, 0.0};
}

extern "C" double2 dd_add(const double2* a, const double2* b) {
    special::dd_real::DoubleDouble dd_a(a->hi, a->lo);
    special::dd_real::DoubleDouble dd_b(b->hi, b->lo);
    special::dd_real::DoubleDouble result = dd_a + dd_b;
    return {result.hi, result.lo};
}

extern "C" double2 dd_mul(const double2* a, const double2* b) {
    special::dd_real::DoubleDouble dd_a(a->hi, a->lo);
    special::dd_real::DoubleDouble dd_b(b->hi, b->lo);
    special::dd_real::DoubleDouble result = dd_a * dd_b;
    return {result.hi, result.lo};
}

extern "C" double2 dd_div(const double2* a, const double2* b) {
    special::dd_real::DoubleDouble dd_a(a->hi, a->lo);
    special::dd_real::DoubleDouble dd_b(b->hi, b->lo);
    special::dd_real::DoubleDouble result = dd_a / dd_b;
    return {result.hi, result.lo};
}

extern "C" double dd_to_double(const double2* a) {
    return a->hi;
}
