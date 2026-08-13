#pragma once

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct double2 {
        double hi;
        double lo;
    } double2;

    double2 dd_create_d(double x);
    double2 dd_create(double x, double y);
    double2 dd_add(const double2* a, const double2* b);
    double2 dd_mul(const double2* a, const double2* b);
    double2 dd_div(const double2* a, const double2* b);
    double2 dd_exp(const double2* x);
    double2 dd_log(const double2* x);
    double dd_to_double(const double2* a);

#ifdef __cplusplus
}
#endif
