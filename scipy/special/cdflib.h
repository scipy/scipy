#ifndef CDFLIB_H
#define CDFLIB_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include <math.h>

struct TupleDD
{
    double d1;
    double d2;
};

struct TupleDI
{
    double d1;
    int i1;
};

struct TupleDDI
{
    double d1;
    double d2;
    int i1;
};

struct TupleDID
{
    double d1;
    int i1;
    double d2;
};

struct TupleDDID
{
    double d1;
    double d2;
    int i1;
    double d3;
};


typedef struct DinvrState DinvrState;
typedef struct DzrorState DzrorState;


struct TupleDDID cdfbet_which1( double, double, double, double);
struct TupleDDID cdfbet_which2(double, double, double, double);
struct TupleDID cdfbet_which3(double, double, double, double, double);
struct TupleDID cdfbet_which4(double, double, double, double, double);
struct TupleDDID cdfbin_which1(double, double, double, double);
struct TupleDID cdfbin_which2(double, double, double, double, double);
struct TupleDID cdfbin_which3(double, double, double, double, double);
struct TupleDDID cdfbin_which4(double, double, double, double);
struct TupleDDID cdfchi_which1(double, double);
struct TupleDID cdfchi_which2(double, double, double);
struct TupleDID cdfchi_which3(double, double, double);
struct TupleDDID cdfchn_which1(double, double, double);
struct TupleDID cdfchn_which2(double, double, double);
struct TupleDID cdfchn_which3(double, double, double);
struct TupleDID cdfchn_which4(double, double, double);
struct TupleDDID cdff_which1(double, double, double);
struct TupleDID cdff_which2(double, double, double, double);
struct TupleDID cdff_which3(double, double, double, double);
struct TupleDID cdff_which4(double, double, double, double);
struct TupleDDID cdffnc_which1(double, double, double, double);
struct TupleDID cdffnc_which2(double, double, double, double, double);
struct TupleDID cdffnc_which3(double, double, double, double, double);
struct TupleDID cdffnc_which4(double, double, double, double, double);
struct TupleDID cdffnc_which5(double, double, double, double, double);
struct TupleDDID cdfgam_which1(double, double, double);
struct TupleDID cdfgam_which2(double, double, double, double);
struct TupleDID cdfgam_which3(double, double, double, double);
struct TupleDID cdfgam_which4(double, double, double, double);
struct TupleDDID cdfnbn_which1(double, double, double, double);
struct TupleDID cdfnbn_which2(double, double, double, double, double);
struct TupleDID cdfnbn_which3(double, double, double, double, double);
struct TupleDDID cdfnbn_which4(double, double, double, double);
struct TupleDDID cdfnor_which1(double, double, double);
struct TupleDID cdfnor_which2(double, double, double, double);
struct TupleDID cdfnor_which3(double, double, double, double);
struct TupleDID cdfnor_which4(double, double, double, double);
struct TupleDDID cdfpoi_which1(double, double);
struct TupleDID cdfpoi_which2(double, double, double);
struct TupleDID cdfpoi_which3(double, double, double);
struct TupleDDID cdft_which1(double, double);
struct TupleDID cdft_which2(double, double, double);
struct TupleDID cdft_which3(double, double, double);
struct TupleDDID cdftnc_which1(double, double, double);
struct TupleDID cdftnc_which2(double, double, double, double);
struct TupleDID cdftnc_which3(double, double, double, double);
struct TupleDID cdftnc_which4(double, double, double, double);

#ifdef __cplusplus
}      /* extern "C" */
#endif /* __cplusplus */

#endif /* ifndef */
