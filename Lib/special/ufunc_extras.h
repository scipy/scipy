#include "Python.h"

typedef int    IntFunc_d_dd(double x, double *y, double *z);
typedef int    IntFunc_d_dddd(double v, double *w, double *x, double *y, double *z);
typedef int    IntFunc_D_DDDD(Py_complex v, Py_complex *w, Py_complex *x, Py_complex *y, Py_complex *z);
typedef int    IntFunc_dd_dddd(double u, double v, double *w, double *x, double *y, double *z);
typedef double DoubleFunc_id_d(int n, double y);
typedef Py_complex CmplxFunc_dD_D(double v, Py_complex z);
typedef double DoubleFunc_ddd_d(double x, double y, double z);
typedef double DoubleFunc_iid_d(int m, int n, double z);
typedef double DoubleFunc_dddd_d(double w, double x, double y, double z);
typedef double DoubleFunc_dddd_dd(double w, double x, double y, double z, double *err);
typedef double DoubleFunc_dddi_dd(double w, double x, double y, int z, double *err);


extern void PyUFunc_f_ff_As_d_dd(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_d_dd(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_dd_dd(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_f_ffff_As_d_dddd(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_d_dddd(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_f_ffff_As_D_DDDD(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_d_dddd_As_D_DDDD(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_F_FFFF_As_D_DDDD(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_D_DDDD(char **args, int *dimensions, int *steps, void *func);

extern void PyUFunc_ff_ffff_As_dd_dddd(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_dd_dddd(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_fff_f_As_ddd_d(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_ddd_d(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_ff_f_As_id_d(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_dd_d_As_id_d(char **args, int *dimensions, int *steps, void *func);

extern void PyUFunc_ff_f_As_dD_D(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_dd_d_As_dD_D(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_fF_F_As_dD_D(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_dD_D(char **args, int *dimensions, int *steps, void *func);

extern void PyUFunc_fff_f_As_iid_d(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_ddd_d_As_iid_d(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_ffff_f_As_dddd_d(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_dddd_d(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_ffff_ff_As_dddd_dd(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_dddd_dd(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_ffff_ff_As_dddi_dd(char **args, int *dimensions, int *steps, void *func);
extern void PyUFunc_dddd_dd_As_dddi_dd(char **args, int *dimensions, int *steps, void *func);


