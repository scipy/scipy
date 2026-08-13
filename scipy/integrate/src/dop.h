#ifndef DOP_H
#define DOP_H

#include <math.h>


/**
 * @brief Callback function type definition for computing the right-hand side of the ODE system.
 *
 * This function evaluates the derivatives for a system of first-order ordinary differential equations
 * at a given point (x, y).
 *
 * @param n   Dimension of the system (number of equations).
 * @param x   Current value of the independent variable.
 * @param y   Array of current values of the dependent variables.
 * @param f   Output array for the computed derivatives (dy/dx).
 * @param rpar User-supplied real parameter array (optional, can be NULL).
 * @param ipar User-supplied integer parameter array (optional, can be NULL).
 *
 * @return None. The computed derivatives should be stored in the array f.
 *
 * @note
 * - The function should fill the array f with the values of the derivatives at (x, y).
 * - Used as a callback by the integrator routines.
 */
typedef void dopri_fcn(int n, double x, double* y, double* f, double* rpar, int* ipar);


/**
 * @brief Callback function type definition for solution output during ODE integration.
 *
 * This function is called by the integrator at each accepted step or for dense output,
 * allowing the user to process or store intermediate results.
 *
 * @param nr      Step number.
 * @param xold    Previous value of the independent variable.
 * @param x       Current value of the independent variable.
 * @param y       Array of current values of the dependent variables.
 * @param n       Dimension of the system (number of equations).
 * @param con     Workspace array for dense output coefficients (if applicable).
 * @param icomp   Array of indices for components to output (if applicable).
 * @param nd      Number of components for dense output.
 * @param rpar    User-supplied real parameter array (optional, can be NULL).
 * @param ipar    User-supplied integer parameter array (optional, can be NULL).
 * @param irtrn   Pointer to integer return flag; set negative to terminate integration.
 *
 * @return Integer status code. Negative value requests termination of integration.
 *
 * @note
 * - This function can be used to store, print, or process the solution at each step.
 * - For dense output, additional information is provided via con and icomp.
 * - Setting *irtrn < 0 will terminate the integration early.
 */
typedef void dopri_solout(int nr, double xold, double x, double* y, int n, double* con, int* icomp, int nd, double* rpar, int* ipar, int*irtrn);


/**
 * @brief Integrates a system of first-order ODEs using the explicit Runge-Kutta DOP853 method.
 *
 * This function numerically solves the initial value problem for a system of ordinary differential equations
 * using the Dormand-Prince 8(5,3) method with adaptive step size control and optional dense output.
 *
 * @param n         Dimension of the system (number of equations).
 * @param fcn       Pointer to the function computing the right-hand side: fcn(n, x, y, f).
 * @param x         Initial value of the independent variable.
 * @param y         Array of initial values for the dependent variables (input/output).
 * @param xend      Pointer to the final value of the independent variable.
 * @param rtol      Array of relative error tolerances.
 * @param atol      Array of absolute error tolerances.
 * @param itol      Tolerance type: 0 for scalar, 1 for vector.
 * @param solout    Pointer to the solution output function (optional, can be NULL).
 * @param iout      Output control flag: 0 (no output), 1 (output at each step), 2 (dense output).
 * @param work      Workspace array for method parameters and internal storage.
 * @param iwork     Integer workspace array for method parameters and internal storage.
 * @param rpar      User-supplied real parameter array (passed to fcn and solout).
 * @param ipar      User-supplied integer parameter array (passed to fcn and solout).
 * @param ierr      Pointer to error code (output).
 *
 * @return None. Results are returned via the y array and ierr.
 *
 * @note
 * - The function supports both scalar and vector error tolerances.
 * - Dense output is available if iout == 2.
 * - The workspace arrays must be properly sized according to the documentation.
 * - Error codes are set in *ierr: 1 (success), negative values (various errors).
 */
void
dopri853(const int n, dopri_fcn* fcn, double* x, double* y, double* xend, double* rtol,
         double* atol, const int itol, dopri_solout* solout, const int iout, double* work,
         int* iwork, double* rpar, int* ipar, int* ierr);


/**
 * @brief Integrates a system of first-order ODEs using the explicit Runge-Kutta DOPRI5 method.
 *
 * This function numerically solves the initial value problem for a system of ordinary differential equations
 * using the Dormand-Prince 5th order method with adaptive step size control.
 *
 * @param n         Dimension of the system (number of equations).
 * @param fcn       Pointer to the function computing the right-hand side: fcn(n, x, y, f).
 * @param x         Initial value of the independent variable.
 * @param y         Array of initial values for the dependent variables (input/output).
 * @param xend      Pointer to the final value of the independent variable.
 * @param rtol      Array of relative error tolerances.
 * @param atol      Array of absolute error tolerances.
 * @param itol      Tolerance type: 0 for scalar, 1 for vector.
 * @param solout    Pointer to the solution output function (optional, can be NULL).
 * @param iout      Output control flag: 0 (no output), 1 (output at each step), 2 (dense output).
 * @param work      Workspace array for method parameters and internal storage.
 * @param iwork     Integer workspace array for method parameters and internal storage.
 * @param ierr      Pointer to error code (output).
 *
 * @return None. Results are returned via the y array and ierr.
 *
 * @note
 * - The function supports both scalar and vector error tolerances.
 * - The workspace arrays must be properly sized according to the documentation.
 * - Error codes are set in *ierr: 1 (success), negative values (various errors).
 */
void dopri5(const int n, dopri_fcn* fcn, double* x, double* y, double* xend, double* rtol,
       double* atol, const int itol, dopri_solout* solout, const int iout, double* work,
       int* iwork, double* rpar, int* ipar, int* ierr);


#endif
