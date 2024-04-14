#pragma once

/*
RNAME = ['S', 'D']
RTYPE = ['float', 'double']
*/



/**
Approximate the steady state of a two-pole filer in polar form ran in backwards
for a step input.
**/
template<typename T>
T _hs(int k, T cs, double rsq, double omega)
{
    T cssq;
    T c0;
    double gamma, rsupk;

    cssq = cs * cs;
    k = abs(k);
    rsupk = pow(rsq, ((double ) k) / 2.0);
    if (omega == 0.0) {
	c0 = (1+rsq)/ ((1-rsq)*(1-rsq)*(1-rsq)) * cssq;
	gamma = (1-rsq) / (1+rsq);
	return c0 * rsupk * (1 + gamma * k);
    }
    if (omega == M_PI) {
	c0 = (1+rsq)/ ((1-rsq)*(1-rsq)*(1-rsq)) * cssq;
	gamma = (1-rsq) / (1+rsq) * (1 - 2 * (k % 2));
	return c0 * rsupk * (1 + gamma * k);
    }
    c0 = cssq * (1.0+rsq)/(1.0-rsq) / (1-2*rsq*cos(2*omega) + rsq*rsq);
    gamma = (1.0 - rsq)/ (1.0+rsq) / tan(omega);
    return c0 * rsupk * (cos(omega*k) + gamma * sin(omega * k));
}



/**
Approximate the steady state of a two-pole filter in polar form for a
step input.
**/
template<typename T>
T _hc(int k, T cs, double r, double omega)
{
    if (k < 0) return 0.0;
    if (omega == 0.0)
	return cs * pow(r, (double )k) * (k+1);
    else if (omega == M_PI)
	return cs * pow(r, (double )k) * (k+1) * (1 - 2*(k % 2));
    return cs * pow(r, (double) k) * sin(omega * (k+1)) / sin(omega);
}



/**
Compute the starting initial conditions for the system (ran in backwards)
cs / (1 - a2 * z - a3 * z^2) against signal x, where:

a2 = 2 * r * cos(omega)
a3 = - r ** 2
cs = 1 - 2 * r * cos(omega) + r ** 2

Arguments
---------
x: double* or float*
    2D strided pointer signal of size (M, N). When M > 1, multiple signals
    will be processed independently.
yp: double* or float*
    Output state condition pointer of size (M, 2), yp[:, 0] will contain the
    initial conditions y[n + 1], whereas yp[:, 1] will contain the initial
    condition y[n + 2].
M: int
    Number of signals to compute initial conditions for.
N: int
    Length of the signals.
precision: double* or float*
    Precision up to which the initial conditions will be computed.
**/
template<typename T>
int _SYM_IIR2_initial_bwd(double r, double omega, T *x, T *yp, int M, int N, T precision)
{
    double rsq = r * r;
    T cs = 1 - 2 * r * cos(omega) + rsq;

    // Fix starting values assuming mirror-symmetric boundary conditions.
    int k = 0;

    T err;
    T diff;

    do {
	    diff = (_hs(k, cs, rsq, omega) + _hs(k+1, cs, rsq, omega));
	    for(int i = 0; i < M; i++) {
            // Compute initial condition y[n + 1]
            yp[2 * i] += diff * x[N * i + N - 1 - k];
        }
	    err = diff * diff;
	    k++;
    } while((err > precision) && (k < N));

    if (k >= N) {return -3;}     // sum did not converge

    k = 0;
    do {
        diff = (_hs(k-1, cs, rsq, omega) + _hs(k+2, cs, rsq, omega));
        for(int i = 0; i < M; i++) {
            // Compute initial condition y[n + 2]
            yp[2 * i + 1] += diff * x[N * i + N - 1 - k];
        }
        err = diff * diff;
        k++;
    } while((err > precision) && (k < N));

    if (k >= N) {return -3;}     // sum did not converge

    return 0;
}


/*
{{for SUB, TYP in zip(RNAME, RTYPE)}}
int {{SUB}}_SYM_IIR2_initial_bwd(double r, double omega,
        {{TYP}} *x, {{TYP}} *yp, int M, int N, {{TYP}} precision) {
    double rsq = r * r;
    {{TYP}} cs = 1 - 2 * r * cos(omega) + rsq;

    // Fix starting values assuming mirror-symmetric boundary conditions.
    int k = 0;

    {{TYP}} err;
    {{TYP}} diff;

    do {
	    diff = ({{SUB}}_hs(k, cs, rsq, omega) + {{SUB}}_hs(k+1, cs, rsq, omega));
	    for(int i = 0; i < M; i++) {
            // Compute initial condition y[n + 1]
            yp[2 * i] += diff * x[N * i + N - 1 - k];
        }
	    err = diff * diff;
	    k++;
    } while((err > precision) && (k < N));

    if (k >= N) {return -3;}     // sum did not converge

    k = 0;
    do {
        diff = ({{SUB}}_hs(k-1, cs, rsq, omega) + {{SUB}}_hs(k+2, cs, rsq, omega));
        for(int i = 0; i < M; i++) {
            // Compute initial condition y[n + 2]
            yp[2 * i + 1] += diff * x[N * i + N - 1 - k];
        }
        err = diff * diff;
        k++;
    } while((err > precision) && (k < N));

    if (k >= N) {return -3;}     // sum did not converge

    return 0;
}
{{endfor}}
*/




/**
Compute the starting initial conditions for the system
cs / (1 - a2 * z^-1 - a3 * z^-2) against signals x, where:

a2 = 2 * r * cos(omega)
a3 = - r ** 2
cs = 1 - 2 * r * cos(omega) + r ** 2

Arguments
---------
x: double* or float*
    2D strided pointer signal of size (M, N). When M > 1, multiple signals
    will be processed independently.
yp: double* or float*
    Strided output state condition pointer of size (M, 2).
    yp[:, 0] will contain the initial conditions y[n - 1],
    whereas yp[:, 1] will contain the initial conditions y[n - 2].
M: int
    Number of signals to compute initial conditions for.
N: int
    Length of the signals.
precision: double* or float*
    Precision up to which the initial conditions will be computed.
**/
template<typename T>
int _SYM_IIR2_initial_fwd(double r, double omega, T *x, T *yp, int M, int N, T precision)
{
    /* Fix starting values assuming mirror-symmetric boundary conditions. */
    T cs = 1 - 2 * r * cos(omega) + r * r;

    for(int i = 0; i < M; i++) {
        // Compute starting condition y[n - 1]
        yp[2 * i] = _hc(0, cs, r, omega) * x[N * i];
    }

    int k = 0;
    precision *= precision;

    T err;
    T diff;

    do {
        diff = _hc(k+1, cs, r, omega);
        for(int i = 0; i < M; i++) {
            // Keep computing starting condition y[n - 1]
            yp[2 * i] += diff * x[N * i + k];
        }
        err = diff * diff;
        k++;
    } while((err > precision) && (k < N));

    if (k >= N) {return -3;}     /* sum did not converge */

    for(int i = 0; i < M; i++) {
        // Compute starting condition y[n - 2]
        yp[2 * i + 1] = _hc(0, cs, r, omega) * x[N * i + 1];
        yp[2 * i + 1] += _hc(1, cs, r, omega) * x[N * i];
    }

    k = 0;
    do {
        diff = _hc(k+2, cs, r, omega);
        for(int i = 0; i < M; i++) {
            // Keep computing starting condition y[n - 2]
            yp[2 * i + 1] += diff * x[N * i + k];
        }
        err = diff * diff;
        k++;
    } while((err > precision) && (k < N));

    if (k >= N) {return -3;}     /* sum did not converge */
    return 0;
}

