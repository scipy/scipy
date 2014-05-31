#include "math.h"
const double PI = 3.141592653589793238462643383279502884;
double
_multivariate_typical(int n, double args[n])
{
    return cos(args[1] * args[0] - args[2] * sin(args[0])) / PI;
}

double
_multivariate_indefinite(int n, double args[n])
{
    return -exp(-args[0]) * log(args[0]);
}

double
_multivariate_sin(int n, double args[n])
{
    return sin(args[0]);
}
 
