
#include <math.h>
extern double Gamma (double x);


#define EPS 1.0e-17

double besselpoly(double a, double lambda, double nu) {

  int m, factor=0;
  double Sm, relerr, Sol;
  double sum=0.0;

  /* Special handling for a = 0.0 */
  if (a == 0.0) {
    if (nu == 0.0) return 1.0/(lambda + 1);
    else return 0.0;
  }
  /* Special handling for negative and integer nu */
  if ((nu < 0) && (floor(nu)==nu)) {
    nu = -nu;
    factor = ((int) nu) % 2;
  }    
  Sm = exp(nu*log(a))/(Gamma(nu+1)*(lambda+nu+1));
  m = 0;
  do {
    sum += Sm;
    Sol = Sm;
    Sm *= -a*a*(lambda+nu+1+2*m)/((nu+m+1)*(m+1)*(lambda+nu+1+2*m+2));
    m++;
    relerr = fabs((Sm-Sol)/Sm);
  } while (relerr > EPS && m < 1000);
  if (!factor)
    return sum;
  else
    return -sum;
}






