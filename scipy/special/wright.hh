#ifndef WRIGHT_HH
#define WRIGHT_HH

#include <complex>

namespace wright {

int wrightomega_ext(std::complex<double> z, std::complex<double> *w,
		    std::complex<double> *cond);
std::complex<double> wrightomega(std::complex<double> z);

};
  
#endif /* wright.hh */
