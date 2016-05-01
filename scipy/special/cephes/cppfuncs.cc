#include "../Faddeeva.hh"

extern "C" double erfc(double x);

double erfc(double x)
{
  return Faddeeva::erfc(x);
}
