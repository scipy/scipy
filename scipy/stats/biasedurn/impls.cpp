// R_BUILD excludes some function implementations;
// patch them in here.

#include <stdexcept>
#include "stocc.h"

double StochasticLib1::Normal(double m, double s) {
  throw std::runtime_error("No normal generator found!");
}

void FatalError(const char* msg) {
  throw std::runtime_error(msg);
}