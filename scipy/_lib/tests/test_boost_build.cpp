// Testing availability of boost headers

#include <boost/math/distributions.hpp>
#include <iostream>

void test() { boost::math::binomial_distribution<double> d(10, 0.5); }
