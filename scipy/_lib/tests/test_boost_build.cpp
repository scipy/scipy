// Testing availability of boost headers

#include <iostream>
#include <boost/math/distributions.hpp>
#include <boost/static_assert.hpp>

void test() {
  BOOST_STATIC_ASSERT(1 == 1);
  boost::math::binomial_distribution<double> d(10, 0.5);
}
