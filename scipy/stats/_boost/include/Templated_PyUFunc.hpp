#ifndef TEMPLATED_PYUFUNC_HPP
#define TEMPLATED_PYUFUNC_HPP

#include <type_traits>
#include <utility>

#include <iostream>

/** \brief Variadic magic to cast to and call a function with NINPUTS */
template <class RealType, std::size_t ... Is>
RealType callfunc(std::index_sequence<Is ...>, void* func, RealType **inputs) {
  return (reinterpret_cast<RealType(*)(...)>(func))(*inputs[Is] ...);
}

/** \brief UFUNC loop function that handles any type and number of inputs.
 *         Always returns a single output. */
template<class RealType = double, std::size_t NINPUTS>
static void PyUFunc_T(char** args, npy_intp const *dimensions, npy_intp const *steps, void *func) {
  static_assert(NINPUTS > 0, "numpy.ufunc demands NINPUT > 0!");
  RealType *inputs[NINPUTS];
  for (std::size_t ii = 0; ii < NINPUTS; ++ii) {
    inputs[ii] = (RealType*)args[ii];
  }
  RealType *output = (RealType*)args[NINPUTS];
  const std::size_t sz = sizeof(RealType);
  for (npy_intp ii = 0; ii < dimensions[0]; ++ii) {
    *output = callfunc<RealType>(std::make_index_sequence<NINPUTS>{}, func, inputs);
    for (std::size_t jj = 0; jj < NINPUTS; ++jj) {
      inputs[jj] += steps[jj]/sz;
    }
    output += steps[NINPUTS]/sz;
  }
}


#endif
