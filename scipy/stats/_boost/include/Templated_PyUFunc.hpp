#ifndef TEMPLATED_PYUFUNC_HPP
#define TEMPLATED_PYUFUNC_HPP

#include <cstddef>

// C++ standard dependent implementation
// TODO: remove C++11 implementation when SciPy supports C++14 on all platforms
#define USE_CPP14 (__cplusplus >= 201402L)

#if USE_CPP14

#include <type_traits>
#include <utility>
#define UFUNC_STATIC_ASSERT(cond, msg) static_assert(cond, msg)
#define UFUNC_CALLFUNC(RealType, NINPUTS, func, inputs)		      \
    *output = callfunc<RealType>(std::make_index_sequence<NINPUTS>{}, \
				 func, inputs);

// Variadic magic to cast to and call a function with NINPUTS
template <class RealType, std::size_t ... Is>
RealType
callfunc(std::index_sequence<Is ...>, void* func, RealType **inputs)
{
    return (reinterpret_cast<RealType(*)(...)>(func))(*inputs[Is] ...);
}

#else // USE_CPP14

#include <stdexcept>
#include <boost/static_assert.hpp>
#define UFUNC_STATIC_ASSERT(cond, msg) BOOST_STATIC_ASSERT_MSG(cond, msg)
#define UFUNC_CALLFUNC(RealType, NINPUTS, func, inputs) \
    *output = callfunc<RealType, NINPUTS>(func, inputs);

// Explicitly handle NINPUTS cases for C++11 compatibility
template<class RealType, std::size_t NINPUTS>
RealType
callfunc(void* func, RealType **inputs)
{
    switch (NINPUTS) {
    case 1:
        return (reinterpret_cast<RealType(*)(RealType)>(func))(*inputs[0]);
    case 2:
        return (reinterpret_cast<RealType(*)(RealType, RealType)>(func))(*inputs[0], *inputs[1]);
    case 3:
        return (reinterpret_cast<RealType(*)(RealType, RealType, RealType)>(func))(*inputs[0], *inputs[1], *inputs[2]);
    case 4:
        return (reinterpret_cast<RealType(*)(RealType, RealType, RealType, RealType)>(func))(*inputs[0], *inputs[1], *inputs[2], *inputs[3]);
    // case n:
    //     add more cases if/when you need them
    default:
        throw std::runtime_error("callfunc() in C++11 compatibility mode "
				 "handles too few args!");
    }
}

#endif // USE_CPP14


/** \brief UFUNC loop function that handles any type and number of inputs.
 *         Always returns a single output. */
template<class RealType = double, std::size_t NINPUTS>
static void
PyUFunc_T(char **args, npy_intp const *dimensions, npy_intp const *steps,
	  void *func)
{
    UFUNC_STATIC_ASSERT(NINPUTS > 0, "numpy.ufunc demands NINPUT > 0!");
    RealType *inputs[NINPUTS];
    for (std::size_t ii = 0; ii < NINPUTS; ++ii) {
        inputs[ii] = (RealType*)args[ii];
    }

    RealType *output = (RealType*)args[NINPUTS];
    const std::size_t sz = sizeof(RealType);
    for (npy_intp ii = 0; ii < dimensions[0]; ++ii) {
        UFUNC_CALLFUNC(RealType, NINPUTS, func, inputs);
        for (std::size_t jj = 0; jj < NINPUTS; ++jj) {
            inputs[jj] += steps[jj]/sz;
        }
        output += steps[NINPUTS]/sz;
    }
}


#endif // TEMPLATED_PYUFUNC_HPP
