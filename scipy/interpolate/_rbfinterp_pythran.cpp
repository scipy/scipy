#include <pythonic/core.hpp>
#include <pythonic/python/core.hpp>
#include <pythonic/types/bool.hpp>
#include <pythonic/types/int.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <pythonic/include/types/ndarray.hpp>
#include <pythonic/include/types/float.hpp>
#include <pythonic/include/types/numpy_texpr.hpp>
#include <pythonic/include/types/str.hpp>
#include <pythonic/include/types/int.hpp>
#include <pythonic/types/str.hpp>
#include <pythonic/types/float.hpp>
#include <pythonic/types/int.hpp>
#include <pythonic/types/numpy_texpr.hpp>
#include <pythonic/types/ndarray.hpp>
#include <pythonic/include/builtins/None.hpp>
#include <pythonic/include/builtins/dict.hpp>
#include <pythonic/include/builtins/float_.hpp>
#include <pythonic/include/builtins/getattr.hpp>
#include <pythonic/include/builtins/pythran/make_shape.hpp>
#include <pythonic/include/builtins/range.hpp>
#include <pythonic/include/builtins/tuple.hpp>
#include <pythonic/include/numpy/empty.hpp>
#include <pythonic/include/numpy/exp.hpp>
#include <pythonic/include/numpy/linalg/norm.hpp>
#include <pythonic/include/numpy/log.hpp>
#include <pythonic/include/numpy/max.hpp>
#include <pythonic/include/numpy/min.hpp>
#include <pythonic/include/numpy/prod.hpp>
#include <pythonic/include/numpy/sqrt.hpp>
#include <pythonic/include/numpy/square.hpp>
#include <pythonic/include/numpy/zeros.hpp>
#include <pythonic/include/operator_/add.hpp>
#include <pythonic/include/operator_/div.hpp>
#include <pythonic/include/operator_/eq.hpp>
#include <pythonic/include/operator_/iadd.hpp>
#include <pythonic/include/operator_/mul.hpp>
#include <pythonic/include/operator_/neg.hpp>
#include <pythonic/include/operator_/pow.hpp>
#include <pythonic/include/operator_/sub.hpp>
#include <pythonic/include/types/slice.hpp>
#include <pythonic/include/types/str.hpp>
#include <pythonic/builtins/None.hpp>
#include <pythonic/builtins/dict.hpp>
#include <pythonic/builtins/float_.hpp>
#include <pythonic/builtins/getattr.hpp>
#include <pythonic/builtins/pythran/make_shape.hpp>
#include <pythonic/builtins/range.hpp>
#include <pythonic/builtins/tuple.hpp>
#include <pythonic/numpy/empty.hpp>
#include <pythonic/numpy/exp.hpp>
#include <pythonic/numpy/linalg/norm.hpp>
#include <pythonic/numpy/log.hpp>
#include <pythonic/numpy/max.hpp>
#include <pythonic/numpy/min.hpp>
#include <pythonic/numpy/prod.hpp>
#include <pythonic/numpy/sqrt.hpp>
#include <pythonic/numpy/square.hpp>
#include <pythonic/numpy/zeros.hpp>
#include <pythonic/operator_/add.hpp>
#include <pythonic/operator_/div.hpp>
#include <pythonic/operator_/eq.hpp>
#include <pythonic/operator_/iadd.hpp>
#include <pythonic/operator_/mul.hpp>
#include <pythonic/operator_/neg.hpp>
#include <pythonic/operator_/pow.hpp>
#include <pythonic/operator_/sub.hpp>
#include <pythonic/types/slice.hpp>
#include <pythonic/types/str.hpp>
namespace __pythran__rbfinterp_pythran
{
  struct polynomial_matrix
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::prod{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type2;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type1>())) __type4;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type4>::type>::type __type5;
      typedef decltype(std::declval<__type2>()(std::declval<__type5>())) __type6;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type6>::type::iterator>::value_type>::type __type7;
      typedef decltype(std::declval<__type1>()[std::declval<__type7>()]) __type8;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type9;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type9>())) __type11;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type11>::type>::type __type12;
      typedef decltype(std::declval<__type2>()(std::declval<__type12>())) __type13;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type13>::type::iterator>::value_type>::type __type14;
      typedef decltype(std::declval<__type9>()[std::declval<__type14>()]) __type15;
      typedef decltype(pythonic::builtins::pow(std::declval<__type8>(), std::declval<__type15>())) __type16;
      typedef decltype(std::declval<__type0>()(std::declval<__type16>())) __type17;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type7>(), std::declval<__type14>())) __type20;
      typedef __type17 __ptype0;
      typedef __type20 __ptype1;
      typedef typename pythonic::returnable<pythonic::types::none_type>::type result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    typename type<argument_type0, argument_type1, argument_type2>::result_type operator()(argument_type0&& x, argument_type1&& powers, argument_type2&& out) const
    ;
  }  ;
  struct kernel_matrix
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::linalg::functor::norm{})>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type3;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type2>())) __type5;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type5>::type>::type __type6;
      typedef decltype(std::declval<__type3>()(std::declval<__type6>())) __type7;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type7>::type::iterator>::value_type>::type __type8;
      typedef decltype(std::declval<__type2>()[std::declval<__type8>()]) __type9;
      typedef long __type12;
      typedef decltype(pythonic::operator_::add(std::declval<__type8>(), std::declval<__type12>())) __type13;
      typedef decltype(std::declval<__type3>()(std::declval<__type13>())) __type14;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type14>::type::iterator>::value_type>::type __type15;
      typedef decltype(std::declval<__type2>()[std::declval<__type15>()]) __type16;
      typedef decltype(pythonic::operator_::sub(std::declval<__type9>(), std::declval<__type16>())) __type17;
      typedef decltype(std::declval<__type1>()(std::declval<__type17>())) __type18;
      typedef decltype(std::declval<__type0>()(std::declval<__type18>())) __type19;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type8>(), std::declval<__type15>())) __type22;
      typedef __type19 __ptype4;
      typedef __type22 __ptype5;
      typedef typename pythonic::returnable<pythonic::types::none_type>::type result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    typename type<argument_type0, argument_type1, argument_type2>::result_type operator()(argument_type0&& x, argument_type1&& kernel_func, argument_type2&& out) const
    ;
  }  ;
  struct polynomial_vector
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::prod{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type2;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type3;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type2>())) __type5;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type5>::type>::type __type6;
      typedef decltype(std::declval<__type3>()(std::declval<__type6>())) __type7;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type7>::type::iterator>::value_type>::type __type8;
      typedef decltype(std::declval<__type2>()[std::declval<__type8>()]) __type9;
      typedef decltype(pythonic::builtins::pow(std::declval<__type1>(), std::declval<__type9>())) __type10;
      typedef decltype(std::declval<__type0>()(std::declval<__type10>())) __type11;
      typedef __type11 __ptype12;
      typedef __type8 __ptype13;
      typedef typename pythonic::returnable<pythonic::types::none_type>::type result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    typename type<argument_type0, argument_type1, argument_type2>::result_type operator()(argument_type0&& x, argument_type1&& powers, argument_type2&& out) const
    ;
  }  ;
  struct kernel_vector
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::linalg::functor::norm{})>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type3;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type4;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type3>())) __type6;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type __type7;
      typedef decltype(std::declval<__type4>()(std::declval<__type7>())) __type8;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type8>::type::iterator>::value_type>::type __type9;
      typedef decltype(std::declval<__type3>()[std::declval<__type9>()]) __type10;
      typedef decltype(pythonic::operator_::sub(std::declval<__type2>(), std::declval<__type10>())) __type11;
      typedef decltype(std::declval<__type1>()(std::declval<__type11>())) __type12;
      typedef decltype(std::declval<__type0>()(std::declval<__type12>())) __type13;
      typedef __type13 __ptype16;
      typedef __type9 __ptype17;
      typedef typename pythonic::returnable<pythonic::types::none_type>::type result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type operator()(argument_type0&& x, argument_type1&& y, argument_type2&& kernel_func, argument_type3&& out) const
    ;
  }  ;
  struct gaussian
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::exp{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::square{})>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
      typedef decltype(std::declval<__type1>()(std::declval<__type2>())) __type3;
      typedef decltype(pythonic::operator_::neg(std::declval<__type3>())) __type4;
      typedef typename pythonic::returnable<decltype(std::declval<__type0>()(std::declval<__type4>()))>::type result_type;
    }  
    ;
    template <typename argument_type0 >
    typename type<argument_type0>::result_type operator()(argument_type0&& r) const
    ;
  }  ;
  struct inverse_quadratic
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef long __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::square{})>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
      typedef decltype(std::declval<__type1>()(std::declval<__type2>())) __type3;
      typedef decltype(pythonic::operator_::add(std::declval<__type3>(), std::declval<__type0>())) __type4;
      typedef typename pythonic::returnable<decltype(pythonic::operator_::div(std::declval<__type0>(), std::declval<__type4>()))>::type result_type;
    }  
    ;
    template <typename argument_type0 >
    typename type<argument_type0>::result_type operator()(argument_type0&& r) const
    ;
  }  ;
  struct inverse_multiquadric
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef long __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::sqrt{})>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::square{})>::type>::type __type2;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type3;
      typedef decltype(std::declval<__type2>()(std::declval<__type3>())) __type4;
      typedef decltype(pythonic::operator_::add(std::declval<__type4>(), std::declval<__type0>())) __type5;
      typedef decltype(std::declval<__type1>()(std::declval<__type5>())) __type6;
      typedef typename pythonic::returnable<decltype(pythonic::operator_::div(std::declval<__type0>(), std::declval<__type6>()))>::type result_type;
    }  
    ;
    template <typename argument_type0 >
    typename type<argument_type0>::result_type operator()(argument_type0&& r) const
    ;
  }  ;
  struct multiquadric
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::sqrt{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::square{})>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
      typedef decltype(std::declval<__type1>()(std::declval<__type2>())) __type3;
      typedef long __type4;
      typedef decltype(pythonic::operator_::add(std::declval<__type3>(), std::declval<__type4>())) __type5;
      typedef decltype(std::declval<__type0>()(std::declval<__type5>())) __type6;
      typedef typename pythonic::returnable<decltype(pythonic::operator_::neg(std::declval<__type6>()))>::type result_type;
    }  
    ;
    template <typename argument_type0 >
    typename type<argument_type0>::result_type operator()(argument_type0&& r) const
    ;
  }  ;
  struct quintic
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::square{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type1;
      typedef decltype(std::declval<__type0>()(std::declval<__type1>())) __type2;
      typedef decltype(std::declval<__type0>()(std::declval<__type2>())) __type3;
      typedef decltype(pythonic::operator_::mul(std::declval<__type3>(), std::declval<__type1>())) __type5;
      typedef typename pythonic::returnable<decltype(pythonic::operator_::neg(std::declval<__type5>()))>::type result_type;
    }  
    ;
    template <typename argument_type0 >
    typename type<argument_type0>::result_type operator()(argument_type0&& r) const
    ;
  }  ;
  struct cubic
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::square{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type1;
      typedef decltype(std::declval<__type0>()(std::declval<__type1>())) __type2;
      typedef typename pythonic::returnable<decltype(pythonic::operator_::mul(std::declval<__type2>(), std::declval<__type1>()))>::type result_type;
    }  
    ;
    template <typename argument_type0 >
    typename type<argument_type0>::result_type operator()(argument_type0&& r) const
    ;
  }  ;
  struct thin_plate_spline
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef double __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::square{})>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
      typedef decltype(std::declval<__type1>()(std::declval<__type2>())) __type3;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::log{})>::type>::type __type4;
      typedef decltype(std::declval<__type4>()(std::declval<__type2>())) __type6;
      typedef decltype(pythonic::operator_::mul(std::declval<__type3>(), std::declval<__type6>())) __type7;
      typedef typename pythonic::returnable<typename __combined<__type0,__type7>::type>::type result_type;
    }  
    ;
    template <typename argument_type0 >
    typename type<argument_type0>::result_type operator()(argument_type0&& r) const
    ;
  }  ;
  struct linear
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
      typedef typename pythonic::returnable<decltype(pythonic::operator_::neg(std::declval<__type0>()))>::type result_type;
    }  
    ;
    template <typename argument_type0 >
    typename type<argument_type0>::result_type operator()(argument_type0&& r) const
    ;
  }  ;
  struct _polynomial_matrix
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::empty{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::pythran::functor::make_shape{})>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type2>())) __type3;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type3>::type>::type __type4;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type5;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type5>())) __type6;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type __type7;
      typedef decltype(std::declval<__type1>()(std::declval<__type4>(), std::declval<__type7>())) __type8;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::float_{})>::type>::type __type9;
      typedef typename pythonic::assignable<decltype(std::declval<__type0>()(std::declval<__type8>(), std::declval<__type9>()))>::type __type10;
      typedef typename polynomial_matrix::type<__type2, __type5, __type10>::__ptype0 __type14;
      typedef container<typename std::remove_reference<__type14>::type> __type15;
      typedef typename __combined<__type10,__type15>::type __type16;
      typedef typename polynomial_matrix::type<__type2, __type5, __type16>::__ptype1 __type18;
      typedef indexable<__type18> __type19;
      typedef typename __combined<__type16,__type19>::type __type20;
      typedef typename pythonic::returnable<typename __combined<__type20,__type20>::type>::type result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 >
    typename type<argument_type0, argument_type1>::result_type operator()(argument_type0&& x, argument_type1&& powers) const
    ;
  }  ;
  struct NAME_TO_FUNC
  {
    typedef void callable;
    typedef void pure;
    struct type
    {
      typedef pythonic::types::str __type0;
      typedef linear __type1;
      typedef pythonic::types::dict<__type0,__type1> __type2;
      typedef thin_plate_spline __type3;
      typedef pythonic::types::dict<__type0,__type3> __type4;
      typedef typename __combined<__type2,__type4>::type __type5;
      typedef cubic __type6;
      typedef pythonic::types::dict<__type0,__type6> __type7;
      typedef typename __combined<__type5,__type7>::type __type8;
      typedef quintic __type9;
      typedef pythonic::types::dict<__type0,__type9> __type10;
      typedef typename __combined<__type8,__type10>::type __type11;
      typedef multiquadric __type12;
      typedef pythonic::types::dict<__type0,__type12> __type13;
      typedef typename __combined<__type11,__type13>::type __type14;
      typedef inverse_multiquadric __type15;
      typedef pythonic::types::dict<__type0,__type15> __type16;
      typedef typename __combined<__type14,__type16>::type __type17;
      typedef inverse_quadratic __type18;
      typedef pythonic::types::dict<__type0,__type18> __type19;
      typedef typename __combined<__type17,__type19>::type __type20;
      typedef gaussian __type21;
      typedef pythonic::types::dict<__type0,__type21> __type22;
      typedef typename pythonic::returnable<typename __combined<__type20,__type22>::type>::type result_type;
    }  ;
    typename type::result_type operator()() const;
    ;
  }  ;
  struct _kernel_matrix
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::empty{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::pythran::functor::make_shape{})>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type2>())) __type3;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type3>::type>::type __type4;
      typedef decltype(std::declval<__type1>()(std::declval<__type4>(), std::declval<__type4>())) __type8;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::float_{})>::type>::type __type9;
      typedef typename pythonic::assignable<decltype(std::declval<__type0>()(std::declval<__type8>(), std::declval<__type9>()))>::type __type10;
      typedef pythonic::types::str __type12;
      typedef linear __type13;
      typedef pythonic::types::dict<__type12,__type13> __type14;
      typedef thin_plate_spline __type15;
      typedef pythonic::types::dict<__type12,__type15> __type16;
      typedef typename __combined<__type14,__type16>::type __type17;
      typedef cubic __type18;
      typedef pythonic::types::dict<__type12,__type18> __type19;
      typedef typename __combined<__type17,__type19>::type __type20;
      typedef quintic __type21;
      typedef pythonic::types::dict<__type12,__type21> __type22;
      typedef typename __combined<__type20,__type22>::type __type23;
      typedef multiquadric __type24;
      typedef pythonic::types::dict<__type12,__type24> __type25;
      typedef typename __combined<__type23,__type25>::type __type26;
      typedef inverse_multiquadric __type27;
      typedef pythonic::types::dict<__type12,__type27> __type28;
      typedef typename __combined<__type26,__type28>::type __type29;
      typedef inverse_quadratic __type30;
      typedef pythonic::types::dict<__type12,__type30> __type31;
      typedef typename __combined<__type29,__type31>::type __type32;
      typedef gaussian __type33;
      typedef pythonic::types::dict<__type12,__type33> __type34;
      typedef typename __combined<__type32,__type34>::type __type35;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type36;
      typedef decltype(std::declval<__type35>()[std::declval<__type36>()]) __type37;
      typedef typename kernel_matrix::type<__type2, __type37, __type10>::__ptype4 __type39;
      typedef container<typename std::remove_reference<__type39>::type> __type40;
      typedef typename __combined<__type10,__type40>::type __type41;
      typedef typename kernel_matrix::type<__type2, __type37, __type41>::__ptype5 __type43;
      typedef indexable<__type43> __type44;
      typedef typename __combined<__type41,__type44>::type __type45;
      typedef typename pythonic::returnable<typename __combined<__type45,__type45>::type>::type result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 >
    typename type<argument_type0, argument_type1>::result_type operator()(argument_type0&& x, argument_type1&& kernel) const
    ;
  }  ;
  struct _evaluate
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 , typename argument_type7 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::zeros{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::pythran::functor::make_shape{})>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type2>())) __type3;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type3>::type>::type __type4;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type7>::type>::type __type5;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type5>())) __type6;
      typedef typename std::tuple_element<1,typename std::remove_reference<__type6>::type>::type __type7;
      typedef decltype(std::declval<__type1>()(std::declval<__type4>(), std::declval<__type7>())) __type8;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::float_{})>::type>::type __type9;
      typedef typename pythonic::assignable<decltype(std::declval<__type0>()(std::declval<__type8>(), std::declval<__type9>()))>::type __type10;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type11;
      typedef decltype(std::declval<__type11>()(std::declval<__type4>())) __type12;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type12>::type::iterator>::value_type>::type __type13;
      typedef decltype(std::declval<__type11>()(std::declval<__type7>())) __type14;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type14>::type::iterator>::value_type>::type __type15;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type13>(), std::declval<__type15>())) __type16;
      typedef indexable<__type16> __type17;
      typedef typename __combined<__type10,__type17>::type __type18;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type20;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type20>())) __type21;
      typedef typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type21>::type>::type>::type __type22;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type23;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type23>())) __type24;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type24>::type>::type __type25;
      typedef decltype(pythonic::operator_::add(std::declval<__type22>(), std::declval<__type25>())) __type26;
      typedef decltype(std::declval<__type11>()(std::declval<__type26>())) __type27;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type27>::type::iterator>::value_type>::type __type28;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type28>(), std::declval<__type15>())) __type30;
      typedef decltype(std::declval<__type5>()[std::declval<__type30>()]) __type31;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::empty{})>::type>::type __type32;
      typedef decltype(std::declval<__type1>()(std::declval<__type26>())) __type35;
      typedef typename pythonic::assignable<decltype(std::declval<__type32>()(std::declval<__type35>(), std::declval<__type9>()))>::type __type36;
      typedef pythonic::types::contiguous_slice __type38;
      typedef decltype(std::declval<__type36>()[std::declval<__type38>()]) __type39;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type41;
      typedef typename pythonic::assignable<decltype(pythonic::operator_::mul(std::declval<__type2>(), std::declval<__type41>()))>::type __type42;
      typedef decltype(std::declval<__type42>()[std::declval<__type13>()]) __type44;
      typedef typename pythonic::assignable<decltype(pythonic::operator_::mul(std::declval<__type20>(), std::declval<__type41>()))>::type __type47;
      typedef pythonic::types::str __type48;
      typedef linear __type49;
      typedef pythonic::types::dict<__type48,__type49> __type50;
      typedef thin_plate_spline __type51;
      typedef pythonic::types::dict<__type48,__type51> __type52;
      typedef typename __combined<__type50,__type52>::type __type53;
      typedef cubic __type54;
      typedef pythonic::types::dict<__type48,__type54> __type55;
      typedef typename __combined<__type53,__type55>::type __type56;
      typedef quintic __type57;
      typedef pythonic::types::dict<__type48,__type57> __type58;
      typedef typename __combined<__type56,__type58>::type __type59;
      typedef multiquadric __type60;
      typedef pythonic::types::dict<__type48,__type60> __type61;
      typedef typename __combined<__type59,__type61>::type __type62;
      typedef inverse_multiquadric __type63;
      typedef pythonic::types::dict<__type48,__type63> __type64;
      typedef typename __combined<__type62,__type64>::type __type65;
      typedef inverse_quadratic __type66;
      typedef pythonic::types::dict<__type48,__type66> __type67;
      typedef typename __combined<__type65,__type67>::type __type68;
      typedef gaussian __type69;
      typedef pythonic::types::dict<__type48,__type69> __type70;
      typedef typename __combined<__type68,__type70>::type __type71;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type72;
      typedef typename pythonic::assignable<decltype(std::declval<__type71>()[std::declval<__type72>()])>::type __type73;
      typedef typename kernel_vector::type<__type44, __type47, __type73, __type39>::__ptype16 __type74;
      typedef container<typename std::remove_reference<__type74>::type> __type75;
      typedef typename __combined<__type39,__type75>::type __type76;
      typedef typename kernel_vector::type<__type44, __type47, __type73, __type76>::__ptype17 __type77;
      typedef indexable<__type77> __type78;
      typedef typename __combined<__type76,__type78>::type __type79;
      typedef typename __combined<__type75,__type78>::type __type80;
      typedef typename __combined<__type36,__type79,__type80>::type __type81;
      typedef decltype(std::declval<__type81>()[std::declval<__type38>()]) __type82;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type __type84;
      typedef decltype(pythonic::operator_::sub(std::declval<__type2>(), std::declval<__type84>())) __type85;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type6>::type>::type __type86;
      typedef typename pythonic::assignable<decltype(pythonic::operator_::div(std::declval<__type85>(), std::declval<__type86>()))>::type __type87;
      typedef decltype(std::declval<__type87>()[std::declval<__type13>()]) __type89;
      typedef typename polynomial_vector::type<__type89, __type23, __type82>::__ptype12 __type91;
      typedef container<typename std::remove_reference<__type91>::type> __type92;
      typedef typename __combined<__type82,__type92>::type __type93;
      typedef typename polynomial_vector::type<__type89, __type23, __type93>::__ptype13 __type94;
      typedef indexable<__type94> __type95;
      typedef typename __combined<__type93,__type95>::type __type96;
      typedef typename __combined<__type92,__type95>::type __type97;
      typedef typename __combined<__type36,__type79,__type80,__type96,__type97>::type __type98;
      typedef decltype(std::declval<__type98>()[std::declval<__type28>()]) __type100;
      typedef decltype(pythonic::operator_::mul(std::declval<__type31>(), std::declval<__type100>())) __type101;
      typedef container<typename std::remove_reference<__type101>::type> __type102;
      typedef typename pythonic::returnable<typename __combined<__type18,__type17,__type102>::type>::type result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 , typename argument_type7 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6, argument_type7>::result_type operator()(argument_type0&& x, argument_type1&& y, argument_type2&& kernel, argument_type3&& epsilon, argument_type4&& powers, argument_type5&& shift, argument_type6&& scale, argument_type7&& coeffs) const
    ;
  }  ;
  struct _build_system
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::empty{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::pythran::functor::make_shape{})>::type>::type __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type2;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type2>())) __type3;
      typedef typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type3>::type>::type>::type __type4;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type __type5;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type5>())) __type6;
      typedef typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type>::type __type7;
      typedef decltype(pythonic::operator_::add(std::declval<__type4>(), std::declval<__type7>())) __type8;
      typedef decltype(std::declval<__type1>()(std::declval<__type8>(), std::declval<__type8>())) __type12;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::float_{})>::type>::type __type13;
      typedef decltype(std::declval<__type0>()(std::declval<__type12>(), std::declval<__type13>())) __type14;
      typedef typename pythonic::assignable<decltype(pythonic::builtins::getattr(pythonic::types::attr::T{}, std::declval<__type14>()))>::type __type15;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type16;
      typedef decltype(std::declval<__type16>()(std::declval<__type4>())) __type18;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type18>::type::iterator>::value_type>::type __type19;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type19>(), std::declval<__type19>())) __type21;
      typedef indexable<__type21> __type22;
      typedef typename __combined<__type15,__type22>::type __type23;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type24;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type25;
      typedef decltype(pythonic::operator_::mul(std::declval<__type24>(), std::declval<__type25>())) __type26;
      typedef pythonic::types::str __type27;
      typedef linear __type28;
      typedef pythonic::types::dict<__type27,__type28> __type29;
      typedef thin_plate_spline __type30;
      typedef pythonic::types::dict<__type27,__type30> __type31;
      typedef typename __combined<__type29,__type31>::type __type32;
      typedef cubic __type33;
      typedef pythonic::types::dict<__type27,__type33> __type34;
      typedef typename __combined<__type32,__type34>::type __type35;
      typedef quintic __type36;
      typedef pythonic::types::dict<__type27,__type36> __type37;
      typedef typename __combined<__type35,__type37>::type __type38;
      typedef multiquadric __type39;
      typedef pythonic::types::dict<__type27,__type39> __type40;
      typedef typename __combined<__type38,__type40>::type __type41;
      typedef inverse_multiquadric __type42;
      typedef pythonic::types::dict<__type27,__type42> __type43;
      typedef typename __combined<__type41,__type43>::type __type44;
      typedef inverse_quadratic __type45;
      typedef pythonic::types::dict<__type27,__type45> __type46;
      typedef typename __combined<__type44,__type46>::type __type47;
      typedef gaussian __type48;
      typedef pythonic::types::dict<__type27,__type48> __type49;
      typedef typename __combined<__type47,__type49>::type __type50;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type51;
      typedef decltype(std::declval<__type50>()[std::declval<__type51>()]) __type52;
      typedef pythonic::types::contiguous_slice __type54;
      typedef decltype(std::declval<__type15>()(std::declval<__type54>(), std::declval<__type54>())) __type55;
      typedef typename kernel_matrix::type<__type26, __type52, __type55>::__ptype4 __type56;
      typedef container<typename std::remove_reference<__type56>::type> __type57;
      typedef container<typename std::remove_reference<__type57>::type> __type58;
      typedef typename kernel_matrix::type<__type26, __type52, __type55>::__ptype5 __type59;
      typedef indexable<__type59> __type60;
      typedef container<typename std::remove_reference<__type60>::type> __type61;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::max{})>::type>::type __type63;
      typedef long __type65;
      typedef typename pythonic::assignable<decltype(std::declval<__type63>()(std::declval<__type24>(), std::declval<__type65>()))>::type __type66;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::min{})>::type>::type __type67;
      typedef typename pythonic::assignable<decltype(std::declval<__type67>()(std::declval<__type24>(), std::declval<__type65>()))>::type __type69;
      typedef decltype(pythonic::operator_::add(std::declval<__type66>(), std::declval<__type69>())) __type70;
      typedef typename pythonic::assignable<decltype(pythonic::operator_::div(std::declval<__type70>(), std::declval<__type65>()))>::type __type71;
      typedef decltype(pythonic::operator_::sub(std::declval<__type24>(), std::declval<__type71>())) __type72;
      typedef decltype(pythonic::operator_::sub(std::declval<__type66>(), std::declval<__type69>())) __type75;
      typedef typename pythonic::assignable<decltype(pythonic::operator_::div(std::declval<__type75>(), std::declval<__type65>()))>::type __type76;
      typedef double __type77;
      typedef container<typename std::remove_reference<__type77>::type> __type78;
      typedef typename __combined<__type76,__type78>::type __type79;
      typedef decltype(pythonic::operator_::eq(std::declval<__type79>(), std::declval<__type77>())) __type80;
      typedef indexable<__type80> __type81;
      typedef typename __combined<__type76,__type81>::type __type82;
      typedef typename __combined<__type82,__type78,__type81,__type78>::type __type83;
      typedef decltype(pythonic::operator_::div(std::declval<__type72>(), std::declval<__type83>())) __type84;
      typedef typename __combined<__type15,__type58,__type61>::type __type86;
      typedef decltype(std::declval<__type86>()(std::declval<__type54>(), std::declval<__type54>())) __type87;
      typedef typename polynomial_matrix::type<__type84, __type5, __type87>::__ptype0 __type88;
      typedef container<typename std::remove_reference<__type88>::type> __type89;
      typedef container<typename std::remove_reference<__type89>::type> __type90;
      typedef typename polynomial_matrix::type<__type84, __type5, __type87>::__ptype1 __type91;
      typedef indexable<__type91> __type92;
      typedef container<typename std::remove_reference<__type92>::type> __type93;
      typedef typename __combined<__type15,__type58,__type61,__type90,__type93>::type __type94;
      typedef decltype(std::declval<__type94>()(std::declval<__type54>(), std::declval<__type54>())) __type95;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::T{}, std::declval<__type95>())) __type96;
      typedef container<typename std::remove_reference<__type96>::type> __type97;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type98;
      typedef decltype(std::declval<__type98>()[std::declval<__type19>()]) __type100;
      typedef container<typename std::remove_reference<__type100>::type> __type101;
      typedef typename __combined<__type23,__type58,__type61,__type90,__type93,__type97,__type78,__type22,__type101>::type __type102;
      typedef typename std::tuple_element<1,typename std::remove_reference<__type3>::type>::type __type105;
      typedef decltype(std::declval<__type1>()(std::declval<__type105>(), std::declval<__type8>())) __type109;
      typedef decltype(std::declval<__type0>()(std::declval<__type109>(), std::declval<__type13>())) __type110;
      typedef typename pythonic::assignable<decltype(pythonic::builtins::getattr(pythonic::types::attr::T{}, std::declval<__type110>()))>::type __type111;
      typedef typename __combined<__type111,__type2,__type77>::type __type113;
      typedef typename __combined<__type82,__type78,__type81>::type __type115;
      typedef typename pythonic::returnable<decltype(pythonic::types::make_tuple(std::declval<__type102>(), std::declval<__type113>(), std::declval<__type71>(), std::declval<__type115>()))>::type result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5>::result_type operator()(argument_type0&& y, argument_type1&& d, argument_type2&& smoothing, argument_type3&& kernel, argument_type4&& epsilon, argument_type5&& powers) const
    ;
  }  ;
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
  typename polynomial_matrix::type<argument_type0, argument_type1, argument_type2>::result_type polynomial_matrix::operator()(argument_type0&& x, argument_type1&& powers, argument_type2&& out) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type0;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type1;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type1>())) __type2;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type2>::type>::type __type3;
    typedef decltype(std::declval<__type0>()(std::declval<__type3>())) __type4;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type5;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type5>())) __type6;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type __type7;
    typedef decltype(std::declval<__type0>()(std::declval<__type7>())) __type8;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type4>::type::iterator>::value_type>::type>::type j;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type8>::type::iterator>::value_type>::type>::type i;
    {
      long  __target6475961920 = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, x));
      for (long  i=0L; i < __target6475961920; i += 1L)
      {
        {
          long  __target6475984656 = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, powers));
          for (long  j=0L; j < __target6475984656; j += 1L)
          {
            out.fast(pythonic::types::make_tuple(i, j)) = pythonic::numpy::functor::prod{}(pythonic::builtins::pow(x.fast(i), powers.fast(j)));
          }
        }
      }
    }
    return pythonic::builtins::None;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
  typename kernel_matrix::type<argument_type0, argument_type1, argument_type2>::result_type kernel_matrix::operator()(argument_type0&& x, argument_type1&& kernel_func, argument_type2&& out) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type0;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type1;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type1>())) __type2;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type2>::type>::type __type3;
    typedef decltype(std::declval<__type0>()(std::declval<__type3>())) __type4;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type4>::type::iterator>::value_type>::type __type5;
    typedef long __type6;
    typedef decltype(pythonic::operator_::add(std::declval<__type5>(), std::declval<__type6>())) __type7;
    typedef decltype(std::declval<__type0>()(std::declval<__type7>())) __type8;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type8>::type::iterator>::value_type>::type>::type j;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type4>::type::iterator>::value_type>::type>::type i;
    {
      long  __target6475950592 = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, x));
      for (long  i=0L; i < __target6475950592; i += 1L)
      {
        {
          long  __target6475962112 = pythonic::operator_::add(i, 1L);
          for (long  j=0L; j < __target6475962112; j += 1L)
          {
            out.fast(pythonic::types::make_tuple(i, j)) = kernel_func(pythonic::numpy::linalg::functor::norm{}(pythonic::operator_::sub(x.fast(i), x.fast(j))));
            out.fast(pythonic::types::make_tuple(j, i)) = out.fast(pythonic::types::make_tuple(i, j));
          }
        }
      }
    }
    return pythonic::builtins::None;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
  typename polynomial_vector::type<argument_type0, argument_type1, argument_type2>::result_type polynomial_vector::operator()(argument_type0&& x, argument_type1&& powers, argument_type2&& out) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type0;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type1;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type1>())) __type2;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type2>::type>::type __type3;
    typedef decltype(std::declval<__type0>()(std::declval<__type3>())) __type4;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type4>::type::iterator>::value_type>::type>::type i;
    {
      long  __target6475934880 = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, powers));
      for (long  i=0L; i < __target6475934880; i += 1L)
      {
        out.fast(i) = pythonic::numpy::functor::prod{}(pythonic::builtins::pow(x, powers.fast(i)));
      }
    }
    return pythonic::builtins::None;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
  typename kernel_vector::type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type kernel_vector::operator()(argument_type0&& x, argument_type1&& y, argument_type2&& kernel_func, argument_type3&& out) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type0;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type1;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type1>())) __type2;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type2>::type>::type __type3;
    typedef decltype(std::declval<__type0>()(std::declval<__type3>())) __type4;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type4>::type::iterator>::value_type>::type>::type i;
    {
      long  __target6475934496 = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, y));
      for (long  i=0L; i < __target6475934496; i += 1L)
      {
        out.fast(i) = kernel_func(pythonic::numpy::linalg::functor::norm{}(pythonic::operator_::sub(x, y.fast(i))));
      }
    }
    return pythonic::builtins::None;
  }
  template <typename argument_type0 >
  typename gaussian::type<argument_type0>::result_type gaussian::operator()(argument_type0&& r) const
  {
    return pythonic::numpy::functor::exp{}(pythonic::operator_::neg(pythonic::numpy::functor::square{}(r)));
  }
  template <typename argument_type0 >
  typename inverse_quadratic::type<argument_type0>::result_type inverse_quadratic::operator()(argument_type0&& r) const
  {
    return pythonic::operator_::div(1L, pythonic::operator_::add(pythonic::numpy::functor::square{}(r), 1L));
  }
  template <typename argument_type0 >
  typename inverse_multiquadric::type<argument_type0>::result_type inverse_multiquadric::operator()(argument_type0&& r) const
  {
    return pythonic::operator_::div(1L, pythonic::numpy::functor::sqrt{}(pythonic::operator_::add(pythonic::numpy::functor::square{}(r), 1L)));
  }
  template <typename argument_type0 >
  typename multiquadric::type<argument_type0>::result_type multiquadric::operator()(argument_type0&& r) const
  {
    return pythonic::operator_::neg(pythonic::numpy::functor::sqrt{}(pythonic::operator_::add(pythonic::numpy::functor::square{}(r), 1L)));
  }
  template <typename argument_type0 >
  typename quintic::type<argument_type0>::result_type quintic::operator()(argument_type0&& r) const
  {
    return pythonic::operator_::neg(pythonic::operator_::mul(pythonic::numpy::functor::square{}(pythonic::numpy::functor::square{}(r)), r));
  }
  template <typename argument_type0 >
  typename cubic::type<argument_type0>::result_type cubic::operator()(argument_type0&& r) const
  {
    return pythonic::operator_::mul(pythonic::numpy::functor::square{}(r), r);
  }
  template <typename argument_type0 >
  typename thin_plate_spline::type<argument_type0>::result_type thin_plate_spline::operator()(argument_type0&& r) const
  {
    if (pythonic::operator_::eq(r, 0L))
    {
      return 0.0;
    }
    else
    {
      return pythonic::operator_::mul(pythonic::numpy::functor::square{}(r), pythonic::numpy::functor::log{}(r));
    }
  }
  template <typename argument_type0 >
  typename linear::type<argument_type0>::result_type linear::operator()(argument_type0&& r) const
  {
    return pythonic::operator_::neg(r);
  }
  template <typename argument_type0 , typename argument_type1 >
  typename _polynomial_matrix::type<argument_type0, argument_type1>::result_type _polynomial_matrix::operator()(argument_type0&& x, argument_type1&& powers) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::empty{})>::type>::type __type0;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::pythran::functor::make_shape{})>::type>::type __type1;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type2>())) __type3;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type3>::type>::type __type4;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type5;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type5>())) __type6;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type __type7;
    typedef decltype(std::declval<__type1>()(std::declval<__type4>(), std::declval<__type7>())) __type8;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::float_{})>::type>::type __type9;
    typedef typename pythonic::assignable<decltype(std::declval<__type0>()(std::declval<__type8>(), std::declval<__type9>()))>::type __type10;
    typedef typename polynomial_matrix::type<__type2, __type5, __type10>::__ptype0 __type14;
    typedef container<typename std::remove_reference<__type14>::type> __type15;
    typedef typename __combined<__type10,__type15>::type __type16;
    typedef typename polynomial_matrix::type<__type2, __type5, __type16>::__ptype1 __type18;
    typedef indexable<__type18> __type19;
    typedef typename __combined<__type16,__type19>::type __type20;
    typename pythonic::assignable<typename __combined<__type20,__type20>::type>::type out = pythonic::numpy::functor::empty{}(pythonic::builtins::pythran::functor::make_shape{}(std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, x)), std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, powers))), pythonic::builtins::functor::float_{});
    polynomial_matrix()(x, powers, out);
    return out;
  }
  typename NAME_TO_FUNC::type::result_type NAME_TO_FUNC::operator()() const
  {
    {
      static typename NAME_TO_FUNC::type::result_type tmp_global = typename pythonic::assignable<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<pythonic::types::dict<pythonic::types::str,linear>,pythonic::types::dict<pythonic::types::str,thin_plate_spline>>::type,pythonic::types::dict<pythonic::types::str,cubic>>::type,pythonic::types::dict<pythonic::types::str,quintic>>::type,pythonic::types::dict<pythonic::types::str,multiquadric>>::type,pythonic::types::dict<pythonic::types::str,inverse_multiquadric>>::type,pythonic::types::dict<pythonic::types::str,inverse_quadratic>>::type,pythonic::types::dict<pythonic::types::str,gaussian>>::type>::type{{{ pythonic::types::str("linear"), linear() }, { pythonic::types::str("thin_plate_spline"), thin_plate_spline() }, { pythonic::types::str("cubic"), cubic() }, { pythonic::types::str("quintic"), quintic() }, { pythonic::types::str("multiquadric"), multiquadric() }, { pythonic::types::str("inverse_multiquadric"), inverse_multiquadric() }, { pythonic::types::str("inverse_quadratic"), inverse_quadratic() }, { pythonic::types::str("gaussian"), gaussian() }}};
      return tmp_global;
    }
  }
  template <typename argument_type0 , typename argument_type1 >
  typename _kernel_matrix::type<argument_type0, argument_type1>::result_type _kernel_matrix::operator()(argument_type0&& x, argument_type1&& kernel) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::empty{})>::type>::type __type0;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::pythran::functor::make_shape{})>::type>::type __type1;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type2>())) __type3;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type3>::type>::type __type4;
    typedef decltype(std::declval<__type1>()(std::declval<__type4>(), std::declval<__type4>())) __type8;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::float_{})>::type>::type __type9;
    typedef typename pythonic::assignable<decltype(std::declval<__type0>()(std::declval<__type8>(), std::declval<__type9>()))>::type __type10;
    typedef pythonic::types::str __type12;
    typedef linear __type13;
    typedef pythonic::types::dict<__type12,__type13> __type14;
    typedef thin_plate_spline __type15;
    typedef pythonic::types::dict<__type12,__type15> __type16;
    typedef typename __combined<__type14,__type16>::type __type17;
    typedef cubic __type18;
    typedef pythonic::types::dict<__type12,__type18> __type19;
    typedef typename __combined<__type17,__type19>::type __type20;
    typedef quintic __type21;
    typedef pythonic::types::dict<__type12,__type21> __type22;
    typedef typename __combined<__type20,__type22>::type __type23;
    typedef multiquadric __type24;
    typedef pythonic::types::dict<__type12,__type24> __type25;
    typedef typename __combined<__type23,__type25>::type __type26;
    typedef inverse_multiquadric __type27;
    typedef pythonic::types::dict<__type12,__type27> __type28;
    typedef typename __combined<__type26,__type28>::type __type29;
    typedef inverse_quadratic __type30;
    typedef pythonic::types::dict<__type12,__type30> __type31;
    typedef typename __combined<__type29,__type31>::type __type32;
    typedef gaussian __type33;
    typedef pythonic::types::dict<__type12,__type33> __type34;
    typedef typename __combined<__type32,__type34>::type __type35;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type36;
    typedef decltype(std::declval<__type35>()[std::declval<__type36>()]) __type37;
    typedef typename kernel_matrix::type<__type2, __type37, __type10>::__ptype4 __type39;
    typedef container<typename std::remove_reference<__type39>::type> __type40;
    typedef typename __combined<__type10,__type40>::type __type41;
    typedef typename kernel_matrix::type<__type2, __type37, __type41>::__ptype5 __type43;
    typedef indexable<__type43> __type44;
    typedef typename __combined<__type41,__type44>::type __type45;
    typename pythonic::assignable<typename __combined<__type45,__type45>::type>::type out = pythonic::numpy::functor::empty{}(pythonic::builtins::pythran::functor::make_shape{}(std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, x)), std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, x))), pythonic::builtins::functor::float_{});
    kernel_matrix()(x, typename pythonic::assignable<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<pythonic::types::dict<pythonic::types::str,linear>,pythonic::types::dict<pythonic::types::str,thin_plate_spline>>::type,pythonic::types::dict<pythonic::types::str,cubic>>::type,pythonic::types::dict<pythonic::types::str,quintic>>::type,pythonic::types::dict<pythonic::types::str,multiquadric>>::type,pythonic::types::dict<pythonic::types::str,inverse_multiquadric>>::type,pythonic::types::dict<pythonic::types::str,inverse_quadratic>>::type,pythonic::types::dict<pythonic::types::str,gaussian>>::type>::type{{{ pythonic::types::str("linear"), linear() }, { pythonic::types::str("thin_plate_spline"), thin_plate_spline() }, { pythonic::types::str("cubic"), cubic() }, { pythonic::types::str("quintic"), quintic() }, { pythonic::types::str("multiquadric"), multiquadric() }, { pythonic::types::str("inverse_multiquadric"), inverse_multiquadric() }, { pythonic::types::str("inverse_quadratic"), inverse_quadratic() }, { pythonic::types::str("gaussian"), gaussian() }}}[kernel], out);
    return out;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 , typename argument_type7 >
  typename _evaluate::type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6, argument_type7>::result_type _evaluate::operator()(argument_type0&& x, argument_type1&& y, argument_type2&& kernel, argument_type3&& epsilon, argument_type4&& powers, argument_type5&& shift, argument_type6&& scale, argument_type7&& coeffs) const
  {
    typedef pythonic::types::str __type0;
    typedef linear __type1;
    typedef pythonic::types::dict<__type0,__type1> __type2;
    typedef thin_plate_spline __type3;
    typedef pythonic::types::dict<__type0,__type3> __type4;
    typedef typename __combined<__type2,__type4>::type __type5;
    typedef cubic __type6;
    typedef pythonic::types::dict<__type0,__type6> __type7;
    typedef typename __combined<__type5,__type7>::type __type8;
    typedef quintic __type9;
    typedef pythonic::types::dict<__type0,__type9> __type10;
    typedef typename __combined<__type8,__type10>::type __type11;
    typedef multiquadric __type12;
    typedef pythonic::types::dict<__type0,__type12> __type13;
    typedef typename __combined<__type11,__type13>::type __type14;
    typedef inverse_multiquadric __type15;
    typedef pythonic::types::dict<__type0,__type15> __type16;
    typedef typename __combined<__type14,__type16>::type __type17;
    typedef inverse_quadratic __type18;
    typedef pythonic::types::dict<__type0,__type18> __type19;
    typedef typename __combined<__type17,__type19>::type __type20;
    typedef gaussian __type21;
    typedef pythonic::types::dict<__type0,__type21> __type22;
    typedef typename __combined<__type20,__type22>::type __type23;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type24;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::zeros{})>::type>::type __type25;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::pythran::functor::make_shape{})>::type>::type __type26;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type27;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type27>())) __type28;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type28>::type>::type __type29;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type7>::type>::type __type30;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type30>())) __type31;
    typedef typename std::tuple_element<1,typename std::remove_reference<__type31>::type>::type __type32;
    typedef decltype(std::declval<__type26>()(std::declval<__type29>(), std::declval<__type32>())) __type33;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::float_{})>::type>::type __type34;
    typedef typename pythonic::assignable<decltype(std::declval<__type25>()(std::declval<__type33>(), std::declval<__type34>()))>::type __type35;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type36;
    typedef decltype(std::declval<__type36>()(std::declval<__type29>())) __type37;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type37>::type::iterator>::value_type>::type __type38;
    typedef decltype(std::declval<__type36>()(std::declval<__type32>())) __type39;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type39>::type::iterator>::value_type>::type __type40;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type38>(), std::declval<__type40>())) __type41;
    typedef indexable<__type41> __type42;
    typedef typename __combined<__type35,__type42>::type __type43;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type45;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type45>())) __type46;
    typedef typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type46>::type>::type>::type __type47;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type48;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type48>())) __type49;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type49>::type>::type __type50;
    typedef decltype(pythonic::operator_::add(std::declval<__type47>(), std::declval<__type50>())) __type51;
    typedef decltype(std::declval<__type36>()(std::declval<__type51>())) __type52;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type52>::type::iterator>::value_type>::type __type53;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type53>(), std::declval<__type40>())) __type55;
    typedef decltype(std::declval<__type30>()[std::declval<__type55>()]) __type56;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::empty{})>::type>::type __type57;
    typedef decltype(std::declval<__type26>()(std::declval<__type51>())) __type60;
    typedef typename pythonic::assignable<decltype(std::declval<__type57>()(std::declval<__type60>(), std::declval<__type34>()))>::type __type61;
    typedef pythonic::types::contiguous_slice __type63;
    typedef decltype(std::declval<__type61>()[std::declval<__type63>()]) __type64;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type66;
    typedef typename pythonic::assignable<decltype(pythonic::operator_::mul(std::declval<__type27>(), std::declval<__type66>()))>::type __type67;
    typedef decltype(std::declval<__type67>()[std::declval<__type38>()]) __type69;
    typedef typename pythonic::assignable<decltype(pythonic::operator_::mul(std::declval<__type45>(), std::declval<__type66>()))>::type __type72;
    typedef typename pythonic::assignable<decltype(std::declval<__type23>()[std::declval<__type24>()])>::type __type73;
    typedef typename kernel_vector::type<__type69, __type72, __type73, __type64>::__ptype16 __type74;
    typedef container<typename std::remove_reference<__type74>::type> __type75;
    typedef typename __combined<__type64,__type75>::type __type76;
    typedef typename kernel_vector::type<__type69, __type72, __type73, __type76>::__ptype17 __type77;
    typedef indexable<__type77> __type78;
    typedef typename __combined<__type76,__type78>::type __type79;
    typedef typename __combined<__type75,__type78>::type __type80;
    typedef typename __combined<__type61,__type79,__type80>::type __type81;
    typedef decltype(std::declval<__type81>()[std::declval<__type63>()]) __type82;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type __type84;
    typedef decltype(pythonic::operator_::sub(std::declval<__type27>(), std::declval<__type84>())) __type85;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type6>::type>::type __type86;
    typedef typename pythonic::assignable<decltype(pythonic::operator_::div(std::declval<__type85>(), std::declval<__type86>()))>::type __type87;
    typedef decltype(std::declval<__type87>()[std::declval<__type38>()]) __type89;
    typedef typename polynomial_vector::type<__type89, __type48, __type82>::__ptype12 __type91;
    typedef container<typename std::remove_reference<__type91>::type> __type92;
    typedef typename __combined<__type82,__type92>::type __type93;
    typedef typename polynomial_vector::type<__type89, __type48, __type93>::__ptype13 __type94;
    typedef indexable<__type94> __type95;
    typedef typename __combined<__type93,__type95>::type __type96;
    typedef typename __combined<__type92,__type95>::type __type97;
    typedef typename __combined<__type61,__type79,__type80,__type96,__type97>::type __type98;
    typedef decltype(std::declval<__type98>()[std::declval<__type53>()]) __type100;
    typedef decltype(pythonic::operator_::mul(std::declval<__type56>(), std::declval<__type100>())) __type101;
    typedef container<typename std::remove_reference<__type101>::type> __type102;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type52>::type::iterator>::value_type>::type>::type k;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type39>::type::iterator>::value_type>::type>::type j;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type37>::type::iterator>::value_type>::type>::type i;
    typename pythonic::assignable_noescape<decltype(std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, y)))>::type p = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, y));
    typename pythonic::assignable<typename pythonic::assignable<decltype(std::declval<__type23>()[std::declval<__type24>()])>::type>::type kernel_func = typename pythonic::assignable<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<pythonic::types::dict<pythonic::types::str,linear>,pythonic::types::dict<pythonic::types::str,thin_plate_spline>>::type,pythonic::types::dict<pythonic::types::str,cubic>>::type,pythonic::types::dict<pythonic::types::str,quintic>>::type,pythonic::types::dict<pythonic::types::str,multiquadric>>::type,pythonic::types::dict<pythonic::types::str,inverse_multiquadric>>::type,pythonic::types::dict<pythonic::types::str,inverse_quadratic>>::type,pythonic::types::dict<pythonic::types::str,gaussian>>::type>::type{{{ pythonic::types::str("linear"), linear() }, { pythonic::types::str("thin_plate_spline"), thin_plate_spline() }, { pythonic::types::str("cubic"), cubic() }, { pythonic::types::str("quintic"), quintic() }, { pythonic::types::str("multiquadric"), multiquadric() }, { pythonic::types::str("inverse_multiquadric"), inverse_multiquadric() }, { pythonic::types::str("inverse_quadratic"), inverse_quadratic() }, { pythonic::types::str("gaussian"), gaussian() }}}[kernel];
    typename pythonic::assignable_noescape<decltype(pythonic::operator_::mul(y, epsilon))>::type yeps = pythonic::operator_::mul(y, epsilon);
    typename pythonic::assignable_noescape<decltype(pythonic::operator_::mul(x, epsilon))>::type xeps = pythonic::operator_::mul(x, epsilon);
    typename pythonic::assignable_noescape<decltype(pythonic::operator_::div(pythonic::operator_::sub(x, shift), scale))>::type xhat = pythonic::operator_::div(pythonic::operator_::sub(x, shift), scale);
    typename pythonic::assignable<typename __combined<__type43,__type42,__type102>::type>::type out = pythonic::numpy::functor::zeros{}(pythonic::builtins::pythran::functor::make_shape{}(std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, x)), std::get<1>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, coeffs))), pythonic::builtins::functor::float_{});
    typename pythonic::assignable<typename __combined<__type61,__type79,__type80,__type96,__type97>::type>::type vec = pythonic::numpy::functor::empty{}(pythonic::builtins::pythran::functor::make_shape{}(pythonic::operator_::add(p, std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, powers)))), pythonic::builtins::functor::float_{});
    {
      long  __target6476126816 = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, x));
      for (long  i=0L; i < __target6476126816; i += 1L)
      {
        kernel_vector()(xeps[i], yeps, kernel_func, vec[pythonic::types::contiguous_slice(pythonic::builtins::None,p)]);
        polynomial_vector()(xhat[i], powers, vec[pythonic::types::contiguous_slice(p,pythonic::builtins::None)]);
        {
          long  __target6476139152 = std::get<1>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, coeffs));
          for (long  j=0L; j < __target6476139152; j += 1L)
          {
            {
              long  __target6476150624 = pythonic::operator_::add(p, std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, powers)));
              for (long  k=0L; k < __target6476150624; k += 1L)
              {
                out.fast(pythonic::types::make_tuple(i, j)) += pythonic::operator_::mul(coeffs.fast(pythonic::types::make_tuple(k, j)), vec.fast(k));
              }
            }
          }
        }
      }
    }
    return out;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 >
  typename _build_system::type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5>::result_type _build_system::operator()(argument_type0&& y, argument_type1&& d, argument_type2&& smoothing, argument_type3&& kernel, argument_type4&& epsilon, argument_type5&& powers) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::max{})>::type>::type __type0;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type1;
    typedef long __type2;
    typedef typename pythonic::assignable<decltype(std::declval<__type0>()(std::declval<__type1>(), std::declval<__type2>()))>::type __type3;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::min{})>::type>::type __type4;
    typedef typename pythonic::assignable<decltype(std::declval<__type4>()(std::declval<__type1>(), std::declval<__type2>()))>::type __type6;
    typedef decltype(pythonic::operator_::sub(std::declval<__type3>(), std::declval<__type6>())) __type7;
    typedef typename pythonic::assignable<decltype(pythonic::operator_::div(std::declval<__type7>(), std::declval<__type2>()))>::type __type8;
    typedef double __type9;
    typedef container<typename std::remove_reference<__type9>::type> __type10;
    typedef typename __combined<__type8,__type10>::type __type11;
    typedef decltype(pythonic::operator_::eq(std::declval<__type11>(), std::declval<__type9>())) __type12;
    typedef indexable<__type12> __type13;
    typedef typename __combined<__type8,__type13>::type __type14;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::empty{})>::type>::type __type15;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::pythran::functor::make_shape{})>::type>::type __type16;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type17;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type17>())) __type18;
    typedef typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type18>::type>::type>::type __type19;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type __type20;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type20>())) __type21;
    typedef typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type21>::type>::type>::type __type22;
    typedef decltype(pythonic::operator_::add(std::declval<__type19>(), std::declval<__type22>())) __type23;
    typedef decltype(std::declval<__type16>()(std::declval<__type23>(), std::declval<__type23>())) __type27;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::float_{})>::type>::type __type28;
    typedef decltype(std::declval<__type15>()(std::declval<__type27>(), std::declval<__type28>())) __type29;
    typedef typename pythonic::assignable<decltype(pythonic::builtins::getattr(pythonic::types::attr::T{}, std::declval<__type29>()))>::type __type30;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type31;
    typedef decltype(std::declval<__type31>()(std::declval<__type19>())) __type33;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type33>::type::iterator>::value_type>::type __type34;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type34>(), std::declval<__type34>())) __type36;
    typedef indexable<__type36> __type37;
    typedef typename __combined<__type30,__type37>::type __type38;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type40;
    typedef decltype(pythonic::operator_::mul(std::declval<__type1>(), std::declval<__type40>())) __type41;
    typedef pythonic::types::str __type42;
    typedef linear __type43;
    typedef pythonic::types::dict<__type42,__type43> __type44;
    typedef thin_plate_spline __type45;
    typedef pythonic::types::dict<__type42,__type45> __type46;
    typedef typename __combined<__type44,__type46>::type __type47;
    typedef cubic __type48;
    typedef pythonic::types::dict<__type42,__type48> __type49;
    typedef typename __combined<__type47,__type49>::type __type50;
    typedef quintic __type51;
    typedef pythonic::types::dict<__type42,__type51> __type52;
    typedef typename __combined<__type50,__type52>::type __type53;
    typedef multiquadric __type54;
    typedef pythonic::types::dict<__type42,__type54> __type55;
    typedef typename __combined<__type53,__type55>::type __type56;
    typedef inverse_multiquadric __type57;
    typedef pythonic::types::dict<__type42,__type57> __type58;
    typedef typename __combined<__type56,__type58>::type __type59;
    typedef inverse_quadratic __type60;
    typedef pythonic::types::dict<__type42,__type60> __type61;
    typedef typename __combined<__type59,__type61>::type __type62;
    typedef gaussian __type63;
    typedef pythonic::types::dict<__type42,__type63> __type64;
    typedef typename __combined<__type62,__type64>::type __type65;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type66;
    typedef decltype(std::declval<__type65>()[std::declval<__type66>()]) __type67;
    typedef pythonic::types::contiguous_slice __type69;
    typedef decltype(std::declval<__type30>()(std::declval<__type69>(), std::declval<__type69>())) __type70;
    typedef typename kernel_matrix::type<__type41, __type67, __type70>::__ptype4 __type71;
    typedef container<typename std::remove_reference<__type71>::type> __type72;
    typedef container<typename std::remove_reference<__type72>::type> __type73;
    typedef typename kernel_matrix::type<__type41, __type67, __type70>::__ptype5 __type74;
    typedef indexable<__type74> __type75;
    typedef container<typename std::remove_reference<__type75>::type> __type76;
    typedef decltype(pythonic::operator_::add(std::declval<__type3>(), std::declval<__type6>())) __type80;
    typedef typename pythonic::assignable<decltype(pythonic::operator_::div(std::declval<__type80>(), std::declval<__type2>()))>::type __type81;
    typedef decltype(pythonic::operator_::sub(std::declval<__type1>(), std::declval<__type81>())) __type82;
    typedef typename __combined<__type14,__type10,__type13>::type __type83;
    typedef decltype(pythonic::operator_::div(std::declval<__type82>(), std::declval<__type83>())) __type84;
    typedef typename __combined<__type30,__type73,__type76>::type __type86;
    typedef decltype(std::declval<__type86>()(std::declval<__type69>(), std::declval<__type69>())) __type87;
    typedef typename polynomial_matrix::type<__type84, __type20, __type87>::__ptype0 __type88;
    typedef container<typename std::remove_reference<__type88>::type> __type89;
    typedef container<typename std::remove_reference<__type89>::type> __type90;
    typedef typename polynomial_matrix::type<__type84, __type20, __type87>::__ptype1 __type91;
    typedef indexable<__type91> __type92;
    typedef container<typename std::remove_reference<__type92>::type> __type93;
    typedef typename __combined<__type30,__type73,__type76,__type90,__type93>::type __type94;
    typedef decltype(std::declval<__type94>()(std::declval<__type69>(), std::declval<__type69>())) __type95;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::T{}, std::declval<__type95>())) __type96;
    typedef container<typename std::remove_reference<__type96>::type> __type97;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type98;
    typedef decltype(std::declval<__type98>()[std::declval<__type34>()]) __type100;
    typedef container<typename std::remove_reference<__type100>::type> __type101;
    typedef typename std::tuple_element<1,typename std::remove_reference<__type18>::type>::type __type104;
    typedef decltype(std::declval<__type16>()(std::declval<__type104>(), std::declval<__type23>())) __type108;
    typedef decltype(std::declval<__type15>()(std::declval<__type108>(), std::declval<__type28>())) __type109;
    typedef typename pythonic::assignable<decltype(pythonic::builtins::getattr(pythonic::types::attr::T{}, std::declval<__type109>()))>::type __type110;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type33>::type::iterator>::value_type>::type>::type i;
    typename pythonic::assignable_noescape<decltype(std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, d)))>::type p = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, d));
    typename pythonic::assignable_noescape<decltype(std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, powers)))>::type r = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, powers));
    typename pythonic::assignable_noescape<decltype(pythonic::numpy::functor::min{}(y, 0L))>::type mins = pythonic::numpy::functor::min{}(y, 0L);
    typename pythonic::assignable_noescape<decltype(pythonic::numpy::functor::max{}(y, 0L))>::type maxs = pythonic::numpy::functor::max{}(y, 0L);
    typename pythonic::assignable_noescape<decltype(pythonic::operator_::div(pythonic::operator_::add(maxs, mins), 2L))>::type shift = pythonic::operator_::div(pythonic::operator_::add(maxs, mins), 2L);
    typename pythonic::assignable<typename __combined<__type14,__type10,__type13>::type>::type scale = pythonic::operator_::div(pythonic::operator_::sub(maxs, mins), 2L);
    scale.fast(pythonic::operator_::eq(scale, 0.0)) = 1.0;
    typename pythonic::assignable<typename __combined<__type38,__type73,__type76,__type90,__type93,__type97,__type10,__type37,__type101>::type>::type lhs = pythonic::builtins::getattr(pythonic::types::attr::T{}, pythonic::numpy::functor::empty{}(pythonic::builtins::pythran::functor::make_shape{}(pythonic::operator_::add(p, r), pythonic::operator_::add(p, r)), pythonic::builtins::functor::float_{}));
    kernel_matrix()(pythonic::operator_::mul(y, epsilon), typename pythonic::assignable<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<typename __combined<pythonic::types::dict<pythonic::types::str,linear>,pythonic::types::dict<pythonic::types::str,thin_plate_spline>>::type,pythonic::types::dict<pythonic::types::str,cubic>>::type,pythonic::types::dict<pythonic::types::str,quintic>>::type,pythonic::types::dict<pythonic::types::str,multiquadric>>::type,pythonic::types::dict<pythonic::types::str,inverse_multiquadric>>::type,pythonic::types::dict<pythonic::types::str,inverse_quadratic>>::type,pythonic::types::dict<pythonic::types::str,gaussian>>::type>::type{{{ pythonic::types::str("linear"), linear() }, { pythonic::types::str("thin_plate_spline"), thin_plate_spline() }, { pythonic::types::str("cubic"), cubic() }, { pythonic::types::str("quintic"), quintic() }, { pythonic::types::str("multiquadric"), multiquadric() }, { pythonic::types::str("inverse_multiquadric"), inverse_multiquadric() }, { pythonic::types::str("inverse_quadratic"), inverse_quadratic() }, { pythonic::types::str("gaussian"), gaussian() }}}[kernel], lhs(pythonic::types::contiguous_slice(pythonic::builtins::None,p),pythonic::types::contiguous_slice(pythonic::builtins::None,p)));
    polynomial_matrix()(pythonic::operator_::div(pythonic::operator_::sub(y, shift), scale), powers, lhs(pythonic::types::contiguous_slice(pythonic::builtins::None,p),pythonic::types::contiguous_slice(p,pythonic::builtins::None)));
    lhs(pythonic::types::contiguous_slice(p,pythonic::builtins::None),pythonic::types::contiguous_slice(pythonic::builtins::None,p)) = pythonic::builtins::getattr(pythonic::types::attr::T{}, lhs(pythonic::types::contiguous_slice(pythonic::builtins::None,p),pythonic::types::contiguous_slice(p,pythonic::builtins::None)));
    lhs(pythonic::types::contiguous_slice(p,pythonic::builtins::None),pythonic::types::contiguous_slice(p,pythonic::builtins::None)) = 0.0;
    {
      long  __target6476080656 = p;
      for (long  i=0L; i < __target6476080656; i += 1L)
      {
        lhs.fast(pythonic::types::make_tuple(i, i)) += smoothing.fast(i);
      }
    }
    typename pythonic::assignable<typename __combined<__type110,__type17,__type9>::type>::type rhs = pythonic::builtins::getattr(pythonic::types::attr::T{}, pythonic::numpy::functor::empty{}(pythonic::builtins::pythran::functor::make_shape{}(std::get<1>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, d)), pythonic::operator_::add(p, r)), pythonic::builtins::functor::float_{}));
    rhs[pythonic::types::contiguous_slice(pythonic::builtins::None,p)] = d;
    rhs[pythonic::types::contiguous_slice(p,pythonic::builtins::None)] = 0.0;
    return pythonic::types::make_tuple(lhs, rhs, shift, scale);
  }
}
#include <pythonic/python/exception_handler.hpp>
#ifdef ENABLE_PYTHON_MODULE
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>::result_type _evaluate0(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& x, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>::result_type _evaluate1(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& x, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>::result_type _evaluate2(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& x, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>::result_type _evaluate3(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& x, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>::result_type _evaluate4(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& x, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>::result_type _evaluate5(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& x, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>::result_type _evaluate6(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& x, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>::result_type _evaluate7(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& x, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>::result_type _evaluate8(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& x, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>::result_type _evaluate9(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& x, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>::result_type _evaluate10(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& x, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>::result_type _evaluate11(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& x, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>::result_type _evaluate12(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& x, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>::result_type _evaluate13(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& x, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>::result_type _evaluate14(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& x, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_evaluate::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>::result_type _evaluate15(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& x, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& shift, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& scale, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& coeffs) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_evaluate()(x, y, kernel, epsilon, powers, shift, scale, coeffs);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_build_system::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>::result_type _build_system0(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& d, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& smoothing, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_build_system()(y, d, smoothing, kernel, epsilon, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_build_system::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>::result_type _build_system1(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& d, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& smoothing, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_build_system()(y, d, smoothing, kernel, epsilon, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_build_system::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>::result_type _build_system2(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& d, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& smoothing, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_build_system()(y, d, smoothing, kernel, epsilon, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_build_system::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>::result_type _build_system3(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& y, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& d, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& smoothing, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_build_system()(y, d, smoothing, kernel, epsilon, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_build_system::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>::result_type _build_system4(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& d, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& smoothing, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_build_system()(y, d, smoothing, kernel, epsilon, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_build_system::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>::result_type _build_system5(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& d, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& smoothing, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_build_system()(y, d, smoothing, kernel, epsilon, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_build_system::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::str, double, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>::result_type _build_system6(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& d, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& smoothing, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_build_system()(y, d, smoothing, kernel, epsilon, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_build_system::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::str, double, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>::result_type _build_system7(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& y, pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& d, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& smoothing, pythonic::types::str&& kernel, double&& epsilon, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_build_system()(y, d, smoothing, kernel, epsilon, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_polynomial_matrix::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>::result_type _polynomial_matrix0(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& x, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_polynomial_matrix()(x, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_polynomial_matrix::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>::result_type _polynomial_matrix1(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& x, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_polynomial_matrix()(x, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_polynomial_matrix::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>::result_type _polynomial_matrix2(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& x, pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_polynomial_matrix()(x, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_polynomial_matrix::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>::result_type _polynomial_matrix3(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& x, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& powers) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_polynomial_matrix()(x, powers);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_kernel_matrix::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, pythonic::types::str>::result_type _kernel_matrix0(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& x, pythonic::types::str&& kernel) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_kernel_matrix()(x, kernel);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
typename __pythran__rbfinterp_pythran::_kernel_matrix::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, pythonic::types::str>::result_type _kernel_matrix1(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& x, pythonic::types::str&& kernel) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__rbfinterp_pythran::_kernel_matrix()(x, kernel);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}

static PyObject *
__pythran_wrap__evaluate0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7]))
        return to_python(_evaluate0(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate1(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7]))
        return to_python(_evaluate1(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate2(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7]))
        return to_python(_evaluate2(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate3(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7]))
        return to_python(_evaluate3(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate4(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7]))
        return to_python(_evaluate4(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate5(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7]))
        return to_python(_evaluate5(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate6(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7]))
        return to_python(_evaluate6(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate7(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7]))
        return to_python(_evaluate7(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate8(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7]))
        return to_python(_evaluate8(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate9(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7]))
        return to_python(_evaluate9(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate10(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7]))
        return to_python(_evaluate10(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate11(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7]))
        return to_python(_evaluate11(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate12(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7]))
        return to_python(_evaluate12(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate13(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7]))
        return to_python(_evaluate13(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate14(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7]))
        return to_python(_evaluate14(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__evaluate15(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[8+1];
    char const* keywords[] = {"x", "y", "kernel", "epsilon", "powers", "shift", "scale", "coeffs",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5], &args_obj[6], &args_obj[7]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]) && is_convertible<double>(args_obj[3]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7]))
        return to_python(_evaluate15(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2]), from_python<double>(args_obj[3]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[4]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[5]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[6]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[7])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__build_system0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[6+1];
    char const* keywords[] = {"y", "d", "smoothing", "kernel", "epsilon", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]) && is_convertible<pythonic::types::str>(args_obj[3]) && is_convertible<double>(args_obj[4]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[5]))
        return to_python(_build_system0(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]), from_python<pythonic::types::str>(args_obj[3]), from_python<double>(args_obj[4]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[5])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__build_system1(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[6+1];
    char const* keywords[] = {"y", "d", "smoothing", "kernel", "epsilon", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]) && is_convertible<pythonic::types::str>(args_obj[3]) && is_convertible<double>(args_obj[4]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[5]))
        return to_python(_build_system1(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]), from_python<pythonic::types::str>(args_obj[3]), from_python<double>(args_obj[4]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[5])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__build_system2(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[6+1];
    char const* keywords[] = {"y", "d", "smoothing", "kernel", "epsilon", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]) && is_convertible<pythonic::types::str>(args_obj[3]) && is_convertible<double>(args_obj[4]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[5]))
        return to_python(_build_system2(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]), from_python<pythonic::types::str>(args_obj[3]), from_python<double>(args_obj[4]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[5])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__build_system3(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[6+1];
    char const* keywords[] = {"y", "d", "smoothing", "kernel", "epsilon", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]) && is_convertible<pythonic::types::str>(args_obj[3]) && is_convertible<double>(args_obj[4]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[5]))
        return to_python(_build_system3(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]), from_python<pythonic::types::str>(args_obj[3]), from_python<double>(args_obj[4]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[5])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__build_system4(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[6+1];
    char const* keywords[] = {"y", "d", "smoothing", "kernel", "epsilon", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]) && is_convertible<pythonic::types::str>(args_obj[3]) && is_convertible<double>(args_obj[4]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[5]))
        return to_python(_build_system4(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]), from_python<pythonic::types::str>(args_obj[3]), from_python<double>(args_obj[4]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[5])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__build_system5(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[6+1];
    char const* keywords[] = {"y", "d", "smoothing", "kernel", "epsilon", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]) && is_convertible<pythonic::types::str>(args_obj[3]) && is_convertible<double>(args_obj[4]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[5]))
        return to_python(_build_system5(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[1]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]), from_python<pythonic::types::str>(args_obj[3]), from_python<double>(args_obj[4]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[5])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__build_system6(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[6+1];
    char const* keywords[] = {"y", "d", "smoothing", "kernel", "epsilon", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]) && is_convertible<pythonic::types::str>(args_obj[3]) && is_convertible<double>(args_obj[4]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[5]))
        return to_python(_build_system6(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]), from_python<pythonic::types::str>(args_obj[3]), from_python<double>(args_obj[4]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[5])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__build_system7(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[6+1];
    char const* keywords[] = {"y", "d", "smoothing", "kernel", "epsilon", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3], &args_obj[4], &args_obj[5]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]) && is_convertible<pythonic::types::str>(args_obj[3]) && is_convertible<double>(args_obj[4]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[5]))
        return to_python(_build_system7(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[1]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[2]), from_python<pythonic::types::str>(args_obj[3]), from_python<double>(args_obj[4]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[5])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__polynomial_matrix0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[2+1];
    char const* keywords[] = {"x", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OO",
                                     (char**)keywords , &args_obj[0], &args_obj[1]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[1]))
        return to_python(_polynomial_matrix0(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[1])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__polynomial_matrix1(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[2+1];
    char const* keywords[] = {"x", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OO",
                                     (char**)keywords , &args_obj[0], &args_obj[1]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[1]))
        return to_python(_polynomial_matrix1(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[1])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__polynomial_matrix2(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[2+1];
    char const* keywords[] = {"x", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OO",
                                     (char**)keywords , &args_obj[0], &args_obj[1]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[1]))
        return to_python(_polynomial_matrix2(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[1])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__polynomial_matrix3(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[2+1];
    char const* keywords[] = {"x", "powers",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OO",
                                     (char**)keywords , &args_obj[0], &args_obj[1]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[1]))
        return to_python(_polynomial_matrix3(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[1])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__kernel_matrix0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[2+1];
    char const* keywords[] = {"x", "kernel",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OO",
                                     (char**)keywords , &args_obj[0], &args_obj[1]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<pythonic::types::str>(args_obj[1]))
        return to_python(_kernel_matrix0(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<pythonic::types::str>(args_obj[1])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__kernel_matrix1(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[2+1];
    char const* keywords[] = {"x", "kernel",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OO",
                                     (char**)keywords , &args_obj[0], &args_obj[1]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<pythonic::types::str>(args_obj[1]))
        return to_python(_kernel_matrix1(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<pythonic::types::str>(args_obj[1])));
    else {
        return nullptr;
    }
}

            static PyObject *
            __pythran_wrapall__evaluate(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__evaluate0(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate1(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate2(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate3(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate4(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate5(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate6(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate7(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate8(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate9(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate10(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate11(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate12(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate13(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate14(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__evaluate15(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_evaluate", "\n""    - _evaluate(float[:,:], float[:,:], str, float, int[:,:], float[:], float[:], float[:,:])", args, kw);
                });
            }


            static PyObject *
            __pythran_wrapall__build_system(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__build_system0(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__build_system1(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__build_system2(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__build_system3(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__build_system4(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__build_system5(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__build_system6(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__build_system7(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_build_system", "\n""    - _build_system(float[:,:], float[:,:], float[:], str, float, int[:,:])", args, kw);
                });
            }


            static PyObject *
            __pythran_wrapall__polynomial_matrix(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__polynomial_matrix0(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__polynomial_matrix1(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__polynomial_matrix2(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__polynomial_matrix3(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_polynomial_matrix", "\n""    - _polynomial_matrix(float[:,:], int[:,:])", args, kw);
                });
            }


            static PyObject *
            __pythran_wrapall__kernel_matrix(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__kernel_matrix0(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__kernel_matrix1(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_kernel_matrix", "\n""    - _kernel_matrix(float[:,:], str)", args, kw);
                });
            }


static PyMethodDef Methods[] = {
    {
    "_evaluate",
    (PyCFunction)__pythran_wrapall__evaluate,
    METH_VARARGS | METH_KEYWORDS,
    "Evaluate the RBF interpolant at `x`.\n""\n""    Supported prototypes:\n""\n""    - _evaluate(float[:,:], float[:,:], str, float, int[:,:], float[:], float[:], float[:,:])\n""\n""    Parameters\n""    ----------\n""    x : (Q, N) float ndarray\n""        Evaluation point coordinates.\n""    y : (P, N) float ndarray\n""        Data point coordinates.\n""    kernel : str\n""        Name of the RBF.\n""    epsilon : float\n""        Shape parameter.\n""    powers : (R, N) int ndarray\n""        The exponents for each monomial in the polynomial.\n""    shift : (N,) float ndarray\n""        Shifts the polynomial domain for numerical stability.\n""    scale : (N,) float ndarray\n""        Scales the polynomial domain for numerical stability.\n""    coeffs : (P + R, S) float ndarray\n""        Coefficients for each RBF and monomial.\n""\n""    Returns\n""    -------\n""    (Q, S) float ndarray\n""\n"""},{
    "_build_system",
    (PyCFunction)__pythran_wrapall__build_system,
    METH_VARARGS | METH_KEYWORDS,
    "Build the system used to solve for the RBF interpolant coefficients.\n""\n""    Supported prototypes:\n""\n""    - _build_system(float[:,:], float[:,:], float[:], str, float, int[:,:])\n""\n""    Parameters\n""    ----------\n""    y : (P, N) float ndarray\n""        Data point coordinates.\n""    d : (P, S) float ndarray\n""        Data values at `y`.\n""    smoothing : (P,) float ndarray\n""        Smoothing parameter for each data point.\n""    kernel : str\n""        Name of the RBF.\n""    epsilon : float\n""        Shape parameter.\n""    powers : (R, N) int ndarray\n""        The exponents for each monomial in the polynomial.\n""\n""    Returns\n""    -------\n""    lhs : (P + R, P + R) float ndarray\n""        Left-hand side matrix.\n""    rhs : (P + R, S) float ndarray\n""        Right-hand side matrix.\n""    shift : (N,) float ndarray\n""        Domain shift used to create the polynomial matrix.\n""    scale : (N,) float ndarray\n""        Domain scaling used to create the polynomial matrix.\n""\n"""},{
    "_polynomial_matrix",
    (PyCFunction)__pythran_wrapall__polynomial_matrix,
    METH_VARARGS | METH_KEYWORDS,
    "Return monomials, with exponents from `powers`, evaluated at `x`.\n""\n""    Supported prototypes:\n""\n""    - _polynomial_matrix(float[:,:], int[:,:])"},{
    "_kernel_matrix",
    (PyCFunction)__pythran_wrapall__kernel_matrix,
    METH_VARARGS | METH_KEYWORDS,
    "Return RBFs, with centers at `x`, evaluated at `x`.\n""\n""    Supported prototypes:\n""\n""    - _kernel_matrix(float[:,:], str)"},
    {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_rbfinterp_pythran",            /* m_name */
    "",         /* m_doc */
    -1,                  /* m_size */
    Methods,             /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL,                /* m_free */
  };
#define PYTHRAN_RETURN return theModule
#define PYTHRAN_MODULE_INIT(s) PyInit_##s
#else
#define PYTHRAN_RETURN return
#define PYTHRAN_MODULE_INIT(s) init##s
#endif
PyMODINIT_FUNC
PYTHRAN_MODULE_INIT(_rbfinterp_pythran)(void)
#ifndef _WIN32
__attribute__ ((visibility("default")))
__attribute__ ((externally_visible))
#endif
;
PyMODINIT_FUNC
PYTHRAN_MODULE_INIT(_rbfinterp_pythran)(void) {
    import_array()
    #if PY_MAJOR_VERSION >= 3
    PyObject* theModule = PyModule_Create(&moduledef);
    #else
    PyObject* theModule = Py_InitModule3("_rbfinterp_pythran",
                                         Methods,
                                         ""
    );
    #endif
    if(! theModule)
        PYTHRAN_RETURN;
    PyObject * theDoc = Py_BuildValue("(sss)",
                                      "0.9.11",
                                      "2021-06-04 16:00:41.042263",
                                      "3c30425550c4548ade4c98d9f66ed93a241515ad72e7efe308ab023945aca246");
    if(! theDoc)
        PYTHRAN_RETURN;
    PyModule_AddObject(theModule,
                       "__pythran__",
                       theDoc);


    PYTHRAN_RETURN;
}

#endif