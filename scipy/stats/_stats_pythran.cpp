#include <pythonic/core.hpp>
#include <pythonic/python/core.hpp>
#include <pythonic/types/bool.hpp>
#include <pythonic/types/int.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <pythonic\include\types\int64.hpp>
#include <pythonic\include\types\float32.hpp>
#include <pythonic\include\types\str.hpp>
#include <pythonic\include\types\float64.hpp>
#include <pythonic\include\types\float.hpp>
#include <pythonic\include\types\ndarray.hpp>
#include <pythonic\include\types\int.hpp>
#include <pythonic\include\types\numpy_texpr.hpp>
#include <pythonic\types\float.hpp>
#include <pythonic\types\float32.hpp>
#include <pythonic\types\numpy_texpr.hpp>
#include <pythonic\types\ndarray.hpp>
#include <pythonic\types\float64.hpp>
#include <pythonic\types\int.hpp>
#include <pythonic\types\int64.hpp>
#include <pythonic\types\str.hpp>
#include <pythonic/include/builtins/getattr.hpp>
#include <pythonic/include/builtins/int_.hpp>
#include <pythonic/include/builtins/len.hpp>
#include <pythonic/include/builtins/list.hpp>
#include <pythonic/include/builtins/list/append.hpp>
#include <pythonic/include/builtins/map.hpp>
#include <pythonic/include/builtins/max.hpp>
#include <pythonic/include/builtins/min.hpp>
#include <pythonic/include/builtins/pythran/static_list.hpp>
#include <pythonic/include/builtins/range.hpp>
#include <pythonic/include/builtins/tuple.hpp>
#include <pythonic/include/functools/partial.hpp>
#include <pythonic/include/numpy/asarray.hpp>
#include <pythonic/include/numpy/ceil.hpp>
#include <pythonic/include/numpy/expand_dims.hpp>
#include <pythonic/include/numpy/float64.hpp>
#include <pythonic/include/numpy/floor.hpp>
#include <pythonic/include/numpy/median.hpp>
#include <pythonic/include/numpy/nonzero.hpp>
#include <pythonic/include/numpy/ones.hpp>
#include <pythonic/include/numpy/square.hpp>
#include <pythonic/include/numpy/sum.hpp>
#include <pythonic/include/numpy/zeros.hpp>
#include <pythonic/include/operator_/add.hpp>
#include <pythonic/include/operator_/div.hpp>
#include <pythonic/include/operator_/eq.hpp>
#include <pythonic/include/operator_/floordiv.hpp>
#include <pythonic/include/operator_/ge.hpp>
#include <pythonic/include/operator_/gt.hpp>
#include <pythonic/include/operator_/iadd.hpp>
#include <pythonic/include/operator_/isub.hpp>
#include <pythonic/include/operator_/le.hpp>
#include <pythonic/include/operator_/lt.hpp>
#include <pythonic/include/operator_/mul.hpp>
#include <pythonic/include/operator_/sub.hpp>
#include <pythonic/include/scipy/special/binom.hpp>
#include <pythonic/include/types/slice.hpp>
#include <pythonic/include/types/str.hpp>
#include <pythonic/builtins/getattr.hpp>
#include <pythonic/builtins/int_.hpp>
#include <pythonic/builtins/len.hpp>
#include <pythonic/builtins/list.hpp>
#include <pythonic/builtins/list/append.hpp>
#include <pythonic/builtins/map.hpp>
#include <pythonic/builtins/max.hpp>
#include <pythonic/builtins/min.hpp>
#include <pythonic/builtins/pythran/static_list.hpp>
#include <pythonic/builtins/range.hpp>
#include <pythonic/builtins/tuple.hpp>
#include <pythonic/functools/partial.hpp>
#include <pythonic/numpy/asarray.hpp>
#include <pythonic/numpy/ceil.hpp>
#include <pythonic/numpy/expand_dims.hpp>
#include <pythonic/numpy/float64.hpp>
#include <pythonic/numpy/floor.hpp>
#include <pythonic/numpy/median.hpp>
#include <pythonic/numpy/nonzero.hpp>
#include <pythonic/numpy/ones.hpp>
#include <pythonic/numpy/square.hpp>
#include <pythonic/numpy/sum.hpp>
#include <pythonic/numpy/zeros.hpp>
#include <pythonic/operator_/add.hpp>
#include <pythonic/operator_/div.hpp>
#include <pythonic/operator_/eq.hpp>
#include <pythonic/operator_/floordiv.hpp>
#include <pythonic/operator_/ge.hpp>
#include <pythonic/operator_/gt.hpp>
#include <pythonic/operator_/iadd.hpp>
#include <pythonic/operator_/isub.hpp>
#include <pythonic/operator_/le.hpp>
#include <pythonic/operator_/lt.hpp>
#include <pythonic/operator_/mul.hpp>
#include <pythonic/operator_/sub.hpp>
#include <pythonic/scipy/special/binom.hpp>
#include <pythonic/types/slice.hpp>
#include <pythonic/types/str.hpp>
namespace __pythran__stats_pythran
{
  struct _count_paths_outside_method_lambda0
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
      typedef __type0 __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type2;
      typedef __type2 __type3;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type4;
      typedef __type4 __type5;
      typedef decltype(pythonic::operator_::mul(std::declval<__type3>(), std::declval<__type5>())) __type6;
      typedef decltype(pythonic::operator_::add(std::declval<__type1>(), std::declval<__type6>())) __type7;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type8;
      typedef __type8 __type9;
      typedef decltype(pythonic::operator_::add(std::declval<__type7>(), std::declval<__type9>())) __type10;
      typedef long __type11;
      typedef decltype(pythonic::operator_::sub(std::declval<__type10>(), std::declval<__type11>())) __type12;
      typedef decltype(pythonic::operator_::functor::floordiv()(std::declval<__type12>(), std::declval<__type9>())) __type14;
      typedef typename pythonic::returnable<__type14>::type __type15;
      typedef __type15 result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    inline
    typename type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type operator()(argument_type0&& h, argument_type1&& mg, argument_type2&& ng, argument_type3&& j) const
    ;
  }  ;
  struct siegelslopes
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::median{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::asarray{})>::type>::type __type1;
      typedef pythonic::types::empty_list __type2;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type2>(), std::declval<__type2>())) __type3;
      typedef typename pythonic::assignable<__type3>::type __type4;
      typedef __type4 __type5;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type5>::type>::type __type6;
      typedef typename pythonic::assignable<__type6>::type __type7;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::expand_dims{})>::type>::type __type8;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type9;
      typedef __type9 __type10;
      typedef long __type11;
      typedef decltype(std::declval<__type8>()(std::declval<__type10>(), std::declval<__type11>())) __type12;
      typedef decltype(pythonic::operator_::sub(std::declval<__type12>(), std::declval<__type10>())) __type14;
      typedef typename pythonic::assignable<__type14>::type __type15;
      typedef __type15 __type16;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type17;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::len{})>::type>::type __type18;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type19;
      typedef __type19 __type20;
      typedef decltype(std::declval<__type18>()(std::declval<__type20>())) __type21;
      typedef decltype(std::declval<__type17>()(std::declval<__type21>())) __type22;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type22>::type::iterator>::value_type>::type __type23;
      typedef __type23 __type24;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::nonzero{})>::type>::type __type25;
      typedef decltype(std::declval<__type8>()(std::declval<__type20>(), std::declval<__type11>())) __type27;
      typedef decltype(pythonic::operator_::sub(std::declval<__type27>(), std::declval<__type20>())) __type29;
      typedef typename pythonic::assignable<__type29>::type __type30;
      typedef __type30 __type31;
      typedef pythonic::types::contiguous_slice __type33;
      typedef decltype(std::declval<__type31>()(std::declval<__type24>(), std::declval<__type33>())) __type34;
      typedef decltype(std::declval<__type25>()(std::declval<__type34>())) __type35;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type35>::type>::type __type36;
      typedef typename pythonic::assignable<__type36>::type __type37;
      typedef __type37 __type38;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type24>(), std::declval<__type38>())) __type39;
      typedef decltype(std::declval<__type16>()[std::declval<__type39>()]) __type40;
      typedef decltype(std::declval<__type31>()[std::declval<__type39>()]) __type45;
      typedef decltype(pythonic::operator_::div(std::declval<__type40>(), std::declval<__type45>())) __type46;
      typedef decltype(std::declval<__type0>()(std::declval<__type46>())) __type47;
      typedef pythonic::types::list<typename std::remove_reference<__type47>::type> __type48;
      typedef __type7 __type49;
      typedef typename __combined<__type7,__type48,__type49>::type __type50;
      typedef __type50 __type51;
      typedef decltype(std::declval<__type1>()(std::declval<__type51>())) __type52;
      typedef decltype(std::declval<__type0>()(std::declval<__type52>())) __type53;
      typedef typename pythonic::assignable<__type53>::type __type54;
      typedef __type54 __type55;
      typedef typename std::tuple_element<1,typename std::remove_reference<__type5>::type>::type __type56;
      typedef typename pythonic::assignable<__type56>::type __type57;
      typedef decltype(std::declval<__type20>()[std::declval<__type24>()]) __type61;
      typedef decltype(pythonic::operator_::mul(std::declval<__type10>(), std::declval<__type61>())) __type62;
      typedef decltype(std::declval<__type10>()[std::declval<__type24>()]) __type65;
      typedef decltype(pythonic::operator_::mul(std::declval<__type65>(), std::declval<__type20>())) __type67;
      typedef decltype(pythonic::operator_::sub(std::declval<__type62>(), std::declval<__type67>())) __type68;
      typedef decltype(std::declval<__type68>()[std::declval<__type38>()]) __type70;
      typedef decltype(pythonic::operator_::div(std::declval<__type70>(), std::declval<__type45>())) __type76;
      typedef decltype(std::declval<__type0>()(std::declval<__type76>())) __type77;
      typedef pythonic::types::list<typename std::remove_reference<__type77>::type> __type78;
      typedef __type57 __type79;
      typedef typename __combined<__type57,__type78,__type79>::type __type80;
      typedef __type80 __type81;
      typedef decltype(std::declval<__type1>()(std::declval<__type81>())) __type82;
      typedef decltype(std::declval<__type0>()(std::declval<__type82>())) __type83;
      typedef typename pythonic::lazy<__type83>::type __type84;
      typedef decltype(pythonic::operator_::mul(std::declval<__type55>(), std::declval<__type20>())) __type88;
      typedef decltype(pythonic::operator_::sub(std::declval<__type10>(), std::declval<__type88>())) __type89;
      typedef decltype(std::declval<__type0>()(std::declval<__type89>())) __type90;
      typedef typename pythonic::lazy<__type90>::type __type91;
      typedef typename __combined<__type84,__type91>::type __type92;
      typedef __type92 __type93;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type55>(), std::declval<__type93>())) __type94;
      typedef typename pythonic::returnable<__type94>::type __type95;
      typedef __type95 result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    inline
    typename type<argument_type0, argument_type1, argument_type2>::result_type operator()(argument_type0&& y, argument_type1&& x, argument_type2&& method) const
    ;
  }  ;
  struct _compute_prob_outside_square
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 >
    struct type
    {
      typedef long __type0;
      typedef double __type1;
      typedef typename pythonic::lazy<__type1>::type __type2;
      typedef typename pythonic::assignable<__type1>::type __type3;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type4;
      typedef __type4 __type5;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::int_{})>::type>::type __type6;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::floor{})>::type>::type __type7;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type9;
      typedef __type9 __type10;
      typedef decltype(pythonic::operator_::div(std::declval<__type5>(), std::declval<__type10>())) __type11;
      typedef decltype(std::declval<__type7>()(std::declval<__type11>())) __type12;
      typedef decltype(std::declval<__type6>()(std::declval<__type12>())) __type13;
      typedef typename pythonic::assignable<__type13>::type __type14;
      typedef __type14 __type15;
      typedef decltype(pythonic::operator_::mul(std::declval<__type15>(), std::declval<__type10>())) __type17;
      typedef decltype(pythonic::operator_::sub(std::declval<__type5>(), std::declval<__type17>())) __type18;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type19;
      typedef decltype(std::declval<__type19>()(std::declval<__type10>())) __type21;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type21>::type::iterator>::value_type>::type __type22;
      typedef __type22 __type23;
      typedef decltype(pythonic::operator_::sub(std::declval<__type18>(), std::declval<__type23>())) __type24;
      typedef __type3 __type25;
      typedef decltype(pythonic::operator_::mul(std::declval<__type24>(), std::declval<__type25>())) __type26;
      typedef decltype(pythonic::operator_::add(std::declval<__type5>(), std::declval<__type17>())) __type31;
      typedef decltype(pythonic::operator_::add(std::declval<__type31>(), std::declval<__type23>())) __type33;
      typedef decltype(pythonic::operator_::add(std::declval<__type33>(), std::declval<__type0>())) __type34;
      typedef decltype(pythonic::operator_::div(std::declval<__type26>(), std::declval<__type34>())) __type35;
      typedef typename pythonic::assignable<__type35>::type __type36;
      typedef typename __combined<__type3,__type36>::type __type37;
      typedef __type37 __type38;
      typedef __type2 __type39;
      typedef decltype(pythonic::operator_::sub(std::declval<__type1>(), std::declval<__type39>())) __type40;
      typedef decltype(pythonic::operator_::mul(std::declval<__type38>(), std::declval<__type40>())) __type41;
      typedef typename pythonic::lazy<__type41>::type __type42;
      typedef typename __combined<__type2,__type42>::type __type43;
      typedef __type43 __type44;
      typedef decltype(pythonic::operator_::mul(std::declval<__type0>(), std::declval<__type44>())) __type45;
      typedef typename pythonic::returnable<__type45>::type __type46;
      typedef __type46 result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 >
    inline
    typename type<argument_type0, argument_type1>::result_type operator()(argument_type0&& n, argument_type1&& h) const
    ;
  }  ;
  struct _compute_outer_prob_inside_method
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type0;
      typedef __type0 __type1;
      typedef __type1 __type2;
      typedef typename pythonic::assignable<__type1>::type __type3;
      typedef __type3 __type4;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type5;
      typedef __type5 __type6;
      typedef typename pythonic::assignable<__type6>::type __type7;
      typedef __type7 __type8;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type4>(), std::declval<__type8>())) __type9;
      typedef typename pythonic::assignable<__type9>::type __type10;
      typedef __type10 __type11;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type11>::type>::type __type12;
      typedef __type12 __type13;
      typedef __type6 __type14;
      typedef typename std::tuple_element<1,typename std::remove_reference<__type11>::type>::type __type15;
      typedef __type15 __type16;
      typedef double __type17;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::ones{})>::type>::type __type18;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::min{})>::type>::type __type19;
      typedef long __type20;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::int_{})>::type>::type __type21;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::ceil{})>::type>::type __type22;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type23;
      typedef __type23 __type24;
      typedef typename pythonic::assignable<__type12>::type __type25;
      typedef typename __combined<__type7,__type25>::type __type26;
      typedef __type26 __type27;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type28;
      typedef __type28 __type29;
      typedef decltype(pythonic::operator_::functor::floordiv()(std::declval<__type27>(), std::declval<__type29>())) __type30;
      typedef typename pythonic::assignable<__type30>::type __type31;
      typedef __type31 __type32;
      typedef decltype(pythonic::operator_::div(std::declval<__type24>(), std::declval<__type32>())) __type33;
      typedef decltype(std::declval<__type22>()(std::declval<__type33>())) __type34;
      typedef decltype(std::declval<__type21>()(std::declval<__type34>())) __type35;
      typedef typename pythonic::assignable<__type15>::type __type36;
      typedef typename __combined<__type3,__type36>::type __type37;
      typedef __type37 __type38;
      typedef decltype(pythonic::operator_::add(std::declval<__type38>(), std::declval<__type20>())) __type39;
      typedef decltype(std::declval<__type19>()(std::declval<__type35>(), std::declval<__type39>())) __type40;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type20>(), std::declval<__type40>())) __type41;
      typedef typename pythonic::assignable<__type41>::type __type42;
      typedef __type42 __type43;
      typedef typename std::tuple_element<1,typename std::remove_reference<__type43>::type>::type __type44;
      typedef typename pythonic::assignable<__type44>::type __type45;
      typedef __type45 __type46;
      typedef decltype(pythonic::operator_::mul(std::declval<__type20>(), std::declval<__type46>())) __type47;
      typedef decltype(pythonic::operator_::add(std::declval<__type47>(), std::declval<__type20>())) __type48;
      typedef decltype(std::declval<__type19>()(std::declval<__type48>(), std::declval<__type39>())) __type51;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::float64{})>::type>::type __type52;
      typedef decltype(std::declval<__type18>()(std::declval<__type51>(), std::declval<__type52>())) __type53;
      typedef typename pythonic::assignable<__type53>::type __type54;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type55;
      typedef decltype(pythonic::operator_::functor::floordiv()(std::declval<__type38>(), std::declval<__type29>())) __type58;
      typedef typename pythonic::assignable<__type58>::type __type59;
      typedef __type59 __type60;
      typedef decltype(pythonic::operator_::add(std::declval<__type27>(), std::declval<__type20>())) __type62;
      typedef decltype(std::declval<__type55>()(std::declval<__type20>(), std::declval<__type62>())) __type63;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type63>::type::iterator>::value_type>::type __type64;
      typedef __type64 __type65;
      typedef decltype(pythonic::operator_::mul(std::declval<__type60>(), std::declval<__type65>())) __type66;
      typedef decltype(pythonic::operator_::add(std::declval<__type66>(), std::declval<__type24>())) __type68;
      typedef decltype(pythonic::operator_::div(std::declval<__type68>(), std::declval<__type32>())) __type70;
      typedef decltype(std::declval<__type22>()(std::declval<__type70>())) __type71;
      typedef decltype(std::declval<__type21>()(std::declval<__type71>())) __type72;
      typedef decltype(std::declval<__type19>()(std::declval<__type72>(), std::declval<__type39>())) __type75;
      typedef typename pythonic::assignable<__type75>::type __type76;
      typedef typename __combined<__type45,__type76>::type __type77;
      typedef __type77 __type78;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type43>::type>::type __type79;
      typedef typename pythonic::assignable<__type79>::type __type80;
      typedef __type80 __type81;
      typedef decltype(pythonic::operator_::sub(std::declval<__type46>(), std::declval<__type81>())) __type84;
      typedef typename pythonic::assignable<__type84>::type __type85;
      typedef __type85 __type86;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type81>(), std::declval<__type86>())) __type87;
      typedef typename pythonic::assignable<__type87>::type __type88;
      typedef __type88 __type89;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type89>::type>::type __type90;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::max{})>::type>::type __type91;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::floor{})>::type>::type __type92;
      typedef decltype(pythonic::operator_::sub(std::declval<__type66>(), std::declval<__type24>())) __type97;
      typedef decltype(pythonic::operator_::div(std::declval<__type97>(), std::declval<__type32>())) __type99;
      typedef decltype(std::declval<__type92>()(std::declval<__type99>())) __type100;
      typedef decltype(std::declval<__type21>()(std::declval<__type100>())) __type101;
      typedef decltype(pythonic::operator_::add(std::declval<__type101>(), std::declval<__type20>())) __type102;
      typedef typename __combined<__type102,__type20>::type __type103;
      typedef decltype(std::declval<__type91>()(std::declval<__type103>(), std::declval<__type20>())) __type104;
      typedef decltype(std::declval<__type19>()(std::declval<__type104>(), std::declval<__type38>())) __type106;
      typedef typename __combined<__type90,__type106>::type __type107;
      typedef typename pythonic::assignable<__type107>::type __type108;
      typedef typename __combined<__type80,__type108>::type __type109;
      typedef __type109 __type110;
      typedef decltype(pythonic::operator_::sub(std::declval<__type78>(), std::declval<__type110>())) __type111;
      typedef decltype(std::declval<__type55>()(std::declval<__type111>())) __type112;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type112>::type::iterator>::value_type>::type __type113;
      typedef __type113 __type114;
      typedef indexable<__type114> __type115;
      typedef typename pythonic::assignable<__type17>::type __type116;
      typedef typename __combined<__type54,__type17>::type __type117;
      typedef __type117 __type118;
      typedef decltype(pythonic::operator_::add(std::declval<__type114>(), std::declval<__type110>())) __type121;
      typedef typename pythonic::assignable<__type90>::type __type122;
      typedef __type122 __type123;
      typedef decltype(pythonic::operator_::sub(std::declval<__type121>(), std::declval<__type123>())) __type124;
      typedef decltype(std::declval<__type118>()[std::declval<__type124>()]) __type125;
      typedef decltype(pythonic::operator_::mul(std::declval<__type125>(), std::declval<__type65>())) __type127;
      typedef __type116 __type128;
      typedef typename pythonic::assignable<__type121>::type __type132;
      typedef __type132 __type133;
      typedef decltype(pythonic::operator_::mul(std::declval<__type128>(), std::declval<__type133>())) __type134;
      typedef decltype(pythonic::operator_::add(std::declval<__type127>(), std::declval<__type134>())) __type135;
      typedef decltype(pythonic::operator_::add(std::declval<__type65>(), std::declval<__type133>())) __type138;
      typedef decltype(pythonic::operator_::div(std::declval<__type135>(), std::declval<__type138>())) __type139;
      typedef typename pythonic::assignable<__type139>::type __type140;
      typedef typename __combined<__type116,__type140>::type __type141;
      typedef __type141 __type142;
      typedef container<typename std::remove_reference<__type142>::type> __type143;
      typedef typename __combined<__type54,__type115,__type17,__type143,__type20>::type __type144;
      typedef __type144 __type145;
      typedef decltype(pythonic::operator_::sub(std::declval<__type111>(), std::declval<__type20>())) __type149;
      typedef decltype(std::declval<__type145>()[std::declval<__type149>()]) __type150;
      typedef typename __combined<__type17,__type150>::type __type151;
      typedef typename pythonic::returnable<__type151>::type __type152;
      typedef __type2 __ptype0;
      typedef __type13 __ptype2;
      typedef __type14 __ptype1;
      typedef __type16 __ptype3;
      typedef __type152 result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    inline
    typename type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type operator()(argument_type0&& m, argument_type1&& n, argument_type2&& g, argument_type3&& h) const
    ;
  }  ;
  struct _a_ij_Aij_Dij2
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
      typedef __type0 __type1;
      typedef __type1 __type2;
      typedef long __type5;
      typedef typename pythonic::assignable<__type5>::type __type6;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type8;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type1>())) __type10;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type10>::type>::type __type11;
      typedef typename pythonic::lazy<__type11>::type __type12;
      typedef __type12 __type13;
      typedef decltype(std::declval<__type8>()(std::declval<__type13>())) __type14;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type14>::type::iterator>::value_type>::type __type15;
      typedef __type15 __type16;
      typedef typename std::tuple_element<1,typename std::remove_reference<__type10>::type>::type __type17;
      typedef typename pythonic::lazy<__type17>::type __type18;
      typedef __type18 __type19;
      typedef decltype(std::declval<__type8>()(std::declval<__type19>())) __type20;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type20>::type::iterator>::value_type>::type __type21;
      typedef __type21 __type22;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type16>(), std::declval<__type22>())) __type23;
      typedef decltype(std::declval<__type1>()[std::declval<__type23>()]) __type24;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::square{})>::type>::type __type25;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::sum{})>::type>::type __type26;
      typedef typename pythonic::assignable<__type1>::type __type27;
      typedef typename __combined<__type27,__type1>::type __type28;
      typedef __type28 __type29;
      typedef pythonic::types::contiguous_slice __type30;
      typedef decltype(std::declval<__type29>()(std::declval<__type30>(), std::declval<__type30>())) __type31;
      typedef decltype(std::declval<__type26>()(std::declval<__type31>())) __type32;
      typedef decltype(pythonic::operator_::add(std::declval<__type32>(), std::declval<__type32>())) __type35;
      typedef typename pythonic::assignable<__type2>::type __type37;
      typedef __type37 __type38;
      typedef decltype(std::declval<__type38>()(std::declval<__type30>(), std::declval<__type30>())) __type39;
      typedef decltype(std::declval<__type26>()(std::declval<__type39>())) __type40;
      typedef decltype(pythonic::operator_::add(std::declval<__type40>(), std::declval<__type40>())) __type43;
      typedef decltype(pythonic::operator_::sub(std::declval<__type35>(), std::declval<__type43>())) __type44;
      typedef decltype(std::declval<__type25>()(std::declval<__type44>())) __type45;
      typedef decltype(pythonic::operator_::mul(std::declval<__type24>(), std::declval<__type45>())) __type46;
      typedef typename __combined<__type6,__type46>::type __type47;
      typedef __type47 __type48;
      typedef typename pythonic::returnable<__type48>::type __type49;
      typedef __type2 __ptype4;
      typedef __type2 __ptype5;
      typedef __type49 result_type;
    }  
    ;
    template <typename argument_type0 >
    inline
    typename type<argument_type0>::result_type operator()(argument_type0&& A) const
    ;
  }  ;
  struct _discordant_pairs
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
      typedef __type0 __type1;
      typedef __type1 __type2;
      typedef long __type3;
      typedef typename pythonic::assignable<__type3>::type __type4;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type6;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type1>())) __type8;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type8>::type>::type __type9;
      typedef typename pythonic::lazy<__type9>::type __type10;
      typedef __type10 __type11;
      typedef decltype(std::declval<__type6>()(std::declval<__type11>())) __type12;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type12>::type::iterator>::value_type>::type __type13;
      typedef __type13 __type14;
      typedef typename std::tuple_element<1,typename std::remove_reference<__type8>::type>::type __type15;
      typedef typename pythonic::lazy<__type15>::type __type16;
      typedef __type16 __type17;
      typedef decltype(std::declval<__type6>()(std::declval<__type17>())) __type18;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type18>::type::iterator>::value_type>::type __type19;
      typedef __type19 __type20;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type14>(), std::declval<__type20>())) __type21;
      typedef decltype(std::declval<__type1>()[std::declval<__type21>()]) __type22;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::sum{})>::type>::type __type23;
      typedef typename pythonic::assignable<__type1>::type __type24;
      typedef __type24 __type25;
      typedef pythonic::types::contiguous_slice __type26;
      typedef decltype(std::declval<__type25>()(std::declval<__type26>(), std::declval<__type26>())) __type27;
      typedef decltype(std::declval<__type23>()(std::declval<__type27>())) __type28;
      typedef decltype(pythonic::operator_::add(std::declval<__type28>(), std::declval<__type28>())) __type31;
      typedef decltype(pythonic::operator_::mul(std::declval<__type22>(), std::declval<__type31>())) __type32;
      typedef typename __combined<__type4,__type32>::type __type33;
      typedef __type33 __type34;
      typedef typename pythonic::returnable<__type34>::type __type35;
      typedef __type2 __ptype6;
      typedef __type35 result_type;
    }  
    ;
    template <typename argument_type0 >
    inline
    typename type<argument_type0>::result_type operator()(argument_type0&& A) const
    ;
  }  ;
  struct _concordant_pairs
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
      typedef __type0 __type1;
      typedef __type1 __type2;
      typedef long __type3;
      typedef typename pythonic::assignable<__type3>::type __type4;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type6;
      typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type1>())) __type8;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type8>::type>::type __type9;
      typedef typename pythonic::lazy<__type9>::type __type10;
      typedef __type10 __type11;
      typedef decltype(std::declval<__type6>()(std::declval<__type11>())) __type12;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type12>::type::iterator>::value_type>::type __type13;
      typedef __type13 __type14;
      typedef typename std::tuple_element<1,typename std::remove_reference<__type8>::type>::type __type15;
      typedef typename pythonic::lazy<__type15>::type __type16;
      typedef __type16 __type17;
      typedef decltype(std::declval<__type6>()(std::declval<__type17>())) __type18;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type18>::type::iterator>::value_type>::type __type19;
      typedef __type19 __type20;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type14>(), std::declval<__type20>())) __type21;
      typedef decltype(std::declval<__type1>()[std::declval<__type21>()]) __type22;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::sum{})>::type>::type __type23;
      typedef typename pythonic::assignable<__type1>::type __type24;
      typedef __type24 __type25;
      typedef pythonic::types::contiguous_slice __type26;
      typedef decltype(std::declval<__type25>()(std::declval<__type26>(), std::declval<__type26>())) __type27;
      typedef decltype(std::declval<__type23>()(std::declval<__type27>())) __type28;
      typedef decltype(pythonic::operator_::add(std::declval<__type28>(), std::declval<__type28>())) __type31;
      typedef decltype(pythonic::operator_::mul(std::declval<__type22>(), std::declval<__type31>())) __type32;
      typedef typename __combined<__type4,__type32>::type __type33;
      typedef __type33 __type34;
      typedef typename pythonic::returnable<__type34>::type __type35;
      typedef __type2 __ptype7;
      typedef __type35 result_type;
    }  
    ;
    template <typename argument_type0 >
    inline
    typename type<argument_type0>::result_type operator()(argument_type0&& A) const
    ;
  }  ;
  struct _Dij
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::sum{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type1;
      typedef __type1 __type2;
      typedef pythonic::types::contiguous_slice __type3;
      typedef decltype(std::declval<__type2>()(std::declval<__type3>(), std::declval<__type3>())) __type4;
      typedef decltype(std::declval<__type0>()(std::declval<__type4>())) __type5;
      typedef decltype(pythonic::operator_::add(std::declval<__type5>(), std::declval<__type5>())) __type9;
      typedef typename pythonic::returnable<__type9>::type __type10;
      typedef __type10 result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    inline
    typename type<argument_type0, argument_type1, argument_type2>::result_type operator()(argument_type0&& A, argument_type1&& i, argument_type2&& j) const
    ;
  }  ;
  struct _Aij
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::sum{})>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type1;
      typedef __type1 __type2;
      typedef pythonic::types::contiguous_slice __type3;
      typedef decltype(std::declval<__type2>()(std::declval<__type3>(), std::declval<__type3>())) __type4;
      typedef decltype(std::declval<__type0>()(std::declval<__type4>())) __type5;
      typedef decltype(pythonic::operator_::add(std::declval<__type5>(), std::declval<__type5>())) __type9;
      typedef typename pythonic::returnable<__type9>::type __type10;
      typedef __type10 result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    inline
    typename type<argument_type0, argument_type1, argument_type2>::result_type operator()(argument_type0&& A, argument_type1&& i, argument_type2&& j) const
    ;
  }  ;
  struct _count_paths_outside_method
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type0;
      typedef __type0 __type1;
      typedef __type1 __type2;
      typedef typename pythonic::assignable<__type1>::type __type3;
      typedef __type3 __type4;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type5;
      typedef __type5 __type6;
      typedef typename pythonic::assignable<__type6>::type __type7;
      typedef __type7 __type8;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type4>(), std::declval<__type8>())) __type9;
      typedef typename pythonic::assignable<__type9>::type __type10;
      typedef __type10 __type11;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type11>::type>::type __type12;
      typedef __type12 __type13;
      typedef __type6 __type14;
      typedef typename std::tuple_element<1,typename std::remove_reference<__type11>::type>::type __type15;
      typedef __type15 __type16;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::scipy::special::functor::binom{})>::type>::type __type17;
      typedef typename pythonic::assignable<__type12>::type __type18;
      typedef typename __combined<__type7,__type18>::type __type19;
      typedef __type19 __type20;
      typedef typename pythonic::assignable<__type15>::type __type21;
      typedef typename __combined<__type3,__type21>::type __type22;
      typedef __type22 __type23;
      typedef decltype(pythonic::operator_::add(std::declval<__type20>(), std::declval<__type23>())) __type24;
      typedef decltype(std::declval<__type17>()(std::declval<__type24>(), std::declval<__type23>())) __type26;
      typedef long __type27;
      typedef typename pythonic::assignable<__type27>::type __type28;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::zeros{})>::type>::type __type29;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type32;
      typedef __type32 __type33;
      typedef decltype(pythonic::operator_::functor::floordiv()(std::declval<__type20>(), std::declval<__type33>())) __type34;
      typedef typename pythonic::assignable<__type34>::type __type35;
      typedef __type35 __type36;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type37;
      typedef __type37 __type38;
      typedef decltype(pythonic::operator_::sub(std::declval<__type36>(), std::declval<__type38>())) __type39;
      typedef decltype(pythonic::operator_::functor::floordiv()(std::declval<__type39>(), std::declval<__type36>())) __type41;
      typedef decltype(pythonic::operator_::add(std::declval<__type23>(), std::declval<__type41>())) __type42;
      typedef typename pythonic::assignable<__type42>::type __type43;
      typedef __type43 __type44;
      typedef decltype(std::declval<__type29>()(std::declval<__type44>())) __type45;
      typedef typename pythonic::assignable<__type45>::type __type46;
      typedef indexable<__type27> __type47;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type48;
      typedef decltype(std::declval<__type48>()(std::declval<__type27>(), std::declval<__type44>())) __type50;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type50>::type::iterator>::value_type>::type __type51;
      typedef __type51 __type52;
      typedef indexable<__type52> __type53;
      typedef std::integral_constant<long,0> __type54;
      typedef indexable_container<__type54, typename std::remove_reference<__type27>::type> __type55;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::pythran::functor::static_list{})>::type>::type __type56;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::map{})>::type>::type __type57;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::functools::functor::partial{})>::type>::type __type58;
      typedef _count_paths_outside_method_lambda0 __type59;
      typedef decltype(pythonic::operator_::functor::floordiv()(std::declval<__type23>(), std::declval<__type33>())) __type64;
      typedef decltype(std::declval<__type58>()(std::declval<__type59>(), std::declval<__type38>(), std::declval<__type36>(), std::declval<__type64>())) __type65;
      typedef decltype(std::declval<__type48>()(std::declval<__type44>())) __type67;
      typedef decltype(std::declval<__type57>()(std::declval<__type65>(), std::declval<__type67>())) __type68;
      typedef decltype(std::declval<__type56>()(std::declval<__type68>())) __type69;
      typedef typename pythonic::assignable<__type69>::type __type70;
      typedef __type70 __type71;
      typedef decltype(std::declval<__type71>()[std::declval<__type52>()]) __type73;
      typedef decltype(pythonic::operator_::add(std::declval<__type73>(), std::declval<__type52>())) __type75;
      typedef decltype(std::declval<__type17>()(std::declval<__type75>(), std::declval<__type52>())) __type77;
      typedef typename pythonic::assignable<__type77>::type __type78;
      typedef decltype(std::declval<__type48>()(std::declval<__type52>())) __type84;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type84>::type::iterator>::value_type>::type __type85;
      typedef __type85 __type86;
      typedef decltype(std::declval<__type71>()[std::declval<__type86>()]) __type87;
      typedef decltype(pythonic::operator_::sub(std::declval<__type73>(), std::declval<__type87>())) __type88;
      typedef decltype(pythonic::operator_::add(std::declval<__type88>(), std::declval<__type52>())) __type90;
      typedef decltype(pythonic::operator_::sub(std::declval<__type90>(), std::declval<__type86>())) __type92;
      typedef decltype(pythonic::operator_::sub(std::declval<__type52>(), std::declval<__type86>())) __type95;
      typedef decltype(std::declval<__type17>()(std::declval<__type92>(), std::declval<__type95>())) __type96;
      typedef typename __combined<__type46,__type47,__type55>::type __type97;
      typedef __type97 __type98;
      typedef decltype(std::declval<__type98>()[std::declval<__type86>()]) __type100;
      typedef decltype(pythonic::operator_::mul(std::declval<__type96>(), std::declval<__type100>())) __type101;
      typedef typename __combined<__type78,__type101>::type __type102;
      typedef __type102 __type103;
      typedef container<typename std::remove_reference<__type103>::type> __type104;
      typedef typename __combined<__type46,__type47,__type53,__type55,__type104>::type __type105;
      typedef __type105 __type106;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type67>::type::iterator>::value_type>::type __type109;
      typedef __type109 __type110;
      typedef decltype(std::declval<__type106>()[std::declval<__type110>()]) __type111;
      typedef decltype(std::declval<__type71>()[std::declval<__type110>()]) __type115;
      typedef decltype(pythonic::operator_::sub(std::declval<__type20>(), std::declval<__type115>())) __type116;
      typedef decltype(pythonic::operator_::sub(std::declval<__type23>(), std::declval<__type110>())) __type119;
      typedef decltype(pythonic::operator_::add(std::declval<__type116>(), std::declval<__type119>())) __type120;
      typedef decltype(std::declval<__type17>()(std::declval<__type120>(), std::declval<__type119>())) __type124;
      typedef decltype(pythonic::operator_::mul(std::declval<__type111>(), std::declval<__type124>())) __type125;
      typedef typename __combined<__type28,__type125>::type __type126;
      typedef __type126 __type127;
      typedef typename __combined<__type26,__type127>::type __type128;
      typedef typename pythonic::returnable<__type128>::type __type129;
      typedef __type2 __ptype8;
      typedef __type13 __ptype10;
      typedef __type14 __ptype9;
      typedef __type16 __ptype11;
      typedef __type129 result_type;
    }  
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    inline
    typename type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type operator()(argument_type0&& m, argument_type1&& n, argument_type2&& g, argument_type3&& h) const
    ;
  }  ;
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
  inline
  typename _count_paths_outside_method_lambda0::type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type _count_paths_outside_method_lambda0::operator()(argument_type0&& h, argument_type1&& mg, argument_type2&& ng, argument_type3&& j) const
  {
    return pythonic::operator_::functor::floordiv()(pythonic::operator_::sub(pythonic::operator_::add(pythonic::operator_::add(h, pythonic::operator_::mul(mg, j)), ng), 1L), ng);
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
  inline
  typename siegelslopes::type<argument_type0, argument_type1, argument_type2>::result_type siegelslopes::operator()(argument_type0&& y, argument_type1&& x, argument_type2&& method) const
  {
    typedef pythonic::types::empty_list __type0;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type0>(), std::declval<__type0>())) __type1;
    typedef typename pythonic::assignable<__type1>::type __type2;
    typedef __type2 __type3;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type3>::type>::type __type4;
    typedef typename pythonic::assignable<__type4>::type __type5;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::median{})>::type>::type __type6;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::expand_dims{})>::type>::type __type7;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type8;
    typedef __type8 __type9;
    typedef long __type10;
    typedef decltype(std::declval<__type7>()(std::declval<__type9>(), std::declval<__type10>())) __type11;
    typedef decltype(pythonic::operator_::sub(std::declval<__type11>(), std::declval<__type9>())) __type13;
    typedef typename pythonic::assignable<__type13>::type __type14;
    typedef __type14 __type15;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type16;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::len{})>::type>::type __type17;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type18;
    typedef __type18 __type19;
    typedef decltype(std::declval<__type17>()(std::declval<__type19>())) __type20;
    typedef decltype(std::declval<__type16>()(std::declval<__type20>())) __type21;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type21>::type::iterator>::value_type>::type __type22;
    typedef __type22 __type23;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::nonzero{})>::type>::type __type24;
    typedef decltype(std::declval<__type7>()(std::declval<__type19>(), std::declval<__type10>())) __type26;
    typedef decltype(pythonic::operator_::sub(std::declval<__type26>(), std::declval<__type19>())) __type28;
    typedef typename pythonic::assignable<__type28>::type __type29;
    typedef __type29 __type30;
    typedef pythonic::types::contiguous_slice __type32;
    typedef decltype(std::declval<__type30>()(std::declval<__type23>(), std::declval<__type32>())) __type33;
    typedef decltype(std::declval<__type24>()(std::declval<__type33>())) __type34;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type34>::type>::type __type35;
    typedef typename pythonic::assignable<__type35>::type __type36;
    typedef __type36 __type37;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type23>(), std::declval<__type37>())) __type38;
    typedef decltype(std::declval<__type15>()[std::declval<__type38>()]) __type39;
    typedef decltype(std::declval<__type30>()[std::declval<__type38>()]) __type44;
    typedef decltype(pythonic::operator_::div(std::declval<__type39>(), std::declval<__type44>())) __type45;
    typedef decltype(std::declval<__type6>()(std::declval<__type45>())) __type46;
    typedef pythonic::types::list<typename std::remove_reference<__type46>::type> __type47;
    typedef __type5 __type48;
    typedef typename __combined<__type5,__type47,__type48>::type __type49;
    typedef typename pythonic::assignable<__type49>::type __type50;
    typedef typename std::tuple_element<1,typename std::remove_reference<__type3>::type>::type __type51;
    typedef typename pythonic::assignable<__type51>::type __type52;
    typedef decltype(std::declval<__type19>()[std::declval<__type23>()]) __type56;
    typedef decltype(pythonic::operator_::mul(std::declval<__type9>(), std::declval<__type56>())) __type57;
    typedef decltype(std::declval<__type9>()[std::declval<__type23>()]) __type60;
    typedef decltype(pythonic::operator_::mul(std::declval<__type60>(), std::declval<__type19>())) __type62;
    typedef decltype(pythonic::operator_::sub(std::declval<__type57>(), std::declval<__type62>())) __type63;
    typedef decltype(std::declval<__type63>()[std::declval<__type37>()]) __type65;
    typedef decltype(pythonic::operator_::div(std::declval<__type65>(), std::declval<__type44>())) __type71;
    typedef decltype(std::declval<__type6>()(std::declval<__type71>())) __type72;
    typedef pythonic::types::list<typename std::remove_reference<__type72>::type> __type73;
    typedef __type52 __type74;
    typedef typename __combined<__type52,__type73,__type74>::type __type75;
    typedef typename pythonic::assignable<__type75>::type __type76;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::asarray{})>::type>::type __type77;
    typedef __type75 __type78;
    typedef decltype(std::declval<__type77>()(std::declval<__type78>())) __type79;
    typedef decltype(std::declval<__type6>()(std::declval<__type79>())) __type80;
    typedef typename pythonic::lazy<__type80>::type __type81;
    typedef __type49 __type83;
    typedef decltype(std::declval<__type77>()(std::declval<__type83>())) __type84;
    typedef decltype(std::declval<__type6>()(std::declval<__type84>())) __type85;
    typedef typename pythonic::assignable<__type85>::type __type86;
    typedef __type86 __type87;
    typedef decltype(pythonic::operator_::mul(std::declval<__type87>(), std::declval<__type19>())) __type89;
    typedef decltype(pythonic::operator_::sub(std::declval<__type9>(), std::declval<__type89>())) __type90;
    typedef decltype(std::declval<__type6>()(std::declval<__type90>())) __type91;
    typedef typename pythonic::lazy<__type91>::type __type92;
    typedef typename __combined<__type81,__type92>::type __type93;
    typedef typename pythonic::lazy<__type93>::type __type94;
    __type94 medinter;
    typename pythonic::assignable_noescape<decltype(pythonic::operator_::sub(pythonic::numpy::functor::expand_dims{}(x, 1L), x))>::type deltax = pythonic::operator_::sub(pythonic::numpy::functor::expand_dims{}(x, 1L), x);
    typename pythonic::assignable_noescape<decltype(pythonic::operator_::sub(pythonic::numpy::functor::expand_dims{}(y, 1L), y))>::type deltay = pythonic::operator_::sub(pythonic::numpy::functor::expand_dims{}(y, 1L), y);
    typename pythonic::assignable_noescape<decltype(pythonic::types::make_tuple(__type0(pythonic::types::empty_list()), __type0(pythonic::types::empty_list())))>::type __tuple0 = pythonic::types::make_tuple(__type0(pythonic::types::empty_list()), __type0(pythonic::types::empty_list()));
    __type50 slopes = std::get<0>(__tuple0);
    __type76 intercepts = std::get<1>(__tuple0);
    {
      long  __target2668469385248 = pythonic::builtins::functor::len{}(x);
      for (long  j=0L; j < __target2668469385248; j += 1L)
      {
        typename pythonic::assignable_noescape<decltype(std::get<0>(pythonic::numpy::functor::nonzero{}(deltax(j,pythonic::types::contiguous_slice(pythonic::builtins::None,pythonic::builtins::None)))))>::type id_nonzero = std::get<0>(pythonic::numpy::functor::nonzero{}(deltax(j,pythonic::types::contiguous_slice(pythonic::builtins::None,pythonic::builtins::None))));
        pythonic::builtins::list::functor::append{}(slopes, pythonic::numpy::functor::median{}(pythonic::operator_::div(deltay[pythonic::types::make_tuple(j, id_nonzero)], deltax[pythonic::types::make_tuple(j, id_nonzero)])));
        if (pythonic::operator_::eq(method, pythonic::types::str("separate")))
        {
          pythonic::builtins::list::functor::append{}(intercepts, pythonic::numpy::functor::median{}(pythonic::operator_::div(pythonic::operator_::sub(pythonic::operator_::mul(y, x[j]), pythonic::operator_::mul(y[j], x))[id_nonzero], deltax[pythonic::types::make_tuple(j, id_nonzero)])));
        }
      }
    }
    typename pythonic::assignable_noescape<decltype(pythonic::numpy::functor::median{}(pythonic::numpy::functor::asarray{}(slopes)))>::type medslope = pythonic::numpy::functor::median{}(pythonic::numpy::functor::asarray{}(slopes));
    if (pythonic::operator_::eq(method, pythonic::types::str("separate")))
    {
      medinter = pythonic::numpy::functor::median{}(pythonic::numpy::functor::asarray{}(intercepts));
    }
    else
    {
      medinter = pythonic::numpy::functor::median{}(pythonic::operator_::sub(y, pythonic::operator_::mul(medslope, x)));
    }
    return pythonic::types::make_tuple(medslope, medinter);
  }
  template <typename argument_type0 , typename argument_type1 >
  inline
  typename _compute_prob_outside_square::type<argument_type0, argument_type1>::result_type _compute_prob_outside_square::operator()(argument_type0&& n, argument_type1&& h) const
  {
    typedef double __type0;
    typedef typename pythonic::lazy<__type0>::type __type1;
    typedef typename pythonic::assignable<__type0>::type __type2;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type3;
    typedef __type3 __type4;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::int_{})>::type>::type __type5;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::floor{})>::type>::type __type6;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type8;
    typedef __type8 __type9;
    typedef decltype(pythonic::operator_::div(std::declval<__type4>(), std::declval<__type9>())) __type10;
    typedef decltype(std::declval<__type6>()(std::declval<__type10>())) __type11;
    typedef decltype(std::declval<__type5>()(std::declval<__type11>())) __type12;
    typedef typename pythonic::assignable<__type12>::type __type13;
    typedef long __type14;
    typedef typename __combined<__type13,__type14>::type __type15;
    typedef __type15 __type16;
    typedef decltype(pythonic::operator_::mul(std::declval<__type16>(), std::declval<__type9>())) __type18;
    typedef decltype(pythonic::operator_::sub(std::declval<__type4>(), std::declval<__type18>())) __type19;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type20;
    typedef decltype(std::declval<__type20>()(std::declval<__type9>())) __type22;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type22>::type::iterator>::value_type>::type __type23;
    typedef __type23 __type24;
    typedef decltype(pythonic::operator_::sub(std::declval<__type19>(), std::declval<__type24>())) __type25;
    typedef __type2 __type26;
    typedef decltype(pythonic::operator_::mul(std::declval<__type25>(), std::declval<__type26>())) __type27;
    typedef decltype(pythonic::operator_::add(std::declval<__type4>(), std::declval<__type18>())) __type32;
    typedef decltype(pythonic::operator_::add(std::declval<__type32>(), std::declval<__type24>())) __type34;
    typedef decltype(pythonic::operator_::add(std::declval<__type34>(), std::declval<__type14>())) __type35;
    typedef decltype(pythonic::operator_::div(std::declval<__type27>(), std::declval<__type35>())) __type36;
    typedef typename pythonic::assignable<__type36>::type __type37;
    typedef typename __combined<__type2,__type37>::type __type38;
    typedef __type38 __type39;
    typedef __type1 __type40;
    typedef decltype(pythonic::operator_::sub(std::declval<__type0>(), std::declval<__type40>())) __type41;
    typedef decltype(pythonic::operator_::mul(std::declval<__type39>(), std::declval<__type41>())) __type42;
    typedef typename pythonic::lazy<__type42>::type __type43;
    typedef typename __combined<__type1,__type43>::type __type44;
    typedef typename pythonic::lazy<__type44>::type __type45;
    typedef typename pythonic::assignable<__type15>::type __type46;
    typedef typename pythonic::assignable<__type38>::type __type47;
    __type45 P = 0.0;
    __type46 k = pythonic::builtins::functor::int_{}(pythonic::numpy::functor::floor{}(pythonic::operator_::div(n, h)));
    while (pythonic::operator_::ge(k, 0L))
    {
      __type47 p1 = 1.0;
      {
        long  __target2668497867760 = h;
        for (long  j=0L; j < __target2668497867760; j += 1L)
        {
          p1 = pythonic::operator_::div(pythonic::operator_::mul(pythonic::operator_::sub(pythonic::operator_::sub(n, pythonic::operator_::mul(k, h)), j), p1), pythonic::operator_::add(pythonic::operator_::add(pythonic::operator_::add(n, pythonic::operator_::mul(k, h)), j), 1L));
        }
      }
      P = pythonic::operator_::mul(p1, pythonic::operator_::sub(1.0, P));
      k -= 1L;
    }
    return pythonic::operator_::mul(2L, P);
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
  inline
  typename _compute_outer_prob_inside_method::type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type _compute_outer_prob_inside_method::operator()(argument_type0&& m, argument_type1&& n, argument_type2&& g, argument_type3&& h) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type0;
    typedef __type0 __type1;
    typedef typename pythonic::assignable<__type1>::type __type2;
    typedef __type2 __type3;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type4;
    typedef __type4 __type5;
    typedef typename pythonic::assignable<__type5>::type __type6;
    typedef __type6 __type7;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type3>(), std::declval<__type7>())) __type8;
    typedef typename pythonic::assignable<__type8>::type __type9;
    typedef __type9 __type10;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type10>::type>::type __type11;
    typedef typename pythonic::assignable<__type11>::type __type12;
    typedef typename __combined<__type6,__type12>::type __type13;
    typedef typename std::tuple_element<1,typename std::remove_reference<__type10>::type>::type __type14;
    typedef typename pythonic::assignable<__type14>::type __type15;
    typedef typename __combined<__type2,__type15>::type __type16;
    typedef typename pythonic::assignable<__type16>::type __type17;
    typedef typename pythonic::assignable<__type13>::type __type18;
    typedef long __type19;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::min{})>::type>::type __type20;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::int_{})>::type>::type __type21;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::ceil{})>::type>::type __type22;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type23;
    typedef __type23 __type24;
    typedef __type13 __type25;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type26;
    typedef __type26 __type27;
    typedef decltype(pythonic::operator_::functor::floordiv()(std::declval<__type25>(), std::declval<__type27>())) __type28;
    typedef typename pythonic::assignable<__type28>::type __type29;
    typedef __type29 __type30;
    typedef decltype(pythonic::operator_::div(std::declval<__type24>(), std::declval<__type30>())) __type31;
    typedef decltype(std::declval<__type22>()(std::declval<__type31>())) __type32;
    typedef decltype(std::declval<__type21>()(std::declval<__type32>())) __type33;
    typedef __type16 __type34;
    typedef decltype(pythonic::operator_::add(std::declval<__type34>(), std::declval<__type19>())) __type35;
    typedef decltype(std::declval<__type20>()(std::declval<__type33>(), std::declval<__type35>())) __type36;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type19>(), std::declval<__type36>())) __type37;
    typedef typename pythonic::assignable<__type37>::type __type38;
    typedef __type38 __type39;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type39>::type>::type __type40;
    typedef typename pythonic::assignable<__type40>::type __type41;
    typedef __type41 __type42;
    typedef typename std::tuple_element<1,typename std::remove_reference<__type39>::type>::type __type43;
    typedef typename pythonic::assignable<__type43>::type __type44;
    typedef decltype(pythonic::operator_::functor::floordiv()(std::declval<__type34>(), std::declval<__type27>())) __type47;
    typedef typename pythonic::assignable<__type47>::type __type48;
    typedef __type48 __type49;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type50;
    typedef decltype(pythonic::operator_::add(std::declval<__type25>(), std::declval<__type19>())) __type52;
    typedef decltype(std::declval<__type50>()(std::declval<__type19>(), std::declval<__type52>())) __type53;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type53>::type::iterator>::value_type>::type __type54;
    typedef __type54 __type55;
    typedef decltype(pythonic::operator_::mul(std::declval<__type49>(), std::declval<__type55>())) __type56;
    typedef decltype(pythonic::operator_::add(std::declval<__type56>(), std::declval<__type24>())) __type58;
    typedef decltype(pythonic::operator_::div(std::declval<__type58>(), std::declval<__type30>())) __type60;
    typedef decltype(std::declval<__type22>()(std::declval<__type60>())) __type61;
    typedef decltype(std::declval<__type21>()(std::declval<__type61>())) __type62;
    typedef decltype(std::declval<__type20>()(std::declval<__type62>(), std::declval<__type35>())) __type65;
    typedef typename pythonic::assignable<__type65>::type __type66;
    typedef typename __combined<__type44,__type66>::type __type67;
    typedef __type67 __type68;
    typedef decltype(pythonic::operator_::sub(std::declval<__type68>(), std::declval<__type42>())) __type70;
    typedef typename pythonic::assignable<__type70>::type __type71;
    typedef __type71 __type72;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type42>(), std::declval<__type72>())) __type73;
    typedef typename pythonic::assignable<__type73>::type __type74;
    typedef __type74 __type75;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type75>::type>::type __type76;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::max{})>::type>::type __type77;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::floor{})>::type>::type __type78;
    typedef decltype(pythonic::operator_::sub(std::declval<__type56>(), std::declval<__type24>())) __type83;
    typedef decltype(pythonic::operator_::div(std::declval<__type83>(), std::declval<__type30>())) __type85;
    typedef decltype(std::declval<__type78>()(std::declval<__type85>())) __type86;
    typedef decltype(std::declval<__type21>()(std::declval<__type86>())) __type87;
    typedef decltype(pythonic::operator_::add(std::declval<__type87>(), std::declval<__type19>())) __type88;
    typedef typename __combined<__type88,__type19>::type __type89;
    typedef decltype(std::declval<__type77>()(std::declval<__type89>(), std::declval<__type19>())) __type90;
    typedef decltype(std::declval<__type20>()(std::declval<__type90>(), std::declval<__type34>())) __type92;
    typedef typename __combined<__type76,__type92>::type __type93;
    typedef typename pythonic::assignable<__type93>::type __type94;
    typedef typename __combined<__type41,__type94>::type __type95;
    typedef __type95 __type97;
    typedef decltype(pythonic::operator_::sub(std::declval<__type68>(), std::declval<__type97>())) __type98;
    typedef typename pythonic::assignable<__type98>::type __type99;
    typedef typename __combined<__type71,__type99>::type __type100;
    typedef typename pythonic::assignable<__type95>::type __type101;
    typedef typename pythonic::assignable<__type67>::type __type102;
    typedef typename pythonic::assignable<__type100>::type __type103;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::ones{})>::type>::type __type104;
    typedef decltype(pythonic::operator_::mul(std::declval<__type19>(), std::declval<__type68>())) __type106;
    typedef decltype(pythonic::operator_::add(std::declval<__type106>(), std::declval<__type19>())) __type107;
    typedef decltype(std::declval<__type20>()(std::declval<__type107>(), std::declval<__type35>())) __type110;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::float64{})>::type>::type __type111;
    typedef decltype(std::declval<__type104>()(std::declval<__type110>(), std::declval<__type111>())) __type112;
    typedef typename pythonic::assignable<__type112>::type __type113;
    typedef decltype(std::declval<__type50>()(std::declval<__type98>())) __type117;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type117>::type::iterator>::value_type>::type __type118;
    typedef __type118 __type119;
    typedef indexable<__type119> __type120;
    typedef double __type121;
    typedef typename pythonic::assignable<__type121>::type __type122;
    typedef typename __combined<__type113,__type121>::type __type123;
    typedef __type123 __type124;
    typedef decltype(pythonic::operator_::add(std::declval<__type119>(), std::declval<__type97>())) __type127;
    typedef typename pythonic::assignable<__type76>::type __type128;
    typedef __type128 __type129;
    typedef decltype(pythonic::operator_::sub(std::declval<__type127>(), std::declval<__type129>())) __type130;
    typedef decltype(std::declval<__type124>()[std::declval<__type130>()]) __type131;
    typedef decltype(pythonic::operator_::mul(std::declval<__type131>(), std::declval<__type55>())) __type133;
    typedef __type122 __type134;
    typedef typename pythonic::assignable<__type127>::type __type138;
    typedef __type138 __type139;
    typedef decltype(pythonic::operator_::mul(std::declval<__type134>(), std::declval<__type139>())) __type140;
    typedef decltype(pythonic::operator_::add(std::declval<__type133>(), std::declval<__type140>())) __type141;
    typedef decltype(pythonic::operator_::add(std::declval<__type55>(), std::declval<__type139>())) __type144;
    typedef decltype(pythonic::operator_::div(std::declval<__type141>(), std::declval<__type144>())) __type145;
    typedef typename pythonic::assignable<__type145>::type __type146;
    typedef typename __combined<__type122,__type146>::type __type147;
    typedef __type147 __type148;
    typedef container<typename std::remove_reference<__type148>::type> __type149;
    typedef typename __combined<__type113,__type120,__type121,__type149,__type19>::type __type150;
    typedef typename pythonic::assignable<__type150>::type __type151;
    typedef typename pythonic::assignable<__type147>::type __type152;
    __type17 n_ = n;
    __type18 m_ = m;
    if (pythonic::operator_::lt(m_, n_))
    {
      typename pythonic::assignable_noescape<decltype(pythonic::types::make_tuple(n_, m_))>::type __tuple0 = pythonic::types::make_tuple(n_, m_);
      m_ = std::get<0>(__tuple0);
      n_ = std::get<1>(__tuple0);
    }
    typename pythonic::assignable_noescape<decltype(pythonic::operator_::functor::floordiv()(m_, g))>::type mg = pythonic::operator_::functor::floordiv()(m_, g);
    typename pythonic::assignable_noescape<decltype(pythonic::operator_::functor::floordiv()(n_, g))>::type ng = pythonic::operator_::functor::floordiv()(n_, g);
    typename pythonic::assignable_noescape<decltype(pythonic::types::make_tuple(0L, pythonic::builtins::functor::min{}(pythonic::builtins::functor::int_{}(pythonic::numpy::functor::ceil{}(pythonic::operator_::div(h, mg))), pythonic::operator_::add(n_, 1L))))>::type __tuple1 = pythonic::types::make_tuple(0L, pythonic::builtins::functor::min{}(pythonic::builtins::functor::int_{}(pythonic::numpy::functor::ceil{}(pythonic::operator_::div(h, mg))), pythonic::operator_::add(n_, 1L)));
    __type101 minj = std::get<0>(__tuple1);
    __type102 maxj = std::get<1>(__tuple1);
    __type103 curlen = pythonic::operator_::sub(maxj, minj);
    __type151 A = pythonic::numpy::functor::ones{}(pythonic::builtins::functor::min{}(pythonic::operator_::add(pythonic::operator_::mul(2L, maxj), 2L), pythonic::operator_::add(n_, 1L)), pythonic::numpy::functor::float64{});
    A[pythonic::types::contiguous_slice(minj,maxj)] = 0.0;
    {
      long  __target2668490266256 = pythonic::operator_::add(m_, 1L);
      for (long  i=1L; i < __target2668490266256; i += 1L)
      {
        typename pythonic::assignable_noescape<decltype(pythonic::types::make_tuple(minj, curlen))>::type __tuple2 = pythonic::types::make_tuple(minj, curlen);
        typename pythonic::assignable_noescape<decltype(std::get<0>(__tuple2))>::type lastminj = std::get<0>(__tuple2);
        minj = pythonic::builtins::functor::min{}(pythonic::builtins::functor::max{}(pythonic::operator_::add(pythonic::builtins::functor::int_{}(pythonic::numpy::functor::floor{}(pythonic::operator_::div(pythonic::operator_::sub(pythonic::operator_::mul(ng, i), h), mg))), 1L), 0L), n_);
        maxj = pythonic::builtins::functor::min{}(pythonic::builtins::functor::int_{}(pythonic::numpy::functor::ceil{}(pythonic::operator_::div(pythonic::operator_::add(pythonic::operator_::mul(ng, i), h), mg))), pythonic::operator_::add(n_, 1L));
        {
          __type152 val;
          if (pythonic::operator_::le(maxj, minj))
          {
            return 1.0;
          }
          else
          {
            val = (((bool)pythonic::operator_::eq(minj, 0L)) ? typename __combined<decltype(0.0), decltype(1.0)>::type(0.0) : typename __combined<decltype(0.0), decltype(1.0)>::type(1.0));
            {
              long  __target2668490270720 = pythonic::operator_::sub(maxj, minj);
              for (long  jj=0L; jj < __target2668490270720; jj += 1L)
              {
                typename pythonic::assignable_noescape<decltype(pythonic::operator_::add(jj, minj))>::type j = pythonic::operator_::add(jj, minj);
                val = pythonic::operator_::div(pythonic::operator_::add(pythonic::operator_::mul(A[pythonic::operator_::sub(pythonic::operator_::add(jj, minj), lastminj)], i), pythonic::operator_::mul(val, j)), pythonic::operator_::add(i, j));
                A[jj] = val;
              }
            }
            curlen = pythonic::operator_::sub(maxj, minj);
            if (pythonic::operator_::gt(std::get<1>(__tuple2), curlen))
            {
              A[pythonic::types::contiguous_slice(pythonic::operator_::sub(maxj, minj),pythonic::operator_::add(pythonic::operator_::sub(maxj, minj), pythonic::operator_::sub(std::get<1>(__tuple2), curlen)))] = 1L;
            }
          }
        }
      }
    }
    return A[pythonic::operator_::sub(pythonic::operator_::sub(maxj, minj), 1L)];
  }
  template <typename argument_type0 >
  inline
  typename _a_ij_Aij_Dij2::type<argument_type0>::result_type _a_ij_Aij_Dij2::operator()(argument_type0&& A) const
  {
    typedef long __type0;
    typedef typename pythonic::assignable<__type0>::type __type1;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
    typedef __type2 __type3;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type4;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type3>())) __type6;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type __type7;
    typedef typename pythonic::lazy<__type7>::type __type8;
    typedef __type8 __type9;
    typedef decltype(std::declval<__type4>()(std::declval<__type9>())) __type10;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type10>::type::iterator>::value_type>::type __type11;
    typedef __type11 __type12;
    typedef typename std::tuple_element<1,typename std::remove_reference<__type6>::type>::type __type13;
    typedef typename pythonic::lazy<__type13>::type __type14;
    typedef __type14 __type15;
    typedef decltype(std::declval<__type4>()(std::declval<__type15>())) __type16;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type16>::type::iterator>::value_type>::type __type17;
    typedef __type17 __type18;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type12>(), std::declval<__type18>())) __type19;
    typedef decltype(std::declval<__type3>()[std::declval<__type19>()]) __type20;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::square{})>::type>::type __type21;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::sum{})>::type>::type __type22;
    typedef typename pythonic::assignable<__type3>::type __type24;
    typedef typename __combined<__type24,__type3>::type __type26;
    typedef __type26 __type27;
    typedef pythonic::types::contiguous_slice __type28;
    typedef decltype(std::declval<__type27>()(std::declval<__type28>(), std::declval<__type28>())) __type29;
    typedef decltype(std::declval<__type22>()(std::declval<__type29>())) __type30;
    typedef decltype(pythonic::operator_::add(std::declval<__type30>(), std::declval<__type30>())) __type33;
    typedef __type3 __type34;
    typedef typename pythonic::assignable<__type34>::type __type35;
    typedef __type35 __type36;
    typedef decltype(std::declval<__type36>()(std::declval<__type28>(), std::declval<__type28>())) __type37;
    typedef decltype(std::declval<__type22>()(std::declval<__type37>())) __type38;
    typedef decltype(pythonic::operator_::add(std::declval<__type38>(), std::declval<__type38>())) __type41;
    typedef decltype(pythonic::operator_::sub(std::declval<__type33>(), std::declval<__type41>())) __type42;
    typedef decltype(std::declval<__type21>()(std::declval<__type42>())) __type43;
    typedef decltype(pythonic::operator_::mul(std::declval<__type20>(), std::declval<__type43>())) __type44;
    typedef typename __combined<__type1,__type44>::type __type45;
    typedef typename pythonic::assignable<__type45>::type __type46;
    typedef typename pythonic::assignable<__type26>::type __type47;
    typedef typename pythonic::assignable<__type35>::type __type48;
    typename pythonic::lazy<decltype(std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A)))>::type m = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A));
    typename pythonic::lazy<decltype(std::get<1>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A)))>::type n = std::get<1>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A));
    __type46 count = 0L;
    {
      long  __target2668487684784 = m;
      for (long  i=0L; i < __target2668487684784; i += 1L)
      {
        {
          long  __target2668487682144 = n;
          for (long  j=0L; j < __target2668487682144; j += 1L)
          {
            __type47 __pythran_inline_AijA2 = A;
            typename pythonic::assignable_noescape<decltype(i)>::type __pythran_inline_Aiji2 = i;
            typename pythonic::assignable_noescape<decltype(j)>::type __pythran_inline_Aijj2 = j;
            __type48 __pythran_inline_DijA3 = A;
            typename pythonic::assignable_noescape<decltype(i)>::type __pythran_inline_Diji3 = i;
            typename pythonic::assignable_noescape<decltype(j)>::type __pythran_inline_Dijj3 = j;
            count += pythonic::operator_::mul(A.fast(pythonic::types::make_tuple(i, j)), pythonic::numpy::functor::square{}(pythonic::operator_::sub(pythonic::operator_::add(pythonic::numpy::functor::sum{}(__pythran_inline_AijA2(pythonic::types::contiguous_slice(pythonic::builtins::None,__pythran_inline_Aiji2),pythonic::types::contiguous_slice(pythonic::builtins::None,__pythran_inline_Aijj2))), pythonic::numpy::functor::sum{}(__pythran_inline_AijA2(pythonic::types::contiguous_slice(pythonic::operator_::add(__pythran_inline_Aiji2, 1L),pythonic::builtins::None),pythonic::types::contiguous_slice(pythonic::operator_::add(__pythran_inline_Aijj2, 1L),pythonic::builtins::None)))), pythonic::operator_::add(pythonic::numpy::functor::sum{}(__pythran_inline_DijA3(pythonic::types::contiguous_slice(pythonic::operator_::add(__pythran_inline_Diji3, 1L),pythonic::builtins::None),pythonic::types::contiguous_slice(pythonic::builtins::None,__pythran_inline_Dijj3))), pythonic::numpy::functor::sum{}(__pythran_inline_DijA3(pythonic::types::contiguous_slice(pythonic::builtins::None,__pythran_inline_Diji3),pythonic::types::contiguous_slice(pythonic::operator_::add(__pythran_inline_Dijj3, 1L),pythonic::builtins::None)))))));
          }
        }
      }
    }
    return count;
  }
  template <typename argument_type0 >
  inline
  typename _discordant_pairs::type<argument_type0>::result_type _discordant_pairs::operator()(argument_type0&& A) const
  {
    typedef long __type0;
    typedef typename pythonic::assignable<__type0>::type __type1;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
    typedef __type2 __type3;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type4;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type3>())) __type6;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type __type7;
    typedef typename pythonic::lazy<__type7>::type __type8;
    typedef __type8 __type9;
    typedef decltype(std::declval<__type4>()(std::declval<__type9>())) __type10;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type10>::type::iterator>::value_type>::type __type11;
    typedef __type11 __type12;
    typedef typename std::tuple_element<1,typename std::remove_reference<__type6>::type>::type __type13;
    typedef typename pythonic::lazy<__type13>::type __type14;
    typedef __type14 __type15;
    typedef decltype(std::declval<__type4>()(std::declval<__type15>())) __type16;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type16>::type::iterator>::value_type>::type __type17;
    typedef __type17 __type18;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type12>(), std::declval<__type18>())) __type19;
    typedef decltype(std::declval<__type3>()[std::declval<__type19>()]) __type20;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::sum{})>::type>::type __type21;
    typedef typename pythonic::assignable<__type3>::type __type23;
    typedef __type23 __type24;
    typedef pythonic::types::contiguous_slice __type25;
    typedef decltype(std::declval<__type24>()(std::declval<__type25>(), std::declval<__type25>())) __type26;
    typedef decltype(std::declval<__type21>()(std::declval<__type26>())) __type27;
    typedef decltype(pythonic::operator_::add(std::declval<__type27>(), std::declval<__type27>())) __type30;
    typedef decltype(pythonic::operator_::mul(std::declval<__type20>(), std::declval<__type30>())) __type31;
    typedef typename __combined<__type1,__type31>::type __type32;
    typedef typename pythonic::assignable<__type32>::type __type33;
    typename pythonic::lazy<decltype(std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A)))>::type m = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A));
    typename pythonic::lazy<decltype(std::get<1>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A)))>::type n = std::get<1>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A));
    __type33 count = 0L;
    {
      long  __target2668490395648 = m;
      for (long  i=0L; i < __target2668490395648; i += 1L)
      {
        {
          long  __target2668490397376 = n;
          for (long  j=0L; j < __target2668490397376; j += 1L)
          {
            typename pythonic::assignable_noescape<decltype(A)>::type __pythran_inline_DijA1 = A;
            typename pythonic::assignable_noescape<decltype(i)>::type __pythran_inline_Diji1 = i;
            typename pythonic::assignable_noescape<decltype(j)>::type __pythran_inline_Dijj1 = j;
            count += pythonic::operator_::mul(A.fast(pythonic::types::make_tuple(i, j)), pythonic::operator_::add(pythonic::numpy::functor::sum{}(__pythran_inline_DijA1(pythonic::types::contiguous_slice(pythonic::operator_::add(__pythran_inline_Diji1, 1L),pythonic::builtins::None),pythonic::types::contiguous_slice(pythonic::builtins::None,__pythran_inline_Dijj1))), pythonic::numpy::functor::sum{}(__pythran_inline_DijA1(pythonic::types::contiguous_slice(pythonic::builtins::None,__pythran_inline_Diji1),pythonic::types::contiguous_slice(pythonic::operator_::add(__pythran_inline_Dijj1, 1L),pythonic::builtins::None)))));
          }
        }
      }
    }
    return count;
  }
  template <typename argument_type0 >
  inline
  typename _concordant_pairs::type<argument_type0>::result_type _concordant_pairs::operator()(argument_type0&& A) const
  {
    typedef long __type0;
    typedef typename pythonic::assignable<__type0>::type __type1;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type2;
    typedef __type2 __type3;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type4;
    typedef decltype(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, std::declval<__type3>())) __type6;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type __type7;
    typedef typename pythonic::lazy<__type7>::type __type8;
    typedef __type8 __type9;
    typedef decltype(std::declval<__type4>()(std::declval<__type9>())) __type10;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type10>::type::iterator>::value_type>::type __type11;
    typedef __type11 __type12;
    typedef typename std::tuple_element<1,typename std::remove_reference<__type6>::type>::type __type13;
    typedef typename pythonic::lazy<__type13>::type __type14;
    typedef __type14 __type15;
    typedef decltype(std::declval<__type4>()(std::declval<__type15>())) __type16;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type16>::type::iterator>::value_type>::type __type17;
    typedef __type17 __type18;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type12>(), std::declval<__type18>())) __type19;
    typedef decltype(std::declval<__type3>()[std::declval<__type19>()]) __type20;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::sum{})>::type>::type __type21;
    typedef typename pythonic::assignable<__type3>::type __type23;
    typedef __type23 __type24;
    typedef pythonic::types::contiguous_slice __type25;
    typedef decltype(std::declval<__type24>()(std::declval<__type25>(), std::declval<__type25>())) __type26;
    typedef decltype(std::declval<__type21>()(std::declval<__type26>())) __type27;
    typedef decltype(pythonic::operator_::add(std::declval<__type27>(), std::declval<__type27>())) __type30;
    typedef decltype(pythonic::operator_::mul(std::declval<__type20>(), std::declval<__type30>())) __type31;
    typedef typename __combined<__type1,__type31>::type __type32;
    typedef typename pythonic::assignable<__type32>::type __type33;
    typename pythonic::lazy<decltype(std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A)))>::type m = std::get<0>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A));
    typename pythonic::lazy<decltype(std::get<1>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A)))>::type n = std::get<1>(pythonic::builtins::getattr(pythonic::types::attr::SHAPE{}, A));
    __type33 count = 0L;
    {
      long  __target2668490392624 = m;
      for (long  i=0L; i < __target2668490392624; i += 1L)
      {
        {
          long  __target2668490394208 = n;
          for (long  j=0L; j < __target2668490394208; j += 1L)
          {
            typename pythonic::assignable_noescape<decltype(A)>::type __pythran_inline_AijA0 = A;
            typename pythonic::assignable_noescape<decltype(i)>::type __pythran_inline_Aiji0 = i;
            typename pythonic::assignable_noescape<decltype(j)>::type __pythran_inline_Aijj0 = j;
            count += pythonic::operator_::mul(A.fast(pythonic::types::make_tuple(i, j)), pythonic::operator_::add(pythonic::numpy::functor::sum{}(__pythran_inline_AijA0(pythonic::types::contiguous_slice(pythonic::builtins::None,__pythran_inline_Aiji0),pythonic::types::contiguous_slice(pythonic::builtins::None,__pythran_inline_Aijj0))), pythonic::numpy::functor::sum{}(__pythran_inline_AijA0(pythonic::types::contiguous_slice(pythonic::operator_::add(__pythran_inline_Aiji0, 1L),pythonic::builtins::None),pythonic::types::contiguous_slice(pythonic::operator_::add(__pythran_inline_Aijj0, 1L),pythonic::builtins::None)))));
          }
        }
      }
    }
    return count;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
  inline
  typename _Dij::type<argument_type0, argument_type1, argument_type2>::result_type _Dij::operator()(argument_type0&& A, argument_type1&& i, argument_type2&& j) const
  {
    return pythonic::operator_::add(pythonic::numpy::functor::sum{}(A(pythonic::types::contiguous_slice(pythonic::operator_::add(i, 1L),pythonic::builtins::None),pythonic::types::contiguous_slice(pythonic::builtins::None,j))), pythonic::numpy::functor::sum{}(A(pythonic::types::contiguous_slice(pythonic::builtins::None,i),pythonic::types::contiguous_slice(pythonic::operator_::add(j, 1L),pythonic::builtins::None))));
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
  inline
  typename _Aij::type<argument_type0, argument_type1, argument_type2>::result_type _Aij::operator()(argument_type0&& A, argument_type1&& i, argument_type2&& j) const
  {
    return pythonic::operator_::add(pythonic::numpy::functor::sum{}(A(pythonic::types::contiguous_slice(pythonic::builtins::None,i),pythonic::types::contiguous_slice(pythonic::builtins::None,j))), pythonic::numpy::functor::sum{}(A(pythonic::types::contiguous_slice(pythonic::operator_::add(i, 1L),pythonic::builtins::None),pythonic::types::contiguous_slice(pythonic::operator_::add(j, 1L),pythonic::builtins::None))));
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
  inline
  typename _count_paths_outside_method::type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type _count_paths_outside_method::operator()(argument_type0&& m, argument_type1&& n, argument_type2&& g, argument_type3&& h) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type0;
    typedef __type0 __type1;
    typedef typename pythonic::assignable<__type1>::type __type2;
    typedef __type2 __type3;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type4;
    typedef __type4 __type5;
    typedef typename pythonic::assignable<__type5>::type __type6;
    typedef __type6 __type7;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type3>(), std::declval<__type7>())) __type8;
    typedef typename pythonic::assignable<__type8>::type __type9;
    typedef __type9 __type10;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type10>::type>::type __type11;
    typedef typename pythonic::assignable<__type11>::type __type12;
    typedef typename __combined<__type6,__type12>::type __type13;
    typedef typename std::tuple_element<1,typename std::remove_reference<__type10>::type>::type __type14;
    typedef typename pythonic::assignable<__type14>::type __type15;
    typedef typename __combined<__type2,__type15>::type __type16;
    typedef typename pythonic::assignable<__type16>::type __type17;
    typedef typename pythonic::assignable<__type13>::type __type18;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::zeros{})>::type>::type __type19;
    typedef __type16 __type20;
    typedef __type13 __type21;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type22;
    typedef __type22 __type23;
    typedef decltype(pythonic::operator_::functor::floordiv()(std::declval<__type21>(), std::declval<__type23>())) __type24;
    typedef typename pythonic::assignable<__type24>::type __type25;
    typedef __type25 __type26;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type27;
    typedef __type27 __type28;
    typedef decltype(pythonic::operator_::sub(std::declval<__type26>(), std::declval<__type28>())) __type29;
    typedef decltype(pythonic::operator_::functor::floordiv()(std::declval<__type29>(), std::declval<__type26>())) __type31;
    typedef decltype(pythonic::operator_::add(std::declval<__type20>(), std::declval<__type31>())) __type32;
    typedef typename pythonic::assignable<__type32>::type __type33;
    typedef __type33 __type34;
    typedef decltype(std::declval<__type19>()(std::declval<__type34>())) __type35;
    typedef typename pythonic::assignable<__type35>::type __type36;
    typedef long __type37;
    typedef indexable<__type37> __type38;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type39;
    typedef decltype(std::declval<__type39>()(std::declval<__type37>(), std::declval<__type34>())) __type41;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type41>::type::iterator>::value_type>::type __type42;
    typedef __type42 __type43;
    typedef indexable<__type43> __type44;
    typedef std::integral_constant<long,0> __type45;
    typedef indexable_container<__type45, typename std::remove_reference<__type37>::type> __type46;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::scipy::special::functor::binom{})>::type>::type __type47;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::pythran::functor::static_list{})>::type>::type __type48;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::map{})>::type>::type __type49;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::functools::functor::partial{})>::type>::type __type50;
    typedef _count_paths_outside_method_lambda0 __type51;
    typedef decltype(pythonic::operator_::functor::floordiv()(std::declval<__type20>(), std::declval<__type23>())) __type56;
    typedef decltype(std::declval<__type50>()(std::declval<__type51>(), std::declval<__type28>(), std::declval<__type26>(), std::declval<__type56>())) __type57;
    typedef decltype(std::declval<__type39>()(std::declval<__type34>())) __type59;
    typedef decltype(std::declval<__type49>()(std::declval<__type57>(), std::declval<__type59>())) __type60;
    typedef decltype(std::declval<__type48>()(std::declval<__type60>())) __type61;
    typedef typename pythonic::assignable<__type61>::type __type62;
    typedef __type62 __type63;
    typedef decltype(std::declval<__type63>()[std::declval<__type43>()]) __type65;
    typedef decltype(pythonic::operator_::add(std::declval<__type65>(), std::declval<__type43>())) __type67;
    typedef decltype(std::declval<__type47>()(std::declval<__type67>(), std::declval<__type43>())) __type69;
    typedef typename pythonic::assignable<__type69>::type __type70;
    typedef decltype(std::declval<__type39>()(std::declval<__type43>())) __type76;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type76>::type::iterator>::value_type>::type __type77;
    typedef __type77 __type78;
    typedef decltype(std::declval<__type63>()[std::declval<__type78>()]) __type79;
    typedef decltype(pythonic::operator_::sub(std::declval<__type65>(), std::declval<__type79>())) __type80;
    typedef decltype(pythonic::operator_::add(std::declval<__type80>(), std::declval<__type43>())) __type82;
    typedef decltype(pythonic::operator_::sub(std::declval<__type82>(), std::declval<__type78>())) __type84;
    typedef decltype(pythonic::operator_::sub(std::declval<__type43>(), std::declval<__type78>())) __type87;
    typedef decltype(std::declval<__type47>()(std::declval<__type84>(), std::declval<__type87>())) __type88;
    typedef typename __combined<__type36,__type38,__type46>::type __type89;
    typedef __type89 __type90;
    typedef decltype(std::declval<__type90>()[std::declval<__type78>()]) __type92;
    typedef decltype(pythonic::operator_::mul(std::declval<__type88>(), std::declval<__type92>())) __type93;
    typedef typename __combined<__type70,__type93>::type __type94;
    typedef __type94 __type95;
    typedef container<typename std::remove_reference<__type95>::type> __type96;
    typedef typename __combined<__type36,__type38,__type44,__type46,__type96>::type __type97;
    typedef typename pythonic::assignable<__type97>::type __type98;
    typedef typename pythonic::assignable<__type94>::type __type99;
    typedef typename pythonic::assignable<__type37>::type __type100;
    typedef __type97 __type101;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type59>::type::iterator>::value_type>::type __type104;
    typedef __type104 __type105;
    typedef decltype(std::declval<__type101>()[std::declval<__type105>()]) __type106;
    typedef decltype(std::declval<__type63>()[std::declval<__type105>()]) __type110;
    typedef decltype(pythonic::operator_::sub(std::declval<__type21>(), std::declval<__type110>())) __type111;
    typedef decltype(pythonic::operator_::sub(std::declval<__type20>(), std::declval<__type105>())) __type114;
    typedef decltype(pythonic::operator_::add(std::declval<__type111>(), std::declval<__type114>())) __type115;
    typedef decltype(std::declval<__type47>()(std::declval<__type115>(), std::declval<__type114>())) __type119;
    typedef decltype(pythonic::operator_::mul(std::declval<__type106>(), std::declval<__type119>())) __type120;
    typedef typename __combined<__type100,__type120>::type __type121;
    typedef typename pythonic::assignable<__type121>::type __type122;
    __type17 n_ = n;
    __type18 m_ = m;
    if (pythonic::operator_::lt(m_, n_))
    {
      typename pythonic::assignable_noescape<decltype(pythonic::types::make_tuple(n_, m_))>::type __tuple0 = pythonic::types::make_tuple(n_, m_);
      m_ = std::get<0>(__tuple0);
      n_ = std::get<1>(__tuple0);
    }
    typename pythonic::assignable_noescape<decltype(pythonic::operator_::functor::floordiv()(m_, g))>::type mg = pythonic::operator_::functor::floordiv()(m_, g);
    typename pythonic::assignable_noescape<decltype(pythonic::operator_::add(n_, pythonic::operator_::functor::floordiv()(pythonic::operator_::sub(mg, h), mg)))>::type lxj = pythonic::operator_::add(n_, pythonic::operator_::functor::floordiv()(pythonic::operator_::sub(mg, h), mg));
    typename pythonic::assignable_noescape<decltype(pythonic::builtins::pythran::functor::static_list{}(pythonic::builtins::functor::map{}(pythonic::functools::functor::partial{}(_count_paths_outside_method_lambda0(), h, mg, pythonic::operator_::functor::floordiv()(n_, g)), pythonic::builtins::functor::range{}(lxj))))>::type xj = pythonic::builtins::pythran::functor::static_list{}(pythonic::builtins::functor::map{}(pythonic::functools::functor::partial{}(_count_paths_outside_method_lambda0(), h, mg, pythonic::operator_::functor::floordiv()(n_, g)), pythonic::builtins::functor::range{}(lxj)));
    if (pythonic::operator_::eq(lxj, 0L))
    {
      return pythonic::scipy::special::functor::binom{}(pythonic::operator_::add(m_, n_), n_);
    }
    else
    {
      __type98 B = pythonic::numpy::functor::zeros{}(lxj);
      std::get<0>(B) = 1L;
      {
        long  __target2668469377136 = lxj;
        for (long  j=1L; j < __target2668469377136; j += 1L)
        {
          __type99 Bj = pythonic::scipy::special::functor::binom{}(pythonic::operator_::add(xj.fast(j), j), j);
          {
            long  __target2668469373248 = j;
            for (long  i=0L; i < __target2668469373248; i += 1L)
            {
              Bj -= pythonic::operator_::mul(pythonic::scipy::special::functor::binom{}(pythonic::operator_::sub(pythonic::operator_::add(pythonic::operator_::sub(xj.fast(j), xj.fast(i)), j), i), pythonic::operator_::sub(j, i)), B.fast(i));
            }
          }
          B.fast(j) = Bj;
        }
      }
      __type122 num_paths = 0L;
      {
        long  __target2668469371472 = lxj;
        for (long  j_=0L; j_ < __target2668469371472; j_ += 1L)
        {
          num_paths += pythonic::operator_::mul(B.fast(j_), pythonic::scipy::special::functor::binom{}(pythonic::operator_::add(pythonic::operator_::sub(m_, xj.fast(j_)), pythonic::operator_::sub(n_, j_)), pythonic::operator_::sub(n_, j_)));
        }
      }
      return num_paths;
    }
  }
}
#include <pythonic/python/exception_handler.hpp>
#ifdef ENABLE_PYTHON_MODULE
inline
typename __pythran__stats_pythran::siegelslopes::type<pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::ndarray<double,pythonic::types::pshape<long>>, pythonic::types::str>::result_type siegelslopes0(pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& y, pythonic::types::ndarray<double,pythonic::types::pshape<long>>&& x, pythonic::types::str&& method) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::siegelslopes()(y, x, method);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::siegelslopes::type<pythonic::types::ndarray<float,pythonic::types::pshape<long>>, pythonic::types::ndarray<float,pythonic::types::pshape<long>>, pythonic::types::str>::result_type siegelslopes1(pythonic::types::ndarray<float,pythonic::types::pshape<long>>&& y, pythonic::types::ndarray<float,pythonic::types::pshape<long>>&& x, pythonic::types::str&& method) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::siegelslopes()(y, x, method);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_count_paths_outside_method::type<npy_int64, npy_int64, npy_int64, npy_int64>::result_type _count_paths_outside_method0(npy_int64&& m, npy_int64&& n, npy_int64&& g, npy_int64&& h) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_count_paths_outside_method()(m, n, g, h);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_compute_prob_outside_square::type<npy_int64, npy_int64>::result_type _compute_prob_outside_square0(npy_int64&& n, npy_int64&& h) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_compute_prob_outside_square()(n, h);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_compute_outer_prob_inside_method::type<npy_int64, npy_int64, npy_int64, npy_int64>::result_type _compute_outer_prob_inside_method0(npy_int64&& m, npy_int64&& n, npy_int64&& g, npy_int64&& h) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_compute_outer_prob_inside_method()(m, n, g, h);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_a_ij_Aij_Dij2::type<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>::result_type _a_ij_Aij_Dij20(pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_a_ij_Aij_Dij2()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_a_ij_Aij_Dij2::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>::result_type _a_ij_Aij_Dij21(pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_a_ij_Aij_Dij2()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_a_ij_Aij_Dij2::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>::result_type _a_ij_Aij_Dij22(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_a_ij_Aij_Dij2()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_a_ij_Aij_Dij2::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>::result_type _a_ij_Aij_Dij23(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_a_ij_Aij_Dij2()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_discordant_pairs::type<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>::result_type _discordant_pairs0(pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_discordant_pairs()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_discordant_pairs::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>::result_type _discordant_pairs1(pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_discordant_pairs()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_discordant_pairs::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>::result_type _discordant_pairs2(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_discordant_pairs()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_discordant_pairs::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>::result_type _discordant_pairs3(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_discordant_pairs()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_concordant_pairs::type<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>::result_type _concordant_pairs0(pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_concordant_pairs()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_concordant_pairs::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>::result_type _concordant_pairs1(pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_concordant_pairs()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_concordant_pairs::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>::result_type _concordant_pairs2(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_concordant_pairs()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_concordant_pairs::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>::result_type _concordant_pairs3(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& A) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_concordant_pairs()(A);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_Dij::type<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>, long, long>::result_type _Dij0(pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& A, long&& i, long&& j) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_Dij()(A, i, j);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_Dij::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>, long, long>::result_type _Dij1(pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& A, long&& i, long&& j) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_Dij()(A, i, j);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_Dij::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, long, long>::result_type _Dij2(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& A, long&& i, long&& j) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_Dij()(A, i, j);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_Dij::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, long, long>::result_type _Dij3(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& A, long&& i, long&& j) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_Dij()(A, i, j);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_Aij::type<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>, long, long>::result_type _Aij0(pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>&& A, long&& i, long&& j) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_Aij()(A, i, j);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_Aij::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>, long, long>::result_type _Aij1(pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>&& A, long&& i, long&& j) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_Aij()(A, i, j);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_Aij::type<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>, long, long>::result_type _Aij2(pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>&& A, long&& i, long&& j) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_Aij()(A, i, j);
                                PyEval_RestoreThread(_save);
                                return res;
                            }
                            catch(...) {
                                PyEval_RestoreThread(_save);
                                throw;
                            }
                            ;
}
inline
typename __pythran__stats_pythran::_Aij::type<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>, long, long>::result_type _Aij3(pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>&& A, long&& i, long&& j) 
{
  
                            PyThreadState *_save = PyEval_SaveThread();
                            try {
                                auto res = __pythran__stats_pythran::_Aij()(A, i, j);
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
__pythran_wrap_siegelslopes0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[3+1];
    
    char const* keywords[] = {"y", "x", "method",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]))
        return to_python(siegelslopes0(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[0]), from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap_siegelslopes1(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[3+1];
    
    char const* keywords[] = {"y", "x", "method",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<float,pythonic::types::pshape<long>>>(args_obj[0]) && is_convertible<pythonic::types::ndarray<float,pythonic::types::pshape<long>>>(args_obj[1]) && is_convertible<pythonic::types::str>(args_obj[2]))
        return to_python(siegelslopes1(from_python<pythonic::types::ndarray<float,pythonic::types::pshape<long>>>(args_obj[0]), from_python<pythonic::types::ndarray<float,pythonic::types::pshape<long>>>(args_obj[1]), from_python<pythonic::types::str>(args_obj[2])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__count_paths_outside_method0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[4+1];
    
    char const* keywords[] = {"m", "n", "g", "h",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3]))
        return nullptr;
    if(is_convertible<npy_int64>(args_obj[0]) && is_convertible<npy_int64>(args_obj[1]) && is_convertible<npy_int64>(args_obj[2]) && is_convertible<npy_int64>(args_obj[3]))
        return to_python(_count_paths_outside_method0(from_python<npy_int64>(args_obj[0]), from_python<npy_int64>(args_obj[1]), from_python<npy_int64>(args_obj[2]), from_python<npy_int64>(args_obj[3])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__compute_prob_outside_square0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[2+1];
    
    char const* keywords[] = {"n", "h",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OO",
                                     (char**)keywords , &args_obj[0], &args_obj[1]))
        return nullptr;
    if(is_convertible<npy_int64>(args_obj[0]) && is_convertible<npy_int64>(args_obj[1]))
        return to_python(_compute_prob_outside_square0(from_python<npy_int64>(args_obj[0]), from_python<npy_int64>(args_obj[1])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__compute_outer_prob_inside_method0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[4+1];
    
    char const* keywords[] = {"m", "n", "g", "h",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2], &args_obj[3]))
        return nullptr;
    if(is_convertible<npy_int64>(args_obj[0]) && is_convertible<npy_int64>(args_obj[1]) && is_convertible<npy_int64>(args_obj[2]) && is_convertible<npy_int64>(args_obj[3]))
        return to_python(_compute_outer_prob_inside_method0(from_python<npy_int64>(args_obj[0]), from_python<npy_int64>(args_obj[1]), from_python<npy_int64>(args_obj[2]), from_python<npy_int64>(args_obj[3])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__a_ij_Aij_Dij20(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[0]))
        return to_python(_a_ij_Aij_Dij20(from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__a_ij_Aij_Dij21(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[0]))
        return to_python(_a_ij_Aij_Dij21(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__a_ij_Aij_Dij22(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]))
        return to_python(_a_ij_Aij_Dij22(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__a_ij_Aij_Dij23(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]))
        return to_python(_a_ij_Aij_Dij23(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__discordant_pairs0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[0]))
        return to_python(_discordant_pairs0(from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__discordant_pairs1(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[0]))
        return to_python(_discordant_pairs1(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__discordant_pairs2(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]))
        return to_python(_discordant_pairs2(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__discordant_pairs3(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]))
        return to_python(_discordant_pairs3(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__concordant_pairs0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[0]))
        return to_python(_concordant_pairs0(from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__concordant_pairs1(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[0]))
        return to_python(_concordant_pairs1(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__concordant_pairs2(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]))
        return to_python(_concordant_pairs2(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__concordant_pairs3(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[1+1];
    
    char const* keywords[] = {"A",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "O",
                                     (char**)keywords , &args_obj[0]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]))
        return to_python(_concordant_pairs3(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__Dij0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[3+1];
    
    char const* keywords[] = {"A", "i", "j",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<long>(args_obj[1]) && is_convertible<long>(args_obj[2]))
        return to_python(_Dij0(from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<long>(args_obj[1]), from_python<long>(args_obj[2])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__Dij1(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[3+1];
    
    char const* keywords[] = {"A", "i", "j",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<long>(args_obj[1]) && is_convertible<long>(args_obj[2]))
        return to_python(_Dij1(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<long>(args_obj[1]), from_python<long>(args_obj[2])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__Dij2(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[3+1];
    
    char const* keywords[] = {"A", "i", "j",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<long>(args_obj[1]) && is_convertible<long>(args_obj[2]))
        return to_python(_Dij2(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<long>(args_obj[1]), from_python<long>(args_obj[2])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__Dij3(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[3+1];
    
    char const* keywords[] = {"A", "i", "j",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<long>(args_obj[1]) && is_convertible<long>(args_obj[2]))
        return to_python(_Dij3(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<long>(args_obj[1]), from_python<long>(args_obj[2])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__Aij0(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[3+1];
    
    char const* keywords[] = {"A", "i", "j",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<long>(args_obj[1]) && is_convertible<long>(args_obj[2]))
        return to_python(_Aij0(from_python<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<long>(args_obj[1]), from_python<long>(args_obj[2])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__Aij1(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[3+1];
    
    char const* keywords[] = {"A", "i", "j",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<long>(args_obj[1]) && is_convertible<long>(args_obj[2]))
        return to_python(_Aij1(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<long>(args_obj[1]), from_python<long>(args_obj[2])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__Aij2(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[3+1];
    
    char const* keywords[] = {"A", "i", "j",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2]))
        return nullptr;
    if(is_convertible<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]) && is_convertible<long>(args_obj[1]) && is_convertible<long>(args_obj[2]))
        return to_python(_Aij2(from_python<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>(args_obj[0]), from_python<long>(args_obj[1]), from_python<long>(args_obj[2])));
    else {
        return nullptr;
    }
}

static PyObject *
__pythran_wrap__Aij3(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject* args_obj[3+1];
    
    char const* keywords[] = {"A", "i", "j",  nullptr};
    if(! PyArg_ParseTupleAndKeywords(args, kw, "OOO",
                                     (char**)keywords , &args_obj[0], &args_obj[1], &args_obj[2]))
        return nullptr;
    if(is_convertible<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]) && is_convertible<long>(args_obj[1]) && is_convertible<long>(args_obj[2]))
        return to_python(_Aij3(from_python<pythonic::types::numpy_texpr<pythonic::types::ndarray<double,pythonic::types::pshape<long,long>>>>(args_obj[0]), from_python<long>(args_obj[1]), from_python<long>(args_obj[2])));
    else {
        return nullptr;
    }
}

            static PyObject *
            __pythran_wrapall_siegelslopes(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap_siegelslopes0(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap_siegelslopes1(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "siegelslopes", "\n""    - siegelslopes(float64[:], float64[:], str)\n""    - siegelslopes(float32[:], float32[:], str)", args, kw);
                });
            }


            static PyObject *
            __pythran_wrapall__count_paths_outside_method(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__count_paths_outside_method0(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_count_paths_outside_method", "\n""    - _count_paths_outside_method(int64, int64, int64, int64)", args, kw);
                });
            }


            static PyObject *
            __pythran_wrapall__compute_prob_outside_square(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__compute_prob_outside_square0(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_compute_prob_outside_square", "\n""    - _compute_prob_outside_square(int64, int64)", args, kw);
                });
            }


            static PyObject *
            __pythran_wrapall__compute_outer_prob_inside_method(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__compute_outer_prob_inside_method0(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_compute_outer_prob_inside_method", "\n""    - _compute_outer_prob_inside_method(int64, int64, int64, int64)", args, kw);
                });
            }


            static PyObject *
            __pythran_wrapall__a_ij_Aij_Dij2(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__a_ij_Aij_Dij20(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__a_ij_Aij_Dij21(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__a_ij_Aij_Dij22(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__a_ij_Aij_Dij23(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_a_ij_Aij_Dij2", "\n""    - _a_ij_Aij_Dij2(int[:,:])\n""    - _a_ij_Aij_Dij2(float[:,:])", args, kw);
                });
            }


            static PyObject *
            __pythran_wrapall__discordant_pairs(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__discordant_pairs0(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__discordant_pairs1(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__discordant_pairs2(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__discordant_pairs3(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_discordant_pairs", "\n""    - _discordant_pairs(int[:,:])\n""    - _discordant_pairs(float[:,:])", args, kw);
                });
            }


            static PyObject *
            __pythran_wrapall__concordant_pairs(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__concordant_pairs0(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__concordant_pairs1(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__concordant_pairs2(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__concordant_pairs3(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_concordant_pairs", "\n""    - _concordant_pairs(int[:,:])\n""    - _concordant_pairs(float[:,:])", args, kw);
                });
            }


            static PyObject *
            __pythran_wrapall__Dij(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__Dij0(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__Dij1(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__Dij2(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__Dij3(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_Dij", "\n""    - _Dij(int[:,:], int, int)\n""    - _Dij(float[:,:], int, int)", args, kw);
                });
            }


            static PyObject *
            __pythran_wrapall__Aij(PyObject *self, PyObject *args, PyObject *kw)
            {
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {

if(PyObject* obj = __pythran_wrap__Aij0(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__Aij1(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__Aij2(self, args, kw))
    return obj;
PyErr_Clear();


if(PyObject* obj = __pythran_wrap__Aij3(self, args, kw))
    return obj;
PyErr_Clear();

                return pythonic::python::raise_invalid_argument(
                               "_Aij", "\n""    - _Aij(int[:,:], int, int)\n""    - _Aij(float[:,:], int, int)", args, kw);
                });
            }


static PyMethodDef Methods[] = {
    {
    "siegelslopes",
    (PyCFunction)__pythran_wrapall_siegelslopes,
    METH_VARARGS | METH_KEYWORDS,
    "Supported prototypes:\n""\n""    - siegelslopes(float64[:], float64[:], str)\n""    - siegelslopes(float32[:], float32[:], str)"},{
    "_count_paths_outside_method",
    (PyCFunction)__pythran_wrapall__count_paths_outside_method,
    METH_VARARGS | METH_KEYWORDS,
    "Count the number of paths that pass outside the specified diagonal.\n""\n""    Supported prototypes:\n""\n""    - _count_paths_outside_method(int64, int64, int64, int64)\n""\n""    Parameters\n""    ----------\n""    m : integer\n""        m > 0\n""    n : integer\n""        n > 0\n""    g : integer\n""        g is greatest common divisor of m and n\n""    h : integer\n""        0 <= h <= lcm(m,n)\n""\n""    Returns\n""    -------\n""    p : float\n""        The number of paths that go low.\n""        The calculation may overflow - check for a finite answer.\n""\n""    Notes\n""    -----\n""    Count the integer lattice paths from (0, 0) to (m, n), which at some\n""    point (x, y) along the path, satisfy:\n""      m*y <= n*x - h*g\n""    The paths make steps of size +1 in either positive x or positive y\n""    directions.\n""\n""    We generally follow Hodges' treatment of Drion/Gnedenko/Korolyuk.\n""    Hodges, J.L. Jr.,\n""    \"The Significance Probability of the Smirnov Two-Sample Test,\"\n""    Arkiv fiur Matematik, 3, No. 43 (1958), 469-86.\n""\n"""},{
    "_compute_prob_outside_square",
    (PyCFunction)__pythran_wrapall__compute_prob_outside_square,
    METH_VARARGS | METH_KEYWORDS,
    "\n""Compute the proportion of paths that pass outside the two diagonal lines.\n""\n""Supported prototypes:\n""\n""- _compute_prob_outside_square(int64, int64)\n""\n""Parameters\n""----------\n""n : integer\n""    n > 0\n""h : integer\n""    0 <= h <= n\n""\n""Returns\n""-------\n""p : float\n""    The proportion of paths that pass outside the lines x-y = +/-h.\n""\n"""},{
    "_compute_outer_prob_inside_method",
    (PyCFunction)__pythran_wrapall__compute_outer_prob_inside_method,
    METH_VARARGS | METH_KEYWORDS,
    "\n""Count the proportion of paths that do not stay strictly inside two\n""diagonal lines.\n""\n""Supported prototypes:\n""\n""- _compute_outer_prob_inside_method(int64, int64, int64, int64)\n""\n""Parameters\n""----------\n""m : integer\n""    m > 0\n""n : integer\n""    n > 0\n""g : integer\n""    g is greatest common divisor of m and n\n""h : integer\n""    0 <= h <= lcm(m,n)\n""\n""Returns\n""-------\n""p : float\n""    The proportion of paths that do not stay inside the two lines.\n""\n""The classical algorithm counts the integer lattice paths from (0, 0)\n""to (m, n) which satisfy |x/m - y/n| < h / lcm(m, n).\n""The paths make steps of size +1 in either positive x or positive y\n""directions.\n""We are, however, interested in 1 - proportion to computes p-values,\n""so we change the recursion to compute 1 - p directly while staying\n""within the \"inside method\" a described by Hodges.\n""\n""We generally follow Hodges' treatment of Drion/Gnedenko/Korolyuk.\n""Hodges, J.L. Jr.,\n""\"The Significance Probability of the Smirnov Two-Sample Test,\"\n""Arkiv fiur Matematik, 3, No. 43 (1958), 469-86.\n""\n""For the recursion for 1-p see\n""Viehmann, T.: \"Numerically more stable computation of the p-values\n""for the two-sample Kolmogorov-Smirnov test,\" arXiv: 2102.08037\n""\n"""},{
    "_a_ij_Aij_Dij2",
    (PyCFunction)__pythran_wrapall__a_ij_Aij_Dij2,
    METH_VARARGS | METH_KEYWORDS,
    "A term that appears in the ASE of Kendall's tau and Somers' D.\n""\n""    Supported prototypes:\n""\n""    - _a_ij_Aij_Dij2(int[:,:])\n""    - _a_ij_Aij_Dij2(float[:,:])"},{
    "_discordant_pairs",
    (PyCFunction)__pythran_wrapall__discordant_pairs,
    METH_VARARGS | METH_KEYWORDS,
    "Twice the number of discordant pairs, excluding ties.\n""\n""    Supported prototypes:\n""\n""    - _discordant_pairs(int[:,:])\n""    - _discordant_pairs(float[:,:])"},{
    "_concordant_pairs",
    (PyCFunction)__pythran_wrapall__concordant_pairs,
    METH_VARARGS | METH_KEYWORDS,
    "Twice the number of concordant pairs, excluding ties.\n""\n""    Supported prototypes:\n""\n""    - _concordant_pairs(int[:,:])\n""    - _concordant_pairs(float[:,:])"},{
    "_Dij",
    (PyCFunction)__pythran_wrapall__Dij,
    METH_VARARGS | METH_KEYWORDS,
    "Sum of lower-left and upper-right blocks of contingency table.\n""\n""    Supported prototypes:\n""\n""    - _Dij(int[:,:], int, int)\n""    - _Dij(float[:,:], int, int)"},{
    "_Aij",
    (PyCFunction)__pythran_wrapall__Aij,
    METH_VARARGS | METH_KEYWORDS,
    "Sum of upper-left and lower right blocks of contingency table.\n""\n""    Supported prototypes:\n""\n""    - _Aij(int[:,:], int, int)\n""    - _Aij(float[:,:], int, int)"},
    {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_stats_pythran",            /* m_name */
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
PYTHRAN_MODULE_INIT(_stats_pythran)(void)
#ifndef _WIN32
__attribute__ ((visibility("default")))
#if defined(GNUC) && !defined(__clang__)
__attribute__ ((externally_visible))
#endif
#endif
;
PyMODINIT_FUNC
PYTHRAN_MODULE_INIT(_stats_pythran)(void) {
    import_array()
    #if PY_MAJOR_VERSION >= 3
    PyObject* theModule = PyModule_Create(&moduledef);
    #else
    PyObject* theModule = Py_InitModule3("_stats_pythran",
                                         Methods,
                                         ""
    );
    #endif
    if(! theModule)
        PYTHRAN_RETURN;
    PyObject * theDoc = Py_BuildValue("(sss)",
                                      "0.12.0",
                                      "2022-10-15 12:07:28.565801",
                                      "7a47ec42678a84c5a3b15e89f802e251ea40d31fa283ebfe78ea6dc6d0b7da4a");
    if(! theDoc)
        PYTHRAN_RETURN;
    PyModule_AddObject(theModule,
                       "__pythran__",
                       theDoc);


    PYTHRAN_RETURN;
}

#endif