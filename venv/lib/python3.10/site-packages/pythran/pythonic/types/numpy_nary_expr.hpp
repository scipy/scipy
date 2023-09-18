#ifndef NUMPY_NARY_FUNC_NAME
#error NUMPY_NARY_FUNC_NAME undefined
#endif

#ifndef NUMPY_NARY_FUNC_SYM
#error NUMPY_NARY_FUNC_SYM undefined
#endif

#ifndef NUMPY_NARY_RESHAPE_MODE
#define NUMPY_NARY_RESHAPE_MODE adapt_type
#endif

#ifndef NUMPY_NARY_EXTRA_METHOD
#define NUMPY_NARY_EXTRA_METHOD
#endif

namespace functor
{

  template <typename... T>
  auto NUMPY_NARY_FUNC_NAME::
  operator()(T &&... args) const -> typename std::enable_if<
      !types::valid_numexpr_parameters<typename std::decay<T>::type...>::value,
      decltype(NUMPY_NARY_FUNC_SYM(std::forward<T>(args)...))>::type
  {
    return NUMPY_NARY_FUNC_SYM(std::forward<T>(args)...);
  }

  template <class... E>
  typename std::enable_if<
      types::valid_numexpr_parameters<typename std::decay<E>::type...>::value,
      types::numpy_expr<
          NUMPY_NARY_FUNC_NAME,
          typename types::NUMPY_NARY_RESHAPE_MODE<E, E...>::type...>>::type
      NUMPY_NARY_FUNC_NAME::
      operator()(E &&... args) const
  {
    return {std::forward<E>(args)...};
  }
}

#undef NUMPY_NARY_FUNC_NAME
#undef NUMPY_NARY_FUNC_SYM
#undef NUMPY_NARY_RESHAPE_MODE
#undef NUMPY_NARY_EXTRA_METHOD
