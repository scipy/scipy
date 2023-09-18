#ifndef NUMPY_BINARY_FUNC_NAME
#error NUMPY_BINARY_FUNC_NAME undefined
#endif
#ifndef NUMPY_BINARY_FUNC_SYM
#error NUMPY_BINARY_FUNC_SYM undefined
#endif

template <class E0, class E1>
typename std::enable_if<
    types::valid_numop_parameters<typename std::decay<E0>::type,
                                  typename std::decay<E1>::type>::value,
    types::numpy_expr<NUMPY_BINARY_FUNC_SYM,
                      typename types::adapt_type<E0, E1>::type,
                      typename types::adapt_type<E1, E0>::type>>::type
NUMPY_BINARY_FUNC_NAME(E0 &&self, E1 &&other);

#undef NUMPY_BINARY_FUNC_NAME
#undef NUMPY_BINARY_FUNC_SYM
