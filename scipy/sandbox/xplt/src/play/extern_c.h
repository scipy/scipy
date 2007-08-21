/* get this ugliness out of important header files */
#ifndef BEGIN_EXTERN_C
# if defined(__cplusplus) || defined(c_plusplus)
#  define BEGIN_EXTERN_C extern "C" {
#  define END_EXTERN_C }
# else
#  define BEGIN_EXTERN_C
#  define END_EXTERN_C
# endif
#endif
BEGIN_EXTERN_C
