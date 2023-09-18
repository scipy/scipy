/* This file selects the right generated file of `__stub_FUNCTION' macros
   based on the architecture being compiled for.  */

#include <bits/wordsize.h>

#if __WORDSIZE == 32
# include <gnu/stubs-32.h>
#elif __WORDSIZE == 64
# include <gnu/stubs-64.h>
#else
# error "unexpected value for __WORDSIZE macro"
#endif
