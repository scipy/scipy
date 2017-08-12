#if !defined(MSVC_9_FMAX) && defined(_MSC_VER) && _MSC_VER <= 1600
#define MSVC_9_FMAX
#include <algorithm>
#include <float.h>
extern double fmax(double, double);
#endif