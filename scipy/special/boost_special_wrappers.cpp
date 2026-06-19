#include <Python.h>

#include "boost_special_functions.h"
#include "boost_special_wrappers.h"

double special_bdtrik(double y, double n, double p) { return bdtrik_double(y, n, p); }
