// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_INFO_CPP_H_
#define IPX_INFO_CPP_H_

#include <ostream>
#include <string>
#include "ipx_internal.h"

namespace ipx {

// Formatted output of Info. For each member one line is printed in the form
//  info.status                                         1000
//  info.status_ipm                                     1
//  info.status_crossover                               1
//  info.errflag                                        0
//  info.num_var                                        89346
//  info.num_constr                                     39855
//  ...
// The result can be parsed using the Julia function ParseInfo(), which recovers
// an object of a Julia type analogue to struct ipx_info.
std::ostream& operator<<(std::ostream&, const Info&);

// Returns a text representation of each status code.
std::string StatusString(Int status);

}  // namespace ipx

#endif  // IPX_INFO_CPP_H_
