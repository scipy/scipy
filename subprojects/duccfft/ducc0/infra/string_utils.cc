/* SPDX-License-Identifier: BSD-3-Clause OR GPL-2.0-or-later */

/*
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.
* Neither the name of the copyright holder nor the names of its contributors may
  be used to endorse or promote products derived from this software without
  specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  This file contains the implementation of various convenience functions
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2002-2026 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cctype>
#include "ducc0/infra/string_utils.h"
#include "ducc0/infra/error_handling.h"

namespace DUCC0_NAMESPACE {

namespace detail_string_utils {

using namespace std;

string trim (const string &orig)
  {
  string::size_type p1=orig.find_first_not_of(" \t");
  if (p1==string::npos) return "";
  string::size_type p2=orig.find_last_not_of(" \t");
  return orig.substr(p1,p2-p1+1);
  }

template<typename T> T stringToData (const string &x)
  {
  istringstream strstrm(x);
  T value;
  strstrm >> value;
  bool ok = bool(strstrm);
  if (ok)
    {
    string rest;
    strstrm >> rest;
    ok = rest.empty();
    }
  MR_assert(ok, "could not convert '", x, "' to desired data type.");
  return value;
  }

template long stringToData (const string &x);

bool equal_nocase (const string &a, const string &b)
  {
  if (a.size()!=b.size()) return false;
  for (size_t m=0; m<a.size(); ++m)
    if (std::tolower(a[m])!=std::tolower(b[m])) return false;
  return true;
  }

}}
