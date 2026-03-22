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
 *  Copyright (C) 2002-2021 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cctype>
#include "ducc0/infra/string_utils.h"
#include "ducc0/infra/error_handling.h"

namespace ducc0 {

namespace detail_string_utils {

using namespace std;

string trim (const string &orig)
  {
  string::size_type p1=orig.find_first_not_of(" \t");
  if (p1==string::npos) return "";
  string::size_type p2=orig.find_last_not_of(" \t");
  return orig.substr(p1,p2-p1+1);
  }

template<typename T> string dataToString (const T &x)
  {
  ostringstream strstrm;
  strstrm << x;
  return trim(strstrm.str());
  }

template<> string dataToString (const bool &x)
  { return x ? "T" : "F"; }
template<> string dataToString (const string &x)
  { return trim(x); }
template<> string dataToString (const float &x)
  {
  ostringstream strstrm;
  strstrm << setprecision(8) << x;
  return trim(strstrm.str());
  }
template<> string dataToString (const double &x)
  {
  ostringstream strstrm;
  strstrm << setprecision(16) << x;
  return trim(strstrm.str());
  }
template<> string dataToString (const long double &x)
  {
  ostringstream strstrm;
  strstrm << setprecision(25) << x;
  return trim(strstrm.str());
  }

template string dataToString (const signed char &x);
template string dataToString (const unsigned char &x);
template string dataToString (const short &x);
template string dataToString (const unsigned short &x);
template string dataToString (const int &x);
template string dataToString (const unsigned int &x);
template string dataToString (const long &x);
template string dataToString (const unsigned long &x);
template string dataToString (const long long &x);
template string dataToString (const unsigned long long &x);

string intToString(int64_t x, size_t width)
  {
  ostringstream strstrm;
  (x>=0) ? strstrm << setw(width) << setfill('0') << x
         : strstrm << "-" << setw(width-1) << setfill('0') << -x;
  string res = strstrm.str();
  MR_assert(res.size()==width,"number too large");
  return trim(res);
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

template<> string stringToData (const string &x)
  { return trim(x); }

template<> bool stringToData (const string &x)
  {
  const char *fval[] = {"f","n","false",".false."};
  const char *tval[] = {"t","y","true",".true."};
  for (size_t i=0; i< sizeof(fval)/sizeof(fval[0]); ++i)
    if (equal_nocase(x,fval[i])) return false;
  for (size_t i=0; i< sizeof(tval)/sizeof(tval[0]); ++i)
    if (equal_nocase(x,tval[i])) return true;
  MR_fail("conversion error in stringToData<bool>(",x,")");
  }

template signed char stringToData (const string &x);
template unsigned char stringToData (const string &x);
template short stringToData (const string &x);
template unsigned short stringToData (const string &x);
template int  stringToData (const string &x);
template unsigned int stringToData (const string &x);
template long stringToData (const string &x);
template unsigned long stringToData (const string &x);
template long long stringToData (const string &x);
template unsigned long long stringToData (const string &x);
template float stringToData (const string &x);
template double stringToData (const string &x);
template long double stringToData (const string &x);

bool equal_nocase (const string &a, const string &b)
  {
  if (a.size()!=b.size()) return false;
  for (size_t m=0; m<a.size(); ++m)
    if (std::tolower(a[m])!=std::tolower(b[m])) return false;
  return true;
  }

string tolower(const string &input)
  {
  string result=input;
  for (size_t m=0; m<result.size(); ++m)
    result[m]=char(std::tolower(result[m]));
  return result;
  }

namespace {

template<typename T> vector<T> split (istream &stream)
  {
  vector<T> list;
  while (stream)
    {
    string word;
    stream >> word;
    MR_assert (stream||stream.eof(),
      "error while splitting stream into components");
    if (stream) list.push_back(stringToData<T>(word));
    }
  return list;
  }

} // unnamed namespace

template<typename T> vector<T> split (const string &inp)
  {
  istringstream is(inp);
  return split<T>(is);
  }

template vector<string> split (const string &inp);
template vector<float> split (const string &inp);
template vector<double> split (const string &inp);
template vector<int> split (const string &inp);
template vector<long> split (const string &inp);

vector<string> tokenize (const string &inp, char delim)
  {
  istringstream stream(inp);
  string token;
  vector<string> list;
  while (getline(stream,token,delim))
    list.push_back(token);
  return list;
  }

vector<string> parse_words_from_file (const string &filename)
  {
  vector<string> words;
  ifstream inp(filename.c_str());
  MR_assert (inp,"Could not open file '", filename, "'.");
  while (inp)
    {
    string word;
    inp>>word;
    word=trim(word);
    if (word.empty()) words.push_back(word);
    }
  return words;
  }

}}
