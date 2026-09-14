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
 *  This file is part of the MR utility library.
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/** \file ducc0/infra/string_utils.h
 *
 *  \copyright Copyright (C) 2019-2021 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef DUCC0_STRING_UTILS_H
#define DUCC0_STRING_UTILS_H

// FIXME: most of this will be superseded by C++20 std::format

#include <string>
#include <vector>
#include <cstdint>
#include <cstddef>

namespace ducc0 {

namespace detail_string_utils {

using namespace std;

/*! \defgroup stringutilsgroup String handling helper functions */
/*! \{ */

/// Returns the string \a orig without leading and trailing whitespace.
string trim (const string &orig);

/// Returns a string containing the text representation of \a x.
/*! Care is taken that no information is lost in the conversion. */
template<typename T> string dataToString(const T &x);
template<> string dataToString (const bool &x);
template<> string dataToString (const string &x);
template<> string dataToString (const float &x);
template<> string dataToString (const double &x);
template<> string dataToString (const long double &x);

/// Returns a string containing the text representation of \a x, padded
/// with leading zeroes to \a width characters.
string intToString(std::int64_t x, std::size_t width);

/// Reads a value of a given datatype from a string.
template<typename T> T stringToData (const string &x);
template<> string stringToData (const string &x);
template<> bool stringToData (const string &x);

/// Case-insensitive string comparison
/*! Returns \a true, if \a a and \a b differ only in capitalisation,
    else \a false. */
bool equal_nocase (const string &a, const string &b);

/// Returns lowercase version of \a input.
string tolower(const string &input);

/// Tries to split \a inp into a white-space separated list of values of
/// type \a T, and appends them to \a list.
template<typename T> inline std::vector<T> split (const string &inp);

/// Breaks the string \a inp into tokens separated by \a delim, and returns them
/// as a vector<string>.
std::vector<string> tokenize (const string &inp, char delim);

/// Breaks the contents of file \a filename into tokens separated by white
/// space, and returns them as a vector<string>.
std::vector<string> parse_words_from_file (const string &filename);

/*! \} */

}

using detail_string_utils::trim;
//using detail_string_utils::intToString;
using detail_string_utils::dataToString;
using detail_string_utils::stringToData;
using detail_string_utils::equal_nocase;
//using detail_string_utils::tolower;
//using detail_string_utils::split;
//using detail_string_utils::tokenize;
//using detail_string_utils::parse_words_from_file;

}

#endif
