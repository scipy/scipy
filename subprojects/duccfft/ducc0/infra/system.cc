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

/* Copyright (C) 2019-2021 Max-Planck-Society
   Author: Martin Reinecke */

#include <regex>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

#include "ducc0/infra/error_handling.h"
#include "ducc0/infra/system.h"
#include "ducc0/infra/string_utils.h"

namespace ducc0 {

namespace detail_system {

using namespace std;

string fileToString(const string &fname)
  {
  ifstream inp(fname);
  stringstream sbuf;
  sbuf << inp.rdbuf();
  return sbuf.str();
  }

template<typename T> T find(const string &s, const string &pattern)
  {
  regex re(pattern);
  sregex_iterator it(s.begin(), s.end(), re);
  sregex_iterator it_end;
  MR_assert (it!=it_end, "did not find pattern '", pattern, "'");
  return stringToData<T>(it->str(1));
  }

size_t getProcessInfo(const string &quantity)
  {
  string text = fileToString("/proc/self/status");
  return find<size_t>(text, quantity + R"(:\s+(\d+))");
  }

size_t getMemInfo(const string &quantity)
  {
  string text = fileToString("/proc/meminfo");
  return find<size_t>(text, quantity + R"(:\s+(\d+) kB)");
  }

size_t usable_memory()
  {
  string text = fileToString("/proc/meminfo");
  size_t MemTotal = find<size_t>(text, R"(MemTotal:\s+(\d+) kB)");
  size_t Committed = find<size_t>(text, R"(Committed_AS:\s+(\d+) kB)");
  return MemTotal-Committed;
  }

}}
