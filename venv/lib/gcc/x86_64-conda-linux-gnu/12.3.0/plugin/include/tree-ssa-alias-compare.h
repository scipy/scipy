/* Comparsion of AO ref.
   Copyright (C) 2020-2022 Free Software Foundation, Inc.

   This file is part of GCC.

   GCC is free software; you can redistribute it and/or modify
   under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   GCC is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with GCC; see the file COPYING3.  If not see
   <http://www.gnu.org/licenses/>.  */

#ifndef TREE_SSA_ALIAS_COMPARE_H
#define TREE_SSA_ALIAS_COMPARE_H

class operand_compare;
/* A class aggregating all connections and semantic equivalents
   for a given pair of semantic function candidates.  */
class ao_compare : public operand_compare
{
  public:
  enum ao_ref_diff
  {
    SEMANTICS = 1,
    BASE_ALIAS_SET = 2,
    REF_ALIAS_SET = 4,
    ACCESS_PATH = 8,
    DEPENDENCE_CLIQUE = 16
  };
  int compare_ao_refs (ao_ref *ref1, ao_ref *ref2, bool lto_streaming_safe,
		       bool tbaa);
  void hash_ao_ref (ao_ref *ref, bool lto_streaming_safe, bool tbaa,
		    inchash::hash &hstate);
};

#endif
