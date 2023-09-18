/* Emit optimization information as JSON files.
   Copyright (C) 2018-2022 Free Software Foundation, Inc.
   Contributed by David Malcolm <dmalcolm@redhat.com>.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3, or (at your option) any later
version.

GCC is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with GCC; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.  */

#ifndef GCC_OPTINFO_EMIT_JSON_H
#define GCC_OPTINFO_EMIT_JSON_H

#include "json.h"

class optinfo;

/* A class for writing out optimization records in JSON format.  */

class optrecord_json_writer
{
public:
  optrecord_json_writer ();
  ~optrecord_json_writer ();
  void write () const;
  void add_record (const optinfo *optinfo);
  void pop_scope ();

  void add_record (json::object *obj);
  json::object *impl_location_to_json (dump_impl_location_t loc);
  json::object *location_to_json (location_t loc);
  json::object *profile_count_to_json (profile_count count);
  json::string *get_id_value_for_pass (opt_pass *pass);
  json::object *pass_to_json (opt_pass *pass);
  json::value *inlining_chain_to_json (location_t loc);
  json::object *optinfo_to_json (const optinfo *optinfo);
  void add_pass_list (json::array *arr, opt_pass *pass);

 private:
  /* The root value for the JSON file.
     Currently the JSON values are stored in memory, and flushed when the
     compiler exits.  It would probably be better to simply write out
     the JSON as we go.  */
  json::array *m_root_tuple;

  /* The currently open scopes, for expressing nested optimization records.  */
  auto_vec<json::array *> m_scopes;
};

#endif /* #ifndef GCC_OPTINFO_EMIT_JSON_H */
