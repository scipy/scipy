/* Optimization information.
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

#ifndef GCC_OPTINFO_H
#define GCC_OPTINFO_H

/* An "optinfo" is a bundle of information describing part of an
   optimization, which can be emitted to zero or more of several
   destinations, such as:

   * saved to a file as an "optimization record"

   They are generated in response to calls to the "dump_*" API in
   dumpfile.h; repeated calls to the "dump_*" API are consolidated
   into a pending optinfo instance, with a "dump_*_loc" starting a new
   optinfo instance.

   The data sent to the dump calls are captured within the pending optinfo
   instance as a sequence of optinfo_items.  For example, given:

      if (dump_enabled_p ())
        {
          dump_printf_loc (MSG_MISSED_OPTIMIZATION, vect_location,
                           "not vectorized: live stmt not supported: ");
          dump_gimple_stmt (MSG_MISSED_OPTIMIZATION, TDF_SLIM, stmt, 0);
        }

   the "dump_printf_loc" call begins a new optinfo containing two items:
   (1) a text item containing "not vectorized: live stmt not supported: "
   (2) a gimple item for "stmt"

   Dump destinations are thus able to access rich metadata about the
   items when the optinfo is emitted to them, rather than just having plain
   text.  For example, when saving the above optinfo to a file as an
   "optimization record", the record could capture the source location of
   "stmt" above, rather than just its textual form.

   The currently pending optinfo is emitted and deleted:
   * each time a "dump_*_loc" call occurs (which starts the next optinfo), or
   * when the dump files are changed (at the end of a pass)

   Dumping to an optinfo instance is non-trivial (due to building optinfo_item
   instances), so all usage should be guarded by

     if (optinfo_enabled_p ())

   which is off by default.  */


/* Forward decls.  */
class opt_pass;
class optinfo_item;

/* Return true if any of the active optinfo destinations make use
   of inlining information.
   (if true, then the information is preserved).  */

extern bool optinfo_wants_inlining_info_p ();

/* The various kinds of optinfo.  */

enum optinfo_kind
{
  OPTINFO_KIND_SUCCESS,
  OPTINFO_KIND_FAILURE,
  OPTINFO_KIND_NOTE,
  OPTINFO_KIND_SCOPE
};

extern const char *optinfo_kind_to_string (enum optinfo_kind kind);

class dump_context;

/* A bundle of information describing part of an optimization.  */

class optinfo
{
  friend class dump_context;

 public:
  optinfo (const dump_location_t &loc,
	   enum optinfo_kind kind,
	   opt_pass *pass)
  : m_loc (loc), m_kind (kind), m_pass (pass), m_items ()
  {}
  ~optinfo ();

  const dump_location_t &
  get_dump_location () const { return m_loc; }

  const dump_user_location_t &
  get_user_location () const { return m_loc.get_user_location (); }

  const dump_impl_location_t &
  get_impl_location () const { return m_loc.get_impl_location (); }

  enum optinfo_kind get_kind () const { return m_kind; }
  opt_pass *get_pass () const { return m_pass; }
  unsigned int num_items () const { return m_items.length (); }
  const optinfo_item *get_item (unsigned int i) const { return m_items[i]; }

  location_t get_location_t () const { return m_loc.get_location_t (); }
  profile_count get_count () const { return m_loc.get_count (); }

  void add_item (optinfo_item *item);

  void emit_for_opt_problem () const;

 private:
  /* Pre-canned ways of manipulating the optinfo, for use by friend class
     dump_context.  */
  void handle_dump_file_kind (dump_flags_t);

 private:
  dump_location_t m_loc;
  enum optinfo_kind m_kind;
  opt_pass *m_pass;
  auto_vec <optinfo_item *> m_items;
};

/* An enum for discriminating between different kinds of optinfo_item.  */

enum optinfo_item_kind
{
  OPTINFO_ITEM_KIND_TEXT,
  OPTINFO_ITEM_KIND_TREE,
  OPTINFO_ITEM_KIND_GIMPLE,
  OPTINFO_ITEM_KIND_SYMTAB_NODE
};

/* An item within an optinfo.  */

class optinfo_item
{
 public:
  optinfo_item (enum optinfo_item_kind kind, location_t location,
		char *text);
  ~optinfo_item ();

  enum optinfo_item_kind get_kind () const { return m_kind; }
  location_t get_location () const { return m_location; }
  const char *get_text () const { return m_text; }

 private:
  /* Metadata (e.g. for optimization records).  */
  enum optinfo_item_kind m_kind;
  location_t m_location;

  /* The textual form of the item, owned by the item.  */
  char *m_text;
};

#endif /* #ifndef GCC_OPTINFO_H */
