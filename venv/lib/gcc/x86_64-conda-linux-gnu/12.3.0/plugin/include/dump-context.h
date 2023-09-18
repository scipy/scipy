/* Support code for handling the various dump_* calls in dumpfile.h
   Copyright (C) 2018-2022 Free Software Foundation, Inc.
   Contributed by David Malcolm <dmalcolm@redhat.com>.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3, or (at your option)
any later version.

GCC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GCC; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.  */


#ifndef GCC_DUMP_CONTEXT_H
#define GCC_DUMP_CONTEXT_H 1

#include "dumpfile.h"
#include "pretty-print.h"
#include "selftest.h"
#include "optinfo.h"

class optrecord_json_writer;
namespace selftest { class temp_dump_context; }
class debug_dump_context;

/* A class for handling the various dump_* calls.

   In particular, this class has responsibility for consolidating
   the "dump_*" calls into optinfo instances (delimited by "dump_*_loc"
   calls), and emitting them.

   Putting this in a class (rather than as global state) allows
   for selftesting of this code.  */

class dump_context
{
  friend class selftest::temp_dump_context;
  friend class debug_dump_context;

 public:
  static dump_context &get () { return *s_current; }

  ~dump_context ();

  void refresh_dumps_are_enabled ();

  void dump_loc (const dump_metadata_t &metadata,
		 const dump_user_location_t &loc);
  void dump_loc_immediate (dump_flags_t dump_kind,
			   const dump_user_location_t &loc);

  void dump_gimple_stmt (const dump_metadata_t &metadata,
			 dump_flags_t extra_dump_flags,
			 gimple *gs, int spc);

  void dump_gimple_stmt_loc (const dump_metadata_t &metadata,
			     const dump_user_location_t &loc,
			     dump_flags_t extra_dump_flags,
			     gimple *gs, int spc);

  void dump_gimple_expr (const dump_metadata_t &metadata,
			 dump_flags_t extra_dump_flags,
			 gimple *gs, int spc);

  void dump_gimple_expr_loc (const dump_metadata_t &metadata,
			     const dump_user_location_t &loc,
			     dump_flags_t extra_dump_flags,
			     gimple *gs,
			     int spc);

  void dump_generic_expr (const dump_metadata_t &metadata,
			  dump_flags_t extra_dump_flags,
			  tree t);

  void dump_generic_expr_loc (const dump_metadata_t &metadata,
			      const dump_user_location_t &loc,
			      dump_flags_t extra_dump_flags,
			      tree t);

  void dump_printf_va (const dump_metadata_t &metadata, const char *format,
		       va_list *ap) ATTRIBUTE_GCC_DUMP_PRINTF (3, 0);

  void dump_printf_loc_va (const dump_metadata_t &metadata,
			   const dump_user_location_t &loc,
			   const char *format, va_list *ap)
    ATTRIBUTE_GCC_DUMP_PRINTF (4, 0);

  template<unsigned int N, typename C>
  void dump_dec (const dump_metadata_t &metadata, const poly_int<N, C> &value);

  void dump_symtab_node (const dump_metadata_t &metadata, symtab_node *node);

  /* Managing nested scopes.  */
  unsigned int get_scope_depth () const;
  void begin_scope (const char *name,
		    const dump_user_location_t &user_location,
		    const dump_impl_location_t &impl_location);
  void end_scope ();

  /* Should optinfo instances be created?
     All creation of optinfos should be guarded by this predicate.
     Return true if any optinfo destinations are active.  */
  bool optinfo_enabled_p () const;

  bool optimization_records_enabled_p () const
  {
    return m_json_writer != NULL;
  }
  void set_json_writer (optrecord_json_writer *writer);
  void finish_any_json_writer ();

  void end_any_optinfo ();

  void emit_optinfo (const optinfo *info);
  void emit_item (optinfo_item *item, dump_flags_t dump_kind);

  bool apply_dump_filter_p (dump_flags_t dump_kind, dump_flags_t filter) const;

 private:
  optinfo &ensure_pending_optinfo (const dump_metadata_t &metadata);
  optinfo &begin_next_optinfo (const dump_metadata_t &metadata,
			       const dump_user_location_t &loc);

  /* The current nesting depth of dump scopes, for showing nesting
     via indentation).  */
  unsigned int m_scope_depth;

  /* The optinfo currently being accumulated since the last dump_*_loc call,
     if any.  */
  optinfo *m_pending;

  /* If -fsave-optimization-record is enabled, the heap-allocated JSON writer
     instance, otherwise NULL.  */
  optrecord_json_writer *m_json_writer;

  /* For use in selftests: if non-NULL, then items are to be printed
     to this, using the given flags.  */
  pretty_printer *m_test_pp;
  dump_flags_t m_test_pp_flags;

  /* The currently active dump_context, for use by the dump_* API calls.  */
  static dump_context *s_current;

  /* The default active context.  */
  static dump_context s_default;
};

/* A subclass of pretty_printer for implementing dump_context::dump_printf_va.
   In particular, the formatted chunks are captured as optinfo_item instances,
   thus retaining metadata about the entities being dumped (e.g. source
   locations), rather than just as plain text.  */

class dump_pretty_printer : public pretty_printer
{
public:
  dump_pretty_printer (dump_context *context, dump_flags_t dump_kind);

  void emit_items (optinfo *dest);

private:
  /* Information on an optinfo_item that was generated during phase 2 of
     formatting.  */
  class stashed_item
  {
  public:
    stashed_item (const char **buffer_ptr_, optinfo_item *item_)
      : buffer_ptr (buffer_ptr_), item (item_) {}
    const char **buffer_ptr;
    optinfo_item *item;
  };

  static bool format_decoder_cb (pretty_printer *pp, text_info *text,
				 const char *spec, int /*precision*/,
				 bool /*wide*/, bool /*set_locus*/,
				 bool /*verbose*/, bool */*quoted*/,
				 const char **buffer_ptr);

  bool decode_format (text_info *text, const char *spec,
		      const char **buffer_ptr);

  void stash_item (const char **buffer_ptr, optinfo_item *item);

  void emit_any_pending_textual_chunks (optinfo *dest);

  void emit_item (optinfo_item *item, optinfo *dest);

  dump_context *m_context;
  dump_flags_t m_dump_kind;
  auto_vec<stashed_item> m_stashed_items;
};

/* An RAII-style class for use in debug dumpers for temporarily using a
   different dump_context.  It enables full details and outputs to
   stderr instead of the currently active dump_file.  */

class debug_dump_context
{
 public:
  debug_dump_context (FILE *f = stderr);
  ~debug_dump_context ();

 private:
  dump_context m_context;
  dump_context *m_saved;
  dump_flags_t m_saved_flags;
  dump_flags_t m_saved_pflags;
  FILE *m_saved_file;
};


#if CHECKING_P

namespace selftest {

/* An RAII-style class for use in selftests for temporarily using a different
   dump_context.  */

class temp_dump_context
{
 public:
  temp_dump_context (bool forcibly_enable_optinfo,
		     bool forcibly_enable_dumping,
		     dump_flags_t test_pp_flags);
  ~temp_dump_context ();

  /* Support for selftests.  */
  optinfo *get_pending_optinfo () const { return m_context.m_pending; }
  const char *get_dumped_text ();

 private:
  pretty_printer m_pp;
  dump_context m_context;
  dump_context *m_saved;
};

/* Implementation detail of ASSERT_DUMPED_TEXT_EQ.  */

extern void verify_dumped_text (const location &loc,
				temp_dump_context *context,
				const char *expected_text);

/* Verify that the text dumped so far in CONTEXT equals
   EXPECTED_TEXT.
   As a side-effect, the internal buffer is 0-terminated.  */

#define ASSERT_DUMPED_TEXT_EQ(CONTEXT, EXPECTED_TEXT)			\
  SELFTEST_BEGIN_STMT							\
    verify_dumped_text (SELFTEST_LOCATION, &(CONTEXT), (EXPECTED_TEXT)); \
  SELFTEST_END_STMT


/* Verify that ITEM has the expected values.  */

void
verify_item (const location &loc,
	     const optinfo_item *item,
	     enum optinfo_item_kind expected_kind,
	     location_t expected_location,
	     const char *expected_text);

/* Verify that ITEM is a text item, with EXPECTED_TEXT.  */

#define ASSERT_IS_TEXT(ITEM, EXPECTED_TEXT) \
  SELFTEST_BEGIN_STMT						    \
    verify_item (SELFTEST_LOCATION, (ITEM), OPTINFO_ITEM_KIND_TEXT, \
		 UNKNOWN_LOCATION, (EXPECTED_TEXT));		    \
  SELFTEST_END_STMT

/* Verify that ITEM is a tree item, with the expected values.  */

#define ASSERT_IS_TREE(ITEM, EXPECTED_LOCATION, EXPECTED_TEXT) \
  SELFTEST_BEGIN_STMT						    \
    verify_item (SELFTEST_LOCATION, (ITEM), OPTINFO_ITEM_KIND_TREE, \
		 (EXPECTED_LOCATION), (EXPECTED_TEXT));	    \
  SELFTEST_END_STMT

/* Verify that ITEM is a gimple item, with the expected values.  */

#define ASSERT_IS_GIMPLE(ITEM, EXPECTED_LOCATION, EXPECTED_TEXT) \
  SELFTEST_BEGIN_STMT						    \
    verify_item (SELFTEST_LOCATION, (ITEM), OPTINFO_ITEM_KIND_GIMPLE, \
		 (EXPECTED_LOCATION), (EXPECTED_TEXT));	    \
  SELFTEST_END_STMT

/* Verify that ITEM is a symtab node, with the expected values.  */

#define ASSERT_IS_SYMTAB_NODE(ITEM, EXPECTED_LOCATION, EXPECTED_TEXT) \
  SELFTEST_BEGIN_STMT						    \
    verify_item (SELFTEST_LOCATION, (ITEM), OPTINFO_ITEM_KIND_SYMTAB_NODE, \
		 (EXPECTED_LOCATION), (EXPECTED_TEXT));	    \
  SELFTEST_END_STMT

} // namespace selftest

#endif /* CHECKING_P */

#endif /* GCC_DUMP_CONTEXT_H */
