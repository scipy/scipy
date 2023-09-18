/* Data structures and declarations used for reading and writing
   GIMPLE to a file stream.

   Copyright (C) 2009-2022 Free Software Foundation, Inc.
   Contributed by Doug Kwan <dougkwan@google.com>

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

#ifndef GCC_LTO_STREAMER_H
#define GCC_LTO_STREAMER_H

#include "plugin-api.h"
#include "gcov-io.h"
#include "diagnostic.h"
#include "version.h"

/* The encoding for a function consists of the following sections:

   1)    The header.
   2)    FIELD_DECLS.
   3)    FUNCTION_DECLS.
   4)    global VAR_DECLS.
   5)    type_decls
   6)    types.
   7)    Names for the labels that have names
   8)    The SSA names.
   9)    The control flow graph.
   10-11)Gimple for local decls.
   12)   Gimple for the function.
   13)   Strings.

   1) THE HEADER.
   2-6) THE GLOBAL DECLS AND TYPES.

      The global decls and types are encoded in the same way.  For each
      entry, there is word with the offset within the section to the
      entry.

   7) THE LABEL NAMES.

      Since most labels do not have names, this section my be of zero
      length.  It consists of an array of string table references, one
      per label.  In the lto code, the labels are given either
      positive or negative indexes.  the positive ones have names and
      the negative ones do not.  The positive index can be used to
      find the name in this array.

   9) THE CFG.

   10) Index into the local decls.  Since local decls can have local
      decls inside them, they must be read in randomly in order to
      properly restore them.

   11-12) GIMPLE FOR THE LOCAL DECLS AND THE FUNCTION BODY.

     The gimple consists of a set of records.

     THE FUNCTION

     At the top level of (8) is the function. It consists of five
     pieces:

     LTO_function     - The tag.
     eh tree          - This is all of the exception handling regions
                        put out in a post order traversial of the
                        tree.  Siblings are output as lists terminated
			by a 0.  The set of fields matches the fields
			defined in except.cc.

     last_basic_block - in uleb128 form.

     basic blocks     - This is the set of basic blocks.

     zero             - The termination of the basic blocks.

     BASIC BLOCKS

     There are two forms of basic blocks depending on if they are
     empty or not.

     The basic block consists of:

     LTO_bb1 or LTO_bb0 - The tag.

     bb->index          - the index in uleb128 form.

     #succs             - The number of successors un uleb128 form.

     the successors     - For each edge, a pair.  The first of the
                          pair is the index of the successor in
                          uleb128 form and the second are the flags in
                          uleb128 form.

     the statements     - A gimple tree, as described above.
                          These are only present for LTO_BB1.
                          Following each statement is an optional
                          exception handling record LTO_eh_region
			  which contains the region number (for
			  regions >= 0).

     zero               - This is only present for LTO_BB1 and is used
			  to terminate the statements and exception
			  regions within this block.

   12) STRINGS

     String are represented in the table as pairs, a length in ULEB128
     form followed by the data for the string.  */

#define LTO_major_version GCC_major_version
#define LTO_minor_version 0

typedef unsigned char	lto_decl_flags_t;

/* Stream additional data to LTO object files to make it easier to debug
   streaming code.  This changes object files.  */
static const bool streamer_debugging = false;

/* Tags representing the various IL objects written to the bytecode file
   (GIMPLE statements, basic blocks, EH regions, tree nodes, etc).

   NOTE, when adding new LTO tags, also update lto_tag_name.  */
enum LTO_tags
{
  LTO_null = 0,

  /* Reference to previously-streamed node.  */
  LTO_tree_pickle_reference,

  /* References to indexable tree nodes.  These objects are stored in
     tables that are written separately from the function bodies
     and variable constructors that reference them.  This way they can be
     instantiated even when the referencing functions aren't (e.g., during WPA)
     and it also allows functions to be copied from one file to another without
     having to unpickle the body first (the references are location
     independent).  */
  LTO_global_stream_ref,

  LTO_ssa_name_ref,

  /* Special for global streamer.  A blob of unnamed tree nodes.  */
  LTO_tree_scc,

  /* Sequence of trees.  */
  LTO_trees,

  /* Shared INTEGER_CST node.  */
  LTO_integer_cst,

  /* Tags of trees are encoded as
     LTO_first_tree_tag + TREE_CODE.  */
  LTO_first_tree_tag,
  /* Tags of gimple typles are encoded as
     LTO_first_gimple_tag + gimple_code.  */
  LTO_first_gimple_tag = LTO_first_tree_tag + MAX_TREE_CODES,

  /* Entry and exit basic blocks.  */
  LTO_bb0 = LTO_first_gimple_tag + LAST_AND_UNUSED_GIMPLE_CODE,
  LTO_bb1,

  /* EH region holding the previous statement.  */
  LTO_eh_region,

  /* Function body.  */
  LTO_function,

  /* EH table.  */
  LTO_eh_table,

  /* EH region types.  These mirror enum eh_region_type.  */
  LTO_ert_cleanup,
  LTO_ert_try,
  LTO_ert_allowed_exceptions,
  LTO_ert_must_not_throw,

  /* EH landing pad.  */
  LTO_eh_landing_pad,

  /* EH try/catch node.  */
  LTO_eh_catch,

  /* This tag must always be last.  */
  LTO_NUM_TAGS
};


/* Set of section types that are in an LTO file.  This list will grow
   as the number of IPA passes grows since each IPA pass will need its
   own section type to store its summary information.

   When adding a new section type, you must also extend the
   LTO_SECTION_NAME array in lto-section-in.cc.  */
enum lto_section_type
{
  LTO_section_decls = 0,
  LTO_section_function_body,
  LTO_section_static_initializer,
  LTO_section_symtab,
  LTO_section_symtab_extension,
  LTO_section_refs,
  LTO_section_asm,
  LTO_section_jump_functions,
  LTO_section_ipa_pure_const,
  LTO_section_ipa_reference,
  LTO_section_ipa_profile,
  LTO_section_symtab_nodes,
  LTO_section_opts,
  LTO_section_cgraph_opt_sum,
  LTO_section_ipa_fn_summary,
  LTO_section_ipcp_transform,
  LTO_section_ipa_icf,
  LTO_section_offload_table,
  LTO_section_mode_table,
  LTO_section_lto,
  LTO_section_ipa_sra,
  LTO_section_odr_types,
  LTO_section_ipa_modref,
  LTO_N_SECTION_TYPES		/* Must be last.  */
};

/* Indices to the various function, type and symbol streams. */
enum lto_decl_stream_e_t
{
  LTO_DECL_STREAM = 0,		/* Must be first.  */
  LTO_N_DECL_STREAMS
};

typedef enum ld_plugin_symbol_resolution ld_plugin_symbol_resolution_t;

/* Return a char pointer to the start of a data stream for an lto pass
   or function.  The first parameter is the file data that contains
   the information.  The second parameter is the type of information
   to be obtained.  The third parameter is the name of the function
   and is only used when finding a function body; otherwise it is
   NULL.  The fourth parameter is the length of the data returned.  */
typedef const char* (lto_get_section_data_f) (struct lto_file_decl_data *,
					      enum lto_section_type,
					      const char *,
					      int,
					      size_t *);

/* Return the data found from the above call.  The first three
   parameters are the same as above.  The fourth parameter is the data
   itself and the fifth is the length of the data. */
typedef void (lto_free_section_data_f) (struct lto_file_decl_data *,
					enum lto_section_type,
					const char *,
					const char *,
					size_t);

/* The location cache holds expanded locations for streamed in trees.
   This is done to reduce memory usage of libcpp linemap that strongly prefers
   locations to be inserted in the source order.  */

class lto_location_cache
{
public:
  /* Apply all changes in location cache.  Add locations into linemap and patch
     trees.  */
  bool apply_location_cache ();
  /* Tree merging did not suceed; mark all changes in the cache as accepted.  */
  void accept_location_cache ();
  /* Tree merging did suceed; throw away recent changes.  */
  void revert_location_cache ();
  void input_location (location_t *loc, struct bitpack_d *bp,
		       class data_in *data_in);
  void input_location_and_block (location_t *loc, struct bitpack_d *bp,
				 class lto_input_block *ib,
				 class data_in *data_in);
  lto_location_cache ()
     : loc_cache (), accepted_length (0), current_file (NULL), current_line (0),
       current_col (0), current_sysp (false), current_loc (UNKNOWN_LOCATION),
       current_block (NULL_TREE)
  {
    gcc_assert (!current_cache);
    current_cache = this;
  }
  ~lto_location_cache ()
  {
    apply_location_cache ();
    gcc_assert (current_cache == this);
    current_cache = NULL;
  }

  /* There can be at most one instance of location cache (combining multiple
     would bring it out of sync with libcpp linemap); point to current
     one.  */
  static lto_location_cache *current_cache;
  
private:
  static int cmp_loc (const void *pa, const void *pb);

  struct cached_location
  {
    const char *file;
    location_t *loc;
    int line, col;
    bool sysp;
    tree block;
  };

  /* The location cache.  */

  auto_vec<cached_location> loc_cache;

  /* Accepted entries are ones used by trees that are known to be not unified
     by tree merging.  */

  int accepted_length;

  /* Bookkeeping to remember state in between calls to lto_apply_location_cache
     When streaming gimple, the location cache is not used and thus
     lto_apply_location_cache happens per location basis.  It is then
     useful to avoid redundant calls of linemap API.  */

  const char *current_file;
  int current_line;
  int current_col;
  bool current_sysp;
  location_t current_loc;
  tree current_block;
};

/* Structure used as buffer for reading an LTO file.  */
class lto_input_block
{
public:
  /* Special constructor for the string table, it abuses this to
     do random access but use the uhwi decoder.  */
  lto_input_block (const char *data_, unsigned int p_, unsigned int len_,
		   const unsigned char *mode_table_)
      : data (data_), mode_table (mode_table_), p (p_), len (len_) {}
  lto_input_block (const char *data_, unsigned int len_,
		   const unsigned char *mode_table_)
      : data (data_), mode_table (mode_table_), p (0), len (len_) {}

  const char *data;
  const unsigned char *mode_table;
  unsigned int p;
  unsigned int len;
};

/* Compression algorithm used for compression of LTO bytecode.  */

enum lto_compression
{
  ZLIB,
  ZSTD
};

/* Structure that represents LTO ELF section with information
   about the format.  */

struct lto_section
{
  int16_t major_version;
  int16_t minor_version;
  unsigned char slim_object;
  unsigned char _padding;

  /* Flags is a private field that is not defined publicly.  */
  uint16_t flags;

  /* Set compression to FLAGS.  */
  inline void set_compression (lto_compression c)
  {
    flags = c;
  }

  /* Get compression from FLAGS.  */
  inline lto_compression get_compression ()
  {
    return (lto_compression) flags;
  }
};

STATIC_ASSERT (sizeof (lto_section) == 8);

/* The is the first part of the record in an LTO file for many of the
   IPA passes.  */
struct lto_simple_header
{
  /* Size of main gimple body of function.  */
  int32_t main_size;
};

struct lto_simple_header_with_strings : lto_simple_header
{
  /* Size of the string table.  */
  int32_t string_size;
};

/* The header for a function body.  */
struct lto_function_header : lto_simple_header_with_strings
{
  /* Size of the cfg.  */
  int32_t cfg_size;
};


/* Structure describing a symbol section.  */
struct lto_decl_header : lto_simple_header_with_strings
{
  /* Size of region for decl state. */
  int32_t decl_state_size;

  /* Number of nodes in globals stream.  */
  int32_t num_nodes;
};


/* Statistics gathered during LTO, WPA and LTRANS.  */
struct lto_stats_d
{
  unsigned HOST_WIDE_INT num_input_cgraph_nodes;
  unsigned HOST_WIDE_INT num_output_symtab_nodes;
  unsigned HOST_WIDE_INT num_input_files;
  unsigned HOST_WIDE_INT num_output_files;
  unsigned HOST_WIDE_INT num_cgraph_partitions;
  unsigned HOST_WIDE_INT section_size[LTO_N_SECTION_TYPES];
  unsigned HOST_WIDE_INT num_function_bodies;
  unsigned HOST_WIDE_INT num_trees[NUM_TREE_CODES];
  unsigned HOST_WIDE_INT num_output_il_bytes;
  unsigned HOST_WIDE_INT num_compressed_il_bytes;
  unsigned HOST_WIDE_INT num_input_il_bytes;
  unsigned HOST_WIDE_INT num_uncompressed_il_bytes;
  unsigned HOST_WIDE_INT num_tree_bodies_output;
  unsigned HOST_WIDE_INT num_pickle_refs_output;
};

/* Entry of LTO symtab encoder.  */
struct lto_encoder_entry
{
  symtab_node *node;
  /* Is the node in this partition (i.e. ltrans of this partition will
     be responsible for outputting it)? */
  unsigned int in_partition:1;
  /* Do we encode body in this partition?  */
  unsigned int body:1;
  /* Do we encode initializer in this partition?
     For example the readonly variable initializers are encoded to aid
     constant folding even if they are not in the partition.  */
  unsigned int initializer:1;
};


/* Encoder data structure used to stream callgraph nodes.  */
struct lto_symtab_encoder_d
{
  vec<lto_encoder_entry> nodes;
  hash_map<symtab_node *, size_t> *map;
};

typedef struct lto_symtab_encoder_d *lto_symtab_encoder_t;

/* Iterator structure for cgraph node sets.  */
struct lto_symtab_encoder_iterator
{
  lto_symtab_encoder_t encoder;
  unsigned index;
};



/* The lto_tree_ref_encoder struct is used to encode trees into indices. */

struct lto_tree_ref_encoder
{
  hash_map<tree, unsigned> *tree_hash_table;	/* Maps pointers to indices. */
  vec<tree> trees;			/* Maps indices to pointers. */
};


/* Structure to hold states of input scope.  */
struct GTY((for_user)) lto_in_decl_state
{
  /* Array of lto_in_decl_buffers to store type and decls streams. */
  vec<tree, va_gc> *streams[LTO_N_DECL_STREAMS];

  /* If this in-decl state is associated with a function. FN_DECL
     point to the FUNCTION_DECL. */
  tree fn_decl;

  /* True if decl state is compressed.  */
  bool compressed;
};

typedef struct lto_in_decl_state *lto_in_decl_state_ptr;

struct decl_state_hasher : ggc_ptr_hash<lto_in_decl_state>
{
  static hashval_t
  hash (lto_in_decl_state *s)
  {
    return htab_hash_pointer (s->fn_decl);
  }

  static bool
  equal (lto_in_decl_state *a, lto_in_decl_state *b)
  {
    return a->fn_decl == b->fn_decl;
  }
};

/* The structure that holds all of the vectors of global types,
   decls and cgraph nodes used in the serialization of this file.  */
struct lto_out_decl_state
{
  /* The buffers contain the sets of decls of various kinds and types we have
     seen so far and the indexes assigned to them.  */
  struct lto_tree_ref_encoder streams[LTO_N_DECL_STREAMS];

  /* Encoder for cgraph nodes.  */
  lto_symtab_encoder_t symtab_node_encoder;

  /* If this out-decl state belongs to a function, fn_decl points to that
     function.  Otherwise, it is NULL. */
  tree fn_decl;

  /* True if decl state is compressed.  */
  bool compressed;
};

typedef struct lto_out_decl_state *lto_out_decl_state_ptr;


/* Compact representation of a index <-> resolution pair. Unpacked to an 
   vector later. */
struct res_pair 
{
  ld_plugin_symbol_resolution_t res;
  unsigned index;
};


/* One of these is allocated for each object file that being compiled
   by lto.  This structure contains the tables that are needed by the
   serialized functions and ipa passes to connect themselves to the
   global types and decls as they are reconstituted.  */
struct GTY(()) lto_file_decl_data
{
  /* Decl state currently used. */
  struct lto_in_decl_state *current_decl_state;

  /* Decl state corresponding to regions outside of any functions
     in the compilation unit. */
  struct lto_in_decl_state *global_decl_state;

  /* Table of cgraph nodes present in this file.  */
  lto_symtab_encoder_t GTY((skip)) symtab_node_encoder;

  /* Hash table maps lto-related section names to location in file.  */
  hash_table<decl_state_hasher> *function_decl_states;

  /* The .o file that these offsets relate to.  */
  const char *GTY((skip)) file_name;

  /* Hash table maps lto-related section names to location in file.  */
  htab_t GTY((skip)) section_hash_table;

  /* Hash new name of renamed global declaration to its original name.  */
  htab_t GTY((skip)) renaming_hash_table;

  /* Linked list used temporarily in reader */
  struct lto_file_decl_data *next;

  /* Order in which the file appears on the command line.  */
  int order;

  /* Sub ID for merged objects. */
  unsigned HOST_WIDE_INT id;

  /* Symbol resolutions for this file */
  vec<res_pair>  GTY((skip)) respairs;
  unsigned max_index;

  gcov_summary GTY((skip)) profile_info;

  /* Map assigning declarations their resolutions.  */
  hash_map<tree, ld_plugin_symbol_resolution> * GTY((skip)) resolution_map;

  /* Mode translation table.  */
  const unsigned char *mode_table;

  /* Read LTO section.  */
  lto_section lto_section_header;

  int order_base;

  int unit_base;
};

typedef struct lto_file_decl_data *lto_file_decl_data_ptr;

struct lto_char_ptr_base
{
  char *ptr;
};

/* An incore byte stream to buffer the various parts of the function.
   The entire structure should be zeroed when created.  The record
   consists of a set of blocks.  The first sizeof (ptr) bytes are used
   as a chain, and the rest store the bytes to be written.  */
struct lto_output_stream
{
  /* The pointer to the first block in the stream.  */
  struct lto_char_ptr_base * first_block;

  /* The pointer to the last and current block in the stream.  */
  struct lto_char_ptr_base * current_block;

  /* The pointer to where the next char should be written.  */
  char * current_pointer;

  /* The number of characters left in the current block.  */
  unsigned int left_in_block;

  /* The block size of the last block allocated.  */
  unsigned int block_size;

  /* The total number of characters written.  */
  unsigned int total_size;
};

/* A simple output block.  This can be used for simple IPA passes that
   do not need more than one stream.  */
struct lto_simple_output_block
{
  enum lto_section_type section_type;
  struct lto_out_decl_state *decl_state;

  /* The stream that the main tree codes are written to.  */
  struct lto_output_stream *main_stream;
};

/* String hashing.  */

struct string_slot
{
  const char *s;
  int len;
  unsigned int slot_num;
};

/* Hashtable helpers.  */

struct string_slot_hasher : nofree_ptr_hash <string_slot>
{
  static inline hashval_t hash (const string_slot *);
  static inline bool equal (const string_slot *, const string_slot *);
};

/* Returns a hash code for DS.  Adapted from libiberty's htab_hash_string
   to support strings that may not end in '\0'.  */

inline hashval_t
string_slot_hasher::hash (const string_slot *ds)
{
  hashval_t r = ds->len;
  int i;

  for (i = 0; i < ds->len; i++)
     r = r * 67 + (unsigned)ds->s[i] - 113;
  return r;
}

/* Returns nonzero if DS1 and DS2 are equal.  */

inline bool
string_slot_hasher::equal (const string_slot *ds1, const string_slot *ds2)
{
  if (ds1->len == ds2->len)
    return memcmp (ds1->s, ds2->s, ds1->len) == 0;

  return 0;
}

/* Data structure holding all the data and descriptors used when writing
   an LTO file.  */
struct output_block
{
  enum lto_section_type section_type;
  struct lto_out_decl_state *decl_state;

  /* The stream that the main tree codes are written to.  */
  struct lto_output_stream *main_stream;

  /* The stream that contains the string table.  */
  struct lto_output_stream *string_stream;

  /* The stream that contains the cfg.  */
  struct lto_output_stream *cfg_stream;

  /* The hash table that contains the set of strings we have seen so
     far and the indexes assigned to them.  */
  hash_table<string_slot_hasher> *string_hash_table;

  /* The current symbol that we are currently serializing.  Null
     if we are serializing something else.  */
  symtab_node *symbol;

  /* These are the last file and line that were seen in the stream.
     If the current node differs from these, it needs to insert
     something into the stream and fix these up.  */
  const char *current_file;
  int current_line;
  int current_col;
  bool current_sysp;
  bool reset_locus;
  bool emit_pwd;
  tree current_block;

  /* Cache of nodes written in this section.  */
  struct streamer_tree_cache_d *writer_cache;

  /* All trees identified as local to the unit streamed.  */
  hash_set<tree> *local_trees;

  /* All data persistent across whole duration of output block
     can go here.  */
  struct obstack obstack;
};


/* Data and descriptors used when reading from an LTO file.  */
class data_in
{
public:
  /* The global decls and types.  */
  struct lto_file_decl_data *file_data;

  /* The string table.  */
  const char *strings;

  /* The length of the string table.  */
  unsigned int strings_len;

  /* Maps each reference number to the resolution done by the linker. */
  vec<ld_plugin_symbol_resolution_t> globals_resolution;

  /* Cache of pickled nodes.  */
  struct streamer_tree_cache_d *reader_cache;

  /* Cache of source code location.  */
  lto_location_cache location_cache;
};


/* In lto-section-in.cc  */
extern class lto_input_block * lto_create_simple_input_block (
			       struct lto_file_decl_data *,
			       enum lto_section_type, const char **, size_t *);
extern void
lto_destroy_simple_input_block (struct lto_file_decl_data *,
				enum lto_section_type,
				class lto_input_block *, const char *, size_t);
extern void lto_set_in_hooks (struct lto_file_decl_data **,
			      lto_get_section_data_f *,
			      lto_free_section_data_f *);
extern struct lto_file_decl_data **lto_get_file_decl_data (void);
extern const char *lto_get_section_data (struct lto_file_decl_data *,
					 enum lto_section_type,
					 const char *, int, size_t *,
					 bool decompress = false);
extern const char *lto_get_summary_section_data (struct lto_file_decl_data *,
						 enum lto_section_type,
						 size_t *);
extern const char *lto_get_raw_section_data (struct lto_file_decl_data *,
					     enum lto_section_type,
					     const char *, int, size_t *);
extern void lto_free_section_data (struct lto_file_decl_data *,
			           enum lto_section_type,
				   const char *, const char *, size_t,
				   bool decompress = false);
extern void lto_free_raw_section_data (struct lto_file_decl_data *,
				       enum lto_section_type,
				       const char *, const char *, size_t);
extern htab_t lto_create_renaming_table (void);
extern void lto_record_renamed_decl (struct lto_file_decl_data *,
				     const char *, const char *);
extern const char *lto_get_decl_name_mapping (struct lto_file_decl_data *,
					      const char *);
extern struct lto_in_decl_state *lto_new_in_decl_state (void);
extern void lto_delete_in_decl_state (struct lto_in_decl_state *);
extern struct lto_in_decl_state *lto_get_function_in_decl_state (
				      struct lto_file_decl_data *, tree);
extern void lto_free_function_in_decl_state (struct lto_in_decl_state *);
extern void lto_free_function_in_decl_state_for_node (symtab_node *);
extern void lto_section_overrun (class lto_input_block *) ATTRIBUTE_NORETURN;
extern void lto_value_range_error (const char *,
				   HOST_WIDE_INT, HOST_WIDE_INT,
				   HOST_WIDE_INT) ATTRIBUTE_NORETURN;

/* In lto-section-out.cc  */
extern void lto_begin_section (const char *, bool);
extern void lto_end_section (void);
extern void lto_write_data (const void *, unsigned int);
extern void lto_write_raw_data (const void *, unsigned int);
extern void lto_write_stream (struct lto_output_stream *);
extern struct lto_simple_output_block *lto_create_simple_output_block (
				enum lto_section_type);
extern void lto_destroy_simple_output_block (struct lto_simple_output_block *);
extern struct lto_out_decl_state *lto_new_out_decl_state (void);
extern void lto_delete_out_decl_state (struct lto_out_decl_state *);
extern struct lto_out_decl_state *lto_get_out_decl_state (void);
extern void lto_push_out_decl_state (struct lto_out_decl_state *);
extern struct lto_out_decl_state *lto_pop_out_decl_state (void);
extern void lto_record_function_out_decl_state (tree,
						struct lto_out_decl_state *);
extern void lto_append_block (struct lto_output_stream *);


/* In lto-streamer.cc.  */

/* Set when streaming LTO for offloading compiler.  */
extern bool lto_stream_offload_p;

extern const char *lto_tag_name (enum LTO_tags);
extern char *lto_get_section_name (int, const char *, int,
				   struct lto_file_decl_data *);
extern void print_lto_report (const char *);
extern void lto_streamer_init (void);
extern bool gate_lto_out (void);
extern void lto_check_version (int, int, const char *);
extern void lto_streamer_hooks_init (void);

/* In lto-streamer-in.cc */
extern void lto_input_cgraph (struct lto_file_decl_data *, const char *);
extern void lto_reader_init (void);
extern void lto_free_file_name_hash (void);
extern void lto_input_function_body (struct lto_file_decl_data *,
				     struct cgraph_node *,
				     const char *);
extern void lto_input_variable_constructor (struct lto_file_decl_data *,
					    struct varpool_node *,
					    const char *);
extern void lto_input_constructors_and_inits (struct lto_file_decl_data *,
					      const char *);
extern void lto_input_toplevel_asms (struct lto_file_decl_data *, int);
extern void lto_input_mode_table (struct lto_file_decl_data *);
extern class data_in *lto_data_in_create (struct lto_file_decl_data *,
				    const char *, unsigned,
				    vec<ld_plugin_symbol_resolution_t> );
extern void lto_data_in_delete (class data_in *);
extern void lto_input_data_block (class lto_input_block *, void *, size_t);
void lto_input_location (location_t *, struct bitpack_d *, class data_in *);
tree lto_input_tree_ref (class lto_input_block *, class data_in *,
			 struct function *, enum LTO_tags);
void lto_tag_check_set (enum LTO_tags, int, ...);
void lto_init_eh (void);
hashval_t lto_input_scc (class lto_input_block *, class data_in *,
			 unsigned *, unsigned *, bool);
tree lto_input_tree_1 (class lto_input_block *, class data_in *,
		       enum LTO_tags, hashval_t hash);
tree lto_input_tree (class lto_input_block *, class data_in *);
tree stream_read_tree_ref (class lto_input_block *, class data_in *);


/* In lto-streamer-out.cc  */
extern void lto_register_decl_definition (tree, struct lto_file_decl_data *);
extern struct output_block *create_output_block (enum lto_section_type);
extern void destroy_output_block (struct output_block *);
extern void lto_output_tree (struct output_block *, tree, bool, bool);
extern void stream_write_tree_ref (struct output_block *, tree);
extern void lto_output_var_decl_ref (struct lto_out_decl_state *,
				     struct lto_output_stream *, tree);
extern void lto_output_fn_decl_ref (struct lto_out_decl_state *,
				    struct lto_output_stream *, tree);
extern tree lto_input_var_decl_ref (lto_input_block *, lto_file_decl_data *);
extern tree lto_input_fn_decl_ref (lto_input_block *, lto_file_decl_data *);
extern void lto_output_toplevel_asms (void);
extern void produce_asm (struct output_block *ob, tree fn);
extern void lto_output ();
extern void produce_asm_for_decls ();
void lto_output_decl_state_streams (struct output_block *,
				    struct lto_out_decl_state *);
void lto_output_decl_state_refs (struct output_block *,
			         struct lto_output_stream *,
			         struct lto_out_decl_state *);
void lto_output_location (struct output_block *, struct bitpack_d *,
			  location_t);
void lto_output_location_and_block (struct output_block *, struct bitpack_d *,
				    location_t);
void lto_output_init_mode_table (void);
void lto_prepare_function_for_streaming (cgraph_node *);


/* In lto-cgraph.cc  */
extern bool asm_nodes_output;
lto_symtab_encoder_t lto_symtab_encoder_new (bool);
int lto_symtab_encoder_encode (lto_symtab_encoder_t, symtab_node *);
void lto_symtab_encoder_delete (lto_symtab_encoder_t);
bool lto_symtab_encoder_delete_node (lto_symtab_encoder_t, symtab_node *);
bool lto_symtab_encoder_encode_body_p (lto_symtab_encoder_t,
				       struct cgraph_node *);
bool lto_symtab_encoder_in_partition_p (lto_symtab_encoder_t,
					symtab_node *);
void lto_set_symtab_encoder_in_partition (lto_symtab_encoder_t,
					  symtab_node *);

bool lto_symtab_encoder_encode_initializer_p (lto_symtab_encoder_t,
					      varpool_node *);
void output_symtab (void);
void input_symtab (void);
void output_offload_tables (void);
void input_offload_tables (bool);
bool referenced_from_other_partition_p (struct ipa_ref_list *,
				        lto_symtab_encoder_t);
bool reachable_from_other_partition_p (struct cgraph_node *,
				       lto_symtab_encoder_t);
bool referenced_from_this_partition_p (symtab_node *,
					lto_symtab_encoder_t);
bool reachable_from_this_partition_p (struct cgraph_node *,
				      lto_symtab_encoder_t);
lto_symtab_encoder_t compute_ltrans_boundary (lto_symtab_encoder_t encoder);
void select_what_to_stream (void);

/* In omp-general.cc.  */
void omp_lto_output_declare_variant_alt (lto_simple_output_block *,
					 cgraph_node *, lto_symtab_encoder_t);
void omp_lto_input_declare_variant_alt (lto_input_block *, cgraph_node *,
					vec<symtab_node *>);

/* In options-save.cc.  */
void cl_target_option_stream_out (struct output_block *, struct bitpack_d *,
				  struct cl_target_option *);

void cl_target_option_stream_in (class data_in *,
				 struct bitpack_d *,
				 struct cl_target_option *);

void cl_optimization_stream_out (struct output_block *,
				 struct bitpack_d *, struct cl_optimization *);

void cl_optimization_stream_in (class data_in *,
				struct bitpack_d *, struct cl_optimization *);



/* In lto-opts.cc.  */
extern void lto_write_options (void);


/* Statistics gathered during LTO, WPA and LTRANS.  */
extern struct lto_stats_d lto_stats;

/* Section names corresponding to the values of enum lto_section_type.  */
extern const char *lto_section_name[];

/* Holds all the out decl states of functions output so far in the
   current output file.  */
extern vec<lto_out_decl_state_ptr> lto_function_decl_states;

/* Return true if LTO tag TAG corresponds to a tree code.  */
static inline bool
lto_tag_is_tree_code_p (enum LTO_tags tag)
{
  return tag > LTO_first_tree_tag && (unsigned) tag <= MAX_TREE_CODES;
}


/* Return true if LTO tag TAG corresponds to a gimple code.  */
static inline bool
lto_tag_is_gimple_code_p (enum LTO_tags tag)
{
  return (unsigned) tag >= LTO_first_gimple_tag
	 && (unsigned) tag
	    < LTO_first_gimple_tag + LAST_AND_UNUSED_GIMPLE_CODE;
}


/* Return the LTO tag corresponding to gimple code CODE.  See enum
   LTO_tags for details on the conversion.  */
static inline enum LTO_tags
lto_gimple_code_to_tag (enum gimple_code code)
{
  return (enum LTO_tags) ((unsigned) code + LTO_first_gimple_tag);
}


/* Return the GIMPLE code corresponding to TAG.  See enum LTO_tags for
   details on the conversion.  */
static inline enum gimple_code
lto_tag_to_gimple_code (enum LTO_tags tag)
{
  gcc_assert (lto_tag_is_gimple_code_p (tag));
  return (enum gimple_code) ((unsigned) tag - LTO_first_gimple_tag);
}


/* Return the LTO tag corresponding to tree code CODE.  See enum
   LTO_tags for details on the conversion.  */
static inline enum LTO_tags
lto_tree_code_to_tag (enum tree_code code)
{
  return (enum LTO_tags) ((unsigned) code + LTO_first_tree_tag);
}


/* Return the tree code corresponding to TAG.  See enum LTO_tags for
   details on the conversion.  */
static inline enum tree_code
lto_tag_to_tree_code (enum LTO_tags tag)
{
  gcc_assert (lto_tag_is_tree_code_p (tag));
  return (enum tree_code) ((unsigned) tag - LTO_first_tree_tag);
}

/* Check that tag ACTUAL == EXPECTED.  */
static inline void
lto_tag_check (enum LTO_tags actual, enum LTO_tags expected)
{
  if (actual != expected)
    internal_error ("bytecode stream: expected tag %s instead of %s",
		    lto_tag_name (expected), lto_tag_name (actual));
}

/* Check that tag ACTUAL is in the range [TAG1, TAG2].  */
static inline void
lto_tag_check_range (enum LTO_tags actual, enum LTO_tags tag1,
		     enum LTO_tags tag2)
{
  if (actual < tag1 || actual > tag2)
    internal_error ("bytecode stream: tag %s is not in the expected range "
		    "[%s, %s]",
		    lto_tag_name (actual),
		    lto_tag_name (tag1),
		    lto_tag_name (tag2));
}

/* Initialize an lto_out_decl_buffer ENCODER.  */
static inline void
lto_init_tree_ref_encoder (struct lto_tree_ref_encoder *encoder)
{
  encoder->tree_hash_table = new hash_map<tree, unsigned> (251);
  encoder->trees.create (0);
}


/* Destroy an lto_tree_ref_encoder ENCODER by freeing its contents.  The
   memory used by ENCODER is not freed by this function.  */
static inline void
lto_destroy_tree_ref_encoder (struct lto_tree_ref_encoder *encoder)
{
  /* Hash table may be delete already.  */
  delete encoder->tree_hash_table;
  encoder->tree_hash_table = NULL;
  encoder->trees.release ();
}

/* Return the number of trees encoded in ENCODER. */
static inline unsigned int
lto_tree_ref_encoder_size (struct lto_tree_ref_encoder *encoder)
{
  return encoder->trees.length ();
}

/* Return the IDX-th tree in ENCODER. */
static inline tree
lto_tree_ref_encoder_get_tree (struct lto_tree_ref_encoder *encoder,
			       unsigned int idx)
{
  return encoder->trees[idx];
}

/* Return number of encoded nodes in ENCODER.  */
static inline int
lto_symtab_encoder_size (lto_symtab_encoder_t encoder)
{
  return encoder->nodes.length ();
}

/* Value used to represent failure of lto_symtab_encoder_lookup.  */
#define LCC_NOT_FOUND	(-1)

/* Look up NODE in encoder.  Return NODE's reference if it has been encoded
   or LCC_NOT_FOUND if it is not there.  */

static inline int
lto_symtab_encoder_lookup (lto_symtab_encoder_t encoder,
			   symtab_node *node)
{
  size_t *slot = encoder->map->get (node);
  return (slot && *slot ? *(slot) - 1 : LCC_NOT_FOUND);
}

/* Return true if iterator LSE points to nothing.  */
static inline bool
lsei_end_p (lto_symtab_encoder_iterator lsei)
{
  return lsei.index >= (unsigned)lto_symtab_encoder_size (lsei.encoder);
}

/* Advance iterator LSE.  */
static inline void
lsei_next (lto_symtab_encoder_iterator *lsei)
{
  lsei->index++;
}

/* Return the node pointed to by LSI.  */
static inline symtab_node *
lsei_node (lto_symtab_encoder_iterator lsei)
{
  return lsei.encoder->nodes[lsei.index].node;
}

/* Return the node pointed to by LSI.  */
static inline struct cgraph_node *
lsei_cgraph_node (lto_symtab_encoder_iterator lsei)
{
  return dyn_cast<cgraph_node *> (lsei.encoder->nodes[lsei.index].node);
}

/* Return the node pointed to by LSI.  */
static inline varpool_node *
lsei_varpool_node (lto_symtab_encoder_iterator lsei)
{
  return dyn_cast<varpool_node *> (lsei.encoder->nodes[lsei.index].node);
}

/* Return the cgraph node corresponding to REF using ENCODER.  */

static inline symtab_node *
lto_symtab_encoder_deref (lto_symtab_encoder_t encoder, int ref)
{
  if (ref == LCC_NOT_FOUND)
    return NULL;

  return encoder->nodes[ref].node;
}

/* Return an iterator to the first node in LSI.  */
static inline lto_symtab_encoder_iterator
lsei_start (lto_symtab_encoder_t encoder)
{
  lto_symtab_encoder_iterator lsei;

  lsei.encoder = encoder;
  lsei.index = 0;
  return lsei;
}

/* Advance iterator LSE.  */
static inline void
lsei_next_in_partition (lto_symtab_encoder_iterator *lsei)
{
  lsei_next (lsei);
  while (!lsei_end_p (*lsei)
	 && !lto_symtab_encoder_in_partition_p (lsei->encoder, lsei_node (*lsei)))
    lsei_next (lsei);
}

/* Return an iterator to the first node in LSI.  */
static inline lto_symtab_encoder_iterator
lsei_start_in_partition (lto_symtab_encoder_t encoder)
{
  lto_symtab_encoder_iterator lsei = lsei_start (encoder);

  if (lsei_end_p (lsei))
    return lsei;
  if (!lto_symtab_encoder_in_partition_p (encoder, lsei_node (lsei)))
    lsei_next_in_partition (&lsei);

  return lsei;
}

/* Advance iterator LSE.  */
static inline void
lsei_next_function_in_partition (lto_symtab_encoder_iterator *lsei)
{
  lsei_next (lsei);
  while (!lsei_end_p (*lsei)
	 && (!is_a <cgraph_node *> (lsei_node (*lsei))
	     || !lto_symtab_encoder_in_partition_p (lsei->encoder, lsei_node (*lsei))))
    lsei_next (lsei);
}

/* Return an iterator to the first node in LSI.  */
static inline lto_symtab_encoder_iterator
lsei_start_function_in_partition (lto_symtab_encoder_t encoder)
{
  lto_symtab_encoder_iterator lsei = lsei_start (encoder);

  if (lsei_end_p (lsei))
    return lsei;
  if (!is_a <cgraph_node *> (lsei_node (lsei))
      || !lto_symtab_encoder_in_partition_p (encoder, lsei_node (lsei)))
    lsei_next_function_in_partition (&lsei);

  return lsei;
}

/* Advance iterator LSE.  */
static inline void
lsei_next_variable_in_partition (lto_symtab_encoder_iterator *lsei)
{
  lsei_next (lsei);
  while (!lsei_end_p (*lsei)
	 && (!is_a <varpool_node *> (lsei_node (*lsei))
	     || !lto_symtab_encoder_in_partition_p (lsei->encoder, lsei_node (*lsei))))
    lsei_next (lsei);
}

/* Return an iterator to the first node in LSI.  */
static inline lto_symtab_encoder_iterator
lsei_start_variable_in_partition (lto_symtab_encoder_t encoder)
{
  lto_symtab_encoder_iterator lsei = lsei_start (encoder);

  if (lsei_end_p (lsei))
    return lsei;
  if (!is_a <varpool_node *> (lsei_node (lsei))
      || !lto_symtab_encoder_in_partition_p (encoder, lsei_node (lsei)))
    lsei_next_variable_in_partition (&lsei);

  return lsei;
}

/* Entry for the delayed registering of decl -> DIE references.  */
struct dref_entry {
    tree decl;
    const char *sym;
    unsigned HOST_WIDE_INT off;
};

extern vec<dref_entry> dref_queue;

extern FILE *streamer_dump_file;

#endif /* GCC_LTO_STREAMER_H  */
