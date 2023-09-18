/* Declarations and definitions dealing with attribute handling.
   Copyright (C) 2013-2022 Free Software Foundation, Inc.

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

#ifndef GCC_ATTRIBS_H
#define GCC_ATTRIBS_H

extern const struct attribute_spec *lookup_attribute_spec (const_tree);
extern void free_attr_data ();
extern void init_attributes (void);

/* Process the attributes listed in ATTRIBUTES and install them in *NODE,
   which is either a DECL (including a TYPE_DECL) or a TYPE.  If a DECL,
   it should be modified in place; if a TYPE, a copy should be created
   unless ATTR_FLAG_TYPE_IN_PLACE is set in FLAGS.  FLAGS gives further
   information, in the form of a bitwise OR of flags in enum attribute_flags
   from tree.h.  Depending on these flags, some attributes may be
   returned to be applied at a later stage (for example, to apply
   a decl attribute to the declaration rather than to its type).  */
extern tree decl_attributes (tree *, tree, int, tree = NULL_TREE);

extern bool cxx11_attribute_p (const_tree);
extern tree get_attribute_name (const_tree);
extern tree get_attribute_namespace (const_tree);
extern void apply_tm_attr (tree, tree);
extern tree make_attribute (const char *, const char *, tree);
extern bool attribute_ignored_p (tree);
extern bool attribute_ignored_p (const attribute_spec *const);

extern struct scoped_attributes* register_scoped_attributes (const struct attribute_spec *,
							     const char *,
							     bool = false);

extern char *sorted_attr_string (tree);
extern bool common_function_versions (tree, tree);
extern tree make_dispatcher_decl (const tree);
extern bool is_function_default_version (const tree);
extern void handle_ignored_attributes_option (vec<char *> *);

/* Return a type like TTYPE except that its TYPE_ATTRIBUTES
   is ATTRIBUTE.

   Such modified types already made are recorded so that duplicates
   are not made.  */

extern tree build_type_attribute_variant (tree, tree);
extern tree build_decl_attribute_variant (tree, tree);
extern tree build_type_attribute_qual_variant (tree, tree, int);

extern bool simple_cst_list_equal (const_tree, const_tree);
extern bool attribute_value_equal (const_tree, const_tree);

/* Return 0 if the attributes for two types are incompatible, 1 if they
   are compatible, and 2 if they are nearly compatible (which causes a
   warning to be generated).  */
extern int comp_type_attributes (const_tree, const_tree);

extern tree affects_type_identity_attributes (tree, bool = true);
extern tree restrict_type_identity_attributes_to (tree, tree);

/* Default versions of target-overridable functions.  */
extern tree merge_decl_attributes (tree, tree);
extern tree merge_type_attributes (tree, tree);

/* Remove any instances of attribute ATTR_NAME in LIST and return the
   modified list.  */

extern tree remove_attribute (const char *, tree);

/* Given two attributes lists, return a list of their union.  */

extern tree merge_attributes (tree, tree);

/* Duplicate all attributes with name NAME in ATTR list to *ATTRS if
   they are missing there.  */

extern void duplicate_one_attribute (tree *, tree, const char *);

/* Duplicate all attributes from user DECL to the corresponding
   builtin that should be propagated.  */

extern void copy_attributes_to_builtin (tree);

/* Given two Windows decl attributes lists, possibly including
   dllimport, return a list of their union .  */
extern tree merge_dllimport_decl_attributes (tree, tree);

/* Handle a "dllimport" or "dllexport" attribute.  */
extern tree handle_dll_attribute (tree *, tree, tree, int, bool *);

extern int attribute_list_equal (const_tree, const_tree);
extern int attribute_list_contained (const_tree, const_tree);

/* The backbone of lookup_attribute().  ATTR_LEN is the string length
   of ATTR_NAME, and LIST is not NULL_TREE.

   The function is called from lookup_attribute in order to optimize
   for size.  */
extern tree private_lookup_attribute (const char *attr_name, size_t attr_len,
				      tree list);

extern unsigned decls_mismatched_attributes (tree, tree, tree,
					     const char* const[],
					     pretty_printer*);

extern void maybe_diag_alias_attributes (tree, tree);

/* For a given string S of length L, strip leading and trailing '_' characters
   so that we have a canonical form of attribute names.  NB: This function may
   change S and L.  */

template <typename T>
inline bool
canonicalize_attr_name (const char *&s, T &l)
{
  if (l > 4 && s[0] == '_' && s[1] == '_' && s[l - 1] == '_' && s[l - 2] == '_')
    {
      s += 2;
      l -= 4;
      return true;
    }
  return false;
}

/* For a given IDENTIFIER_NODE, strip leading and trailing '_' characters
   so that we have a canonical form of attribute names.  */

static inline tree
canonicalize_attr_name (tree attr_name)
{
  size_t l = IDENTIFIER_LENGTH (attr_name);
  const char *s = IDENTIFIER_POINTER (attr_name);

  if (canonicalize_attr_name (s, l))
    return get_identifier_with_length (s, l);

  return attr_name;
}

/* Compare attribute identifiers ATTR1 and ATTR2 with length ATTR1_LEN and
   ATTR2_LEN.  */

static inline bool
cmp_attribs (const char *attr1, size_t attr1_len,
	     const char *attr2, size_t attr2_len)
{
  return attr1_len == attr2_len && strncmp (attr1, attr2, attr1_len) == 0;
}

/* Compare attribute identifiers ATTR1 and ATTR2.  */

static inline bool
cmp_attribs (const char *attr1, const char *attr2)
{
  return cmp_attribs (attr1, strlen (attr1), attr2, strlen (attr2));
}

/* Given an identifier node IDENT and a string ATTR_NAME, return true
   if the identifier node is a valid attribute name for the string.  */

static inline bool
is_attribute_p (const char *attr_name, const_tree ident)
{
  return cmp_attribs (attr_name, strlen (attr_name),
		      IDENTIFIER_POINTER (ident), IDENTIFIER_LENGTH (ident));
}

/* Given an attribute name ATTR_NAME and a list of attributes LIST,
   return a pointer to the attribute's list element if the attribute
   is part of the list, or NULL_TREE if not found.  If the attribute
   appears more than once, this only returns the first occurrence; the
   TREE_CHAIN of the return value should be passed back in if further
   occurrences are wanted.  ATTR_NAME must be in the form 'text' (not
   '__text__').  */

static inline tree
lookup_attribute (const char *attr_name, tree list)
{
  if (CHECKING_P && attr_name[0] != '_')
    {
      size_t attr_len = strlen (attr_name);
      gcc_checking_assert (!canonicalize_attr_name (attr_name, attr_len));
    }
  /* In most cases, list is NULL_TREE.  */
  if (list == NULL_TREE)
    return NULL_TREE;
  else
    {
      size_t attr_len = strlen (attr_name);
      /* Do the strlen() before calling the out-of-line implementation.
	 In most cases attr_name is a string constant, and the compiler
	 will optimize the strlen() away.  */
      return private_lookup_attribute (attr_name, attr_len, list);
    }
}

/* Given an attribute name ATTR_NAME and a list of attributes LIST,
   return a pointer to the attribute's list first element if the attribute
   starts with ATTR_NAME.  ATTR_NAME must be in the form 'text' (not
   '__text__').  */

static inline tree
lookup_attribute_by_prefix (const char *attr_name, tree list)
{
  gcc_checking_assert (attr_name[0] != '_');
  /* In most cases, list is NULL_TREE.  */
  if (list == NULL_TREE)
    return NULL_TREE;
  else
    {
      size_t attr_len = strlen (attr_name);
      while (list)
	{
	  tree name = get_attribute_name (list);
	  size_t ident_len = IDENTIFIER_LENGTH (name);

	  if (attr_len > ident_len)
	    {
	      list = TREE_CHAIN (list);
	      continue;
	    }

	  const char *p = IDENTIFIER_POINTER (name);
	  gcc_checking_assert (attr_len == 0 || p[0] != '_');

	  if (strncmp (attr_name, p, attr_len) == 0)
	    break;

	  list = TREE_CHAIN (list);
	}

      return list;
    }
}

/* Description of a function argument declared with attribute access.
   Used as an "iterator" over all such arguments in a function declaration
   or call.  */

struct attr_access
{
  /* The beginning and end of the internal string representation.  */
  const char *str, *end;
  /* The attribute pointer argument.  */
  tree ptr;
  /* For a declaration, a TREE_CHAIN of VLA bound expressions stored
     in TREE_VALUE and their positions in the argument list (stored
     in TREE_PURPOSE).  Each expression may be a PARM_DECL or some
     other DECL (for ordinary variables), or an EXPR for other
     expressions (e.g., funcion calls).  */
  tree size;

  /* The zero-based position of each of the formal function arguments.
     For the optional SIZARG, UINT_MAX when not specified.  For VLAs
     with multiple variable bounds, SIZARG is the position corresponding
     to the most significant bound in the argument list.  Positions of
     subsequent bounds are in the TREE_PURPOSE field of the SIZE chain.  */
  unsigned ptrarg;
  unsigned sizarg;
  /* For internal specifications only, the constant minimum size of
     the array, zero if not specified, and HWI_M1U for the unspecified
     VLA [*] notation.  Meaningless for external (explicit) access
     specifications.  */
  unsigned HOST_WIDE_INT minsize;

  /* The access mode.  */
  access_mode mode;

  /* Set for an attribute added internally rather than by an explicit
     declaration. */
  bool internal_p;
  /* Set for the T[static MINSIZE] array notation for nonzero MINSIZE
     less than HWI_M1U.  */
  bool static_p;

  /* Return the number of specified VLA bounds.  */
  unsigned vla_bounds (unsigned *) const;

  /* Return internal representation as STRING_CST.  */
  tree to_internal_string () const;

  /* Return the human-readable representation of the external attribute
     specification (as it might appear in the source code) as STRING_CST.  */
  tree to_external_string () const;

  /* Return argument of array type formatted as a readable string.  */
  std::string array_as_string (tree) const;

  /* Return the access mode corresponding to the character code.  */
  static access_mode from_mode_char (char);

  /* Reset front end-specific attribute access data from attributes.  */
  static void free_lang_data (tree);

  /* The character codes corresponding to all the access modes.  */
  static constexpr char mode_chars[5] = { '-', 'r', 'w', 'x', '^' };

  /* The strings corresponding to just the external access modes.  */
  static constexpr char mode_names[4][11] =
    {
     "none", "read_only", "write_only", "read_write"
    };
};

inline access_mode
attr_access::from_mode_char (char c)
{
  switch (c)
    {
    case mode_chars[access_none]: return access_none;
    case mode_chars[access_read_only]: return access_read_only;
    case mode_chars[access_write_only]: return access_write_only;
    case mode_chars[access_read_write]: return access_read_write;
    case mode_chars[access_deferred]: return access_deferred;
    }
  gcc_unreachable ();
}

/* Used to define rdwr_map below.  */
struct rdwr_access_hash: int_hash<int, -1> { };

/* A mapping between argument number corresponding to attribute access
   mode (read_only, write_only, or read_write) and operands.  */
struct attr_access;
typedef hash_map<rdwr_access_hash, attr_access> rdwr_map;

extern void init_attr_rdwr_indices (rdwr_map *, tree);
extern attr_access *get_parm_access (rdwr_map &, tree,
				     tree = current_function_decl);

#endif // GCC_ATTRIBS_H
