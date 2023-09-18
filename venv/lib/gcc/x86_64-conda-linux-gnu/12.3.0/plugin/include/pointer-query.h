/* Definitions of the pointer_query and related classes.

   Copyright (C) 2020-2022 Free Software Foundation, Inc.

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

#ifndef GCC_POINTER_QUERY_H
#define GCC_POINTER_QUERY_H

/* Describes recursion limits used by functions that follow use-def
   chains of SSA_NAMEs.  */

class ssa_name_limit_t
{
  bitmap visited;         /* Bitmap of visited SSA_NAMEs.  */
  unsigned ssa_def_max;   /* Longest chain of SSA_NAMEs to follow.  */

  /* Not copyable or assignable.  */
  DISABLE_COPY_AND_ASSIGN (ssa_name_limit_t);

public:

  ssa_name_limit_t ()
    : visited (),
      ssa_def_max (param_ssa_name_def_chain_limit) { }

  /* Set a bit for the PHI in VISITED and return true if it wasn't
     already set.  */
  bool visit_phi (tree);
  /* Clear a bit for the PHI in VISITED.  */
  void leave_phi (tree);
  /* Return false if the SSA_NAME chain length counter has reached
     the limit, otherwise increment the counter and return true.  */
  bool next ();

  /* If the SSA_NAME has already been "seen" return a positive value.
     Otherwise add it to VISITED.  If the SSA_NAME limit has been
     reached, return a negative value.  Otherwise return zero.  */
  int next_phi (tree);

  ~ssa_name_limit_t ();
};

class pointer_query;

/* Describes a reference to an object used in an access.  */
struct access_ref
{
  /* Set the bounds of the reference.  */
  access_ref ();

  /* Return the PHI node REF refers to or null if it doesn't.  */
  gphi *phi () const;

  /* Merge the result for a pointer with *THIS.  */
  void merge_ref (vec<access_ref> *all_refs, tree, gimple *, int, bool,
		  ssa_name_limit_t &, pointer_query &);

  /* Return the object to which REF refers.  */
  tree get_ref (vec<access_ref> *, access_ref * = nullptr, int = 1,
		ssa_name_limit_t * = nullptr, pointer_query * = nullptr) const;

  /* Return true if OFFRNG is the constant zero.  */
  bool offset_zero () const
  {
    return offrng[0] == 0 && offrng[1] == 0;
  }

  /* Return true if OFFRNG is bounded to a subrange of offset values
     valid for the largest possible object.  */
  bool offset_bounded () const;

  /* Return the maximum amount of space remaining and if non-null, set
     argument to the minimum.  */
  offset_int size_remaining (offset_int * = nullptr) const;

/* Return true if the offset and object size are in range for SIZE.  */
  bool offset_in_range (const offset_int &) const;

  /* Return true if *THIS is an access to a declared object.  */
  bool ref_declared () const
  {
    return DECL_P (ref) && base0 && deref < 1;
  }

  /* Set the size range to the maximum.  */
  void set_max_size_range ()
  {
    sizrng[0] = 0;
    sizrng[1] = wi::to_offset (max_object_size ());
  }

  /* Add OFF to the offset range.  */
  void add_offset (const offset_int &off)
  {
    add_offset (off, off);
  }

  /* Add the range [MIN, MAX] to the offset range.  */
  void add_offset (const offset_int &, const offset_int &);

  /* Add the maximum representable offset to the offset range.  */
  void add_max_offset ()
  {
    offset_int maxoff = wi::to_offset (TYPE_MAX_VALUE (ptrdiff_type_node));
    add_offset (-maxoff - 1, maxoff);
  }

  /* Issue an informational message describing the target of an access
     with the given mode.  */
  void inform_access (access_mode, int = 1) const;

  /* Dump *THIS to a file.  */
  void dump (FILE *) const;

  /* Reference to the accessed object(s).  */
  tree ref;

  /* Range of byte offsets into and sizes of the object(s).  */
  offset_int offrng[2];
  offset_int sizrng[2];
  /* The minimum and maximum offset computed.  */
  offset_int offmax[2];

  /* Used to fold integer expressions when called from front ends.  */
  tree (*eval)(tree);
  /* Positive when REF is dereferenced, negative when its address is
     taken.  */
  int deref;
  /* Set if trailing one-element arrays should be treated as flexible
     array members.  */
  bool trail1special;
  /* Set if valid offsets must start at zero (for declared and allocated
     objects but not for others referenced by pointers).  */
  bool base0;
  /* Set if REF refers to a function array parameter not declared
     static.  */
  bool parmarray;
};

class range_query;

/* Queries and caches compute_objsize results.  */
class pointer_query
{
  DISABLE_COPY_AND_ASSIGN (pointer_query);

  /* Type of the two-level cache object defined by clients of the class
     to have pointer SSA_NAMEs cached for speedy access.  */
  struct cache_type
  {
    /* 1-based indices into cache.  */
    auto_vec<unsigned> indices;
    /* The cache itself.  */
    auto_vec<access_ref> access_refs;
  };

public:
  /* Construct an object with the given Ranger instance.  */
  explicit pointer_query (range_query * = nullptr);

  /* Retrieve the access_ref for a variable from cache if it's there.  */
  const access_ref* get_ref (tree, int = 1) const;

  /* Retrieve the access_ref for a variable from cache or compute it.  */
  bool get_ref (tree, gimple *, access_ref*, int = 1);

  /* Add an access_ref for the SSA_NAME to the cache.  */
  void put_ref (tree, const access_ref&, int = 1);

  /* Flush the cache.  */
  void flush_cache ();

  /* Dump statistics and optionally cache contents to DUMP_FILE.  */
  void dump (FILE *, bool = false);

  /* A Ranger instance.  May be null to use global ranges.  */
  range_query *rvals;

  /* Cache performance counters.  */
  mutable unsigned hits;
  mutable unsigned misses;
  mutable unsigned failures;
  mutable unsigned depth;
  mutable unsigned max_depth;

private:
  /* Cache of SSA_NAMEs.  May be null to disable caching.  */
  cache_type var_cache;
};

/* Describes a pair of references used in an access by built-in
   functions like memcpy.  */
struct access_data
{
  /* Set the access to at most MAXWRITE and MAXREAD bytes, and
     at least 1 when MINWRITE or MINREAD, respectively, is set.  */
  access_data (range_query *, gimple *, access_mode,
	       tree = NULL_TREE, bool = false,
	       tree = NULL_TREE, bool = false);

  /* Set the access to at most MAXWRITE and MAXREAD bytes, and
     at least 1 when MINWRITE or MINREAD, respectively, is set.  */
  access_data (range_query *, tree, access_mode,
	       tree = NULL_TREE, bool = false,
	       tree = NULL_TREE, bool = false);

  /* Constructor helper.  */
  static void set_bound (offset_int[2], tree, bool, range_query *, gimple *);

  /* Access statement.  */
  gimple *stmt;
  /* Built-in function call.  */
  tree call;
  /* Destination and source of the access.  */
  access_ref dst, src;

  /* Range of the bound of the access: denotes that the access is at
     least XXX_BNDRNG[0] bytes but no more than XXX_BNDRNG[1].  For
     string functions the size of the actual access is further
     constrained by the length of the string.  */
  offset_int dst_bndrng[2];
  offset_int src_bndrng[2];

  /* Read-only for functions like memcmp or strlen, write-only
     for memset, read-write for memcpy or strcat.  */
  access_mode mode;
  /* The object size type.  */
  int ostype;
};

enum size_range_flags
  {
   /* Set to consider zero a valid range.  */
   SR_ALLOW_ZERO = 1,
   /* Set to use the largest subrange of a set of ranges as opposed
      to the smallest.  */
   SR_USE_LARGEST = 2
  };
extern bool get_size_range (tree, tree[2], int = 0);
extern bool get_size_range (range_query *, tree, gimple *, tree[2], int = 0);

class range_query;
extern tree gimple_call_alloc_size (gimple *, wide_int[2] = nullptr,
				    range_query * = nullptr);

/* Compute the size of an object referenced by the first argument in
   a statement given by second argument, using Object Size Type given
   by third argument.  Store result in an access_ref.  */
extern tree compute_objsize (tree, gimple *, int, access_ref *,
			     range_query * = nullptr);
extern tree compute_objsize (tree, gimple *, int, access_ref *,
			     pointer_query *);
inline tree compute_objsize (tree ptr, int ostype, access_ref *pref)
{
  return compute_objsize (ptr, nullptr, ostype, pref, (range_query *)nullptr);
}

/* Legacy/transitional API.  Should not be used in new code.  */
extern tree compute_objsize (tree, gimple *, int, tree * = nullptr,
			     tree * = nullptr, range_query * = nullptr);
inline tree compute_objsize (tree ptr, int ostype, tree *pdecl = nullptr,
			     tree *poff = nullptr, range_query *rvals = nullptr)
{
  return compute_objsize (ptr, nullptr, ostype, pdecl, poff, rvals);
}

/* Return the field at the constant offset.  */
extern tree field_at_offset (tree, tree, HOST_WIDE_INT,
			     HOST_WIDE_INT * = nullptr,
			     HOST_WIDE_INT * = nullptr);
/* Return the array at the constant offset.  */
extern tree array_elt_at_offset (tree, HOST_WIDE_INT,
				 HOST_WIDE_INT * = nullptr,
				 HOST_WIDE_INT * = nullptr);

/* Helper to build an array type that can be printed.  */
extern tree build_printable_array_type (tree, unsigned HOST_WIDE_INT);

#endif   // GCC_POINTER_QUERY_H
