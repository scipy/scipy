/* Data structure for the modref pass.
   Copyright (C) 2020-2022 Free Software Foundation, Inc.
   Contributed by David Cepelik and Jan Hubicka

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

/* modref_tree represent a decision tree that can be used by alias analysis
   oracle to determine whether given memory access can be affected by a function
   call.  For every function we collect two trees, one for loads and other
   for stores.  Tree consist of following levels:

   1) Base: this level represent base alias set of the access and refers
      to sons (ref nodes). Flag all_refs means that all possible references
      are aliasing.

      Because for LTO streaming we need to stream types rather than alias sets
      modref_base_node is implemented as a template.
   2) Ref: this level represent ref alias set and links to accesses unless
      all_refs flag is set.
      Again ref is an template to allow LTO streaming.
   3) Access: this level represent info about individual accesses.  Presently
      we record whether access is through a dereference of a function parameter
      and if so we record the access range.
*/

#ifndef GCC_MODREF_TREE_H
#define GCC_MODREF_TREE_H

struct ipa_modref_summary;

/* parm indexes greater than 0 are normal parms.
   Some negative values have special meaning.  */
enum modref_special_parms {
  MODREF_UNKNOWN_PARM = -1,
  MODREF_STATIC_CHAIN_PARM = -2,
  MODREF_RETSLOT_PARM = -3,
  /* Used for bases that points to memory that escapes from function.  */
  MODREF_GLOBAL_MEMORY_PARM = -4,
  /* Used in modref_parm_map to take references which can be removed
     from the summary during summary update since they now points to local
     memory.  */
  MODREF_LOCAL_MEMORY_PARM = -5
};

/* Modref record accesses relative to function parameters.
   This is entry for single access specifying its base and access range.

   Accesses can be collected to boundedly sized arrays using
   modref_access_node::insert.  */
struct GTY(()) modref_access_node
{
  /* Access range information (in bits).  */
  poly_int64 offset;
  poly_int64 size;
  poly_int64 max_size;

  /* Offset from parameter pointer to the base of the access (in bytes).  */
  poly_int64 parm_offset;

  /* Index of parameter which specifies the base of access. -1 if base is not
     a function parameter.  */
  int parm_index;
  bool parm_offset_known;
  /* Number of times interval was extended during dataflow.
     This has to be limited in order to keep dataflow finite.  */
  unsigned char adjustments;

  /* Return true if access node holds some useful info.  */
  bool useful_p () const
    {
      return parm_index != MODREF_UNKNOWN_PARM;
    }
  /* Return true if access can be used to determine a kill.  */
  bool useful_for_kill_p () const
    {
      return parm_offset_known && parm_index != MODREF_UNKNOWN_PARM
	     && parm_index != MODREF_GLOBAL_MEMORY_PARM
	     && parm_index != MODREF_RETSLOT_PARM && known_size_p (size)
	     && known_eq (max_size, size)
	     && known_gt (size, 0);
    }
  /* Dump range to debug OUT.  */
  void dump (FILE *out);
  /* Return true if both accesses are the same.  */
  bool operator == (modref_access_node &a) const;
  /* Return true if range info is useful.  */
  bool range_info_useful_p () const;
  /* Return tree corresponding to parameter of the range in STMT.  */
  tree get_call_arg (const gcall *stmt) const;
  /* Build ao_ref corresponding to the access and return true if successful.  */
  bool get_ao_ref (const gcall *stmt, class ao_ref *ref) const;
  /* Stream access to OB.  */
  void stream_out (struct output_block *ob) const;
  /* Stream access in from IB.  */
  static modref_access_node stream_in (struct lto_input_block *ib);
  /* Insert A into vector ACCESSES.  Limit size of vector to MAX_ACCESSES and
     if RECORD_ADJUSTMENT is true keep track of adjustment counts.
     Return 0 if nothing changed, 1 is insertion succeeded and -1 if failed.  */
  static int insert (vec <modref_access_node, va_gc> *&accesses,
		     modref_access_node a, size_t max_accesses,
		     bool record_adjustments);
  /* Same as insert but for kills where we are conservative the other way
     around: if information is lost, the kill is lost.  */
  static bool insert_kill (vec<modref_access_node> &kills,
			   modref_access_node &a, bool record_adjustments);
private:
  bool contains (const modref_access_node &) const;
  bool contains_for_kills (const modref_access_node &) const;
  void update (poly_int64, poly_int64, poly_int64, poly_int64, bool);
  bool update_for_kills (poly_int64, poly_int64, poly_int64,
			 poly_int64, poly_int64, bool);
  bool merge (const modref_access_node &, bool);
  bool merge_for_kills (const modref_access_node &, bool);
  static bool closer_pair_p (const modref_access_node &,
			     const modref_access_node &,
			     const modref_access_node &,
			     const modref_access_node &);
  void forced_merge (const modref_access_node &, bool);
  void update2 (poly_int64, poly_int64, poly_int64, poly_int64,
		poly_int64, poly_int64, poly_int64, bool);
  bool combined_offsets (const modref_access_node &,
			 poly_int64 *, poly_int64 *, poly_int64 *) const;
  static void try_merge_with (vec <modref_access_node, va_gc> *&, size_t);
};

/* Access node specifying no useful info.  */
const modref_access_node unspecified_modref_access_node
		 = {0, -1, -1, 0, MODREF_UNKNOWN_PARM, false, 0};

template <typename T>
struct GTY((user)) modref_ref_node
{
  T ref;
  bool every_access;
  vec <modref_access_node, va_gc> *accesses;

  modref_ref_node (T ref):
    ref (ref),
    every_access (false),
    accesses (NULL)
  {}

  /* Collapse the tree.  */
  void collapse ()
  {
    vec_free (accesses);
    accesses = NULL;
    every_access = true;
  }

  /* Insert access with OFFSET and SIZE.
     Collapse tree if it has more than MAX_ACCESSES entries.
     If RECORD_ADJUSTMENTs is true avoid too many interval extensions.
     Return true if record was changed.  */
  bool insert_access (modref_access_node a, size_t max_accesses,
		      bool record_adjustments)
  {
    /* If this base->ref pair has no access information, bail out.  */
    if (every_access)
      return false;

    /* Only the following kind of parameters needs to be tracked.
       We do not track return slots because they are seen as a direct store
       in the caller.  */
    gcc_checking_assert (a.parm_index >= 0
			 || a.parm_index == MODREF_STATIC_CHAIN_PARM
			 || a.parm_index == MODREF_GLOBAL_MEMORY_PARM
			 || a.parm_index == MODREF_UNKNOWN_PARM);

    if (!a.useful_p ())
      {
	if (!every_access)
	  {
	    collapse ();
	    return true;
	  }
	return false;
      }

    int ret = modref_access_node::insert (accesses, a, max_accesses,
					  record_adjustments);
    if (ret == -1)
      {
	if (dump_file)
	  fprintf (dump_file,
		   "--param modref-max-accesses limit reached; collapsing\n");
	collapse ();
      }
    return ret != 0;
  }
};

/* Base of an access.  */
template <typename T>
struct GTY((user)) modref_base_node
{
  T base;
  vec <modref_ref_node <T> *, va_gc> *refs;
  bool every_ref;

  modref_base_node (T base):
    base (base),
    refs (NULL),
    every_ref (false) {}

  /* Search REF; return NULL if failed.  */
  modref_ref_node <T> *search (T ref)
  {
    size_t i;
    modref_ref_node <T> *n;
    FOR_EACH_VEC_SAFE_ELT (refs, i, n)
      if (n->ref == ref)
	return n;
    return NULL;
  }

  /* Insert REF; collapse tree if there are more than MAX_REFS.
     Return inserted ref and if CHANGED is non-null set it to true if
     something changed.  */
  modref_ref_node <T> *insert_ref (T ref, size_t max_refs,
				   bool *changed = NULL)
  {
    modref_ref_node <T> *ref_node;

    /* If the node is collapsed, don't do anything.  */
    if (every_ref)
      return NULL;

    /* Otherwise, insert a node for the ref of the access under the base.  */
    ref_node = search (ref);
    if (ref_node)
      return ref_node;

    /* We always allow inserting ref 0.  For non-0 refs there is upper
       limit on number of entries and if exceeded,
       drop ref conservatively to 0.  */
    if (ref && refs && refs->length () >= max_refs)
      {
	if (dump_file)
	  fprintf (dump_file, "--param modref-max-refs limit reached;"
		   " using 0\n");
	ref = 0;
	ref_node = search (ref);
	if (ref_node)
	  return ref_node;
      }

    if (changed)
      *changed = true;

    ref_node = new (ggc_alloc <modref_ref_node <T> > ())modref_ref_node <T>
								 (ref);
    vec_safe_push (refs, ref_node);
    return ref_node;
  }

  void collapse ()
  {
    size_t i;
    modref_ref_node <T> *r;

    if (refs)
      {
	FOR_EACH_VEC_SAFE_ELT (refs, i, r)
	  {
	    r->collapse ();
	    ggc_free (r);
	  }
	vec_free (refs);
      }
    refs = NULL;
    every_ref = true;
  }
};

/* Map translating parameters across function call.  */

struct modref_parm_map
{
  /* Default constructor.  */
  modref_parm_map ()
  : parm_index (MODREF_UNKNOWN_PARM), parm_offset_known (false), parm_offset ()
  {}

  /* Index of parameter we translate to.
     Values from special_params enum are permitted too.  */
  int parm_index;
  bool parm_offset_known;
  poly_int64 parm_offset;
};

/* Access tree for a single function.  */
template <typename T>
struct GTY((user)) modref_tree
{
  vec <modref_base_node <T> *, va_gc> *bases;
  bool every_base;

  modref_tree ():
    bases (NULL),
    every_base (false) {}

  /* Insert BASE; collapse tree if there are more than MAX_REFS.
     Return inserted base and if CHANGED is non-null set it to true if
     something changed.
     If table gets full, try to insert REF instead.  */

  modref_base_node <T> *insert_base (T base, T ref,
				     unsigned int max_bases,
				     bool *changed = NULL)
  {
    modref_base_node <T> *base_node;

    /* If the node is collapsed, don't do anything.  */
    if (every_base)
      return NULL;

    /* Otherwise, insert a node for the base of the access into the tree.  */
    base_node = search (base);
    if (base_node)
      return base_node;

    /* We always allow inserting base 0.  For non-0 base there is upper
       limit on number of entries and if exceeded,
       drop base conservatively to ref and if it still does not fit to 0.  */
    if (base && bases && bases->length () >= max_bases)
      {
	base_node = search (ref);
	if (base_node)
	  {
	    if (dump_file)
	      fprintf (dump_file, "--param modref-max-bases"
		       " limit reached; using ref\n");
	    return base_node;
	  }
	if (dump_file)
	  fprintf (dump_file, "--param modref-max-bases"
		   " limit reached; using 0\n");
	base = 0;
	base_node = search (base);
	if (base_node)
	  return base_node;
      }

    if (changed)
      *changed = true;

    base_node = new (ggc_alloc <modref_base_node <T> > ())
			 modref_base_node <T> (base);
    vec_safe_push (bases, base_node);
    return base_node;
  }

  /* Insert memory access to the tree.
     Return true if something changed.  */
  bool insert (unsigned int max_bases,
	       unsigned int max_refs,
	       unsigned int max_accesses,
	       T base, T ref, modref_access_node a,
	       bool record_adjustments)
  {
    if (every_base)
      return false;

    bool changed = false;

    /* We may end up with max_size being less than size for accesses past the
       end of array.  Those are undefined and safe to ignore.  */
    if (a.range_info_useful_p ()
	&& known_size_p (a.size) && known_size_p (a.max_size)
	&& known_lt (a.max_size, a.size))
      {
	if (dump_file)
	  fprintf (dump_file,
		   "   - Paradoxical range. Ignoring\n");
	return false;
      }
    if (known_size_p (a.size)
	&& known_eq (a.size, 0))
      {
	if (dump_file)
	  fprintf (dump_file,
		   "   - Zero size. Ignoring\n");
	return false;
      }
    if (known_size_p (a.max_size)
	&& known_eq (a.max_size, 0))
      {
	if (dump_file)
	  fprintf (dump_file,
		   "   - Zero max_size. Ignoring\n");
	return false;
      }
    gcc_checking_assert (!known_size_p (a.max_size)
			 || !known_le (a.max_size, 0));

    /* No useful information tracked; collapse everything.  */
    if (!base && !ref && !a.useful_p ())
      {
	collapse ();
	return true;
      }

    modref_base_node <T> *base_node
      = insert_base (base, ref, max_bases, &changed);
    base = base_node->base;
    /* If table got full we may end up with useless base.  */
    if (!base && !ref && !a.useful_p ())
      {
	collapse ();
	return true;
      }
    if (base_node->every_ref)
      return changed;
    gcc_checking_assert (search (base) != NULL);

    /* No useful ref info tracked; collapse base.  */
    if (!ref && !a.useful_p ())
      {
	base_node->collapse ();
	return true;
      }

    modref_ref_node <T> *ref_node
	    = base_node->insert_ref (ref, max_refs, &changed);
    ref = ref_node->ref;

    if (ref_node->every_access)
      return changed;
    changed |= ref_node->insert_access (a, max_accesses,
					record_adjustments);
    /* See if we failed to add useful access.  */
    if (ref_node->every_access)
      {
	/* Collapse everything if there is no useful base and ref.  */
	if (!base && !ref)
	  {
	    collapse ();
	    gcc_checking_assert (changed);
	  }
	/* Collapse base if there is no useful ref.  */
	else if (!ref)
	  {
	    base_node->collapse ();
	    gcc_checking_assert (changed);
	  }
      }
    return changed;
  }

  /* Insert memory access to the tree.
     Return true if something changed.  */
  bool insert (tree fndecl,
	       T base, T ref, const modref_access_node &a,
	       bool record_adjustments)
  {
     return insert (opt_for_fn (fndecl, param_modref_max_bases),
		    opt_for_fn (fndecl, param_modref_max_refs),
		    opt_for_fn (fndecl, param_modref_max_accesses),
		    base, ref, a, record_adjustments);
  }

 /* Remove tree branches that are not useful (i.e. they will always pass).  */

 void cleanup ()
 {
   size_t i, j;
   modref_base_node <T> *base_node;
   modref_ref_node <T> *ref_node;

   if (!bases)
     return;

   for (i = 0; vec_safe_iterate (bases, i, &base_node);)
     {
       if (base_node->refs)
	 for (j = 0; vec_safe_iterate (base_node->refs, j, &ref_node);)
	   {
	     if (!ref_node->every_access
		 && (!ref_node->accesses
		     || !ref_node->accesses->length ()))
	       {
		 base_node->refs->unordered_remove (j);
		 vec_free (ref_node->accesses);
		 ggc_delete (ref_node);
	       }
	     else
	       j++;
	   }
       if (!base_node->every_ref
	   && (!base_node->refs || !base_node->refs->length ()))
	 {
	   bases->unordered_remove (i);
	   vec_free (base_node->refs);
	   ggc_delete (base_node);
	 }
       else
	 i++;
     }
   if (bases && !bases->length ())
     {
       vec_free (bases);
       bases = NULL;
     }
 }

  /* Merge OTHER into the tree.
     PARM_MAP, if non-NULL, maps parm indexes of callee to caller.
     Similar CHAIN_MAP, if non-NULL, maps static chain of callee to caller.
     Return true if something has changed.  */
  bool merge (unsigned int max_bases,
	      unsigned int max_refs,
	      unsigned int max_accesses,
	      modref_tree <T> *other, vec <modref_parm_map> *parm_map,
	      modref_parm_map *static_chain_map,
	      bool record_accesses,
	      bool promote_unknown_to_global = false)
  {
    if (!other || every_base)
      return false;
    if (other->every_base)
      {
	collapse ();
	return true;
      }

    bool changed = false;
    size_t i, j, k;
    modref_base_node <T> *base_node, *my_base_node;
    modref_ref_node <T> *ref_node;
    modref_access_node *access_node;
    bool release = false;

    /* For self-recursive functions we may end up merging summary into itself;
       produce copy first so we do not modify summary under our own hands.  */
    if (other == this)
      {
	release = true;
	other = modref_tree<T>::create_ggc ();
	other->copy_from (this);
      }

    FOR_EACH_VEC_SAFE_ELT (other->bases, i, base_node)
      {
	if (base_node->every_ref)
	  {
	    my_base_node = insert_base (base_node->base, 0,
					max_bases, &changed);
	    if (my_base_node && !my_base_node->every_ref)
	      {
		my_base_node->collapse ();
		cleanup ();
		changed = true;
	      }
	  }
	else
	  FOR_EACH_VEC_SAFE_ELT (base_node->refs, j, ref_node)
	    {
	      if (ref_node->every_access)
		{
		  changed |= insert (max_bases, max_refs, max_accesses,
				     base_node->base,
				     ref_node->ref,
				     unspecified_modref_access_node,
				     record_accesses);
		}
	      else
		FOR_EACH_VEC_SAFE_ELT (ref_node->accesses, k, access_node)
		  {
		    modref_access_node a = *access_node;

		    if (a.parm_index != MODREF_UNKNOWN_PARM
			&& a.parm_index != MODREF_GLOBAL_MEMORY_PARM
			&& parm_map)
		      {
			if (a.parm_index >= (int)parm_map->length ())
			  a.parm_index = MODREF_UNKNOWN_PARM;
			else
			  {
			    modref_parm_map &m
				    = a.parm_index == MODREF_STATIC_CHAIN_PARM
				      ? *static_chain_map
				      : (*parm_map) [a.parm_index];
			    if (m.parm_index == MODREF_LOCAL_MEMORY_PARM)
			      continue;
			    a.parm_offset += m.parm_offset;
			    a.parm_offset_known &= m.parm_offset_known;
			    a.parm_index = m.parm_index;
			  }
		      }
		    if (a.parm_index == MODREF_UNKNOWN_PARM
			&& promote_unknown_to_global)
		      a.parm_index = MODREF_GLOBAL_MEMORY_PARM;
		    changed |= insert (max_bases, max_refs, max_accesses,
				       base_node->base, ref_node->ref,
				       a, record_accesses);
		  }
	    }
      }
    if (release)
      ggc_delete (other);
    return changed;
  }

  /* Merge OTHER into the tree.
     PARM_MAP, if non-NULL, maps parm indexes of callee to caller.
     Similar CHAIN_MAP, if non-NULL, maps static chain of callee to caller.
     Return true if something has changed.  */
  bool merge (tree fndecl,
	      modref_tree <T> *other, vec <modref_parm_map> *parm_map,
	      modref_parm_map *static_chain_map,
	      bool record_accesses,
	      bool promote_unknown_to_global = false)
  {
     return merge (opt_for_fn (fndecl, param_modref_max_bases),
		   opt_for_fn (fndecl, param_modref_max_refs),
		   opt_for_fn (fndecl, param_modref_max_accesses),
		   other, parm_map, static_chain_map, record_accesses,
		   promote_unknown_to_global);
  }

  /* Copy OTHER to THIS.  */
  void copy_from (modref_tree <T> *other)
  {
    merge (INT_MAX, INT_MAX, INT_MAX, other, NULL, NULL, false);
  }

  /* Search BASE in tree; return NULL if failed.  */
  modref_base_node <T> *search (T base)
  {
    size_t i;
    modref_base_node <T> *n;
    FOR_EACH_VEC_SAFE_ELT (bases, i, n)
      if (n->base == base)
	return n;
    return NULL;
  }

  /* Return true if tree contains access to global memory.  */
  bool global_access_p ()
  {
    size_t i, j, k;
    modref_base_node <T> *base_node;
    modref_ref_node <T> *ref_node;
    modref_access_node *access_node;
    if (every_base)
      return true;
    FOR_EACH_VEC_SAFE_ELT (bases, i, base_node)
      {
	if (base_node->every_ref)
	  return true;
	FOR_EACH_VEC_SAFE_ELT (base_node->refs, j, ref_node)
	  {
	    if (ref_node->every_access)
	      return true;
	    FOR_EACH_VEC_SAFE_ELT (ref_node->accesses, k, access_node)
	      if (access_node->parm_index == MODREF_UNKNOWN_PARM
		  || access_node->parm_index == MODREF_GLOBAL_MEMORY_PARM)
		return true;
	  }
      }
    return false;
  }

  /* Return ggc allocated instance.  We explicitly call destructors via
     ggc_delete and do not want finalizers to be registered and
     called at the garbage collection time.  */
  static modref_tree<T> *create_ggc ()
  {
    return new (ggc_alloc_no_dtor<modref_tree<T>> ())
	 modref_tree<T> ();
  }

  /* Remove all records and mark tree to alias with everything.  */
  void collapse ()
  {
    size_t i;
    modref_base_node <T> *n;

    if (bases)
      {
	FOR_EACH_VEC_SAFE_ELT (bases, i, n)
	  {
	    n->collapse ();
	    ggc_free (n);
	  }
	vec_free (bases);
      }
    bases = NULL;
    every_base = true;
  }

  /* Release memory.  */
  ~modref_tree ()
  {
    collapse ();
  }

  /* Update parameter indexes in TT according to MAP.  */
  void
  remap_params (vec <int> *map)
  {
    size_t i;
    modref_base_node <T> *base_node;
    FOR_EACH_VEC_SAFE_ELT (bases, i, base_node)
      {
	size_t j;
	modref_ref_node <T> *ref_node;
	FOR_EACH_VEC_SAFE_ELT (base_node->refs, j, ref_node)
	  {
	    size_t k;
	    modref_access_node *access_node;
	    FOR_EACH_VEC_SAFE_ELT (ref_node->accesses, k, access_node)
	      if (access_node->parm_index >= 0)
		{
		  if (access_node->parm_index < (int)map->length ())
		    access_node->parm_index = (*map)[access_node->parm_index];
		  else
		    access_node->parm_index = MODREF_UNKNOWN_PARM;
		}
	  }
      }
  }
};

void gt_ggc_mx (modref_tree <int>* const&);
void gt_ggc_mx (modref_tree <tree_node*>* const&);
void gt_pch_nx (modref_tree <int>* const&);
void gt_pch_nx (modref_tree <tree_node*>* const&);
void gt_pch_nx (modref_tree <int>* const&, gt_pointer_operator op, void *cookie);
void gt_pch_nx (modref_tree <tree_node*>* const&, gt_pointer_operator op,
		void *cookie);

void gt_ggc_mx (modref_base_node <int>*);
void gt_ggc_mx (modref_base_node <tree_node*>* &);
void gt_pch_nx (modref_base_node <int>* const&);
void gt_pch_nx (modref_base_node <tree_node*>* const&);
void gt_pch_nx (modref_base_node <int>* const&, gt_pointer_operator op,
		void *cookie);
void gt_pch_nx (modref_base_node <tree_node*>* const&, gt_pointer_operator op,
		void *cookie);

void gt_ggc_mx (modref_ref_node <int>*);
void gt_ggc_mx (modref_ref_node <tree_node*>* &);
void gt_pch_nx (modref_ref_node <int>* const&);
void gt_pch_nx (modref_ref_node <tree_node*>* const&);
void gt_pch_nx (modref_ref_node <int>* const&, gt_pointer_operator op,
		void *cookie);
void gt_pch_nx (modref_ref_node <tree_node*>* const&, gt_pointer_operator op,
		void *cookie);

#endif
