/* A type-safe hash map.
   Copyright (C) 2014-2022 Free Software Foundation, Inc.

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


#ifndef hash_map_h
#define hash_map_h

/* Class hash_map is a hash-value based container mapping objects of
   KeyId type to those of the Value type.
   Both KeyId and Value may be non-trivial (non-POD) types provided
   a suitabe Traits class.  A few default Traits specializations are
   provided for basic types such as integers, pointers, and std::pair.
   Inserted elements are value-initialized either to zero for POD types
   or by invoking their default ctor.  Removed elements are destroyed
   by invoking their dtor.   On hash_map destruction all elements are
   removed.  Objects of hash_map type are copy-constructible but not
   assignable.  */

const size_t default_hash_map_size = 13;
template<typename KeyId, typename Value,
	 typename Traits /* = simple_hashmap_traits<default_hash_traits<Key>,
			                            Value> */>
class GTY((user)) hash_map
{
  typedef typename Traits::key_type Key;
  struct hash_entry
  {
    Key m_key;
    Value m_value;

    typedef hash_entry value_type;
    typedef Key compare_type;

    static hashval_t hash (const hash_entry &e)
      {
       	return Traits::hash (e.m_key);
      }

    static bool equal (const hash_entry &a, const Key &b)
       	{
	  return Traits::equal_keys (a.m_key, b);
       	}

    static void remove (hash_entry &e) { Traits::remove (e); }

    static void mark_deleted (hash_entry &e) { Traits::mark_deleted (e); }

    static bool is_deleted (const hash_entry &e)
      {
       	return Traits::is_deleted (e);
      }

    static const bool empty_zero_p = Traits::empty_zero_p;
    static void mark_empty (hash_entry &e) { Traits::mark_empty (e); }
    static bool is_empty (const hash_entry &e) { return Traits::is_empty (e); }

    static void ggc_mx (hash_entry &e)
      {
	gt_ggc_mx (e.m_key);
	gt_ggc_mx (e.m_value);
      }

    static void ggc_maybe_mx (hash_entry &e)
      {
	if (Traits::maybe_mx)
	  ggc_mx (e);
      }

    static void pch_nx (hash_entry &e)
      {
	gt_pch_nx (e.m_key);
	gt_pch_nx (e.m_value);
      }

    static void pch_nx (hash_entry &e, gt_pointer_operator op, void *c)
      {
	pch_nx_helper (e.m_key, op, c);
	pch_nx_helper (e.m_value, op, c);
      }

    static int keep_cache_entry (hash_entry &e)
      {
	return ggc_marked_p (e.m_key);
      }

  private:
    template<typename T>
    static void
      pch_nx_helper (T &x, gt_pointer_operator op, void *cookie)
	{
	  gt_pch_nx (&x, op, cookie);
	}

    template<typename T>
      static void
      pch_nx_helper (T *&x, gt_pointer_operator op, void *cookie)
	{
	  op (&x, NULL, cookie);
	}

    /* The overloads below should match those in ggc.h.  */
#define DEFINE_PCH_HELPER(T)			\
    static void pch_nx_helper (T, gt_pointer_operator, void *) { }

    DEFINE_PCH_HELPER (bool);
    DEFINE_PCH_HELPER (char);
    DEFINE_PCH_HELPER (signed char);
    DEFINE_PCH_HELPER (unsigned char);
    DEFINE_PCH_HELPER (short);
    DEFINE_PCH_HELPER (unsigned short);
    DEFINE_PCH_HELPER (int);
    DEFINE_PCH_HELPER (unsigned int);
    DEFINE_PCH_HELPER (long);
    DEFINE_PCH_HELPER (unsigned long);
    DEFINE_PCH_HELPER (long long);
    DEFINE_PCH_HELPER (unsigned long long);

#undef DEFINE_PCH_HELPER
  };

public:
  explicit hash_map (size_t n = default_hash_map_size, bool ggc = false,
		     bool sanitize_eq_and_hash = true,
		     bool gather_mem_stats = GATHER_STATISTICS
		     CXX_MEM_STAT_INFO)
    : m_table (n, ggc, sanitize_eq_and_hash, gather_mem_stats,
	       HASH_MAP_ORIGIN PASS_MEM_STAT)
  {
  }

  explicit hash_map (const hash_map &h, bool ggc = false,
		     bool sanitize_eq_and_hash = true,
		     bool gather_mem_stats = GATHER_STATISTICS
		     CXX_MEM_STAT_INFO)
    : m_table (h.m_table, ggc, sanitize_eq_and_hash, gather_mem_stats,
	       HASH_MAP_ORIGIN PASS_MEM_STAT) {}

  /* Create a hash_map in ggc memory.  */
  static hash_map *create_ggc (size_t size = default_hash_map_size,
			       bool gather_mem_stats = GATHER_STATISTICS
			       CXX_MEM_STAT_INFO)
    {
      hash_map *map = ggc_alloc<hash_map> ();
      new (map) hash_map (size, true, true, gather_mem_stats PASS_MEM_STAT);
      return map;
    }

  /* If key k isn't already in the map add key k with value v to the map, and
     return false.  Otherwise set the value of the entry for key k to be v and
     return true.  */

  bool put (const Key &k, const Value &v)
    {
      hash_entry *e = m_table.find_slot_with_hash (k, Traits::hash (k),
						   INSERT);
      bool ins = hash_entry::is_empty (*e);
      if (ins)
	{
	  e->m_key = k;
	  new ((void *) &e->m_value) Value (v);
	}
      else
	e->m_value = v;

      return !ins;
    }

  /* If the passed in key is in the map return pointer to its value
     otherwise NULL.  */

  Value *get (const Key &k)
    {
      hash_entry &e = m_table.find_with_hash (k, Traits::hash (k));
      return Traits::is_empty (e) ? NULL : &e.m_value;
    }

  /* Return a reference to the value for the passed in key, creating the entry
     if it doesn't already exist.  If existed is not NULL then it is set to
     false if the key was not previously in the map, and true otherwise.  */

  Value &get_or_insert (const Key &k, bool *existed = NULL)
    {
      hash_entry *e = m_table.find_slot_with_hash (k, Traits::hash (k),
						   INSERT);
      bool ins = Traits::is_empty (*e);
      if (ins)
	{
	  e->m_key = k;
	  new ((void *)&e->m_value) Value ();
	}

      if (existed != NULL)
	*existed = !ins;

      return e->m_value;
    }

  void remove (const Key &k)
    {
      m_table.remove_elt_with_hash (k, Traits::hash (k));
    }

  /* Call the call back on each pair of key and value with the passed in
     arg until either the call back returns false or all pairs have been seen.
     The traversal is unordered.  */

  template<typename Arg, bool (*f)(const typename Traits::key_type &,
				   const Value &, Arg)>
  void traverse (Arg a) const
    {
      for (typename hash_table<hash_entry>::iterator iter = m_table.begin ();
	   iter != m_table.end (); ++iter)
	if (!f ((*iter).m_key, (*iter).m_value, a))
	  break;
    }

  template<typename Arg, bool (*f)(const typename Traits::key_type &,
				   Value *, Arg)>
  void traverse (Arg a) const
    {
      for (typename hash_table<hash_entry>::iterator iter = m_table.begin ();
	   iter != m_table.end (); ++iter)
	if (!f ((*iter).m_key, &(*iter).m_value, a))
	  break;
    }

  size_t elements () const { return m_table.elements (); }

  void empty () { m_table.empty(); }

  /* Return true when there are no elements in this hash map.  */
  bool is_empty () const { return m_table.is_empty (); }

  class iterator
  {
  public:
    explicit iterator (const typename hash_table<hash_entry>::iterator &iter) :
      m_iter (iter) {}

    iterator &operator++ ()
    {
      ++m_iter;
      return *this;
    }

    /* Can't use std::pair here, because GCC before 4.3 don't handle
       std::pair where template parameters are references well.
       See PR86739.  */
    class reference_pair {
    public:
      const Key &first;
      Value &second;

      reference_pair (const Key &key, Value &value) : first (key), second (value) {}

      template <typename K, typename V>
      operator std::pair<K, V> () const { return std::pair<K, V> (first, second); }
    };

    reference_pair operator* ()
    {
      hash_entry &e = *m_iter;
      return reference_pair (e.m_key, e.m_value);
    }

    bool operator== (const iterator &other) const
    {
      return m_iter == other.m_iter;
    }

    bool operator != (const iterator &other) const
    {
      return m_iter != other.m_iter;
    }

  private:
    typename hash_table<hash_entry>::iterator m_iter;
  };

  /* Standard iterator retrieval methods.  */

  iterator  begin () const { return iterator (m_table.begin ()); }
  iterator end () const { return iterator (m_table.end ()); }

private:

  template<typename T, typename U, typename V> friend void gt_ggc_mx (hash_map<T, U, V> *);
  template<typename T, typename U, typename V> friend void gt_pch_nx (hash_map<T, U, V> *);
  template<typename T, typename U, typename V> friend void gt_pch_nx (hash_map<T, U, V> *, gt_pointer_operator, void *);
  template<typename T, typename U, typename V> friend void gt_cleare_cache (hash_map<T, U, V> *);

  hash_table<hash_entry> m_table;
};

/* ggc marking routines.  */

template<typename K, typename V, typename H>
static inline void
gt_ggc_mx (hash_map<K, V, H> *h)
{
  gt_ggc_mx (&h->m_table);
}

template<typename K, typename V, typename H>
static inline void
gt_pch_nx (hash_map<K, V, H> *h)
{
  gt_pch_nx (&h->m_table);
}

template<typename K, typename V, typename H>
static inline void
gt_cleare_cache (hash_map<K, V, H> *h)
{
  if (h)
    gt_cleare_cache (&h->m_table);
}

template<typename K, typename V, typename H>
static inline void
gt_pch_nx (hash_map<K, V, H> *h, gt_pointer_operator op, void *cookie)
{
  op (&h->m_table.m_entries, NULL, cookie);
}

enum hm_alloc { hm_heap = false, hm_ggc = true };
template<bool ggc, typename K, typename V, typename H>
inline hash_map<K,V,H> *
hash_map_maybe_create (hash_map<K,V,H> *&h,
		       size_t size = default_hash_map_size)
{
  if (!h)
    {
      if (ggc)
	h = hash_map<K,V,H>::create_ggc (size);
      else
	h = new hash_map<K,V,H> (size);
    }
  return h;
}

/* Like h->get, but handles null h.  */
template<typename K, typename V, typename H>
inline V*
hash_map_safe_get (hash_map<K,V,H> *h, const K& k)
{
  return h ? h->get (k) : NULL;
}

/* Like h->get, but handles null h.  */
template<bool ggc, typename K, typename V, typename H>
inline V&
hash_map_safe_get_or_insert (hash_map<K,V,H> *&h, const K& k, bool *e = NULL,
			     size_t size = default_hash_map_size)
{
  return hash_map_maybe_create<ggc> (h, size)->get_or_insert (k, e);
}

/* Like h->put, but handles null h.  */
template<bool ggc, typename K, typename V, typename H>
inline bool
hash_map_safe_put (hash_map<K,V,H> *&h, const K& k, const V& v,
		   size_t size = default_hash_map_size)
{
  return hash_map_maybe_create<ggc> (h, size)->put (k, v);
}

#endif
