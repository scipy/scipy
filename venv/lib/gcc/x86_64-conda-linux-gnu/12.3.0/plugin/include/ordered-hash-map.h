/* A type-safe hash map that retains the insertion order of keys.
   Copyright (C) 2019-2022 Free Software Foundation, Inc.

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


#ifndef GCC_ORDERED_HASH_MAP_H
#define GCC_ORDERED_HASH_MAP_H

/* Notes:
   - The keys must be PODs, since vec<> uses assignment to populate slots
     without properly initializing them.
   - doesn't have GTY support.
   - supports removal, but retains order of original insertion.
     (Removal might be better handled by using a doubly-linked list
     of nodes, holding the values).  */

template<typename KeyId, typename Value,
	 typename Traits>
class ordered_hash_map
{
  typedef typename Traits::key_type Key;

public:
  ordered_hash_map () {}

  ordered_hash_map (const ordered_hash_map &other)
  : m_map (other.m_map),
    m_keys (other.m_keys.length ()),
    m_key_index (other.m_key_index)
  {
     unsigned i;
     Key key;
     FOR_EACH_VEC_ELT (other.m_keys, i, key)
       m_keys.quick_push (key);
  }

  /* If key K isn't already in the map add key K with value V to the map, and
     return false.  Otherwise set the value of the entry for key K to be V and
     return true.  */

  bool put (const Key &k, const Value &v)
  {
    bool existed = m_map.put (k, v);
    if (!existed)
      {
        bool key_present;
        int &slot = m_key_index.get_or_insert (k, &key_present);
        if (!key_present)
          {
             slot = m_keys.length ();
             m_keys.safe_push (k);
          }
      }
    return existed;
  }

  /* If the passed in key is in the map return its value otherwise NULL.  */

  Value *get (const Key &k)
  {
    return m_map.get (k);
  }

  /* Removing a key removes it from the map, but retains the insertion
     order.  */

  void remove (const Key &k)
  {
     m_map.remove (k);
  }

  size_t elements () const { return m_map.elements (); }

  class iterator
  {
  public:
    explicit iterator (const ordered_hash_map &map, unsigned idx) :
      m_ordered_hash_map (map), m_idx (idx) {}

    iterator &operator++ ()
    {
       /* Increment m_idx until we find a non-deleted element, or go beyond
	  the end.  */
       while (1)
	 {
	   ++m_idx;
	   if (valid_index_p ())
	     break;
	}
      return *this;
    }

    /* Can't use std::pair here, because GCC before 4.3 don't handle
       std::pair where template parameters are references well.
       See PR86739.  */
    struct reference_pair {
      const Key &first;
      Value &second;

      reference_pair (const Key &key, Value &value)
      : first (key), second (value) {}

      template <typename K, typename V>
      operator std::pair<K, V> () const { return std::pair<K, V> (first, second); }
    };

    reference_pair operator* ()
    {
      const Key &k = m_ordered_hash_map.m_keys[m_idx];
      Value *slot
        = const_cast<ordered_hash_map &> (m_ordered_hash_map).get (k);
      gcc_assert (slot);
      return reference_pair (k, *slot);
    }

    bool
    operator != (const iterator &other) const
    {
      return m_idx != other.m_idx;
    }

    /* Treat one-beyond-the-end as valid, for handling the "end" case.  */

    bool valid_index_p () const
    {
      if (m_idx > m_ordered_hash_map.m_keys.length ())
	return false;
      if (m_idx == m_ordered_hash_map.m_keys.length ())
	return true;
      const Key &k = m_ordered_hash_map.m_keys[m_idx];
      Value *slot
	= const_cast<ordered_hash_map &> (m_ordered_hash_map).get (k);
      return slot != NULL;
    }

    const ordered_hash_map &m_ordered_hash_map;
    unsigned m_idx;
  };

  /* Standard iterator retrieval methods.  */

  iterator begin () const
  {
    iterator i = iterator (*this, 0);
    while (!i.valid_index_p () && i != end ())
      ++i;
    return i;
  }
  iterator end () const { return iterator (*this, m_keys.length ()); }

private:
  /* The assignment operator is not yet implemented; prevent erroneous
     usage of unsafe compiler-generated one.  */
  void operator= (const ordered_hash_map &);

  /* The underlying map.  */
  hash_map<KeyId, Value, Traits> m_map;

  /* The ordering of the keys.  */
  auto_vec<Key> m_keys;

  /* For each key that's ever been in the map, its index within m_keys.  */
  hash_map<KeyId, int> m_key_index;
};

/* Two-argument form.  */

template<typename Key, typename Value,
	 typename Traits = simple_hashmap_traits<default_hash_traits<Key>,
						 Value> >
class ordered_hash_map;

#endif /* GCC_ORDERED_HASH_MAP_H */
