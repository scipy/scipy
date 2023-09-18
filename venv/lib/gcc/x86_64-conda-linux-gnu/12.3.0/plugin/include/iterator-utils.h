// Iterator-related utilities.
// Copyright (C) 2002-2022 Free Software Foundation, Inc.
//
// This file is part of GCC.
//
// GCC is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 3, or (at your option) any later
// version.
//
// GCC is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// along with GCC; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.

#ifndef GCC_ITERATOR_UTILS_H
#define GCC_ITERATOR_UTILS_H 1

// A half-open [begin, end) range of iterators.
template<typename T>
struct iterator_range
{
public:
  using const_iterator = T;

  iterator_range () = default;
  iterator_range (const T &begin, const T &end)
    : m_begin (begin), m_end (end) {}

  T begin () const { return m_begin; }
  T end () const { return m_end; }

  explicit operator bool () const { return m_begin != m_end; }

private:
  T m_begin;
  T m_end;
};

// Provide an iterator like BaseIT, except that it yields values of type T,
// which is derived from the type that BaseIT normally yields.
//
// The class doesn't inherit from BaseIT for two reasons:
// - using inheritance would stop the class working with plain pointers
// - not using inheritance increases type-safety for writable iterators
//
// Constructing this class from a BaseIT involves an assertion that all
// contents really do have type T.  The constructor is therefore explicit.
template<typename T, typename BaseIT>
class derived_iterator
{
public:
  using value_type = T;

  derived_iterator () = default;

  template<typename... Ts>
  explicit derived_iterator (Ts... args)
    : m_base (std::forward<Ts> (args)...) {}

  derived_iterator &operator++ () { ++m_base; return *this; }
  derived_iterator operator++ (int);

  T operator* () const { return static_cast<T> (*m_base); }
  T *operator-> () const { return static_cast<T *> (m_base.operator-> ()); }

  bool operator== (const derived_iterator &other) const;
  bool operator!= (const derived_iterator &other) const;

protected:
  BaseIT m_base;
};

template<typename T, typename BaseIT>
inline derived_iterator<T, BaseIT>
derived_iterator<T, BaseIT>::operator++ (int)
{
  derived_iterator ret = *this;
  ++m_base;
  return ret;
}

template<typename T, typename BaseIT>
inline bool
derived_iterator<T, BaseIT>::operator== (const derived_iterator &other) const
{
  return m_base == other.m_base;
}

template<typename T, typename BaseIT>
inline bool
derived_iterator<T, BaseIT>::operator!= (const derived_iterator &other) const
{
  return m_base != other.m_base;
}

// Provide a constant view of a BaseCT in which every value is known to
// have type T, which is derived from the type that BaseCT normally presents.
//
// Constructing this class from a BaseCT involves an assertion that all
// contents really do have type T.  The constructor is therefore explicit.
template<typename T, typename BaseCT>
class const_derived_container : public BaseCT
{
  using base_const_iterator = typename BaseCT::const_iterator;

public:
  using value_type = T;
  using const_iterator = derived_iterator<T, base_const_iterator>;

  const_derived_container () = default;

  template<typename... Ts>
  explicit const_derived_container (Ts... args)
    : BaseCT (std::forward<Ts> (args)...) {}

  const_iterator begin () const { return const_iterator (BaseCT::begin ()); }
  const_iterator end () const { return const_iterator (BaseCT::end ()); }

  T front () const { return static_cast<T> (BaseCT::front ()); }
  T back () const { return static_cast<T> (BaseCT::back ()); }
  T operator[] (unsigned int i) const;
};

template<typename T, typename BaseCT>
inline T
const_derived_container<T, BaseCT>::operator[] (unsigned int i) const
{
  return static_cast<T> (BaseCT::operator[] (i));
}

// A base class for iterators whose contents consist of a StoredT and that
// when dereferenced yield those StoredT contents as a T.  Derived classes
// should implement at least operator++ or operator--.
template<typename T, typename StoredT = T>
class wrapper_iterator
{
public:
  using value_type = T;

  wrapper_iterator () = default;

  template<typename... Ts>
  wrapper_iterator (Ts... args) : m_contents (std::forward<Ts> (args)...) {}

  T operator* () const { return static_cast<T> (m_contents); }
  bool operator== (const wrapper_iterator &) const;
  bool operator!= (const wrapper_iterator &) const;

protected:
  StoredT m_contents;
};

template<typename T, typename StoredT>
inline bool
wrapper_iterator<T, StoredT>::operator== (const wrapper_iterator &other) const
{
  return m_contents == other.m_contents;
}

template<typename T, typename StoredT>
inline bool
wrapper_iterator<T, StoredT>::operator!= (const wrapper_iterator &other) const
{
  return m_contents != other.m_contents;
}

// A forward iterator for a linked list whose nodes are referenced using
// type T.  Given a node "T N", the next element is given by (N->*Next) ().
template<typename T, T *(T::*Next) () const>
class list_iterator : public wrapper_iterator<T *>
{
private:
  using parent = wrapper_iterator<T *>;

public:
  using parent::parent;
  list_iterator &operator++ ();
  list_iterator operator++ (int);
};

template<typename T, T *(T::*Next) () const>
inline list_iterator<T, Next> &
list_iterator<T, Next>::operator++ ()
{
  this->m_contents = (this->m_contents->*Next) ();
  return *this;
}

template<typename T, T *(T::*Next) () const>
inline list_iterator<T, Next>
list_iterator<T, Next>::operator++ (int)
{
  list_iterator ret = *this;
  this->m_contents = (this->m_contents->*Next) ();
  return ret;
}

#endif
