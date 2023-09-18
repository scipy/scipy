// Multiplexer utilities
// Copyright (C) 2020-2022 Free Software Foundation, Inc.
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

#ifndef GCC_MUX_UTILS_H
#define GCC_MUX_UTILS_H 1

// A class that stores a choice "A or B", where A has type T1 * and B has
// type T2 *.  Both T1 and T2 must have an alignment greater than 1, since
// the low bit is used to identify B over A.  T1 and T2 can be the same.
//
// A can be a null pointer but B cannot.
//
// Barring the requirement that B must be nonnull, using the class is
// equivalent to using:
//
//     union { T1 *A; T2 *B; };
//
// and having a separate tag bit to indicate which alternative is active.
// However, using this class can have two advantages over a union:
//
// - It avoides the need to find somewhere to store the tag bit.
//
// - The compiler is aware that B cannot be null, which can make checks
//   of the form:
//
//       if (auto *B = mux.dyn_cast<T2 *> ())
//
//   more efficient.  With a union-based representation, the dyn_cast
//   check could fail either because MUX is an A or because MUX is a
//   null B, both of which require a run-time test.  With a pointer_mux,
//   only a check for MUX being A is needed.
template<typename T1, typename T2 = T1>
class pointer_mux
{
public:
  // Return an A pointer with the given value.
  static pointer_mux first (T1 *);

  // Return a B pointer with the given (nonnull) value.
  static pointer_mux second (T2 *);

  pointer_mux () = default;

  // Create a null A pointer.
  pointer_mux (std::nullptr_t) : m_ptr (nullptr) {}

  // Create an A or B pointer with the given value.  This is only valid
  // if T1 and T2 are distinct and if T can be resolved to exactly one
  // of them.
  template<typename T,
	   typename Enable = typename
	     std::enable_if<std::is_convertible<T *, T1 *>::value
			    != std::is_convertible<T *, T2 *>::value>::type>
  pointer_mux (T *ptr);

  // Return true unless the pointer is a null A pointer.
  explicit operator bool () const { return m_ptr; }

  // Assign A and B pointers respectively.
  void set_first (T1 *ptr) { *this = first (ptr); }
  void set_second (T2 *ptr) { *this = second (ptr); }

  // Return true if the pointer is an A pointer.
  bool is_first () const { return !(uintptr_t (m_ptr) & 1); }

  // Return true if the pointer is a B pointer.
  bool is_second () const { return uintptr_t (m_ptr) & 1; }

  // Return the contents of the pointer, given that it is known to be
  // an A pointer.
  T1 *known_first () const { return reinterpret_cast<T1 *> (m_ptr); }

  // Return the contents of the pointer, given that it is known to be
  // a B pointer.
  T2 *known_second () const { return reinterpret_cast<T2 *> (m_ptr - 1); }

  // If the pointer is an A pointer, return its contents, otherwise
  // return null.  Thus a null return can mean that the pointer is
  // either a null A pointer or a B pointer.
  //
  // If all A pointers are nonnull, it is more efficient to use:
  //
  //    if (ptr.is_first ())
  //      ...use ptr.known_first ()...
  //
  // over:
  //
  //    if (T1 *a = ptr.first_or_null ())
  //      ...use a...
  T1 *first_or_null () const;

  // If the pointer is a B pointer, return its contents, otherwise
  // return null.  Using:
  //
  //    if (T1 *b = ptr.second_or_null ())
  //      ...use b...
  //
  // should be at least as efficient as:
  //
  //    if (ptr.is_second ())
  //      ...use ptr.known_second ()...
  T2 *second_or_null () const;

  // Return true if the pointer is a T.
  //
  // This is only valid if T1 and T2 are distinct and if T can be
  // resolved to exactly one of them.  The condition is checked using
  // a static assertion rather than SFINAE because it gives a clearer
  // error message.
  template<typename T>
  bool is_a () const;

  // Assert that the pointer is a T and return it as such.  See is_a
  // for the restrictions on T.
  template<typename T>
  T as_a () const;

  // If the pointer is a T, return it as such, otherwise return null.
  // See is_a for the restrictions on T.
  template<typename T>
  T dyn_cast () const;

private:
  pointer_mux (char *ptr) : m_ptr (ptr) {}

  // Points to the first byte of an object for A pointers or the second
  // byte of an object for B pointers.  Using a pointer rather than a
  // uintptr_t tells the compiler that second () can never return null,
  // and that second_or_null () is only null if is_first ().
  char *m_ptr;
};

template<typename T1, typename T2>
inline pointer_mux<T1, T2>
pointer_mux<T1, T2>::first (T1 *ptr)
{
  gcc_checking_assert (!(uintptr_t (ptr) & 1));
  return reinterpret_cast<char *> (ptr);
}

template<typename T1, typename T2>
inline pointer_mux<T1, T2>
pointer_mux<T1, T2>::second (T2 *ptr)
{
  gcc_checking_assert (ptr && !(uintptr_t (ptr) & 1));
  return reinterpret_cast<char *> (ptr) + 1;
}

template<typename T1, typename T2>
template<typename T, typename Enable>
inline pointer_mux<T1, T2>::pointer_mux (T *ptr)
  : m_ptr (reinterpret_cast<char *> (ptr))
{
  if (std::is_convertible<T *, T2 *>::value)
    {
      gcc_checking_assert (m_ptr);
      m_ptr += 1;
    }
}

template<typename T1, typename T2>
inline T1 *
pointer_mux<T1, T2>::first_or_null () const
{
  return is_first () ? known_first () : nullptr;
}

template<typename T1, typename T2>
inline T2 *
pointer_mux<T1, T2>::second_or_null () const
{
  // Micro optimization that's effective as of GCC 11: compute the value
  // of the second pointer as an integer and test that, so that the integer
  // result can be reused as the pointer and so that all computation can
  // happen before a branch on null.  This reduces the number of branches
  // needed for loops.
  return (uintptr_t (m_ptr) - 1) & 1 ? nullptr : known_second ();
}

template<typename T1, typename T2>
template<typename T>
inline bool
pointer_mux<T1, T2>::is_a () const
{
  static_assert (std::is_convertible<T1 *, T>::value
		 != std::is_convertible<T2 *, T>::value,
		 "Ambiguous pointer type");
  if (std::is_convertible<T2 *, T>::value)
    return is_second ();
  else
    return is_first ();
}

template<typename T1, typename T2>
template<typename T>
inline T
pointer_mux<T1, T2>::as_a () const
{
  static_assert (std::is_convertible<T1 *, T>::value
		 != std::is_convertible<T2 *, T>::value,
		 "Ambiguous pointer type");
  if (std::is_convertible<T2 *, T>::value)
    {
      gcc_checking_assert (is_second ());
      return reinterpret_cast<T> (m_ptr - 1);
    }
  else
    {
      gcc_checking_assert (is_first ());
      return reinterpret_cast<T> (m_ptr);
    }
}

template<typename T1, typename T2>
template<typename T>
inline T
pointer_mux<T1, T2>::dyn_cast () const
{
  static_assert (std::is_convertible<T1 *, T>::value
		 != std::is_convertible<T2 *, T>::value,
		 "Ambiguous pointer type");
  if (std::is_convertible<T2 *, T>::value)
    {
      if (is_second ())
	return reinterpret_cast<T> (m_ptr - 1);
    }
  else
    {
      if (is_first ())
	return reinterpret_cast<T> (m_ptr);
    }
  return nullptr;
}

#endif
