/* Rich information on why an optimization wasn't possible.
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

#ifndef GCC_OPT_PROBLEM_H
#define GCC_OPT_PROBLEM_H

#include "diagnostic-core.h" /* for ATTRIBUTE_GCC_DIAG.  */
#include "optinfo.h" /* for optinfo.  */

/* This header declares a family of wrapper classes for tracking a
   success/failure value, while optionally supporting propagating an
   opt_problem * describing any failure back up the call stack.

   For instance, at the deepest point of the callstack where the failure
   happens, rather than:

     if (!check_something ())
       {
         if (dump_enabled_p ())
           dump_printf_loc (MSG_MISSED_OPTIMIZATION, vect_location,
                            "foo is unsupported.\n");
         return false;
       }
     // [...more checks...]

     // All checks passed:
     return true;

   we can capture the cause of the failure via:

     if (!check_something ())
       return opt_result::failure_at (stmt, "foo is unsupported");
     // [...more checks...]

     // All checks passed:
     return opt_result::success ();

   which effectively returns true or false, whilst recording any problem.

   opt_result::success and opt_result::failure return opt_result values
   which "looks like" true/false respectively, via operator bool().
   If dump_enabled_p, then opt_result::failure also creates an opt_problem *,
   capturing the pertinent data (here, "foo is unsupported " and "stmt").
   If dumps are disabled, then opt_problem instances aren't
   created, and it's equivalent to just returning a bool.

   The opt_problem can be propagated via opt_result values back up
   the call stack to where it makes most sense to the user.
   For instance, rather than:

     bool ok = try_something_that_might_fail ();
     if (!ok)
       {
         if (dump_enabled_p ())
           dump_printf_loc (MSG_MISSED_OPTIMIZATION, vect_location,
                            "some message.\n");
         return false;
       }

   we can replace the bool with an opt_result, so if dump_enabled_p, we
   assume that if try_something_that_might_fail, an opt_problem * will be
   created, and we can propagate it up the call chain:

     opt_result ok = try_something_that_might_fail ();
     if (!ok)
       {
         if (dump_enabled_p ())
           dump_printf_loc (MSG_MISSED_OPTIMIZATION, vect_location,
                            "some message.\n");
         return ok; // propagating the opt_result
       }

   opt_result is an opt_wrapper<bool>, where opt_wrapper<T> is a base
   class for wrapping a T, optionally propagating an opt_problem in
   case of failure_at (when dumps are enabled).  Similarly,
   opt_pointer_wrapper<T> can be used to wrap pointer types (where non-NULL
   signifies success, NULL signifies failure).

   In all cases, opt_wrapper<T> acts as if the opt_problem were one of its
   fields, but the opt_problem is actually stored in a global, so that when
   compiled, an opt_wrapper<T> is effectively just a T, so that we're
   still just passing e.g. a bool around; the opt_wrapper<T> classes
   simply provide type-checking and an API to ensure that we provide
   error-messages deep in the callstack at the places where problems
   occur, and that we propagate them.  This also avoids having
   to manage the ownership of the opt_problem instances.

   Using opt_result and opt_wrapper<T> documents the intent of the code
   for the places where we represent success values, and allows the C++ type
   system to track where the deepest points in the callstack are where we
   need to emit the failure messages from.  */

/* A bundle of information about why an optimization failed (e.g.
   vectorization), and the location in both the user's code and
   in GCC itself where the problem occurred.

   Instances are created by static member functions in opt_wrapper
   subclasses, such as opt_result::failure.

   Instances are only created when dump_enabled_p ().  */

class opt_problem
{
 public:
  static opt_problem *get_singleton () { return s_the_problem; }

  opt_problem (const dump_location_t &loc,
	       const char *fmt, va_list *ap)
    ATTRIBUTE_GCC_DUMP_PRINTF (3, 0);

  const dump_location_t &
  get_dump_location () const { return m_optinfo.get_dump_location (); }

  const optinfo & get_optinfo () const { return m_optinfo; }

  void emit_and_clear ();

 private:
  optinfo m_optinfo;

  static opt_problem *s_the_problem;
};

/* A base class for wrapper classes that track a success/failure value, while
   optionally supporting propagating an opt_problem * describing any
   failure back up the call stack.  */

template <typename T>
class opt_wrapper
{
 public:
  typedef T wrapped_t;

  /* Be accessible as the wrapped type.  */
  operator wrapped_t () const { return m_result; }

  /* No public ctor.  */

  wrapped_t get_result () const { return m_result; }
  opt_problem *get_problem () const { return opt_problem::get_singleton (); }

 protected:
  opt_wrapper (wrapped_t result, opt_problem */*problem*/)
  : m_result (result)
  {
    /* "problem" is ignored: although it looks like a field, we
       actually just use the opt_problem singleton, so that
       opt_wrapper<T> in memory is just a T.  */
  }

 private:
  wrapped_t m_result;
};

/* Subclass of opt_wrapper<T> for bool, where
   - true signifies "success", and
   - false signifies "failure"
   whilst effectively propagating an opt_problem * describing any failure
   back up the call stack.  */

class opt_result : public opt_wrapper <bool>
{
 public:
  /* Generate a "success" value: a wrapper around "true".  */

  static opt_result success () { return opt_result (true, NULL); }

  /* Generate a "failure" value: a wrapper around "false", and,
     if dump_enabled_p, an opt_problem.  */

  static opt_result failure_at (const dump_location_t &loc,
				const char *fmt, ...)
	  ATTRIBUTE_GCC_DUMP_PRINTF (2, 3)
  {
    opt_problem *problem = NULL;
    if (dump_enabled_p ())
      {
	va_list ap;
	va_start (ap, fmt);
	problem = new opt_problem (loc, fmt, &ap);
	va_end (ap);
      }
    return opt_result (false, problem);
  }

  /* Given a failure wrapper of some other kind, make an opt_result failure
     object, for propagating the opt_problem up the call stack.  */

  template <typename S>
  static opt_result
  propagate_failure (opt_wrapper <S> other)
  {
    return opt_result (false, other.get_problem ());
  }

 private:
  /* Private ctor.  Instances should be created by the success and failure
     static member functions.  */
  opt_result (wrapped_t result, opt_problem *problem)
  : opt_wrapper <bool> (result, problem)
  {}
};

/* Subclass of opt_wrapper<T> where T is a pointer type, for tracking
   success/failure, where:
   - a non-NULL value signifies "success", and
   - a NULL value signifies "failure",
   whilst effectively propagating an opt_problem * describing any failure
   back up the call stack.  */

template <typename PtrType_t>
class opt_pointer_wrapper : public opt_wrapper <PtrType_t>
{
 public:
  typedef PtrType_t wrapped_pointer_t;

  /* Given a non-NULL pointer, make a success object wrapping it.  */

  static opt_pointer_wrapper <wrapped_pointer_t>
  success (wrapped_pointer_t ptr)
  {
    return opt_pointer_wrapper <wrapped_pointer_t> (ptr, NULL);
  }

  /* Make a NULL pointer failure object, with the given message
     (if dump_enabled_p).  */

  static opt_pointer_wrapper <wrapped_pointer_t>
  failure_at (const dump_location_t &loc,
	      const char *fmt, ...)
    ATTRIBUTE_GCC_DUMP_PRINTF (2, 3)
  {
    opt_problem *problem = NULL;
    if (dump_enabled_p ())
      {
	va_list ap;
	va_start (ap, fmt);
	problem = new opt_problem (loc, fmt, &ap);
	va_end (ap);
      }
    return opt_pointer_wrapper <wrapped_pointer_t> (NULL, problem);
  }

  /* Given a failure wrapper of some other kind, make a NULL pointer
     failure object, propagating the problem.  */

  template <typename S>
  static opt_pointer_wrapper <wrapped_pointer_t>
  propagate_failure (opt_wrapper <S> other)
  {
    return opt_pointer_wrapper <wrapped_pointer_t> (NULL,
						    other.get_problem ());
  }

  /* Support accessing the underlying pointer via ->.  */

  wrapped_pointer_t operator-> () const { return this->get_result (); }

 private:
  /* Private ctor.  Instances should be built using the static member
     functions "success" and "failure".  */
  opt_pointer_wrapper (wrapped_pointer_t result, opt_problem *problem)
  : opt_wrapper<PtrType_t> (result, problem)
  {}
};

/* A typedef for wrapping "tree" so that NULL_TREE can carry an
   opt_problem describing the failure (if dump_enabled_p).  */

typedef opt_pointer_wrapper<tree> opt_tree;

#endif /* #ifndef GCC_OPT_PROBLEM_H */
