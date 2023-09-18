/* Handling of fnspec attribute specifiers
   Copyright (C) 2008-2022 Free Software Foundation, Inc.
   Contributed by Richard Guenther  <rguenther@suse.de>

   This file is part of GCC.

   GCC is free software; you can redistribute it and/or modify
   under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   GCC is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with GCC; see the file COPYING3.  If not see
   <http://www.gnu.org/licenses/>.  */

/* Parse string of attribute "fn spec".  This is an internal attribute
   describing side effects of a function as follows:

   character 0  specifies properties of return values as follows:
     '1'...'4'  specifies number of argument function returns (as in memset)
     'm'	specifies that returned value is noalias (as in malloc)
     '.'	specifies that nothing is known.
   character 1  specifies additional function properties
     ' '        specifies that nothing is known
     'p' or 'P' specifies that function is pure except for described side
		effects.
     'c' or 'C' specifies that function is const except for described side
		effects.
   The uppercase letter in addition specifies that function clobbers errno.

   character 2+2i specifies properties of argument number i as follows:
     'x' or 'X' specifies that parameter is unused.
     'r' or 'R' specifies that the memory pointed to by the parameter is only
		read and does not escape
     'o' or 'O' specifies that the memory pointed to by the parameter is only
		written and does not escape
     'w' or 'W' specifies that the memory pointed to by the parameter does not
		escape
     '1'....'9' specifies that the memory pointed to by the parameter is
		copied to memory pointed to by different parameter
		(as in memcpy).
     '.'	specifies that nothing is known.
   The uppercase letter in addition specifies that the memory pointed to
   by the parameter is not dereferenced.  For 'r' only read applies
   transitively to pointers read from the pointed-to memory.

   character 3+2i specifies additional properties of argument number i
   as follows:
     ' '        nothing is known
     't'	the size of value written/read corresponds to the size of
		of the pointed-to type of the argument type
     '1'...'9'  specifies the size of value written/read is bound by the
		specified argument
 */

#ifndef ATTR_FNSPEC_H
#define ATTR_FNSPEC_H

class attr_fnspec
{
private:
  /* fn spec attribute string.  */
  const char *str;
  /* length of the fn spec string.  */
  const unsigned len;
  /* Number of characters specifying return value.  */
  const unsigned int return_desc_size = 2;
  /* Number of characters specifying size.  */
  const unsigned int arg_desc_size = 2;

  /* Return start of specifier of arg i.  */
  unsigned int arg_idx (int i)
  {
    return return_desc_size + arg_desc_size * i;
  }

public:
  attr_fnspec (const char *str, unsigned len)
  : str (str), len (len)
  {
    if (flag_checking)
      verify ();
  }
  attr_fnspec (const char *str)
  : str (str), len (strlen (str))
  {
    if (flag_checking)
      verify ();
  }
  attr_fnspec (const_tree identifier)
  : str (TREE_STRING_POINTER (identifier)),
    len (TREE_STRING_LENGTH (identifier))
  {
    if (flag_checking)
      verify ();
  }
  attr_fnspec ()
  : str (NULL), len (0)
  {
  }

  /* Return true if fn spec is known.  */
  bool
  known_p ()
  {
    return len;
  }

  /* Return true if arg I is specified.  */
  bool
  arg_specified_p (unsigned int i)
  {
    return len >= arg_idx (i + 1);
  }

  /* True if the argument is not dereferenced recursively, thus only
     directly reachable memory is read or written.  */
  bool
  arg_direct_p (unsigned int i)
  {
    unsigned int idx = arg_idx (i);
    gcc_checking_assert (arg_specified_p (i));
    return str[idx] == 'R' || str[idx] == 'O'
	   || str[idx] == 'W' || (str[idx] >= '1' && str[idx] <= '9');
  }

  /* True if argument is used.  */
  bool
  arg_used_p (unsigned int i)
  {
    unsigned int idx = arg_idx (i);
    gcc_checking_assert (arg_specified_p (i));
    return str[idx] != 'x' && str[idx] != 'X';
  }

  /* True if memory reached by the argument is readonly (not clobbered).  */
  bool
  arg_readonly_p (unsigned int i)
  {
    unsigned int idx = arg_idx (i);
    gcc_checking_assert (arg_specified_p (i));
    return str[idx] == 'r' || str[idx] == 'R' || (str[idx] >= '1' && str[idx] <= '9');
  }

  /* True if memory reached by the argument is read (directly or indirectly)  */
  bool
  arg_maybe_read_p (unsigned int i)
  {
    unsigned int idx = arg_idx (i);
    gcc_checking_assert (arg_specified_p (i));
    return str[idx] != 'o' && str[idx] != 'O'
	   && str[idx] != 'x' && str[idx] != 'X';
  }

  /* True if memory reached by the argument is written.
     (directly or indirectly)  */
  bool
  arg_maybe_written_p (unsigned int i)
  {
    unsigned int idx = arg_idx (i);
    gcc_checking_assert (arg_specified_p (i));
    return str[idx] != 'r' && str[idx] != 'R'
	   && (str[idx] < '1' || str[idx] > '9')
	   && str[idx] != 'x' && str[idx] != 'X';
  }

  /* Return true if load of memory pointed to by argument I is bound
     by another argument.  In this case set ARG.  */
  bool
  arg_max_access_size_given_by_arg_p (unsigned int i, unsigned int *arg)
  {
    unsigned int idx = arg_idx (i);
    gcc_checking_assert (arg_specified_p (i));
    if (str[idx + 1] >= '1' && str[idx + 1] <= '9')
      {
	*arg = str[idx + 1] - '1';
	return true;
      }
    else
      return false;
  }

  /* Return true if the pointed-to type of the argument correspond to the
     size of the memory acccess.  */
  bool
  arg_access_size_given_by_type_p (unsigned int i)
  {
    unsigned int idx = arg_idx (i);
    gcc_checking_assert (arg_specified_p (i));
    return str[idx + 1] == 't';
  }

  /* Return true if memory pointer to by argument is copied to a memory
     pointed to by a different argument (as in memcpy).
     In this case set ARG.  */
  bool
  arg_copied_to_arg_p (unsigned int i, unsigned int *arg)
  {
    unsigned int idx = arg_idx (i);
    gcc_checking_assert (arg_specified_p (i));
    if (str[idx] < '1' || str[idx] > '9')
      return false;
    *arg = str[idx] - '1';
    return true;
  }


  /* True if the argument does not escape.  */
  bool
  arg_noescape_p (unsigned int i)
  {
    unsigned int idx = arg_idx (i);
    gcc_checking_assert (arg_specified_p (i));
    return str[idx] == 'w' || str[idx] == 'W'
	   || str[idx] == 'r' || str[idx] == 'R'
	   || str[idx] == 'o' || str[idx] == 'O';
  }

  /* Return true if function returns value of its parameter.  If ARG_NO is
     non-NULL return initialize it to the argument returned.  */
  bool
  returns_arg (unsigned int *arg_no)
  {
    if (str[0] >= '1' && str[0] <= '4')
      {
	if (arg_no)
	  *arg_no = str[0] - '1';
	return true;
      }
    return false;
  }

  /* Nonzero if the return value does not alias with anything.  Functions
     with the malloc attribute have this set on their return value.  */
  bool
  returns_noalias_p ()
  {
    return str[0] == 'm';
  }

  /* Return true if all memory read by the function is specified by fnspec.  */
  bool
  global_memory_read_p ()
  {
    return str[1] != 'c' && str[1] != 'C';
  }

  /* Return true if all memory written by the function
     is specified by fnspec.  */
  bool
  global_memory_written_p ()
  {
    return str[1] != 'c' && str[1] != 'C' && str[1] != 'p' && str[1] != 'P';
  }

  bool
  errno_maybe_written_p ()
  {
    return str[1] == 'C' || str[1] == 'P';
  }

  /* Return EAF flags for arg I.  */
  int
  arg_eaf_flags (unsigned int i)
  {
    int flags = 0;

    if (!arg_specified_p (i))
      ;
    else if (!arg_used_p (i))
      flags = EAF_UNUSED;
    else
      {
	if (arg_direct_p (i))
	  flags |= EAF_NO_INDIRECT_READ | EAF_NO_INDIRECT_ESCAPE
		   | EAF_NOT_RETURNED_INDIRECTLY | EAF_NO_INDIRECT_CLOBBER;
	if (arg_noescape_p (i))
	  flags |= EAF_NO_DIRECT_ESCAPE | EAF_NO_INDIRECT_ESCAPE;
	if (arg_readonly_p (i))
	  flags |= EAF_NO_DIRECT_CLOBBER | EAF_NO_INDIRECT_CLOBBER;
      }
    return flags;
  }

  /* Check validity of the string.  */
  void verify ();

  /* Return the fnspec string.  */
  const char *
  get_str ()
  {
    return str;
  }
};

extern attr_fnspec gimple_call_fnspec (const gcall *stmt);
extern attr_fnspec builtin_fnspec (tree);

#endif /* ATTR_FNSPEC_H  */
