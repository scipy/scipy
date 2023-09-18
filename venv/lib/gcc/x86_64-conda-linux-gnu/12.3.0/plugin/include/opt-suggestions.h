/* Provide suggestions to handle misspelled options, and implement the
   --complete option for auto-completing options from a prefix.
   Copyright (C) 2016-2022 Free Software Foundation, Inc.

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

#ifndef GCC_OPT_PROPOSER_H
#define GCC_OPT_PROPOSER_H

/* Option proposer is class used by driver in order to provide hints
   for wrong options provided.  And it's used by --complete option that's
   intended to be invoked by BASH in order to provide better option
   completion support.  */

class option_proposer
{
public:
  /* Default constructor.  */
  option_proposer (): m_option_suggestions (NULL)
  {}

  /* Default destructor.  */
  ~option_proposer ();

  /* Helper function for driver::handle_unrecognized_options.

     Given an unrecognized option BAD_OPT (without the leading dash),
     locate the closest reasonable matching option (again, without the
     leading dash), or NULL.

     The returned string is owned by the option_proposer instance.  */
  const char *suggest_option (const char *bad_opt);

  /* Print on stdout a list of valid options that begin with OPTION_PREFIX,
     one per line, suitable for use by Bash completion.

     Implementation of the "-completion=" option.  */
  void suggest_completion (const char *option_prefix);

  /* Populate RESULTS with valid completions of options that begin
     with OPTION_PREFIX.  */
  void get_completions (const char *option_prefix, auto_string_vec &results);

private:
  /* Helper function for option_proposer::suggest_option.  Populate
     m_option_suggestions with candidate strings for misspelled options.
     The strings will be freed by the option_proposer's dtor.
     PREFIX is used for bash completion suggestions, otherwise
     it's set to NULL.  */
  void build_option_suggestions (const char *prefix);

private:
  /* Cache with all suggestions.  */
  auto_string_vec *m_option_suggestions;
};

#endif  /* GCC_OPT_PROPOSER_H */
