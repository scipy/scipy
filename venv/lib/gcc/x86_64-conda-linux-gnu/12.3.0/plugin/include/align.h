/* Alignment-related classes.
   Copyright (C) 2018-2022 Free Software Foundation, Inc.

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

/* Align flags tuple with alignment in log form and with a maximum skip.  */

struct align_flags_tuple
{
  /* Values of the -falign-* flags: how much to align labels in code.
     log is "align to 2^log" (so 0 means no alignment).
     maxskip is the maximum allowed amount of padding to insert.  */
  int log;
  int maxskip;

  /* Normalize filled values so that maxskip is not bigger than 1 << log.  */
  void normalize ()
  {
    int n = (1 << log);
    if (maxskip > n)
      maxskip = n - 1;
  }

  /* Return original value of an alignment flag.  */
  int get_value ()
  {
    return maxskip + 1;
  }
};

/* Alignment flags is structure used as value of -align-* options.
   It's used in target-dependant code.  */

class align_flags
{
public:
  /* Default constructor.  */
  align_flags (int log0 = 0, int maxskip0 = 0, int log1 = 0, int maxskip1 = 0)
  {
    levels[0].log = log0;
    levels[0].maxskip = maxskip0;
    levels[1].log = log1;
    levels[1].maxskip = maxskip1;
    normalize ();
  }

  /* Normalize both components of align_flags.  */
  void normalize ()
  {
    for (unsigned i = 0; i < 2; i++)
      levels[i].normalize ();
  }

  /* Get alignment that is common bigger alignment of alignments F0 and F1.  */
  static align_flags max (const align_flags f0, const align_flags f1)
    {
      int log0 = MAX (f0.levels[0].log, f1.levels[0].log);
      int maxskip0 = MAX (f0.levels[0].maxskip, f1.levels[0].maxskip);
      int log1 = MAX (f0.levels[1].log, f1.levels[1].log);
      int maxskip1 = MAX (f0.levels[1].maxskip, f1.levels[1].maxskip);
      return align_flags (log0, maxskip0, log1, maxskip1);
    }

  align_flags_tuple levels[2];
};

/* Define maximum supported code alignment.  */
#define MAX_CODE_ALIGN 16
#define MAX_CODE_ALIGN_VALUE (1 << MAX_CODE_ALIGN)
