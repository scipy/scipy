/* Example illustrating how to use mpfr_can_round. */

/*
Copyright 2016-2023 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
https://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <stdio.h>
#include <mpfr.h>

int
main (void)
{
  mpfr_t x, y;
  mpfr_prec_t px = 53, py = 50;
  mpfr_rnd_t r1, r2;
  int ok;

  /* Given an approximation of Pi to px bits computed with rounding mode r1,
     we call mpfr_can_round() to see if we can deduced the correct rounding
     of Pi to py bits with rounding mode r2.
     The error is at most 1 = 2^0 ulp. This translates into err = prec(x). */
  mpfr_init2 (x, px);
  mpfr_init2 (y, py);
  for (r1 = 0; r1 < 4; r1++)
    {
      mpfr_const_pi (x, r1);
      printf ("r1=%s approx=", mpfr_print_rnd_mode (r1));
      mpfr_out_str (stdout, 2, 0, x, MPFR_RNDN);
      printf ("\n");
      for (r2 = 0; r2 < 4; r2++)
        {
          ok = mpfr_can_round (x, mpfr_get_prec (x), r1, r2, py);
          printf ("r2=%s ok=%d", mpfr_print_rnd_mode (r2), ok);
          if (ok)
            {
              mpfr_set (y, x, r2);
              printf ("   ");
              mpfr_out_str (stdout, 2, 0, y, MPFR_RNDN);
            }
          printf ("\n");
        }
    }
  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_free_cache ();
  return 0;
}
