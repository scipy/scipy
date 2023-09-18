/* Multithreading test to detect scaling issues with MPFR.

Define:
  * the function F;
  * the precision PREC;
  * the value V as an expression that will have the type double
    (it may depend on the thread number i).

Example:
  gcc threads.c -lmpfr -lgmp -lpthread -DF=mpfr_sin -DPREC=200 -DV=100

Copyright 2018-2023 Free Software Foundation, Inc.
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
#include <stdlib.h>
#include <pthread.h>

#include <mpfr.h>

#define MAX_THREADS 256

static int m;

static void *start_routine (void *arg)
{
  mpfr_t x, y;
  int i = *(int *) arg, j;

  (void) i;  /* avoid a warning if i is not used by V */

  mpfr_inits2 (PREC, x, y, (mpfr_ptr) 0);
  mpfr_set_d (x, (V), MPFR_RNDN);

  for (j = 0; j < m; j++)
    F (y, x, MPFR_RNDN);

  mpfr_clears (x, y, (mpfr_ptr) 0);
  pthread_exit (NULL);
}

int main (int argc, char *argv[])
{
  int i, n;
  pthread_t tid[MAX_THREADS];

  if (argc != 3 ||
      (m = atoi (argv[1]), m < 1) ||
      (n = atoi (argv[2]), n < 1 || n > MAX_THREADS))
    {
      fprintf (stderr, "Usage: %s <#iterations> <#threads>\n", argv[0]);
      exit (1);
    }

  printf ("%d iteration(s), %d thread(s).\n", m, n);

  for (i = 0; i < n; i++)
    if (pthread_create (&tid[i], NULL, start_routine, &i) != 0)
      {
        fprintf (stderr, "%s: failed to create thread %d\n", argv[0], i);
        exit (1);
      }

  for (i = 0; i < n; i++)
    if (pthread_join (tid[i], NULL) != 0)
      {
        fprintf (stderr, "%s: failed to join thread %d\n", argv[0], i);
        exit (1);
      }

  return 0;
}
