/**
  This file contains two algorithms written in C to generate random two-way
  tables. The algorithms rcont1 and rcont2 originate from Boyett and Patefield,
  respectively. For more information, see the docs for each function.

  If you wonder about the spelling, rcont is short for random contingency
  table. Random contingency table is a bit of a misnomer. The tables generated
  by these algorithms have no contingency/association between the two
  variables.

  Author: Hans Dembinski
*/
#include <math.h>
#include <stdbool.h>

#include "logfactorial.h"
#include "_rcont.h"

// helper function to access a 1D array like a C-style 2D array
tab_t *ptr(tab_t *m, int nr, int nc, int ir, int ic)
{
  return m + nc * ir + ic;
}

/*
  Call this once to initialize workspace for rcont1.

  The work space must have size N, where N is the total number of entries.
*/
void rcont1_init(tab_t *work, int nc, const tab_t *c)
{
  for (int i = 0; i < nc; ++i)
  {
    tab_t ci = c[i];
    while (ci--)
      *work++ = i;
  }
}

/*
  Generate random two-way table with given marginal totals.

  Boyett's shuffling algorithm adapted from AS 144 Appl. Statist. (1979)
  329-332. The algorithm has O(N) complexity in space and time for
  a table with N entries in total. The algorithm performs poorly for large N,
  but is insensitive to the number K of table cells.

  This function uses a work space of size N which must be preallocated and
  initialized with rcont1_init.
*/
void rcont1(tab_t *table, int nr, const tab_t *r, int nc, const tab_t *c,
            const tab_t ntot, tab_t *work, bitgen_t *rstate)
{
  // nothing to do
  if (ntot == 0)
    return;

  // shuffle work with Knuth's algorithm
  for (tab_t i = ntot - 1; i > 0; --i)
  {
    tab_t j = random_interval(rstate, i);
    tab_t tmp = work[j];
    work[j] = work[i];
    work[i] = tmp;
  }

  // clear table
  for (int i = 0, nrc = (nr * nc); i < nrc; ++i)
    table[i] = 0;

  // fill table
  for (int ir = 0; ir < nr; ++ir)
  {
    tab_t ri = r[ir];
    while (ri--)
      *ptr(table, nr, nc, ir, *work++) += 1;
  }
}

/*
  Generate random two-way table with given marginal totals.

  Patefield's algorithm adapted from AS 159 Appl. Statist. (1981) 91-97. This
  algorithm has O(K log(N)) complexity in time for a table with K cells and N
  entries in total. It requires only a small constant stack space.

  The original FORTRAN code was hand-translated to C. Changes to the original:

  - The computation of a look-up table of log-factorials was replaced with
    logfactorial function from numpy (which does something similar).
  - The original implementation allocated a column vector JWORK, but this is
    not necessary. The vector can be folded into the last column of the output
    table.
  - The function uses Numpy's random number generator and distribution library.
  - The algorithm now handles zero entries in row or column vector. When a
    zero is encountered, the output table is filled with zeros along that row
    or column and the algorithm proceeds to the next entry.

  The argument ntot is used to detect whether the function has been run before
  and has to be zero initialised.
*/
void rcont2(tab_t *table, int nr, const tab_t *r, int nc, const tab_t *c,
            const tab_t ntot, bitgen_t *rstate)
{
  // nothing to do
  if (ntot == 0)
    return;

  // jwork is folded into table using last row
  tab_t *jwork = ptr(table, nr, nc, nr - 1, 0);
  // last entry of jwork is never used
  for (int i = 0; i < nc - 1; ++i)
  {
    jwork[i] = c[i];
  }

  tab_t jc = ntot;
  tab_t ib = 0;
  // last row is not random due to constraint
  for (int l = 0; l < nr - 1; ++l)
  {
    tab_t ia = r[l]; // first term
    if (ia == 0)
    {
      for (int i = 0; i < nc; ++i)
        *ptr(table, nr, nc, l, i) = 0;
      continue;
    }
    tab_t ic = jc; // second term
    jc -= r[l];
    // last column is not random due to constraint
    for (int m = 0; m < nc - 1; ++m)
    {
      const tab_t id = jwork[m]; // third term
      const tab_t ie = ic;       // eight term
      ic -= id;
      ib = ie - ia;
      // must be after ib calculation, which is used at the end
      if (c[m] == 0)
      {
        for (int i = 0; i < nr; ++i)
          *ptr(table, nr, nc, i, m) = 0;
        continue;
      }
      const tab_t ii = ib - id; // forth term
      if (ie == 0)
      {
        for (int j = m; j < nc - 1; ++j)
          *ptr(table, nr, nc, l, j) = 0;
        ia = 0;
        break;
      }
      double z = random_standard_uniform(rstate);
      tab_t nlm;
    l131:
      nlm = (tab_t)floor((double)(ia * id) / ie + 0.5);
      double x = exp(
          logfactorial(ia) + logfactorial(ib) + logfactorial(ic) + logfactorial(id) - logfactorial(ie) - logfactorial(nlm) - logfactorial(id - nlm) - logfactorial(ia - nlm) - logfactorial(ii + nlm));
      if (x >= z)
        goto l160;
      double sumprb = x;
      double y = x;
      tab_t nll = nlm;
      bool lsp = false;
      bool lsm = false;
      // increment entry at (l,m)
      tab_t j;
    l140:
      j = (id - nlm) * (ia - nlm);
      if (j == 0)
        goto l156;
      nlm += 1;
      x *= (double)j / (nlm * (ii + nlm));
      sumprb += x;
      if (sumprb >= z)
        goto l160;
    l150:
      if (lsm)
        goto l155;
      // decrement entry at (l,m)
      j = nll * (ii + nll);
      if (j == 0)
        goto l154;
      nll -= 1;
      y *= (double)j / ((id - nll) * (ia - nll));
      sumprb += y;
      if (sumprb >= z)
        goto l159;
      if (!lsp)
        goto l140;
      goto l150;
    l154:
      lsm = true;
    l155:
      if (!lsp)
        goto l140;
      z = random_standard_uniform(rstate) * sumprb;
      goto l131;
    l156:
      lsp = true;
      goto l150;
    l159:
      nlm = nll;
    l160:
      *ptr(table, nr, nc, l, m) = nlm;
      ia -= nlm;
      jwork[m] -= nlm;
    }
    // compute entry in last column of table
    *ptr(table, nr, nc, l, nc - 1) = ia;
  }
  // compute entries in last row of table
  // jwork is already last row of table, so nothing to be done up to nc - 2
  *ptr(table, nr, nc, nr - 1, nc - 1) = ib - *ptr(table, nr, nc, nr - 1, nc - 2);
}
