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

#include "_rcont.h"
#include <math.h>

// helper function to access a 1D array like a C-style 2D array
double *ptr(double *m, int nr, int nc, int ir, int ic)
{
  return m + nc * ir + ic;
}

/*
  Call this once to initialize workspace for rcont1.

  The work space must have size N, where N is the total number of entries.
*/
void rcont1_init(int *work, int nc, const double *c)
{
  for (int i = 0; i < nc; ++i)
  {
    int ci = (int)c[i];
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
void rcont1(double *matrix, int nr, const double *r, int nc, const double *c,
            double ntot, int *work, bitgen_t *rstate)
{
  // nothing to do
  if (ntot == 0)
    return;

  // shuffle work with Knuth's algorithm
  for (int i = (int)ntot - 1; i > 0; --i)
  {
    int j = random_interval(rstate, i);
    int tmp = work[j];
    work[j] = work[i];
    work[i] = tmp;
  }

  // clear table
  for (int i = 0, nrc = (nr * nc); i < nrc; ++i)
    matrix[i] = 0;

  // fill table
  for (int ir = 0; ir < nr; ++ir)
  {
    int ri = (int)r[ir];
    while (ri--)
      *ptr(matrix, nr, nc, ir, *work++) += 1;
  }
}

/*
  Generate random two-way table with given marginal totals.

  Patefield's algorithm adapted from AS 159 Appl. Statist. (1981) 91-97. This
  algorithm has O(K log(N)) complexity in time for a table with K cells and N
  entries in total. It requires only a small constant stack space.

  The original FORTRAN code was hand-translated to C. Changes to the original:

  - The computation of a look-up table of log-factorials was replaced with
    calls to lgamma, which lifts the limitation that the code only works for
    tables with less then 5000 entries.
  - The data type of input and output arrays was changed to double to minimize
    type conversions. Users are responsible for passing only integral numbers.
  - The original implementation allocated a column vector JWORK, but this is
    not necessary. The vector can be folded into the last column of the output
    matrix.
  - The function uses Numpy's random number generator and distribution library.
  - The algorithm now handles zero entries in row or column vector. When a
    zero is encountered, the output matrix is filled with zeros along that row
    or column and the algorithm proceeds to the next entry.

  The argument ntot is used to detect whether the function has been run before
  and has to be zero initialised.
*/
void rcont2(double *matrix, int nr, const double *r, int nc, const double *c,
            double ntot, bitgen_t *rstate)
{
  // nothing to do
  if (ntot == 0.0)
    return;

  // jwork is folded into matrix using last row
  double *jwork = ptr(matrix, nr, nc, nr - 1, 0);
  // last entry of jwork is never used
  for (int i = 0; i < nc - 1; ++i)
  {
    jwork[i] = c[i];
  }

  double jc = ntot;
  double ib = 0;
  // last row is not random due to constraint
  for (int l = 0; l < nr - 1; ++l)
  {
    double ia = r[l]; // first term
    if (ia == 0)
    {
      for (int i = 0; i < nc; ++i)
        *ptr(matrix, nr, nc, l, i) = 0;
      continue;
    }
    double ic = jc; // second term
    jc -= r[l];
    // last column is not random due to constraint
    for (int m = 0; m < nc - 1; ++m)
    {
      const double id = jwork[m]; // third term
      const double ie = ic;       // eight term
      ic -= id;
      ib = ie - ia;
      // must be after ib calculation, which is used at the end
      if (c[m] == 0)
      {
        for (int i = 0; i < nr; ++i)
          *ptr(matrix, nr, nc, i, m) = 0;
        continue;
      }
      const double ii = ib - id; // forth term
      if (ie == 0)
      {
        for (int j = m; j < nc - 1; ++j)
          *ptr(matrix, nr, nc, l, j) = 0;
        ia = 0;
        break;
      }
      double z = random_standard_uniform(rstate);
      double nlm;
    l131:
      nlm = floor(ia * id / ie + 0.5);
      double x = exp(
          lgamma(ia + 1) + lgamma(ib + 1) + lgamma(ic + 1) + lgamma(id + 1) - lgamma(ie + 1) - lgamma(nlm + 1) - lgamma(id - nlm + 1) - lgamma(ia - nlm + 1) - lgamma(ii + nlm + 1));
      if (x >= z)
        goto l160;
      double sumprb = x;
      double y = x;
      double nll = nlm;
      int lsp = 0;
      int lsm = 0;
      // increment entry at (l,m)
      double j;
    l140:
      j = (id - nlm) * (ia - nlm);
      if (j == 0)
        goto l156;
      nlm += 1;
      x *= j / (nlm * (ii + nlm));
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
      y *= j / ((id - nll) * (ia - nll));
      sumprb += y;
      if (sumprb >= z)
        goto l159;
      if (!lsp)
        goto l140;
      goto l150;
    l154:
      lsm = 1;
    l155:
      if (!lsp)
        goto l140;
      z = random_standard_uniform(rstate) * sumprb;
      goto l131;
    l156:
      lsp = 1;
      goto l150;
    l159:
      nlm = nll;
    l160:
      *ptr(matrix, nr, nc, l, m) = nlm;
      ia -= nlm;
      jwork[m] -= nlm;
    }
    // compute entry in last column of matrix
    *ptr(matrix, nr, nc, l, nc - 1) = ia;
  }
  // compute entries in last row of matrix
  // jwork is already last row of matrix, so nothing to be done up to nc - 2
  *ptr(matrix, nr, nc, nr - 1, nc - 1) = ib - *ptr(matrix, nr, nc, nr - 1, nc - 2);
}
