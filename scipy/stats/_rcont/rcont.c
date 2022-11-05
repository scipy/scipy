/**
  This file contains two algorithms written in C to generate random two-way
  tables. The algorithm rcont1 and rcont2 originate from Boyett and Patefield,
  respectively. For more information, see the docs for each function.

  If you wonder about the spelling, rcont is short for random contingency
  table. Random contingency table is a bit of a misnomer. The tables generated
  by these algorithms have no contingency/association between the two
  variables.
*/

#include <math.h>
#include <numpy/random/distributions.h>

// helper function to access a 1D array like a C-style 2D array
double *ptr(double *m, int nr, int nc, int ir, int ic)
{
  return m + nc * ir + ic;
}

// sanity checks for inputs to rcont1 and rcont2,
// also sets total number of entries
int rcont_check(double *n, const double *m, int nr, const double *r, int nc, const double *c)
{
  if (m == 0 || r == 0 || c == 0 || n == 0)
    return 1;

  if (nr < 2 || nc < 2)
    return 2;

  // check sum(r) == sum(c); r[i] >= 0, c[i] >= 0; sum(r) > 0
  *n = 0;
  for (int i = 0; i < nc; ++i)
  {
    if (!(c[i] >= 0))
      return 3;
    *n += c[i];
  }
  double n2 = 0;
  for (int i = 0; i < nr; ++i)
  {
    if (!(r[i] >= 0))
      return 3;
    n2 += r[i];
  }
  if (*n != n2)
    return 4;
  if (!(*n > 0))
    return 5;

  return 0;
}

/*
  Generate random two-way table with given marginal totals.

  Boyett's shuffling algorithm adapted from AS 144 Appl. Statist. (1979)
  329-332. The algorithm has O(N) complexity in space and time for
  a table with N entries in total. The algorithm performs poorly for large N,
  but is insensitive to the number K of table cells.

  This function uses a work space that is allocated into the argument work
  (which must be zero initialised) and has to freed by the user.
*/
int rcont1(double *matrix, int nr, const double *r, int nc, const double *c,
           int **work, bitgen_t *rstate)
{
  int status = 0;
  if (*work == 0)
  {
    double nd = 0;
    status = rcont_check(&nd, matrix, nr, r, nc, c);
    if (status != 0)
      return status;

    int n = (int)nd;
    *work = (int *)malloc(sizeof(int) * (n + 1));
    if (!*work)
    {
      status = 1;
      return status;
    }
    *work[0] = n;
    int *ymap = *work + 1;
    for (int i = 0; i < nc; ++i)
    {
      int ci = (int)c[i];
      while (ci--)
        *ymap++ = i;
    }
  }

  int n = *work[0];
  int *ymap = *work + 1;

  // shuffle ymap
  for (int i = n - 1; i > 0; --i)
  {
    int j = random_interval(rstate, i);
    int tmp = ymap[j];
    ymap[j] = ymap[i];
    ymap[i] = tmp;
  }

  // clear table
  for (int i = 0, nrc = (nr * nc); i < nrc; ++i)
    matrix[i] = 0;

  // fill table
  for (int ir = 0; ir < nr; ++ir)
  {
    int ri = (int)r[ir];
    while (ri--)
      *ptr(matrix, nr, nc, ir, *ymap++) += 1;
  }

  return 0;
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
int rcont2(double *matrix, int nr, const double *r, int nc, const double *c,
           double *ntot, bitgen_t *rstate)
{
  int status = 0;
  if (*ntot == 0) // perform checks only once
    status = rcont_check(ntot, matrix, nr, r, nc, c);
  if (status != 0)
    return status;

  // jwork is folded into matrix using last row
  double *jwork = ptr(matrix, nr, nc, nr - 1, 0);
  for (int i = 0; i < nc; ++i)
  {
    jwork[i] = c[i];
  }

  double jc = *ntot;
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

  return 0;
}
