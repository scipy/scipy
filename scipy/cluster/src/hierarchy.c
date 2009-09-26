/**
 * hierarchy.c
 *
 * Author: Damian Eads
 * Date:   September 22, 2007
 *
 * Copyright (c) 2007, 2008, Damian Eads. All rights reserved.
 * Adapted for incorporation into Scipy, April 9, 2008.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   - Redistributions of source code must retain the above
 *     copyright notice, this list of conditions and the
 *     following disclaimer.
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer
 *     in the documentation and/or other materials provided with the
 *     distribution.
 *   - Neither the name of the author nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <Python.h>
#include <numpy/ndarrayobject.h>

#include "common.h"

#define ISCLUSTER(_nd) ((_nd)->id >= n)
#define GETCLUSTER(_id) ((lists + _id - n))

/** The number of link stats (for the inconsistency computation) for each
    cluster. */

#define CPY_NIS 4

/** The column offsets for the different link stats for the inconsistency
    computation. */
#define CPY_INS_MEAN 0
#define CPY_INS_STD 1
#define CPY_INS_N 2
#define CPY_INS_INS 3

/** The number of linkage stats for each cluster. */
#define CPY_LIS 4

/** The column offsets for the different link stats for the linkage matrix. */
#define CPY_LIN_LEFT 0
#define CPY_LIN_RIGHT 1
#define CPY_LIN_DIST 2
#define CPY_LIN_CNT 3

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "hierarchy.h"

static NPY_INLINE double euclidean_distance(const double *u, const double *v, int n) {
  int i = 0;
  double s = 0.0, d;
  for (i = 0; i < n; i++) {
    d = u[i] - v[i];
    s = s + d * d;
  }
  return sqrt(s);
}

void chopmins(int *ind, int mini, int minj, int np) {
  int i;
  for (i = mini; i < minj - 1; i++) {
    ind[i] = ind[i + 1];
  }
  for (i = minj - 1; i < np - 2; i++) {
    ind[i] = ind[i + 2];
  }
  /**  CPY_DEBUG_MSG("[Remove mini=%d minj=%d]\n", mini, minj);**/
}

void chopmin(int *ind, int minj, int np) {
  int i;
  for (i = minj; i < np - 1; i++) {
    ind[i] = ind[i + 1];
  }
  /**  CPY_DEBUG_MSG("[Remove mini=%d minj=%d]\n", mini, minj);**/
}

void chopmins_ns_ij(double *ind, int mini, int minj, int np) {
  int i;
  for (i = mini; i < minj - 1; i++) {
    ind[i] = ind[i + 1];
  }
  for (i = minj - 1; i < np - 2; i++) {
    ind[i] = ind[i + 2];
  }
}

void chopmins_ns_i(double *ind, int mini, int np) {
  int i;
  for (i = mini; i < np - 1; i++) {
    ind[i] = ind[i + 1];
  }
}

void dist_single(cinfo *info, int mini, int minj, int np, int n) {
  double **rows = info->rows;
  double *buf = info->buf;
  double *bit;
  int i;
  bit = buf;
  for (i = 0; i < mini; i++, bit++) {
    *bit = CPY_MIN(*(rows[i] + mini - i - 1), *(rows[i] + minj - i - 1));
  }
  for (i = mini + 1; i < minj; i++, bit++) {
    *bit = CPY_MIN(*(rows[mini] + i - mini - 1), *(rows[i] + minj - i - 1));
  }
  for (i = minj + 1; i < np; i++, bit++) {
    *bit = CPY_MIN(*(rows[mini] + i - mini - 1), *(rows[minj] + i - minj - 1));
  }
}

void dist_complete(cinfo *info, int mini, int minj, int np, int n) {
  double **rows = info->rows;
  double *buf = info->buf;
  double *bit;
  int i;
  bit = buf;
  for (i = 0; i < mini; i++, bit++) {
    *bit = CPY_MAX(*(rows[i] + mini - i - 1), *(rows[i] + minj - i - 1));
  }
  for (i = mini + 1; i < minj; i++, bit++) {
    *bit = CPY_MAX(*(rows[mini] + i - mini - 1), *(rows[i] + minj - i - 1));
  }
  for (i = minj + 1; i < np; i++, bit++) {
    *bit = CPY_MAX(*(rows[mini] + i - mini - 1), *(rows[minj] + i - minj - 1));
  }
}

void dist_average(cinfo *info, int mini, int minj, int np, int n) {
  double **rows = info->rows, *buf = info->buf, *bit;
  int *inds = info->ind;
  double drx, dsx, mply, rscnt, rc, sc;
  int i, xi, xn;
  cnode *rn = info->nodes + inds[mini];
  cnode *sn = info->nodes + inds[minj];
  cnode *xnd;
  bit = buf;
  rc = (double)rn->n;
  sc = (double)sn->n;
  rscnt = rc + sc;

  for (i = 0; i < mini; i++, bit++) {
    /** d(r,x) **/
    drx = *(rows[i] + mini - i - 1);
    dsx = *(rows[i] + minj - i - 1);
    xi = inds[i];
    xnd = info->nodes + xi;
    xn = xnd->n;
    mply = (double)1.0 / (((double)xn) * rscnt);
    *bit = mply * ((drx * (rc * xn)) + (dsx * (sc * xn)));
  }
  for (i = mini + 1; i < minj; i++, bit++) {
    drx = *(rows[mini] + i - mini - 1);
    dsx = *(rows[i] + minj - i - 1);
    xi = inds[i];
    xnd = info->nodes + xi;
    xn = xnd->n;
    mply = (double)1.0 / (((double)xn) * rscnt);
    *bit = mply * ((drx * (rc * xn)) + (dsx * (sc * xn)));
  }
  for (i = minj + 1; i < np; i++, bit++) {
    drx = *(rows[mini] + i - mini - 1);
    dsx = *(rows[minj] + i - minj - 1);
    xi = inds[i];
    xnd = info->nodes + xi;
    xn = xnd->n;
    mply = (double)1.0 / (((double)xn) * rscnt);
    *bit = mply * ((drx * (rc * xn)) + (dsx * (sc * xn)));
  }
}

void dist_centroid(cinfo *info, int mini, int minj, int np, int n) {
  double *buf = info->buf, *bit;
  int *inds = info->ind;
  const double *centroid_tq;
  int i, m, xi;
  centroid_tq = info->centroids[info->nid];
  bit = buf;
  m = info->m;
  for (i = 0; i < np; i++, bit++) {
    /** d(r,x) **/
    if (i == mini || i == minj) {
      bit--;
      continue;
    }
    xi = inds[i];
    *bit = euclidean_distance(info->centroids[xi], centroid_tq, m);
    /**    CPY_DEBUG_MSG("%5.5f ", *bit);**/
  }
  /**  CPY_DEBUG_MSG("\n");**/
}

void combine_centroids(double *centroidResult,
		       const double *centroidA, const double *centroidB,
		       double na, double nb, int n) {
  int i;
  double nr = (double)na + (double)nb;
  for (i = 0; i < n; i++) {
    centroidResult[i] = ((centroidA[i] * na) + (centroidB[i] * nb)) / nr;
  }
}

void dist_weighted(cinfo *info, int mini, int minj, int np, int n) {
  double **rows = info->rows, *buf = info->buf, *bit;
  int i;
  double drx, dsx;

  bit = buf;

  for (i = 0; i < mini; i++, bit++) {
    /** d(r,x) **/
    drx = *(rows[i] + mini - i - 1);
    dsx = *(rows[i] + minj - i - 1);
    *bit = (drx + dsx) / 2;
  }
  for (i = mini + 1; i < minj; i++, bit++) {
    drx = *(rows[mini] + i - mini - 1);
    dsx = *(rows[i] + minj - i - 1);
    *bit = (drx + dsx) / 2;
  }
  for (i = minj + 1; i < np; i++, bit++) {
    drx = *(rows[mini] + i - mini - 1);
    dsx = *(rows[minj] + i - minj - 1);
    *bit = (drx + dsx) / 2;
  }
  /**  CPY_DEBUG_MSG("\n");**/
}

void dist_ward(cinfo *info, int mini, int minj, int np, int n) {
  double **rows = info->rows, *buf = info->buf, *bit;
  int *inds = info->ind;
  const double *centroid_tq;
  int i, m, xi, rind, sind;
  double drx, dsx, rf, sf, xf, xn, rn, sn, drsSq;
  cnode *newNode;
  cnode *xnd;

  rind = inds[mini];
  sind = inds[minj];
  rn = (double)info->nodes[rind].n;
  sn = (double)info->nodes[sind].n;
  newNode = info->nodes + info->nid;
  drsSq = newNode->d;
  drsSq = drsSq * drsSq;
  centroid_tq = info->centroids[info->nid];
  bit = buf;
  m = info->m;

  for (i = 0; i < mini; i++, bit++) {
    /** d(r,x) **/
    drx = *(rows[i] + mini - i - 1);
    dsx = *(rows[i] + minj - i - 1);
    xi = inds[i];
    xnd = info->nodes + xi;
    xn = xnd->n;
    rf = (rn + xn) / (rn + sn + xn);
    sf = (sn + xn) / (rn + sn + xn);
    xf = -xn / (rn + sn + xn);
    *bit = sqrt(rf * (drx * drx) +
		sf * (dsx * dsx) +
		xf * drsSq);
		
  }
  for (i = mini + 1; i < minj; i++, bit++) {
    drx = *(rows[mini] + i - mini - 1);
    dsx = *(rows[i] + minj - i - 1);
    xi = inds[i];
    xnd = info->nodes + xi;
    xn = xnd->n;
    rf = (rn + xn) / (rn + sn + xn);
    sf = (sn + xn) / (rn + sn + xn);
    xf = -xn / (rn + sn + xn);
    *bit = sqrt(rf * (drx * drx) +
		sf * (dsx * dsx) +
		xf * drsSq);
  }
  for (i = minj + 1; i < np; i++, bit++) {
    drx = *(rows[mini] + i - mini - 1);
    dsx = *(rows[minj] + i - minj - 1);
    xi = inds[i];
    xnd = info->nodes + xi;
    xn = xnd->n;
    rf = (rn + xn) / (rn + sn + xn);
    sf = (sn + xn) / (rn + sn + xn);
    xf = -xn / (rn + sn + xn);
    *bit = sqrt(rf * (drx * drx) +
		sf * (dsx * dsx) +
		xf * drsSq);
  }
  /**  CPY_DEBUG_MSG("\n");**/
}


void print_dm(const double **rows, int np) {
  int i, j, k;
  const double *row;
  CPY_DEBUG_MSG("[DM, np=%d\n", np);
  for (i = 0; i < np - 1; i++) {
    row = rows[i];
    for (j = 0; j <= i; j++) {
      CPY_DEBUG_MSG("%5.5f ", 0.0);
    }

    for (k = 0, j = i + 1; j < np; j++, k++) {
      CPY_DEBUG_MSG("%5.5f ", *(row + k));
    }
    CPY_DEBUG_MSG("|j=%d|\n", i + 1);
  }
}

void print_ind(const int *inds, int np) {
  int i;
  CPY_DEBUG_MSG("[IND, np=%d || ", np);
  for (i = 0; i < np; i++) {
    CPY_DEBUG_MSG("%d ", inds[i]);
  }
  CPY_DEBUG_MSG("]\n");
}

void print_vec(const double *d, int n) {
  int i;
  CPY_DEBUG_MSG("[");
  for (i = 0; i < n; i++) {
    CPY_DEBUG_MSG("%5.5f ", d[i]);
  }
  CPY_DEBUG_MSG("]");
}

/**
 * notes to self:
 * dm:    The distance matrix.
 * Z:     The result of the linkage, a (n-1) x 3 matrix.
 * X:     The original observations as row vectors (=NULL if not needed).
 * n:     The number of objects.
 * ml:    A boolean indicating whether a list of objects in the forest
 *        clusters should be maintained.
 * kc:    Keep track of the centroids.
 */
void linkage(double *dm, double *Z, double *X,
	     int m, int n, int ml, int kc, distfunc dfunc,
	     int method) {
  int i, j, k, t, np, nid, mini, minj, npc2;
  double min, ln, rn, qn;
  int *ind;
  /** An iterator through the distance matrix. */
  double *dmit, *buf;

  int *rowsize;

  /** Temporary array to store modified distance matrix. */
  double *dmt, **rows, *Zrow;
  double *centroidsData;
  double **centroids;
  const double *centroidL, *centroidR;
  double *centroid;
  clist *lists, *listL, *listR, *listC;
  clnode *lnodes;
  cnode *nodes, *node;

  cinfo info;

  /** The next two are only necessary for euclidean distance methods. */
  if (ml) {
    lists = (clist*)malloc(sizeof(clist) * (n-1));
    lnodes = (clnode*)malloc(sizeof(clnode) * n);
  }
  else {
    lists = 0;
    lnodes = 0;
  }
  if (kc) {
    centroids = (double**)malloc(sizeof(double*) * (2 * n));
    centroidsData = (double*)malloc(sizeof(double) * n * m);
    for (i = 0; i < n; i++) {
      centroids[i] = X + i * m;
    }
    for (i = 0; i < n; i++) {
      centroids[i+n] = centroidsData + i * m;
    }
  }
  else {
    centroids = 0;
    centroidsData = 0;
  }

  nodes = (cnode*)malloc(sizeof(cnode) * (n * 2) - 1);
  ind = (int*)malloc(sizeof(int) * n);
  dmt = (double*)malloc(sizeof(double) * NCHOOSE2(n));
  buf = (double*)malloc(sizeof(double) * n);
  rows = (double**)malloc(sizeof(double*) * n);
  rowsize = (int*)malloc(sizeof(int) * n);
  memcpy(dmt, dm, sizeof(double) * NCHOOSE2(n));

  info.X = X;
  info.m = m;
  info.n = n;
  info.nodes = nodes;
  info.ind = ind;
  info.dmt = dmt;
  info.buf = buf;
  info.rows = rows;
  info.rowsize = rowsize;
  info.dm = dm;
  info.centroids = centroids;
  if (kc) {
    info.centroidBuffer = centroids[2*n - 1];
  }
  else {
    info.centroidBuffer = 0;
  }
  info.lists = lists;
  for (i = 0; i < n; i++) {
    ind[i] = i;
    node = nodes + i;
    node->left = 0;
    node->right = 0;
    node->id = i;
    node->n = 1;
    node->d = 0.0;
    rowsize[i] = n - 1 - i;
  }
  rows[0] = dmt;
  for (i = 1; i < n; i++) {
    rows[i] = rows[i-1] + n - i;
  }
  
  if (ml) {
    for (i = 0; i < n; i++) {
      (lnodes + i)->val = nodes + i;
      (lnodes + i)->next = 0;
    }
  }

  for (k = 0, nid = n; k < n - 1; k++, nid++) {
    info.nid = nid;
    np = n - k;
    npc2 = NCHOOSE2(np);
    /**    CPY_DEBUG_MSG("k=%d, nid=%d, n=%d np=%d\n", k, nid, n, np);**/
    min = dmt[0];
    mini = 0;
    minj = 1;
    /** Note that mini < minj since j > i is always true. */
    for (i = 0; i < np - 1; i++) {
      dmit = rows[i];
      for (j = i + 1; j < np; j++, dmit++) {
	if (*dmit <= min) {
	  min = *dmit;
	  mini = i;
	  minj = j;
	}
      }
    }

    node = nodes + nid;
    node->left = nodes + ind[mini];
    node->right = nodes + ind[minj];
    ln = (double)node->left->n;
    rn = (double)node->right->n;
    qn = ln + rn;
    node->n = node->left->n + node->right->n;
    node->d = min;
    node->id = nid;

    Zrow = Z + (k * CPY_LIS);
    Zrow[CPY_LIN_LEFT] = node->left->id;
    Zrow[CPY_LIN_RIGHT] = node->right->id;
    Zrow[CPY_LIN_DIST] = min;
    Zrow[CPY_LIN_CNT] = node->n;

    /**    fCPY_DEBUG_MSG(stderr,
	    "[lid=%d, rid=%d, llid=%d, rrid=%d m=%5.8f]",
	    node->left->id, node->right->id, ind[mini], ind[minj], min);**/

    if (ml) {
      listC = GETCLUSTER(nid);
      if (ISCLUSTER(node->left) != 0) {
	listL = GETCLUSTER(node->left->id);
	if (ISCLUSTER(node->right) != 0) {
	  listR = GETCLUSTER(node->right->id);
	  listL->tail->next = listR->head;
	  listC->tail = listR->tail;
	  listR->tail->next = 0;
	}
	else {
	  listC->tail = lnodes + node->right->id;
	  listL->tail->next = listC->tail;
	  listC->tail->next = 0;
	}
	listC->head = listL->head;
      }
      else {
	listC->head = lnodes + node->left->id;
	if (ISCLUSTER(node->right)) {
	  listR = GETCLUSTER(node->right->id);
	  listC->head->next = listR->head;
	  listC->tail = listR->tail;
	  listC->tail->next = 0;
	}
	else {
	  listC->tail = lnodes + node->right->id;
	  listC->tail->next = 0;
	  listC->head->next = listC->tail;
	}
      }
    }
    if (kc) {
      centroidL = centroids[ind[mini]];
      centroidR = centroids[ind[minj]];
      centroid = centroids[nid];
      switch(method) {
      case CPY_LINKAGE_MEDIAN:
	for (t = 0; t < m; t++) {
	  centroid[t] = (centroidL[t] * 0.5 + centroidR[t] * 0.5);
	}
	break;
      case CPY_LINKAGE_CENTROID:
      case CPY_LINKAGE_WARD:
      default:
	for (t = 0; t < m; t++) {
	  centroid[t] = (centroidL[t] * ln + centroidR[t] * rn) / qn;
	}
	break;
      }
      /**      CPY_DEBUG_MSG("L: ");
      print_vec(centroidL, m);
      CPY_DEBUG_MSG("\nR: ");
      print_vec(centroidR, m);
      CPY_DEBUG_MSG("\nT: ");
      print_vec(centroid, m);**/
    }

    /**    print_dm(rows, np);**/
    /**    dfunc(buf, rows, mini, minj, np, dm, n, ind, nodes);**/
    dfunc(&info, mini, minj, np, n);

    /** For these rows, we must remove, i and j but leave all unused space
        at the end. This reduces their size by two.*/
    for (i = 0; i < mini; i++) {
      chopmins_ns_ij(rows[i], mini - i - 1, minj - i - 1, rowsize[i]);
    }

    /** We skip the i'th row. For rows i+1 up to j-1, we just remove j. */
    for (i = mini + 1; i < minj; i++) {
      chopmins_ns_i(rows[i], minj - i - 1, rowsize[i]);
    }

    /** For rows 0 to mini - 1, we move them down the matrix, leaving the
	first row free. */
    /**    for (i = mini; i > 0; i--) {
      memcpy(rows[i], rows[i-1], sizeof(double) * rowsize[i]-k);
      }**/

    for (i = mini; i < minj - 1; i++) {
      memcpy(rows[i], rows[i+1], sizeof(double) * (rowsize[i+1]));
    }

    /** For rows mini+1 to minj-1, we do nothing since they are in the
	right place for the next iteration. For rows minj+1 onward,
	we move them to the right. */
	
    for (i = minj - 1; i < np - 2; i++) {
      memcpy(rows[i], rows[i+2], sizeof(double) * (rowsize[i+2]));
    }

    /** Rows i+1 to j-1 lose one unit of space, so we move them up. */
    /** Rows j to np-1 lose no space. We do nothing to them. */

    /**    memcpy(rows[0], buf, sizeof(double) * rowsize[0] - k);*/

    for (i = 0; i < np - 2; i++) {
      *(rows[i] + np - 3 - i) = buf[i];
    }

    /**    print_dm(rows, np - 1);
	   print_ind(ind, np);**/
    chopmins(ind, mini, minj, np);
    ind[np - 2] = nid;
    /**    print_ind(ind, np - 1);**/
  }
  free(lists);
  free(lnodes);
  free(nodes);
  free(ind);
  free(dmt);
  free(buf);
  free(rows);
  free(rowsize);
  free(centroidsData);
  free(centroids);
}

/** Trying to reimplement so that output is consistent with MATLAB's in
    cases where there are is than one correct choice to make at each
    iteration of the algorithm. This implementation is not active. */

void linkage_alt(double *dm, double *Z, double *X,
	     int m, int n, int ml, int kc, distfunc dfunc,
	     int method) {
  int i, j, k, t, np, nid, mini, minj, npc2;
  double min, ln, rn, qn;
  int *ind;
  /** An iterator through the distance matrix. */
  double *dmit, *buf;

  int *rowsize;

  /** Temporary array to store modified distance matrix. */
  double *dmt, **rows, *Zrow;
  double *centroidsData;
  double **centroids;
  const double *centroidL, *centroidR;
  double *centroid;
  clist *lists, *listL, *listR, *listC;
  clnode *lnodes;
  cnode *nodes, *node;

  cinfo info;

  /** The next two are only necessary for euclidean distance methods. */
  if (ml) {
    lists = (clist*)malloc(sizeof(clist) * (n-1));
    lnodes = (clnode*)malloc(sizeof(clnode) * n);
  }
  else {
    lists = 0;
    lnodes = 0;
  }
  if (kc) {
    centroids = (double**)malloc(sizeof(double*) * (2 * n));
    centroidsData = (double*)malloc(sizeof(double) * n * m);
    for (i = 0; i < n; i++) {
      centroids[i] = X + i * m;
    }
    for (i = 0; i < n; i++) {
      centroids[i+n] = centroidsData + i * m;
    }
  }
  else {
    centroids = 0;
    centroidsData = 0;
  }

  nodes = (cnode*)malloc(sizeof(cnode) * (n * 2) - 1);
  ind = (int*)malloc(sizeof(int) * n);
  dmt = (double*)malloc(sizeof(double) * NCHOOSE2(n));
  buf = (double*)malloc(sizeof(double) * n);
  rows = (double**)malloc(sizeof(double*) * n);
  rowsize = (int*)malloc(sizeof(int) * n);
  memcpy(dmt, dm, sizeof(double) * NCHOOSE2(n));

  info.X = X;
  info.m = m;
  info.n = n;
  info.nodes = nodes;
  info.ind = ind;
  info.dmt = dmt;
  info.buf = buf;
  info.rows = rows;
  info.rowsize = rowsize;
  info.dm = dm;
  info.centroids = centroids;
  if (kc) {
    info.centroidBuffer = centroids[2*n - 1];
  }
  else {
    info.centroidBuffer = 0;
  }
  info.lists = lists;
  for (i = 0; i < n; i++) {
    ind[i] = i;
    node = nodes + i;
    node->left = 0;
    node->right = 0;
    node->id = i;
    node->n = 1;
    node->d = 0.0;
    rowsize[i] = n - 1 - i;
  }
  rows[0] = dmt;
  for (i = 1; i < n; i++) {
    rows[i] = rows[i-1] + n - i;
  }
  
  if (ml) {
    for (i = 0; i < n; i++) {
      (lnodes + i)->val = nodes + i;
      (lnodes + i)->next = 0;
    }
  }

  for (k = 0, nid = n; k < n - 1; k++, nid++) {
    info.nid = nid;
    np = n - k;
    npc2 = NCHOOSE2(np);
    /**    CPY_DEBUG_MSG("k=%d, nid=%d, n=%d np=%d\n", k, nid, n, np);**/
    min = dmt[0];
    mini = 0;
    minj = 1;
    /** Note that mini < minj since j > i is always true. */
    /** BEGIN NEW CODE **/
    for (i = 0; i < np - 1; i++) {
      dmit = rows[i];
      for (j = i + 1; j < np; j++, dmit++) {
	if (*dmit < min) {
	  min = *dmit;
	  mini = i;
	  minj = j;
	}
      }
    }

    node = nodes + nid;
    node->left = nodes + ind[mini];
    node->right = nodes + ind[minj];
    ln = (double)node->left->n;
    rn = (double)node->right->n;
    qn = ln + rn;
    node->n = node->left->n + node->right->n;
    node->d = min;
    node->id = nid;

    Zrow = Z + (k * CPY_LIS);
    Zrow[CPY_LIN_LEFT] = node->left->id;
    Zrow[CPY_LIN_RIGHT] = node->right->id;
    Zrow[CPY_LIN_DIST] = min;
    Zrow[CPY_LIN_CNT] = node->n;

    /**    fprintf(stderr,
	    "[lid=%d, rid=%d, llid=%d, rrid=%d m=%5.8f]",
	    node->left->id, node->right->id, ind[mini], ind[minj], min);**/

    if (ml) {
      listC = GETCLUSTER(nid);
      if (ISCLUSTER(node->left) != 0) {
	listL = GETCLUSTER(node->left->id);
	if (ISCLUSTER(node->right) != 0) {
	  listR = GETCLUSTER(node->right->id);
	  listL->tail->next = listR->head;
	  listC->tail = listR->tail;
	  listR->tail->next = 0;
	}
	else {
	  listC->tail = lnodes + node->right->id;
	  listL->tail->next = listC->tail;
	  listC->tail->next = 0;
	}
	listC->head = listL->head;
      }
      else {
	listC->head = lnodes + node->left->id;
	if (ISCLUSTER(node->right)) {
	  listR = GETCLUSTER(node->right->id);
	  listC->head->next = listR->head;
	  listC->tail = listR->tail;
	  listC->tail->next = 0;
	}
	else {
	  listC->tail = lnodes + node->right->id;
	  listC->tail->next = 0;
	  listC->head->next = listC->tail;
	}
      }
    }
    if (kc) {
      centroidL = centroids[ind[mini]];
      centroidR = centroids[ind[minj]];
      centroid = centroids[nid];
      switch(method) {
      case CPY_LINKAGE_MEDIAN:
	for (t = 0; t < m; t++) {
	  centroid[t] = (centroidL[t] * 0.5 + centroidR[t] * 0.5);
	}
	break;
      case CPY_LINKAGE_CENTROID:
      case CPY_LINKAGE_WARD:
      default:
	for (t = 0; t < m; t++) {
	  centroid[t] = (centroidL[t] * ln + centroidR[t] * rn) / qn;
	}
	break;
      }
      /**      CPY_DEBUG_MSG("L: ");
      print_vec(centroidL, m);
      CPY_DEBUG_MSG("\nR: ");
      print_vec(centroidR, m);
      CPY_DEBUG_MSG("\nT: ");
      print_vec(centroid, m);**/
    }

    /**    print_dm(rows, np);**/
    /**    dfunc(buf, rows, mini, minj, np, dm, n, ind, nodes);**/
    dfunc(&info, mini, minj, np, n);

    /** For these rows, we must remove, i and j but leave all unused space
        at the end. This reduces their size by two.*/
    for (i = 0; i < minj; i++) {
      chopmins_ns_i(rows[i], minj - i - 1, rowsize[i]);
    }

    /** We skip the i'th row. For rows i+1 up to j-1, we just remove j. */
    /**for (i = mini + 1; i < minj; i++) {
      chopmins_ns_i(rows[i], minj - i - 1, rowsize[i]);
      }**/

    /** For rows 0 to mini - 1, we move them down the matrix, leaving the
	first row free. */
    /**for (i = mini; i > 0; i--) {
      memcpy(rows[i], rows[i-1], sizeof(double) * rowsize[i]-k);
    }

    for (i = mini; i < minj - 1; i++) {
      memcpy(rows[i], rows[i+1], sizeof(double) * (rowsize[i+1]));
      }**/

    /** For rows mini+1 to minj-1, we do nothing since they are in the
	right place for the next iteration. For rows minj+1 onward,
	we move them to the right. */
	
    for (i = minj; i < np - 1; i++) {
      memcpy(rows[i], rows[i+1], sizeof(double) * (rowsize[i+1]));
    }

    /** Rows i+1 to j-1 lose one unit of space, so we move them up. */
    /** Rows j to np-1 lose no space. We do nothing to them. */
    /**    memcpy(rows[0], buf, sizeof(double) * rowsize[0] - k);*/

    for (i = 0; i < mini; i++) {
      *(rows[i] + mini - i - 1) = buf[i];
    }

    for (i = mini + 1; i < np - 2; i++) {
      *(rows[mini] + i - mini - 1) = buf[i-1];
    }

    /**    print_dm(rows, np - 1);
	   print_ind(ind, np);**/
    chopmin(ind, minj, np);
    ind[mini] = nid;
    /**    print_ind(ind, np - 1);**/
  }
  free(lists);
  free(lnodes);
  free(nodes);
  free(ind);
  free(dmt);
  free(buf);
  free(rows);
  free(rowsize);
  free(centroidsData);
  free(centroids);
}

void cpy_to_tree(const double *Z, cnode **tnodes, int n) {
  const double *row;
  cnode *node;
  cnode *nodes;
  int i;
  nodes = (cnode*)malloc(sizeof(cnode) * (n * 2) - 1);
  *tnodes = nodes;
  for (i = 0; i < n; i++) {
    node = nodes + i;
    node->left = 0;
    node->right = 0;
    node->id = i;
    node->n = 1;
    node->d = 0.0;
  }
  for (i = 0; i < n - 1; i++) {
    node = nodes + i + n;
    row = Z + (i * CPY_LIS);
    node->id = i + n;
    node->left = nodes + (int)row[CPY_LIN_LEFT];
    node->right = nodes + (int)row[CPY_LIN_RIGHT];
    node->d = row[CPY_LIN_DIST];
    node->n = (int)row[CPY_LIN_CNT];
    /**    CPY_DEBUG_MSG("l: %d r: %d d: %5.5f n: %d\n", (int)row[0],
	   (int)row[1], row[2], (int)row[3]);**/
  }
}

NPY_INLINE void set_dist_entry(double *d, double val, int i, int j, int n) {
  if (i < j) {
    *(d + (NCHOOSE2(n)-NCHOOSE2(n - i)) + j) = val;
  }
  if (j < i) {
    *(d + (NCHOOSE2(n)-NCHOOSE2(n - j)) + i) = val;
  }
}

void cophenetic_distances(const double *Z, double *d, int n) {
  int *curNode, *left;
  int ndid, lid, rid, i, j, k, t = 0, ln, rn, ii, jj, nc2;
  unsigned char *lvisited, *rvisited;
  const double *Zrow;
  int *members = (int*)malloc(n * sizeof(int));
  const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(n);
  k = 0;
  curNode = (int*)malloc(n * sizeof(int));
  left = (int*)malloc(n * sizeof(int));
  lvisited = (unsigned char*)malloc(bff);
  rvisited = (unsigned char*)malloc(bff);
  curNode[k] = (n * 2) - 2;
  left[k] = 0;
  nc2 = NCHOOSE2(n);
  memset(lvisited, 0, bff);
  memset(rvisited, 0, bff);

  while (k >= 0) {
    ndid = curNode[k];
    Zrow = Z + ((ndid-n) * CPY_LIS);
    lid = (int)Zrow[CPY_LIN_LEFT];
    rid = (int)Zrow[CPY_LIN_RIGHT];
    if (lid >= n) {
      ln = (int)*(Z + (CPY_LIS * (lid-n)) + CPY_LIN_CNT);
    }
    else {
      ln = 1;
    }
    if (rid >= n) {
      rn = (int)*(Z + (CPY_LIS * (rid-n)) + CPY_LIN_CNT);
    }
    else {
      rn = 1;
    }
    
    /**    CPY_DEBUG_MSG("[fp] ndid=%d, ndid-n=%d, k=%d, lid=%d, rid=%d\n",
	   ndid, ndid-n, k, lid, rid);**/

    if (lid >= n && !CPY_GET_BIT(lvisited, ndid-n)) {
      CPY_SET_BIT(lvisited, ndid-n);
      curNode[k+1] = lid;
      left[k+1] = left[k];
      k++;
      continue;
    }
    else if (lid < n) {
      members[left[k]] = lid;
    }
    if (rid >= n && !CPY_GET_BIT(rvisited, ndid-n)) {
      CPY_SET_BIT(rvisited, ndid-n);
      curNode[k+1] = rid;
      left[k+1] = left[k] + ln;
      k++;
      continue;
    }
    else if (rid < n) {
      members[left[k]+ln] = rid;
    }

    /** If it's not a leaf node, and we've visited both children,
	record the final mean in the table. */
    if (ndid >= n) {
      for (ii = 0; ii < ln; ii++) {
	i = *(members + left[k] + ii);
	for (jj = 0; jj < rn; jj++) {
	  j = *(members + left[k] + ln + jj);
	  if (i < j) {
	    t = nc2 - NCHOOSE2(n - i) + (j - i - 1);
	  }
	  if (j < i) {
	    t = nc2 - NCHOOSE2(n - j) + (i - j - 1);
	  }
	  d[t] = Zrow[CPY_LIN_DIST];
	  /**	CPY_DEBUG_MSG("i=%d j=%d k=%d d=%5.5f \n", i, j, k, dist);**/
	}
      }
    }
    k--;
  }
  free(members);
  free(left);
  free(curNode);
  free(lvisited);
  free(rvisited);
}

void inconsistency_calculation_alt(const double *Z, double *R, int n, int d) {
  int *curNode;
  int ndid, lid, rid, i, k;
  unsigned char *lvisited, *rvisited;
  const double *Zrow;
  double *Rrow;
  double levelSum, levelStdSum;
  int levelCnt;
  const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(n);
  k = 0;
  curNode = (int*)malloc(n * sizeof(int));
  lvisited = (unsigned char*)malloc(bff);
  rvisited = (unsigned char*)malloc(bff);
  /** for each node in the original linkage matrix. */
  for (i = 0; i < n - 1; i++) {
    /** the current depth j */
    k = 0;
    levelSum = 0.0;
    levelCnt = 0;
    levelStdSum = 0.0;
    memset(lvisited, 0, bff);
    memset(rvisited, 0, bff);
    curNode[0] = i;
    for (k = 0; k >= 0;) {
      ndid = curNode[k];
      Zrow = Z + ((ndid) * CPY_LIS);
      lid = (int)Zrow[CPY_LIN_LEFT];
      rid = (int)Zrow[CPY_LIN_RIGHT];
      /** CPY_DEBUG_MSG("[fp] ndid=%d, ndid-n=%d, k=%d, lid=%d, rid=%d\n",
	          ndid, ndid, k, lid, rid);**/
      if (k < d - 1) {
	if (lid >= n && !CPY_GET_BIT(lvisited, ndid)) {
	  CPY_SET_BIT(lvisited, ndid);
	  k++;
	  curNode[k] = lid-n;
	  continue;
	}
	if (rid >= n && !CPY_GET_BIT(rvisited, ndid)) {
	  CPY_SET_BIT(rvisited, ndid);
	  k++;
	  curNode[k] = rid-n;
	  continue;
	}
      }
      levelCnt++;
      levelSum += Zrow[CPY_LIN_DIST];
      levelStdSum += Zrow[CPY_LIN_DIST] * Zrow[CPY_LIN_DIST];
	/**CPY_DEBUG_MSG("  Using range %d to %d, levelCnt[k]=%d\n", lb, ub, levelCnt[k]);**/
      /** Let the count and sum slots be used for the next newly visited
	  node. */
      k--;
    }
    Rrow = R + (CPY_NIS * i);
    Rrow[CPY_INS_N] = (double)levelCnt;
    Rrow[CPY_INS_MEAN] = levelSum / levelCnt;
    if (levelCnt < 2) {
      Rrow[CPY_INS_STD] = (levelStdSum - (levelSum * levelSum)) / levelCnt;
    }
    else {
      Rrow[CPY_INS_STD] = (levelStdSum - ((levelSum * levelSum) / levelCnt)) / (levelCnt - 1);
    }
    Rrow[CPY_INS_STD] = sqrt(CPY_MAX(0, Rrow[CPY_INS_STD]));
    if (Rrow[CPY_INS_STD] > 0) {
      Rrow[CPY_INS_INS] = (Zrow[CPY_LIN_DIST] - Rrow[CPY_INS_MEAN]) / Rrow[CPY_INS_STD];
    }
  }
  
  free(curNode);
  free(lvisited);
  free(rvisited);
}

void calculate_cluster_sizes(const double *Z, double *cs, int n) {
  int i, j, k, q;
  const double *row;
  for (k = 0; k < n - 1; k++) {
    row = Z + (k * 3);
    i = (int)row[CPY_LIN_LEFT];
    j = (int)row[CPY_LIN_RIGHT];
    /** If the left node is a non-singleton, add its count. */
    if (i >= n) {
      q = i - n;
      cs[k] += cs[q];
    }
    /** Otherwise just add 1 for the leaf. */
    else {
      cs[k] += 1.0;
    }
    /** If the right node is a non-singleton, add its count. */
    if (j >= n) {
      q = j - n;
      cs[k] += cs[q];
    }
    /** Otherwise just add 1 for the leaf. */
    else {
      cs[k] += 1.0;
    }
    CPY_DEBUG_MSG("i=%d, j=%d, cs[%d]=%d\n", i, j, k, (int)cs[k]);
  }
}

/** Returns an array of original observation indices (pre-order traversal). */
void form_member_list(const double *Z, int *members, int n) {
  int *curNode, *left;
  int ndid, lid, rid, k, ln, rn;
  unsigned char *lvisited, *rvisited;
  const double *Zrow;
  const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(n);

  k = 0;
  curNode = (int*)malloc(n * sizeof(int));
  left = (int*)malloc(n * sizeof(int));
  lvisited = (unsigned char*)malloc(bff);
  rvisited = (unsigned char*)malloc(bff);
  curNode[k] = (n * 2) - 2;
  left[k] = 0;
  memset(lvisited, 0, bff);
  memset(rvisited, 0, bff);

  while (k >= 0) {
    ndid = curNode[k];
    Zrow = Z + ((ndid-n) * CPY_LIS);
    lid = (int)Zrow[CPY_LIN_LEFT];
    rid = (int)Zrow[CPY_LIN_RIGHT];
    if (lid >= n) {
      ln = (int)*(Z + (CPY_LIS * (lid-n)) + CPY_LIN_CNT);
    }
    else {
      ln = 1;
    }
    if (rid >= n) {
      rn = (int)*(Z + (CPY_LIS * (rid-n)) + CPY_LIN_CNT);
    }
    else {
      rn = 1;
    }
    
    /**    CPY_DEBUG_MSG("[fp] ndid=%d, ndid-n=%d, k=%d, lid=%d, rid=%d\n",
	   ndid, ndid-n, k, lid, rid);**/

    if (lid >= n && !CPY_GET_BIT(lvisited, ndid-n)) {
      CPY_SET_BIT(lvisited, ndid-n);
      curNode[k+1] = lid;
      left[k+1] = left[k];
      k++;
      continue;
    }
    else if (lid < n) {
      members[left[k]] = lid;
    }
    if (rid >= n && !CPY_GET_BIT(rvisited, ndid-n)) {
      CPY_SET_BIT(rvisited, ndid-n);
      curNode[k+1] = rid;
      left[k+1] = left[k] + ln;
      k++;
      continue;
    }
    else if (rid < n) {
      members[left[k]+ln] = rid;
    }
    k--;
  }
  free(curNode);
  free(left);
  free(lvisited);
  free(rvisited);
}

void form_flat_clusters_from_in(const double *Z, const double *R, int *T,
				double cutoff, int n) {
  double *max_inconsists = (double*)malloc(sizeof(double) * n);
  get_max_Rfield_for_each_cluster(Z, R, max_inconsists, n, 3);
  form_flat_clusters_from_monotonic_criterion(Z, max_inconsists, T, cutoff, n);
  free(max_inconsists);
}

void form_flat_clusters_from_dist(const double *Z, int *T,
				  double cutoff, int n) {
  double *max_dists = (double*)malloc(sizeof(double) * n);
  get_max_dist_for_each_cluster(Z, max_dists, n);
  CPY_DEBUG_MSG("cupid: n=%d cutoff=%5.5f MD[0]=%5.5f MD[n-1]=%5.5f\n", n, cutoff, max_dists[0], max_dists[n-2]);
  form_flat_clusters_from_monotonic_criterion(Z, max_dists, T, cutoff, n);
  free(max_dists);
}

void form_flat_clusters_maxclust_dist(const double *Z, int *T, int n, int mc) {
  
  double *MD = (double*)malloc(sizeof(double) * n);
  get_max_dist_for_each_cluster(Z, MD, n);
  CPY_DEBUG_MSG("fumble: n=%d mc=%d MD[0]=%5.5f MD[n-1]=%5.5f\n", n, mc, MD[0], MD[n-2]);
  form_flat_clusters_maxclust_monocrit(Z, MD, T, n, mc);
  free(MD);
}
						 
/** form flat clusters by thresholding a monotonic criterion. */
void form_flat_clusters_from_monotonic_criterion(const double *Z,
						 const double *mono_crit,
						 int *T, double cutoff, int n) {
  int *curNode;
  int ndid, lid, rid, k, ms, nc;
  unsigned char *lvisited, *rvisited;
  double max_crit;
  const double *Zrow;
  const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(n);

  curNode = (int*)malloc(n * sizeof(int));
  lvisited = (unsigned char*)malloc(bff);
  rvisited = (unsigned char*)malloc(bff);

  /** number of clusters formed so far. */
  nc = 0;
  /** are we in part of a tree below the cutoff? .*/
  ms = -1;
  k = 0;
  curNode[k] = (n * 2) - 2;
  memset(lvisited, 0, bff);
  memset(rvisited, 0, bff);
  ms = -1;
  while (k >= 0) {
    ndid = curNode[k];
    Zrow = Z + ((ndid-n) * CPY_LIS);
    lid = (int)Zrow[CPY_LIN_LEFT];
    rid = (int)Zrow[CPY_LIN_RIGHT];
    max_crit = mono_crit[ndid-n];
    CPY_DEBUG_MSG("cutoff: %5.5f maxc: %5.5f nc: %d\n", cutoff, max_crit, nc);
    if (ms == -1 && max_crit <= cutoff) {
      CPY_DEBUG_MSG("leader: i=%d\n", ndid);
      ms = k;
      nc++;
    }
    if (lid >= n && !CPY_GET_BIT(lvisited, ndid-n)) {
      CPY_SET_BIT(lvisited, ndid-n);
      curNode[k+1] = lid;
      k++;
      continue;
    }
    if (rid >= n && !CPY_GET_BIT(rvisited, ndid-n)) {
      CPY_SET_BIT(rvisited, ndid-n);
      curNode[k+1] = rid;
      k++;
      continue;
    }
    if (ndid >= n) {
      if (lid < n) {
	if (ms == -1) {
	  nc++;
	  T[lid] = nc;
	}
	else {
	  T[lid] = nc;
	}
      }
      if (rid < n) {
	if (ms == -1) {
	  nc++;
	  T[rid] = nc;
	}
	else {
	  T[rid] = nc;
	}
      }
      if (ms == k) {
	ms = -1;
      }
    }
    k--;
  }

  free(curNode);
  free(lvisited);
  free(rvisited);  
}

void form_flat_clusters_maxclust_monocrit(const double *Z,
					  const double *mono_crit,
					  int *T, int n, int mc) {
  int *curNode;
  int ndid, lid, rid, k, nc, g, ms;
  unsigned char *lvisited, *rvisited;
  const double *Zrow;
  double thresh, maxmono_crit;
  /** The maximum unsuccessful distance is initially -1.0 (hack). */
  double max_illegal = -1.0;
  double min_legal = 0.0;
  int min_legal_nc = 1;
  const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(n);
  k = 0;

  min_legal = mono_crit[n-2];
  curNode = (int*)malloc(n * sizeof(int));
  lvisited = (unsigned char*)malloc(bff);
  rvisited = (unsigned char*)malloc(bff);
  curNode[k] = (n * 2) - 2;
  memset(lvisited, 0, bff);
  memset(rvisited, 0, bff);

  /** number of clusters formed so far. */
  nc = 0;

  CPY_DEBUG_MSG("[BEGIN] min legal: %5.5f nc: %d mc: %d\n", min_legal, min_legal_nc, mc);

  for (g = n - 2; g >= 0; g--) {
    thresh = mono_crit[g];
    /** 1. If the threshold is <= the minimum threshold we've tried
        unsuccessfully, skip the threshold. (or)

        2. If the threshold is > the minimum legal threshold, it is
	less optimal so skip it. */
    if (thresh > min_legal) { /** commented out : && thresh <= max_illegal **/
      continue;
    }
    k = 0;
    curNode[k] = (n * 2) - 2;
    memset(lvisited, 0, bff);
    memset(rvisited, 0, bff);
    nc = 0;
    ms = -1;
    /** See if the threshold MD[g] works. **/
    while (k >= 0) {
      ndid = curNode[k];
      Zrow = Z + ((ndid-n) * CPY_LIS);
      lid = (int)Zrow[CPY_LIN_LEFT];
      rid = (int)Zrow[CPY_LIN_RIGHT];
      maxmono_crit = mono_crit[ndid-n];
      /**      CPY_DEBUG_MSG("cutoff: %5.5f maxi: %5.5f nc: %d\n", cutoff, max_mono_crit, nc);**/

      /** If the current nodes maxmono_crit is <= the threshold, stop exploring
	  deeper in the tree. The node and its descendent leaves will be their
	  own cluster. */
      if (maxmono_crit <= thresh) {
	nc++;
	k--;
	CPY_SET_BIT(lvisited, ndid-n);
	CPY_SET_BIT(rvisited, ndid-n);
	continue;
      }
      /** Otherwise, the node is above the threshold, so we need to explore
	  it's children. */
      if (!CPY_GET_BIT(lvisited, ndid-n)) {
	CPY_SET_BIT(lvisited, ndid-n);
	if (lid >= n) {
	  curNode[k+1] = lid;
	  k++;
	  continue;
	}
	else if (lid < n) {
	  nc++;
	}
      }
      if (!CPY_GET_BIT(rvisited, ndid-n)) {
	if (rid >= n) {
	  CPY_SET_BIT(rvisited, ndid-n);
	  curNode[k+1] = rid;
	  k++;
	  continue;
	}
	else if (rid < n) {
	  nc++;
	}
      }
      k--;
    }

    if (thresh > max_illegal && nc > mc) {
      CPY_DEBUG_MSG("max illegal: %5.5f mc: %d", max_illegal, mc);
      max_illegal = thresh;
      continue;
    }
    /** If the threshold is less than the current minimum legal threshold
	but has a legal number of clusters, set the new legal minimum. */
    if (thresh < min_legal && nc <= mc) {
      min_legal = thresh;
      min_legal_nc = nc;
      CPY_DEBUG_MSG("min legal: %5.5f nc: %d mc: %d\n", min_legal, min_legal_nc, mc);
    }
  }

  form_flat_clusters_from_monotonic_criterion(Z, mono_crit, T, min_legal, n);

  free(curNode);
  free(lvisited);
  free(rvisited);
}

void get_max_dist_for_each_cluster(const double *Z, double *max_dists, int n) {
  int *curNode;
  int ndid, lid, rid, k;
  unsigned char *lvisited, *rvisited;
  const double *Zrow;
  double max_dist;
  const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(n);

  k = 0;
  curNode = (int*)malloc(n * sizeof(int));
  lvisited = (unsigned char*)malloc(bff);
  rvisited = (unsigned char*)malloc(bff);
  curNode[k] = (n * 2) - 2;
  memset(lvisited, 0, bff);
  memset(rvisited, 0, bff);
  while (k >= 0) {
    ndid = curNode[k];
    Zrow = Z + ((ndid-n) * CPY_LIS);
    lid = (int)Zrow[CPY_LIN_LEFT];
    rid = (int)Zrow[CPY_LIN_RIGHT];
    if (lid >= n && !CPY_GET_BIT(lvisited, ndid-n)) {
      CPY_SET_BIT(lvisited, ndid-n);
      curNode[k+1] = lid;
      k++;
      continue;
    }
    if (rid >= n && !CPY_GET_BIT(rvisited, ndid-n)) {
      CPY_SET_BIT(rvisited, ndid-n);
      curNode[k+1] = rid;
      k++;
      continue;
    }
    max_dist = Zrow[CPY_LIN_DIST];
    if (lid >= n) {
      max_dist = CPY_MAX(max_dist, max_dists[lid-n]);
    }
    if (rid >= n) {
      max_dist = CPY_MAX(max_dist, max_dists[rid-n]);
    }
    max_dists[ndid-n] = max_dist;
    CPY_DEBUG_MSG("i=%d maxdist[i]=%5.5f verif=%5.5f\n",
		  ndid-n, max_dist, max_dists[ndid-n]);
    k--;
  }
  free(curNode);
  free(lvisited);
  free(rvisited);
}

/**
   Returns the maximum Rrow[rf] field for each cluster node where
   0 <= rf < 3. */

void get_max_Rfield_for_each_cluster(const double *Z, const double *R,
				     double *max_rfs, int n, int rf) {
  int *curNode;
  int ndid, lid, rid, k;
  unsigned char *lvisited, *rvisited;
  const double *Zrow, *Rrow;
  double max_rf;
  const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(n);
  k = 0;
  curNode = (int*)malloc(n * sizeof(int));
  lvisited = (unsigned char*)malloc(bff);
  rvisited = (unsigned char*)malloc(bff);
  curNode[k] = (n * 2) - 2;
  memset(lvisited, 0, bff);
  memset(rvisited, 0, bff);
  while (k >= 0) {
    ndid = curNode[k];
    Zrow = Z + ((ndid-n) * CPY_LIS);
    Rrow = R + ((ndid-n) * CPY_NIS);
    lid = (int)Zrow[CPY_LIN_LEFT];
    rid = (int)Zrow[CPY_LIN_RIGHT];
    if (lid >= n && !CPY_GET_BIT(lvisited, ndid-n)) {
      CPY_SET_BIT(lvisited, ndid-n);
      curNode[k+1] = lid;
      k++;
      continue;
    }
    if (rid >= n && !CPY_GET_BIT(rvisited, ndid-n)) {
      CPY_SET_BIT(rvisited, ndid-n);
      curNode[k+1] = rid;
      k++;
      continue;
    }
    max_rf = Rrow[rf];
    if (lid >= n) {
      max_rf = CPY_MAX(max_rf, max_rfs[lid-n]);
    }
    if (rid >= n) {
      max_rf = CPY_MAX(max_rf, max_rfs[rid-n]);
    }
    max_rfs[ndid-n] = max_rf;
    k--;
  }
  free(curNode);
  free(lvisited);
  free(rvisited);
}

/** find the leaders. report an error if found. */
int leaders(const double *Z, const int *T, int *L, int *M, int kk, int n) {
  int *curNode;
  int ndid, lid, rid, k, nc;
  unsigned char *lvisited, *rvisited;
  const double *Zrow;
  const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(n);
  int *fid; /** done vector, flat cluster ids **/
  int lfid = 0, rfid = 0, errid = -1;

  curNode = (int*)malloc(n * sizeof(int));
  lvisited = (unsigned char*)malloc(bff);
  rvisited = (unsigned char*)malloc(bff);
  fid = (int*)malloc((2 * n - 1) * sizeof(int));

  for (k = 0; k < n; k++) {
    fid[k] = T[k];
  }
  for (k = n; k < 2 * n - 1; k++) {
    fid[k] = -1;
  }

  /** number of clusters formed so far. */
  nc = 0;
  k = 0;
  curNode[k] = (n * 2) - 2;
  memset(lvisited, 0, bff);
  memset(rvisited, 0, bff);
  while (k >= 0) {
    ndid = curNode[k];
    Zrow = Z + ((ndid-n) * CPY_LIS);
    lid = (int)Zrow[CPY_LIN_LEFT];
    rid = (int)Zrow[CPY_LIN_RIGHT];
    CPY_DEBUG_MSG("ndid=%d lid=%d rid=%d\n", ndid, lid, rid);
    if (lid >= n && !CPY_GET_BIT(lvisited, ndid-n)) {
      CPY_SET_BIT(lvisited, ndid-n);
      curNode[k+1] = lid;
      k++;
      continue;
    }
    if (rid >= n && !CPY_GET_BIT(rvisited, ndid-n)) {
      CPY_SET_BIT(rvisited, ndid-n);
      curNode[k+1] = rid;
      k++;
      continue;
    }
    lfid = fid[lid];
    rfid = fid[rid];
    CPY_DEBUG_MSG("[Q] ndid=%d lid=%d lfid=%d rid=%d rfid=%d\n", ndid, lid, lfid, rid, rfid);

    /** If the left and right have the same id, neither can be a leader,
        and their parent takes on their flat cluster id. **/
    if (lfid == rfid) {
      fid[ndid] = lfid;
    }
    /** Otherwise, they are both leaders. */
    else {
      if (lfid != -1) {
	/** If there isn't more room in the result vectors,
	    something is wrong. Condition (2) in help(hcluster.leader)
	    is violated. */
	if (nc >= kk) {
	  errid = ndid;
	  break;
	}
	CPY_DEBUG_MSG("[L] new leader i=%d nc=%d, M[nc]=%d kk=%d n=%d\n", lid, nc, lfid, kk, n);
	L[nc] = lid;
	M[nc] = lfid;
	nc++;
      }
      if (rfid != -1) {
	if (nc >= kk) {
	  errid = ndid;
	  break;
	}
	CPY_DEBUG_MSG("[R] new leader i=%d nc=%d, M[nc]=%d kk=%d n=%d\n", rid, nc, rfid, kk, n);
	L[nc] = rid;
	M[nc] = rfid;
	nc++;
      }
      /** Want to make sure this guy doesn't become a leader since
	  it's children are both leaders. **/
      fid[ndid] = -1;
      
    }
    k--;
  }
  /** For the root node, if its too children have the same flat cluster id,
      neither is negative, the root becomes the leader. */
  Zrow = Z + ((n-2) * CPY_LIS);
  lid = (int)Zrow[CPY_LIN_LEFT];
  rid = (int)Zrow[CPY_LIN_RIGHT];
  lfid = fid[lid];
  rfid = fid[rid];
  if (lfid == rfid && lfid != -1 && errid == -1) {
    if (nc >= kk) {
      errid = (n * 2) - 2;
      /** I know, I know, this looks bad! First time in a good 10 years that I've used one of
	  these. Don't want to copy the free statements. I don't think this detracts from
          the code's readability.*/
      goto leaders_free;
    }
    L[nc] = (n * 2) - 2;
    M[nc] = lfid;
    nc++;
  }
 leaders_free:
  free(curNode);
  free(lvisited);
  free(rvisited);
  free(fid);
  return errid;
}
