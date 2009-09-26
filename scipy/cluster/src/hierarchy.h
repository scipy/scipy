/**
 * hierarchy.h
 *
 * Author: Damian Eads
 * Date:   September 22, 2007
 * Adapted for incorporation into Scipy, April 9, 2008.
 *
 * Copyright (c) 2007, 2008, Damian Eads. All rights reserved.
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

#ifndef _CPY_HIERARCHY_H
#define _CPY_HIERARCHY_H

#define CPY_LINKAGE_SINGLE 0
#define CPY_LINKAGE_COMPLETE 1
#define CPY_LINKAGE_AVERAGE 2
#define CPY_LINKAGE_CENTROID 3
#define CPY_LINKAGE_MEDIAN 4
#define CPY_LINKAGE_WARD 5
#define CPY_LINKAGE_WEIGHTED 6

#define CPY_CRIT_INCONSISTENT 0
#define CPY_CRIT_DISTANCE 1
#define CPY_CRIT_MAXCLUST 2

typedef struct cnode {
  int n;
  int id;
  double d;
  struct cnode *left;
  struct cnode *right;
} cnode;

typedef struct clnode {
  struct clnode *next;
  struct cnode *val;
} clnode;

typedef struct clist {
  struct clnode *head;
  struct clnode *tail;
} clist;

typedef struct cinfo {
  struct cnode *nodes;
  struct clist *lists;
  int *ind;
  double *dmt;
  double *dm;
  double *buf;
  double **rows;
  double **centroids;
  double *centroidBuffer;
  const double *X;
  int *rowsize;
  int m;
  int n;
  int nid;
} cinfo;

typedef void (distfunc) (cinfo *info, int mini, int minj, int np, int n); 

void inconsistency_calculation(const double *Z, double *R, int n, int d);
void inconsistency_calculation_alt(const double *Z, double *R, int n, int d);

void chopmins(int *ind, int mini, int minj, int np);
void chopmins_ns_i(double *ind, int mini, int np);
void chopmins_ns_ij(double *ind, int mini, int minj, int np);

void dist_single(cinfo *info, int mini, int minj, int np, int n);
void dist_average(cinfo *info, int mini, int minj, int np, int n);
void dist_complete(cinfo *info, int mini, int minj, int np, int n);
void dist_centroid(cinfo *info, int mini, int minj, int np, int n);
void dist_ward(cinfo *info, int mini, int minj, int np, int n);
void dist_weighted(cinfo *info, int mini, int minj, int np, int n);

int leaders(const double *Z, const int *T, int *L, int *M, int kk, int n);

void linkage(double *dm, double *Z, double *X, int m, int n, int ml, int kc, distfunc dfunc, int method);
void linkage_alt(double *dm, double *Z, double *X, int m, int n, int ml, int kc, distfunc dfunc, int method);

void cophenetic_distances(const double *Z, double *d, int n);
void cpy_to_tree(const double *Z, cnode **tnodes, int n);
void calculate_cluster_sizes(const double *Z, double *cs, int n);

void form_member_list(const double *Z, int *members, int n);
void form_flat_clusters_from_in(const double *Z, const double *R, int *T,
				double cutoff, int n);
void form_flat_clusters_from_dist(const double *Z, int *T,
				  double cutoff, int n);
void form_flat_clusters_from_monotonic_criterion(const double *Z,
						 const double *mono_crit,
						 int *T, double cutoff, int n);

void form_flat_clusters_maxclust_dist(const double *Z, int *T, int n, int mc);

void form_flat_clusters_maxclust_monocrit(const double *Z,
					  const double *mono_crit,
					  int *T, int n, int mc);

void get_max_dist_for_each_cluster(const double *Z, double *max_dists, int n);
void get_max_Rfield_for_each_cluster(const double *Z, const double *R,
				     double *max_rfs, int n, int rf);
#endif
