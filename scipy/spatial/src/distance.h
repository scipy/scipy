/**
 * distance.h
 *
 * Author: Damian Eads
 * Date:   September 22, 2007 (moved to new file on June 8, 2008)
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

#ifndef _CPY_DISTANCE_H
#define _CPY_DISTANCE_H

void dist_to_squareform_from_vector(double *M, const double *v, int n);
void dist_to_vector_from_squareform(const double *M, double *v, int n);
void pdist_euclidean(const double *X, double *dm, int m, int n);
void pdist_seuclidean(const double *X,
		      const double *var, double *dm, int m, int n);
void pdist_mahalanobis(const double *X, const double *covinv,
		       double *dm, int m, int n);
void pdist_bray_curtis(const double *X, double *dm, int m, int n);
void pdist_canberra(const double *X, double *dm, int m, int n);
void pdist_hamming(const double *X, double *dm, int m, int n);
void pdist_hamming_bool(const char *X, double *dm, int m, int n);
void pdist_city_block(const double *X, double *dm, int m, int n);
void pdist_cosine(const double *X, double *dm, int m, int n, const double *norms);
void pdist_chebyshev(const double *X, double *dm, int m, int n);
void pdist_jaccard(const double *X, double *dm, int m, int n);
void pdist_jaccard_bool(const char *X, double *dm, int m, int n);
void pdist_kulsinski_bool(const char *X, double *dm, int m, int n);
void pdist_minkowski(const double *X, double *dm, int m, int n, double p);
void pdist_weighted_minkowski(const double *X, double *dm, int m, int n, double p, const double *w);
void pdist_yule_bool(const char *X, double *dm, int m, int n);
void pdist_matching_bool(const char *X, double *dm, int m, int n);
void pdist_dice_bool(const char *X, double *dm, int m, int n);
void pdist_rogerstanimoto_bool(const char *X, double *dm, int m, int n);
void pdist_russellrao_bool(const char *X, double *dm, int m, int n);
void pdist_sokalmichener_bool(const char *X, double *dm, int m, int n);
void pdist_sokalsneath_bool(const char *X, double *dm, int m, int n);

void cdist_euclidean(const double *XA, const double *XB, double *dm, int mA, int mB, int n);
void cdist_mahalanobis(const double *XA, const double *XB,
		       const double *covinv,
		       double *dm, int mA, int mB, int n);
void cdist_bray_curtis(const double *XA, const double *XB,
		       double *dm, int mA, int mB, int n);
void cdist_canberra(const double *XA,
		    const double *XB, double *dm, int mA, int mB, int n);
void cdist_hamming(const double *XA,
		   const double *XB, double *dm, int mA, int mB, int n);
void cdist_hamming_bool(const char *XA,
			const char *XB, double *dm,
			int mA, int mB, int n);
void cdist_jaccard(const double *XA,
		   const double *XB, double *dm, int mA, int mB, int n);
void cdist_jaccard_bool(const char *XA,
			const char *XB, double *dm, int mA, int mB, int n);
void cdist_chebyshev(const double *XA,
		     const double *XB, double *dm, int mA, int mB, int n);
void cdist_cosine(const double *XA,
		  const double *XB, double *dm, int mA, int mB, int n,
		  const double *normsA, const double *normsB);
void cdist_seuclidean(const double *XA,
		      const double *XB,
		      const double *var,
		      double *dm, int mA, int mB, int n);
void cdist_city_block(const double *XA, const double *XB, double *dm,
		      int mA, int mB, int n);
void cdist_minkowski(const double *XA, const double *XB, double *dm,
		     int mA, int mB, int n, double p);
void cdist_weighted_minkowski(const double *XA, const double *XB, double *dm,
			      int mA, int mB, int n, double p, const double *w);
void cdist_yule_bool(const char *XA, const char *XB, double *dm,
		     int mA, int mB, int n);
void cdist_matching_bool(const char *XA, const char *XB, double *dm,
			 int mA, int mB, int n);
void cdist_dice_bool(const char *XA, const char *XB, double *dm,
		     int mA, int mB, int n);
void cdist_rogerstanimoto_bool(const char *XA, const char *XB, double *dm,
			       int mA, int mB, int n);
void cdist_russellrao_bool(const char *XA, const char *XB, double *dm,
			   int mA, int mB, int n);
void cdist_kulsinski_bool(const char *XA, const char *XB, double *dm,
			  int mA, int mB, int n);
void cdist_sokalsneath_bool(const char *XA, const char *XB, double *dm,
			    int mA, int mB, int n);
void cdist_sokalmichener_bool(const char *XA, const char *XB, double *dm,
			      int mA, int mB, int n);

#endif
