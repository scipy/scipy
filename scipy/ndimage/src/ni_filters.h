/* Copyright (C) 2003-2005 Peter J. Verveer
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote
 *    products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef NI_FILTERS_H
#define NI_FILTERS_H

int NI_Correlate1D(PyArrayObject*, PyArrayObject*, int, PyArrayObject*,
                   NI_ExtendMode, double, npy_intp);
int NI_Correlate(PyArrayObject*, PyArrayObject*, PyArrayObject*,
                 NI_ExtendMode, double, npy_intp*);
int NI_UniformFilter1D(PyArrayObject*, npy_intp, int, PyArrayObject*,
                       NI_ExtendMode, double, npy_intp);
int NI_MinOrMaxFilter1D(PyArrayObject*, npy_intp, int, PyArrayObject*,
                        NI_ExtendMode, double, npy_intp, int);
int NI_MinOrMaxFilter(PyArrayObject*, PyArrayObject*, PyArrayObject*,
                      PyArrayObject*, NI_ExtendMode, double, npy_intp*,
                                            int);
int NI_RankFilter(PyArrayObject*, int, PyArrayObject*, PyArrayObject*,
                                    NI_ExtendMode, double, npy_intp*);
int NI_GenericFilter1D(PyArrayObject*, int (*)(double*, npy_intp,
                       double*, npy_intp, void*), void*, npy_intp, int,
                       PyArrayObject*, NI_ExtendMode, double, npy_intp);
int NI_GenericFilter(PyArrayObject*, int (*)(double*, npy_intp, double*,
                                         void*), void*, PyArrayObject*, PyArrayObject*,
                     NI_ExtendMode, double, npy_intp*);
#endif
