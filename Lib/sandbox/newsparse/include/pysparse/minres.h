/**************************************************************************
*                                                                         *
*               Swiss Federal Institute of Technology (ETH),              *
*                       CH-8092 Zuerich, Switzerland                      *
*                                                                         *
*                       (C) 1999 All Rights Reserved                      *
*                                                                         *
*                                NOTICE                                   *
*                                                                         *
*  Permission to use, copy, modify, and distribute this software and      *
*  its documentation for any purpose and without fee is hereby granted    *
*  provided that the above copyright notice appear in all copies and      *
*  that both the copyright notice and this permission notice appear in    *
*  supporting documentation.                                              *
*                                                                         *
*  Neither the Swiss Federal Institute of Technology nor the author make  *
*  any representations about the suitability of this software for any     *
*  purpose.  This software is provided ``as is'' without express or       *
*  implied warranty.                                                      *
*                                                                         *
**************************************************************************/

#ifndef MINRES_H
#define MINRES_H

#include "Python.h"

int Itsolvers_minres_kernel(int n, double errtol, int it_max,
			    int *it, double *nrm_res, int clvl,
			    double *x, double *b, double *work,
			    PyObject *mat_obj,
			    PyObject *prec_obj);

#endif
