/*
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*
 -- SuperLU routine (version 4.4) --
 Lawrence Berkeley National Lab, Univ. of California Berkeley.
 November 20, 2012
*/

#include <stdio.h>
#include "slu_Cnames.h"

/*! @file input_error.c
 * \brief Error handler for input parameters.
 */

/*! \brief Error handler for input parameters.
 *
 * It is called if an input parameter has an
 * invalid value.
 * \param[in] srname Six characters, the name of the routine which called input_error.
 * \param[in] info   The position of the invalid parameter in the parameter list
 *                   of the calling routine.
 */
void input_error(char *srname, int *info)
{
    printf("** On entry to %6s, parameter number %2d had an illegal value\n",
		srname, *info);
}
