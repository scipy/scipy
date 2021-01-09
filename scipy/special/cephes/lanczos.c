/*  (C) Copyright John Maddock 2006.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)
 */

/* Scipy changes:
 * - 06-22-2016: Removed all code not related to double precision and
 *   ported to c for use in Cephes
 */

#include "mconf.h"
#include "lanczos.h"


static double lanczos_sum(double x)
{
    return ratevl(x, lanczos_num,
		  sizeof(lanczos_num) / sizeof(lanczos_num[0]) - 1,
		  lanczos_denom,
		  sizeof(lanczos_denom) / sizeof(lanczos_denom[0]) - 1);
}


double lanczos_sum_expg_scaled(double x)
{
    return ratevl(x, lanczos_sum_expg_scaled_num,
		  sizeof(lanczos_sum_expg_scaled_num) / sizeof(lanczos_sum_expg_scaled_num[0]) - 1,
		  lanczos_sum_expg_scaled_denom,
		  sizeof(lanczos_sum_expg_scaled_denom) / sizeof(lanczos_sum_expg_scaled_denom[0]) - 1);
}


static double lanczos_sum_near_1(double dx)
{
    double result = 0;
    unsigned k;
    
    for (k = 1; k <= sizeof(lanczos_sum_near_1_d)/sizeof(lanczos_sum_near_1_d[0]); ++k) {
	result += (-lanczos_sum_near_1_d[k-1]*dx)/(k*dx + k*k);
    }
    return result;
}


static double lanczos_sum_near_2(double dx)
{
    double result = 0;
    double x = dx + 2;
    unsigned k;
    
    for(k = 1; k <= sizeof(lanczos_sum_near_2_d)/sizeof(lanczos_sum_near_2_d[0]); ++k) {
	result += (-lanczos_sum_near_2_d[k-1]*dx)/(x + k*x + k*k - 1);
    }
    return result;
}
