/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file 
 * \brief Read a matrix stored in Harwell-Boeing format.
 * Contributed by Francois-Henry Rouet.
 *
 */
#include <ctype.h>
#include "slu_ddefs.h"

#undef EXPAND_SYM

/*! brief
 *
 * <pre>
 * Output parameters
 * =================
 *   (nzval, rowind, colptr): (*rowind)[*] contains the row subscripts of
 *      nonzeros in columns of matrix A; (*nzval)[*] the numerical values;
 *	column i of A is given by (*nzval)[k], k = (*rowind)[i],...,
 *      (*rowind)[i+1]-1.
 * </pre>
 */

void
dreadMM(FILE *fp, int *m, int *n, int_t *nonz,
	    double **nzval, int_t **rowind, int_t **colptr)
{
    int_t    j, k, jsize, nnz, nz, new_nonz;
    double *a, *val;
    int_t    *asub, *xa;
    int      *row, *col;
    int    zero_base = 0;
    char *p, line[512], banner[64], mtx[64], crd[64], arith[64], sym[64];
    int expand;

    /* 	File format:
     *    %%MatrixMarket matrix coordinate real general/symmetric/...
     *    % ...
     *    % (optional comments)
     *    % ...
     *    #rows    #non-zero
     *    Triplet in the rest of lines: row    col    value
     */

     /* 1/ read header */ 
     fgets(line,512,fp);
     for (p=line; *p!='\0'; *p=tolower(*p),p++);

     if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, arith, sym) != 5) {
       printf("Invalid header (first line does not contain 5 tokens)\n");
       exit(-1);
     }
 
     if(strcmp(banner,"%%matrixmarket")) {
       printf("Invalid header (first token is not \"%%%%MatrixMarket\")\n");
       exit(-1);
     }

     if(strcmp(mtx,"matrix")) {
       printf("Not a matrix; this driver cannot handle that.\n");
       exit(-1);
     }

     if(strcmp(crd,"coordinate")) {
       printf("Not in coordinate format; this driver cannot handle that.\n");
       exit(-1);
     }

     if(strcmp(arith,"real")) {
       if(!strcmp(arith,"complex")) {
         printf("Complex matrix; use zreadMM instead!\n");
         exit(-1);
       }
       else if(!strcmp(arith, "pattern")) {
         printf("Pattern matrix; values are needed!\n");
         exit(-1);
       }
       else {
         printf("Unknown arithmetic\n");
         exit(-1);
       }
     }

     if(strcmp(sym,"general")) {
       printf("Symmetric matrix: will be expanded\n");
       expand=1;
     } else expand=0;

     /* 2/ Skip comments */
     while(banner[0]=='%') {
       fgets(line,512,fp);
       sscanf(line,"%s",banner);
     }

     /* 3/ Read n and nnz */
#ifdef _LONGINT
    sscanf(line, "%d%d%lld",m, n, nonz);
#else
    sscanf(line, "%d%d%d",m, n, nonz);
#endif
    printf("m %lld, n %lld, nonz %lld\n", (long long) *m, (long long) *n, (long long) *nonz);
    if(*m!=*n) {
      printf("Rectangular matrix!. Abort\n");
      exit(-1);
   }

    if(expand)
      new_nonz = 2 * *nonz - *n;
    else
      new_nonz = *nonz;


    dallocateA(*n, new_nonz, nzval, rowind, colptr); /* Allocate storage */
    a    = *nzval;
    asub = *rowind;
    xa   = *colptr;

    if ( !(val = (double *) SUPERLU_MALLOC(new_nonz * sizeof(double))) )
        ABORT("Malloc fails for val[]");
    if ( !(row = int32Malloc(new_nonz)) )
        ABORT("Malloc fails for row[]");
    if ( !(col = int32Malloc(new_nonz)) )
        ABORT("Malloc fails for col[]");

    for (j = 0; j < *n; ++j) xa[j] = 0;

    /* 4/ Read triplets of values */
    for (nnz = 0, nz = 0; nnz < *nonz; ++nnz) {
	fscanf(fp, "%d%d%lf\n", &row[nz], &col[nz], &val[nz]);

	if ( nnz == 0 ) { /* first nonzero */
	    if ( row[0] == 0 || col[0] == 0 ) {
		zero_base = 1;
		printf("triplet file: row/col indices are zero-based.\n");
	    } else {
		printf("triplet file: row/col indices are one-based.\n");
            }
	}

	if ( !zero_base ) {
	    /* Change to 0-based indexing. */
	    --row[nz];
	    --col[nz];
	}

	if (row[nz] < 0 || row[nz] >= *m || col[nz] < 0 || col[nz] >= *n
	    /*|| val[nz] == 0.*/) {
	    fprintf(stderr, "nz %d, (%d, %d) = %e out of bound, removed\n", 
                    (int) nz, row[nz], col[nz], val[nz]);
	    exit(-1);
	} else {
	    ++xa[col[nz]];
            if(expand) {
	        if ( row[nz] != col[nz] ) { /* Excluding diagonal */
	          ++nz;
	          row[nz] = col[nz-1];
	          col[nz] = row[nz-1];
	          val[nz] = val[nz-1];
	          ++xa[col[nz]];
	        }
            }	
	    ++nz;
	}
    }

    *nonz = nz;
    if(expand) {
      printf("new_nonz after symmetric expansion:\t%lld\n", (long long)*nonz);
    }
    

    /* Initialize the array of column pointers */
    k = 0;
    jsize = xa[0];
    xa[0] = 0;
    for (j = 1; j < *n; ++j) {
	k += jsize;
	jsize = xa[j];
	xa[j] = k;
    }
    
    /* Copy the triplets into the column oriented storage */
    for (nz = 0; nz < *nonz; ++nz) {
	j = col[nz];
	k = xa[j];
	asub[k] = row[nz];
	a[k] = val[nz];
	++xa[j];
    }

    /* Reset the column pointers to the beginning of each column */
    for (j = *n; j > 0; --j)
	xa[j] = xa[j-1];
    xa[0] = 0;

    SUPERLU_FREE(val);
    SUPERLU_FREE(row);
    SUPERLU_FREE(col);

#ifdef CHK_INPUT
    int i;
    for (i = 0; i < *n; i++) {
	printf("Col %d, xa %d\n", i, xa[i]);
	for (k = xa[i]; k < xa[i+1]; k++)
	    printf("%d\t%16.10f\n", asub[k], a[k]);
    }
#endif

}


static void dreadrhs(int m, double *b)
{
    FILE *fp = fopen("b.dat", "r");
    int i;

    if ( !fp ) {
        fprintf(stderr, "dreadrhs: file does not exist\n");
	exit(-1);
    }
    for (i = 0; i < m; ++i)
      fscanf(fp, "%lf\n", &b[i]);

    fclose(fp);
}
