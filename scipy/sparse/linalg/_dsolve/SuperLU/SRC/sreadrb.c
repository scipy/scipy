/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file sreadrb.c
 * \brief Read a matrix stored in Rutherford-Boeing format
 *
 * <pre>
 * -- SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Laboratory.
 * June 30, 2009
 * </pre>
 *
 * Purpose
 * =======
 *
 * Read a FLOAT PRECISION matrix stored in Rutherford-Boeing format
 * as described below.
 *
 * Line 1 (A72, A8)
 *      Col. 1 - 72   Title (TITLE)
 *      Col. 73 - 80  Matrix name / identifier (MTRXID)
 *
 * Line 2 (I14, 3(1X, I13))
 *      Col. 1 - 14   Total number of lines excluding header (TOTCRD)
 *      Col. 16 - 28  Number of lines for pointers (PTRCRD)
 *      Col. 30 - 42  Number of lines for row (or variable) indices (INDCRD)
 *      Col. 44 - 56  Number of lines for numerical values (VALCRD)
 *
 * Line 3 (A3, 11X, 4(1X, I13))
 *      Col. 1 - 3    Matrix type (see below) (MXTYPE)
 *      Col. 15 - 28  Compressed Column: Number of rows (NROW)
 *                    Elemental: Largest integer used to index variable (MVAR)
 *      Col. 30 - 42  Compressed Column: Number of columns (NCOL)
 *                    Elemental: Number of element matrices (NELT)
 *      Col. 44 - 56  Compressed Column: Number of entries (NNZERO)
 *                    Elemental: Number of variable indeces (NVARIX)
 *      Col. 58 - 70  Compressed Column: Unused, explicitly zero
 *                    Elemental: Number of elemental matrix entries (NELTVL)
 *
 * Line 4 (2A16, A20)
 *      Col. 1 - 16   Fortran format for pointers (PTRFMT)
 *      Col. 17 - 32  Fortran format for row (or variable) indices (INDFMT)
 *      Col. 33 - 52  Fortran format for numerical values of coefficient matrix
 *                    (VALFMT)
 *                    (blank in the case of matrix patterns)
 *
 * The three character type field on line 3 describes the matrix type.
 * The following table lists the permitted values for each of the three
 * characters. As an example of the type field, RSA denotes that the matrix
 * is real, symmetric, and assembled.
 *
 * First Character:
 *      R Real matrix
 *      C Complex matrix
 *      I integer matrix
 *      P Pattern only (no numerical values supplied)
 *      Q Pattern only (numerical values supplied in associated auxiliary value
 *        file)
 *
 * Second Character:
 *      S Symmetric
 *      U Unsymmetric
 *      H Hermitian
 *      Z Skew symmetric
 *      R Rectangular
 *
 * Third Character:
 *      A Compressed column form
 *      E Elemental form
 *
 * </pre>
 */

#include <stdio.h>
#include <stdlib.h>
#include "slu_sdefs.h"


/*! \brief Eat up the rest of the current line */
static int sDumpLine(FILE *fp)
{
    register int c;
    while ((c = fgetc(fp)) != '\n') ;
    return 0;
}

static int sParseIntFormat(char *buf, int *num, int *size)
{
    char *tmp;

    tmp = buf;
    while (*tmp++ != '(') ;
    sscanf(tmp, "%d", num);
    while (*tmp != 'I' && *tmp != 'i') ++tmp;
    ++tmp;
    sscanf(tmp, "%d", size);
    return 0;
}

static int sParseFloatFormat(char *buf, int *num, int *size)
{
    char *tmp, *period;

    tmp = buf;
    while (*tmp++ != '(') ;
    *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
    while (*tmp != 'E' && *tmp != 'e' && *tmp != 'D' && *tmp != 'd'
           && *tmp != 'F' && *tmp != 'f') {
        /* May find kP before nE/nD/nF, like (1P6F13.6). In this case the
           num picked up refers to P, which should be skipped. */
        if (*tmp=='p' || *tmp=='P') {
           ++tmp;
           *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
        } else {
           ++tmp;
        }
    }
    ++tmp;
    period = tmp;
    while (*period != '.' && *period != ')') ++period ;
    *period = '\0';
    *size = atoi(tmp); /*sscanf(tmp, "%2d", size);*/

    return 0;
}

static int ReadVector(FILE *fp, int n, int *where, int perline, int persize)
{
    register int i, j, item;
    char tmp, buf[100], *dummy;

    i = 0;
    while (i < n) {
        dummy = fgets(buf, 100, fp);    /* read a line at a time */
        if(dummy == NULL) {
            ABORT("Unable to read from the file");
        }
        for (j=0; j<perline && i<n; j++) {
            tmp = buf[(j+1)*persize];     /* save the char at that place */
            buf[(j+1)*persize] = 0;       /* null terminate */
            item = atoi(&buf[j*persize]);
            buf[(j+1)*persize] = tmp;     /* recover the char at that place */
            where[i++] = item - 1;
        }
    }

    return 0;
}

static int sReadValues(FILE *fp, int n, float *destination, int perline,
        int persize)
{
    register int i, j, k, s;
    char tmp, buf[100], *dummy;

    i = 0;
    while (i < n) {
        dummy = fgets(buf, 100, fp);    /* read a line at a time */
        if(dummy == NULL) {
            ABORT("Unable to read from the file");
        }
        for (j=0; j<perline && i<n; j++) {
            tmp = buf[(j+1)*persize];     /* save the char at that place */
            buf[(j+1)*persize] = 0;       /* null terminate */
            s = j*persize;
            for (k = 0; k < persize; ++k) /* No D_ format in C */
                if ( buf[s+k] == 'D' || buf[s+k] == 'd' ) buf[s+k] = 'E';
            destination[i++] = atof(&buf[s]);
            buf[(j+1)*persize] = tmp;     /* recover the char at that place */
        }
    }

    return 0;
}



/*! \brief
 *
 * <pre>
 * On input, nonz/nzval/rowind/colptr represents lower part of a symmetric
 * matrix. On exit, it represents the full matrix with lower and upper parts.
 * </pre>
 */
static void
FormFullA(int n, int *nonz, float **nzval, int **rowind, int **colptr)
{
    register int i, j, k, col, new_nnz;
    int *t_rowind, *t_colptr, *al_rowind, *al_colptr, *a_rowind, *a_colptr;
    int *marker;
    float *t_val, *al_val, *a_val;

    al_rowind = *rowind;
    al_colptr = *colptr;
    al_val = *nzval;

    if ( !(marker =(int *) SUPERLU_MALLOC( (n+1) * sizeof(int)) ) )
	ABORT("SUPERLU_MALLOC fails for marker[]");
    if ( !(t_colptr = (int *) SUPERLU_MALLOC( (n+1) * sizeof(int)) ) )
	ABORT("SUPERLU_MALLOC t_colptr[]");
    if ( !(t_rowind = (int *) SUPERLU_MALLOC( *nonz * sizeof(int)) ) )
	ABORT("SUPERLU_MALLOC fails for t_rowind[]");
    if ( !(t_val = (float*) SUPERLU_MALLOC( *nonz * sizeof(float)) ) )
	ABORT("SUPERLU_MALLOC fails for t_val[]");

    /* Get counts of each column of T, and set up column pointers */
    for (i = 0; i < n; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = al_colptr[j]; i < al_colptr[j+1]; ++i)
	    ++marker[al_rowind[i]];
    }
    t_colptr[0] = 0;
    for (i = 0; i < n; ++i) {
	t_colptr[i+1] = t_colptr[i] + marker[i];
	marker[i] = t_colptr[i];
    }

    /* Transpose matrix A to T */
    for (j = 0; j < n; ++j)
	for (i = al_colptr[j]; i < al_colptr[j+1]; ++i) {
	    col = al_rowind[i];
	    t_rowind[marker[col]] = j;
	    t_val[marker[col]] = al_val[i];
	    ++marker[col];
	}

    new_nnz = *nonz * 2 - n;
    if ( !(a_colptr = (int *) SUPERLU_MALLOC( (n+1) * sizeof(int)) ) )
	ABORT("SUPERLU_MALLOC a_colptr[]");
    if ( !(a_rowind = (int *) SUPERLU_MALLOC( new_nnz * sizeof(int)) ) )
	ABORT("SUPERLU_MALLOC fails for a_rowind[]");
    if ( !(a_val = (float*) SUPERLU_MALLOC( new_nnz * sizeof(float)) ) )
	ABORT("SUPERLU_MALLOC fails for a_val[]");

    a_colptr[0] = 0;
    k = 0;
    for (j = 0; j < n; ++j) {
      for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
	if ( t_rowind[i] != j ) { /* not diagonal */
	  a_rowind[k] = t_rowind[i];
	  a_val[k] = t_val[i];
#ifdef DEBUG
	  if ( fabs(a_val[k]) < 4.047e-300 )
	      printf("%5d: %e\n", k, a_val[k]);
#endif
	  ++k;
	}
      }

      for (i = al_colptr[j]; i < al_colptr[j+1]; ++i) {
	a_rowind[k] = al_rowind[i];
	a_val[k] = al_val[i];
#ifdef DEBUG
	if ( fabs(a_val[k]) < 4.047e-300 )
	    printf("%5d: %e\n", k, a_val[k]);
#endif
	++k;
      }

      a_colptr[j+1] = k;
    }

    printf("FormFullA: new_nnz = %d, k = %d\n", new_nnz, k);

    SUPERLU_FREE(al_val);
    SUPERLU_FREE(al_rowind);
    SUPERLU_FREE(al_colptr);
    SUPERLU_FREE(marker);
    SUPERLU_FREE(t_val);
    SUPERLU_FREE(t_rowind);
    SUPERLU_FREE(t_colptr);

    *nzval = a_val;
    *rowind = a_rowind;
    *colptr = a_colptr;
    *nonz = new_nnz;
}


void
sreadrb(int *nrow, int *ncol, int *nonz,
        float **nzval, int **rowind, int **colptr)
{

    register int i, numer_lines = 0;
    int tmp, colnum, colsize, rownum, rowsize, valnum, valsize;
    char buf[100], type[4], *dummy;
    int sym, f_count = 0, s_count = 0;
    FILE *fp;

    fp = stdin;

    /* Line 1 */
    dummy = fgets(buf, 100, fp);
    if(dummy == NULL) {
        ABORT("Unable to read from the file");
    }
    fputs(buf, stdout);

    /* Line 2 */
    for (i=0; i<4; i++) {
        f_count = fscanf(fp, "%14c", buf); buf[14] = 0;
        check_read(f_count);
        s_count = sscanf(buf, "%d", &tmp);
        check_read(s_count);
        if (i == 3) numer_lines = tmp;
    }
    sDumpLine(fp);

    /* Line 3 */
    f_count = fscanf(fp, "%3c", type);
    check_read(f_count);
    f_count = fscanf(fp, "%11c", buf); /* pad */
    check_read(f_count);
    type[3] = 0;
#ifdef DEBUG
    printf("Matrix type %s\n", type);
#endif

    f_count = fscanf(fp, "%14c", buf);
    check_read(f_count);
    s_count = sscanf(buf, "%d", nrow);
    check_read(s_count);
    f_count = fscanf(fp, "%14c", buf);
    check_read(f_count);
    s_count = sscanf(buf, "%d", ncol);
    check_read(s_count);
    f_count = fscanf(fp, "%14c", buf);
    check_read(f_count);
    s_count = sscanf(buf, "%d", nonz);
    check_read(s_count);
    f_count = fscanf(fp, "%14c", buf);
    check_read(f_count);
    s_count = sscanf(buf, "%d", &tmp);
    check_read(s_count);

    if (tmp != 0)
        printf("This is not an assembled matrix!\n");
    if (*nrow != *ncol)
        printf("Matrix is not square.\n");
    sDumpLine(fp);

    /* Allocate storage for the three arrays ( nzval, rowind, colptr ) */
    sallocateA(*ncol, *nonz, nzval, rowind, colptr);

    /* Line 4: format statement */
    f_count = fscanf(fp, "%16c", buf);
    check_read(f_count);
    sParseIntFormat(buf, &colnum, &colsize);
    f_count = fscanf(fp, "%16c", buf);
    check_read(f_count);
    sParseIntFormat(buf, &rownum, &rowsize);
    f_count = fscanf(fp, "%20c", buf);
    check_read(f_count);
    sParseFloatFormat(buf, &valnum, &valsize);
    sDumpLine(fp);

#ifdef DEBUG
    printf("%d rows, %d nonzeros\n", *nrow, *nonz);
    printf("colnum %d, colsize %d\n", colnum, colsize);
    printf("rownum %d, rowsize %d\n", rownum, rowsize);
    printf("valnum %d, valsize %d\n", valnum, valsize);
#endif

    ReadVector(fp, *ncol+1, *colptr, colnum, colsize);
    ReadVector(fp, *nonz, *rowind, rownum, rowsize);
    if ( numer_lines ) {
        sReadValues(fp, *nonz, *nzval, valnum, valsize);
    }

    sym = (type[1] == 'S' || type[1] == 's');
    if ( sym ) {
	FormFullA(*ncol, nonz, nzval, rowind, colptr);
    }

    fclose(fp);
}
