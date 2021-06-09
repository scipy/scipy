/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: matrix_source.h                                                   *
 *                                                                           *
 *   Routines for computations with square matrices.                         *
 *   Matrices are stored as double array (i.e.: double *matrix).             *
 *   The rows of the matrix have to be stored consecutively in this array,   *
 *   i.e. the entry with index [i,j] can be entered via (i*dim+j), where     *
 *   dim is the dimension of the dim x dim - matrix.                         *
 *                                                                           *
 *   Routines are mainly adapted from the GSL (GNU Scientifiy Library)       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   adapted by Wolfgang Hoermann and Josef Leydold                          *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/
				                                                                                    
/*---------------------------------------------------------------------------*/

int _unur_matrix_transform_diagonal (int dim, const double *M, const double *D, double *res);
/* Computes the transformation M^t . D . M of diagonal matrix D              */
/* and stores it in matrix 'res'.                                            */

int _unur_matrix_cholesky_decomposition (int dim, const double *S, double *L );
/* The Cholesky factor L of a variance-covariance matrix S is computed:      */
/*   S = LL'                                                                 */

int _unur_matrix_invert_matrix (int dim, const double *A, double *Ainv, double *det );
/* Calculates the inverse matrix (by means of LU decomposition).             */
/* As a side effect det(A) is computed.                                      */

double _unur_matrix_determinant ( int dim, const double *A );
/* Calculates the determinant of the matrix A                                */

int _unur_matrix_multiplication(int dim, const double *A, const double *B, double *AB);
/* Calculates the matrix multiplication of two matrices A and B              */

double _unur_matrix_qf (int dim, double *x, double *A );
/* Compute quadratic form x'Ax.                                              */

int _unur_matrix_eigensystem (int dim, const double *M, double *values, double *vectors );
/* Calculates eigenvalues and eigenvectors of real symmetric matrix M.       */
/* The eigenvectors are normalized and (almost) orthognal.                   */
/* The eigenvectors are stored consecutively in the array vectors.           */

void _unur_matrix_print_vector ( int dim, const double *vec, const char *info,
				 FILE *LOG, const char *genid, const char *indent );
/* Print elements of vector in a single row enclosed by parenthesis into     */
/* LOG file. The line starts with <genid>: <indent>                          */
/* <info> is printed in a (first) separate line.                             */
/* A blank line is inserted after the printed vector.                        */
/* If the NULL pointer is given, the string "[unknown]" is printed.          */

void _unur_matrix_print_matrix ( int dim, const double *mat, const char *info,
				 FILE *LOG, const char *genid, const char *indent );
/* Print elements of the given <dim>x<dim> square matrix into LOG file.      */
/* The matrix is stored row-wise in <mat>.                                   */
/* The lines start with <genid>: <indent>                                    */
/* <info> is printed in a (first) separate line.                             */
/* A blank line is inserted after the printed matrix.                        */
/* If the NULL pointer is given, the string "[unknown]" is printed.          */

/*---------------------------------------------------------------------------*/
