#include "Python.h"
#include <math.h>
#include "newsparse/mmio.h"
#include "newsparse/ll_mat.h"
#include "newsparse/csr_mat.h"
#include "newsparse/sss_mat.h"

#define SPMATRIX_MODULE
#include "newsparse/spmatrix.h"

#define PY_ARRAY_UNIQUE_SYMBOL spmatrix
/* Was: #include "Numeric/arrayobject.h" */
#include "scipy/arrayobject.h"

#define INCREASE_FACTOR   1.5	/* increase rate for memory reallocation of ll_mat arrays */
#define PPRINT_ROW_THRESH 500	/* row threshold for choosing between print formats */
#define PPRINT_COL_THRESH 20	/* column threshold for choosing between print formats */
#define OPT_MATMATMUL     1	/* enable optimised matrix-matrix-mult */

#define MAX(A,B) ( (A) > (B) ? (A) : (B) )


/*************************************************************************************/
/*  R o u t i n e s   f o r   b u i l d i n g   d a t a   s t r u c t u r e   f o r  */
/*  c o l u m n - w i s e   t r a v e r s a l   o f   l l _ m a t   o b j e c t s    */
/*************************************************************************************/

static int 
SpMatrix_LLMatBuildColIndex(struct llColIndex **idx, LLMatObject *self, int includeDiagonal) {
/*
 *  build data structure for column-wise traversal
 *
 *    build a (root,row,link) linked-list data structure which links the
 *    entries of each column of "self".
 *
 *    if "includeDiagonal" is zero the diagonal elements of "self" are
 *    not included in the linked-list data structure.
 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */
  int i, j, k;

  if(!(*idx = (struct llColIndex *)malloc(sizeof(struct llColIndex))))
    goto fail;
  
  /* allocate (link,row,root) arrays
   *
   *   can we do better than nalloc ???
   */
  (*idx)->link = PyMem_New(int, self->nalloc);
  if ((*idx)->link == NULL)
    goto fail;
  (*idx)->row = PyMem_New(int, self->nalloc);
  if ((*idx)->row == NULL)
    goto fail;
  (*idx)->root = PyMem_New(int, self->dim[1]);
  if ((*idx)->root == NULL)
    goto fail;

  /* initialize root arrays
   */
  for (i = 0; i < self->dim[1]; i ++)
    (*idx)->root[i] = -1;

  /* scan matrix from bottom up
     so the resulting lists will be sorted by ascending row*/
  (*idx)->nzLo = 0; (*idx)->nzDiag = 0; (*idx)->nzUp = 0;
  for (i = self->dim[0] - 1; i >= 0; i --) {
    k = self->root[i];
    while (k != -1) {
      j = self->col[k];
      if (i > j)
	(*idx)->nzLo ++;
      else if (i == j)
	(*idx)->nzDiag ++;
      else
	(*idx)->nzUp ++;
      if (includeDiagonal || i != j) {
	(*idx)->link[k] = (*idx)->root[j];
	(*idx)->root[j] = k;
	(*idx)->row[k] = i;
      }
      k = self->link[k];
    }
  }
  return 0;

 fail:
  if (*idx != NULL) {
    PyMem_Del((*idx)->link);    
    PyMem_Del((*idx)->row);    
    PyMem_Del((*idx)->root);
    free(*idx);
    *idx = NULL;
  }
  PyErr_NoMemory();
  return 1;
}

static void
SpMatrix_LLMatDestroyColIndex(struct llColIndex **idx) {
/*
 *  free memory of (root,row,link) linked-list data
 *
* * * * * * * * * * * * * * * * * * * * * * * * * *
 */
  if (*idx != NULL) {
    PyMem_Del((*idx)->link);    
    PyMem_Del((*idx)->row); 
    PyMem_Del((*idx)->root);
    free(*idx);
    *idx = NULL;
  }
}  


/*******************************************************************************/
/*  R o u t i n e s   f o r   s e t t i n g ,   u p d a t i n g                */
/*  a n d   r e a d i n g   e n t r i e s   o f   l l _ m a t   o b j e c t s  */
/*******************************************************************************/

static double
SpMatrix_LLMatGetItem(LLMatObject *a, int i, int j) {
/*
 *  return matrix entry a[i,j] as a double value
 *
* * * * * * * * * * * * * * * * * * * * * * * * *
 */
  int k, t;

  if (a->issym && i < j) {
    t = i; i = j; j = t;
  }

  k = a->root[i];
  while (k != -1) {
    if (a->col[k] == j)
      return a->val[k];
    k = a->link[k];
  }
  return 0.0;
}

static int 
SpMatrix_LLMatSetItem(LLMatObject *a, int i, int j, double x) {
/*
 *  set matrix entry: a[i,j] = x
 *
* * * * * * * * * * * * * * * * *
 */
  void *temp;
  int k, new_elem, last, col;

  if (a->issym && i < j) {
    PyErr_SetString(SpMatrix_ErrorObject, 
		    "write operation to upper triangle of symmetric matrix");
    return -1;
  } 

  /* find element to be set (or removed) */
  col = last = -1;
  k = a->root[i]; 
  while (k != -1) {
    col = a->col[k];
    if (col >= j)
      break;
    last = k;
    k = a->link[k];
  }

  if (x != 0.0) {
    
    if (col == j) {

      /* element already exist */
      a->val[k] = x;

    } else {
      /* new element */

      /* find location for new element */
      if (a->free != -1) {
	
	/* use element from the free chain */
	new_elem = a->free;
	a->free = a->link[new_elem];
	
      } else {
	
	/* append new element to the end */
	new_elem = a->nnz;
	
	/* test if there is space for a new element */
	if (a->nnz == a->nalloc) {
	  int nalloc_new;
	  
	  /* increase size of idx, val and link arrays */
	  nalloc_new = (int)((double)INCREASE_FACTOR * a->nalloc) + 1;
	  if ((temp = PyMem_Resize(a->col, int, nalloc_new)) == NULL)
	    goto fail;
	  else
	    a->col = temp;
	  if ((temp = PyMem_Resize(a->link, int, nalloc_new)) == NULL)
	    goto fail;
	  else
	    a->link = temp;
	  if ((temp = PyMem_Resize(a->val, double, nalloc_new)) == NULL)
	    goto fail;
	  else
	    a->val = temp;
	  a->nalloc = nalloc_new;
	}
      }

      a->val[new_elem] = x;
      a->col[new_elem] = j;
      a->link[new_elem] = k;
      if (last == -1)
	a->root[i] = new_elem;
      else
	a->link[last] = new_elem;
      
      a->nnz ++;
    }

  } else { /* x == 0.0 */
    
    if (col == j) {
      /* relink row i */
      if (last == -1)
	a->root[i] = a->link[k];
      else
	a->link[last] = a->link[k];
      /* add element to free list */
      a->link[k] = a->free;
      a->free = k;
      
      a->nnz --;
    }
  }
  return 0;
    
 fail:
  PyErr_NoMemory();
  return -1;
}

static int 
SpMatrix_LLMatUpdateItemAdd(LLMatObject *a, int i, int j, double x) {
/*
 *  Update-add matrix entry: a[i,j] += x
 *
* * * * * * * * * * * * * * * * * * * * *
 */
  void *temp;
  int k, new_elem, col, last;

  if (a->issym && i < j) {
    PyErr_SetString(SpMatrix_ErrorObject, 
		    "write operation to upper triangle of symmetric matrix");
    return -1;
  } 

  if (x == 0.0)
    return 0;

  /* find element to be updated */
  col = last = -1;
  k = a->root[i]; 
  while (k != -1) {
    col = a->col[k];
    if (col >= j)
      break;
    last = k;
    k = a->link[k];    
  }
  if (col == j) {
    /* element already exists: compute updated value */
    x += a->val[k];

    if (x == 0.0) {
      /* the updated element is zero and must be removed */
    
      /* relink row i */
      if (last == -1)
	a->root[i] = a->link[k];
      else
	a->link[last] = a->link[k];
      /* add element to free list */
      a->link[k] = a->free;
      a->free = k;
	
      a->nnz --;
    } else {
      a->val[k] = x;
    }

  } else {
    /* new item */
    if (a->free != -1) {

      /* use element from the free chain */
      new_elem = a->free;
      a->free = a->link[new_elem];

    } else {

      /* append new element to the end */
      new_elem = a->nnz;

      /* test if there is space for a new element */
      if (a->nnz == a->nalloc) {
	int nalloc_new;

	/* increase size of idx, val and link arrays */
	nalloc_new = (int)((double)INCREASE_FACTOR * a->nalloc) + 1;
	if ((temp = PyMem_Resize(a->col, int, nalloc_new)) == NULL)
	  goto fail;
	else
	  a->col = temp;
	if ((temp = PyMem_Resize(a->link, int, nalloc_new)) == NULL)
	  goto fail;
	else
	  a->link = temp;
	if ((temp = PyMem_Resize(a->val, double, nalloc_new)) == NULL)
	  goto fail;
	else
	  a->val = temp;
	a->nalloc = nalloc_new;
      }
    }
  
    a->val[new_elem] = x;
    a->col[new_elem] = j;
    a->link[new_elem] = k;
    if (last == -1)
      a->root[i] = new_elem;
    else
      a->link[last] = new_elem;
    a->nnz ++;
  }
  return 0;
  
 fail:
  PyErr_NoMemory();
  return -1;
}

static int
LLMat_Compress(LLMatObject *self, int *nofFreed) {
/*
 *  shrink val, col and link arrays of matrix to minimal size
 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */
  double *val;
  int *col, *link;
  int i, k, k_next, k_last, k_new, nalloc_new;

  /* new size for val, col and link arrays */
  nalloc_new = self->nnz;
  /* remove entries with k >= nalloc_new from free list */
  k_last = -1;
  k = self->free;
  while (k != -1) {
    k_next =  self->link[k];
    if (k >= nalloc_new) {
      if (k_last == -1)
	self->free = k_next;
      else
	self->link[k_last] = k_next;
    } else
      k_last = k;
    k = k_next;
  }
  /* reposition matrix entries with k >= nalloc_new */
  for (i = 0; i < self->dim[0]; i ++) {
    k_last = -1;
    for (k = self->root[i]; k != -1; k = self->link[k]) {
      if (k >= nalloc_new) {
	k_new = self->free;
	if (k_last == -1)
	  self->root[i] = k_new;
	else
	  self->link[k_last] = k_new;
	self->free = self->link[k_new];
	self->val[k_new] = self->val[k];
	self->col[k_new] = self->col[k];
	self->link[k_new] = self->link[k];
	k_last = k_new;
      } else 
	k_last = k;
    }
  }
  /* shrink arrays */
  col = PyMem_Resize(self->col, int, nalloc_new);
  if (col == NULL) goto fail;
  self->col = col;
  link = PyMem_Resize(self->link, int, nalloc_new);
  if (link == NULL) goto fail;
  self->link = link;
  val = PyMem_Resize(self->val, double, nalloc_new);
  if (val == NULL) goto fail;
  self->val = val;

  *nofFreed = self->nalloc - nalloc_new;
  self->nalloc = nalloc_new;
  return 0;

 fail:
  PyErr_NoMemory();
  return -1;
}

/***********************************************************************/
/*  R o u t i n e s   f o r   h a n d l i n g   s u b m a t r i c e s  */
/***********************************************************************/

static int
copySubMatrix(LLMatObject *src, LLMatObject *dst, 
	      int start0, int stop0, 
	      int start1, int stop1, 
	      int transpose) {
/*
 *  Copy submatrix src[start0:stop0,start1:stop1] into a newly allocated matrix.
 *  If transpose is non-zero, the transposed submatrix is copied.
 *
 *  Helper routine for get_submatrix
 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */
  int i, j, k, res;

  for (i = start0; i < stop0; i ++) {
    k = src->root[i];
    while (k != -1) {
      j = src->col[k];
      if (start1 <= j && j < stop1) {
	if (!transpose)
	  res = SpMatrix_LLMatSetItem(dst, i - start0, j - start1, src->val[k]);
	else
	  res = SpMatrix_LLMatSetItem(dst, j - start1, i - start0, src->val[k]);
	if (res) {
	  PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
	  return -1;
	}
      }
      k = src->link[k];
    }
  }
  return 0;
}


static PyObject *
get_submatrix(LLMatObject *self, int start0, int stop0, int start1, int stop1) {
/*
 *  return a new LLMatObject constructed from the slice
 *  self[start0:stop0, start1:stop1].
 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */

  PyObject *op;
  int dim[2];
  int symmetric, transpose, t, res;
  
  dim[0] = stop0 - start0;
  if (dim[0] < 0) dim[0] = 0;
  dim[1] = stop1 - start1;
  if (dim[1] < 0) dim[1] = 0;

  if (self->issym) {

    /* symmetric matrix - 4 cases
     */
    if (start0 == start1 && stop0 == stop1) {
      
      /* case 1: diagonal block -- return symmetric matrix */
      symmetric = 1; transpose = 0;
      
    } else if (start0 >= stop1){
      
      /* case 2: block in lower triangle -- return general matrix */
      symmetric = 0; transpose = 0;

    } else if (stop0 <= start1){
    
      /* case 3: block in upper triangle -- return general matrix */
      t = start0; start0 = start1; start1 = t;
      t = stop0; stop0 = stop1; stop1 = t;
      symmetric = 0; transpose = 1;
      
    } else {

      /* case 4: block overlaps diagonal -- too complicated */
      PyErr_SetString(PyExc_NotImplementedError, 
		      "submatrix overlaps diagonal");
      return NULL;
    }

  } else {

    /* general matrix
     */
    symmetric = 0; transpose = 0;
    
  }

  op = SpMatrix_NewLLMatObject(dim, symmetric, 1000); 
  if (op == NULL)
    return NULL;

  res = copySubMatrix(self, (LLMatObject *)op, start0, stop0, start1, stop1, transpose);
  if (res) {
    Py_DECREF(op);
    return NULL;
  }
  return op;

}

static int
clear_submatrix(LLMatObject *self, int start0, int stop0, int start1, int stop1) {
/*
 *  Delete all non-zero entries from slice self[start0:stop0, start1:stop1]
 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */
  int i, j, k, next, last;

  for (i = start0; i < stop0; i ++) {
    last = -1;
    k = self->root[i];
    while (k != -1) {
      j = self->col[k];
      next = self->link[k];
      if (start1 <= j && j < stop1) {
	/* remove element */
	if (last == -1)
	  self->root[i] = next;
	else
	  self->link[last] = next;
	/* add element to free list */
	self->link[k] = self->free;
	self->free = k;
	self->nnz --;
      } else
	last = k;
      k = next;
    }
  }
  return 0;
}

static int
set_submatrix(LLMatObject *self, int start0, int stop0, int start1, int stop1, LLMatObject *mat) {
/*
 *  assign slice self[start0:stop0, start1:stop1] to matrix mat.
 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */
  double v;
  int dim[2];
  int i, j, k, mirror, res;

  /* validate argument "mat" */
  if (mat->ob_type != &LLMatType) {
    PyErr_SetString(PyExc_ValueError, "Matrix must be of type ""ll_mat"".");
    return -1;
  } 
  dim[0] = stop0 - start0;
  if (dim[0] < 0) dim[0] = 0;
  dim[1] = stop1 - start1;
  if (dim[1] < 0) dim[1] = 0;
  if (mat->dim[0] != dim[0] || mat->dim[1] != dim[1]) {
    PyErr_SetString(PyExc_ValueError, "Matrix shapes are different.");
    return -1;
  }
  
  if (self->issym && !((start0 >= stop1) ||
		       (start0 == start1 && stop0 == stop1 && mat->issym)
		       )) {
    PyErr_SetString(SpMatrix_ErrorObject, 
		    "write operation to upper triangle of symmetric matrix");
    return -1;
  }

  /* clear elements from slice [start0:stop0, start1:stop1] */
  clear_submatrix(self, start0, stop0, start1, stop1);
  
  /* if mat is symmetric and mat is not assigned to a diagonal block
   * of symmetric matrix, then the non-diagonal elements of mat will
   * have to be mirrored
   */
  mirror = mat->issym && !(start0 == start1 && stop0 == stop1 && self->issym);
    
  /* copy mat, element by element */
  for (i = 0; i < mat->dim[0]; i ++) {
    k = mat->root[i];
    while (k != -1) {
      j = mat->col[k];
      v = mat->val[k];
      res = SpMatrix_LLMatSetItem(self, start0 + i, start1 + j, v);
      if (mirror && i != j)
	res = res || SpMatrix_LLMatSetItem(self, start0 + j, start1 +i, v);
      if (res) {
	PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
	return -1;
      }
      k = mat->link[k];
    }
  }
  
  return 0;
}

/*****************************************************************************/
/*  M a t r i x - v e c t o r   m u l t i p l i c a t i o n   k e r n e l s  */
/*****************************************************************************/

static void
ll_matvec_kernel(int m, double *x, double *y,
		 double *val, int *col, int *link, int *root) {
  double s;
  int i, k;
  
  for (i = 0; i < m; i ++) {
    s = 0.0;
    k = root[i];
    while (k != -1) {
      s += val[k] * x[col[k]];
      k = link[k];
    }
    y[i] = s;
  }
}

static void
ll_matvec_kernel_stride(int m, 
			double *x, int incx, 
			double *y, int incy,
			double *val, int *col, int *link, int *root) {
  double s;
  int i, k;
  
  for (i = 0; i < m; i ++) {
    s = 0.0;
    k = root[i];
    while (k != -1) {
      s += val[k] * x[col[k]*incx];
      k = link[k];
    }
    y[i*incy] = s;
  }
}

static void
ll_matvec_kernel_sym(int m, double *x, double *y,
		     double *val, int *col, int *link, int *root) {
  double s, v, xi;
  int i, j, k;
  
  for (i = 0; i < m; i ++) {
    xi = x[i];
    s = 0.0;
    k = root[i];
    while (k != -1) {
      j = col[k];
      v = val[k];
      s += v * x[j];
      if (i != j)
	y[j] += v * xi;
      k = link[k];
    }
    y[i] = s;
  }
}

static void
ll_matvec_kernel_stride_sym(int m, 
			    double *x, int incx, 
			    double *y, int incy,
			    double *val, int *col, int *link, int *root) {
  double s, v, xi;
  int i, j, k;
  
  for (i = 0; i < m; i ++) {
    xi = x[i*incx];
    s = 0.0;
    k = root[i];
    while (k != -1) {
      j = col[k];
      v = val[k];
      s += v * x[j*incx];
      if (i != j)
	y[j*incy] += v * xi;
      k = link[k];
    }
    y[i*incy] = s;
  }
}

static void
ll_matvec_transp_kernel(int m, int n, double *x, double *y,
			    double *val, int *col, int *link, int *root) {
  double xi;
  int i, k;
  
  for (i = 0; i < n; i ++)
    y[i] = 0.0;
  
  for (i = 0; i < m; i ++) {
    xi = x[i];
    k = root[i];
    while (k != -1) {
      y[col[k]] += val[k] * xi;
      k = link[k];
    }
  }
}

static void
ll_matvec_transp_kernel_stride(int m, int n, 
			       double *x, int incx, 
			       double *y, int incy,
			       double *val, int *col, int *link, int *root) {
  double xi;
  int i, k;
  
  for (i = 0; i < n; i ++)
    y[i*incy] = 0.0;
  
  for (i = 0; i < m; i ++) {
    xi = x[i*incx];
    k = root[i];
    while (k != -1) {
      y[col[k]*incy] += val[k] * xi;
      k = link[k];
    }
  }
}

/*********************************************/
/*  L L M a t   o b j e c t   m e t h o d s  */
/*********************************************/

static char LLMat_matvec_transp_doc[] = "a.matvec_transp(x, y)\n\
\n\
compute the sparse matrix-vector product y := a^T * x. \n\
a^T is the transpose of a, which is a d1 by d2 sparse matrix.\n\
x and y are two 1-dimensional Numeric arrays of appropriate size.";

static PyObject *
LLMat_matvec_transp(LLMatObject *self, PyObject *args)
{
  PyArrayObject *xp, *yp;

  SPMATRIX_PARSE_ARGS_ARR_ARR_STRIDE(args, xp, yp, self->dim[0], self->dim[1]);
  
  if (xp->flags & CONTIGUOUS &&  yp->flags & CONTIGUOUS)
    if (self->issym)
      ll_matvec_kernel_sym(self->dim[0], (double *)(xp->data), (double *)(yp->data), 
			   self->val, self->col, self->link, self->root);
    else
      ll_matvec_transp_kernel(self->dim[0], self->dim[1], (double *)(xp->data), (double *)(yp->data), 
			      self->val, self->col, self->link, self->root);
  else {
    if (self->issym)
      ll_matvec_kernel_stride_sym(self->dim[0], 
				  (double *)(xp->data), xp->strides[0] / sizeof(double),
				  (double *)(yp->data), yp->strides[0] / sizeof(double),
				  self->val, self->col, self->link, self->root);
    else
      ll_matvec_transp_kernel_stride(self->dim[0], self->dim[1], 
				     (double *)(xp->data), xp->strides[0] / sizeof(double),
				     (double *)(yp->data), yp->strides[0] / sizeof(double),
				     self->val, self->col, self->link, self->root);
  }

  Py_INCREF(Py_None); 
  return Py_None;
}

static char LLMat_matvec_doc[] = "a.matvec(x, y)\n\
\n\
compute the sparse matrix-vector product y := a * x. \n\
a is a d1 by d2 sparse matrix.\n\
x and y are two 1-dimensional Numeric arrays of appropriate size.";

static PyObject *
LLMat_matvec(LLMatObject *self, PyObject *args)
{
  PyArrayObject *xp, *yp;

  SPMATRIX_PARSE_ARGS_ARR_ARR_STRIDE(args, xp, yp, self->dim[1], self->dim[0]);
     
  if (xp->flags & CONTIGUOUS &&  yp->flags & CONTIGUOUS)
    if (self->issym)
      ll_matvec_kernel_sym(self->dim[0], (double *)(xp->data), (double *)(yp->data), 
			   self->val, self->col, self->link, self->root);
    else
      ll_matvec_kernel(self->dim[0], (double *)(xp->data), (double *)(yp->data), 
		       self->val, self->col, self->link, self->root);
  else {
    if (self->issym)
      ll_matvec_kernel_stride_sym(self->dim[0], 
				  (double *)(xp->data), xp->strides[0] / sizeof(double),
				  (double *)(yp->data), yp->strides[0] / sizeof(double),
				  self->val, self->col, self->link, self->root);
    else
      ll_matvec_kernel_stride(self->dim[0], 
			      (double *)(xp->data), xp->strides[0] / sizeof(double),
			      (double *)(yp->data), yp->strides[0] / sizeof(double),
			      self->val, self->col, self->link, self->root);
  }

  Py_INCREF(Py_None); 
  return Py_None;
}

static char to_csr_doc[] = "a.to_csr()\n\
\n\
return a new CSRMatObject constructed from data of a";

static PyObject *
LLMat_to_csr(LLMatObject *self, PyObject *args)
{
  CSRMatObject *op;
  int i, j, k, r;
  
  if (!PyArg_ParseTuple(args, ""))
    return NULL;

  if (self->issym) {

    /* symmetric case
     */
    struct llColIndex *colIdx;

    if (SpMatrix_LLMatBuildColIndex(&colIdx, self, 0))
      return NULL;
    assert(colIdx->nzUp == 0);

    op = (CSRMatObject *)newCSRMatObject(self->dim, 2*colIdx->nzLo + colIdx->nzDiag);
    if (op == NULL) {
      SpMatrix_LLMatDestroyColIndex(&colIdx);
      return NULL;
    }

    r = 0;
    op->ind[0] = 0;
    for (i = 0; i < self->dim[0]; i ++) {

      /* store self[0:i+1,i] in op[0:i+1,i] */
      k = self->root[i];
      while (k != -1) {
	op->val[r] = self->val[k];
	op->col[r] = self->col[k];
	r ++;
	k = self->link[k];
      }

      /* store self[i,i+1:n] in op[i+1:n,i] */
      k = colIdx->root[i];
      while (k != -1) {
	j = colIdx->row[k];
	op->val[r] = self->val[k];
	op->col[r] = j;
	r ++;
	k = colIdx->link[k];
      }
      
      op->ind[i+1] = r;
    }
  
    SpMatrix_LLMatDestroyColIndex(&colIdx);

  } else {
    
    /* unsymmetric case
     */

    op = (CSRMatObject *)newCSRMatObject(self->dim, self->nnz);
    if (op == NULL)
      return NULL;

    r = 0;
    op->ind[0] = 0;
    for (i = 0; i < self->dim[0]; i ++) {
      k = self->root[i];
      while (k != -1) {
	op->val[r] = self->val[k];
	op->col[r] = self->col[k];
	r ++;
	k = self->link[k];
      }
      op->ind[i+1] = r;
    }
  }

  return (PyObject *)op;
}

static char to_sss_doc[] = "a.to_sss()\n\
\n\
return a new SSSMatObject constructed from the lower triangle of a";

static PyObject *
LLMat_to_sss(LLMatObject *self, PyObject *args)
{
  SSSMatObject *op;
  int i, j, k, r, n, nnz;
  
  if (!PyArg_ParseTuple(args, ""))
    return NULL;

  /* test for square matrix */
  n = self->dim[0];
  if (n != self->dim[1]) {
    PyErr_SetString(PyExc_ValueError, "Matrix must be square");
    return NULL;
  }
  
  /* 1st pass: compute number of non-zeros in lower triangle */
  nnz = 0;
  for (i = 0; i < n; i ++) {
    k = self->root[i];
    while (k != -1) {
      j = self->col[k];
      if (i > j)
	nnz ++;
      k = self->link[k];
    }
  }
  
  /* allocate new SSS matrix */
  op = (SSSMatObject *)newSSSMatObject(n, nnz);
  if (op == NULL)
    return NULL;

  /* 2nd pass: fill SSSMatObject */
  for (i = 0; i < n; i ++)
    op->diag[i] = 0.0;
  r = 0;
  op->ind[0] = 0;
  for (i = 0; i < n; i ++) {
    k = self->root[i];
    while (k != -1) {
      j = self->col[k];
      if (i > j) {
	op->val[r] = self->val[k];
	op->col[r] = j;
	r ++;
      } else if (i == j)
	op->diag[i] = self->val[k];
      k = self->link[k];
    }
    op->ind[i+1] = r;
  }

  return (PyObject *)op;
}

static char LLMat_generalize_doc[] = "convert ll_mat object from symmetric to non-symmetric form (in-place).";

static PyObject *
LLMat_generalize(LLMatObject *self, PyObject *args) {
  int i, j, k;

  if (!PyArg_ParseTuple(args, ""))
    return NULL;

  if (self->issym) {
    self->issym = 0;
    for (i = 0; i < self->dim[0]; i ++) {
      /* New elements are inserted into the matrix while it is being traversed.
	 However, this should not be a problem */
      for (k = self->root[i]; k != -1; k = self->link[k]) {
	j = self->col[k];
	if (i > j)
	  if (SpMatrix_LLMatSetItem(self, j, i, self->val[k]))
	    return NULL;
      }
    }
  }

  Py_INCREF(Py_None); 
  return Py_None;
}

static char LLMat_compress_doc[] = "a.compress() frees memory by reclaiming unused space in the data structures of a. returns number of elements freed";

static PyObject *
LLMat_compress(LLMatObject *self, PyObject *args) {
  int nofFreed;

  if (!PyArg_ParseTuple(args, ""))
    return NULL;
  
  if (LLMat_Compress(self, &nofFreed))
    return NULL;

  return PyInt_FromLong(nofFreed);
}

static char export_mtx_doc[] = "a.export_mtx(fileName, precision=16)\n\
\n\
write matrix in Matrix Market format to file f.\n\
\n\
parameters:\n\
\n\
fileName:  string, name of the file to be created\n\
precision: number of significant digits to be written for double values";

static PyObject *
LLMat_export_mtx(LLMatObject *self, PyObject *args)
{
  char *fileName;
  int precision = 16;
  MM_typecode matcode;
  FILE *f;
  int ret, i, k;

  if (!PyArg_ParseTuple(args, "s|i", &fileName, &precision))
    return NULL;

  f = fopen(fileName, "w");
  if (f == NULL)
    return PyErr_SetFromErrno(PyExc_IOError);
  
  mm_set_matrix(matcode); mm_set_sparse(matcode); mm_set_real(matcode); 
  if (self->issym)
    mm_set_symmetric(matcode);
  else
    mm_set_general(matcode);
    
  ret = mm_write_banner(f, matcode);
  if (ret) {
    PyErr_SetString(SpMatrix_ErrorObject, "error writing file header");    
    return NULL;
  }
  
  ret = fprintf(f, "%% file created by pysparse module\n");
  if (ret < 0) {
    PyErr_SetString(PyExc_IOError, "error writing file header");    
    return NULL;
  }

  ret = mm_write_mtx_crd_size(f, self->dim[0], self->dim[1], self->nnz);
  if (ret) {
    PyErr_SetString(SpMatrix_ErrorObject, "error writing file header");    
    return NULL;
  }
  
  for (i = 0; i < self->dim[0]; i ++) {
    k = self->root[i];
    while (k != -1) {
      ret = fprintf(f, "%d %d %.*e\n", i+1, self->col[k]+1, precision-1, self->val[k]);
      if (ret < 0) {
	PyErr_SetString(PyExc_IOError, "error writing matrix data");    
	return NULL;
      }
      k = self->link[k];
    }
  }
  
  ret = fclose(f);
  if (ret)
    return PyErr_SetFromErrno(PyExc_IOError);

  Py_INCREF(Py_None); 
  return Py_None;  
}

static char copy_doc[] = "a.copy()\n\
\n\
return a (deep) copy of the matrix a.";

static PyObject *
LLMat_copy(LLMatObject *self, PyObject *args)
{
  LLMatObject *new;
  int i, k;

  if (!PyArg_ParseTuple(args, ""))
    return NULL;
  
  new = (LLMatObject *)SpMatrix_NewLLMatObject(self->dim, self->issym, self->nnz);
  if (new == NULL)
    return NULL;

  for (i = 0; i < self->dim[0]; i ++) {
    k = self->root[i];
    while (k != -1) {
      if (SpMatrix_LLMatSetItem(new, i, self->col[k], self->val[k]) == -1) {
	Py_DECREF(new);
	return NULL;
      }
      k = self->link[k];
    }
  }
  return (PyObject *)new;
}

static char update_add_at_doc[] = "a.update_add_at(b,id1,id2)\n\
\n\
for i in range(len(b)):\n\
    a[id1[i],id2[i]] += b[i]";

static PyObject *
LLMat_update_add_at(LLMatObject *self, PyObject *args) {
  PyObject *bIn;
  PyObject *id1in;
  PyObject *id2in;
  PyArrayObject *b = NULL;
  PyArrayObject *id1 = NULL;
  PyArrayObject *id2 = NULL;
  double v;
  int lenb,i;

  if (self->issym) 
    {
      PyErr_SetString(SpMatrix_ErrorObject, "method not allowed for symmetric matrices");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "OOO", &bIn,&id1in,&id2in))
    return NULL;

  b = (PyArrayObject *)PyArray_ContiguousFromObject(bIn, PyArray_DOUBLE, 1, 1);
  if (b == NULL)
    goto fail;

  lenb = b->dimensions[0];

  id1 = (PyArrayObject *)PyArray_ContiguousFromObject(id1in, PyArray_LONG, 1, 1);
  if (id1 == NULL)
    goto fail;
  
  id2 = (PyArrayObject *)PyArray_ContiguousFromObject(id2in, PyArray_LONG, 1, 1);
  if (id2 == NULL)
    goto fail;
  
  if(self->dim[0]!=self->dim[1]){
    PyErr_SetString(PyExc_IndexError, "dim[0] and dim[1] are different sizes");
    goto fail;}

  if (lenb < 0 ) {
    PyErr_SetString(PyExc_IndexError, "vector b is a negative size");
    goto fail;}

  if (id1->dimensions[0] != lenb ) {
    PyErr_SetString(PyExc_IndexError, "id1 is not the same size as b");
    goto fail;}

  if (id2->dimensions[0] != lenb ) {
    PyErr_SetString(PyExc_IndexError, "id2 is not the same size as b");
    goto fail;}

  /* perform update add operation */

  for (i = 0; i < lenb; i ++)
    {
      v = ((double *)b->data)[i];
      if (SpMatrix_LLMatUpdateItemAdd(self, ((long *) id1->data)[i], ((long *) id2->data)[i], v) == -1)
	goto fail;
    }

  Py_DECREF(b);
  Py_DECREF(id1);
  Py_DECREF(id2);
  Py_INCREF(Py_None); 
  return Py_None;

 fail:
  if (b) {
      Py_XDECREF(b);
  }  
  if (id1) {
      Py_XDECREF(id1);
  }
  if (id2) {
      Py_XDECREF(id2);
  }
  
  return NULL;
}

static char LLMat_norm_doc[] = "return p-norm of matrix\n\
\n\
A.norm(p) returns the p-norm of matrix A\n\
\n\
p ist a string identifying the type of norm to be returned\n\
\n\
  '1'   -- return the 1-norm of A\n\
  'inf' -- return the infinity norm of A\n\
  'fro' -- return the frobenius norm of A";

static PyObject *
LLMat_norm(LLMatObject *self, PyObject *args)
{
  char *normType;

  struct llColIndex *colIdx;
  double norm, s, sMax, v;
  int i, k;
  
  if (!PyArg_ParseTuple(args, "s", &normType))
    return NULL;
  
  if (strcmp(normType, "1") == 0) {
    if (self->issym) {
      PyErr_SetString(PyExc_NotImplementedError, "not implemented for symmetric matrices");
      return NULL;
    } else {
      if (SpMatrix_LLMatBuildColIndex(&colIdx, self, 1))
	return NULL;
      sMax = 0.0;
      for (i = 0; i < self->dim[1]; i ++) {
	s = 0.0;
	for (k = colIdx->root[i]; k != -1; k = colIdx->link[k]) {
	  s += fabs(self->val[k]);
	}
	if (s > sMax)
	  sMax = s;
      }
      norm = sMax;
      SpMatrix_LLMatDestroyColIndex(&colIdx);
    }
  } else if (strcmp(normType, "inf") == 0) {
    if (self->issym) {
      PyErr_SetString(PyExc_NotImplementedError, "not implemented for symmetric matrices");
      return NULL;
    } else {
      sMax = 0.0;
      for (i = 0; i < self->dim[0]; i ++) {
	s = 0.0;
	for (k = self->root[i]; k != -1; k = self->link[k]) {
	  s += fabs(self->val[k]);
	}
	if (s > sMax)
	  sMax = s;
      }
      norm = sMax;
    }
  } else if (strcmp(normType, "fro") == 0) {
    s = 0.0;
    for (i = 0; i < self->dim[0]; i ++) {
      for (k = self->root[i]; k != -1; k = self->link[k]) {
	v = self->val[k];
	s += v*v;
	if (self->issym && self->col[k] != i)
	  s += v*v;
      }
    }
    norm = sqrt(s);
  } else {
    PyErr_SetString(PyExc_ValueError, "unknown norm type");
    return NULL;
  }

  return Py_BuildValue("d", norm);
}

static char shift_doc[] = "a.shift(sigma, b)\n\
\n\
shift the matrix:\n\
compute a = a + sigma * b";

static PyObject *
LLMat_shift(LLMatObject *self, PyObject *args)
{
  LLMatObject *mat;
  double sigma, v;
  int i, j, k;
  
  if (!PyArg_ParseTuple(args, "dO!", &sigma, &LLMatType, &mat))
    return NULL;
  if (self->dim[0] != mat->dim[0] || self->dim[1] != mat->dim[1]) {
    PyErr_SetString(PyExc_ValueError, "matrix shapes do not match");
    return NULL;
  }

  if (self->issym == mat->issym) {
    for (i = 0; i < mat->dim[0]; i ++) {
      k = mat->root[i];
      while (k != -1) {
	if (SpMatrix_LLMatUpdateItemAdd(self, i, mat->col[k], sigma * mat->val[k]) == -1)
	  return NULL;
	k = mat->link[k];
      }
    }
  } else if (mat->issym) {
    for (i = 0; i < mat->dim[0]; i ++) {
      k = mat->root[i];
      while (k != -1) {
	j = mat->col[k];
	v = sigma * mat->val[k];
	if (SpMatrix_LLMatUpdateItemAdd(self, i, j, v) == -1)
	  return NULL;
	if (i != j)
	  if (SpMatrix_LLMatUpdateItemAdd(self, j, i, v) == -1)
	    return NULL;
	
	k = mat->link[k];
      }
    }
  } else {
    PyErr_SetString(PyExc_NotImplementedError, 
		    "shift of symmetric matrix by non-symmetric matrix not supported");
    return NULL;
    
  }

  Py_INCREF(Py_None); 
  return Py_None;
}


static char keys_doc[] = "A.keys()\n\
\n\
Return a list of tuples (i,j) of non-zero matrix entries.";

static PyObject *
LLMat_keys(LLMatObject *a, PyObject *args)
{
    PyObject *list;             /* the list that will hold the keys */
    int i, j, k;
    int pos = 0;                /* position in list */
    
    if (!PyArg_ParseTuple(args, ""))
        return NULL;
    
    if (!a->issym) {

        if ((list = PyList_New(a->nnz)) == NULL)
            return NULL;
        
        for (i = 0; i < a->dim[0]; i ++) {
	    k = a->root[i];
	    while (k != -1) {
                j = a->col[k];
                PyList_SET_ITEM(list, pos++, Py_BuildValue("ii", i, j));
                k = a->link[k];
	    }
        }
        return list;
        
    } else {
        PyErr_SetString(PyExc_NotImplementedError, 
                        "keys() doesn't yet support symmetric matrices");
        return NULL;
    }
}

static char values_doc[] = "A.values()\n\
\n\
Return a list of the non-zero matrix entries as floats.";

static PyObject *
LLMat_values(LLMatObject *a, PyObject *args)
{
    PyObject *list;           /* the list that will hold the values */
    int i, k;
    int pos = 0;                /* position in list */

    if (!PyArg_ParseTuple(args, ""))
        return NULL;
    
    if (!a->issym) {
        
        if ((list = PyList_New(a->nnz)) == NULL)
            return NULL;
        
        for (i = 0; i < a->dim[0]; i ++) {
	    k = a->root[i];
	    while (k != -1) {
                PyList_SET_ITEM(list, pos++, PyFloat_FromDouble(a->val[k]));
                k = a->link[k];
	    }
        }
        return list;

    } else {
        PyErr_SetString(PyExc_NotImplementedError, 
                        "values() doesn't yet support symmetric matrices");
        return NULL;
    }
}



static char items_doc[] = "A.items()\n\
\n\
Return a list of tuples (indices, value) of\n\
the non-zero matrix entries' keys and values.\n\
\n\
The indices are themselves tuples (i,j) of row\n\
and column values.";

static PyObject *
LLMat_items(LLMatObject *a, PyObject *args)
{
    PyObject *list;           /* the list that will hold the values */
    int i, j, k;
    int pos = 0;                /* position in list */
    double val;
    
    if (!PyArg_ParseTuple(args, ""))
        return NULL;
    
    if (!a->issym) {
        
        if ((list = PyList_New(a->nnz)) == NULL)
            return NULL;
        
        for (i = 0; i < a->dim[0]; i ++) {
            k = a->root[i];
            while (k != -1) {
                j = a->col[k];
                val = a->val[k];
                PyList_SET_ITEM(list, pos++, Py_BuildValue("((ii)d)", i, j, val));
                k = a->link[k];
            }
        }
        return list;
        
    } else {
        PyErr_SetString(PyExc_NotImplementedError, 
                        "items() doesn't yet support symmetric matrices");
        return NULL;
    }
}

static char scale_doc[] = "a.scale(sigma)\n\
\n\
Scale each element in the matrix by the constant sigma.\n";

static PyObject *
LLMat_scale(LLMatObject *self, PyObject *args)
{
  double sigma;
  int i, k;
  
  if (!PyArg_ParseTuple(args, "d", &sigma))
      return NULL;
  
  for (i = 0; i < self->dim[0]; i ++) {
      k = self->root[i];
      while (k != -1) {
          self->val[k] *= sigma;
          k = self->link[k];
      }
  }
  
  Py_INCREF(Py_None); 
  return Py_None;
}


static char update_add_mask_doc[] = "a.update_add_mask(b, ind0, ind1, mask0, mask1)\n\
\n\
Update of global FEM matrix. Equivalent to:\n\
\n\
for i in range(len(ind0)):\n\
    for j in range(len(ind1)):\n\
        if mask0[i] and mask1[j]:\n\
            a[ind0[i],ind1[j]] += b[i,j]";

static PyObject *
LLMat_update_add_mask(LLMatObject *self, PyObject *args) {
  PyObject *bIn, *ind0In, *ind1In, *mask0In, *mask1In; 
  PyArrayObject *b, *ind0, *ind1, *mask0, *mask1;
  double v;
  int len0, len1, i, j, i1, j1, ldb;

  if (self->issym) {
    PyErr_SetString(SpMatrix_ErrorObject, "method not allowed for symmetric matrices");
    return NULL;
  }

  if (!PyArg_ParseTuple(args, "OOOOO", &bIn, &ind0In, &ind1In, &mask0In, &mask1In)) 
    return NULL;

  b = (PyArrayObject *)PyArray_ContiguousFromObject(bIn, PyArray_DOUBLE, 2, 2);
  ind0 = (PyArrayObject *)PyArray_ContiguousFromObject(ind0In, PyArray_LONG, 1, 1);
  ind1 = (PyArrayObject *)PyArray_ContiguousFromObject(ind1In, PyArray_LONG, 1, 1);
  mask0 = (PyArrayObject *)PyArray_ContiguousFromObject(mask0In, PyArray_LONG, 1, 1);
  mask1 = (PyArrayObject *)PyArray_ContiguousFromObject(mask1In, PyArray_LONG, 1, 1);

  if (b == NULL || ind0 == NULL || ind1 == NULL || mask0 == NULL || mask1 == NULL)
    goto fail;

  len0 = ind0->dimensions[0];
  len1 = ind1->dimensions[0];

  /* validate array shapes */
  if (mask0->dimensions[0] != len0 || mask1->dimensions[0] != len1) {
    PyErr_SetString(PyExc_ValueError, "shapes of index and mask array do not match");
    goto fail;
  }
  if (b->dimensions[0] != len0 || b->dimensions[1] != len1) {
    PyErr_SetString(PyExc_ValueError, "shapes of input matrix and index arrays do not match");
    goto fail;
  }
  
  /* perform update add operation */
  ldb = b->dimensions[0];
  for (i = 0; i < len0; i ++) {
    if (((long *)mask0->data)[i]) {

      i1 = ((long *)ind0->data)[i];
      if (i1 < 0)
	i1 += self->dim[0];
      if (i1 < 0 || i1 >= self->dim[0]) {
	PyErr_SetString(PyExc_IndexError, "element of arg 2 out of range");
	goto fail;
      }
   
      for (j = 0; j < len1; j ++) {
	if (((long *)mask1->data)[j]) {

	  j1 = ((long *)ind1->data)[j];
	  if (j1 < 0)
	    j1 += self->dim[1];
	  if (j1 < 0 || j1 >= self->dim[1]) {
	    PyErr_SetString(PyExc_IndexError, "element of arg 3 out of range");
	    goto fail;
	  }

	  v = ((double *)b->data)[i + ldb*j];
	  if (SpMatrix_LLMatUpdateItemAdd(self, i1, j1, v) == -1)
	    goto fail;
	}
      }
    }
  }
  
  Py_DECREF(b);
  Py_DECREF(ind0);
  Py_DECREF(ind1);
  Py_DECREF(mask0);
  Py_DECREF(mask1);
  Py_INCREF(Py_None); 
  return Py_None;

 fail:
  Py_XDECREF(b);
  Py_XDECREF(ind0);
  Py_XDECREF(ind1);
  Py_XDECREF(mask0);
  Py_XDECREF(mask1);
  return NULL;
}

static char update_add_mask_sym_doc[] = "a.update_add_mask(b, ind, mask)\n\
\n\
Symmetric update of global FEM matrix. Equivalent to:\n\
\n\
for i in range(len(ind)):\n\
    for j in range(len(ind)):\n\
        if mask[i] and mask[i]:\n\
            i1 = ind[i]; j1 = ind[j]\n\
            if i >= j:\n\
                a[i1,j1] += b[i,j]\n\
\n\
Only operates on the lower triangle of a. Used for symmetric matrices.";

static PyObject *
LLMat_update_add_mask_sym(LLMatObject *self, PyObject *args) {
  PyObject *bIn, *indIn, *maskIn; 
  PyArrayObject *b, *ind, *mask;
  double v;
  int len, i, j, i1, j1, ldb;

  if (!PyArg_ParseTuple(args, "OOO", &bIn, &indIn, &maskIn)) 
    return NULL;

  b = (PyArrayObject *)PyArray_ContiguousFromObject(bIn, PyArray_DOUBLE, 2, 2);
  ind = (PyArrayObject *)PyArray_ContiguousFromObject(indIn, PyArray_LONG, 1, 1);
  mask = (PyArrayObject *)PyArray_ContiguousFromObject(maskIn, PyArray_LONG, 1, 1);

  if (b == NULL || ind == NULL || mask == NULL)
    goto fail;

  len = ind->dimensions[0];

  /* validate array shapes */
  if (mask->dimensions[0] != len) {
    PyErr_SetString(PyExc_ValueError, "shapes of index and mask array do not match");
    goto fail;
  }
  if (b->dimensions[0] != len || b->dimensions[1] != len) {
    PyErr_SetString(PyExc_ValueError, "shapes of input matrix and index arrays do not match");
    goto fail;
  }
  
  /* perform update add operation */
  ldb = b->dimensions[0];
  for (i = 0; i < len; i ++) {
    if (((long *)mask->data)[i]) {

      i1 = ((long *)ind->data)[i];
      if (i1 < 0)
	i1 += self->dim[0];
      if (i1 < 0 || i1 >= self->dim[0]) {
	PyErr_SetString(PyExc_IndexError, "element of arg 2 out of range");
	goto fail;
      }
   
      for (j = 0; j <= i; j ++) {
	if (((long *)mask->data)[j]) {

	  j1 = ((long *)ind->data)[j]; /* index check not necessary here */
	  if (j1 < 0)
	    j1 += self->dim[1];
	  v = ((double *)b->data)[i + ldb*j];

	  if (self->issym) {
	    /* symmetric matrix: update entry in lower triangle */
	    if (i1 > j1) {
	      if (SpMatrix_LLMatUpdateItemAdd(self, i1, j1, v) == -1)
		goto fail;
	    } else {
	      if (SpMatrix_LLMatUpdateItemAdd(self, j1, i1, v) == -1)
		goto fail;
	    }
	  } else {
	    /* non-symmetric matrix: update two entries if not on diagonal */
	    if (SpMatrix_LLMatUpdateItemAdd(self, i1, j1, v) == -1)
	      goto fail;
	    if (i1 != j1) {
	      if (SpMatrix_LLMatUpdateItemAdd(self, j1, i1, v) == -1)
		goto fail;
	    }
	  }

	}
      }
    }
  }
  
  Py_DECREF(b);
  Py_DECREF(ind);
  Py_DECREF(mask);
  Py_INCREF(Py_None); 
  return Py_None;

 fail:
  Py_XDECREF(b);
  Py_XDECREF(ind);
  Py_XDECREF(mask);
  return NULL;
}

static char LLMat_take_doc[] = "a.take(b,id1,id2)\n\
\n\
for i in range(len(b)):\n\
    b[i] = a[id1[i],id2[i]]";

static PyObject *
LLMat_take(LLMatObject *self, PyObject *args) {
  PyObject *bIn;
  PyObject *id1in = NULL;
  PyObject *id2in = NULL;
  PyArrayObject *b;
  PyArrayObject *id1 = NULL;
  PyArrayObject *id2 = NULL;
  int lenb,i;

  if (!PyArg_ParseTuple(args, "O|OO", &bIn,&id1in,&id2in))
    return NULL;

  b = (PyArrayObject *)PyArray_ContiguousFromObject(bIn, PyArray_DOUBLE, 1, 1);
  if (b == NULL)
    goto fail;

  lenb = b->dimensions[0];

  if (id1in) {
      id1 = (PyArrayObject *)PyArray_ContiguousFromObject(id1in, PyArray_LONG, 1, 1);
      if (id1 == NULL)
	goto fail;
  }
  
  if (id2in) {
      id2 = (PyArrayObject *)PyArray_ContiguousFromObject(id2in, PyArray_LONG, 1, 1);
      if (id2 == NULL)
	goto fail;
  }
  
  if (self->dim[0] != self->dim[1]){
      PyErr_SetString(PyExc_IndexError, "dim[0] and dim[1] are different sizes");
      goto fail;
  }
  
  if (lenb < 0 ) {
      PyErr_SetString(PyExc_IndexError, "vector b is a negative size");
      goto fail;
  }
  
  if (id1 && id1->dimensions[0] != lenb ) {
      PyErr_SetString(PyExc_IndexError, "id1 is not the same size as b");
      goto fail;
  }

  if (id2 && id2->dimensions[0] != lenb ) {
      PyErr_SetString(PyExc_IndexError, "id2 is not the same size as b");
      goto fail;
  }
  
  if (id1 != id2 && self->issym) {
      PyErr_SetString(SpMatrix_ErrorObject, "symmetric matrices require identical sets of indices");
      goto fail;
  }
    
  /* perform take operation */

    for (i = 0; i < lenb; i ++)
      {
	  int	i1, j1;
	  
	  if (id1) {
	      i1 = ((long *) id1->data)[i];
	  } else {
	      i1 = i;
	  }
	  if (id2) {
	      j1 = ((long *) id2->data)[i];
	  } else {
	      j1 = i1;
	  }

	  if (i1 > j1 || !self->issym) {
	      /* get entries as given */
	      ((double *)b->data)[i] = SpMatrix_LLMatGetItem(self, i1, j1);
	  } else {
	      /* symmetric matrix: get entries from lower triangle */
	      ((double *)b->data)[i] = SpMatrix_LLMatGetItem(self, j1, i1);
	  }
      }
      
    Py_DECREF(b);
    if (id1) {
	Py_DECREF(id1);
    }
    if (id2) {
	Py_DECREF(id2);
    }
    Py_INCREF(Py_None); 
    return Py_None;

   fail:
    Py_XDECREF(b);
    if (id1) {
	Py_XDECREF(id1);
    }
    if (id2) {
	Py_XDECREF(id2);
    }
    return NULL;
}

static char LLMat_put_doc[] = "a.put(b,id1,id2)\n\
\n\
for i in range(len(b)):\n\
    a[id1[i],id2[i]] = b[i]";

static PyObject *
LLMat_put(LLMatObject *self, PyObject *args) {
    PyObject *bIn;
    PyObject *id1in = NULL;
    PyObject *id2in = NULL;
    PyArrayObject *b;
    PyArrayObject *id1 = NULL;
    PyArrayObject *id2 = NULL;
    int lenb,i;

    if (!PyArg_ParseTuple(args, "O|OO", &bIn,&id1in,&id2in))
      return NULL;

    b = (PyArrayObject *)PyArray_ContiguousFromObject(bIn, PyArray_DOUBLE, 1, 1);
    if (b == NULL)
      goto fail;

    lenb = b->dimensions[0];

    if (id1in) {
	id1 = (PyArrayObject *)PyArray_ContiguousFromObject(id1in, PyArray_LONG, 1, 1);
	if (id1 == NULL)
	  goto fail;
    }
    
    if (id2in) {
	id2 = (PyArrayObject *)PyArray_ContiguousFromObject(id2in, PyArray_LONG, 1, 1);
	if (id2 == NULL)
	  goto fail;
    }
    
    if (self->dim[0] != self->dim[1]){
	PyErr_SetString(PyExc_IndexError, "dim[0] and dim[1] are different sizes");
	goto fail;
    }
    
    if (lenb < 0 ) {
	PyErr_SetString(PyExc_IndexError, "vector b is a negative size");
	goto fail;
    }
    
    if (id1 && id1->dimensions[0] != lenb ) {
	PyErr_SetString(PyExc_IndexError, "id1 is not the same size as b");
	goto fail;
    }

    if (id2 && id2->dimensions[0] != lenb ) {
	PyErr_SetString(PyExc_IndexError, "id2 is not the same size as b");
	goto fail;
    }
    
  /* perform put operation */

   for (i = 0; i < lenb; i ++)
     {
	 int	i1, j1;
	 
	 if (id1) {
	     i1 = ((long *) id1->data)[i];
	 } else {
	     i1 = i;
	 }
	 if (id2) {
	     j1 = ((long *) id2->data)[i];
	 } else {
	     j1 = i1;
	 }

	 if (i1 > j1 || !self->issym) {
	     /* update entries as given */
	     if (SpMatrix_LLMatSetItem(self, i1, j1, ((double *)b->data)[i]) == -1) {
		 goto fail;
	     }
	 } else {
	     /* symmetric matrix: update entries in lower triangle */
	     if (SpMatrix_LLMatSetItem(self, j1, i1, ((double *)b->data)[i]) == -1) {
		 goto fail;
	     }
	 }
     }
     
    Py_DECREF(b);
    if (id1) {
	Py_DECREF(id1);
    }
    if (id2) {
	Py_DECREF(id2);
    }
    Py_INCREF(Py_None); 
    return Py_None;

   fail:
    Py_XDECREF(b);
    if (id1) {
	Py_XDECREF(id1);
    }
    if (id2) {
	Py_XDECREF(id2);
    }
    return NULL;
}

static char LLMat_delete_rows_doc[] = 
"Delete rows from matrix (inplace). The rows to be deleted are specified by the mask array.\n\
\n\
Arguments:\n\
\n\
  mask: A 1D integer NumPy array. If mask[i] == 0, then row i is deleted,\n\
        otherwise row i is kept.\n\
\n\
This method may not be applied to a matrix in symmetric format.";

static PyObject*
LLMat_delete_rows(LLMatObject *self, PyObject* args){
  PyArrayObject *maskObj;
  int newm, newnnz;
  int act, row;
  
  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &maskObj)) 
    return NULL;
  if (maskObj->nd != 1 || maskObj->descr->type_num != PyArray_LONG || maskObj->dimensions[0] != self->dim[0]) {
    PyErr_SetString(PyExc_ValueError, "mask must be a 1D integer NumPy array of appropriate length");
    return NULL;
  }
  if (self->issym) {
    PyErr_SetString(SpMatrix_ErrorObject, "method not allowed for symmetric matrices");
    return NULL;
  }

  /* Delete the rows to be cancelled by rearranging the row */
  /* array. After having done so, newdim is the new matrix dim. */
  newm = 0;
  newnnz = self->nnz;
  for(row = 0; row < self->dim[0]; row ++){
    if (*(int *)(maskObj->data + row*maskObj->strides[0]) != 0){ /* This row has to be kept */
      self->root[newm] = self->root[row];
      newm ++;
    } else {			/* row let out; update free list */
      act = self->root[row];
      if(act != -1){		/* only do smth. for non-empty rows */
	newnnz --;
	while(self->link[act] != -1) { /* Walk to the end of the list */
	  act = self->link[act];
	  newnnz --;
	}
	self->link[act] = self->free;	/* Attach end of row to free list */
	self->free = self->root[row];	/* Start free list where row began */
      }
    }
  }

  /* Set the new values */
  self->dim[0] = newm;
  self->nnz = newnnz;

  Py_INCREF(Py_None); 
  return Py_None;
}

static char LLMat_delete_cols_doc[] = 
"Delete columns from matrix (inplace). The columns to be deleted are\n\
specified by the mask array.\n\
\n\
Arguments:\n\
\n\
  mask: A 1D integer NumPy array. If mask[i] == 0, then column i is deleted,\n\
        otherwise column i is kept.";

static PyObject*
LLMat_delete_cols(LLMatObject *self, PyObject* args){
  PyArrayObject *maskObj;
  int newn, newnnz;
  int act, old, col, row;
  int *shift;
  
  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &maskObj)) 
    return NULL;
  if (maskObj->nd != 1 || maskObj->descr->type_num != PyArray_LONG || maskObj->dimensions[0] != self->dim[1]) {
    PyErr_SetString(PyExc_ValueError, "mask must be a 1D integer NumPy array of appropriate length");
    return NULL;
  }
  if (self->issym) {
    PyErr_SetString(SpMatrix_ErrorObject, "method not allowed for symmetric matrices");
    return NULL;
  }

#define MASK(i) *(long *)(maskObj->data + (i)*maskObj->strides[0])

  /* Allocate column shift vector (after deletion col[i] is at */
  /* col[i] - shift[i]). */
  shift = (int*)malloc((self->dim[1])*sizeof(int));
  newn = self->dim[1];
  if (MASK(0)) shift[0] = 0; else {shift[0] = 1; newn --;}
  for (col = 1; col < self->dim[1]; col++){
    if (MASK(col)) 
      shift[col] = shift[col-1]; 
    else 
      {shift[col] = shift[col-1]+1; newn --; }
  }
    
  /* Deleteting columns in the remainig rows */
  newnnz = self->nnz;
  for(row = 0; row < self->dim[0]; row ++) {
    old = -1; act = self->root[row];
    while (act != -1){
      if (MASK(self->col[act])) {	       /* Keep this column */
	self->col[act] -= shift[self->col[act]];
	old = act; act = self->link[act];
      } else {				       /* Drop the column */
	newnnz--;                              
	if (self->root[row] == act) {	       /* Special case: first row element */
	  self->root[row] = self->link[act];
	  old = act; act = self->link[act];
	  self->link[old] = self->free;	       /* Append element into freelist */
	  self->free = old;
	} else {			       /* Regular case: element inbetween */
	  act = self->link[act];
	  self->link[self->link[old]] = self->free;
	  self->free = self->link[old];
	  self->link[old] = act;	       /* Append element into freelist */
	}
      }
    }
  }

  /* Set the new values */
  self->dim[1] = newn;
  self->nnz = newnnz;

  /* clean up */
  free(shift);

  Py_INCREF(Py_None); 
  return Py_None;

#undef MASK
}

static char LLMat_delete_rowcols_doc[] = 
"Delete rows and columns from matrix (inplace). The rows and columns to be deleted are\n\
specified by the mask array.\n\
\n\
Arguments:\n\
\n\
  mask: A 1D integer NumPy array. If mask[i] == 0, then row and column i are deleted,\n\
        otherwise row and column i are kept.";

static PyObject*
LLMat_delete_rowcols(LLMatObject *self, PyObject* args){
  PyArrayObject *maskObj;
  int newn, newm, newnnz;
  int act, old, col, row;
  int *shift;
  
  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &maskObj)) 
    return NULL;
  if (maskObj->nd != 1 || maskObj->descr->type_num != PyArray_LONG || maskObj->dimensions[0] != self->dim[0]) {
    PyErr_SetString(PyExc_ValueError, "mask must be a 1D integer NumPy array of appropriate length");
    return NULL;
  }
  if (self->dim[0] != self->dim[1]) {
    PyErr_SetString(SpMatrix_ErrorObject, "method only allowed for square matrices");
    return NULL;
  }

#define MASK(i) *(long *)(maskObj->data + (i)*maskObj->strides[0])

  /* Delete the rows to be cancelled by rearranging the row */
  /* array. After having done so, newdim is the new matrix dim. */
  newm = 0;
  newnnz = self->nnz;
  for(row = 0; row < self->dim[0]; row ++){
    if (MASK(row)){			       /* This row has to be kept */
      self->root[newm] = self->root[row];
      newm ++;
    } else {				       /* row let out; update free list */
      act = self->root[row];
      if(act != -1){			       /* only do smth. for non-empty rows */
	newnnz --;
	while(self->link[act] != -1) {	       /* Walk to the end of the list */
	  act = self->link[act];
	  newnnz --;
	}
	self->link[act] = self->free;	       /* Attach end of row to free list */
	self->free = self->root[row];	       /* Start free list where row began */
      }
    }
  }

  /* Set the new values */
  self->dim[0] = newm;
  self->nnz = newnnz;

  /* Allocate column shift vector (after deletion col[i] is at */
  /* col[i] - shift[i]). */
  shift = (int*)malloc((self->dim[1])*sizeof(int));
  newn = self->dim[1];
  if (MASK(0)) shift[0] = 0; else {shift[0] = 1; newn --;}
  for (col = 1; col < self->dim[1]; col++){
    if (MASK(col)) 
      shift[col] = shift[col-1]; 
    else 
      {shift[col] = shift[col-1]+1; newn --; }
  }
    
  /* Deleteting columns in the remainig rows */
  newnnz = self->nnz;
  for(row = 0; row < self->dim[0]; row ++) {
    old = -1; act = self->root[row];
    while (act != -1){
      if (MASK(self->col[act])) {	       /* Keep this column */
	self->col[act] -= shift[self->col[act]];
	old = act; act = self->link[act];
      } else {				       /* Drop the column */
	newnnz--;                              
	if (self->root[row] == act) {	       /* Special case: first row element */
	  self->root[row] = self->link[act];
	  old = act; act = self->link[act];
	  self->link[old] = self->free;	       /* Append element into freelist */
	  self->free = old;
	} else {			       /* Regular case: element inbetween */
	  act = self->link[act];
	  self->link[self->link[old]] = self->free;
	  self->free = self->link[old];
	  self->link[old] = act;	       /* Append element into freelist */
	}
      }
    }
  }

  /* Set the new values */
  self->dim[1] = newn;
  self->nnz = newnnz;

  /* clean up */
  free(shift);

  Py_INCREF(Py_None); 
  return Py_None;

#undef MASK
}

/*********************************/
/*  O b j e c t   m e t h o d s  */
/*********************************/

PyMethodDef LLMat_methods[] = {
  {"matvec",          (PyCFunction)LLMat_matvec,          METH_VARARGS, LLMat_matvec_doc},
  {"matvec_transp",   (PyCFunction)LLMat_matvec_transp,   METH_VARARGS, LLMat_matvec_transp_doc},
  {"to_csr",          (PyCFunction)LLMat_to_csr,          METH_VARARGS, to_csr_doc},
  {"to_sss",          (PyCFunction)LLMat_to_sss,          METH_VARARGS, to_sss_doc},
  {"generalize",      (PyCFunction)LLMat_generalize,      METH_VARARGS, LLMat_generalize_doc},
  {"compress",        (PyCFunction)LLMat_compress,        METH_VARARGS, LLMat_compress_doc},
  {"export_mtx",      (PyCFunction)LLMat_export_mtx,      METH_VARARGS, export_mtx_doc},
  {"copy",            (PyCFunction)LLMat_copy,            METH_VARARGS, copy_doc},
  {"norm",            (PyCFunction)LLMat_norm,            METH_VARARGS, LLMat_norm_doc},
  {"shift",           (PyCFunction)LLMat_shift,           METH_VARARGS, shift_doc},
  {"scale",           (PyCFunction)LLMat_scale,           METH_VARARGS, scale_doc},
  {"keys",            (PyCFunction)LLMat_keys,            METH_VARARGS, keys_doc},
  {"values",          (PyCFunction)LLMat_values,          METH_VARARGS, values_doc},
  {"items",           (PyCFunction)LLMat_items,           METH_VARARGS, items_doc},
  {"update_add_mask", (PyCFunction)LLMat_update_add_mask, METH_VARARGS, update_add_mask_doc},
  {"update_add_mask_sym", (PyCFunction)LLMat_update_add_mask_sym, METH_VARARGS, update_add_mask_sym_doc},
  {"delete_rows",     (PyCFunction)LLMat_delete_rows,     METH_VARARGS, LLMat_delete_rows_doc},
  {"delete_cols",     (PyCFunction)LLMat_delete_cols,     METH_VARARGS, LLMat_delete_cols_doc},
  {"delete_rowcols",  (PyCFunction)LLMat_delete_rowcols,  METH_VARARGS, LLMat_delete_rowcols_doc},
  {"update_add_at",  (PyCFunction)LLMat_update_add_at,  METH_VARARGS, update_add_at_doc},
  {"put",             (PyCFunction)LLMat_put,             METH_VARARGS, LLMat_put_doc},
  {"take",            (PyCFunction)LLMat_take,            METH_VARARGS, LLMat_take_doc},
  {NULL, NULL}			/* sentinel */
};

/*****************************************/
/*  L L M a t   t y p e   m e t h o d s  */
/*****************************************/

static void
LLMatType_dealloc(LLMatObject *a)
{
  PyMem_DEL(a->root);
  PyMem_DEL(a->val);
  PyMem_DEL(a->col);
  PyMem_DEL(a->link);
  PyObject_Del(a);
}

static int
LLMatType_print(LLMatObject *a, FILE *fp, int flags)
{
  int i, k, first = 1;
  char *symStr;

  if (a->issym)
    symStr = "symmetric";
  else
    symStr = "general";

  if (a->dim[1] <= PPRINT_COL_THRESH && a->dim[0] <= PPRINT_ROW_THRESH) {
    double *mat;
    int j;
    double val;
    mat = (double *)malloc(a->dim[0]*a->dim[1] * sizeof(double));
    if (mat == NULL) {
      PyErr_NoMemory();
      return -1;
    }
    fprintf(fp, "ll_mat(%s, [%d,%d]):\n", symStr, a->dim[0], a->dim[1]);
    for (i = 0; i < a->dim[0]; i ++) {
      for (j = 0; j < a->dim[1]; j ++)
	mat[i*a->dim[1] + j] = 0.0;
      k = a->root[i];
      while (k != -1) {
	mat[(i*a->dim[1])+a->col[k]] = a->val[k];
	k = a->link[k];
      }
    }

    for (i = 0; i < a->dim[0]; i ++) {
      for (j = 0; j < a->dim[1]; j ++) {
	val = mat[(i*a->dim[1])+j];
	if (val != 0.0) {
	  int exp = (int)log10(fabs(val));
	  if (abs(exp) <= 4) {
	    if (exp < 0)
	      fprintf(fp, "%9.*f ", 6, val);
	    else
	      fprintf(fp, "%9.*f ", 6-exp, val);
	  } else
	    fprintf(fp, "%9.1e ", val);
	}
	else
	  if (!(a->issym) || i > j)
	    fprintf(fp, " -------- ");
      }
      fprintf(fp, "\n");
    }
    free(mat);

  } else {

  if (a->nnz == 0) {
    fprintf(fp, "ll_mat(%s, [%d,%d])", symStr, a->dim[0], a->dim[1]);
    return 0;
  }
  fprintf(fp, "ll_mat(%s, [%d,%d], [", symStr, a->dim[0], a->dim[1]);
  for (i = 0; i < a->dim[0]; i ++) {
    k = a->root[i];
    while (k != -1) {
      if (!first)
	fprintf(fp, ", ");
      first = 0;
      fprintf(fp, "(%d,%d): %g", i, a->col[k], a->val[k]);
      k = a->link[k];
    }
  }
  fprintf(fp, "])");
  }
  return 0;
}

static PyObject *
LLMatType_getattr(LLMatObject *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->dim[0], self->dim[1]);
  if (strcmp(name, "nnz") == 0)
    return PyInt_FromLong(self->nnz);
  if (strcmp(name, "issym") == 0)
    return PyInt_FromLong(self->issym);
  if (strcmp(name, "__members__") == 0) {
    char *members[] = {"shape", "nnz", "issym"};
    int i;

    PyObject *list = PyList_New(sizeof(members)/sizeof(char *));
    if (list != NULL) {
      for (i = 0; i < sizeof(members)/sizeof(char *); i ++)
	PyList_SetItem(list, i, PyString_FromString(members[i]));
      if (PyErr_Occurred()) {
	Py_DECREF(list);
	list = NULL;
      }
    }
    return list;
  }
  return Py_FindMethod(LLMat_methods, (PyObject *)self, name);
}
/***********************************************************************
 * mapping functions
 */

/** LLMat_length - number of items in mapping
 *    == number of matrix entries
 */
static int LLMat_length(LLMatObject *self) {
  return self->dim[0] * self->dim[1];
}

static int 
slice_GetIndices(PySliceObject *r, int length, 
		 int *start, int *stop, int *step) {
  if (r->step == Py_None) {
    *step = 1;
  } else {
    if (!PyInt_Check(r->step)) return -1;
    *step = PyInt_AsLong(r->step);
  }
  if (r->start == Py_None) {
    *start = *step < 0 ? length-1 : 0;
  } else {
    if (!PyInt_Check(r->start)) return -1;
    *start = PyInt_AsLong(r->start);
    if (*start < 0) *start += length;
  }
  if (r->stop == Py_None) {
    *stop = *step < 0 ? -1 : length;
  } else {
    if (!PyInt_Check(r->stop)) return -1;
    *stop = PyInt_AsLong(r->stop);
    if (*stop < 0) *stop += length;
  }
  
  if (*start > (length-1)) *start = length;
  if (*start < 0) *start = 0;
  if (*stop < -1) *stop = -1;
  else if (*stop > length) *stop = length;
  return 0;
}

static int 
LLMat_parse_subindex(PyObject *sidx, int max, int *start, int *stop, int *step) {
  int tmp;

  if (PyInt_Check(sidx)) {
    /* parse integer */
    tmp = PyInt_AsLong(sidx);
    if (tmp < 0) tmp += max;
    if (tmp >= max || tmp < 0) {
      return -1;
    }
    *start = tmp;
    return 1;
  }
  if (PySlice_Check(sidx)) {
    /* parse slice */
    if (slice_GetIndices((PySliceObject *)sidx, max, 
			 start, stop, step) == -1) return -1;
    return 2;
  }
  /* invalid subindex */
  return -1;
}

static int 
LLMat_parse_index(PyObject *op, int dim[],
		  int *start0, int *stop0, 
		  int *start1, int *stop1) {
  PyObject *op0, *op1;
  int type0, step0, type1, step1;

  if (!PySequence_Check(op)) {
    PyErr_SetString(PyExc_IndexError, "index must be a sequence");
    return -1;
  }
  if (PySequence_Length(op) != 2) {
    PyErr_SetString(PyExc_IndexError, "There must be exactly two indices");
    return -1;
  }

  /* parse first index */
  if (!(op0 = PySequence_GetItem(op, 0))) {
    PyErr_SetString(PyExc_IndexError, "first index is invalid");
    return -1;
  }
  type0 = LLMat_parse_subindex(op0, dim[0], start0, stop0, &step0);
  Py_DECREF(op0);
  if (type0 == -1) {
    PyErr_SetString(PyExc_IndexError, "first index is invalid");
    return -1;
  }
  
  /* parse second index */
  if (!(op1 = PySequence_GetItem(op, 1))) {
    PyErr_SetString(PyExc_IndexError, "second index is invalid");
    return -1;
  }
  type1 = LLMat_parse_subindex(op1, dim[1], start1, stop1, &step1);
  Py_DECREF(op1);
  if (type1 == -1) {
    PyErr_SetString(PyExc_IndexError, "second index is invalid");
    return -1;
  }

  if (type0 == 1 && type1 == 1)
    return 1;

  /* disallow strides other than one */
  if ((type0 == 2 && step0 != 1) || (type1 == 2 && step1 != 1)) {
    PyErr_SetString(PyExc_IndexError, "strides other than one not allowed");
    return -1;
  }
  
  /* treat single row or column as submatrix */
  if (type0 == 1)
    *stop0 = *start0 + 1;
  if (type1 == 1)
    *stop1 = *start1 + 1;

  return 2;
}

/** LLMat_subscript
 *    Called when treating array object like a mapping. This is used
 *    implement two-dimensional idices, e.g. A[i,j] or A[i1:i2,j1:j2]
 */
static PyObject *
LLMat_subscript(LLMatObject *self, PyObject *idx) {
  int type, start0, stop0, start1, stop1;

  if ((type = LLMat_parse_index(idx, self->dim, &start0, &stop0, &start1, &stop1)) == -1)
    return NULL;
  if (type == 1)
    return PyFloat_FromDouble(SpMatrix_LLMatGetItem(self, start0, start1));
  else
    return get_submatrix(self, start0, stop0, start1, stop1);
}

/** LLMat_ass_subscript
 *    Called when treating array object like a mapping. This is used
 *    implement two-dimensional indices, e.g. A[i,j]
 */
static int
LLMat_ass_subscript(LLMatObject *self, PyObject *idx, PyObject *v) {
  double x;
  int type, start0, stop0, start1, stop1;
  
  if ((type = LLMat_parse_index(idx, self->dim, &start0, &stop0, &start1, &stop1)) == -1)
    return -1;
  if (type == 1) {
    if (!PyArg_Parse(v, "d;array item must be float", &x))
      return -1;
    else
      return SpMatrix_LLMatSetItem(self, start0, start1, x);
  }
  else {
    return set_submatrix(self, start0, stop0, start1, stop1, (LLMatObject *)v);
  }
}

static PyMappingMethods LLMat_as_mapping = {
    (inquiry)LLMat_length,	/*mp_length*/
    (binaryfunc)LLMat_subscript, /*mp_subscript*/
    (objobjargproc)LLMat_ass_subscript, /*mp_ass_subscript*/
};

/*************************************/
/*  L L M a t T y p e   o b j e c t  */
/*************************************/

static PyTypeObject LLMatType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "ll_mat",
  sizeof(LLMatObject),
  0,
  (destructor)LLMatType_dealloc, /* tp_dealloc */
  (printfunc)LLMatType_print,	/* tp_print */
  (getattrfunc)LLMatType_getattr, /* tp_getattr */
  0,				/* tp_setattr */
  0,				/* tp_compare */
  0,				/* tp_repr */
  0,				/* tp_as_number*/
  0,				/* tp_as_sequence*/
  &LLMat_as_mapping,		/* tp_as_mapping*/
  0,				/* tp_hash */
};

/*************************************************/
/*  M i s c .   h e l p e r   f u n c t i o n s  */
/*************************************************/

static PyObject *
SpMatrix_NewLLMatObject(int dim[], int sym, int sizeHint) {
  int i;
  LLMatObject *op;

  if (dim[0] < 0 || dim[1] < 0) {
    PyErr_SetString(PyExc_ValueError, "matrix dimension must be non-negative");
    return NULL;
  }
  if (sizeHint < 1)
    sizeHint = 1;

  /* create new SparseArrayt object */
  op = PyObject_New(LLMatObject, &LLMatType);
  if (op == NULL)
    return PyErr_NoMemory();

  op->root = NULL;
  op->val = NULL;
  op->col = NULL;
  op->link = NULL;

  /* allocate ob_val and on_idx arrays */
  op->root = PyMem_New(int, dim[0]);
  if (op->root == NULL)
    goto fail;
  op->val = PyMem_New(double, sizeHint);
  if (op->val == NULL)
    goto fail;
  op->col = PyMem_New(int, sizeHint);
  if (op->col == NULL)
    goto fail;
  op->link = PyMem_New(int, sizeHint);
  if (op->link == NULL)
    goto fail;

  /* initialize rest of fields */
  for (i = 0; i < dim[0]; i ++)
    op->root[i] = -1;
  op->dim[0] = dim[0];
  op->dim[1] = dim[1];
  op->issym = sym;
  op->nnz = 0;
  op->nalloc = sizeHint;
  op->free = -1;

  return (PyObject *) op;

 fail:
    PyMem_Del(op->link);    
    PyMem_Del(op->col);    
    PyMem_Del(op->val);    
    PyMem_Del(op->root);    
    PyObject_Del(op);
    return PyErr_NoMemory();
}

static PyObject*
LLMat_from_mtx(PyObject *module, PyObject *args) {
  LLMatObject *self = NULL;
  char *fileName;
  MM_typecode matcode;
  int dim[2], nz;
  FILE *f;
  int ret, i;
  double val;
  int row, col;

  if (!PyArg_ParseTuple(args, "s", &fileName))
    return NULL;
  
  /* open file */
  f = fopen(fileName, "r");
  if (f == NULL)
    return PyErr_SetFromErrno(PyExc_IOError);

  /* read MTX header */
  ret = mm_read_banner(f, matcode);
  if (ret != 0) {
    PyErr_SetString(PyExc_IOError, "error reading MTX file header");
    goto fail;
  }
  if (!(mm_is_real(matcode) && mm_is_matrix(matcode) &&
	mm_is_sparse(matcode))) {
    PyErr_SetString(SpMatrix_ErrorObject, "must be real, sparse matrix");
    goto fail;
  }
  ret = mm_read_mtx_crd_size(f, dim, dim+1, &nz);
  if (ret != 0) {
    PyErr_SetString(PyExc_IOError, "error reading MTX file size information");
    goto fail;
  }

  /* allocate matrix object */
  self = (LLMatObject *)SpMatrix_NewLLMatObject(dim, mm_is_symmetric(matcode), nz);
  if (self == NULL)
    goto fail;

  for (i = 0; i < nz; i ++) {
    ret = fscanf(f, "%d %d %lg\n", &row, &col, &val);
    if (ret != 3) {
      PyErr_SetString(PyExc_IOError, "error reading MTX file data");
      goto fail;
    }
    row --; col --;
    if (!(0 <= row && row < dim[0] && 0 <= col && col < dim[1])) {
      PyErr_SetString(SpMatrix_ErrorObject, "matrix indices out of range");
      fclose(f);
      return NULL;
    }
    ret = SpMatrix_LLMatSetItem(self, row, col, val);
    if (ret)
      goto fail;
  }
  fclose(f);
  return (PyObject *)self;

 fail:
  fclose(f);
  Py_XDECREF(self);
  return NULL;
}

char LLMat_matrixmultiply_doc[] = "matrixmultiply(A, B)\n\
\n\
Returns a new ll_mat object representing the matrix A*B";

static PyObject *
LLMat_matrixmultiply(PyObject *self, PyObject *args)
{
  int sizeHint = 1000;
  LLMatObject *matA, *matB, *matC;
  int dimC[2];
  int symCode, ret;
  
  if (!PyArg_ParseTuple(args, "O!O!", &LLMatType, &matA, &LLMatType, &matB))
    return NULL;

  /* matrix dimensions
   */
  dimC[0] = matA->dim[0];
  dimC[1] = matB->dim[1];

  if (matA->dim[1] != matB->dim[0]) {
    PyErr_SetString(PyExc_ValueError, "matrix dimensions must agree");
    return NULL;
  }

  /* create result object
   */
  matC = (LLMatObject *)SpMatrix_NewLLMatObject(dimC, 0, sizeHint);
  if (matC == NULL)
    return NULL;

  symCode = matB->issym << 1 | matA->issym;
  if (symCode == 0) {
    /* unsym-unsym multiplication
     */

#if !OPT_MATMATMUL
    double valA;
    int iA, jA, kA, kB; 

    for (iA = 0; iA < matA->dim[0]; iA ++) {
      kA = matA->root[iA];
      while (kA != -1) {
	valA = matA->val[kA];
	jA = matA->col[kA];
	kA = matA->link[kA];
	
	/* add jA-th row of B to iA-th row of C */
	kB = matB->root[jA];
	while (kB != -1) {
	  ret = SpMatrix_LLMatUpdateItemAdd(matC, iA, matB->col[kB], valA*matB->val[kB]);
	  if (ret == -1)
	    goto fail;
	  kB = matB->link[kB];
	}
      }

    }
#else
    int *tmpage = NULL;
    int *tmpind = NULL;
    int *tmpcol = NULL;
    double *tmpval = NULL;
    int tmpsize;
    int nxttmp;
    int row;
    int indA, colA, colB, dummy, indB;
    double valA;
    
    tmpsize = 5; nxttmp = -1;
    tmpage = (int*)malloc(matB->dim[1] * sizeof(int));
    tmpind = (int*)malloc(matB->dim[1] * sizeof(int));
    tmpcol = (int*)malloc(tmpsize * sizeof(int));
    tmpval = (double*)malloc(tmpsize * sizeof(double));
    if (tmpage == NULL || tmpind == NULL || tmpcol == NULL ||tmpval == NULL) {
      PyErr_NoMemory();
      goto fail_unsym_unsym;
    }

    /* main loop */
    
    for(row=0; row < matB->dim[1]; row++){ tmpage[row] = -1;}

    /* Go through each row of A and perform necessary computations */
    for(row=0; row < matA->dim[0]; row++) {
      indA = matA->root[row];         /* Pick first entry of A[row,:] */
      while(indA != -1){              /* As long as there is an element in A[row,:] ... */
	colA = matA->col[indA];       /* ... get its column number ... */
	valA = matA->val[indA];       /* ... and value ... */
	
	indB = matB->root[colA];       /* colA is equivalent to rowB! */
	while(indB != -1){
	  colB = matB->col[indB];
	  
	  if(tmpage[colB] != row){         /* This column never appeared so far */
	    nxttmp++;
	    tmpage[colB]  = row;   
	    tmpind[colB]  = nxttmp;
	    
	    if(nxttmp >= tmpsize){          /* If tmp storage is too small, realloc */
	      tmpsize = (int)((tmpsize*12)/10)+1;
	      tmpcol = (int*)realloc(tmpcol, tmpsize * sizeof(int));
	      tmpval = (double*)realloc(tmpval, tmpsize * sizeof(double));
	      if (tmpcol == NULL ||tmpval == NULL) {
		PyErr_NoMemory();
		goto fail_unsym_unsym;
	      }
	    }
	    
	    tmpcol[nxttmp] = colB;
	    tmpval[nxttmp] = valA * matB->val[indB];
	  }else{                   /* This column appeared at least once already */
	    dummy = tmpind[colB];
	    tmpval[dummy] += valA * matB->val[indB];
	  }
	  
	  indB = matB->link[indB];
	}
	indA = matA->link[indA];
      }
	
      /* All the new values for rowC = rowA have now to be filled in */
      /* into the matrix C */
      for(dummy=0; dummy<=nxttmp; dummy++) {
	if (SpMatrix_LLMatSetItem(matC,row,tmpcol[dummy],tmpval[dummy]))
	  goto fail_unsym_unsym;
      }
      
      nxttmp=-1; /* For the next row of A we need a "fresh" tmp storage */
    }
    /* Get the memory back ... */
    free(tmpage);
    free(tmpind);
    free(tmpcol);
    free(tmpval);
    return (PyObject *)matC;

  fail_unsym_unsym:
    free(tmpage);
    free(tmpind);
    free(tmpcol);
    free(tmpval);
    goto fail;
#endif

  } else if (symCode == 1) {
    
    /* sym-unsym multiplication
     */
    double valA;
    int iA, jA, kA, kB;
    
    for (iA = 0; iA < matA->dim[0]; iA ++) {
      kA = matA->root[iA];
      while (kA != -1) {
	valA = matA->val[kA];
	jA = matA->col[kA];
	kA = matA->link[kA];

	/* add jA-th row of B to iA-th row of C */
	kB = matB->root[jA];
	while (kB != -1) {
	  ret = SpMatrix_LLMatUpdateItemAdd(matC, iA, matB->col[kB], valA*matB->val[kB]);
	  if (ret == -1)
	    goto fail;
	  kB = matB->link[kB];
	}
	
	if (iA == jA)
	  continue;
	
	/* add iA-th row of B to jA-th row of C */
	kB = matB->root[iA];
	while (kB != -1) {
	  ret = SpMatrix_LLMatUpdateItemAdd(matC, jA, matB->col[kB], valA*matB->val[kB]);
	  if (ret == -1)
	    goto fail;
	  kB = matB->link[kB];
	}	
      }
    }

  } else if (symCode == 2) {

    /* unsym-sym multiplication
     */

    PyErr_SetString(PyExc_NotImplementedError, "multiply of an unsymmetric and a symmetric matrix not supported");
    goto fail;

  } else {

    /* sym-sym multiplication
     */

    PyErr_SetString(PyExc_NotImplementedError, "multiply of two symmetric matrices not supported");
    goto fail;

  }
  
  return (PyObject *)matC;

 fail:
  Py_DECREF(matC);
  return NULL;
}

static char LLMat_dot_doc[] = "dot(A, B)\n\
\n\
Returns a new ll_mat object representing the matrix transpose(A)*B";

static PyObject *
LLMat_dot(PyObject *self, PyObject *args)
{
  int sizeHint = 1000;
  LLMatObject *matA, *matB, *matC;
  int dimC[2];
  double valA;
  int iA, kA, iC, kB, ret;
  
  if (!PyArg_ParseTuple(args, "O!O!", &LLMatType, &matA, &LLMatType, &matB))
    return NULL;

  dimC[0] = matA->dim[1];
  dimC[1] = matB->dim[1];

  if (matA->dim[0] != matB->dim[0]) {
    PyErr_SetString(PyExc_ValueError, "matrix dimensions must agree");
    return NULL;
  }

  if (matA->issym || matB->issym) {
    PyErr_SetString(PyExc_NotImplementedError, "ddot operation with symmetric matrices not supported");
    return NULL;
  }

  matC = (LLMatObject *)SpMatrix_NewLLMatObject(dimC, 0, sizeHint);
  if (matC == NULL)
    return NULL;

  for (iA = 0; iA < matA->dim[0]; iA ++) {
    kA = matA->root[iA];
    while (kA != -1) {
      valA = matA->val[kA];
      iC = matA->col[kA];
      kB = matB->root[iA];
      while (kB != -1) {
	ret = SpMatrix_LLMatUpdateItemAdd(matC, iC, matB->col[kB], valA*matB->val[kB]);
	if (ret == -1)
	  goto fail;
	kB = matB->link[kB];
      }
      kA = matA->link[kA];
    }
  }
  return (PyObject *)matC;

 fail:
  Py_DECREF(matC);
  return NULL;
}
