/* numpyio.c -- Version 0.9.9
 * 
 * Author:  Travis E. Oliphant
 * Date  :  March 1999
 * 
 * This file is a module for python that defines basically two functions for
 * reading from and writing to a binary file.  It also has some functions
 * for byteswapping data and packing and unpacking bits.
 *
 *  The data goes into a NumPy array object (multiarray)
 *
 * It is basically an implemetation of read and write with the data
 *  going directly into a NumPy array
 * 
 * Permission is granted to use this program however you see fit, but I give 
 *   no guarantees as to its usefulness or reliability.  You assume full
 *   responsibility for using this program.
 *
 *  Thanks to Michael A. Miller <miller5@uiuc.edu> 
 *    whose TableIO packages helped me learn how
 *    to write an extension package.  I've adapted his Makefile as well.
 */

#include "Python.h"             /* Python header files */
#include "Numeric/arrayobject.h"
/* #include <math.h> */
#include <stdio.h>

void rbo(char *, int, int);
void packbits(char *, int, char *, int, int);
void unpackbits(char *, int, char *, int, int, int);
int is_little_endian(void);

static PyObject *ErrorObject;     /* locally-raised exception */

#define PYERR(message) do {PyErr_SetString(PyExc_ValueError, message); goto fail;} while(0)
#define DATA(arr) ((arr)->data)
#define DIMS(arr) ((arr)->dimensions)
#define STRIDES(arr) ((arr)->strides)
#define ELSIZE(arr) ((arr)->descr->elsize)
#define OBJECTTYPE(arr) ((arr)->descr->type_num)
#define BASEOBJ(arr) ((PyArrayObject *)((arr)->base))
#define RANK(arr) ((arr)->nd)
#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)
#define MIN(a,b) (((a) > (b)) ? (b) : (a))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define PYSETERROR(message) \
{ PyErr_SetString(ErrorObject, message); goto fail; }

#define INCREMENT(ret_ind, nd, max_ind) \
{ \
  int k; \
  k = (nd) - 1; \
  if (++(ret_ind)[k] >= (max_ind)[k]) { \
    while (k >= 0 && ((ret_ind)[k] >= (max_ind)[k]-1)) \
      (ret_ind)[k--] = 0; \
    if (k >= 0) (ret_ind)[k]++; \
    else (ret_ind)[0] = (max_ind)[0]; \
  }  \
}

#define CALCINDEX(indx, nd_index, strides, ndim) \
{ \
  int i; \
 \
  indx = 0; \
  for (i=0; i < (ndim); i++)  \
    indx += nd_index[i]*strides[i]; \
} 

static PyObject *
 numpyio_fromfile(PyObject *self, PyObject *args)  /* args: number of bytes and type */
{
  PyObject *file;
  PyArrayObject *arr=NULL;
  PyArray_Descr *indescr=NULL;
  void     *ibuff=NULL;
  int      myelsize;
  int      ibuff_cleared = 1;
  long      n,nread;
  char      read_type;
  FILE     *fp;
  char      dobyteswap = 0;
  int      swap_factor;
  char     out_type = 124;    /* set to unused value */

  if (!PyArg_ParseTuple( args, "Olc|cb" , &file, &n, &read_type, &out_type, &dobyteswap ))
    return NULL;

  if (out_type == 124)
    out_type = read_type;

  fp = PyFile_AsFile(file);

  if (fp == NULL) {
    PYSETERROR("First argument must be an open file");
  } 

  if (n <= 0) {
    PYSETERROR("Second argument (number of bytes to read) must be positive.");
  }
  /* Make a 1-D NumPy array of type read_type with n elements */ 

  if ((arr = (PyArrayObject *)PyArray_FromDims(1,(int*)&n,out_type)) == NULL)
    return NULL;
  
      /* Read the data into the array from the file */
  if (out_type == read_type) {
    ibuff = arr -> data;
    myelsize = arr -> descr -> elsize;
  }
  else {                    /* Alocate a storage buffer for data read in */
    indescr = PyArray_DescrFromType((int ) read_type);
    if (indescr == NULL) goto fail;
    myelsize = indescr -> elsize;
    ibuff = malloc(myelsize*n);
    if (ibuff == NULL)
      PYSETERROR("Could not allocate memory for type casting")
	ibuff_cleared = 0;
  }
  
  nread = fread(ibuff,myelsize,n,fp);
  if (ferror(fp)) {
    clearerr(fp);
    PYSETERROR("There was an error reading from the file");
  }
  
  /* Check to see correct number of bytes were read.  If not, then
     resize the array to the number of bytes actually read in.
  */

  if (nread < n) {
    fprintf(stderr,"Warning: %ld bytes requested, %ld bytes read.\n", n, nread);
    arr->dimensions[0] = nread;
    arr->data = realloc(arr->data,arr->descr->elsize*nread);
  }
  
  if (dobyteswap) {
    swap_factor = ((read_type=='F' || read_type=='D') ? 2 : 1);
    rbo(ibuff,myelsize/swap_factor,nread*swap_factor);
  }
  
  if (out_type != read_type) {    /* We need to type_cast it */
    (indescr->cast[arr->descr->type_num])(ibuff, 1, arr->data, 1, nread );
    free(ibuff);
    ibuff_cleared = 1;
  }

  return PyArray_Return(arr);

 fail:
  if (!ibuff_cleared) free(ibuff);
  Py_XDECREF(arr);
  return NULL;

}

static int write_buffered_output(FILE *fp, PyArrayObject *arr, PyArray_Descr* outdescr, char *buffer, int buffer_size, int bswap) {

  /* INITIALIZE N-D index */

  /* Loop over the N-D index filling the buffer with the data in arr
        (indexed correctly using strides)
     Each time dimension subdim is about to roll
     write the buffer to disk and fill it again. */

  char  *buff_ptr, *output_ptr;
  int nwrite, *nd_index, indx;
  int buffer_size_bytes, elsize;

  buff_ptr = buffer;
  nd_index = (int *)calloc(arr->nd,sizeof(int));
  if (NULL == nd_index) {
     PyErr_SetString(ErrorObject,"Could not allocate memory for index array.");
     return -1;
  }
  buffer_size_bytes = buffer_size * arr->descr->elsize;
  while(nd_index[0] != arr->dimensions[0]) {
    CALCINDEX(indx,nd_index,arr->strides,arr->nd);
    memcpy(buff_ptr, arr->data+indx, arr->descr->elsize);
    buff_ptr += arr->descr->elsize;
    INCREMENT(nd_index,arr->nd,arr->dimensions);
    if ((buff_ptr - buffer) >= buffer_size_bytes) {
      buff_ptr = buffer;

      if (outdescr->type != arr->descr->type) {  /* Cast to new type before writing */
	output_ptr = buffer + buffer_size_bytes;
        (arr->descr->cast[outdescr->type_num])(buffer, 1, output_ptr, 1, buffer_size);
	elsize = outdescr->elsize;
      }
      else {
	output_ptr = buffer;
	elsize = arr->descr->elsize;
      }
      if (bswap) {
	rbo((char *)output_ptr, elsize, buffer_size);
      }

      nwrite = fwrite(output_ptr, elsize, buffer_size, fp);

      if (ferror(fp)) {
	clearerr(fp);
	PyErr_SetString(ErrorObject,"There was an error writing to the file");
	return -1;
      }
      if (nwrite < buffer_size) {
	fprintf(stderr,"Warning: %d of %d specified bytes written.\n",nwrite, buffer_size);
      }
    }    

  }
  return 0;
}

static PyObject *
 numpyio_tofile(PyObject *self, PyObject *args)  /* args: number of bytes and type */
{
  PyObject *file;
  PyArrayObject *arr = NULL;
  PyObject *obj;
  PyArray_Descr *outdescr;
  void     *obuff = NULL;
  long      n, k, nwrite, maxN, elsize_bytes;
  int      myelsize, buffer_size;
  FILE     *fp;
  char     *buffer = NULL;
  char     dobyteswap = 0;
  int      swap_factor;
  char     ownalloc = 0;
  char     write_type = 124;

  if (!PyArg_ParseTuple( args, "OlO|cb" , &file, &n, &obj, &write_type, &dobyteswap))
    return NULL;
  
  fp = PyFile_AsFile(file);

  if (fp == NULL) {
    PYSETERROR("First argument must be an open file");
  }

  if (!PyArray_Check(obj)) {
    PYSETERROR("Third argument must be a NumPy array.");
  }

  maxN = PyArray_SIZE((PyArrayObject *)obj);
  if (n > maxN)
    PYSETERROR("The NumPy array does not have that many elements.");

  if (((PyArrayObject *)obj)->descr->type_num == PyArray_OBJECT)
    PYSETERROR("Cannot write an object array.");

  if (!PyArray_ISCONTIGUOUS((PyArrayObject *)obj)) {
    arr = (PyArrayObject *)PyArray_CopyFromObject(obj,((PyArrayObject *)obj) -> descr -> type_num, 0, 0); 
    if (NULL == arr) { /* Memory allocation failed 
			 Write out buffered data using strides info */
      arr = (PyArrayObject *)obj;
      Py_INCREF(arr);
      if (write_type == 124)
	write_type = arr -> descr -> type;
      
      if (write_type != arr -> descr -> type) {
	outdescr = PyArray_DescrFromType((int) write_type);
	if (outdescr == NULL) goto fail;
	elsize_bytes = (outdescr->elsize + arr->descr->elsize); /* allocate space for buffer and casted buffer */
      }
      else {
	outdescr = arr->descr;
	elsize_bytes = (arr->descr->elsize);
      }
      k = 0;
      do {
	k++;
	buffer_size = _PyArray_multiply_list(arr->dimensions + k, arr->nd - k);
	buffer = (char *)malloc(elsize_bytes*buffer_size);
      }
      while ((NULL == buffer) && (k < arr->nd - 1));

      if (NULL == buffer)  /* Still NULL no size was small enough */
	PYSETERROR("Could not allocate memory for any attempted output buffer size.");

      /* Write a buffered output */

      if (write_buffered_output(fp, (PyArrayObject *)obj, outdescr, buffer, buffer_size, dobyteswap) < 0) {
	free(buffer);
	goto fail;
      }
      free(buffer);
      Py_DECREF(arr);
      Py_INCREF(Py_None);
      return Py_None;
    }
  }
  else {
    arr = (PyArrayObject *)obj;
    Py_INCREF(arr);
  }

  /* Write the array to file (low-level data transfer) */
  if (n > 0) {
    
    if (write_type == 124)   /* Wasn't specified:  use input type */
      write_type = arr -> descr -> type;    
    
    if (write_type == arr -> descr -> type) {  /* point output buffer to data */
      obuff = arr -> data;
      myelsize = arr -> descr -> elsize;
    }
    else {
      if ((outdescr = PyArray_DescrFromType((int ) write_type)) == NULL) goto fail;
      myelsize = outdescr -> elsize;
      obuff = malloc(n*myelsize);
      if (obuff == NULL)
	PYSETERROR("Could not allocate memory for type-casting");
      ownalloc = 1;
      (arr->descr->cast[(int)outdescr->type_num])(arr->data,1,obuff,1,n);
    }      
    /* Write the data from the array to the file */
    if (dobyteswap) {
      swap_factor = ((write_type=='F' || write_type=='D') ? 2 : 1);
      rbo((char *)obuff,myelsize/swap_factor,n*swap_factor); 
    }
    
    nwrite = fwrite(obuff,myelsize,n,fp);
    
    if (dobyteswap) {       /* Swap data in memory back if allocated obuff */
      if (write_type == arr -> descr -> type)  /* otherwise we changed obuff only */
	rbo(arr->data,arr->descr->elsize/swap_factor,PyArray_SIZE(arr)*swap_factor);
    }
    
    if (ferror(fp)) {
      clearerr(fp);
      PYSETERROR("There was an error writing to the file");
    }
    if (nwrite < n) {
      fprintf(stderr,"Warning: %ld of %ld specified bytes written.\n",nwrite,n);
    }
  }

  if (ownalloc == 1) {
    free(obuff); 
  }

  Py_DECREF(arr);
  Py_INCREF(Py_None);
  return Py_None;

 fail:
  if (ownalloc == 1) free(obuff);
  Py_XDECREF(arr);
  return NULL;

}

static PyObject *
 numpyio_byteswap(PyObject *self, PyObject *args)  /* args: number of bytes and type */
{
  PyArrayObject *arr = NULL;
  PyObject *obj;
  int type;

  if (!PyArg_ParseTuple( args, "O" , &obj))
    return NULL;
  
  type = PyArray_ObjectType(obj,0);
  if ((arr = (PyArrayObject *)PyArray_ContiguousFromObject(obj,type,0,0)) == NULL)
    return NULL;

  rbo(arr->data,arr->descr->elsize,PyArray_SIZE(arr));

  return PyArray_Return(arr);
}

static PyObject *
 numpyio_pack(PyObject *self, PyObject *args)  /* args: in */
{
  PyArrayObject *arr = NULL, *out = NULL;
  PyObject *obj;
  int      els_per_slice;
  int      out_size;
  int      type;

  if (!PyArg_ParseTuple( args, "O" , &obj))
    return NULL;
  
  type = PyArray_ObjectType(obj,0);
  if ((arr = (PyArrayObject *)PyArray_ContiguousFromObject(obj,type,0,0)) == NULL)
    return NULL;

  if (arr->descr->type_num > PyArray_LONG)
    PYSETERROR("Expecting an input array of integer type (no floats).");

  /* Get size information from input array and make a 1-D output array of bytes */

  els_per_slice = arr->dimensions[arr->nd - 1];
  if (arr->nd > 1) 
     els_per_slice =  els_per_slice * arr->dimensions[arr->nd - 2]; 

  out_size = (PyArray_SIZE(arr)/els_per_slice)*ceil ( (float) els_per_slice / 8);

  if ((out = (PyArrayObject *)PyArray_FromDims(1,&out_size,PyArray_UBYTE))==NULL) {
      goto fail;
  }
  
  packbits(arr->data,arr->descr->elsize,out->data,PyArray_SIZE(arr),els_per_slice);

  Py_DECREF(arr);
  return PyArray_Return(out);

 fail:
  Py_XDECREF(arr);
  return NULL;

}

static PyObject *
 numpyio_unpack(PyObject *self, PyObject *args)  /* args: in, out_type */
{
  PyArrayObject *arr = NULL, *out=NULL;
  PyObject *obj;
  int      els_per_slice, arrsize;
  int      out_size, type;
  char     out_type = 'b';

  if (!PyArg_ParseTuple( args, "Oi|c" , &obj, &els_per_slice, &out_type))
    return NULL;
  
  if (els_per_slice < 1)
    PYSETERROR("Second argument is elements_per_slice and it must be >= 1.");

  type = PyArray_ObjectType(obj,0);
  if ((arr = (PyArrayObject *)PyArray_ContiguousFromObject(obj,type,0,0)) == NULL)
    return NULL;

  arrsize = PyArray_SIZE(arr);

  if ((arrsize % (int) (ceil( (float) els_per_slice / 8))) != 0)
    PYSETERROR("That cannot be the number of elements per slice for this array size.");

  if (arr->descr->type_num > PyArray_LONG)
    PYSETERROR("Can only unpack arrays that are of integer type.");

  /* Make an 1-D output array of type out_type */

  out_size = els_per_slice * arrsize / ceil( (float) els_per_slice / 8);

  if ((out = (PyArrayObject *)PyArray_FromDims(1,&out_size,out_type))==NULL)
      goto fail;

  if (out->descr->type_num > PyArray_LONG) {
    PYSETERROR("Can only unpack bits into integer type.");
  }
  
  unpackbits(arr->data,arr->descr->elsize,out->data,out->descr->elsize,out_size,els_per_slice);

  Py_DECREF(arr);
  return PyArray_Return(out);

 fail:
  Py_XDECREF(out);
  Py_XDECREF(arr);
  return NULL;
}


static char fread_doc[] = 
"g = numpyio.fread( fid, Num, read_type { mem_type, byteswap})\n\n"
"     fid =       open file pointer object (i.e. from fid = open('filename') )\n"
"     Num =       number of elements to read of type read_type\n"
"     read_type = a character in 'cb1silfdFD' (PyArray types)\n"
"                 describing how to interpret bytes on disk.\nOPTIONAL\n"
"     mem_type =  a character (PyArray type) describing what kind of\n"
"                 PyArray to return in g.   Default = read_type\n"
"     byteswap =  0 for no byteswapping or a 1 to byteswap (to handle\n"
"                 different endianness).    Default = 0.";

static char fwrite_doc[] = 
"numpyio.fwrite( fid, Num, myarray { write_type, byteswap} )\n\n"
"     fid =       open file stream\n"
"     Num =       number of elements to write\n"
"     myarray =   NumPy array holding the data to write (will be\n"
"                 written as if ravel(myarray) was passed)\nOPTIONAL\n"
"     write_type = character ('cb1silfdFD') describing how to write the\n"
"                  data (what datatype to use)  Default = type of\n"
"                  myarray.\n"
"     byteswap =   0 or 1 to determine if byteswapping occurs on write.\n"
"                  Default = 0.";

static char bswap_doc[] = 
"     out = numpyio.bswap(myarray)\n\n"
"     myarray = an array whose elements you want to byteswap.\n"
"     out     = a reference to byteswapped myarray.\n\n"
"     This does an inplace byte-swap so that myarray is changed in\n"
"     memory.";

static char packbits_doc[] = 
"out = numpyio.packbits(myarray)\n\n"
"  myarray = an array whose (assumed binary) elements you want to\n"
"             pack into bits (must be of integer type, 'cb1sl')\n\n"
"   This routine packs the elements of a binary-valued dataset into a\n"
"   1-D NumPy array of type PyArray_UBYTE ('b') whose bits correspond to\n"
"   the logical (0 or nonzero) value of the input elements. \n\n"
"   If myarray has more dimensions than 2 it packs each slice (rows*columns)\n"
"   separately.  The number of elements per slice (rows*columns) is\n"
"   important to know to be able to unpack the data later.\n\n"
"     Example:\n"
"     >>> a = array([[[1,0,1],\n"
"     ...             [0,1,0]],\n"
"     ...            [[1,1,0],\n"
"     ...             [0,0,1]]])\n"
"     >>> b = numpyio.packbits(a)\n"
"     >>> b\n"
"     array([168, 196], 'b')\n\n"
"     Note that 168 = 128 + 32 + 8\n"
"               196 = 128 + 64 + 4";

static char unpackbits_doc[] = 
"out = numpyio.unpackbits(myarray, elements_per_slice {, out_type} )\n\n"
"     myarray =        Array of integer type ('cb1sl') whose least\n"
"                      significant byte is a bit-field for the\n"
"                      resulting output array.\n\n"
"     elements_per_slice = Necessary for interpretation of myarray.\n"
"                          This is how many elements in the\n "
"                         rows*columns of original packed structure.\n\nOPTIONAL\n"
"     out_type =       The type of output array to populate with 1's\n"
"                      and 0's.  Must be an integer type.\n\n\nThe output array\n"
"                      will be a 1-D array of 1's and zero's";


#define BUFSIZE 256
/* Convert a Python string object to a complex number */
static int convert_from_object(PyObject *obj, Py_complex *cnum)
{
  PyObject *res=NULL, *elobj=NULL;
  PyObject *newstr=NULL, *finalobj=NULL, *valobj=NULL;
  char strbuffer[2*BUFSIZE];
  char *xptr, *elptr;
  char *newstrbuff, thischar;
  char buffer[BUFSIZE];
  char validnum[] = "0123456789.eE+-";
  int validlen = 15;
  int inegflag = 1;
  int rnegflag = 1;
  int n, k, m, i, elN, size, state, count;
  double val;

  if (!PyString_Check(obj)) return -1;

  /* strip string */
  newstr = PyObject_CallMethod(obj, "strip", NULL);
  if (newstr == NULL) goto fail;

  /* Replace any 'e+' or 'e-' */
  size = PyString_GET_SIZE(newstr);
  newstrbuff = PyString_AsString(newstr);
  if (newstrbuff == NULL) goto fail;
  if (size > 2*BUFSIZE) PYERR("String too large.");

  state = 0;
  count = 0;
  for (k=0; k<size; k++) {
    thischar = newstrbuff[k];
    if (state == 1) {
      if (thischar == '+') {
	thischar = '\254';
      }
      else if (thischar == '-') {
	thischar = '\253';
      }
    }
    if ((thischar == 'e') || (thischar == 'E')) state = 1;
    else state = 0;
    strbuffer[count] = thischar;
    count++;
  }
  Py_DECREF(newstr);
  newstr = PyString_FromStringAndSize(strbuffer, count);
  if (newstr == NULL) goto fail;
  xptr = strbuffer;
    
  /* Split the string into two substrings first on a ',' then on a '+'
     or '-' */
  res = PyObject_CallMethod(newstr, "split", "s", ",");
  if (res == NULL) goto fail;
  if (PySequence_Size(res) < 2) {
    Py_DECREF(res);
    res = PyObject_CallMethod(newstr, "split", "s", "+");
    if (res == NULL) goto fail;
  }
  if (PySequence_Size(res) < 2) {
    if ((strbuffer[0] == '(') || (strbuffer[0] == '[') || 
	(strbuffer[0] == '{')) {
      xptr++;
      count--;
      /* strip leading whitespaces */
      while (isspace(*xptr)) {xptr++; count--;}
    }
    if (xptr[0] == '-') {
      rnegflag = -1;
      xptr++;
      count--;
    }
    Py_DECREF(newstr);
    newstr = PyString_FromStringAndSize(xptr, count);
    if (newstr == NULL) goto fail;
    Py_DECREF(res);
    res = PyObject_CallMethod(newstr, "split", "s", "-");
    if (res == NULL) goto fail;
    inegflag = -1;
  }

  size = PySequence_Size(res);
  for (k=0; k < MIN(size,2); k++) {
    elobj = PySequence_GetItem(res, k);
    if (elobj == NULL) goto fail;
    elN = PyString_Size(elobj);
    if ((elN > BUFSIZE))
	PYSETERROR("String too large.");

    /* Replace back the + and - and strip away invalid characters */
    elptr = PyString_AsString(elobj);
    m = 0;
    for (n=0; n < elN; n++) {
      thischar = elptr[n];
      if (thischar == '\254')
	buffer[m++] = '+';
      else if (thischar == '\253')
	buffer[m++] = '-';
      else {
	for (i=0; i< validlen; i++) {
	  if (thischar == validnum[i]) break;
	}
	if (i < validlen) buffer[m++] = thischar;
      }
    }
    finalobj = PyString_FromStringAndSize(buffer, m);
    if (finalobj == NULL) goto fail;
    valobj = PyFloat_FromString(finalobj, NULL);  /* Try to make a float */
    if (valobj == NULL) goto fail;
    val = PyFloat_AsDouble(valobj);
    if (PyErr_Occurred()) goto fail;
    Py_DECREF(finalobj);
    Py_DECREF(valobj);
    Py_DECREF(elobj);
    if (k==0) {
      cnum->real = val*rnegflag;
    }
    else {
      cnum->imag = val*inegflag;
    }
    
  }
  Py_DECREF(newstr);
  Py_DECREF(res);
  return 0;
  
 fail:
  Py_XDECREF(res);  
  Py_XDECREF(elobj);
  Py_XDECREF(newstr);
  Py_XDECREF(finalobj);
  Py_XDECREF(valobj);
  return -1;
}



static int PyTypeFromChar(char ctype)
{
  switch(ctype) {
  case 'c': return PyArray_CHAR;
  case 'b': return PyArray_UBYTE;
  case '1': return PyArray_SBYTE;
  case 's': return PyArray_SHORT;
  case 'i': return PyArray_INT;
#ifdef PyArray_UNSIGNED_TYPES
  case 'u': return PyArray_UINT;
  case 'w': return PyArray_USHORT;
#endif
  case 'l': return PyArray_LONG;
  case 'f': return PyArray_FLOAT;
  case 'd': return PyArray_DOUBLE;
  case 'F': return PyArray_CFLOAT;
  case 'D': return PyArray_CDOUBLE;
  case 'O': return PyArray_OBJECT; 
  }
  return PyArray_NOTYPE;
}


static PyObject *
 numpyio_convert_objects(PyObject *self, PyObject *args)
{
  PyObject *obj = NULL, *missing_val = NULL;
  PyArrayObject *arr = NULL, *out=NULL;
  PyArrayObject *missing_arr = NULL;
  PyArray_Descr *descr;
  PyObject *builtins, *dict;
  char out_type;
  int int_type, i, err;
  char *outptr;
  PyObject **arrptr;
  PyObject *numobj=NULL;
  PyObject *comp_obj;
  Py_complex numc;
  PyArray_VectorUnaryFunc *funcptr;

  if (!PyArg_ParseTuple( args, "Oc|O" , &obj, &out_type, &missing_val))
    return NULL;

  if (missing_val == NULL) {
    missing_val = PyInt_FromLong(0);
  }
  else {
    Py_INCREF(missing_val);  /* Increment missing_val for later DECREF */
  }

  int_type = PyTypeFromChar(out_type);
  if ((int_type == PyArray_NOTYPE) || (int_type == PyArray_OBJECT))
    PYERR("Invalid output type.");

  missing_arr = (PyArrayObject *)PyArray_ContiguousFromObject(missing_val, 
							      int_type, 0, 0);
  Py_DECREF(missing_val);
  missing_val = NULL;  /* So later later failures don't decrement it */

  if ((missing_arr == NULL)) goto fail;
  if ((RANK(missing_arr) > 0)) PYERR("Missing value must be as scalar");

  arr = (PyArrayObject *)PyArray_ContiguousFromObject(obj, PyArray_OBJECT, 
						      0, 0);
  if (arr == NULL) goto fail;

  out = (PyArrayObject *)PyArray_FromDims(RANK(arr), DIMS(arr), int_type);
  if (out == NULL) goto fail;

  /* Get the builtin_functions from the builtin module */
  builtins = PyImport_AddModule("__builtin__");
  if (builtins == NULL) goto fail;

  dict = PyModule_GetDict(builtins);
  comp_obj = PyDict_GetItemString(dict, "complex");
  if (comp_obj == NULL) goto fail;
  
  /*  get_complex = PyDict_GetItemString(dict, "complex");
  get_float = PyDict_GetItemString(dict, "float");
  get_int = PyDict_GetItemString(dict, "int");
  if ((get_complex == NULL) || (get_float == NULL) || (get_int == NULL) ) goto fail;
  */
  /* 
  get_complex_self = PyCFunction_GetSelf(PyDict_GetItemString(dict, "complex"));
  get_float_self = PyCFunction_GetSelf(PyDict_GetItemString(dict, "float"));
  get_int_self = PyCFunction_GetSelf(PyDict_GetItemString(dict, "int"));
  */

  /* Loop through arr and convert each element and place in out */
  i = PyArray_Size((PyObject *)arr);
  arrptr = ((PyObject **)DATA(arr)) - 1;
  outptr = (DATA(out)) - ELSIZE(out);

  descr = PyArray_DescrFromType(PyArray_CDOUBLE);
  funcptr = descr->cast[int_type];

  while (i--) {
    outptr += ELSIZE(out);
    arrptr  += 1;
    numc.real = 0;
    numc.imag = 0;
    numobj = PyObject_CallFunction(comp_obj, "O", *arrptr);
    if (numobj != NULL) {
      numc = PyComplex_AsCComplex(numobj);
      Py_DECREF(numobj);
    }
    if (PyErr_Occurred()) {   /* Use our own homegrown converter... */
	PyErr_Clear();
	err = convert_from_object(*arrptr, &numc);
	if (PyErr_Occurred()) PyErr_Clear();
	if (err < 0) {     /* Nothing works fill with missing value... */
	    memcpy(outptr, DATA(missing_arr), ELSIZE(out));
	}
    }
    /* Place numc into the array */
    funcptr((void *)&(numc.real), 1, (void *)outptr, 1, 1);
  }

  Py_DECREF(missing_arr);
  Py_DECREF(arr);
  return PyArray_Return(out);
  
 fail:
  Py_XDECREF(out);
  Py_XDECREF(arr); 
  Py_XDECREF(missing_arr);
  Py_XDECREF(missing_val);
  return NULL;
}


static char convert_objects_doc[] = 
"convert_objectarray(myarray, arraytype{, missing_value} ) -> out \n\n"
"    myarray = Sequence of strings.\n"
"    arraytype = Type of output array.\n"
"    missing_value = Value to insert when conversion fails.";

/* *************************************************************************** */
/* Method registration table: name-string -> function-pointer */

static struct PyMethodDef numpyio_methods[] = {
  {"fread",     numpyio_fromfile,   1, fread_doc},
  {"fwrite",    numpyio_tofile,     1, fwrite_doc},
  {"bswap",     numpyio_byteswap,   1, bswap_doc},
  {"packbits",  numpyio_pack,       1, packbits_doc},
  {"unpackbits", numpyio_unpack,     1, unpackbits_doc},
  {"convert_objectarray", numpyio_convert_objects, 1, convert_objects_doc},
  {NULL,         NULL}
};

DL_EXPORT(void) initnumpyio(void)
{
  PyObject *m, *d;

  import_array();   /* allows multiarray to be a shared library (I think) */
  /* Should be defined in arrayobject.h */

  /* create the module and add the functions */
  m = Py_InitModule("numpyio", numpyio_methods);        /* registration hook */
  
  /* add symbolic constants to the module */
  d = PyModule_GetDict(m);
  ErrorObject = Py_BuildValue("s", "numpyio.error");   /* export exception */
  PyDict_SetItemString(d, "error", ErrorObject);       /* add more if need */

}

/**********************************************************/
/*                                                        */
/*   SYNOPSIS: rbo(data, bpe, nel) ;                      */
/*             where:                                     */
/* 	    nel..... number of array elements             */
/* 	    data.... pointer to the first byte in the     */
/* 	             array                                */
/* 	    bpe..... bytes per array element              */
/*                                                        */
/*   PURPOSE: convert data from little to big endian (and */
/*            visa-versa)                                 */
/*                                                        */
/**********************************************************/

void rbo(char * data, int bpe, int nel) 
{
	int nswaps, i,j;		/* number of swaps to make per element */
	char tmp;			/* temporary storage for swapping      */
	long int p1, p2;		/* indexes for elements to be swapped  */
	
	nswaps = bpe / 2;  		/* divide element size by two          */
	if (nswaps == 0) return;	/* return if it is a byte array        */

	p1 = 0;
	for ( i=0; i<nel; i++) {
		p1 = i*bpe;
		p2 = p1 + bpe - 1;
		for (j=0; j<nswaps; j++) {
			tmp      = data[p1];
			data[p1] = data[p2];
			data[p2] = tmp;
			p1++;
			p2--;
		}
	}
	return;
}

/*  PACKBITS


    This function packs binary (0 or 1) 1-bit per pixel images
        into bytes for writing to disk. 

*/

void packbits(
	      char	In[],
              int       element_size,  /* in bytes */
	      char	Out[],
              int       total_elements,
              int       els_per_slice
	     )
{
  char          build;
  int           i,index,slice,slices,out_bytes;
  int           maxi, remain, nonzero, j;
  char          *outptr,*inptr;

  outptr = Out;                          /* pointer to output buffer */
  inptr  = In;                           /* pointer to input buffer */
  slices = total_elements/els_per_slice;
  out_bytes = ceil( (float) els_per_slice / 8);     /* number of bytes in each slice */
  remain = els_per_slice % 8;                      /* uneven bits */
  if (remain == 0) remain = 8;           /* */
  /*  printf("Start: %d %d %d %d %d\n",inM,MN,slices,out_bytes,remain);
   */
  for (slice = 0; slice < slices; slice++) {
    for (index = 0; index < out_bytes; index++) {
      build = 0;
      maxi = (index != out_bytes - 1 ? 8 : remain);
      for (i = 0; i < maxi ; i++) {
        build <<= 1;                 /* shift bits left one bit */
        nonzero = 0;
        for (j = 0; j < element_size; j++)  /* determine if this number is non-zero */
          nonzero += (*(inptr++) != 0);
        build += (nonzero > 0);                   /* add to this bit if the input value is non-zero */
      }
      if (index == out_bytes - 1) build <<= (8-remain);
      /*      printf("Here: %d %d %d %d\n",build,slice,index,maxi); 
       */
      *(outptr++) = build;
    }
  }
  return;
}


void unpackbits(
		char    In[],
		int     in_element_size,
	        char    Out[],
                int     element_size,
	        int     total_elements,
                int     els_per_slice
               )
{
  unsigned char mask;
  int           i,index,slice,slices,out_bytes;
  int           maxi, remain;
  char          *outptr,*inptr;

  outptr = Out;
  inptr  = In;
  if (is_little_endian()) {
     fprintf(stderr,"This is a little-endian machine.\n"); 
  }
  else {
     fprintf(stderr,"This is a big-endian machine.\n");
     outptr += (element_size - 1);
     inptr  += (in_element_size - 1);
  }
  slices = total_elements / els_per_slice;
  out_bytes = ceil( (float) els_per_slice / 8);
  remain = els_per_slice % 8;
  if (remain == 0) remain = 8;
  /*  printf("Start: %d %d %d %d %d\n",inM,MN,slices,out_bytes,remain);
   */
  for (slice = 0; slice < slices; slice++) {
    for (index = 0; index < out_bytes; index++) {
      maxi = (index != out_bytes - 1 ? 8 : remain);
      mask = 128;
      for (i = 0; i < maxi ; i++) {
        *outptr = ((mask & (unsigned char)(*inptr)) > 0);
        outptr += element_size;
        mask >>= 1;
      }
      /*      printf("Here: %d %d %d %d\n",build,slice,index,maxi); 
       */
      inptr += in_element_size;
    }
  }
  return;
}

int is_little_endian()
{                      /*                             high low  */
  short testnum = 1;   /* If little endian it will be 0x00 0x01 */
                       /* If big endian it will be    0x01 0x00 */
  void *testptr;
  char *myptr;

  testptr = (void *)(&testnum);  /* Assumes address gives low-byte in memory */
  myptr = (char*)testptr;
                              
  return (*(myptr) == 1);
  
}


