/* Copyright (c) 1996, 1997, The Regents of the University of California.
 * All rights reserved.  See Legal.htm for full text and disclaimer. */

/* Gistmodule provides glue between Python and the Gist library of
 * graphics routines.  Quite a bit of the code in this module was adapted
 * from similar code in Dave Munro's Yorick interpreter.  Thanks, Dave.
 */

/* Modified in December, 1997 to coexist with the readline library used
 * in Python's interpreter. William Magro, Cornell Theory Center.
 */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*
 * TODO:
 *   Fri May 24 15:54:16 PDT 1996
 *   There is a significant potential problem with memory leaks in the
 *   current version of the Gist C module. There are many places where
 *   allocated objects don't get cleaned up correctly in the event of an
 *   untimely error. Typical code looks like this:
 * 
 *     PyObject *op = foo (x);
 *     TRY (something_which_might_fail (y));
 *     // Calculate, display, etc.
 *     Py_DECREF (op);
 *     Py_INCREF (Py_None);
 *     return Py_None;
 *   
 *   The problem, of course, is the "return 0" in the usual definition of
 *   the TRY() macro. If its function fails, then the premature return is
 *   taken, "op" never gets DECREF'ed, and its memory leaks. This idiom
 *   shows up all over the place in gistCmodule, and I don't know any
 *   really nice ways to fix it up.
 * 
 *   The best approach I've thought of so far is to redefine TRY() so that
 *   it does a "goto errexit" instead of returning immediately.  The code
 *   fragment above changes to:
 * 
 *     PyObject *op = foo (x);
 *     TRY (something_which_might_fail (y));
 *     // Calculate, display, etc.
 *     Py_DECREF (op);
 *     Py_INCREF (Py_None);
 *     return Py_None;
 *   errexit:
 *     Py_DECREF (op);
 *     return 0;
 *   
 *   The general case is complicated by the possibility of freeing the same
 *   pointer twice. The mechanism of choice for handling this circumstance
 *   would seem to be to adopt a strict convention of zeroing all pointers
 *   at initialization time and later when they are freed, and (of course)
 *   never free a nil pointer. Py_XDECREF() does part of the job (it won't
 *   free a nil pointer), but it needs a little help to finish the job by
 *   zeroing its argument. One other fine point is that not every function
 *   in the Python API sets an exception upon failure. This leads to the
 *   strange-looking error in the interpreter "Zero return without apparent
 *   error" (I'm quoting from memory, sorry.) Therefore the careful errexit
 *   code should check for a legitimate error with PyErr_Occurred(), and
 *   set a default error if there was none.
 * 
 *   All this is straightforward enough, but somewhat tedious and verbose
 *   to apply. And the payback is usually small, because this sort of thing
 *   only happens when there is a previous error to provoke it. And I wish I
 *   could think of a nicer way to do it than with goto and a bunch of even
 *   fancier macros.  So I will wait awhile....
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "Python.h"
#include "Numeric/arrayobject.h"
#include "gist/hlevel.h"
#ifndef NO_XLIB
#  include "gist/dispas.h"
#endif
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #include <dmalloc.h> */

/*#ifndef WANT_SIGFPE_HANDLER*/
#undef PyFPE_START_PROTECT
#undef PyFPE_END_PROTECT
#define PyFPE_START_PROTECT(err_string, leave_stmt)
#define PyFPE_END_PROTECT
/*#endif*/


/* primitive allowance for other non-X windows systems (untested) */
#ifdef DISPATCH_FILE
#  include DISPATCH_FILE
#endif
#ifndef DISPLAY_ENGINE
#  define DISPLAY_ENGINE GpFXEngine
#endif
#ifndef DISPLAY_SET_HANDLER
#  define DISPLAY_SET_HANDLER GhSetXHandler
#endif
#ifdef NO_XLIB
#  define NO_MOUSE
#else
#  ifndef NO_MOUSE
#    ifndef DISPLAY_MOUSE
#      define DISPLAY_MOUSE GxPointClick
#    endif
#  endif
#endif
#ifdef NO_MOUSE
#  ifndef DISPLAY_ZOOM_FACTOR
static double gxZoomFactor = 1.0;
#  endif
#endif
#ifndef DISPLAY_ZOOM_FACTOR
#  define DISPLAY_ZOOM_FACTOR gxZoomFactor
#endif

extern int G_poll(long, unsigned long *, long);

/* We add a component to the default Gist search path, for style and
   palette files in our Python distribution.
*/
static char *gistpath = 0, *oldgistpath = 0;
#define OUR_SPECIAL_DIR "/graphics/gist"

/* Mouse() related stuff */
#ifndef NO_MOUSE
static int MouseCallBack (Engine * engine, int system,
			  int release, GpReal x, GpReal y,
			  int butmod, GpReal xn, GpReal yn);
static double mouseX0, mouseY0, mouseX1, mouseY1;
static double mouseX0ndc, mouseY0ndc, mouseX1ndc, mouseY1ndc;
static int mouseButton, mouseModifier, mouseSystem, mouseError;
static char *defaultPrompts[2] = {
  "<Click mouse at point>", "<Press, drag, and release mouse>" };
#endif

/* The Gist library uses SINGLE_P in order to typedef GpReal, so use the
 * same setting for SINGLE_P as in your copy of libgist.a. (Note that
 * anything other than GpReal == double is untested.)
 */
#ifndef SINGLE_P
#  define Py_GpReal PyArray_DOUBLE
#else
#  define Py_GpReal PyArray_FLOAT
#endif
#define Py_GpColor PyArray_UBYTE

static PyObject *GistError;

#define ERRSS(s) ((PyObject *)(PyErr_SetString(GistError,s),0))
#define SETJMP0 if(setjmp(jmpbuf))return(0)
#define SETKW(ob, target, func, s) \
  if(ob && ob != Py_None && !func(ob, &target, s)) return 0
#define BUILD_KWT(kd,K,kwt) if (-1 == build_kwt (kd, K, kwt)) return 0;
#define NELT(a) (sizeof(a)/sizeof(a[0]))
#define SAFE_FREE(p) (p = (p)?free(p),(void *)0:(void *)0)
/* Old definition of TRY: */
/* #define TRY(e) do{if(!(e))return 0;}while(0) */
/* New definition of TRY: (second argument m = 1 for memerr) */
#define TRY(e, m) do{if(!(e)){clearArrayList();clearFreeList(0);\
               clearMemList();return m;}} while(0)
#define DECL_ZERO(type, var) type var = 0
#define DECREF_AND_ZERO(p) do{Py_XDECREF(p);p=0;}while(0)
#define SETERR(s) if(!PyErr_Occurred()) ERRSS(errstr ? errstr : s)



/* %%%%%%%% */
/* Macros defining most of my uses of a PyArrayObject.  */

/* The number of dimensions of the array. */
#define A_NDIM(a) (((PyArrayObject *)a)->nd)

/* The length of the ith dimension of the array. */
#define A_DIM(a,i) (((PyArrayObject *)a)->dimensions[i])

/* The total number of elements in the array. */
#define A_SIZE(a) PyArray_Size((PyObject *)a)

/* The type number of the array. */
#define A_TYPE(a) (int)(((PyArrayObject *)a)->descr->type_num)

/* The pointer to the array data. */
#define A_DATA(a) (((PyArrayObject *)a)->data)

/* Object is non-null and a PyArrayObject */
#define isARRAY(a) ((a) && ( (PyObject *)a != Py_None) && PyArray_Check((PyArrayObject *)a))

/* Build an array from an object, possibly upcasting the type and copying
   the data to make it contiguous. */
/* Old version: */
/* #define GET_ARR(ap,op,type,dim) \
  TRY(ap=(PyArrayObject *)PyArray_ContiguousFromObject(op,type,dim,dim)) */
/* The new versions below try to add the new object to an appropriate
   list of objects. If they fail, routines which delete all objects
   created so far are called prior to return (see new definition of TRY). */
/* New version: */
#define GET_ARR(ap,op,type,dim,cast) \
  TRY(addToArrayList((PyObject *)(ap=\
  (PyArrayObject *)PyArray_ContiguousFromObject(op,type,dim,dim))), \
  (cast)PyErr_NoMemory ())
/* Build a new array */
#define NEW_ARR(ap,n,dims,type,cast) \
  TRY(addToArrayList((PyObject *)(ap=\
  (PyArrayObject *)PyArray_FromDims(n,dims,type))), \
  (cast)PyErr_NoMemory ())
/* Build an array from existing data */
#define RET_ARR(op,ndim,dim,type,data,cast)\
  TRY(addToArrayList(op=\
  PyArray_FromDimsAndData(ndim,dim,type,data)), \
  (cast)PyErr_NoMemory ())
/* Array owns its data so if DECREF'ed, its data will be freed */
#define SET_OWN(op) ( (PyArrayObject *) op)->flags |= OWN_DATA
/* Array does not own its data so if DECREF'ed, its data will not be freed */
#define UNSET_OWN(op) ( (PyArrayObject *) op)->flags &= ~OWN_DATA
/* Create a block of memory */
#define NEW_MEM(mem,n,type,cast) \
  TRY(addToMemList((void *)(mem=(type *)malloc(n*sizeof(type)))), \
  (cast)PyErr_NoMemory ())

/* %%%%%%%% */

/* Routines needed from outside (but not declared in a header file). */
extern int Py_AtExit (void (*func) (void));
extern int (*PyOS_InputHook)(void);

/********************M E M O R Y   M A N A G E M E N T*****************/
typedef unsigned char Uchar;

/* struct arrayobject is used in slice2 routines */
typedef struct arrayobject {
   void * data ;
   int size ;
   char typecode ;
   } ArrayObject;

#define SAVE -1
#define FREE0 0
#define FREE1 1

/* We have three lists of things to be released upon error return:        */
/* freeList = ArrayObjects (list 0 is slice2 list, list1 is _slice2_part) */
/* PyArrayList = list of NumPy arrays created                             */
/* PyMemList = list of heap pointers created by malloc                    */
/* Some of these items should not be purged upon normal return. The caller*/
/* is responsible for purging only temporaries upon return, and setting   */
/* list lengths to zero.                                                  */

#define MAX_NO_LISTS 2
#define MAX_LIST_SIZE 30

static ArrayObject * freeList [MAX_NO_LISTS] [MAX_LIST_SIZE];
static int freeListLen [MAX_NO_LISTS] = {0, 0};

#define ARRAY_LIST_SIZE 30

static PyObject * PyArrayList [ARRAY_LIST_SIZE];
static int array_list_length = 0;

#define MEM_LIST_SIZE 15

static void * PyMemList [MEM_LIST_SIZE];
static int mem_list_length = 0;

static void clearFreeList (int n);
static void clearArrayList (void);
static void clearMemList (void);

/*****************************************************************/
/*           ArrayObject creation                                */
/* allocateArray creates a complete new ArrayObject including a  */
/* newly malloc'ed data area set to zeros.                       */
/* arrayFromPointer creates a new ArrayObject pointing to a data */
/* area supplied by the client.                                  */
/* Each routine places the new structure on the specified        */
/* FreeList.                                                     */
/* Each returns NULL if it fails for any reason.                 */
/*****************************************************************/

static int addToFreeList (ArrayObject * x, int n);
static void freeArray (ArrayObject * a, int n);
static void removeFromFreeList (ArrayObject * x, int n);
static void dumpFreeList (int n);

static ArrayObject * allocateArray (int size, char tc, int nlist) {
   /* allocate an array object containing size components of type tc */
   ArrayObject * res;
   if (size <= 0)
      return (ArrayObject *) NULL;
   TRY (res = calloc (1, sizeof (ArrayObject)), 
      (ArrayObject *) PyErr_NoMemory ());
   res->size = size ;
   res->typecode = tc ;
   if (size != 0) {
      switch (tc) {
         case 'i' :
            if ( !(res->data = (void *) calloc (size, sizeof (int)))) {
               free (res) ;
               return (ArrayObject *) PyErr_NoMemory ();
               }
            break;
         case 'd' :
            if ( !(res->data = (void *) calloc (size, sizeof (double)))) {
               free (res) ;
               return (ArrayObject *) PyErr_NoMemory ();
               }
            break;
         case 'b' :
            if ( !(res->data = (void *) calloc (size, sizeof (Uchar)))) {
               free (res) ;
               return (ArrayObject *) PyErr_NoMemory ();
               }
            break;
         default :
            free (res) ;
            return (ArrayObject *) NULL ;
         }
      }
   else
      res->data = (void *)NULL;
   if (addToFreeList (res, nlist) != 0) {
      freeArray (res, nlist);
      res = (ArrayObject *) NULL ;
      }
   return res;
   }

static ArrayObject * copyArray (ArrayObject * a) {
   /* returns a copy of the argument array. Does not put it on */
   /* the free list.                                           */
   ArrayObject * res;
   int size;
   if ( a == (ArrayObject *) NULL || a->size <= 0) {
      return (ArrayObject *) NULL;
      }
   switch (a->typecode) {
      case 'b':
         size = sizeof (Uchar);
         break;
      case 'i':
         size = sizeof (int);
         break;
      case 'd':
         size = sizeof (double);
         break;
      default:
         return (ArrayObject *) NULL;
      }
   TRY (res = calloc (1, sizeof (ArrayObject)), 
      (ArrayObject *) PyErr_NoMemory ());
   TRY (res->data = calloc (a->size, size), 
      (ArrayObject *) PyErr_NoMemory ());
   TRY (memcpy (res->data, a->data, a->size * size),
      (ArrayObject *) ERRSS ("memcpy failed in copyArray."));
   res->size = a->size;
   res->typecode = a->typecode;
   return res;
   }

static ArrayObject * arrayFromPointer (int size, char tc, void * data,
   int nlist) {
   /* allocate an array object containing size components of type tc. */
   /* Caller supplies address of data and takes responsibility for    */
   /* its being the right size.                                       */
   ArrayObject * res;
   if (size <= 0)
      return (ArrayObject *) NULL;
   TRY (res = calloc (1, sizeof (ArrayObject)), 
      (ArrayObject *) PyErr_NoMemory ());
   res->size = size ;
   res->typecode = tc ;
   res->data = data ;
   if (addToFreeList (res, nlist) != 0) {
      freeArray (res, nlist);
      res = (ArrayObject *) NULL ;
      }
   return res ;
   }

/*****************************************************************/
/*               getting rid of ArrayObjects                     */
/* freeArray frees the array and its data space, and removes it  */
/* from the FreeList, if it's on it.                             */
/*****************************************************************/

static void freeArray (ArrayObject * a, int n) {
   /* deallocate an array object */
   if (a == (ArrayObject *) NULL)
      return;
   removeFromFreeList (a, n);
   SAFE_FREE (a->data);
   SAFE_FREE (a);
   }

/***********************************************************************
 * FREE LIST MAINTENANCE
 * Routines to create ArrayObjects may return a NULL object, which
 * sometimes is a memory error. Testing for this, and keeping track
 * of created objects in routines like _slice2_part so that they can
 * be removed before error return is a real hassle.
 * I have semi-automated the process with the following routines and
 * conventions.
 * Every time an ArrayObject is successfully created and returned
 * by a routine, it will be added to the freeList. Every time freeArray
 * frees a real object, that object will be removed from the freeList.
 * Prior to return from any routine, a call to clearFreeList will
 * then get rid of all unFreed objects.
 * Note the existence of two freeLists, 0 is for slice2 and 1 is
 * for _slice2_part.
 ************************************************************************/

static void clearFreeList (int n) {
   int i;
   if (n < 0 || n >= MAX_NO_LISTS) return;
   for (i = 0; i < freeListLen [n]; i++) {
      if (freeList [n] [i]) {
         SAFE_FREE (freeList [n] [i]->data);
         }
      SAFE_FREE (freeList [n] [i]);
      }
   freeListLen [n] = 0;
   return;
   }

static void dumpFreeList (int n) {
   /* Useful for debugging ??? */
   int i;
   printf ("-----------start-%d-----------\n", n); fflush (stdout);
   for (i =0; i < freeListLen [n]; i++) {
      printf ("entry %x points to %c data (%d) at %x.\n",
         freeList [n] [i], freeList [n] [i]->typecode, freeList [n] [i]->size,
         freeList [n] [i]->data); fflush (stdout);
      }
   printf ("----------finish-------------\n"); fflush (stdout);
   }

static int addToFreeList (ArrayObject * x, int n) {
   /* A new element is always added at the end. */
   if (n < 0 || n >= MAX_NO_LISTS || freeListLen [n] >= MAX_LIST_SIZE) {
      return -1;
      }
   freeList [n] [ freeListLen [n] ++] = x;
   return 0;
   }

static void removeArrayOnly (ArrayObject * x, int n) {
   /* If the specified list exists, and the item is on it, */
   /* then remove the item, free the item only, and        */
   /* compress the list together.                          */
   int i;
   int found = 0;
   if (n < 0 || n >= MAX_NO_LISTS || x == NULL)
      return;
   for (i =0; i < freeListLen [n]; i++)
      if (!found && freeList [n] [i] == x) {
         SAFE_FREE (freeList [n] [i]);
         found = 1;
         }
      else if (found)
         freeList [n] [i - 1] = freeList [n] [i];
   if (found)
      freeListLen [n] --;
   return;
   }

static void removeFromFreeList (ArrayObject * x, int n) {
   /* If the specified list exists, and the item is on it, */
   /* then remove the item and compress the list together.  */
   int i;
   int found = 0;
   if (n < 0 || n >= MAX_NO_LISTS || x == NULL)
      return;
   for (i =0; i < freeListLen [n]; i++)
      if (!found && freeList [n] [i] == x) {
         found = 1;
         }
      else if (found)
         freeList [n] [i - 1] = freeList [n] [i];
   if (found) {
      freeListLen [n] --;
      }
   return;
   }

/*************************************************************************/
/*              ArrayList maintenance                                    */
/*************************************************************************/

static int addToArrayList (PyObject * obj) {
   if (obj == (PyObject *) NULL || array_list_length > ARRAY_LIST_SIZE)
      return 0;
   PyArrayList [array_list_length ++] = obj;
   return 1;
   }

static void clearArrayList (void) {
   /* DECREF's everything on the ArrayList; needs to be done if there */
   /* is an error.                                                    */
   int i;
   for (i = 0; i < array_list_length; i++)
      Py_DECREF (PyArrayList [i]);
   array_list_length = 0;
   return;
   }

static void removeFromArrayList (PyObject * obj) {
   int i;
   int found = 0;
   if (obj == (PyObject *) NULL)
      return;
   for (i = 0; i < array_list_length; i++)
      if (! found && PyArrayList [i] == obj) {
         Py_DECREF (obj);
         found = 1;
         }
      else if (found)
         PyArrayList [i - 1] = PyArrayList [i];
   if (found)
      array_list_length --;
   }

static void takeOffArrayList (PyObject * obj) {
   /* removes object from list without DECREFing. */
   int i;
   int found = 0;
   for (i = 0; i < array_list_length; i++)
      if (! found && PyArrayList [i] == obj) {
         found = 1;
         }
      else if (found)
         PyArrayList [i - 1] = PyArrayList [i];
   if (found)
      array_list_length --;
   }

/*************************************************************************/
/*              MemList maintenance                                    */
/*************************************************************************/
 
static int addToMemList (void * addr) {
   if (addr == (void *) NULL || mem_list_length >= MEM_LIST_SIZE)
      return 0;
   PyMemList [mem_list_length ++] = addr;
   return 1;
   }
   
static void clearMemList (void) {
   int i;
   for (i = 0; i < mem_list_length; i++) {
      SAFE_FREE (PyMemList [i]);
      }
   mem_list_length = 0;
   return;
   }

/******************************E N D***********************************/

/* This module's initialization routine */
void initgistC (void);

/* Routines in the Gist user interface */
static PyObject *animate (PyObject * self, PyObject * args);
static PyObject *bytscl (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *contour (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *current_window (PyObject * self, PyObject * args);
static PyObject *debug_array (PyObject * self, PyObject * args);
static PyObject *g_fma (PyObject * self, PyObject * args);
static PyObject *get_slice2_precision (PyObject * self, PyObject * args);
static PyObject *gridxy (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *hcp (PyObject * self, PyObject * args);
static PyObject *hcp_file (PyObject * self, PyObject * args, PyObject *kd);
static PyObject *hcp_finish (PyObject * self, PyObject * args);
static PyObject *hcpoff (PyObject * self, PyObject * args);
static PyObject *hcpon (PyObject * self, PyObject * args);
static PyObject *limits (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *logxy (PyObject * self, PyObject * args);
static PyObject *mesh_loc (PyObject * self, PyObject * args);
static PyObject *mfit (PyObject * self, PyObject * args);
static PyObject *mouse (PyObject * self, PyObject * args);
static PyObject *palette (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *pyg_pause (PyObject * self, PyObject * args);
static PyObject *plc (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *pldefault (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *pldj (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *pledit (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *plf (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *plfp (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *plg (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *pli (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *plm (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *plmesh (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *plq (PyObject * self, PyObject * args);
static PyObject *plsys (PyObject * self, PyObject * args);
static PyObject *plt (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *plv (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *redraw (PyObject * self, PyObject * args);
static PyObject *set_slice2_precision (PyObject * self, PyObject * args);
static PyObject *setdpi (PyObject * self, PyObject * args);
static PyObject *slice2 (PyObject * self, PyObject * args);
static PyObject *unzoom (PyObject * self, PyObject * args);
static PyObject *viewport (PyObject * self, PyObject * args);
static PyObject *window (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *zoom_factor (PyObject * self, PyObject * args);

/* Utility routines */
static GpColor *PushColors(double *z, long len, double zmin,
  double zmax, double scale, double offset);
static char *GetHCPname (int n);
static char *SetHCPname (int n, char *name);
static char *expand_pathname (const char *name);
static double Safe_dbl (double x);
static int CheckDefaultWindow (void);
static int CheckPalette (void);
static int GetTypeface (char *s, int *f);
static int GrabByteScale ( PyObject **kwt, char **keywrds, double *scale,
  double *offset, double *zn, double *zx, double *z, int *reg, 
    int region, long iMax, long jMax, int zCompressed);
static int KeyboardInput (FILE * file, void *context);
static int SetHCPDefault (void);
static int YPrompt(const char *s);
static int build_kwt (PyObject *kd, char *Keys[], PyObject * kwt[]);
static int set_def_reg (int nr, int nc);
static int set_limit (PyObject * ob, double *lim, int *flags, int fval);
static int set_pyMsh(PyObject *args, char *errstr, PyObject *tri);
static int set_reg (PyObject *rop);
static int set_tri (PyObject *top);
static int set_yx (PyObject *yop, PyObject *xop);
static int setkw_boolean (PyObject * v, int *t, char *kw);
static int setkw_color (PyObject * v, int *t, char *kw);
static int setkw_double (PyObject * v, double *t, char *kw);
static int setkw_fonttype (PyObject * v, int *t, char *kw);
static int setkw_integer (PyObject * v, int *t, char *kw);
static int setkw_justify (PyObject * v, int *t, char *kw);
static int setkw_linetype (PyObject * v, int *t, char *kw);
static int setkw_string (PyObject * v, char **t, char *kw);
static int setkw_xinteger (PyObject * v, int *t, char *kw);
static int setvu_mesh( PyObject *args, PyObject **vop, 
  PyObject **uop, char *errstr);
static int setz_mesh( PyObject *args, PyObject **zop, 
  char *errstr, PyObject *tri);
static int unpack_limit_tuple (PyObject * ob, double limits[], int *flags);
static int verify_kw (char *keyword, char * kwlist[]);
static long FindMeshZone(double xx, double yy, double *x, double *y, 
  int *reg, long ix, long jx);
static void CheckDefaultPalette (void);
static void CleanUpGraphics (void);
static void ForceNewline (void);
static void GetPCrange (double *zmn, double *zmx, double *z, 
  int *reg, int region, long iMax, long jMax);
static void GetZCrange(double *zmn, double *zmx, double *z, 
  int *reg, int region, long iMax, long jMax, int zCompressed);
static void PermitNewline (int nSpaces);
static void PrintColor (char *line, int color, int suffix);
static void PrintFunc (const char *s);
static void PrintHideLegend (char *line, int type);
static void PrintInit (int (*puts) (const char *));
static void PrintMarks (char *line, int suffix);
static void PrintRegion (char *line, int suffix);
static void PrintSuffix (int suffix);
static void PrintTypeWidth (char *line, int suffix);
static void Xerror_longjmp (char *message);
static void YGDispatch (void);
static int YGEventHandler(void);
static void clear_pyMsh(void);
static void get_mesh(GaQuadMesh *m);

/* Global variables */
/*
 * All the mesh routines share this pyMsh struct. (Why? Because the mesh
 * arrays are saved between calls to mesh routines, and in general, mesh
 * routines allow the mesh array arguments to be optional, defaulting to
 * the currently saved mesh.) PyMsh mirrors the Gist GaQuadMesh structure.
 * Gist documentation/code uses the variables iMax and jMax to designate
 * the number of elements in each dimension of the two 2D arrays X and Y
 * whose respective elements give the coordinates of the quadrilateral mesh
 * nodes. I sometimes use iMax, jMax in this file, but more often use "nr",
 * "nc" to designate the number of rows and columns in the arrays. The
 * following correspondence generally is true:
 *
 *      iMax == nc == A_DIM(array, 1)
 *      jMax == nr == A_DIM(array, 0)
 */
static struct {
  PyArrayObject *y, *x, *reg, *triangle;
} pyMsh = { 0, 0, 0, 0 };

static int yPendingIn = 0;
static jmp_buf jmpbuf;
static int already_initialized = 0;
static int curPlotter = -1;
static int curElement = -1;
static int paletteSize = 0;
static int maxColors = 200; /* maximum number of colors for GpReadPalette */
static int hcpDump = 1;
static int hcpPSdefault = 0;
static int hcpOnFMA = 0;
static int defaultDPI = 100;
static int defaultLegends = 1;
static char *defaultStyle = 0;
static char *defaultPalette = 0;
static char *hcpNames[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static char *windowNames[8] = {
  "Pygist 0", "Pygist 1", "Pygist 2", "Pygist 3",
  "Pygist 4", "Pygist 5", "Pygist 6", "Pygist 7"
};

/* Next few variables are used by plq(), which does some fancy printing. */
/*
static long prop3sizes[10] = {0, 8, 2, 5, 5, 3, 3, 7, 0, 0};
static long prop4sizes[10] = {0, 8, 1, 3, 1, 1, 3, 4, 4, 0};
static long prop5sizes[10] = {0, 3, 5, 2, 5, 6, 7, 9, 3, 5};
*/
static int curIX = -1, curIXc = -1;
static char specialMarkers[5] = ".+*ox";
static int (*RawPrinter) (const char *s);
static int printLength = 79;	/* maximum number of characters on a line */
  /*
static int lenPrintBuf = 79;
  */
static char printBuf[80];
static long maxPrintLines = 5000;
static int printNow, permitNow;
static long printLines;
static double _slice2_precision = 0.0;

static char gist_module_documentation[] =
"Gist Graphics Package, version1.3"
;

#define PYCF (PyCFunction)
/*#define PYCFWK (PyCFunctionWithKeywords)*/
#define PYCFWK (PyCFunction) /* Make the compiler shut up. */
#define KWFLG (METH_VARARGS | METH_KEYWORDS)

static struct PyMethodDef gist_methods[] =
{ 
  { "animate",        PYCF   animate,        1,     0 },
  { "bytscl",         PYCFWK bytscl,         KWFLG, 0 },
  { "contour",        PYCFWK contour,        KWFLG, 0 },
  { "current_window", PYCF   current_window, 1,     0 },
  { "debug_array",    PYCF   debug_array,    1,     0 },
  { "fma",            PYCF   g_fma,          1,     0 },
  { "gridxy",         PYCFWK gridxy,         KWFLG, 0 },
  { "get_slice2_precision", PYCF get_slice2_precision, 1, 0},
  { "hcp",            PYCF   hcp,            1,     0 },
  { "hcp_file",       PYCFWK hcp_file,       KWFLG, 0 },
  { "hcp_finish",     PYCF   hcp_finish,     1,     0 },
  { "hcpoff",         PYCF   hcpoff,         1,     0 },
  { "hcpon",          PYCF   hcpon,          1,     0 },
  { "limits",         PYCFWK limits,         KWFLG, 0 },
  { "logxy",          PYCF   logxy,          1,     0 },
  { "mesh_loc",       PYCF   mesh_loc,       1,     0 },
  { "mfit",           PYCF   mfit,           1,     0 },
  { "mouse",          PYCF   mouse,          1,     0 },
  { "palette",        PYCFWK palette,        KWFLG, 0 },
  { "pause",          PYCF   pyg_pause,      1,     0 },
  { "plc",            PYCFWK plc,            KWFLG, 0 },
  { "pldefault",      PYCFWK pldefault,      KWFLG, 0 },
  { "pldj",           PYCFWK pldj,           KWFLG, 0 },
  { "pledit",         PYCFWK pledit,         KWFLG, 0 },
  { "plf",            PYCFWK plf,            KWFLG, 0 },
  { "plfp",           PYCFWK plfp,           KWFLG, 0 },
  { "plg",            PYCFWK plg,            KWFLG, 0 },
  { "pli",            PYCFWK pli,            KWFLG, 0 },
  { "plm",            PYCFWK plm,            KWFLG, 0 },
  { "plmesh",         PYCFWK plmesh,         KWFLG, 0 },
  { "plq",            PYCF   plq,            1,     0 },
  { "plsys",          PYCF   plsys,          1,     0 },
  { "plt",            PYCFWK plt,            KWFLG, 0 },
  { "plv",            PYCFWK plv,            KWFLG, 0 },
  { "redraw",         PYCF   redraw,         1,     0 },
  { "set_slice2_precision", PYCF set_slice2_precision, 1, 0},
  { "slice2",         PYCF   slice2,         1,     0 },
  { "set_default_dpi", PYCF setdpi,          1,     0 },
  { "unzoom",         PYCF   unzoom,         1,     0 },
  { "viewport",       PYCF   viewport,       1,     0 },
  { "window",         PYCFWK window,         KWFLG, 0 },
  { "zoom_factor",    PYCF   zoom_factor,    1,     0 },

  { 0, 0 }
};

/*#######################################################################*/
/*           Auxiliary routines used by the slice2 suite                 */
/*#######################################################################*/

static ArrayObject * Add1 ( ArrayObject * i, int freei, int n) {
   /* add 1 to an integer array i and return the value */
   ArrayObject * res ;
   int * isrc ,
       * itar;
   if (i == (ArrayObject *) NULL)
      return (ArrayObject *) NULL;
   if (freei >= 0 && freei == n) {
      /* We don't need i any more and it's in the right list, so *
       * just increment it.                                      */
      for (isrc = (int *) (i->data);
           isrc < (int *) (i->data) + i->size; isrc ++)
         (* isrc) ++;
      return i;
      }
   else {
      TRY (res = allocateArray (i->size, 'i', n),
         (ArrayObject *) PyErr_NoMemory ());
      for (isrc = (int *) (i->data), itar = (int *) (res->data);
           isrc < (int *) (i->data) + i->size; isrc ++, itar ++)
         * itar = (* isrc) + 1;
      if (freei >= 0)
         freeArray (i, freei);
      return res;
      }
   }
 
static ArrayObject * concatenate (ArrayObject * x, ArrayObject * y,
   int freex, int freey, int n) {
   /* concatenates array objects x and y and returns the result */
   /* returns NULL if memory error or if both objects are NULL. */
   /* Note that the new array is put onto the free list.        */
   ArrayObject * result;
   int datasize;
   Uchar tc;
   if (x == (ArrayObject *) NULL && y == (ArrayObject *) NULL) {
      return (ArrayObject *) NULL;
      }
   tc = (y == (ArrayObject *) NULL) ? x->typecode : y->typecode;
   switch (tc) {
      case 'b':
         datasize = sizeof (Uchar);
         break;
      case 'i':
         datasize = sizeof (int);
         break;
      case 'd':
         datasize = sizeof (double);
         break;
      default:
         return (ArrayObject *) NULL;
      }
   if (y == (ArrayObject *) NULL) {
      TRY (result = allocateArray (x->size, x->typecode, n), 
         (ArrayObject *) PyErr_NoMemory ());
      TRY (memcpy (result->data, x->data, x->size * datasize),
         (ArrayObject *) ERRSS ("memcpy failed in concatenate."));
      if (freex >= 0)
         freeArray (x, freex);
      return result;
      }
   if (x == (ArrayObject *) NULL) {
      TRY (result = allocateArray (y->size, y->typecode, n), 
         (ArrayObject *) PyErr_NoMemory ());
      TRY (memcpy (result->data, y->data, y->size * datasize),
         (ArrayObject *) ERRSS ("memcpy failed in concatenate."));
      if (freey >= 0)
         freeArray (y, freey);
      return result;
      }
   if (x->typecode != y->typecode)
      return (ArrayObject *) NULL;
   TRY (result = allocateArray (x->size + y->size, tc, n), 
      (ArrayObject *) PyErr_NoMemory ());
   TRY (memcpy (result->data, x->data, x->size * datasize),
      (ArrayObject *) ERRSS ("memcpy failed in concatenate."));
   TRY (memcpy ( (void *) ( (Uchar *) (result->data) + x->size * datasize),
      y->data, y->size * datasize),
      (ArrayObject *) ERRSS ("memcpy failed in concatenate."));
   if (freex >= 0)
      freeArray (x, freex);
   if (freey >= 0)
      freeArray (y, freey);
   return result;
   }

static ArrayObject * cumsum (ArrayObject * i, int freei, int n) {
   /* compute the cumulative sum of an integer array object */
   ArrayObject * res ;
   int * src,
       * tar,
       pre;
   if (i == (ArrayObject *) NULL)
      return (ArrayObject *) NULL;
   src = (int *) (i->data);
   TRY (res = allocateArray (i->size, 'i', n), 
      (ArrayObject *) PyErr_NoMemory ());
   for (src = (int *) (i->data), tar = (int *) (res->data), pre = 0;
        src < (int *) (i->data) + i->size; src ++, tar++) {
      * tar += * src + pre;
      pre = * tar;
      }
   if (freei >= 0)
      freeArray (i, freei);
   return res;
   }

static ArrayObject * equal (ArrayObject * i, ArrayObject * j,
   int freei, int freej, int n) {
   /* Compare two integer arrays for equality, return an unsigned *
    * character array of results.                                 */
   int * src1,
       * src2;
   ArrayObject * res;
   Uchar * tar;
   if (i == (ArrayObject *)NULL || j == (ArrayObject *)NULL ||
      i->size != j->size)
      return (ArrayObject *) NULL;
   TRY (res = allocateArray (i->size, 'b', n), 
      (ArrayObject *) PyErr_NoMemory ());
   for (src1 = (int *) (i->data), src2 = (int *) (j->data),
        tar = (Uchar *) (res->data); src1 < (int *) (i->data) + i->size;
        src1 ++, src2 ++, tar ++)
      * tar = (Uchar) (* src1 == * src2);
   if (freej >= 0)
      freeArray (j, freej);
   if (freei >= 0)
      freeArray (i, freei);
   return res ;
   }

static ArrayObject * greater (ArrayObject * d, double v, int freed, int n) {
   /* compare an array object of type 'd' with a double v; return an
    * array  of Uchar the same size with 1's where the compare succeeds. */
   double * dr ;
   Uchar * c ;
   ArrayObject * res ;
   if (d == (ArrayObject *) NULL)
      return (ArrayObject *) NULL;
   TRY (res = allocateArray (d->size, 'b', n), 
      (ArrayObject *) PyErr_NoMemory ());
   for (dr = (double *) (d->data), c = (Uchar *) (res->data);
        dr < (double *) (d->data) + d->size; dr ++, c++)
      * c = (Uchar) (* dr > v);
   if (freed >= 0)
      freeArray (d, freed);
   return res ;
   }

static ArrayObject * greater_equal (ArrayObject * d, double v,
   int freed, int n) {
   /* compare an array object of type 'd' with a double v; return an
    * array  of Uchar the same size with 1's where the compare succeeds. */
   double * dr ;
   Uchar * c ;
   ArrayObject * res ;
   if (d == (ArrayObject *) NULL)
      return (ArrayObject *) NULL;
   TRY (res = allocateArray (d->size, 'b', n), 
      (ArrayObject *) PyErr_NoMemory ());
   for (dr = (double *) (d->data), c = (Uchar *) (res->data);
        dr < (double *) (d->data) + d->size; dr ++, c++)
      * c = (Uchar) (* dr >= v);
   if (freed >= 0)
      freeArray (d, freed);
   return res ;
   }

static ArrayObject * less (ArrayObject * i, ArrayObject * j,
   int freei, int freej, int n) {
   int * src1,
       * src2;
   ArrayObject * res;
   Uchar * tar;
   if (i == (ArrayObject *) NULL || j == (ArrayObject *) NULL ||
      i->size != j->size)
      return (ArrayObject *) NULL;
   TRY (res = allocateArray (i->size, 'b', n), 
      (ArrayObject *) PyErr_NoMemory ());
   for (src1 = (int *) (i->data), src2 = (int *) (j->data),
        tar = (Uchar *) (res->data); src1 < (int *) (i->data) + i->size;
        src1 ++, src2 ++, tar ++)
      * tar = (Uchar) (* src1 < * src2);
   if (freej >= 0)
      freeArray (j, freej);
   if (freei >= 0)
      freeArray (i, freei);
   return res ;
   }

static ArrayObject * logical_and (ArrayObject * a, ArrayObject * b,
   int freea, int freeb, int n) {
   /* Take the logical and of two Uchar arrays */
   ArrayObject * res;
   Uchar * src1,
         * src2,
         * tar;
   if (a == (ArrayObject *) NULL || b == (ArrayObject *) NULL ||
      a->size != b->size || a->typecode != 'b' || b->typecode != 'b')
      return (ArrayObject *) NULL;
   src1 = (Uchar *) (a->data);
   src2 = (Uchar *) (b->data);
   if (freea == n) /* can use a as target */ {
      tar = (Uchar *) (a->data);
      freea = -1;
      }
   else if (freeb == n) /* can use b as target */ {
      tar = (Uchar *) (b->data);
      freeb = -1;
      }
   else { /* need new result as target */
      TRY (res = allocateArray (a->size, 'b', n),
         (ArrayObject *) PyErr_NoMemory ());
      tar = (Uchar *) (res->data);
      }
   for (; tar < (Uchar *) (res->data) + res->size; src1 ++, src2 ++, tar ++)
      * tar = * src1 && * src2;
   if (freea >= 0)
      freeArray (a, freea);
   if (freeb >= 0)
      freeArray (b, freea);
   return (res);
   }

static ArrayObject * logical_not (ArrayObject * b, int freeb, int n) {
   /* Take the logical not of an int or Uchar array and return Uchar */
   ArrayObject * res;
   Uchar * usrc,
         * tar;
   int * isrc;
   if (b == (ArrayObject *) NULL)
      return (ArrayObject *) NULL;
   if (freeb >= 0 && freeb == n) {
      /* do it in place */
      if (b->typecode == 'b') {
         for (usrc = (Uchar *) (b->data);
              usrc < (Uchar *) (b->data) + b->size; usrc ++)
            * usrc = (Uchar) ! * usrc;
         return b;
         }
      else if (b->typecode == 'i') {
         for (tar = (Uchar *) (b->data), isrc = (int *) (b->data);
              tar < (Uchar *) (b->data) + b->size; isrc ++, tar ++)
            * tar = (Uchar) ! * isrc;
         b->typecode = 'b';
         return b;
         }
      else
         return (ArrayObject *) NULL;
      }
   else {
      TRY (res = allocateArray (b->size, 'b', n), 
         (ArrayObject *) PyErr_NoMemory ());
      if (b->typecode == 'b')
         for (usrc = (Uchar *) (b->data), tar = (Uchar *) (res->data);
              usrc < (Uchar *) (b->data) + b->size; usrc ++, tar++)
            * tar = (Uchar) ! * usrc;
      else if (b->typecode == 'i')
         for (isrc = (int *) (b->data), tar = (Uchar *) (res->data);
              isrc < (int *) (b->data) + b->size; isrc ++, tar++)
            * tar = (Uchar) ! * isrc;
      else {
         freeArray (res, n);
         return (ArrayObject *) NULL;
         }
      if (freeb >= 0)
         freeArray (b, freeb);
      return res ;
      }
   }

static ArrayObject * not_equal (ArrayObject * i, int j, int freei, int n) {
   /* compare integer array against integer for not equal */
   int * src;
   ArrayObject * res;
   Uchar * tar;
   if (i == (ArrayObject *) NULL)
      return (ArrayObject *) NULL;
   TRY (res = allocateArray (i->size, 'b', n), 
      (ArrayObject *) PyErr_NoMemory ());
   for (src = (int *) (i->data), tar = (Uchar *) (res->data);
        src < (int *) (i->data) + i->size;
        src ++, tar ++)
      * tar = (Uchar) (* src != j);
   if (freei >= 0)
      freeArray (i, freei);
   return res ;
   }

static ArrayObject * SimpleHist (ArrayObject * i, int freei, int n) {
   /* returns an array which counts the number of occurrences of
    * each element of i, which must consist of non-negative
    * integers.                                             */
   ArrayObject * res ;
   int * src,
       * tar;
   int max;
   if (i == (ArrayObject *) NULL) {
      return (ArrayObject *) NULL;
      }
   for (src = (int *) (i->data), max = * (src ++);
        src < (int *) (i->data) + i->size; src++)
      if (* src > max)
         max = * src;
      else if (* src < 0) {
         return (ArrayObject *) NULL;
         }
   TRY (res = allocateArray (max + 1, 'i', n), 
      (ArrayObject *) PyErr_NoMemory ());
   for (src = (int *) (i->data), tar = (int *) (res->data);
        src < (int *) (i->data) + i->size; src++)
      tar [* src] += 1;
   if (freei >= 0)
      freeArray (i, freei);
   return res;
   }

static ArrayObject * take (ArrayObject * a, ArrayObject * i, int nelts,
   int freea, int freei, int n) {
   /* Return an array the same size as i, which consists of the
    * elements of a subscripted by the elements of i. if nelts is
    * not equal to 1, this acts as if a is an array of things
    * nelts long, and does the extraction accordingly. Currently
    * this is only implemented for type 'd'.
    * returns NULL in the event of any error.                      */
   ArrayObject * res;
   double * dtar,
          * dsrc;
   int * itar,
       * isrc;
   int * ind;
   Uchar * utar,
         * usrc;
   int j;
   if (a == (ArrayObject *) NULL || nelts <= 0 ||
       i == (ArrayObject *) NULL || i->size == 0) {
      return (ArrayObject *) NULL ;
      }
   TRY (res = allocateArray (nelts * i->size, a->typecode, n), 
      (ArrayObject *) PyErr_NoMemory ());
   ind = (int *) (i->data);
   switch (a->typecode) {
      case 'b' :
         for (usrc = (Uchar *) (a->data), utar = (Uchar *) (res->data);
              ind < (int *) (i->data) + i->size; utar ++, ind ++)
            if (* ind >= 0 && * ind < a->size)
               * utar = usrc [* ind];
            else {
               freeArray (res, n);
               return (ArrayObject *) NULL;
               }
         break ;
      case 'i' :
         for (isrc = (int *) (a->data), itar = (int *) (res->data);
              ind < (int *) (i->data) + i->size; itar ++, ind ++)
            if (* ind >= 0 && * ind < a->size)
               * itar = isrc [* ind];
            else {
               freeArray (res, n);
               return (ArrayObject *) NULL;
               }
         break ;
      case 'd' :
         if (nelts == 1)
            for (dsrc = (double *) (a->data), dtar = (double *) (res->data);
                 ind < (int *) (i->data) + i->size; dtar ++, ind ++)
               if (* ind >= 0 && * ind < a->size)
                  * dtar = dsrc [* ind];
               else {
                  freeArray (res, n);
                  return (ArrayObject *) NULL;
                  }
         else if (nelts > 0)
            for (dsrc = (double *) (a->data), dtar = (double *) (res->data);
                 ind < (int *) (i->data) + i->size; dtar += nelts, ind ++)
               if (* ind >= 0 && * ind * nelts < a->size)
                  for (j = 0; j < nelts; j++)
                     dtar [j] = dsrc [* ind * nelts + j];
               else {
                  freeArray (res, n);
                  return (ArrayObject *) NULL;
                  }
         else {
            freeArray (res, n);
            return (ArrayObject *) NULL;
            }
         break ;
      default :
         freeArray (res, n);
         return (ArrayObject *) NULL ;
      }
   if (freea >= 0)
      freeArray (a, freea);
   if (freei >= 0)
      freeArray (i, freea);
   return res ;
   }

static ArrayObject * WeightedHist (ArrayObject *i, ArrayObject *w,
   int freei, int freew, int n) {
   /* weighted histogram: the kth entry of the result contains the     *
    * sum of all entries in w corresponding to entries in i containing *
    * the value k. (SimpleHist consists of all weights being 1.)       */
   ArrayObject * res ;
   int * src,
       * tar;
   Uchar * wgt;
   int max;
   if (i == (ArrayObject *) NULL || w == (ArrayObject *) NULL ||
      i->size > w->size)
      return (ArrayObject *) NULL;
   for (src = (int *) (i->data), max = * (src ++);
        src < (int *) (i->data) + i->size; src++)
      if (* src > max)
         max = * src;
      else if (* src < 0)
         return (ArrayObject *) NULL;
   TRY (res = allocateArray (max + 1, 'i', n), 
      (ArrayObject *) PyErr_NoMemory ());
   for (src = (int *) (i->data), tar = (int *) (res->data),
        wgt = (Uchar *) (w->data);
        src < (int *) (i->data) + i->size; src++, wgt++)
      tar [* src] += (int) (* wgt);
   if (freew >= 0)
      freeArray (w, freew);
   if (freei >= 0)
      freeArray (i, freei);
   return res;
   }

/*           _slice2_part is a big helper routine for slice2              */

/************************************************************************/

static int _slice2_part (ArrayObject * xyzc, ArrayObject * keep,
   ArrayObject * next, ArrayObject * dp, ArrayObject * prev,
   ArrayObject * last, ArrayObject * valc, ArrayObject ** xyzc_new,
   ArrayObject ** nvertc, ArrayObject ** valc_new, int freexyzc, 
   int freevalc, char atype) {
   /* auxiliary routine for slice2. returns 0 if there was any kind
    * of error, 1 otherwise.                                        */
   int nlist0, nlist1, i, j, i0, i1, listtmp;
   double dpl, dpu, dpu_dpl;
   ArrayObject * list = (ArrayObject *) NULL,
               * list0 = (ArrayObject *) NULL,
               * list1 = (ArrayObject *) NULL,
               * mask = (ArrayObject *) NULL,
               * ndxs = (ArrayObject *) NULL,
               * valc0 = (ArrayObject *) NULL,
               * valc1 = (ArrayObject *) NULL,
               * xold = (ArrayObject *) NULL,
               * xyz0 = (ArrayObject *) NULL,
               * xyz1 = (ArrayObject *) NULL;
   Uchar * keepd,
         * maskd;
   int * nextd,
       * nvertcd,
       * prevd,
       * list0d = NULL,
       * list1d = NULL,
       * listd,
       * xoldd;
   /*#######*/int * ndxsd;
   double * xyz0d = NULL,
          * valcd = (double *) NULL,
          * valc0d = (double *) NULL,
          * valc1d = (double *) NULL,
          * dpd,
          * xyzcd,
          * xyz1d = NULL,
          * xyzc_newd,
          * valc_newd = (double *) NULL;
   Uchar * valcc = (Uchar *) NULL,
         * valc0c = (Uchar *) NULL,
         * valc1c = (Uchar *) NULL,
         * valc_newc = (Uchar *) NULL;

   /* i = dmalloc_verify (0); */
   * xyzc_new = (ArrayObject *) NULL;
   * nvertc = (ArrayObject *) NULL;
   if (valc) {
      /* Don't touch the following if no values are coming in!! */
      * valc_new = (ArrayObject *) NULL;
      if (atype == 'd')
         valcd = (double *) (valc->data);
      else if (atype == 'b')
         valcc = (Uchar *) (valc->data);
      }
   dpd = (double *) (dp->data);
   xyzcd = (double *) (xyzc->data);
   TRY (mask = allocateArray (keep->size, 'b', 1), (int) NULL);
   maskd = (Uchar *) (mask->data);
   keepd = (Uchar *) (keep->data);
   nextd = (int *) (next->data);
   prevd = (int *) (prev->data);
   for (nlist0 = 0, nlist1 = 0, i = 0; i < keep->size; i ++) {
      i0 = ! keepd [i] && keepd [nextd [i]];
      i1 = ! keepd [i] && keepd [prevd [i]];
      nlist0 += i0;
      nlist1 += i1;
      maskd [i] = keepd [i] + i0 + i1;
      }
   if (nlist0 != 0) {
      TRY (list0 = allocateArray (nlist0, 'i', 1), (int) NULL);
      list0d = (int *) (list0->data);
      for (j = 0, i = 0; i < keep->size; i++)
         if ( ! keepd [i] && keepd [nextd [i]])
            list0d [j ++] = i;
      TRY (xyz0 = allocateArray (3 * nlist0, 'd', 1), (int) NULL);
      xyz0d = (double *) (xyz0->data);
      if (valc != (ArrayObject *)NULL) {
         TRY (valc0 = allocateArray (nlist0, atype, 1), (int) NULL);
         if (atype == 'd')
            valc0d = (double *) (valc0->data);
         else if (atype == 'b')
            valc0c = (Uchar *) (valc0->data);
         }
      for (i = 0; i < nlist0; i++) {
         dpl = dpd [list0d [i]];
         listtmp = nextd [list0d [i]];
         dpu = dpd [listtmp];
         dpu_dpl = dpu - dpl;
         xyz0d [3 * i] = (xyzcd [3 * list0d [i]] * dpu -
            xyzcd [3 * listtmp] * dpl) / dpu_dpl;
         xyz0d [3 * i + 1] = (xyzcd [3 * list0d [i] + 1] * dpu -
            xyzcd [3 * listtmp + 1] * dpl) / dpu_dpl;
         xyz0d [3 * i + 2] = (xyzcd [3 * list0d [i] + 2] * dpu -
            xyzcd [3 * listtmp + 2] * dpl) / dpu_dpl;
         if (valc != (ArrayObject *)NULL) {
            if (atype == 'd')
               valc0d [i] = (valcd [list0d [i]] * dpu -
                  valcd [listtmp] * dpl) / dpu_dpl;
            else if (atype == 'b')
               valc0c [i] = (valcc [list0d [i]] * dpu -
                  valcc [listtmp] * dpl) / dpu_dpl;
            }
         }
      }
   if (nlist1 != 0) {
      TRY (list1 = allocateArray (nlist1, 'i', 1), (int) NULL);
      list1d = (int *) (list1->data);
      for (j = 0, i = 0; i < keep->size; i++)
         if (! keepd [i] && keepd [prevd [i]])
            list1d [j ++] = i;
      TRY (xyz1 = allocateArray (3 * nlist1, 'd', 1), (int) NULL);
      xyz1d = (double *) (xyz1->data);
      if (valc != (ArrayObject *)NULL) {
         TRY (valc1 = allocateArray (nlist1, atype, 1), (int) NULL);
         if (atype == 'd')
            valc1d = (double *) (valc1->data);
         else if (atype == 'b')
            valc1c = (Uchar *) (valc1->data);
         }
      for (i = 0; i < nlist1; i++) {
         dpl = dpd [list1d [i]];
         listtmp = prevd [list1d [i]];
         dpu = dpd [listtmp];
         dpu_dpl = dpu - dpl;
         xyz1d [3 * i] = (xyzcd [3 * list1d [i]] * dpu -
            xyzcd [3 * listtmp] * dpl) / dpu_dpl;
         xyz1d [3 * i + 1] = (xyzcd [3 * list1d [i] + 1] * dpu -
            xyzcd [3 * listtmp + 1] * dpl) / dpu_dpl;
         xyz1d [3 * i + 2] = (xyzcd [3 * list1d [i] + 2] * dpu -
            xyzcd [3 * listtmp + 2] * dpl) / dpu_dpl;
         if (valc != (ArrayObject *)NULL) {
            if (atype == 'd')
               valc1d [i] = (valcd [list1d [i]]  * dpu -
                  valcd [listtmp] * dpl) / dpu_dpl;
            else if (atype == 'b')
               valc1c [i] = (valcc [list1d [i]]  * dpu -
                  valcc [listtmp] * dpl) / dpu_dpl;
            }
         }
      }
   TRY (list = allocateArray (mask->size, 'i', 1), (int) PyErr_NoMemory ());
   listd = (int *) (list->data);
   for (i = 0, i0 = 0; i < mask->size; i ++) {
      i0 += maskd [i];
      listd [i] = i0;
      }
   TRY (xold = allocateArray (listd [list->size - 1], 'i', 1), (int) NULL);
   xoldd = (int *) (xold->data);
   for (i = 0; i < mask->size; i ++) {
      if ( maskd [i] != 0) 
         xoldd [listd [i] - 1] = i;
      if ( maskd [i] == 2)
         xoldd [listd [i] - 2] = i;
      }
   /* Note: allocate return values in caller's list. */
   TRY (* xyzc_new = allocateArray (3 * xold->size, 'd', 0), (int) NULL);
   xyzc_newd = (double *) ( (* xyzc_new)->data);
   if (valc != (ArrayObject *)NULL) {
      TRY (* valc_new = allocateArray (xold->size, atype, 0), (int) NULL);
      if (atype == 'd')
         valc_newd = (double *) ( (* valc_new)->data);
      else if (atype == 'b')
         valc_newc = (Uchar *) ( (* valc_new)->data);
      }
   for (i = 0; i < xold->size; i++) {
      xyzc_newd [3 * i] = xyzcd [3 * xoldd [i]];
      xyzc_newd [3 * i + 1] = xyzcd [3 * xoldd [i] + 1];
      xyzc_newd [3 * i + 2] = xyzcd [3 * xoldd [i] + 2];
      if (valc != (ArrayObject *)NULL) {
         if (atype == 'd')
            valc_newd [i] = valcd [xoldd [i]];
         else if (atype == 'b')
            valc_newc [i] = valcc [xoldd [i]];
      }
      }
   for (i = 0; i < nlist0; i ++) {
      xyzc_newd [3 * listd [list0d [i]] - 3] = xyz0d [3 * i];
      xyzc_newd [3 * listd [list0d [i]] - 2] = xyz0d [3 * i + 1];
      xyzc_newd [3 * listd [list0d [i]] - 1] = xyz0d [3 * i + 2];
      if (valc != (ArrayObject *)NULL) {
         if (atype == 'd')
            valc_newd [listd [list0d [i]] - 1] = valc0d [i];
         else if (atype == 'b')
            valc_newc [listd [list0d [i]] - 1] = valc0c [i];
      }
      }
   freeArray (list0, 1);
   freeArray (xyz0, 1);
   if (valc != (ArrayObject *)NULL)
      freeArray (valc0, 1);
   for (i = 0; i < nlist1; i ++) {
      xyzc_newd [3 * (listd [list1d [i]] - maskd [list1d [i]])] =
         xyz1d [3 * i];
      xyzc_newd [3 * (listd [list1d [i]] - maskd [list1d [i]]) + 1] =
         xyz1d [3 * i + 1];
      xyzc_newd [3 * (listd [list1d [i]] - maskd [list1d [i]]) + 2] =
         xyz1d [3 * i + 2];
      if (valc != (ArrayObject *)NULL) {
         if (atype == 'd')
            valc_newd [listd [list1d [i]] - maskd [list1d [i]]] = valc1d [i];
         else if (atype == 'b')
            valc_newc [listd [list1d [i]] - maskd [list1d [i]]] = valc1c [i];
      }
      }
   freeArray (mask, 1);
   freeArray (list1, 1);
   freeArray (xyz1, 1);
   if (valc != (ArrayObject *)NULL)
      freeArray (valc1, 1);
   TRY (ndxs = allocateArray ( ( (int *) (last->data)) [last->size - 1],
      'i', 1), (int) NULL);
   ndxsd = (int *) (ndxs->data);
   for (i = 0, i0 = 0; i < last->size; i0 = ( (int *) (last->data)) [i ++])
      for (j = 0; i0 + j < ( (int *) (last->data)) [i]; j++)
         ndxsd [i0 + j] = i;
   /* N. B. the following removes ndxs and xold. */
   /* Note: allocate return values in caller's list. */
   TRY (* nvertc = SimpleHist (take (ndxs, xold, 1, FREE1, FREE1, 1), FREE1, 0),
      (int) NULL);
   nvertcd = (int *) ( (* nvertc)->data);

   if (freexyzc >= 0)
      freeArray (xyzc, freexyzc);
   if (freevalc >= 0)
      freeArray (valc, freevalc);
   clearFreeList (1); /* just in case; should be OK */
   return 1;
   }

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Initialize the module.  This should be the only symbol with
   external linkage. */
void initgistC (void)
{
  PyObject *m, *d, *sys_path;
  int i, n;

  if (already_initialized)
    return;
  m = Py_InitModule4 ("gistC", gist_methods,
		      gist_module_documentation,
		      (PyObject *) 0, PYTHON_API_VERSION);
  d = PyModule_GetDict (m);
  GistError = PyString_FromString ("gist.error");
  PyDict_SetItemString (d, "error", GistError);
  if (PyErr_Occurred ()) {
    Py_FatalError ("Cannot initialize module gist");
  }
#ifdef import_array
  import_array();
#endif
  PyOS_InputHook = YGEventHandler;
  AddFDispatcher (stdin, &KeyboardInput, 0);
  DISPLAY_SET_HANDLER (Xerror_longjmp);
  if (0 != Py_AtExit (CleanUpGraphics)) {
    fprintf (stderr, "Gist: Warning: Exit procedure not registered\n");
  }
  /* Default is to put occasional markers on curves. */
  GhGetLines ();
  gistA.dl.marks = 1;
  GhSetLines ();

  /* Default text is 14 point Helvetica. */
  GhGetText ();
  gistA.t.font = T_HELVETICA;
  gistA.t.height = 14.0 * ONE_POINT;
  GhSetText ();

/* Find that component of sys.path which ends in "/graphics/gist", and
   add it to gistPathDefault.
*/
  m = (PyObject *) PyImport_AddModule ("sys");
  d = PyModule_GetDict (m);
  sys_path = PyDict_GetItemString (d, "path");
  n = PySequence_Length(sys_path); /* sys.path is a list of strings. */
  for(i=0; i<n; i++){
    PyObject *op;
    char *s;
    op = PySequence_GetItem( sys_path, i );
    s = PyString_AsString (op);
    if( strstr(s, OUR_SPECIAL_DIR)){
      gistpath = (char *) malloc(strlen(gistPathDefault) + strlen(s) + 2);
      if( gistpath ){
	oldgistpath = gistPathDefault;
        strcpy(gistpath, gistPathDefault);
        strcat(gistpath, ":");
        strcat(gistpath, s);
        gistPathDefault = gistpath;
      }
    break;
    }
  }

  already_initialized = 1;
}

static void CheckDefaultPalette (void)
{
  GpColorCell *palette;
  GhGetPalette (curPlotter, &palette);
  if (!palette)
    paletteSize = GhReadPalette (curPlotter,
       defaultPalette ? defaultPalette : "earth.gp", &palette, maxColors);
}

static int CheckDefaultWindow (void)
{
  int i;
  for (i = 0; i < 8; i++)
    if (ghDevices[i].drawing) {
      if (!ghDevices[i].display && !ghDevices[i].hcp) {
	Drawing *drawing = ghDevices[i].drawing;
	ghDevices[i].drawing = 0;
	GdKillDrawing (drawing);
	curElement = -1;
      }
    }
  if (curPlotter < 0) {
    for (i = 0; i < 8; i++)
      if (ghDevices[i].drawing) {
	return (int) ERRSS ("graphics window killed -- "
			    "use window command to re-select");
      }
    ghDevices[0].drawing =
	GdNewDrawing (defaultStyle ? defaultStyle : "work.gs");
    curElement = -1;
    if (!ghDevices[0].drawing) {
      return (int) ERRSS ("failed to create drawing -- "
			   "Gist work.gs style sheet missing");
    }
    ghDevices[0].doLegends = defaultLegends;

#ifndef NO_XLIB
    ghDevices[0].display =
	DISPLAY_ENGINE (windowNames[0], 0, defaultDPI, (char *) 0);
    if (!ghDevices[0].display) {
      return (int) ERRSS ("failed to open X display or create X window");
    }
#else
    ghDevices[0].display = 0;
    ghDevices[0].hcp = hcpDefault;
    hcpDefault = 0;
#endif

    curPlotter = 0;
    GhSetPlotter (0);
  }
  return 1;
}

static int CheckPalette (void)
{
  int n = curPlotter;
  if (n >= 0 && !ghDevices[n].hcp) {
    if (!hcpDefault) {
      if (!SetHCPDefault ())
	return 0; /* Failure */
    }
    SetHCPPalette ();
  }
  return 1; /* Success */
}

/* Shutdown routine */
static void CleanUpGraphics (void)
{
  int n;
  if (hcpDefault) {
    GpKillEngine (hcpDefault);
  }
  for (n = 7; n >= 0; n--) {
    if (ghDevices[n].display)
      GpKillEngine (ghDevices[n].display);
    if (ghDevices[n].hcp)
      GpKillEngine (ghDevices[n].hcp);
  }
  if ( gistpath ){
    gistPathDefault = oldgistpath;
    free(gistpath);
  }
  PyOS_InputHook = 0;
}

/* Search for 0-origin zone index of (xx,yy) in mesh (x,y).
 * reg==0 is allowed, and reg need not have overreach zones.
 * Return -1 if (xx,yy) not within mesh, zone index otherwise.
 */
static long FindMeshZone (double xx, double yy, double *x, double *y,
			  int *reg, long ix, long jx)
{
  int wind = 0;		/* winding number of zone around (xx,yy) */
  long i, i0, ijx = ix * jx;
  /* coordinates of zone corners relative to (xx,yy) */
  double x00, y00, x01, y01, x10, y10, x11, y11;

  i0 = 0;
  for (i = ix + 1; i < ijx; i++) {
    if ((++i0) == ix) {
      i0 = 1;
      i++;
    }			/* skip 1st point of each row */
    if (reg && !reg[i])
      continue;			/* skip non-existent zones */
    /* initialize 00 and 10 points at beginning of a row or
       after a non-existent zone */
    x00 = x[i - ix - 1] - xx;
    y00 = y[i - ix - 1] - yy;
    x10 = x[i - 1] - xx;
    y10 = y[i - 1] - yy;

    x01 = x[i - ix] - xx;
    y01 = y[i - ix] - yy;
    x11 = x[i] - xx;
    y11 = y[i] - yy;

    /* hopefully spend most of the time in these fast scan loops
       with all four points in one quadrant relative to (xx,yy) */
    if (ix < 8)
      goto test;		/* pointless if rows too short */
    if (x11 < 0.0) {
      if (x01 < 0.0 && x10 < 0.0 && x00 < 0.0) {
      left:
	for (i++;; i++) {
	  if (i >= ijx)
	    goto done;
	  if ((++i0) == ix) {
	    i0 = 0;
	    goto abort;
	  }
	  if (reg && !reg[i])
	    goto abort;
	  x01 = x[i - ix] - xx;
	  x11 = x[i] - xx;
	  if (x01 >= 0.0 || x11 >= 0.0)
	    break;
	}
	y00 = y[i - ix - 1] - yy;
	y10 = y[i - 1] - yy;
	y01 = y[i - ix] - yy;
	y11 = y[i] - yy;
	if (y11 < 0.0) {
	  if (y01 < 0.0 && y10 < 0.0 && y00 < 0.0)
	    goto low;
	} else {
	  if (y01 >= 0.0 && y10 >= 0.0 && y00 >= 0.0)
	    goto high;
	}
	x00 = x[i - ix - 1] - xx;
	x10 = x[i - 1] - xx;
	goto test;
      }
    } else {
      if (x01 >= 0.0 && x10 >= 0.0 && x00 >= 0.0) {
      right:
	for (i++;; i++) {
	  if (i >= ijx)
	    goto done;
	  if ((++i0) == ix) {
	    i0 = 0;
	    goto abort;
	  }
	  if (reg && !reg[i])
	    goto abort;
	  x01 = x[i - ix] - xx;
	  x11 = x[i] - xx;
	  if (x01 < 0.0 || x11 < 0.0)
	    break;
	}
	y00 = y[i - ix - 1] - yy;
	y10 = y[i - 1] - yy;
	y01 = y[i - ix] - yy;
	y11 = y[i] - yy;
	if (y11 < 0.0) {
	  if (y01 < 0.0 && y10 < 0.0 && y00 < 0.0)
	    goto low;
	} else {
	  if (y01 >= 0.0 && y10 >= 0.0 && y00 >= 0.0)
	    goto high;
	}
	x00 = x[i - ix - 1] - xx;
	x10 = x[i - 1] - xx;
	goto test;
      }
    }

    if (y11 < 0.0) {
      if (y01 < 0.0 && y10 < 0.0 && y00 < 0.0) {
      low:
	for (i++;; i++) {
	  if (i >= ijx)
	    goto done;
	  if ((++i0) == ix) {
	    i0 = 0;
	    goto abort;
	  }
	  if (reg && !reg[i])
	    goto abort;
	  y01 = y[i - ix] - yy;
	  y11 = y[i] - yy;
	  if (y01 >= 0.0 || y11 >= 0.0)
	    break;
	}
	x00 = x[i - ix - 1] - xx;
	x10 = x[i - 1] - xx;
	x01 = x[i - ix] - xx;
	x11 = x[i] - xx;
	if (x11 < 0.0) {
	  if (x01 < 0.0 && x10 < 0.0 && x00 < 0.0)
	    goto left;
	} else {
	  if (x01 >= 0.0 && x10 >= 0.0 && x00 >= 0.0)
	    goto right;
	}
	y00 = y[i - ix - 1] - yy;
	y10 = y[i - 1] - yy;
	goto test;
      }
    } else {
      if (y01 >= 0.0 && y10 >= 0.0 && y00 >= 0.0) {
      high:
	for (i++;; i++) {
	  if (i >= ijx)
	    goto done;
	  if ((++i0) == ix) {
	    i0 = 0;
	    goto abort;
	  }
	  if (reg && !reg[i])
	    goto abort;
	  y01 = y[i - ix] - yy;
	  y11 = y[i] - yy;
	  if (y01 < 0.0 || y11 < 0.0)
	    break;
	}
	x00 = x[i - ix - 1] - xx;
	x10 = x[i - 1] - xx;
	x01 = x[i - ix] - xx;
	x11 = x[i] - xx;
	if (x11 < 0.0) {
	  if (x01 < 0.0 && x10 < 0.0 && x00 < 0.0)
	    goto left;
	} else {
	  if (x01 >= 0.0 && x10 >= 0.0 && x00 >= 0.0)
	    goto right;
	}
	y00 = y[i - ix - 1] - yy;
	y10 = y[i - 1] - yy;
	goto test;
      }
    }

  test:
    /* compute counterclockwise crossings of x==xx for each of the
       four edges of the zone */
    wind = 0;
    if ((x00 < 0.0) ^ (x01 < 0.0)) {
      if (x00 * y01 > x01 * y00)
	wind++;
      else
	wind--;
    }
    if ((x01 < 0.0) ^ (x11 < 0.0)) {
      if (x01 * y11 > x11 * y01)
	wind++;
      else
	wind--;
    }
    if ((x11 < 0.0) ^ (x10 < 0.0)) {
      if (x11 * y10 > x10 * y11)
	wind++;
      else
	wind--;
    }
    if ((x10 < 0.0) ^ (x00 < 0.0)) {
      if (x10 * y00 > x00 * y10)
	wind++;
      else
	wind--;
    }
    if (wind)
      break;			/* this zone winds around (xx,yy) */

  abort:
    continue;
  }
done:

  return wind ? i : -1;
}

/* Used only by plq() */
static void ForceNewline (void)
{
  if (printNow) {
    if (printLines++ < maxPrintLines)
      RawPrinter (printBuf);
    printNow = permitNow = 0;
    printBuf[0] = '\0';
  }
}

/* Return name of current hardcopy file. */
static char *GetHCPname (int n)
{
  if (n >= 0 && n < 8)
    return ghDevices[n].hcp ? hcpNames[n] : hcpNames[8];
  else
    return hcpNames[8];
}

static void GetPCrange (double *zmn, double *zmx, double *z, int *reg,
			int region, long iMax, long jMax)
{
  double zmin = 0.0, zmax = 0.0;
  long i, j, k/* , len = iMax * jMax */;
  int have_min_max = 0;
  
  for (i = 0, k = 0; i < iMax; i ++)  {
     for (j = 0; j < jMax; j ++) {
        if (reg ? (region ?
                   i != iMax - 1 && j != jMax - 1 &&
                   (reg [k] == region || reg [k+1] == region ||
                    reg [k + jMax] == region || reg [k + jMax + 1] == region) :
                   (reg [k] || (i != iMax - 1 && j != jMax - 1 &&
                    (reg [k+1] || reg [k + jMax] || reg [k + jMax + 1])))) : 1) {
           if (! have_min_max) {
              have_min_max = 1;
              zmin = zmax = z [k];
           }
           else {
              if (z [k] < zmin) zmin = z [k];
              else if (z [k] > zmax) zmax = z [k];
           }
        }
        k ++;
     }
  }

  if (!have_min_max)
     ERRSS (
      "Unable to find maximum and minimum of data??");
  *zmn = zmin;
  *zmx = zmax;
}

/* Used by setkw_fonttype() */
static int GetTypeface (char *s, int *f)
{
  int face = 0;
  while (*s) {
    if (*s == 'B' && !(face & T_BOLD))
      face |= T_BOLD;
    else if (*s == 'I' && !(face & T_ITALIC))
      face |= T_ITALIC;
    else
      return (int) ERRSS (
      "illegal font keyword suffix -- B is bold, I is italic");
    s++;
  }
  *f = face;
  return 1;
}

static void GetZCrange(double *zmn, double *zmx, double *z, int *reg,
 int region, long iMax, long jMax, int zCompressed)
{
  double zmin= 0.0, zmax= 0.0;
  long i, j= iMax-1;
  long len= (zCompressed? j : iMax)*(jMax-1);

  if (zCompressed) {
    long len= (iMax-1)*(jMax-1);
    if (reg) reg+= iMax+1;
    for (i=0 ; i<len ; i++) {	/* first loop finds first z */
      if (reg? (region? (*reg==region) : (*reg!=0)) : 1) {
	zmin= zmax= z[i];
	break;
      }
      if (reg) {
	if (!(--j)) { reg+= 2; j= iMax-1; }
	else reg++;
      }
    }
    if (reg) {
      if (!(--j)) { reg+= 2; j= iMax-1; }
      else reg++;
    }
    for (i++ ; i<len ; i++) {	/* second loop judges extreme values */
      if (reg? (region? (*reg==region) : (*reg!=0)) : 1) {
	if (zmin>z[i]) zmin= z[i];
	else if (zmax<z[i]) zmax= z[i];
      }
      if (reg) {
	if (!(--j)) { reg+= 2; j= iMax-1; }
	else reg++;
      }
    }

  } else {
    z+= iMax+1;			/* set_yx guarantees at least 2-by-2 */
    if (reg) reg+= iMax+1;
    for (i=1 ; i<len ; i++) {	/* first loop finds first z */
      if (--j) {
	if (reg? (region? (*reg==region) : (*reg!=0)) : 1) {
	  zmin= zmax= z[i];
	  break;
	}
      } else {
	j= iMax;
      }
    }
    for (i++ ; i<len ; i++) {	/* second loop judges extreme values */
      if (--j) {
	if (reg? (region? (*reg==region) : (*reg!=0)) : 1) {
	  if (zmin>z[i]) zmin= z[i];
	  else if (zmax<z[i]) zmax= z[i];
	}
      } else {
	j= iMax;
      }
    }
  }

  *zmn= zmin;
  *zmx= zmax;
}

/* Called from bytscl, pli, plf, and plfp */
static int GrabByteScale (
  PyObject *kwt[], char *keywrds[], double *scale, double *offset,
  double *zn, double *zx, double *z, int *reg,
  int region, long iMax, long jMax, int zCompressed)
{
  int top;
  double zmin = 0.0, zmax = 0.0;
  int minGiven = 0, maxGiven = 0;

  /* get any parameters specified as keywords */
  if (kwt[0]) {
    SETKW(kwt[0],  top,    setkw_integer, keywrds[0]);
  } else {
    top = paletteSize - 1;
  }

  if (kwt[1] && kwt[1] != Py_None) {
    minGiven = 1;
    SETKW(kwt[1],  zmin,    setkw_double, keywrds[1]);
  }
  if (kwt[2] && kwt[2] != Py_None) {
    maxGiven = 1;
    SETKW(kwt[2],  zmax,    setkw_double, keywrds[2]);
  }

  /* fill in zmin and zmax from data if not specified */
  if (!minGiven || !maxGiven) {
    double zmn, zmx;
    GetZCrange(&zmn, &zmx, z, reg, region, iMax, jMax, zCompressed);
    if (!minGiven) zmin = zmn;
    if (!maxGiven) zmax = zmx;
  }

  /* adjust zmin and zmax to avert numerical catastrophes */
  if (zmin > zmax) { double tmp = zmin; zmin = zmax; zmax = tmp; }
  else if (zmin == zmax) {
    if (zmin > 0.0) { zmin = 0.9999*zmin; zmax = 1.0001*zmax; }
    if (zmin < 0.0) { zmin = 1.0001*zmin; zmax = 0.9999*zmax; }
    else { zmin = -0.0001; zmax = 0.0001; }
  }
  *zn = zmin;
  *zx = zmax;

  /* adjust top value if it is silly */
  if (top < 0 || top > 255) top = 255;

  /* (byte value)= scale*(z cut off at zmin, zmax)+offset
     maps from z to interval [0, top] */
  *scale = (double)top/(zmax-zmin);
  *offset = zmin-(0.4999/(*scale));	  /* zmin->0.5, zmax->top+0.5 */
  return 1;
}

static PyObject *contour (PyObject * self, PyObject * args, PyObject * kd)
{
  GaQuadMesh mesh;         /* object which gets sent to Gist      */
  int n,
      npt,
      i,
      nlevels,             /* number of contour levels (1 or 2)   */
      dims [2],            /* dimensions of mesh                  */
      ntotal,              /* number of coordinates sent back     */
      own_triangle = 0,    /* set to 1 if I create temp tri array */
      region = 0;
  long nparts,
       * np;
  static char * contourKeys [] = {"region", 0};
  PyObject * olevels,
           * zop,
           * kwt [NELT (contourKeys) - 1],
           * retval;
  PyArrayObject * zap,
                * alevels;
  PyObject      * anp,
                * axcp,
                * aycp;
  GpReal levels [2],
         * lev,
         * xcp,
         * ycp;
  double * z;

  if (!pyMsh.y)
    return ERRSS("contour: no current mesh - use plmesh(y, x) to initialize");

  n = PyTuple_Size (args);
  /* contour (levels, z [, region = num]) */
  if (n != 2)
     return ERRSS ("contour requires 2 positional parameters (levels and z).");
  BUILD_KWT (kd, contourKeys, kwt);
  TRY (PyArg_ParseTuple (args, "OO", &olevels, &zop),
     ERRSS ("contour: unable to parse arguments."));
  GET_ARR (zap, zop, PyArray_DOUBLE, 2, PyObject *);
  dims [0] = A_DIM (zap, 0);
  dims [1] = A_DIM (zap, 1);
  if (dims [0] != A_DIM (pyMsh.y, 0) || dims [1] != A_DIM (pyMsh.y, 1)) {
     ERRSS ("z array must have same dimensions as mesh in contour.");
     return (PyObject *) NULL; }
  /* Provide a triangle if none supplied */
  if ( !pyMsh.triangle )
     TRY (pyMsh.triangle = (PyArrayObject *) PyArray_FromDims (2, dims, PyArray_SHORT),
        ERRSS ("contour: unable to create triangle array."));
  gistD.region = 0;
  SETKW (kwt [0], gistD.region, setkw_integer, contourKeys [0]);
  get_mesh (&mesh);

  /* Figure out the contour levels */
  if (isARRAY (olevels)) {
     GET_ARR (alevels, olevels, Py_GpReal, 1, PyObject *);
     lev = (GpReal *) A_DATA (alevels);
     nlevels = A_SIZE (alevels);
     if (nlevels > 2) {
        clearArrayList ();
        return ERRSS ("contour: only 1 or 2 levels allowed."); }
     for (i = 0; i < nlevels; i++)
        levels [i] = lev [i];
     removeFromArrayList ( (PyObject *) alevels);
     }
  /* levels argument can be scalar--allow Int or Float */
  else if (PyFloat_Check (olevels) || PyInt_Check (olevels)) {
     nlevels = 1;
     if (PyFloat_Check (olevels))
        levels [0] = (GpReal) PyFloat_AsDouble (olevels);
     else
        levels [0] = (GpReal) PyInt_AsLong (olevels);
     }
  else {
     clearArrayList ();
     return ERRSS ("contour: levels argument is wrong type."); }

  z = (double *) A_DATA (zap);
  ntotal = (nlevels == 2) ?
     GcInit2 (&mesh, gistD.region, z, levels, 30L, &nparts):
     GcInit1 (&mesh, region, z, levels [0], &nparts);

  /* The following is necessary because I must own references to objects */
  /* that go in the list to be returned. Py_None will be INCREF'ed 3     */
  /* times when it is put on this list.                                  */
  if ( !(retval = Py_BuildValue ("[O,O,O]", Py_None, Py_None, Py_None))) {
     clearArrayList ();
     return ERRSS ("contour: unable to create return value list."); }
  if (ntotal == 0) {
     clearArrayList ();
     /* return a list [None, None, None] */
     return retval; }
  
  /* the tedium of all the SAFE_FREEs and Py_DECREF's has been */
  /* subsumed by the memory management routines.               */
  /* create three arrays and their data, make sure DECREF'able */
  npt = (int) nparts;
  NEW_MEM (np, npt, long, PyObject *);
  RET_ARR ( anp, 1, &npt, PyArray_INT, (char *) np, PyObject *);
  SET_OWN (anp);
  NEW_MEM (xcp, ntotal, GpReal, PyObject *);
  RET_ARR ( axcp, 1, &ntotal, Py_GpReal, (char *) xcp, PyObject *);
  SET_OWN (axcp);
  NEW_MEM (ycp, ntotal, GpReal, PyObject *);
  RET_ARR ( aycp, 1, &ntotal, Py_GpReal, (char *) ycp, PyObject *);
  SET_OWN (aycp);
  i = GcTrace (np, xcp, ycp);
  if ( i != ntotal) {
     clearArrayList ();
     clearMemList ();
     return ERRSS ("contour: GcTrace has failed."); }
  /* For some reason, if PyList_SetItem fails, it returns -1. */
  if (own_triangle) {
    Py_DECREF (kwt [0]);
  }
  if (PyList_SetItem (retval, 0, anp) < 0 ||
      PyList_SetItem (retval, 1, aycp) < 0 ||
      PyList_SetItem (retval, 2, axcp) < 0) {
     clearArrayList ();
     clearMemList ();
     return ERRSS ("contour was unable to build return list."); }
  removeFromArrayList ( (PyObject *) zap);
  mem_list_length = 0;
  array_list_length = 0;
  return retval;
}

/* Called from the event dispatcher when keyboard input is detected. */
/* See also YGDispatch. */
static int KeyboardInput (FILE * file, void *context)
{
  yPendingIn = 1;		/* Cause DispatchEvents to return */
  return 1;			/* immediately when keyboard input is
				 * available.  */
}

#ifndef NO_MOUSE
static int MouseCallBack (Engine * engine, int system,
			  int release, GpReal x, GpReal y,
			  int butmod, GpReal xn, GpReal yn)
{
  int n = curPlotter;
  if (n < 0 || ghDevices[n].display != engine) {
    mouseError = 1;
    return 1;
  }
  if (!release) {
    mouseSystem = system;
    mouseButton = butmod;
    mouseX0 = x;
    mouseY0 = y;
    mouseX0ndc = xn;
    mouseY0ndc = yn;
  } else {
    mouseModifier = butmod;
    mouseX1 = x;
    mouseY1 = y;
    mouseX1ndc = xn;
    mouseY1ndc = yn;
  }
  return 0;
}
#endif

/* Used only by plq() */
static void PermitNewline (int nSpaces)
{
  if (printNow + nSpaces > printLength)
    ForceNewline ();
  else
    while (nSpaces--)
      printBuf[printNow++] = ' ';
  printBuf[printNow] = '\0';
  permitNow = printNow;
}

/* Used only by plq() */
static void PrintColor (char *line, int color, int suffix)
{
  if (color >= 0) {
    sprintf (line, "color= %d,", color);
    PrintFunc (line);
  } else if (color == FG_COLOR)
    PrintFunc ("color= \"fg\"");
  else if (color == BG_COLOR)
    PrintFunc ("color= \"bg\"");
  else if (color == RED_COLOR)
    PrintFunc ("color= \"red\"");
  else if (color == GREEN_COLOR)
    PrintFunc ("color= \"green\"");
  else if (color == BLUE_COLOR)
    PrintFunc ("color= \"blue\"");
  else if (color == CYAN_COLOR)
    PrintFunc ("color= \"cyan\"");
  else if (color == MAGENTA_COLOR)
    PrintFunc ("color= \"magenta\"");
  else if (color == YELLOW_COLOR)
    PrintFunc ("color= \"yellow\"");
  else if (color == GREEN_COLOR)
    PrintFunc ("color= \"green\"");
  else
    PrintFunc ("color= <unknown>");
  PrintSuffix (suffix);
}

/* Used only by plq() */
static void PrintFunc (const char *s)
{
  long len = strlen (s);
  while (printNow + len > printLength) {
    if (permitNow) {
      char savec[2];
      int i = permitNow, j = 1;
      savec[0] = printBuf[i];
      printBuf[i++] = '\0';
      if (printLines++ < maxPrintLines)
	RawPrinter (printBuf);
      printBuf[0] = savec[0];
      while (i <= printNow)
	printBuf[j++] = printBuf[i++];
      printNow -= permitNow;
      permitNow = 0;
    } else {
      long nhere = printLength - printNow - 1;
      char movec = '\0';
      if (nhere > 0) {
	strncpy (&printBuf[printNow], s, nhere);
	s += nhere;
	len -= nhere;
      } else if (nhere < 0) {	/* only -1 is possible */
	movec = printBuf[printLength - 1];
      }
      strcpy (&printBuf[printLength - 1], "\\");
      if (printLines++ < maxPrintLines)
	RawPrinter (printBuf);
      if (nhere >= 0) {
	printNow = 0;
	printBuf[0] = '\0';
      } else {
	printNow = 1;
	printBuf[0] = movec;
	printBuf[1] = '\0';
      }
    }
  }
  strcpy (&printBuf[printNow], s);
  printNow += len;
}

/* Used only by plq() */
static void PrintHideLegend (char *line, int type)
{
  int offset = 0;
  char marker[5];
  marker[0] = '\0';
  sprintf (line, "hide= %d,", gistD.hidden);
  PrintFunc (line);
  ForceNewline ();
  if ((type == 1 || type == 7) && gistD.legend && gistD.legend[0] == '\001') {
    marker[0] = '\\';
    marker[1] = marker[2] = '0';
    marker[3] = '1';
    marker[4] = '\0';
    offset = 1;
  }
  sprintf (line, "legend= \"%s%.104s\",", marker,
	   gistD.legend ? gistD.legend + offset : "");
  PrintFunc (line);
  ForceNewline ();
}

/* Used only by plq() */
static void PrintInit (int (*puts) (const char *))
{
  RawPrinter = puts;
  printNow = permitNow = 0;
  printLines = 0;
  printBuf[0] = '\0';
}

/* Used only by plq() */
static void PrintMarks (char *line, int suffix)
{
  sprintf (line, "marks= %d,  mcolor= %d,  ", gistA.dl.marks, gistA.m.color);
  PrintFunc (line);
  if (gistA.m.type <= ' ' || gistA.m.type >= 0xff)
    sprintf (line, "marker= '\\%o',", gistA.m.type);
  else
    sprintf (line, "marker= '%c',", gistA.m.type);
  PrintFunc (line);
  ForceNewline ();
  sprintf (line,
	   "  msize= %.2f, mspace= %.5f, mphase= %.5f",
	   Safe_dbl (gistA.m.size),
	   Safe_dbl (gistA.dl.mSpace), Safe_dbl (gistA.dl.mPhase));
  PrintFunc (line);
  PrintSuffix (suffix);
}

/* Used only by plq() */
static void PrintRegion (char *line, int suffix)
{
  sprintf (line, "region= %d", gistD.region);
  PrintFunc (line);
  PrintSuffix (suffix);
}

/* Used only by plq() */
static void PrintSuffix (int suffix)
{
  if (suffix == 1)
    PrintFunc (",  ");
  else if (suffix == 3)
    PrintFunc (",");
  if (suffix & 2)
    ForceNewline ();
}

/* Used only by plq() */
static void PrintTypeWidth (char *line, int suffix)
{
  if (gistA.l.type == L_NONE)
    PrintFunc ("type= \"none\"");
  else if (gistA.l.type == L_SOLID)
    PrintFunc ("type= \"solid\"");
  else if (gistA.l.type == L_DASH)
    PrintFunc ("type= \"dash\"");
  else if (gistA.l.type == L_DOT)
    PrintFunc ("type= \"dot\"");
  else if (gistA.l.type == L_DASHDOT)
    PrintFunc ("type= \"dashdot\"");
  else if (gistA.l.type == L_DASHDOTDOT)
    PrintFunc ("type= \"dashdotdot\"");
  else
    PrintFunc ("type= <unknown>");
  sprintf (line, ",  width= %.2f", Safe_dbl (gistA.l.width));
  PrintFunc (line);
  PrintSuffix (suffix);
}

static GpColor *PushColors(double *z, long len, double zmin, double zmax,
 double scale, double offset)
{
  long i;
  double zz;
  GpColor *zc = (GpColor *) malloc (len * sizeof(GpColor));
  if (!zc) return (GpColor *) PyErr_NoMemory();

  for (i = 0 ; i < len ; i++) {
    zz = z[i];
    if (zz < zmin) zz = zmin;
    else if (zz > zmax) zz = zmax;
    zc[i] = (int) ((zz - offset) * scale);
  }
  return zc;
}

/* Used only by plq() */
static double Safe_dbl (double x)
{
  if (x > 1000.0)
    return 1000.0;
  else if (x < -1000.0)
    return -1000.0;
  else
    return x;
}

static int SetHCPDefault (void)
{
  int i, j;
  FILE *f;
  char hcpName[16];

  if (!hcpPSdefault)
    strcpy(hcpName, "Aa00.cgm");
  else
    strcpy(hcpName, "Aa00.ps");

  for (j = 'A'; j <= 'Z'; j++) {
    hcpName[0] = j;
    for (i = 'a'; i <= 'z'; i++) {
      hcpName[1] = i;
      if ((f = fopen (hcpName, "rb")))
	fclose (f);
      else
	goto got1;
    }
  }
  return (int) ERRSS (
		"you appear to have Aa00 through Zz00 hcp files -- clean up");

got1:
  if (!hcpPSdefault)
    hcpDefault = GpCGMEngine ("Pygist default", 0, hcpDump,
			    SetHCPname (-1, hcpName));
  else
    hcpDefault = GpPSEngine ("Pygist default", 0, hcpDump,
			    SetHCPname (-1, hcpName));

  if (!hcpDefault)
    return (int) ERRSS ("failed to create default hcp file");

  return 1;
}

static char *SetHCPname (int n, char *name)
{
  char *now;
  if (n < 0 || n > 7)
    n = 8;
  now = hcpNames[n];
  hcpNames[n] = expand_pathname (name);
  if (now)
    free (now);
  return hcpNames[n];
}

/* Xlib errors call back this function, which hopefully will have someplace
 * to longjmp back to.  So if a plot routine will risk an Xlib error,
 * it's important that it should SETJMP0 as it begins.
 */
static void Xerror_longjmp (char *message)
{
  PyErr_SetString (GistError, message);
  longjmp (jmpbuf, 1);
}

/* Call Gist Event Dispatcher. */
static void YGDispatch (void)
{
  GhBeforeWait ();		/* be sure X displays are updated */
  yPendingIn = 0;
  DispatchEvents ();		/* use Gist dispatcher */
}

static int YGEventHandler(void) { YGDispatch(); return 0; }

/* Used only by mouse() */
static int YPrompt(const char *s)
{
  int val = fputs(s, stdout);
  fflush (stdout);
  return val;
}

static PyObject *animate (PyObject * self, PyObject * args)
{
  int i = 3;			/* default is to toggle */

  SETJMP0;
  if (!PyArg_ParseTuple (args, "|i", &i))
    return ERRSS ("Animate takes zero or one argument.");
  PyFPE_START_PROTECT("animate", return 0)
  CheckDefaultWindow ();
  GhFMAMode (2, i);
  PyFPE_END_PROTECT
  Py_INCREF (Py_None);
  return Py_None;
}

/* Build an array of pointers to the keyword values in kwt[],
 * or nil, if that keyword wasn't given in the command arguments.
 * Return -1 on failure, number of keywords set otherwise.
 */
static int build_kwt (PyObject *kd, char *kwlist[], PyObject * kwt[])
{
  int i, n, nkw_set = 0;
  char *kw;
  PyObject *kob, *keylist;
  char *kword, errstr[256];

  for (i = 0; (kw = kwlist[i]); i++) { kwt[i] = 0; }
  if (!PyMapping_Check (kd))
    return 0; /* No keywords were passed. */

  /* Check that all keywords passed are legal for this command. */
  keylist = PyMapping_Keys (kd);
  n = PyList_Size (keylist);
  for (i = 0; i < n; i++) {
    kob = PySequence_GetItem (keylist, i);
    kword = PyString_AsString (kob);
    if (!verify_kw (kword, kwlist)) {
      sprintf (errstr, "Unrecognized keyword: %s", kword);
      (void) ERRSS (errstr);
      return -1;
    }
  }
  Py_DECREF(keylist);

  /* Ok, all keywords were legal.  Now store pointers to their value.
   * Note that PyDict_GetItemString() puts 0 in kwt[i] if
     that key isn't found. */
  for (i = 0; (kw = kwlist[i]); i++) {
    if((kwt[i] = PyDict_GetItemString (kd, kw)))
      ++nkw_set;
    /* I tried PyMapping_GetItemString() above, but kept getting
     * "KeyError: wait" messages back from Python.
     */
  }
  return nkw_set;
}

static PyObject *bytscl (PyObject * self, PyObject * args, PyObject * kd)
{
  static char *bsKeys[] = { "top", "cmin", "cmax", 0 };
  PyObject *zop, *kwt[NELT (bsKeys) - 1];
  PyArrayObject *zap, *zcap;
  double *z, zmin, zmax, scale, offset;
  GpColor *zc, *zc1;
  int i, len;

  if (!PyArg_ParseTuple (args, "O", &zop))
    return ERRSS ("bytscl requires exactly one non-keyword argument");

  TRY (addToArrayList((PyObject *)(zap = (PyArrayObject *)
      PyArray_ContiguousFromObject (zop, Py_GpReal, 1, 0))),
      (PyObject *)PyErr_NoMemory ());
  z = (GpReal *) A_DATA (zap);
  len = A_SIZE (zap);

  BUILD_KWT(kd, bsKeys, kwt);
  TRY (GrabByteScale(&kwt[0], &bsKeys[0], &scale, &offset, &zmin, &zmax,
     z, (int *) 0, 0, len + 1, 2L, 1), (PyObject *) NULL);
  TRY (zc = PushColors(z, len, zmin, zmax, scale, offset), (PyObject *) NULL);
  NEW_ARR (zcap, zap->nd, zap->dimensions, Py_GpColor, PyObject *);
  Py_DECREF (zap);
  zc1 = (GpColor *) A_DATA (zcap);
  for (i = 0; i < len; i++) zc1[i] = zc[i];

  array_list_length = 0;
  free(zc);
  return (PyObject *) zcap;
}

/* Zero out the global pyMsh struct, and free any mesh arrays. */
static void clear_pyMsh(void)
{
  Py_XDECREF (pyMsh.y);
  Py_XDECREF (pyMsh.x);
  Py_XDECREF (pyMsh.reg);
  Py_XDECREF (pyMsh.triangle);
  pyMsh.y = pyMsh.x = pyMsh.reg = pyMsh.triangle = 0;
}

static PyObject *current_window (PyObject * self, PyObject * args)
{
  return PyInt_FromLong (curPlotter);
}

/* The following routine has been added to check the integrity of what */
/* we think might be a NumPy array, including looking at addresses.    */
static PyObject *debug_array (PyObject * self, PyObject * args)
{
 PyObject *oarray;
 PyArrayObject * aarray;
 int i;
 int max;
 long mmax;
 TRY (PyArg_ParseTuple (args, "O", &oarray),
    ERRSS ("debug_array: argument should be one PyObject*."));
 printf("Value of input pointer is %x.", oarray); fflush (stdout);
 printf(" Reference count %d, size %d.\n", oarray->ob_refcnt, oarray->ob_type->ob_size);
 fflush (stdout);
 if (! isARRAY (oarray)) {
    return (PyObject *) (ERRSS ("debug_array: argument should be a NumPy array."));
    }
 aarray = (PyArrayObject *) oarray;
 printf("Data pointer: %x; nd %d; dim1 %d; type %c.\n", aarray->data,
   aarray->nd, aarray->dimensions [0], aarray->descr->type); fflush (stdout);
 if (aarray->descr->type == 'i') {
    printf ("%d ", ( (int *)(aarray->data)) [0]); fflush (stdout);
    for (i = 1, max = ( (int *)(aarray->data)) [0]; i < aarray->dimensions [0]; i ++){
       if ( ( (int *)(aarray->data)) [i] > max) max = ( (int *)(aarray->data)) [i];
       printf ("%d ", ( (int *)(aarray->data)) [i]);
       if (i % 10 == 0) printf ("\n");
       fflush (stdout);
       }
    printf ("maximum value is %d.\n", max); fflush (stdout);
    }
 else if (aarray->descr->type == 'l') {
    printf ("%ld ", ( (long *)(aarray->data)) [0]); fflush (stdout);
    for (i = 1, mmax = ( (long *)(aarray->data)) [0]; i < aarray->dimensions [0]; i ++){
       if ( ( (long *)(aarray->data)) [i] > mmax) mmax = ( (long *)(aarray->data)) [i];
       printf ("%ld ", ( (long *)(aarray->data)) [i]);
       if (i % 10 == 0) printf ("\n");
       fflush (stdout);
       }
    printf ("maximum value is %ld.\n", mmax); fflush (stdout);
    }
 Py_INCREF (Py_None);
 return Py_None;
}

/* Expand tildes and environment variables in pathnames. */
static char *expand_pathname (const char *name)
{
  PyObject *m, *d, *xpnduser, *xpndvars;
  DECL_ZERO (PyObject *, p1);
  DECL_ZERO (PyObject *, p2);
  DECL_ZERO (PyObject *, p3);
  DECL_ZERO (PyObject *, p4);
  char *path, *errstr = (char *) NULL, *name2, *name3;

  if (!name) return 0;

  /* Next four object refs are borrowed, and should not be DECREF'ed.
   * I know that module "os" has already been imported by gist.py, thus
   * can safely call PyImport_AddModule.  Otherwise, would need to call
   * PyImport_ImportModule to get posixmodule initialized.
   */
  TRY (m = (PyObject *) PyImport_AddModule ("posixpath"), (char *) NULL);
  TRY (d = PyModule_GetDict (m), (char *) NULL);
  TRY (xpnduser = PyDict_GetItemString (d, "expanduser"), (char *) NULL);
  TRY (xpndvars = PyDict_GetItemString (d, "expandvars"), (char *) NULL);

  /*
   * Here's a scorecard to keep track of the variables which follow:
   * "p1" is the PyObject (tuple of length 1) built from the string "name".
   * "p2" is the PyObject after expansion of any leading ~.
   * "name2" is the C string corresponding to "p2".
   * "p3" is the PyObject (tuple of length 1) built from the string "name2".
   * "p4" is the PyObject after expansion of any environment variables.
   * "name3" is the C string corresponding to "p4".
   * Finally, "path" is a copy of "name3" in memory that we own.
  */

  if (!(p1 = Py_BuildValue ("(s)", name))) goto errexit;
  if (!(p2 = PyEval_CallObject (xpnduser, p1))) goto errexit;
  name2 = PyString_AsString (p2);
  if (!(p3 = Py_BuildValue ("(s)", name2))) goto errexit;
  if (!(p4 = PyEval_CallObject (xpndvars, p3))) goto errexit;
  name3 = PyString_AsString (p4);

  path = (char *) malloc (1 + strlen (name3));
  if (path) (void)strcpy (path, name3);

  DECREF_AND_ZERO(p1);
  DECREF_AND_ZERO(p2);
  DECREF_AND_ZERO(p3);
  DECREF_AND_ZERO(p4);
  return path;

errexit:
  SETERR ("error in expand_path");
  DECREF_AND_ZERO(p1);
  DECREF_AND_ZERO(p2);
  DECREF_AND_ZERO(p3);
  DECREF_AND_ZERO(p4);
  return 0;
}

static PyObject *g_fma (PyObject * self, PyObject * args)
{
  SETJMP0;
  if (!CheckDefaultWindow ())
    return 0;
  if (hcpOnFMA) {
    if (!CheckPalette ())
      return 0;
  }
  curElement = -1;
  GhFMA ();
  Py_INCREF (Py_None);
  return Py_None;
}

/* Set pointers in the GaQuadMesh struct from values in the current
 * pyMsh struct.  Naturally, pyMsh must be fully initialized before
 * this is called.
 */
static void get_mesh(GaQuadMesh *m)
{
  m->iMax = A_DIM (pyMsh.y, 1);
  m->jMax = A_DIM (pyMsh.y, 0);
  m->y = (GpReal *) A_DATA (pyMsh.y);
  m->x = (GpReal *) A_DATA (pyMsh.x);
  m->reg = (int *) A_DATA (pyMsh.reg);
  if (isARRAY (pyMsh.triangle))
    m->triangle = (short *) A_DATA (pyMsh.triangle);
  else
    m->triangle = 0; /* Gist will provide a default in this case. */
}

static PyObject* get_slice2_precision (PyObject * self, PyObject * args)
{
 if (PyTuple_Size (args) > 0)
    return ERRSS ("get_slice2_precision takes no arguments.") ;
 return Py_BuildValue ( "d", _slice2_precision) ;
}

static PyObject *gridxy (PyObject * self, PyObject * args, PyObject * kd)
{
  int xgrid = 0, ygrid = 0, narg;
  static char *gridKeys[]= { "color", "type", "width", 0 };
  PyObject * kwt[NELT(gridKeys) - 1];

  SETJMP0;
  if (!PyArg_ParseTuple (args, "|ii", &xgrid, &ygrid)) {
    return ERRSS ("gridxy takes zero, one or two non-keyword arguments.");
  }
  /* If one argument is given, use it for both x and y. */
  if((narg = PyTuple_Size(args)) == 1)
    ygrid = xgrid;

  CheckDefaultWindow ();

  BUILD_KWT(kd, gridKeys, kwt);
  SETKW(kwt[0], gistD.ticks.horiz.gridStyle.color, setkw_color,    gridKeys[0]);
  SETKW(kwt[0], gistD.ticks.vert.gridStyle.color,  setkw_color,    gridKeys[0]);
  SETKW(kwt[1], gistD.ticks.horiz.gridStyle.type,  setkw_linetype, gridKeys[1]);
  SETKW(kwt[1], gistD.ticks.vert.gridStyle.type,   setkw_linetype, gridKeys[1]);
  SETKW(kwt[2], gistD.ticks.horiz.gridStyle.width, setkw_double,   gridKeys[2]);
  SETKW(kwt[2], gistD.ticks.vert.gridStyle.width,  setkw_double,   gridKeys[2]);

  if(narg > 0){
    gistD.ticks.horiz.flags &= ~(GRID_F | GRID_O);
    if (xgrid == 1)
      gistD.ticks.horiz.flags |= GRID_F;
    else if (xgrid == 2)
      gistD.ticks.horiz.flags |= GRID_O;
    if (xgrid & 0x200) {
      gistD.ticks.horiz.flags = (xgrid & 0x1ff);
      gistD.ticks.frame = (xgrid & 0x400) != 0;
    }
    gistD.ticks.vert.flags &= ~(GRID_F | GRID_O);
    if (ygrid & 1)
      gistD.ticks.vert.flags |= GRID_F;
    else if (ygrid & 2)
      gistD.ticks.vert.flags |= GRID_O;
    if (ygrid & 0x200) {
      gistD.ticks.vert.flags = (ygrid & 0x1ff);
      gistD.ticks.frame = (ygrid & 0x400) != 0;
    }
  }
  GdSetPort ();
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *hcp (PyObject * self, PyObject * args)
{
  SETJMP0;
  PyFPE_START_PROTECT("hcp", return 0)
  CheckDefaultWindow ();
  CheckPalette ();
  GhHCP ();
  PyFPE_END_PROTECT
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *hcp_file (PyObject * self, PyObject * args, PyObject *kd)
{
  Engine *engine = hcpDefault;
  char *hcp = 0;
  int gotDump = 0;
  static char *hcpKeys[]= { "dump", "ps", 0 };
  PyObject * kwt[NELT(hcpKeys) - 1];

  if (!PyArg_ParseTuple (args, "|s", &hcp)) {
    return ERRSS ("Bad arguments for hcp_file.");
  }

  BUILD_KWT(kd, hcpKeys, kwt);
  gotDump = (kwt[0] != 0);
  SETKW(kwt[0], hcpDump,      setkw_boolean, hcpKeys[0]);
  SETKW(kwt[1], hcpPSdefault, setkw_boolean, hcpKeys[1]);

  if (hcp) {
    long len = strlen (hcp);

    if (engine) {
      hcpDefault = 0;
      GpKillEngine (engine);
      SetHCPname (-1, (char *) 0);
      engine = 0;
    }
    if (len > 3 && strcmp (&hcp[len - 3], ".ps") == 0) {
      engine = GpPSEngine ("Pygist default", 0, hcpDump, SetHCPname (-1, hcp));
      if (!engine) {
	return ERRSS ("failed to create PostScript file");
      }
    } else if (len > 0) {
      engine = GpCGMEngine ("Pygist default", 0, hcpDump, SetHCPname (-1, hcp));
      if (!engine) {
	return ERRSS ("failed to create binary CGM file");
      }
    }
    hcpDefault = engine;
  } else if (gotDump) {
    GhDumpColors (-1, 1, hcpDump);
  }
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *hcp_finish (PyObject * self, PyObject * args)
{
  /* just return name of current hcp file */
  int n = curPlotter;
  Engine *engine;
  PyObject *p;

  if (!PyArg_ParseTuple (args, "|i", &n)) {
    return ERRSS ("Bad argument for hcp_finish.");
  }
  if (n < -1 || n > 7) {
    return ERRSS ("hcp_finish argument must be -1 through 7 inclusive");
  }
  p = PyString_FromString (GetHCPname (n));

  if (n >= 0)
    engine = ghDevices[n].hcp ? ghDevices[n].hcp : hcpDefault;
  else
    engine = hcpDefault;
  if (engine) {
    if (engine == hcpDefault) {
      hcpDefault = 0;
    } else {
      ghDevices[n].hcp = 0;
    }
    GpKillEngine (engine);
    SetHCPname (n, (char *) 0);
  }
  return p;
}

static PyObject *hcpoff (PyObject * self, PyObject * args)
{
  SETJMP0;
  PyFPE_START_PROTECT("hcpoff", return 0)
  CheckDefaultWindow ();
  hcpOnFMA = 0;
  GhFMAMode (0, 2);
  PyFPE_END_PROTECT
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *hcpon (PyObject * self, PyObject * args)
{
  SETJMP0;
  PyFPE_START_PROTECT("hcpon", return 0)
  CheckDefaultWindow ();
  hcpOnFMA = 1;
  GhFMAMode (1, 2);
  PyFPE_END_PROTECT
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *limits (PyObject * self, PyObject * args, PyObject * kd)
{
  /* NB-- If the plot has not been displayed yet, this will not retrieve
     the latest extreme values calculated by GdScan.  Nevertheless,
     it DOES retrieve the precise state of the limits at the time
     of this call, and retoring them will work correctly.  */

  double old_limits[4], new_limits[4];
  PyObject *xmin_ob = 0, *xmax_ob = 0;
  PyObject *ymin_ob = 0, *ymax_ob = 0;
  int nkw, old_flags, new_flags = 0, changed = 0, j;
  int square, nice, g_restrict;
  static char *limKeys[]= { "square", "nice", "restrict", 0 };
  PyObject * kwt[NELT(limKeys) - 1];

  SETJMP0;
  if (!PyArg_ParseTuple (args, "|OOOO",
			 &xmin_ob, &xmax_ob, &ymin_ob, &ymax_ob)) {
    return ERRSS ("limits takes at most 4 non-keyword arguments.");
  }
  /* retrieve current limits and flags */
  GdGetLimits ();
  old_limits[0] = gistD.limits.xmin;
  old_limits[1] = gistD.limits.xmax;
  old_limits[2] = gistD.limits.ymin;
  old_limits[3] = gistD.limits.ymax;
  old_flags = gistD.flags;

  if (xmin_ob && PyTuple_Check (xmin_ob)) {  /* Restore saved limits. */
    if (PyMapping_Check (kd)) {
      return ERRSS ("Keywords not allowed when restoring saved limits.");
    }
    if (!unpack_limit_tuple (xmin_ob, new_limits, &new_flags)) {
      return 0;
    }
    gistD.limits.xmin = new_limits[0];
    gistD.limits.xmax = new_limits[1];
    gistD.limits.ymin = new_limits[2];
    gistD.limits.ymax = new_limits[3];
    gistD.flags = new_flags;
    GdSetLimits ();
    return Py_BuildValue ("ddddi",
     old_limits[0], old_limits[1], old_limits[2], old_limits[3], old_flags);
  }

  if ((nkw = build_kwt (kd, limKeys, kwt)) > 0) { /* At least one keyword */
    SETKW(kwt[0], square,   setkw_boolean, limKeys[0]);
    SETKW(kwt[1], nice,     setkw_boolean, limKeys[1]);
    SETKW(kwt[2], g_restrict, setkw_boolean, limKeys[2]);

    if (kwt[0]) {
      if(square) gistD.flags |= D_SQUARE;
      else gistD.flags &= ~D_SQUARE;
    }
    if (kwt[1]) {
      if(nice) gistD.flags |= D_NICE;
      else gistD.flags &= ~D_NICE;
    }
    if (kwt[2]) {
      if(g_restrict) gistD.flags |= D_RESTRICT;
      else gistD.flags &= ~D_RESTRICT;
    }
    ++changed;

  } else if (-1 == nkw) { /* Error unpacking keyword dictionary */
    return 0;
  } else if (!xmin_ob) { /* No arguments nor keywords were passed. */
    gistD.flags = (D_XMIN | D_XMAX | D_YMIN | D_YMAX);
    ++changed;
  }

  if (xmin_ob) {
    j = set_limit (xmin_ob, &gistD.limits.xmin, &gistD.flags, D_XMIN);
    if(0 == j) /* Error */
      return ERRSS ("bad xmin argument: Use float or 'e'");
    else if(1 == j) /* Xmin changed or set to extreme. */
      ++changed;
  }
  if (xmax_ob) {
    j = set_limit (xmax_ob, &gistD.limits.xmax, &gistD.flags, D_XMAX);
    if(0 == j) /* Error */
      return ERRSS ("bad xmax argument: Use float or 'e'");
    else if(1 == j) /* Xmax changed or set to extreme. */
      ++changed;
  }
  if (ymin_ob) {
    j = set_limit (ymin_ob, &gistD.limits.ymin, &gistD.flags, D_YMIN);
    if(0 == j) /* Error */
      return ERRSS ("bad ymin argument: Use float or 'e'");
    else if(1 == j) /* Ymin changed or set to extreme. */
      ++changed;
  }
  if (ymax_ob) {
    j = set_limit (ymax_ob, &gistD.limits.ymax, &gistD.flags, D_YMAX);
    if(0 == j) /* Error */
      return ERRSS ("bad ymax argument: Use float or 'e'");
    else if(1 == j) /* Ymax changed or set to extreme. */
      ++changed;
  }

  if (changed) GdSetLimits ();

  return Py_BuildValue ("ddddi",
     old_limits[0], old_limits[1], old_limits[2], old_limits[3], old_flags);
}

static PyObject *logxy (PyObject * self, PyObject * args)
{
  int xflag = -1, yflag = -1, changed = 0;

  SETJMP0;
  if (!PyArg_ParseTuple (args, "|ii", &xflag, &yflag)) {
    return ERRSS ("Bad arguments for logxy.");
  }
  if (-1 != xflag)
    changed |= 1;
  if (-1 != yflag)
    changed |= 2;

  if (changed) {
    GdGetLimits ();
    if (changed & 1) {
      if (1 == xflag)
	gistD.flags |= D_LOGX;
      else
	gistD.flags &= ~D_LOGX;
    }
    if (changed & 2) {
      if (1 == yflag)
	gistD.flags |= D_LOGY;
      else
	gistD.flags &= ~D_LOGY;
    }
    GdSetLimits ();
  }
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *mesh_loc (PyObject * self, PyObject * args)
{
  long *zone;
  int i, n;
  GaQuadMesh mesh;
  double *x0 = 0, *y0 = 0;
  PyObject *y0op, *x0op;
  PyArrayObject *y0ap = 0, *x0ap = 0, *rap = 0;
  char *errstr = "mesh_loc requires arguments (y0, x0 [ , y, x [ ,ireg ] ])";
  struct { double x, y; } dbuf;

  if (PyTuple_Size (args) < 2)
    return ERRSS("mesh_loc requires at least two arguments");
  TRY (setvu_mesh (args, &y0op, &x0op, errstr), (PyObject *) NULL);
  if (!pyMsh.y)
    return ERRSS("No current mesh - set (y, x) first");
  get_mesh (&mesh);

  if (isARRAY (y0op)) {
    TRY (addToArrayList((PyObject *)(y0ap = (PyArrayObject *)
        PyArray_ContiguousFromObject (y0op, Py_GpReal, 1, 0))),
        (PyObject *)PyErr_NoMemory ());
    n = A_SIZE ( y0ap );
    TRY (addToArrayList((PyObject *)(x0ap = (PyArrayObject *)
        PyArray_ContiguousFromObject (x0op, Py_GpReal, 1, 0))),
        (PyObject *)PyErr_NoMemory ());
    if (n != A_SIZE ( x0ap )) {
      clearArrayList();
      return ERRSS ("(y0, x0) must be same size");
    }
    y0 = (GpReal *) A_DATA (y0ap);
    x0 = (GpReal *) A_DATA (x0ap);
  } else if (PyFloat_Check (y0op)) {
    y0 = &dbuf.y;
    x0 = &dbuf.x;
    y0[0] = PyFloat_AsDouble (y0op);
    x0[0] = PyFloat_AsDouble (x0op);
    n = 1;
  } else {
    return ERRSS ("(y0, x0) must be floating point scalars or arrays.");
  }

  NEW_ARR (rap, 1, &n, PyArray_LONG, PyObject *);
  zone = (long *) A_DATA (rap);

  for (i = 0; i < n; i++)
    zone[i] = 1 + FindMeshZone (x0[i], y0[i], mesh.x, mesh.y, mesh.reg,
      mesh.iMax, mesh.jMax);

  if (isARRAY (y0op)) {
    Py_DECREF (y0ap);
    Py_DECREF (x0ap);
  }
  array_list_length = 0;
  return PyArray_Return (rap);
}

static PyObject *mfit (PyObject * self, PyObject * args)
{
 /* Computes multiquadric fit to data; used for contour plotting
  * of random data. Calling sequence from Python:
  * zcplot = mfit (alpha, x, xcplot, y, ycplot, rsqmqd)
  * where alpha are the interpolation coefficients, x and y
  * are the original randomly distributed coordinates
  * (alpha, x, and y are all the same length, say nzcplot).
  * xcplot (nxcplot long) and ycplot (nycplot long) specify
  * an overlying rectangular mesh. rsqmod is a scalar peculiar
  * to the problem.                                              */
 int nxcplot,
     nycplot,
     nzcplot;
 double *x,
        *y,
        *alpha,
        *xcplot,
        *ycplot,
        *zcplot,
        rsqmqd;
 PyObject *oalpha,
          *ox,
          *oy,
          *oxcplot,
          *oycplot/* , */
/*           *ozcplot */;
 PyArrayObject *aalpha,
               *ax,
               *ay,
               *axcplot,
               *aycplot,
               *azcplot;
 int i,
     j,
     k,
     l,
     dims [2];
 double sum;
 
 TRY (PyArg_ParseTuple (args, "OOOOOd", &oalpha, &ox, &oxcplot, &oy,
   &oycplot, &rsqmqd), ERRSS ("mfit: unable to parse arguments."));
 GET_ARR (aalpha, oalpha, PyArray_DOUBLE, 1, PyObject *);
 GET_ARR (ax, ox, PyArray_DOUBLE, 1, PyObject *);
 GET_ARR (axcplot, oxcplot, PyArray_DOUBLE, 1, PyObject *);
 GET_ARR (ay, oy, PyArray_DOUBLE, 1, PyObject *);
 GET_ARR (aycplot, oycplot, PyArray_DOUBLE, 1, PyObject *);
 /* There is no other error checking, really. It is intended that
  * this routine be called only from Python code, not by the user. */
 nzcplot = A_DIM (aalpha, 0);
 nxcplot = A_DIM (axcplot, 0);
 dims [0] = nxcplot;
 nycplot = A_DIM (aycplot, 0);
 dims [1] = nycplot;
 x = (double *) A_DATA (ax);
 y = (double *) A_DATA (ay);
 xcplot = (double *) A_DATA (axcplot);
 ycplot = (double *) A_DATA (aycplot);
 alpha = (double *) A_DATA (aalpha);
 TRY (azcplot = (PyArrayObject *) PyArray_FromDims (2, dims, PyArray_DOUBLE),
    ERRSS ("mfit: unable to create return array."));
 zcplot = (double *) A_DATA (azcplot);
 l = 0;
 for (i = 0; i < nxcplot; i++)
    for (j = 0; j < nycplot; j ++) {
       sum = 0.0;
       for (k = 0; k < nzcplot; k ++)
          sum += alpha [k] *
             sqrt ( (x [k] - xcplot [i]) * (x [k] - xcplot [i]) +
                    (y [k] - ycplot [j]) * (y [k] - ycplot [j]) + rsqmqd);
       zcplot [l ++] = sum;
    }
 clearArrayList ();
 return PyArray_Return (azcplot);
}

static PyObject *mouse (PyObject * self, PyObject * args)
{
#ifdef DISPLAY_MOUSE
  char *prompt = 0;
  int system = -1, style = 0;
  int n = curPlotter;

  SETJMP0;
  if (n < 0 || !ghDevices[n].display)
    return ERRSS ("no current graphics window for mouse function");

  if (!PyArg_ParseTuple (args, "|iis", &system, &style, &prompt))
    return ERRSS ("call with (system, style, prompt)");

  GhWaitDisplay ();		/* otherwise can lock up */
  GhBeforeWait ();		/* be sure display is current */
  if (!prompt)
    YPrompt (defaultPrompts[style != 0]);
  else if (prompt[0])
    YPrompt (prompt);
  mouseError = 0;
  mouseError |= DISPLAY_MOUSE (ghDevices[n].display, style, system,
			       &MouseCallBack);
  if (!prompt || prompt[0])
    YPrompt ("\n");

  if (mouseError) {
    Py_INCREF (Py_None);
    return Py_None;
  } else {
    return Py_BuildValue ("ddddddddiii",
      mouseX0,    mouseY0,    mouseX1,    mouseY1,
      mouseX0ndc, mouseY0ndc, mouseX1ndc, mouseY1ndc,
      mouseSystem, mouseButton, mouseModifier);
  }
#else
  return ERRSS ("no mouse function in this version of Python/Gist");
#endif
}

static PyObject *palette (PyObject * self, PyObject * args, PyObject * kd)
{
  GpColorCell *palette = 0;
  static char *paletteKeys[] = { "ntsc", "query", 0 };
  PyObject * kwt[NELT(paletteKeys) - 1];
  unsigned char *red = 0, *green = 0, *blue = 0, *gray = 0;
  PyObject *rop, *gop, *bop, *grayop;
  PyArrayObject *rap = 0, *gap = 0, *bap = 0, *grayap = 0;
  int nred = 0, ngreen = 0, nblue = 0, ngray = 0;
  int i, nColors=0, nDevice, query = 0, ntsc = 0, len_match;
  Engine *engine;
  int sourceDevice = -2;
  char *filename = 0;
  char *errstr = "palette takes a filename, or source window number, or\n"
    "red, blue, green [, gray] color arrays to specify or query the palette";

  SETJMP0;			/* See Xerror_longjmp() */

  BUILD_KWT(kd, paletteKeys, kwt);
  SETKW(kwt[0], ntsc,  setkw_boolean, paletteKeys[0]);
  SETKW(kwt[1], query, setkw_boolean, paletteKeys[1]);

  switch (PyTuple_Size (args)) {
  case 4: /* (red, green, blue, gray) given */
    TRY (grayop = PyTuple_GetItem (args, 3), (PyObject *) NULL);
    GET_ARR (grayap, grayop, Py_GpColor, 1, PyObject *);
    ngray = A_SIZE (grayap);
    gray = (GpColor *) A_DATA (grayap);
    /* Fall through. */
  case 3: /* (red, green, blue) given */
    TRY (PyArg_ParseTuple (args, "OOO", &rop, &gop, &bop), (PyObject *) NULL);
    GET_ARR (rap, rop, Py_GpColor, 1, PyObject *);
    nred = A_SIZE (rap);
    red = (GpColor *) A_DATA (rap);

    GET_ARR (gap, gop, Py_GpColor, 1, PyObject *);
    ngreen = A_SIZE (gap);
    green = (GpColor *) A_DATA (gap);

    GET_ARR (bap, bop, Py_GpColor, 1, PyObject *);
    nblue = A_SIZE (bap);
    blue = (GpColor *) A_DATA (bap);

    /* Check for matched array sizes and set nColors. */
    len_match = (nred == ngreen && nred == nblue);
    if(ngray) { len_match = (len_match && nred == ngray); }
    if(len_match){
      nColors = nred;
    }else{
      clearArrayList ();
      return ERRSS ("red, green, blue, and gray arguments must be same length");
    }

    break;
  case 1: /* (filename) or (source_window_number) given */
    if (query) return ERRSS ("query requires (r,g,b) arrays as arguments");

    if (PyArg_ParseTuple (args, "s", &filename)) {
      break; /* call was palette (filename) */
    } else if (PyArg_ParseTuple (args, "i", &sourceDevice)) {
      if (sourceDevice < 0 || sourceDevice >= 7 ||
	  (!(engine = ghDevices[sourceDevice].display) &&
	   !(engine = ghDevices[sourceDevice].hcp)))
        return ERRSS ("specified palette source window does not exist");
      break;
    } else {
      return ERRSS (errstr);
    }
  default:
    return ERRSS (errstr);
  }

  TRY (CheckDefaultWindow(), (PyObject *) NULL);
  nDevice = curPlotter;
  engine = ghDevices[nDevice].display;
  if (!engine) engine = ghDevices[nDevice].hcp;

  if (query) { /* Retrieve the current palette. */
    nColors = GpGetPalette(engine, &palette);
    if (nColors > 256) {
      clearArrayList ();
      return ERRSS ("Gist palettes can never have more than 256 colors");
    }
    if (nColors > nred  || nColors > ngreen ||
        nColors > nblue || (ngray && nColors > ngray)) {
      clearArrayList ();
      return ERRSS ("arrays passed are too small to hold all the colors");
    }
    if (nColors > 0) {
      for (i = 0 ; i < nColors ; i++) {
	red[i]   = palette[i].red;
	green[i] = palette[i].green;
	blue[i]  = palette[i].blue;
      }
      if (ngray)
	for (i = 0 ; i < nColors ; i++)
	  gray[i]  = palette[i].gray;

    }
  } else { /* Set a new palette. */

    if (sourceDevice != nDevice) {
      /* be sure to preserve dump = 1 setting even if hcp palette
         is deleted */
      int dump;
      if (hcpDefault) dump = GhGetColorMode(hcpDefault);
      else dump = 0;
      GhDeletePalette(nDevice);
      paletteSize = 0;
      if (hcpDefault) GhDumpColors(-1, 1, dump);
    }
    if (red) {
      /* palette is unprotected against asynchronous interrupts...
	 fix this someday */
      palette = malloc (sizeof(GpColorCell) * nColors);
      for (i = 0 ; i < nColors ; i++) {
	palette[i].red = red[i];
	palette[i].green = green[i];
	palette[i].blue = blue[i];
	if (gray) palette[i].gray = gray[i];
      }
      if (!gray) {
	if (ntsc)
	  GpPutNTSC(nColors, palette);
	else
	  GpPutGray(nColors, palette);
      }
      GhSetPalette (nDevice, palette, nColors);
      paletteSize = nColors;

    } else if (filename) {
      nColors = GhReadPalette(nDevice, filename, &palette, maxColors);
      if (nColors <= 0) {
        char s[1024];
	sprintf(s, "%s: Gist palette not found", filename);
        clearArrayList ();
	return ERRSS (s);
      }
    }
    paletteSize = nColors;
  }
  clearArrayList ();
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *pyg_pause (PyObject * self, PyObject * args)
{
  long timeout;
  unsigned long mask = (unsigned long) (-1L);

  if (!PyArg_ParseTuple (args, "i", &timeout)) {
    return ERRSS ("Pause requires one integer argument.");
  }
  if (timeout < 0)
    timeout = 0;
  G_poll (2L, &mask, timeout);
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *plc (PyObject * self, PyObject * args, PyObject * kd)
{
  GaQuadMesh mesh;
  PyArrayObject *zap;
  PyObject *zop;
  int i, nr, nc, nLevels;
  GpReal *z = 0, *levels;
  static char *plcKeys[]= {
    "legend", "hide", "region", "color", "type", "width",
    "marks", "mcolor", "marker", "msize", "mspace", "mphase",
    "smooth", "triangle", "levs", 0 };
  PyObject * kwt[NELT(plcKeys) - 1];
  char *errstr =
    "plc requires 2D arguments (z [ , y, x, ireg, levs = levels ] )";

  SETJMP0;			/* See Xerror_longjmp() */

  if (PyTuple_Size (args) == 0)
    return ERRSS("plc requires at least one argument");
  BUILD_KWT(kd, plcKeys, kwt);
  TRY (setz_mesh (args, &zop, errstr, kwt[13]), (PyObject *) NULL);
  if (!pyMsh.y)
    return ERRSS("No current mesh - set (y, x) first");
  GET_ARR (zap, zop, Py_GpReal, 2, PyObject *);
  nr = A_DIM(zap, 0);
  nc = A_DIM(zap, 1);
  if (A_DIM (pyMsh.y, 0) != nr || A_DIM (pyMsh.y, 1) != nc) {
    clearArrayList ();
    return ERRSS("Z array must match (y, x) mesh arrays in shape");
  }
  z = (GpReal *) A_DATA (zap);
  get_mesh (&mesh);

  TRY(CheckDefaultWindow (), (PyObject *) NULL);

  GhGetLines();	/* Properties start from defaults for decorated polylines. */
  gistD.region = 0;

  SETKW(kwt[0],  gistD.legend,    setkw_string,   plcKeys[0]);
  SETKW(kwt[1],  gistD.hidden,    setkw_boolean,  plcKeys[1]);
  SETKW(kwt[2],  gistD.region,    setkw_integer,  plcKeys[2]);
  SETKW(kwt[3],  gistA.l.color,   setkw_color,    plcKeys[3]);
  if(kwt[3]) gistA.m.color = gistA.l.color;
  SETKW(kwt[4],  gistA.l.type,    setkw_linetype, plcKeys[4]);

  SETKW(kwt[5],  gistA.l.width,   setkw_double,   plcKeys[5]);
  SETKW(kwt[6],  gistA.dl.marks,  setkw_boolean,  plcKeys[6]);
  SETKW(kwt[7],  gistA.m.color,   setkw_color,    plcKeys[7]);
  SETKW(kwt[8],  gistA.m.type,    setkw_xinteger, plcKeys[8]);
  SETKW(kwt[9],  gistA.m.size,    setkw_double,   plcKeys[9]);

  SETKW(kwt[10], gistA.dl.mSpace, setkw_double,   plcKeys[10]);
  SETKW(kwt[11], gistA.dl.mPhase, setkw_double,   plcKeys[11]);
  SETKW(kwt[12], gistA.dl.smooth, setkw_boolean,  plcKeys[12]);
  /* kwt[13] ("triangle=") was handled by setz_mesh. */
  if(kwt[14]) { /* levs= keyword */
    PyArrayObject *lap;
    GpReal *lev;

    GET_ARR (lap, kwt[14], Py_GpReal, 1, PyObject *);
    lev = (GpReal *) A_DATA (lap);
    nLevels = A_SIZE (lap);
    levels = malloc (sizeof(GpReal) * nLevels);
    if(!levels) return PyErr_NoMemory();
    for(i = 0; i < nLevels; i++)
      levels[i] = lev[i];
    removeFromArrayList ( (PyObject *) lap);
  } else {
    double zmin, zmax, step;

    nLevels = 9;
    NEW_MEM (levels, nLevels, GpReal, PyObject *);
    PyFPE_START_PROTECT("plc", return 0)
    GetPCrange (&zmin, &zmax, z, mesh.reg, gistD.region, nc, nr);
    step = (zmax - zmin)/9.0;
    levels[0] = zmin/*  + 0.5 * step */; /*########TMP#########*/
    for(i = 1; i < nLevels; i++)
      levels[i] = levels[i-1] + step;
    PyFPE_END_PROTECT
  }

  curElement = -1;
  PyFPE_START_PROTECT("plc", return 0)
  curElement =
    GdContours (NOCOPY_MESH, &mesh, gistD.region, z, levels, (int)nLevels);
  PyFPE_END_PROTECT
  Py_DECREF (zap);
  SAFE_FREE (levels);
  array_list_length = 0;
  mem_list_length = 0;
  if (curElement < 0)
    return ERRSS ("Gist GdContour plotter failed");

  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *pldefault (PyObject * self, PyObject * args, PyObject * kd)
{
  static char *dfltKeys[]= {
    "color", "type", "width",
    "marks", "mcolor", "marker", "msize", "mspace", "mphase",
    "rays", "arrowl", "arroww", "rspace", "rphase",
    "font", "height", "orient", "justify", "opaque",
    "hollow", "aspect", "dpi", "style", "legends", "palette", "maxcolors",
    "edges", "ecolor", "ewidth", 0 };
  PyObject * kwt[NELT(dfltKeys) - 1];
  char *errstr = "pldefault takes no non-keyword arguments";
  int dummy, dpi, type;

  if(PyTuple_Size(args) > 0) return ERRSS (errstr);
  
  /* retrieve all default settings */
  GhGetLines();
  GhGetMesh();
  GhGetVectors();
  GhGetText();

  BUILD_KWT(kd, dfltKeys, kwt);
  SETKW(kwt[0],  gistA.l.color,     setkw_color,    dfltKeys[0]);
  SETKW(kwt[1],  gistA.l.type,      setkw_linetype, dfltKeys[1]);
  SETKW(kwt[2],  gistA.l.width,     setkw_double,   dfltKeys[2]);
  SETKW(kwt[3],  gistA.dl.marks,    setkw_boolean,  dfltKeys[3]);
  SETKW(kwt[4],  gistA.m.color,     setkw_color,    dfltKeys[4]);
  SETKW(kwt[5],  gistA.m.type,      setkw_xinteger, dfltKeys[5]);

  SETKW(kwt[6],  gistA.m.size,      setkw_double,   dfltKeys[6]);
  SETKW(kwt[7],  gistA.dl.mSpace,   setkw_double,   dfltKeys[7]);
  SETKW(kwt[8],  gistA.dl.mPhase,   setkw_double,   dfltKeys[8]);
  SETKW(kwt[9],  gistA.dl.rays,     setkw_boolean,  dfltKeys[9]);
  SETKW(kwt[10], gistA.dl.arrowL,   setkw_double,   dfltKeys[10]);

  SETKW(kwt[11], gistA.dl.arrowW,   setkw_double,   dfltKeys[11]);
  SETKW(kwt[12], gistA.dl.rSpace,   setkw_double,   dfltKeys[12]);
  SETKW(kwt[13], gistA.dl.rPhase,   setkw_double,   dfltKeys[13]);
  SETKW(kwt[14], gistA.t.font,      setkw_fonttype, dfltKeys[14]);
  if(kwt[15]) {
    SETKW(kwt[15], gistA.t.height,  setkw_double,   dfltKeys[15]);
    gistA.t.height *= ONE_POINT;
  }

  SETKW(kwt[16], gistA.t.orient,    setkw_integer,  dfltKeys[16]);
  if (!gistA.t.orient) {
    gistA.t.orient = TX_RIGHT;
  } else {
    if (gistA.t.orient == 1) gistA.t.orient = TX_UP;
    else if (gistA.t.orient == 2) gistA.t.orient = TX_LEFT;
    else if (gistA.t.orient == 3) gistA.t.orient = TX_DOWN;
    else {
      gistA.t.orient= TX_RIGHT;
      return ERRSS("orient= keyword must be 0, 1, 2, or 3");
    }
  }

  SETKW(kwt[17], dummy,             setkw_justify,  dfltKeys[17]);
  SETKW(kwt[18], gistA.t.opaque,    setkw_boolean,  dfltKeys[18]);
  SETKW(kwt[19], gistA.vect.hollow, setkw_boolean,  dfltKeys[19]);
  SETKW(kwt[20], gistA.vect.aspect, setkw_double,   dfltKeys[20]);

  if(kwt[21]) {
    SETKW(kwt[21], dpi,             setkw_integer,  dfltKeys[21]);
    if (dpi != 100 && dpi != 75)
      return ERRSS ("dpi=100 or dpi=75 are only legal values");
    else
      defaultDPI = dpi;
  }

  if(kwt[22]) {
    char *style;
    SAFE_FREE(defaultStyle);
    SETKW(kwt[22], style,           setkw_string,   dfltKeys[22]);
    if(style && style[0]){
      NEW_MEM (defaultStyle, strlen(style) + 1, char, PyObject *);
      strcpy(defaultStyle, style);
    }
  }

  SETKW(kwt[23], defaultLegends,    setkw_boolean,  dfltKeys[23]);

  if(kwt[24]) {
    char *name;
    SAFE_FREE(defaultPalette);
    SETKW(kwt[24], name,            setkw_string,   dfltKeys[24]);
    if(name && name[0]){
      NEW_MEM (defaultPalette, strlen(name) + 1, char, PyObject *);
      strcpy(defaultPalette, name);
    }
  }

  SETKW(kwt[25], maxColors,         setkw_integer,  dfltKeys[25]);
  if(kwt[26]) {
    SETKW(kwt[26], type,            setkw_boolean,  dfltKeys[26]);
    gistA.e.type = type ? L_SOLID : L_NONE;
  }
  SETKW(kwt[27], gistA.e.color,     setkw_color,    dfltKeys[27]);
  SETKW(kwt[28], gistA.e.width,     setkw_double,   dfltKeys[28]);

  /* store all default settings */
  GhSetLines();
  GhSetMesh();
  GhSetVectors();
  GhSetText();
  GhSetFill();

  mem_list_length = 0;
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *pldj (PyObject * self, PyObject * args, PyObject * kd)
{
  PyObject *op[4];
  PyArrayObject *ap[4];
  GpReal *d[4];
  int i, n;
  static char *pldjKeys[]= { "legend", "hide", "color", "type", "width", 0 };
  PyObject * kwt[NELT(pldjKeys) - 1];
  char *errstr = "pldj requires exactly four non-keyword arguments";

  SETJMP0;

  if (!PyArg_ParseTuple (args, "OOOO", &op[0], &op[1], &op[2], &op[3]))
    return ERRSS (errstr);
  
  for (i=0; i<4; i++)
    TRY (addToArrayList ((PyObject *)(ap[i] = (PyArrayObject *)
        PyArray_ContiguousFromObject (op[i], Py_GpReal, 1, 0))),
        (PyObject *)PyErr_NoMemory ());

  n = A_SIZE ( ap[0] );
  for (i=1; i<4; i++)
    if ( A_SIZE (ap[i]) != n) {
      clearArrayList ();
      return ERRSS ("pldj arguments must all be the same size");
    }

  TRY (CheckDefaultWindow (), (PyObject *) NULL);

  GhGetMesh();	/* Properties start from defaults for simple polylines. */
  BUILD_KWT(kd, pldjKeys, kwt);
  SETKW(kwt[0],  gistD.legend,    setkw_string,   pldjKeys[0]);
  SETKW(kwt[1],  gistD.hidden,    setkw_boolean,  pldjKeys[1]);
  SETKW(kwt[2],  gistA.l.color,   setkw_color,    pldjKeys[2]);
  SETKW(kwt[3],  gistA.l.type,    setkw_linetype, pldjKeys[3]);
  SETKW(kwt[4],  gistA.l.width,   setkw_double,   pldjKeys[4]);

  for (i=0; i<4; i++)
    d[i] = (GpReal *) A_DATA (ap[i]);

  curElement = -1;
  PyFPE_START_PROTECT("pldj", return 0)
  curElement = GdDisjoint (n, d[0], d[1], d[2], d[3]);
  PyFPE_END_PROTECT
  clearArrayList ();
  if (curElement < 0)
    return ERRSS ("Gist GdDisjoint plotter failed");

  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *pledit (PyObject * self, PyObject * args, PyObject * kd)
{
  int dummy, type = 0, n_element = 0, n_contour = 0;
  int changes = 0;
  static char *editKeys[] = {
    "legend", "hide",
    "color", "type", "width",
    "marks", "mcolor", "marker", "msize", "mspace", "mphase",
    "rays", "arrowl", "arroww", "rspace", "rphase", "closed", "smooth",
    "font", "height", "orient", "justify", "opaque",
    "hollow", "aspect", "region", "boundary", "levs", "scale", "scalem",
    "dx", "dy", "edges", "ecolor", "ewidth", "inhibit", 0 };
  PyObject *kwt[NELT (editKeys) - 1];
  char *legend = 0;

  switch (PyTuple_Size (args)) {
  case 2: /* (n_element, n_contour) given */
    TRY (PyArg_ParseTuple (args, "ii", &n_element, &n_contour),
       (PyObject *) NULL);
    break;
  case 1: /* (n_element) given */
    TRY (PyArg_ParseTuple (args, "i", &n_element), (PyObject *) NULL);
    break;
  case 0: /* () given */
    break;
  default:
    return ERRSS ("pledit function takes no more than two arguments");
  }

  /* Pygist uses 1-origin element numbering, Gist uses 0-origin */
  n_element--;
  n_contour--;

  if (n_element < 0) {
    if (curElement >= 0) {
      n_element = GdFindIndex (curElement);
      if (n_element < 0) {
	curElement = -1;
	return ERRSS ("lost current graphical element for pledit (BUG?)");
      }
    } else if (curElement == -6666) {
      n_element = curIX;
      n_contour = curIXc;
    } else {
      return ERRSS ("no current graphical element for pledit");
    }
  }
  if (n_element >= 0 || n_contour >= 0) {
    /* retrieve specified element */
    if (n_element >= 0)
      type = GdSetElement (n_element);
    if (n_contour >= 0) {
      if (type != E_CONTOURS)
	return ERRSS ("current graphical element is not contours in pledit");
      type = GdSetContour (n_contour);
    }
    curElement = -6666;		/* differs from -1 to allow pledit after plq */
    curIX = n_element;		/* need these, too */
    curIXc = n_contour;
    if (type == E_LINES) type = 1;
    else if (type == E_DISJOINT) type = 2;
    else if (type == E_TEXT) type = 3;
    else if (type == E_MESH) type = 4;
    else if (type == E_FILLED) type = 5;
    else if (type == E_VECTORS) type = 6;
    else if (type == E_CONTOURS) type = 7;
    else if (type == E_CELLS) type = 8;
    else type = 0;
    if (type == 0)
      return ERRSS ("no such graphical element for pledit");
  }

  BUILD_KWT(kd, editKeys, kwt);
  SETKW(kwt[0],   legend,          setkw_string,   editKeys[0]);
  SETKW(kwt[1],   gistD.hidden,    setkw_boolean,  editKeys[1]);
  if(kwt[2]){
    SETKW(kwt[2], gistA.l.color,   setkw_color,    editKeys[2]);
    gistA.m.color = gistA.f.color = gistA.t.color = gistA.l.color;
  }
  SETKW(kwt[3],   gistA.l.type,    setkw_linetype, editKeys[3]);
  SETKW(kwt[4],   gistA.l.width,   setkw_double,   editKeys[4]);
  SETKW(kwt[5],   gistA.dl.marks,  setkw_boolean,  editKeys[5]);
  SETKW(kwt[6],   gistA.m.color,   setkw_color,    editKeys[6]);
  SETKW(kwt[7],   gistA.m.type,    setkw_xinteger, editKeys[7]);
  SETKW(kwt[8],   gistA.m.size,    setkw_double,   editKeys[8]);
  SETKW(kwt[9],   gistA.dl.mSpace, setkw_double,   editKeys[9]);
  SETKW(kwt[10],  gistA.dl.mPhase, setkw_double,   editKeys[10]);
  SETKW(kwt[11],  gistA.dl.rays,   setkw_boolean,  editKeys[11]);
  SETKW(kwt[12],  gistA.dl.arrowL, setkw_double,   editKeys[12]);
  SETKW(kwt[13],  gistA.dl.arrowW, setkw_double,   editKeys[13]);
  SETKW(kwt[14],  gistA.dl.rSpace, setkw_double,   editKeys[14]);
  SETKW(kwt[15],  gistA.dl.rPhase, setkw_double,   editKeys[15]);
  SETKW(kwt[16],  gistA.dl.closed, setkw_boolean,  editKeys[16]);
  SETKW(kwt[17],  gistA.dl.smooth, setkw_boolean,  editKeys[17]);
  SETKW(kwt[18],  gistA.t.font,    setkw_fonttype, editKeys[18]);
  if(kwt[19]) {
    SETKW(kwt[19], gistA.t.height, setkw_double,   editKeys[19]);
    gistA.t.height *= ONE_POINT;
  }
  SETKW(kwt[20], gistA.t.orient,   setkw_integer,  editKeys[20]);
  if (!gistA.t.orient) {
    gistA.t.orient = TX_RIGHT;
  } else {
    if (gistA.t.orient == 1) gistA.t.orient = TX_UP;
    else if (gistA.t.orient == 2) gistA.t.orient = TX_LEFT;
    else if (gistA.t.orient == 3) gistA.t.orient = TX_DOWN;
    else {
      gistA.t.orient= TX_RIGHT;
      return ERRSS("orient= keyword must be 0, 1, 2, or 3");
    }
  }

  SETKW(kwt[21], dummy,            setkw_justify,  editKeys[21]);
  SETKW(kwt[22], gistA.t.opaque,   setkw_boolean,  editKeys[22]);
  SETKW(kwt[23], gistA.vect.hollow, setkw_boolean,  editKeys[23]);
  SETKW(kwt[24], gistA.vect.aspect, setkw_double,   editKeys[24]);

  if (kwt[25]) {	/* region */
    if (type < 4 || type > 7)
      return ERRSS ("region = in pledit allowed only for plm, plf, plv, plc");
    SETKW(kwt[25],  gistD.region,   setkw_integer,  editKeys[25]);
  }
  if (kwt[26]) {	/* boundary */
    if (type != 4)
      return ERRSS ("boundary = in pledit allowed only for plm");
    SETKW(kwt[26],  gistD.boundary, setkw_boolean,  editKeys[26]);
  }

  if (kwt[27]) {	/* levs */
    double *levels;
    long nLevels = 0;
    PyArrayObject *lap;
    GpReal *lev;
    int i;

    if (type != 7)
      return ERRSS ("levs = in pledit allowed only for plc");

    GET_ARR (lap, kwt[27], Py_GpReal, 1, PyObject *);
    lev = (GpReal *) A_DATA (lap);
    nLevels = A_SIZE (lap);
    if (0 == nLevels) {
      clearArrayList ();
      return ERRSS ("pledit cannot recompute default contour levels");
      }
    levels = malloc (sizeof(GpReal) * nLevels);
    if(!levels) return PyErr_NoMemory();
    for(i = 0; i < nLevels; i++)
      levels[i] = lev[i];
    removeFromArrayList ( (PyObject *) lap);
    /* WARNING --
       this is a critical code section, since until GdEdit successfully
       completes, Gist owns a pointer to the freed levels -- no way to
       gracefully avoid this without "knowing" more about guts of Gist's
       data structures than seem reasonable here... */
    GmFree (gistD.levels);
    gistD.levels = levels;
    gistD.nLevels = nLevels;
    changes |= CHANGE_Z;
  }
  if (kwt[28]) {	/* scale */
    if (type != 6)
      return ERRSS ("scale = in pledit allowed only for plv");
    SETKW(kwt[28],  gistD.scale,  setkw_double,  editKeys[28]);
  }
  if (kwt[29]) {	/* scalem */
    double scalem;
    if (type != 6)
      return ERRSS ("scalem = in pledit allowed only for plv");
    SETKW(kwt[29],  scalem,       setkw_double,  editKeys[29]);
    PyFPE_START_PROTECT("pledit", return 0)
    gistD.scale *= scalem;
    PyFPE_END_PROTECT
  }
  if (kwt[30]) {	/* dx */
    double x0;
    if (type != 3)
      return ERRSS ("dx = in pledit allowed only for plt");
    SETKW(kwt[30],  x0,           setkw_double,  editKeys[30]);
    PyFPE_START_PROTECT("pledit", return 0)
    gistD.x0 += x0;
    PyFPE_END_PROTECT
  }
  if (kwt[31]) {	/* dy */
    double y0;
    if (type != 3)
      return ERRSS ("dy = in pledit allowed only for plt");
    SETKW(kwt[31],  y0,           setkw_double,  editKeys[31]);
    PyFPE_START_PROTECT("pledit", return 0)
    gistD.y0 += y0;
    PyFPE_END_PROTECT
  }
  if (kwt[32]) {
    int edgetype = 0;
    SETKW(kwt[32],  edgetype,     setkw_boolean, editKeys[32]);
    gistA.e.type = edgetype ? L_SOLID : L_NONE;
  }
  SETKW(kwt[33],  gistA.e.color,  setkw_color,   editKeys[33]);
  SETKW(kwt[34],  gistA.e.width,  setkw_double,  editKeys[34]);

  if (kwt[35]) {	/* inhibit */
    if (type != 4)
      return ERRSS ("inhibit = in pledit allowed only for plm");
    SETKW(kwt[35],  gistD.inhibit, setkw_integer, editKeys[35]);
  }
  if (legend) {
    /* Some jiggery-pokery necessary to get the old legend deleted properly,
       and the new legend allocated properly, so that Gist will delete it
       correctly when the graphical element is deleted.  */
    char *oldleg = gistD.legend;
    if (!(gistD.legend = GmMalloc (strlen (legend) + 1)))
       return PyErr_NoMemory();
    strcpy (gistD.legend, legend);
    legend = oldleg;
  }
  GdEdit (changes);
  if (legend)
    GmFree (legend);
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *plf (PyObject * self, PyObject * args, PyObject * kd)
{
  GaQuadMesh mesh;
  PyArrayObject *zap;
  PyObject *zop = 0;
  int nr, nc, convertedZ = 0;
  GpReal *z = 0;
  GpColor *zc = 0, *zc1;
  static char *plfKeys[]= {
    "legend", "hide", "region", "top", "cmin", "cmax",
    "edges", "ecolor", "ewidth", 0 };
  PyObject * kwt[NELT(plfKeys) - 1];
  char *errstr = "plf requires 2D arguments (z [ , y, x, ireg ] )";

  SETJMP0;			/* See Xerror_longjmp() */

  if (PyTuple_Size (args) == 0)
    return ERRSS("plf requires at least one argument");
  BUILD_KWT(kd, plfKeys, kwt);
  TRY (setz_mesh (args, &zop, errstr, 0), (PyObject *) NULL);
  if (!pyMsh.y)
    return ERRSS("No current mesh - set (y, x) first");
  get_mesh (&mesh);

  if (isARRAY(zop) && (A_TYPE(zop) == Py_GpColor)) {
    GET_ARR (zap, zop, Py_GpColor, 2, PyObject *);
    zc = (GpColor *) A_DATA (zap);
  } else {
    if (isARRAY(zop) && (A_TYPE(zop) == Py_GpReal)) {
      GET_ARR (zap, zop, Py_GpReal, 2, PyObject *);
      z = (GpReal *) A_DATA (zap);
    } else {
      z = 0;
      zc = 0;
      zap = 0;
    }
  }

  if (zap) {
    nr = A_DIM(zap, 0);
    nc = A_DIM(zap, 1);
  } else {
    nr = nc = 0;
  }
  if ((z || zc) && ((mesh.iMax != nc   || mesh.jMax != nr) &&
		    (mesh.iMax != nc+1 || mesh.jMax != nr+1))) {
    removeFromArrayList ( (PyObject *) zap);
    return ERRSS (
      "z array must have same or 1 smaller dimensions as mesh in plf");
  }

  TRY (CheckDefaultWindow (), (PyObject *) NULL);
  CheckDefaultPalette ();

  gistD.region = 0;
  SETKW(kwt[2],  gistD.region,    setkw_integer,  plfKeys[2]);

  if (!zc && z) {
    double zmin, zmax, scale, offset;

    TRY (GrabByteScale(&kwt[3], &plfKeys[3], &scale, &offset, &zmin, &zmax,
       z, mesh.reg, gistD.region, mesh.iMax, mesh.jMax,
       mesh.iMax != nc), (PyObject *) NULL);
    TRY (zc = PushColors(z, nc*nr, zmin, zmax, scale, offset),
       (PyObject *) NULL);
    convertedZ= 1;
  }

  GhGetFill();

  SETKW(kwt[0],  gistD.legend,    setkw_string,   plfKeys[0]);
  SETKW(kwt[1],  gistD.hidden,    setkw_boolean,  plfKeys[1]);
  if (kwt[6]) {
    int edgetype = 0;
    SETKW(kwt[6],  edgetype,  setkw_boolean,  plfKeys[6]);
    gistA.e.type = edgetype ? L_SOLID : L_NONE;
  }
  SETKW(kwt[7],  gistA.e.color,   setkw_color,    plfKeys[7]);
  SETKW(kwt[8],  gistA.e.width,    setkw_double, plfKeys[8]);

  zc1 = (mesh.iMax == nc) ? zc + (nc + 1) : zc;
  curElement = -1;
  PyFPE_START_PROTECT("plf", return 0)
  curElement = GdFillMesh(NOCOPY_MESH, &mesh, gistD.region, zc1, nc);
  PyFPE_END_PROTECT
  clearArrayList ();
  if (convertedZ && zc) free (zc);
  if (curElement < 0)
    return ERRSS ("Gist GdFillMesh plotter failed");

  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *plfp (PyObject * self, PyObject * args, PyObject * kd)
{
  PyArrayObject *zap = 0, *yap, *xap, *nap;
  PyObject *zop, *yop, *xop, *nop;
  int i, convertedZ = 0;
  long nz, ny, nx, nn, np, *pn = 0;
  double *z = 0, *x, *y;
  GpColor *zc = 0;
  static char *plfpKeys[]= {
    "legend", "hide", "top", "cmin", "cmax", "edges", "ecolor", "ewidth", 0 };
  PyObject * kwt[NELT(plfpKeys) - 1];
  char *errstr = "plfp requires arguments (z, y, x, n)";

  SETJMP0;			/* See Xerror_longjmp() */

  if (!PyArg_ParseTuple (args, "OOOO", &zop, &yop, &xop, &nop))
    return ERRSS (errstr);

  if (isARRAY(zop) && (A_TYPE(zop) == Py_GpColor)) {
    GET_ARR (zap, zop, Py_GpColor, 1, PyObject *);
    zc = (GpColor *) A_DATA (zap);
  } else if (isARRAY(zop) && (A_TYPE(zop) == Py_GpReal)) {
    GET_ARR (zap, zop, Py_GpReal, 1, PyObject *);
    z = (GpReal *) A_DATA (zap);
  }
  GET_ARR (yap, yop, Py_GpReal, 1, PyObject *);
  GET_ARR (xap, xop, Py_GpReal, 1, PyObject *);
  GET_ARR (nap, nop, PyArray_LONG, 1, PyObject *);
  nn = A_SIZE (nap);
  nx = A_SIZE (xap);
  ny = A_SIZE (yap);
  nz = (zap) ? A_SIZE (zap) : nn;
  y = (GpReal *) A_DATA (yap);
  x = (GpReal *) A_DATA (xap);
  pn = (long *) A_DATA (nap);

  /* Error checking is complicated by required DECREF's on failure. */
  {
    char *es = 0;

    if (nx != ny) es = "len(x) != len(y)";
    if (nz && (nz != nn)) es = "len(n) != len(z)";
    for (np = i = 0; i < nn; i++) np += pn[i];
    if (np != ny) es = "sum(n) != len(y)";
    if (es) {
      clearArrayList ();
      return ERRSS (es);
    }
  }

  BUILD_KWT(kd, plfpKeys, kwt);

  TRY (CheckDefaultWindow (), (PyObject *) NULL);
  CheckDefaultPalette ();

  if (!zc && z) {
    double zmin, zmax, scale, offset;

    TRY (GrabByteScale(&kwt[2], &plfpKeys[2], &scale, &offset, &zmin, &zmax,
       z, (int *)0, 0, nz + 1, 2L, 1), (PyObject *) NULL);
    TRY (zc = PushColors(z, nz, zmin, zmax, scale, offset), (PyObject *) NULL);
    convertedZ= 1;
  }

  GhGetFill();

  SETKW(kwt[0],  gistD.legend,    setkw_string,   plfpKeys[0]);
  SETKW(kwt[1],  gistD.hidden,    setkw_boolean,  plfpKeys[1]);
  if (kwt[5]) {
    int edgetype = 0;
    SETKW(kwt[5],  edgetype,  setkw_boolean,  plfpKeys[5]);
    gistA.e.type = edgetype ? L_SOLID : L_NONE;
  }
  SETKW(kwt[6],  gistA.e.color,   setkw_color,    plfpKeys[6]);
  SETKW(kwt[7],  gistA.e.width,   setkw_double,   plfpKeys[7]);

  curElement = -1;
  PyFPE_START_PROTECT("plfp", return 0)
  curElement = GdFill (nz, zc, x, y, pn);
  PyFPE_END_PROTECT
  clearArrayList ();
  if (convertedZ) free (zc);
  if (curElement < 0)
    return ERRSS ("Gist GdFill plotter failed");

  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *plg (PyObject * self, PyObject * args, PyObject * kd)
{
  PyObject *xop = 0, *yop;
  PyArrayObject *xap, *yap;
  GpReal *x, *y;
  int i;
  long length;
  static char *plgKeys[] = {
    "legend", "hide", "color", "type", "width",
    "marks", "mcolor", "marker", "msize", "mspace", "mphase",
    "rays", "arrowl", "arroww", "rspace", "rphase",
    "closed", "smooth", 0 };
  PyObject * kwt[NELT(plgKeys) - 1];
  char *errstr =
    "plg requires one or two 1-D double arrays, of the same length";

  SETJMP0;			/* See Xerror_longjmp() */

  if (!PyArg_ParseTuple (args, "O|O", &yop, &xop)) {
    return ERRSS (errstr);
  }
  GET_ARR(yap, yop, Py_GpReal, 1, PyObject *);
  length = A_SIZE(yap);
  y = (GpReal *) A_DATA(yap);

  TRY(CheckDefaultWindow (), (PyObject *) NULL);

  GhGetLines();	/* Properties start from defaults for decorated polylines. */
  BUILD_KWT(kd, plgKeys, kwt);
  SETKW(kwt[0],  gistD.legend,    setkw_string,   plgKeys[0]);
  SETKW(kwt[1],  gistD.hidden,    setkw_boolean,  plgKeys[1]);
  SETKW(kwt[2],  gistA.l.color,   setkw_color,    plgKeys[2]);
  SETKW(kwt[2],  gistA.m.color,   setkw_color,    plgKeys[2]);
  SETKW(kwt[3],  gistA.l.type,    setkw_linetype, plgKeys[3]);
  SETKW(kwt[4],  gistA.l.width,   setkw_double,   plgKeys[4]);
  SETKW(kwt[5],  gistA.dl.marks,  setkw_boolean,  plgKeys[5]);
  SETKW(kwt[6],  gistA.m.color,   setkw_color,    plgKeys[6]);
  SETKW(kwt[7],  gistA.m.type,    setkw_xinteger, plgKeys[7]);
  SETKW(kwt[8],  gistA.m.size,    setkw_double,   plgKeys[8]);
  SETKW(kwt[9],  gistA.dl.mSpace, setkw_double,   plgKeys[9]);
  SETKW(kwt[10], gistA.dl.mPhase, setkw_double,   plgKeys[10]);
  SETKW(kwt[11], gistA.dl.rays,   setkw_boolean,  plgKeys[11]);
  SETKW(kwt[12], gistA.dl.arrowL, setkw_double,   plgKeys[12]);
  SETKW(kwt[13], gistA.dl.arrowW, setkw_double,   plgKeys[13]);
  SETKW(kwt[14], gistA.dl.rSpace, setkw_double,   plgKeys[14]);
  SETKW(kwt[15], gistA.dl.rPhase, setkw_double,   plgKeys[15]);
  SETKW(kwt[16], gistA.dl.closed, setkw_boolean,  plgKeys[16]);
  SETKW(kwt[17], gistA.dl.smooth, setkw_boolean,  plgKeys[17]);

  if (xop) {
    GET_ARR(xap, xop, Py_GpReal, 1, PyObject *);
    if(A_SIZE(xap) != length) {
      clearArrayList ();
      return ERRSS (errstr);
    }
    x = (GpReal *) A_DATA(xap);
  } else {
    NEW_MEM (x, length, GpReal, PyObject *);
    for (i = 0; i < length; i++)
	x[i] = (GpReal) (1+i);
  }

  curElement = -1;
  PyFPE_START_PROTECT("plg", return 0)
  curElement = GdLines (length, x, y);
  PyFPE_END_PROTECT

  clearArrayList ();
  clearMemList ();

  if (curElement < 0)
    return ERRSS ("Gist GdLines plotter failed");

  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *pli (PyObject * self, PyObject * args, PyObject * kd)
{
  PyArrayObject *zap;
  PyObject *zop = 0;
  int nr, nc, convertedZ = 0, nargs;
  double *z = 0, x0, y0, x1, y1;
  GpColor *zc = 0;
  static char *pliKeys[]= { "legend", "hide", "top", "cmin", "cmax", 0 };
  PyObject * kwt[NELT(pliKeys) - 1];
  char *errstr = "pli requires arguments (z [ , [ x0, y0, ] x1, y1 ] )";

  SETJMP0;			/* See Xerror_longjmp() */

  switch (nargs = PyTuple_Size (args)) {
  case 5: /* (z, x0, y0, x1, y1) given */
    TRY (PyArg_ParseTuple (args, "Odddd", &zop, &x0, &y0, &x1, &y1), 
       (PyObject *) NULL);
    break;
  case 3: /* (z, x1, y1) given */
    TRY (PyArg_ParseTuple (args, "Odd", &zop, &x1, &y1), 
       (PyObject *) NULL);
    x0 = y0 = 0.0;
    break;
  case 1: /* (z) only given */
    TRY (PyArg_ParseTuple (args, "O", &zop), 
       (PyObject *) NULL);
    x0 = y0 = 0.0;
    break;
  default:
    return ERRSS (errstr);
  }

  if (isARRAY(zop) && (A_TYPE(zop) == Py_GpColor)) {
    GET_ARR (zap, zop, Py_GpColor, 2, PyObject *);
    zc = (GpColor *) A_DATA (zap);
  } else {
    GET_ARR (zap, zop, Py_GpReal, 2, PyObject *);
    z = (GpReal *) A_DATA (zap);
  }
  nr = A_DIM(zap, 0);
  nc = A_DIM(zap, 1);
  if (1 == nargs) {
    x1 = (double) nc;
    y1 = (double) nr;
  }

  BUILD_KWT(kd, pliKeys, kwt);

  TRY (CheckDefaultWindow (), 
     (PyObject *) NULL);
  CheckDefaultPalette ();

  if (!zc) {
    double zmin, zmax, scale, offset;

    TRY (GrabByteScale(&kwt[2], &pliKeys[2], &scale, &offset, &zmin, &zmax,
       z, (int *)0, 0, nc + 1, nr + 1, 1), 
       (PyObject *) NULL);
    TRY (zc = PushColors(z, nc*nr, zmin, zmax, scale, offset), 
       (PyObject *) NULL);
    convertedZ= 1;
  }

  SETKW(kwt[0],  gistD.legend,    setkw_string,   pliKeys[0]);
  SETKW(kwt[1],  gistD.hidden,    setkw_boolean,  pliKeys[1]);

  curElement = -1;
  PyFPE_START_PROTECT("pli", return 0)
  curElement = GdCells (x0, y0, x1, y1, nc, nr, nc, zc);
  PyFPE_END_PROTECT
  removeFromArrayList ( (PyObject *) zap);
  if (convertedZ) free (zc);
  if (curElement < 0)
    return ERRSS ("Gist GdCells plotter failed");

  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *plm (PyObject * self, PyObject * args, PyObject * kd)
{
  GaQuadMesh mesh;
  static char *plmKeys[]= {
    "legend", "hide", "color", "type", "width", "region", "boundary",
    "inhibit", 0 };
  PyObject *kwt[NELT(plmKeys) - 1];
  char *errstr = "plm takes 1-3 non-keyword arguments: (y, x, ireg).";

  SETJMP0;
  if (PyTuple_Size (args) > 0)
    TRY (set_pyMsh (args, errstr, 0), 
       (PyObject *) NULL);

  get_mesh(&mesh);

  GhGetMesh();
  gistD.region = 0;
  gistD.boundary = 0;
  gistD.inhibit = 0;

  BUILD_KWT(kd, plmKeys, kwt);
  SETKW(kwt[0],  gistD.legend,    setkw_string,   plmKeys[0]);
  SETKW(kwt[1],  gistD.hidden,    setkw_boolean,  plmKeys[1]);
  SETKW(kwt[2],  gistA.l.color,   setkw_color,    plmKeys[2]);
  SETKW(kwt[3],  gistA.l.type,    setkw_linetype, plmKeys[3]);
  SETKW(kwt[4],  gistA.l.width,   setkw_double,   plmKeys[4]);
  SETKW(kwt[5],  gistD.region,    setkw_integer,  plmKeys[5]);
  SETKW(kwt[6],  gistD.boundary,  setkw_boolean,  plmKeys[6]);
  SETKW(kwt[7],  gistD.inhibit,   setkw_integer,  plmKeys[7]);

  if (!pyMsh.y)
    return ERRSS("no current mesh - use plmesh(y, x) to initialize");

  TRY(CheckDefaultWindow(), 
    (PyObject *) NULL);
  curElement = -1;
  PyFPE_START_PROTECT("plm", return 0)
  curElement = GdMesh(NOCOPY_MESH, &mesh, gistD.region, gistD.boundary,
		     gistD.inhibit);
  PyFPE_END_PROTECT

  if (curElement < 0)
    return ERRSS ("Gist GdMesh plotter failed");

  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *plmesh (PyObject * self, PyObject * args, PyObject * kd)
{
  static char *meshKeys[] = { "triangle", 0 };
  PyObject * kwt[NELT(meshKeys) - 1];
  char *errstr = "plmesh takes 0-3 non-keyword arguments: (y, x, ireg).";

  BUILD_KWT(kd, meshKeys, kwt);
  TRY (set_pyMsh (args, errstr, kwt[0]), 
     (PyObject *) NULL);

  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *plq (PyObject * self, PyObject * args)
{
  int type, n_element = 0, n_contour = 0;
  extern int puts (const char *);

  switch (PyTuple_Size (args)) {
  case 2: /* (n_element, n_contour) given */
    TRY (PyArg_ParseTuple (args, "ii", &n_element, &n_contour), 
       (PyObject *) NULL);
    break;
  case 1: /* (n_element) given */
    TRY (PyArg_ParseTuple (args, "i", &n_element), 
       (PyObject *) NULL);
    break;
  case 0: /* () given */
    break;
  default:
    return ERRSS ("plq function takes no more than two arguments");
  }

  /* Pygist uses 1-origin element numbering, Gist uses 0-origin */
  n_element--;
  n_contour--;

  if (n_element >= 0) {
    /* retrieve specified element */
    type = GdSetElement (n_element);
    if (n_contour >= 0) {
      if (type != E_CONTOURS)
	return ERRSS ("current graphical element is not contours in pledit");
      type = GdSetContour (n_contour);
    }
    curElement = -6666;		/* differs from -1 to allow pledit after plq */
    curIX = n_element;		/* need these, too */
    curIXc = n_contour;
    if (type == E_LINES) type = 1;
    else if (type == E_DISJOINT) type = 2;
    else if (type == E_TEXT) type = 3;
    else if (type == E_MESH) type = 4;
    else if (type == E_FILLED) type = 5;
    else if (type == E_VECTORS) type = 6;
    else if (type == E_CONTOURS) type = 7;
    else if (type == E_CELLS) type = 8;
    else if (type == E_POLYS) type = 9;
    else type = 0;

    if (1) {
      /* return printed summary of keyword values */
      char line[120];
      PrintInit (puts);

      if (type == 0) {
	sprintf (line, "<no such object>  element# %d", n_element + 1);
	PrintFunc (line);
	if (n_contour >= 0) {
	  sprintf (line, "  contour# %d", n_contour + 1);
	  PrintFunc (line);
	}
	ForceNewline ();

      } else if (type == 1) {
	sprintf (line, "plg  element# %d", n_element + 1);
	PrintFunc (line);
	if (n_contour >= 0) {
	  sprintf (line, "  contour# %d", n_contour + 1);
	  PrintFunc (line);
	  ForceNewline ();
	  sprintf (line, "  at level value %g", gistD.levels[n_contour]);
	  PrintFunc (line);
	}
	ForceNewline ();
	PrintHideLegend (line, type);
	PrintColor (line, gistA.l.color, 1);
	PrintTypeWidth (line, 3);
	PrintMarks (line, 3);
	sprintf (line, "rays= %d,", gistA.dl.rays);
	PrintFunc (line);
	ForceNewline ();
	sprintf (line,
		 "  arrowl= %.2f, arroww= %.2f, rspace= %.5f, rphase= %.5f,",
		 Safe_dbl (gistA.dl.arrowL), Safe_dbl (gistA.dl.arrowW),
		 Safe_dbl (gistA.dl.rSpace), Safe_dbl (gistA.dl.rPhase));
	PrintFunc (line);
	ForceNewline ();
	sprintf (line, "smooth= %d,  closed= %d",
		 gistA.dl.smooth, gistA.dl.closed);
	PrintFunc (line);
	ForceNewline ();

      } else if (type == 2) {
	sprintf (line, "pldj  element# %d", n_element + 1);
	PrintFunc (line);
	ForceNewline ();
	PrintHideLegend (line, type);
	PrintColor (line, gistA.l.color, 1);
	PrintTypeWidth (line, 2);

      } else if (type == 3) {
	sprintf (line, "plt  element# %d", n_element + 1);
	PrintFunc (line);
	ForceNewline ();
	PrintHideLegend (line, type);
	PrintColor (line, gistA.t.color, 3);
	sprintf (line, "text= %.80s", gistD.text);
	PrintFunc (line);
	ForceNewline ();

      } else if (type == 4) {
	sprintf (line, "plm  element# %d", n_element + 1);
	PrintFunc (line);
	ForceNewline ();
	PrintHideLegend (line, type);
	PrintColor (line, gistA.l.color, 1);
	PrintTypeWidth (line, 2);
	PrintRegion (line, 1);
	sprintf (line, "boundary= %d, inhibit= %d", gistD.boundary,
		 gistD.inhibit);
	PrintFunc (line);
	ForceNewline ();

      } else if (type == 5) {
	sprintf (line, "plf  element# %d", n_element + 1);
	PrintFunc (line);
	ForceNewline ();
	PrintHideLegend (line, type);
	sprintf (line, "edges= %d, e", gistA.e.type != L_NONE);
	PrintFunc (line);
	PrintColor (line, gistA.e.color, 1);
	sprintf (line, "ewidth= %.2f", Safe_dbl (gistA.e.width));
	PrintFunc (line);
	ForceNewline ();
	PrintRegion (line, 2);

      } else if (type == 6) {
	sprintf (line, "plv  element# %d", n_element + 1);
	PrintFunc (line);
	ForceNewline ();
	PrintHideLegend (line, type);
	PrintColor (line, gistA.l.color, 1);
	sprintf (line, "width= %.2f,", Safe_dbl (gistA.l.width));
	PrintFunc (line);
	ForceNewline ();
	sprintf (line, "hollow= %d,  aspect= %.4f,", gistA.vect.hollow,
		 Safe_dbl (gistA.vect.aspect));
	PrintFunc (line);
	ForceNewline ();
	PrintRegion (line, 3);
	sprintf (line, "scale= %g", gistD.scale);
	PrintFunc (line);
	ForceNewline ();

      } else if (type == 7) {
	int i;
	sprintf (line, "plc  element# %d", n_element + 1);
	PrintFunc (line);
	ForceNewline ();
	PrintHideLegend (line, type);
	PrintColor (line, gistA.l.color, 1);
	PrintTypeWidth (line, 3);
	PrintMarks (line, 3);
	sprintf (line, "smooth= %d,", gistA.dl.smooth);
	PrintFunc (line);
	ForceNewline ();
	PrintRegion (line, 2);
	sprintf (line, "%d contour levels, levs=", gistD.nLevels);
	PrintFunc (line);
	ForceNewline ();
	PrintFunc ("[");
	if (gistD.nLevels > 0) {
	  for (i = 0;; i++) {
	    sprintf (line, "%g", gistD.levels[i]);
	    PrintFunc (line);
	    if (i == gistD.nLevels - 1)
	      break;
	    PrintFunc (",");
	    PermitNewline (0);
	  }
	}
	PrintFunc ("]");
	ForceNewline ();

      } else if (type == 8) {
	sprintf (line, "pli  element# %d", n_element + 1);
	PrintFunc (line);
	ForceNewline ();
	PrintHideLegend (line, type);
	sprintf (line, "x0= %g,  y0= %g,  x1= %g,  y1= %g",
		 gistD.px, gistD.py, gistD.qx, gistD.qy);
	PrintFunc (line);
	ForceNewline ();

      } else if (type == 9) {
	sprintf (line, "plfp  element# %d", n_element + 1);
	PrintFunc (line);
	ForceNewline ();
	PrintHideLegend (line, type);
	sprintf (line, "%d polygons", gistD.n);
	PrintFunc (line);
	ForceNewline ();
      }
    } else {
      /* Future implementation of array return */
    }

  } else if (n_contour >= 0) {
    return ERRSS ("contour number cannot be specified without element number");

  } else {
    char line[16];
    int i, offset;
    /* print list of legends... */
    /* if (CalledAsSubroutine()) FUTURE */
    if (1) {
      /* ...to terminal */
      PrintInit (puts);
    } else {
     /* Future Implementation */
    }

    curElement = -1;
    for (i = 0; (type = GdSetElement (i)) != E_NONE; i++) {
      sprintf (line, "%s%2d: ", gistD.hidden ? "(H)" : "", i + 1);
      PrintFunc (line);
      offset = 0;
      if ((type == E_LINES || type == E_CONTOURS) && gistD.legend &&
	  gistD.legend[0] == '\001') {
	char marker[2];
	marker[1] = '\0';
	if (gistA.m.type >= ' ' && gistA.m.type < '\177')
	  marker[0] = (char) gistA.m.type;
	else if (gistA.m.type >= 1 && gistA.m.type <= 5)
	  marker[0] = specialMarkers[gistA.m.type - 1];
	else
	  marker[0] = '?';
	PrintFunc (marker);
	offset = 1;
      }
      if (gistD.legend)
	PrintFunc (gistD.legend + offset);
      ForceNewline ();
    }
  }
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *plsys (PyObject * self, PyObject * args)
{
  int n = -9999, n0;
  char *errstr = "Error: plsys takes zero or one integer argument.";

  SETJMP0;
  if (!PyArg_ParseTuple (args, "|i", &n)) {
    return ERRSS (errstr);
  }

  CheckDefaultWindow ();
  n0 = GdGetSystem();

  if (n != -9999){
    if (GdSetSystem (n) != E_SYSTEM && n != 0) {
      return ERRSS (
       "No such coordinate system exists in current graphics window.");
    }
  }
  return Py_BuildValue ("i",n0);
}

static PyObject *plt (PyObject * self, PyObject * args, PyObject * kd)
{
  char *text = 0;
  double x = 0.0, y = 0.0;
  int toSys = 0, dummy;
  static char *pltKeys[] = {
    "legend", "hide",
    "color", "font", "height", "orient", "justify", "opaque", "tosys", 0 };
  PyObject *kwt[NELT (pltKeys) - 1];

  SETJMP0;			/* See Xerror_longjmp() */

  if (!PyArg_ParseTuple (args, "sdd", &text, &x, &y)) {
    return ERRSS ("plt requires exactly three non-keyword arguments");
  }

  CheckDefaultWindow ();

  /* set properties, starting from defaults for vectors */
  GhGetText ();
  BUILD_KWT(kd, pltKeys, kwt);
  SETKW(kwt[0],  gistD.legend,    setkw_string,   pltKeys[0]);
  SETKW(kwt[1],  gistD.hidden,    setkw_boolean,  pltKeys[1]);
  SETKW(kwt[2],  gistA.t.color,   setkw_color,    pltKeys[2]);
  SETKW(kwt[3],  gistA.t.font,    setkw_fonttype, pltKeys[3]);
  SETKW(kwt[4],  gistA.t.height,  setkw_double,   pltKeys[4]);
  if(kwt[4])  gistA.t.height *= ONE_POINT;
  SETKW(kwt[5],  gistA.t.orient,  setkw_integer,  pltKeys[5]);
  if (!gistA.t.orient) {
    gistA.t.orient = TX_RIGHT;
  } else {
    if (gistA.t.orient == 1) gistA.t.orient = TX_UP;
    else if (gistA.t.orient == 2) gistA.t.orient = TX_LEFT;
    else if (gistA.t.orient == 3) gistA.t.orient = TX_DOWN;
    else {
      gistA.t.orient= TX_RIGHT;
      return ERRSS("orient= keyword must be 0, 1, 2, or 3");
    }
  }
  SETKW(kwt[6],  dummy,           setkw_justify,  pltKeys[6]);
  SETKW(kwt[7],  gistA.t.opaque,  setkw_boolean,  pltKeys[7]);
  SETKW(kwt[8],  toSys,           setkw_boolean,  pltKeys[8]);


  if (!text)
    text = "";
  curElement = -1;
  PyFPE_START_PROTECT("plt", return 0)
  curElement = GdText (x, y, text, toSys);
  PyFPE_END_PROTECT
  if (curElement < 0)
    return ERRSS ("Gist GdText plotter failed");

  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *plv (PyObject * self, PyObject * args, PyObject * kd)
{
  GaQuadMesh mesh;
  PyArrayObject *uap, *vap;
  PyObject *uop, *vop;
  int nr, nc;
  GpReal *u = 0, *v = 0, scale;
  static char *plvKeys[] = {
    "legend", "hide", "region",
    "color", "hollow", "width", "aspect", "scale", 0 };
  PyObject * kwt[NELT(plvKeys) - 1];
  char *errstr =
    "plv requires 2D arguments (v, u [ , y, x, ireg, scale = dt ] )";

  SETJMP0;			/* See Xerror_longjmp() */
  TRY(CheckDefaultWindow (), 
     (PyObject *) NULL);

  if (PyTuple_Size (args) < 2)
    return ERRSS("plv requires at least two arguments");
  BUILD_KWT(kd, plvKeys, kwt);
  TRY (setvu_mesh (args, &vop, &uop, errstr), 
     (PyObject *) NULL);
  if (!pyMsh.y)
    return ERRSS("No current mesh - set (y, x) first");
  GET_ARR (vap, vop, Py_GpReal, 2, PyObject *);
  GET_ARR (uap, uop, Py_GpReal, 2, PyObject *);
  nr = (A_DIM(vap, 0) == A_DIM(uap, 0)) ? A_DIM(vap, 0) : 0;
  nc = (A_DIM(vap, 1) == A_DIM(uap, 1)) ? A_DIM(vap, 1) : 0;
  if (A_DIM (pyMsh.y, 0) != nr || A_DIM (pyMsh.y, 1) != nc) {
    clearArrayList ();
    return ERRSS("(v, u) arrays must match (y, x) mesh arrays in shape");
  }
  v = (GpReal *) A_DATA (vap);
  u = (GpReal *) A_DATA (uap);
  get_mesh (&mesh);

  GhGetVectors();	/* Properties start from defaults for vectors. */
  gistD.region = 0;

  SETKW(kwt[0],   gistD.legend,      setkw_string,  plvKeys[0]);
  SETKW(kwt[1],   gistD.hidden,      setkw_boolean, plvKeys[1]);
  SETKW(kwt[2],   gistD.region,      setkw_integer, plvKeys[2]);
  SETKW(kwt[3],   gistA.l.color,     setkw_color,   plvKeys[3]);
  if(kwt[3])      gistA.f.color = gistA.l.color;
  SETKW(kwt[4],   gistA.vect.hollow, setkw_boolean, plvKeys[4]);

  SETKW(kwt[5],   gistA.l.width,     setkw_double,  plvKeys[5]);
  SETKW(kwt[6],   gistA.vect.aspect, setkw_double,  plvKeys[6]);
  if(kwt[7]) { /* scale= keyword */
    SETKW(kwt[7], scale,             setkw_double,  plvKeys[7]);
  } else {
    /* set vector scale factor to make maximum vector length a
       "typical" zone dimension */
    double umin, umax, vmin, vmax, xmin, xmax, ymin, ymax;

    PyFPE_START_PROTECT("plv", return 0)
    GetPCrange(&xmin, &xmax, mesh.x, mesh.reg, gistD.region, nc, nr);
    GetPCrange(&ymin, &ymax, mesh.y, mesh.reg, gistD.region, nc, nr);
    GetPCrange(&umin, &umax, u, mesh.reg, gistD.region, nc, nr);
    GetPCrange(&vmin, &vmax, v, mesh.reg, gistD.region, nc, nr);

    umax -= umin;
    vmax -= vmin;
    if (vmax > umax) umax = vmax;
    xmax = (xmax - xmin) + (ymax - ymin);
    xmax /= (nc + nr);

    if (umax > 0.0) scale = xmax / umax;
    else scale = 1.0;
    PyFPE_END_PROTECT
  }

  curElement = -1;
  PyFPE_START_PROTECT("plv", return 0)
  curElement = GdVectors(NOCOPY_MESH, &mesh, gistD.region, u, v, scale);
  PyFPE_END_PROTECT
  clearArrayList ();
  if (curElement < 0)
    return ERRSS ("Gist GdVectors plotter failed");

  Py_INCREF (Py_None);
  return Py_None;
}

#if 0
static void print_array_stats(PyArrayObject *op)
{
  int i,ne;
  double *dp;
  fprintf(stderr,"Data pointer: %p Base pointer: %p\n", op->data, op->base);
  fprintf(stderr,"Num dims: %d Flags: %d\n", op->nd, op->flags);
  fprintf(stderr,"Dims & strides:");
  for(i=0; i<op->nd; i++)
    fprintf(stderr," i: %d dim: %d stride: %d",i,op->dimensions[i], op->strides[i]);
  fprintf(stderr,"\n");
  ne = op->dimensions[0];
  for(i=1; i<op->nd; i++)
    ne *= op->dimensions[i];
  fprintf(stderr,"Data: (ne = %d)", ne);
  for(i=0,dp = (double *)op->data; i < ne; i++, dp++)
    fprintf(stderr," %.1g", *dp);
  fprintf(stderr,"\n\n");
}
#endif

static PyObject *redraw (PyObject * self, PyObject * args)
{
  SETJMP0;
  PyFPE_START_PROTECT("redraw", return 0)
  CheckDefaultWindow ();
  GhRedraw ();
  PyFPE_END_PROTECT
  Py_INCREF (Py_None);
  return Py_None;
}

/* Create a default region array for the current mesh. */
static int set_def_reg (int nr, int nc)
{
  int i, ne, newlen, *p1;
  PyArrayObject *ra1;

  ne = nr * nc;
  newlen = ne + nc + 1;
  TRY (ra1 = (PyArrayObject *) PyArray_FromDims (1, &newlen, PyArray_INT), 0);
  p1 = (int *) A_DATA (ra1);

  /* Fill in the data part of the new region array. */
  for (i = 0;      i <= nc; i++)    p1[i]      = 0;
  for (i = nc+1;   i <  ne; i++)    p1[i]      = 1;
  for (i = 0;      i <  nc; i++)    p1[ne + i] = 0;
  for (i = 2*nc;   i <  ne; i+=nc)  p1[i]      = 0;

  Py_XDECREF (pyMsh.reg);
  pyMsh.reg = ra1;
  return 1;
}

/* Set xmin, xmax, ymin, or ymax.  Return 2 if limit was unchanged
 * (e.g., user passed 'u' argument), return 1 if limit was reset,
 * return 0 if error. */
static int set_limit (PyObject * ob, double *lim, int *flags, int fval)
{
  if (PyString_Check (ob)) {
    char *s = PyString_AsString (ob);
    if (*s == 'e' || *s == 'E') {
      *flags |= fval;
    } else if (*s == 'u' || *s == 'U') {
      return 2; /* No change in this limit. */
    } else {
      return 0; /* Error */
    }
  } else if (PyFloat_Check (ob)) {
    *lim = PyFloat_AsDouble (ob);
    *flags &= ~fval;
  } else if (PyInt_Check (ob)) {
    *lim = (double) PyInt_AsLong (ob);
    *flags &= ~fval;
  } else {
    return 0; /* Error */
  }
  return 1; /* Limit was changed, either to specific value or extreme. */
}

/* Called from plmesh and plm, and indirectly from plc, plf, and plv,
 * which all allow the arguments (y, x, ireg).  Plc and plmesh also
 * allow the "triangle=" keyword, so handle that also.
 */
static int set_pyMsh(PyObject *args, char *errstr, PyObject *tri)
{
  PyObject *op1, *op2, *op3;

  if (!PyArg_ParseTuple (args, "|OOO", &op1, &op2, &op3))
    return (int) ERRSS (errstr);

  switch (PyTuple_Size(args)) {
  case 3: /* Arguments were (y, x, ireg). */
    TRY (set_yx (op1, op2), 0);
    TRY (set_reg (op3), 0);
    break;

  case 2: /* Arguments were (y, x). */
    TRY(set_yx(op1, op2), 0);
    TRY(set_def_reg(A_DIM(op1, 0), A_DIM(op1, 1)), 0); /* Default region array. */
    break;

  case 1: /* Arguments were (ireg). */
    TRY(set_reg(op1), 0);
    break;

  case 0: /* Clear the mesh, unless "triangle" keyword was given. */
    if(tri) TRY (set_tri(tri), 0);
    else clear_pyMsh();
    break;

  default: /* This **REALLY** shouldn't happen. */
    return (int) ERRSS (errstr);
  }
  return 1;
}

/* Create a non-default region (mesh) array from the passed-in object. */
static int set_reg (PyObject *op)
{
  int i, ok, nr, nc, ne, newlen, *p2, *p1;
  PyArrayObject *ra2, *ra1;

  ok = (isARRAY(op) && (A_NDIM(op) == 2) && ((A_TYPE(op) == PyArray_INT)
     || (A_TYPE(op) == PyArray_LONG)));
  if (!ok) return (int) ERRSS("(ireg) must be a 2-D int array");

  if (!pyMsh.y)
    return (int) ERRSS("No current mesh - ireg not set - set (y, x) first");
  nr = A_DIM (op, 0);
  nc = A_DIM (op, 1);
  if (A_DIM (pyMsh.y, 0) != nr || A_DIM (pyMsh.y, 1) != nc)
    return (int) ERRSS("(ireg) must match (y, x) in shape");

  ne = nr * nc;
  newlen = ne + nc + 1;
  NEW_ARR (ra1, 1, &newlen, PyArray_INT, int);
  p1 = (int *) A_DATA (ra1);
  GET_ARR (ra2, op, PyArray_INT, 2, int);
  p2 = (int *) A_DATA (ra2);

  /* Fill in the data part of the new region array. */
  for (i = 0;      i <= nc; i++)    p1[i]      = 0;
  for (i = nc+1;   i <  ne; i++)    p1[i]      = p2[i];
  for (i = 0;      i <  nc; i++)    p1[ne + i] = 0;
  for (i = 2*nc;   i <  ne; i+=nc)  p1[i]      = 0;

  Py_DECREF (ra2);
  Py_XDECREF (pyMsh.reg);
  array_list_length = 0;
  pyMsh.reg = ra1;
  takeOffArrayList ( (PyObject *) ra1);
  takeOffArrayList ( (PyObject *) ra2);
  return 1;
}

static PyObject *set_slice2_precision (PyObject * self, PyObject * args)
{
 if ( ! PyArg_ParseTuple (args, "d", &_slice2_precision))
    return (PyObject *) ERRSS ("set_slice2_precision: bad value.");
 Py_INCREF (Py_None);
 return Py_None;
}

static PyObject *setdpi (PyObject * self, PyObject * args)
{
  int temp_dpi;
  if ( ! PyArg_ParseTuple (args, "i", &temp_dpi))
    return (PyObject *) ERRSS ("set_default_dpi: bad value.");
  if ((temp_dpi != 75) && (temp_dpi != 100))
    return (PyObject *) ERRSS ("set_default_dpi: value must be 75 or 100.");
  defaultDPI = temp_dpi;
  Py_INCREF (Py_None);
  return Py_None;
}

/* Create a triangulation (mesh) array. */
static int set_tri (PyObject *top)
{
  int nr, nc;

  if (!pyMsh.y)
    return (int) ERRSS("No current mesh - triangle not set - set (y, x) first");
  nr = A_DIM (pyMsh.y, 0);
  nc = A_DIM (pyMsh.y, 1);

  Py_XDECREF (pyMsh.triangle);
  GET_ARR (pyMsh.triangle, top, PyArray_SHORT, 2, int);

  if (A_DIM (pyMsh.triangle, 0) != nr || A_DIM (pyMsh.triangle, 1) != nc) {
    removeFromArrayList ((PyObject *)pyMsh.triangle);
    return (int) ERRSS("triangle array must match shape of (y, x).");
  }
  array_list_length = 0;
  return 1;
}

static int set_yx (PyObject *yop, PyObject *xop)
{
  int nr, nc;

  clear_pyMsh();
  GET_ARR (pyMsh.y, yop, Py_GpReal, 2, int);
  nr = A_DIM (pyMsh.y, 0);
  nc = A_DIM (pyMsh.y, 1);
  if (nr < 2 || nc < 2) {
    clearArrayList ();
    return (int) ERRSS("(y, x) arrays must be at least 2X2");
  }
  GET_ARR (pyMsh.x, xop, Py_GpReal, 2, int);
  if (A_DIM (pyMsh.x, 0) != nr || A_DIM (pyMsh.x, 1) != nc) {
    clearArrayList ();
    return (int) ERRSS("x array must match shape of y");
  }

#if 0
print_array_stats(xop);
print_array_stats(yop);
#endif

  array_list_length = 0;
  return 1;
}

static int setkw_boolean (PyObject * v, int *t, char *kw)
{
  char buf[256];
  char *format = "%s keyword requires argument of 0 (False) or 1 (True)";

  if (PyInt_Check (v)) {
    *t = (PyInt_AsLong (v) != 0);
    return 1;
  }
  (void)sprintf(buf, format, kw);
  return (int) ERRSS (buf);
}

/* Set value for "color=" keyword.  Value passed can be either a string
 * or an integer.  All these setkw_*() functions return 0 on error,
 * non-zero otherwise. */
static int setkw_color (PyObject * v, int *t, char *kw)
{
  int color = FG_COLOR;

  if (PyString_Check (v)) {
    char *s = PyString_AsString (v);
    if (strcmp (s, "bg") == 0)
      color = BG_COLOR;
    else if (strcmp (s, "fg") == 0)
      color = FG_COLOR;
    else if (strcmp (s, "black") == 0)
      color = BLACK_COLOR;
    else if (strcmp (s, "white") == 0)
      color = WHITE_COLOR;
    else if (strcmp (s, "red") == 0)
      color = RED_COLOR;
    else if (strcmp (s, "green") == 0)
      color = GREEN_COLOR;
    else if (strcmp (s, "blue") == 0)
      color = BLUE_COLOR;
    else if (strcmp (s, "cyan") == 0)
      color = CYAN_COLOR;
    else if (strcmp (s, "magenta") == 0)
      color = MAGENTA_COLOR;
    else if (strcmp (s, "yellow") == 0)
      color = YELLOW_COLOR;
    else {
      char errstr[256];
      sprintf (errstr, "Unrecognized color keyword: %s: "
	       "Use fg, bg, or 8 primaries only", s);
      return (int) ERRSS (errstr);
    }
  } else if (PyInt_Check (v)) {
    color = PyInt_AsLong (v);
  } else {
    return (int) ERRSS ("Color keyword value must be string or integer");
  }
  *t = color;
  return 1;
}

static int setkw_double (PyObject * v, double *t, char *kw)
{
  char buf[256];
  char *format = "%s keyword requires floating point argument";

  if (PyFloat_Check (v)) {
    *t = PyFloat_AsDouble (v);
    return 1;
  } else if (PyInt_Check (v)) {
    *t = (double) PyInt_AsLong (v);
    return 1;
  }
  (void)sprintf(buf, format, kw);
  return (int) ERRSS (buf);
}

static int setkw_fonttype (PyObject * v, int *t, char *kw)
{
  char buf[256];
  char *format = "%s keyword requires string or integer argument";
  char *s;
  int font, face;

#undef SETFONT
#define SETFONT(s, n, font, face, type) \
if (GetTypeface (&s[n], &face)) font = type|face; else return 0

  if (PyString_Check (v)) {
    s = PyString_AsString (v);

    if (strncmp (s, "courier", 7) == 0) {
      SETFONT (s, 7, font, face, T_COURIER);
    } else if (strncmp (s, "times", 5) == 0) {
      SETFONT (s, 5, font, face, T_TIMES);
    } else if (strncmp (s, "helvetica", 9) == 0) {
      SETFONT (s, 9, font, face, T_HELVETICA);
    } else if (strncmp (s, "symbol", 6) == 0) {
      SETFONT (s, 6, font, face, T_SYMBOL);
    } else if (strncmp (s, "schoolbook", 10) == 0) {
      SETFONT (s, 10, font, face, T_NEWCENTURY);
    } else {
      return (int) ERRSS ("unrecognized font keyword");
    }
  } else if (PyInt_Check (v)) {
    font = PyInt_AsLong (v);
  } else {
    (void)sprintf(buf, format, kw);
    return (int) ERRSS (buf);
  }
  *t = font;
  return 1;
#undef SETFONT
}

static int setkw_integer (PyObject * v, int *t, char *kw)
{
  char buf[256];
  char *format = "%s keyword requires integer argument";

  if (PyInt_Check (v)) {
    *t = PyInt_AsLong (v);
    return 1;
  }
  (void)sprintf(buf, format, kw);
  return (int) ERRSS (buf);
}

static int setkw_justify (PyObject * v, int *t, char *kw)
{
  char *s;
  char buf[256];
  char *format = "%s keyword requires string or integer argument";

  if (PyString_Check (v)) {
    s = PyString_AsString (v);

    /* Inspect the first character. */
    if (*s == 'N') {
      gistA.t.alignH = TH_NORMAL;
      s++;
    } else if (*s == 'L') {
      gistA.t.alignH = TH_LEFT;
      s++;
    } else if (*s == 'C') {
      gistA.t.alignH = TH_CENTER;
      s++;
    } else if (*s == 'R') {
      gistA.t.alignH = TH_RIGHT;
      s++;
    } else {
      while (*s)
	s++;
    }

    /* Inspect the second character. */
    if (*s == 'N') {
      gistA.t.alignV = TV_NORMAL;
    } else if (*s == 'T') {
      gistA.t.alignV = TV_TOP;
    } else if (*s == 'C') {
      gistA.t.alignV = TV_CAP;
    } else if (*s == 'H') {
      gistA.t.alignV = TV_HALF;
    } else if (*s == 'A') {
      gistA.t.alignV = TV_BASE;
    } else if (*s == 'B') {
      gistA.t.alignV = TV_BOTTOM;
    } else {
      return (int) ERRSS ("unrecognized justify keyword");
    }

  } else if (PyInt_Check (v)) {
    int justify = PyInt_AsLong (v);
    gistA.t.alignH = justify & 3;
    gistA.t.alignV = justify >> 2;

  } else {
    (void)sprintf(buf, format, kw);
    return (int) ERRSS (buf);
  }
  return 1;
}

static int setkw_linetype (PyObject * v, int *t, char *kw)
{
  int type = 0;
  char *errstr =
    "bad linetype: Use \"none\", \"solid\", \"dash\", \"dot\",\n"
    "\"dashdot\", \"dashdotdot\", or 0-5, respectively.";

  if (PyString_Check (v)) {
    char *s = PyString_AsString (v);
    if (strcmp (s, "none") == 0)
      type = L_NONE;
    else if (strcmp (s, "solid") == 0)
      type = L_SOLID;
    else if (strcmp (s, "dash") == 0)
      type = L_DASH;
    else if (strcmp (s, "dot") == 0)
      type = L_DOT;
    else if (strcmp (s, "dashdot") == 0)
      type = L_DASHDOT;
    else if (strcmp (s, "dashdotdot") == 0)
      type = L_DASHDOTDOT;
    else
      return (int) ERRSS (errstr);
  } else if (PyInt_Check (v)) {
    type = PyInt_AsLong (v);
  } else {
    return (int) ERRSS (errstr);
  }
  if (type < 0)
    type = 0;
  else if (type > 5)
    type = 1 + (type - 1) % 5;

  *t = type;
  return 1;
}

static int setkw_string (PyObject * v, char **t, char *kw)
{
  char buf[256];
  char *format = "%s keyword requires string argument";

  if (PyString_Check (v)) {
    *t = PyString_AsString (v);
    /* Should I Py_INCREF(v)?  PyString_AsString() just returns the pointer
       to its internal string value, which I'm saving away.  Hopefully,
       Gist strcpy's the string into its own space, but I dunno.... */
    return 1;
  }
  (void)sprintf(buf, format, kw);
  return (int) ERRSS (buf);
}

/* "Extended" integer: Allow the user to input a character. */
static int setkw_xinteger (PyObject * v, int *t, char *kw)
{
  char buf[256];
  char *format = "%s keyword requires integer or single character argument";

  if (PyInt_Check (v)) {
    *t = PyInt_AsLong (v);
    return 1;
  } else if (PyString_Check (v)) {
    char *s = PyString_AsString (v);
    *t = (int) s[0];
    return 1;
  }
  (void)sprintf(buf, format, kw);
  return (int) ERRSS (buf);
}

/* Set v, u, and the (y, x, ireg) mesh variables.
 * Called from plv and from mesh_loc, which happens to take the same args.
 * Note that PyObject references returned in vop and uop are
 * borrowed, so should usually NOT be DECREF'ed.
 * Returns 0 on failure, 1 otherwise.
 */
static int setvu_mesh(
  PyObject *args, PyObject **vop, PyObject **uop, char *errstr)
{
  int n;
  PyObject *newargs;

  switch (n = PyTuple_Size(args)) {
  case 5: /* (v, u, y, x, ireg) given */
  case 4: /* (v, u, y, x) given */
    TRY (newargs = PyTuple_GetSlice (args, 2, n), 0);
    TRY (set_pyMsh (newargs, errstr, 0), 0);
    Py_DECREF (newargs);
    /* (Fall through.) */
  case 2: /* (v, u) only given */
    TRY ( *vop = PyTuple_GetItem (args, 0), 0); /* Borrowed reference returned */
    TRY ( *uop = PyTuple_GetItem (args, 1), 0);
    break;
  default:
    return (int) ERRSS(errstr);
  }
  return 1;
}

/* Set z and the (y, x, ireg) mesh variables.  Called from plc and plf.
 * Note that PyObject reference returned in zop is borrowed, so should
 * usually NOT be DECREF'ed.
 * Returns 0 on failure, 1 otherwise.
 */
static int setz_mesh (
  PyObject *args, PyObject **zop, char *errstr, PyObject *tri)
{
  int n;
  PyObject *newargs;

  switch (n = PyTuple_Size(args)) {
  case 4: /* (z, y, x, ireg) given */
  case 3: /* (z, y, x) given */
    TRY (newargs = PyTuple_GetSlice (args, 1, n), 0);
    TRY (set_pyMsh (newargs, errstr, tri), 0);
    Py_DECREF (newargs);
    /* (Fall through.) */
  case 1: /* (z) only given */
    TRY ( *zop = PyTuple_GetItem (args, 0), 0);
    break;
  default:
    return (int) ERRSS(errstr);
  }
  return 1;
}

/* Used only by limits() */
static int unpack_limit_tuple (PyObject * ob, double limits[], int *flags)
{
  int i, size = PyTuple_Size (ob);
  PyObject *item;
  if (5 != size) {
    return (int) ERRSS ("Old limits must have four doubles + 1 integer");
  }
  for (i = 0; i < 4; i++) {
    if ((item = PyTuple_GetItem (ob, i)) == 0) {
      return (int) ERRSS ("Error unpacking limit tuple.");
    }
    if (PyFloat_Check (item)) {
      limits[i] = PyFloat_AsDouble (item);
    } else if (PyInt_Check (item)) {
      limits[i] = (double) PyInt_AsLong (item);
    } else {
      return (int) ERRSS ("Expected floating point value");
    }
  }
  if ((item = PyTuple_GetItem (ob, 4)) == 0) {
    return (int) ERRSS ("Error unpacking flags in limit tuple.");
  }
  if (PyInt_Check (item)) {
    *flags = (int) PyInt_AsLong (item);
  } else {
    return (int) ERRSS ("Expected integer value");
  }
  return 1;
}

static PyObject * slice2 (PyObject * self, PyObject * args)
/* slice2 (plane, nverts, xyzverts, values = None, _slice2x = 0)
   returns either a triple [nverts, xyzverts, values] or a sextuple
   [nverts, xyzverts, values, nvertb, xyzvertb, valueb]            */
{
 PyObject * oplane,
          * onverts,
          * oxyzverts,
          * ovalues = (PyObject *) NULL,
          * oreturn_value,
          /* return values */
          * ornverts,
          * orxyzverts,
          * orvalues = (PyObject *) NULL,
          /* additional slice2x return values */
          * ornvertb,
          * orxyzvertb,
          * orvalueb = (PyObject *) NULL;
 PyArrayObject * aplane = (PyArrayObject *) NULL,
               * anverts,
               * axyzverts,
               * avalues = (PyArrayObject *) NULL;
 ArrayObject /* * plane = (ArrayObject *) NULL, */
/*              * nverts, */
/*              * xyzverts, */
/*              * values = (ArrayObject *) NULL, */
             * rnverts = (ArrayObject *) NULL,
             * rxyzverts = (ArrayObject *) NULL,
             * rvalues = (ArrayObject *) NULL,
             * rnvertb = (ArrayObject *) NULL,
             * rxyzvertb = (ArrayObject *) NULL,
             * rvalueb = (ArrayObject *) NULL,
/*              * ndxs, */
             * dp = (ArrayObject *) NULL,
             * ndp = (ArrayObject *) NULL,
             * nvertc = (ArrayObject *) NULL,
             * xyzc = (ArrayObject *) NULL,
             * valuec = (ArrayObject *) NULL,
             /* values returned from _slice2_part */
             * rnvertc,
             * rxyzc,
             * rvaluec = (ArrayObject *) NULL,
             * nvertc0 = (ArrayObject *) NULL,
             * xyzc0 = (ArrayObject *) NULL,
             * valuec0 = (ArrayObject *) NULL,
             * keep = (ArrayObject *) NULL,
             * prev = (ArrayObject *) NULL,
             * next = (ArrayObject *) NULL,
             * last = (ArrayObject *) NULL,
             * nkeep,
             * nkeep2,
   /*        * mask0 = (ArrayObject *) NULL, */
             * mask2,
             * list = (ArrayObject *) NULL,
             * listc = (ArrayObject *) NULL;
 double * planed = NULL,
        * xyzvertsd,
        * valuesd = (double *) NULL,
        * rxyzvertsd,
        * rvaluecd = (double *) NULL,
        * rvaluesd = (double *) NULL,
        * rxyzvertbd = NULL,
        * rvaluebd = (double *) NULL,
        * dpd = (double *) NULL,
        * ndpd,
        * xyzcd,
        * xyzc0d,
        * valuecd = (double *) NULL,
        * valuec0d = (double *) NULL;
 int * nvertsd,
     * rnvertsd,
     * rnvertbd = NULL,
   /*     * ndxsd, */
     * nvertcd = NULL,
     * nvertc0d = NULL,
     * nkeepd,
     * nkeep2d,
     * prevd,
     * nextd,
     * lastd,
     * listd = NULL,
     * listcd = NULL;
 Uchar * keepd = NULL,
       * mask2d,
       * valuesc = (Uchar *) NULL,
       * rvaluecc = (Uchar *) NULL,
       * rvaluesc = (Uchar *) NULL,
       * rvaluebc = (Uchar *) NULL,
       * valuecc = (Uchar *) NULL,
       * valuec0c = (Uchar *) NULL;
 char atype = '\0';
 int _slice2x = 0,
     plane_is_scalar,
     isum,
     i, j, k,
     list_length = 0,
     list0_length = 0,
     list1_length = 0,
     list2_length = 0,
     listc_length = 0,
     listc_length1 = 0,
     listc0_length = 0,
     node,
     sumt,
     sumv,
     xdims [2];
 double dplane = 0.0;
 
 if (!PyArg_ParseTuple (args, "OOO|Oi", &oplane, &onverts, &oxyzverts,
    &ovalues, &_slice2x))
    return ERRSS ("slice2: unable to parse arguments.");
 plane_is_scalar = PyFloat_Check (oplane);
 if (plane_is_scalar)
    dplane = PyFloat_AsDouble (oplane);
 else {
    GET_ARR (aplane, oplane, Py_GpReal, 1, PyObject *);
    }
 /* convert arguments to arrays */
 GET_ARR (anverts, onverts, PyArray_INT, 1, PyObject *);
 GET_ARR (axyzverts, oxyzverts, Py_GpReal, 2, PyObject *);
 if (isARRAY (ovalues)) {
    if (A_TYPE (ovalues) == Py_GpReal) {
       GET_ARR (avalues, ovalues, Py_GpReal, 1, PyObject *);
       valuesd = (double *) A_DATA (avalues);
       atype = 'd';
       }
    else if (A_TYPE (ovalues) == Py_GpColor) {
       GET_ARR (avalues, ovalues, Py_GpColor, 1, PyObject *);
       valuesc = (Uchar *) A_DATA (avalues);
       atype = 'b';
       }
    else {
       ERRSS ("Data type for values must be 'b' or 'd'.");
       clearFreeList (0);
       return (PyObject *) NULL;
       }
    }
 /* convert to our type of array for ease of handling */
 if (! plane_is_scalar) {
    planed = (double *) A_DATA (aplane);
    }
 nvertsd = (int *) A_DATA (anverts);
 xyzvertsd = (double *) A_DATA (axyzverts);


 for (isum = 0, i = 0; i < A_SIZE (anverts); i++)
    isum += nvertsd [i];
 if (avalues) {
    if (A_SIZE (avalues) == A_SIZE (anverts))
       node = 0;
    else if (A_SIZE (avalues) == isum)
       node = 1;
    else {
       ERRSS (
       "Number of data values must equal number of nodes or number of cells.");
       clearFreeList (0);
       return (PyObject *) NULL;
       }
    }
 else
    node = 0;
 if (! plane_is_scalar) {
    TRY (dp = allocateArray (A_SIZE (axyzverts) / 3, 'd', 0), 
       PyErr_NoMemory());
    dpd = (double *) (dp->data);
    for (i = 0; i < dp->size; i ++) {
       dpd [i] = xyzvertsd [3 * i] * planed [0] +
          xyzvertsd [3 * i + 1] * planed [1] +
          xyzvertsd [3 * i + 2] * planed [2] - planed [3];
       }
    /* At this point we are done with plane */
    removeFromArrayList ( (PyObject *) aplane);
    }
 else if (plane_is_scalar && avalues && node == 1) {
    TRY (dp = allocateArray (A_SIZE (avalues), 'd', 0), 
       PyErr_NoMemory());
    dpd = (double *) (dp->data);
    if (atype == 'd')
       for (i = 0; i < dp->size; i ++) {
          dpd [i] = valuesd [i] - dplane;
          }
    else
       for (i = 0; i < dp->size; i ++) {
          dpd [i] = (double) valuesc [i] - dplane;
       }
    }
 else {
    ERRSS ("Not sure what kind of slice you're asking for.");
    clearFreeList (0);
    return (PyObject *) NULL;
    }
 /* nkeep is an integer array whose ith entry tells you the number */
 /* of vertices for polygon i which are "above" the isosurface or  */
 /* plane doing the slicing.                                       */
 TRY (nkeep = allocateArray (A_SIZE (anverts), 'i', 0), PyErr_NoMemory());
 nkeepd = (int *) (nkeep->data);
 k = 0;
 for (i = 0; i < nkeep->size; i++) {
    for (j = 0; j < nvertsd [i]; j ++) {
       nkeepd [i] += (int) (dpd [k] >= _slice2_precision) ;
       k ++;
       }
    }
 for (i = 0; i < nkeep->size; i++) {
    /* Compute length of vertex lists */
    /* list1: indices of polygons cut by the plane or isosurface */
    /* list: indices of the vertices of the above polygons       */
    if (nkeepd [i] != 0 && nkeepd [i] != nvertsd [i]) {
       list1_length ++;
       list_length += nvertsd [i];
       }
    }
 if (list1_length != 0) {
   /* nvertc: The number of vertices in each cut polygon.            */
   /* xyzc: The vertex coordinates of each cut polygon.              */
   /* valuec (if given): The values of the function on each vertex   */
   /*    (if node == 1) or on each polygon (if node == 0).           */
   TRY (nvertc = allocateArray (list1_length, 'i', 0), PyErr_NoMemory());
   nvertcd = (int *) (nvertc ->data);
   TRY (list = allocateArray (list_length, 'i', 0), PyErr_NoMemory());
   listd = (int *) (list->data);
   TRY (xyzc = allocateArray (list_length * 3, 'd', 0), PyErr_NoMemory());
   xyzcd = (double *) (xyzc->data);
   if (avalues && node == 0) {
      TRY (valuec = allocateArray (list1_length, atype, 0), 
         PyErr_NoMemory());
      }
   else if (avalues && node == 1) {
      TRY (valuec = allocateArray (list_length, atype, 0), PyErr_NoMemory());
      }
   if (avalues) {
      if (atype == 'd')
         valuecd = (double *) (valuec->data);
      else
         valuecc = (Uchar *) (valuec->data);
   }
   for (i = 0, k = 0, sumv = 0, sumt = 0; i < nkeep->size; i++) {
      if (nkeepd [i] != 0 && nkeepd [i] != nvertsd [i]) {
         nvertcd [k] = nvertsd [i];
         if (avalues && node == 0) {
            if (atype == 'd')
               valuecd [k] = valuesd [i];
            else
               valuecc [k] = valuesc [i];
	 }
         for (j = 0; j < nvertsd [i]; j ++) {
            listd [sumv + j] = sumt + j;
            xyzcd [3 * (sumv + j)] = xyzvertsd [3 * (sumt + j)];
            xyzcd [3 * (sumv + j) + 1] = xyzvertsd [3 * (sumt + j) + 1];
            xyzcd [3 * (sumv + j) + 2] = xyzvertsd [3 * (sumt + j) + 2];
            if (avalues && node == 1) {
               if (atype == 'd')
                  valuecd [sumv + j] = valuesd [sumt + j];
               else
                  valuecc [sumv + j] = valuesc [sumt + j];
	    }
            }
         k ++;
         sumv += nvertsd [i];
         }
      sumt += nvertsd [i];
      }
   }

 if (_slice2x) {
   /* doing a double slice */
   if (! _slice2_precision) 
      if (list1_length == 0) {
         mask2 = (ArrayObject *) NULL;
         mask2d = (void *) NULL;
         nvertc0 = (ArrayObject *) NULL;
         nvertc0d = (void *) NULL;
         valuec0 = (ArrayObject *) NULL;
         valuec0d = (void *) NULL;
         valuec0c = (void *) NULL;
         xyzc0 = (ArrayObject *) NULL;
         xyzc0d = (void *) NULL;
         }
      else {
         TRY (mask2 = logical_not (nkeep, SAVE, 0), PyErr_NoMemory());
         mask2d = (Uchar *) (mask2->data);
         TRY (addToFreeList (nvertc0 = copyArray (nvertc), 0) != -1,
            PyErr_NoMemory ());
         nvertc0d = (int *) (nvertc0->data);
         if (avalues) {
            TRY (addToFreeList (valuec0 = copyArray (valuec), 0) != -1,
               PyErr_NoMemory ());
            if (atype == 'd')
               valuec0d = (double *) (valuec0->data);
            else
               valuec0c = (Uchar *) (valuec0->data);
            }
         TRY (addToFreeList (xyzc0 = copyArray (xyzc), 0) != -1,
            PyErr_NoMemory ());
         xyzc0d = (double *) (xyzc0->data);
         }
   else {
      TRY (nkeep2 = allocateArray (A_SIZE (anverts), 'i', 0), PyErr_NoMemory());
      nkeep2d = (int *) (nkeep->data);
      TRY (mask2 = allocateArray (A_SIZE (anverts), 'b', 0), PyErr_NoMemory());
      mask2d = (Uchar *) (mask2->data);
      k = 0;
      for (i = 0; i < nkeep2->size; i++) {
         for (j = 0; j < nvertsd [i]; j ++) {
            nkeep2d [i] += (int) (dpd [k] >= _slice2_precision) ;
            k ++;
            }
         }
      for (i = 0; i < nkeep2->size; i++) {
         /* Compute length of vertex lists */
         /* list2: indices of polygons cut _slice2_precision below surface */
         /* listc: indices of the vertices of the above polygons           */
         mask2d [i] = ! nkeep2d [i];
         if (nkeep2d [i] != 0 && nkeep2d [i] < nvertsd [i]) {
            list2_length ++;
            listc_length += nvertsd [i];
            }
         }
      if (list2_length != 0) {
         TRY (nvertc0 = allocateArray (list2_length, 'i', 0), 
            PyErr_NoMemory());
         nvertc0d = (int *) (nvertc0->data);
         TRY (listc = allocateArray (listc_length, 'i', 0), 
            PyErr_NoMemory());
         listcd = (int *) (listc->data);
         TRY (xyzc0 = allocateArray (listc_length * 3, 'd', 0), 
            PyErr_NoMemory());
         xyzc0d = (double *) (xyzc0->data);
         if (avalues && node == 0) {
            TRY (valuec0 = allocateArray (list2_length, atype, 0), 
               PyErr_NoMemory());
            if (atype == 'd')
               valuec0d = (double *) (valuec0->data);
            else
               valuec0c = (Uchar *) (valuec0->data);
            }
         else if (avalues && node == 1) {
            TRY (valuec0 = allocateArray (listc_length, atype, 0),
               PyErr_NoMemory());
            if (atype == 'd')
               valuec0d = (double *) (valuec0->data);
            else
               valuec0c = (Uchar *) (valuec0->data);
            }
         for (i = 0, k = 0, sumv = 0; i < nkeep2->size; i++) {
            if (nkeep2d [i] != 0 && nkeep2d [i] != nvertsd [i]) {
               nvertc0d [k] = nvertsd [i];
               if (avalues && node == 0) {
                  if (atype == 'd')
                     valuec0d [k] = valuesd [i];
                  else
                     valuec0c [k] = valuesc [i];
	       }
               for (j = 0; j < nvertsd [i]; j ++) {
                  listcd [sumv + k] = i + j;
                  xyzc0d [3 * (sumv + j)] = xyzvertsd [3 * (i + j)];
                  xyzc0d [3 * (sumv + j) + 1] = xyzvertsd [3 * (i + j) + 1];
                  xyzc0d [3 * (sumv + j) + 2] = xyzvertsd [3 * (i + j) + 2];
                  if (avalues && node == 1) {
                     if (atype == 'd')
                        valuec0d [sumv + j] = valuesd [i + j];
                     else
                        valuec0c [sumv + j] = valuesc [i + j];
		  }
                  }
               k ++;
               sumv += nvertsd [i];
               }
            }
         }
      }
    /* N. B. It's a little confusing, but list2 is being reused here, */
    /* but listc has to stay around for awhile.                       */
    list2_length = 0;
    if (mask2 != (ArrayObject *) NULL)
       for (i = 0, list2_length = 0, listc_length1 = 0; i < mask2->size; i++) {
          list2_length += mask2d [i];
          listc_length1 += mask2d [i] * nvertsd [i];
          }
    if (list2_length != 0) {
      TRY (rnvertb = allocateArray (list2_length, 'i', 0), 
         PyErr_NoMemory());
      rnvertbd = (int *) (rnvertb->data);
      if (avalues && node == 0) {
         TRY (rvalueb = allocateArray (list2_length, atype, 0), 
            PyErr_NoMemory());
         if (atype == 'd')
            rvaluebd = (double *) (rvalueb->data);
         else
            rvaluebc = (Uchar *) (rvalueb->data);
         }
      else if (avalues && node != 0) {
         TRY (rvalueb = allocateArray (listc_length1, atype, 0),
            PyErr_NoMemory());
         if (atype == 'd')
            rvaluebd = (double *) (rvalueb->data);
         else
            rvaluebc = (Uchar *) (rvalueb->data);
         }
      TRY (rxyzvertb = allocateArray (3 * listc_length1, 'd', 0), 
         PyErr_NoMemory());
      rxyzvertbd = (double *) (rxyzvertb->data);
      for (i = 0, k = 0, sumv = 0, sumt = 0; i < A_SIZE (anverts); i ++) {
         if (mask2d [i] != 0) {
            rnvertbd [k] = nvertsd [i];
            if (avalues && node == 0) {
               if (atype == 'd')
                  rvaluebd [k] = valuesd [i];
               else 
                  rvaluebc [k] = valuesc [i];
	    }
            for (j = 0; j < nvertsd [i]; j ++) {
               rxyzvertbd [3 * (sumv + j)] = xyzvertsd [3 * (sumt + j)];
               rxyzvertbd [3 * (sumv + j) + 1] = xyzvertsd [3 * (sumt + j) + 1];
               rxyzvertbd [3 * (sumv + j) + 2] = xyzvertsd [3 * (sumt + j) + 2];
               if (avalues && node != 0) {
                  if (atype == 'd')
                     rvaluebd [(sumv + j)] = valuesd [(sumt + j)];
                  else
                     rvaluebc [(sumv + j)] = valuesc [(sumt + j)];
	       }
               }
            k ++;
            sumv += nvertsd [i];
            }
         sumt += nvertsd [i];
         }
      }
    freeArray (mask2, 0);
   }
/* freeArray (nkeep, 0); Need nkeep and nverts instead of mask0 */

 for (i = 0; i < nkeep->size; i++) {
    /* list0: the uncut polygons.                               */
    /* listc0: the vertices of the uncut polygons.              */
    list0_length += (int) (nkeepd [i] == nvertsd [i]);
    listc0_length += (nkeepd [i] == nvertsd [i]) ? nvertsd [i] : 0;
    }
 if (list0_length < A_SIZE (anverts) ) {
    if (list0_length == 0) {
       rnverts = (ArrayObject *) NULL;
       rxyzverts = (ArrayObject *) NULL;
       rvalues = (ArrayObject *) NULL;
       rnvertsd = (void *) NULL;
       rxyzvertsd = (void *) NULL;
       rvaluesd = (void *) NULL;
       }
    else {
       /* Extract the uncut data. */
       TRY (rnverts = allocateArray (list0_length, 'i', 0), 
          PyErr_NoMemory());
       rnvertsd = (int *) (rnverts->data);
       TRY (rxyzverts = allocateArray (listc0_length * 3, 'd', 0), 
          PyErr_NoMemory());
       rxyzvertsd = (double *) (rxyzverts->data);
       if (avalues && node != 0) {
          TRY (rvalues = allocateArray (listc0_length, atype, 0), 
             PyErr_NoMemory());
          if (atype == 'd')
             rvaluesd = (double *) (rvalues->data);
          else
             rvaluesc = (Uchar *) (rvalues->data);
          }
       else if (avalues) {
          TRY (rvalues = allocateArray (list0_length, atype, 0),
             PyErr_NoMemory());
          if (atype == 'd')
             rvaluesd = (double *) (rvalues->data);
          else
             rvaluesc = (Uchar *) (rvalues->data);
          }
       for (i = 0, k = 0, sumv = 0, sumt = 0; i < nkeep->size; i++) {
          if (nkeepd [i] == nvertsd [i]) {
             rnvertsd [k] = nvertsd [i];
             if (avalues && node == 0) {
                if (atype == 'd')
                   rvaluesd [k] = valuesd [i];
                else
                   rvaluesc [k] = valuesc [i];
	     }
             for (j = 0; j < nvertsd [i]; j++) {
                rxyzvertsd [3 * (sumv + j)] = xyzvertsd [3 * (sumt + j)];
                rxyzvertsd [3 * (sumv + j) + 1] = xyzvertsd [3 * (sumt + j) + 1];
                rxyzvertsd [3 * (sumv + j) + 2] = xyzvertsd [3 * (sumt + j) + 2];
                if (avalues && node != 0) {
                   if (atype == 'd')
                      rvaluesd [(sumv + j)] = valuesd [(sumt + j)];
                   else
                      rvaluesc [(sumv + j)] = valuesc [(sumt + j)];
		}
                }
             k ++;
             sumv += nvertsd [i];
             }
          sumt += nvertsd [i];
          }
       }
    }
 else {
    /* inputs unchanged. But copy them. */
    TRY (rnverts = arrayFromPointer (A_SIZE (anverts), 'i', A_DATA (anverts), 0),
       PyErr_NoMemory());
    rnvertsd = (int *) (rnverts->data);
    TRY (rxyzverts = arrayFromPointer (A_SIZE (axyzverts), 'd', A_DATA (axyzverts), 0),
       PyErr_NoMemory());
    rxyzvertsd = (double *) (rxyzverts->data);
    /* We've given the data pointers of these two or three arrays to others, so they
       no longer own their data. Clear the OWN_DATA flags, so that the
       DECREF applied when removed from the array list will not free their data. */
       UNSET_OWN (anverts);
       UNSET_OWN (axyzverts);
    if (avalues == (PyArrayObject *) NULL)
       rvalues = (ArrayObject *) NULL;
    else {
       TRY (rvalues = arrayFromPointer (A_SIZE (avalues), atype, A_DATA (avalues), 0),
          PyErr_NoMemory());
       if (atype == 'd')
          rvaluesd = (double *) (rvalues->data);
       else
          rvaluesc = (Uchar *) (rvalues->data);
       UNSET_OWN (avalues);
       }
    }
 /* Free these now to avoid leaking both the data memory and the */
 /* memory for each ArrayObject.                                 */
 freeArray (nkeep, 0);
 removeFromArrayList ( (PyObject *) anverts);
 removeFromArrayList ( (PyObject *) axyzverts);
 removeFromArrayList ( (PyObject *) avalues);

 /* done if no partially clipped polys */
 if (list1_length == 0 && listc_length == 0) {
    if (rnverts) {
       RET_ARR (ornverts, 1, & (rnverts->size), PyArray_INT, (char *) rnvertsd,
          PyObject *);
       SET_OWN (ornverts);
       }
    else {
       Py_INCREF (Py_None);
       ornverts = Py_None;
       }
    if (rxyzverts) {
       xdims [0] = rxyzverts->size / 3;
       xdims [1] = 3;
       RET_ARR (orxyzverts, 2, xdims, PyArray_DOUBLE, (char *) rxyzvertsd,
          PyObject *);
       SET_OWN (orxyzverts);
       }
    else {
       Py_INCREF (Py_None);
       orxyzverts = Py_None;
       }
    if (rvalues) {
       if (atype == 'd') {
          RET_ARR (orvalues, 1, & (rvalues->size), atype,
             (char *) rvaluesd, PyObject *);
          }
       else {
          RET_ARR (orvalues, 1, & (rvalues->size), atype,
             (char *) rvaluesc, PyObject *);
          }
       SET_OWN (orvalues);
       }
    else {
       Py_INCREF (Py_None);
       orvalues = Py_None;
       }
    if (rxyzvertb == NULL) {
       /* We remove the following three objects only because addresses */
       /* are in use.                                                  */
       removeArrayOnly (rnverts, 0);
       removeArrayOnly (rxyzverts, 0);
       removeArrayOnly (rvalues, 0);
       clearFreeList (0);
       if (_slice2x) {
          TRY (addToArrayList (oreturn_value = PyList_New (6)),
             PyErr_NoMemory());
          if (PyList_SetItem (oreturn_value, 3, Py_None) < 0 ||
              PyList_SetItem (oreturn_value, 4, Py_None) < 0 ||
              PyList_SetItem (oreturn_value, 5, Py_None) < 0) {
             clearArrayList ();
             return ERRSS ("slice2: unable to assign to return list.");
             }
          else {
             Py_INCREF (Py_None);
             Py_INCREF (Py_None);
             Py_INCREF (Py_None);
             }
          }
       else {
          TRY (addToArrayList (oreturn_value = PyList_New (3)), 
             PyErr_NoMemory());
          }
       if (PyList_SetItem (oreturn_value, 0, ornverts) < 0 ||
           PyList_SetItem (oreturn_value, 1, orxyzverts) < 0 ||
           PyList_SetItem (oreturn_value, 2, orvalues) < 0) {
          clearArrayList ();
          return ERRSS ("slice2: unable to assign to return list.");
          }
       else {
          /* don't want to decref returning values, just clear the list */
          array_list_length = 0;
          return oreturn_value;
          }
       }
    else {
       /* Build the rest */
       RET_ARR (ornvertb, 1, & (rnvertb->size), PyArray_INT, (char *) rnvertbd,
          PyObject *);
       SET_OWN (ornvertb);
       xdims [0] = rxyzvertb->size / 3;
       xdims [1] = 3;
       RET_ARR (orxyzvertb, 2, xdims, PyArray_DOUBLE, (char *) rxyzvertbd,
          PyObject *);
       SET_OWN (orxyzvertb);
       if (avalues) {
          if (atype == 'd') {
             RET_ARR (orvalueb, 1, & (rvalueb->size), atype,
                (char *) rvaluebd, PyObject *);
             }
          else {
             RET_ARR (orvalueb, 1, & (rvalueb->size), atype,
                (char *) rvaluebc, PyObject *);
             }
          SET_OWN (orvalueb);
          }
       else {
          Py_INCREF (Py_None);
          orvalueb = Py_None;
          }
       /* We remove the following six objects only because addresses */
       /* are in use.                                                */
       removeArrayOnly (rnverts, 0);
       removeArrayOnly (rxyzverts, 0);
       removeArrayOnly (rvalues, 0);
       removeArrayOnly (rnvertb, 0);
       removeArrayOnly (rxyzvertb, 0);
       removeArrayOnly (rvalueb, 0);
       clearFreeList (0);
       TRY (addToArrayList (oreturn_value = PyList_New (6)), 
          PyErr_NoMemory());
       if (PyList_SetItem (oreturn_value, 0, ornverts) < 0 ||
           PyList_SetItem (oreturn_value, 1, orxyzverts) < 0 ||
           PyList_SetItem (oreturn_value, 2, orvalues) < 0 ||
           PyList_SetItem (oreturn_value, 3, ornvertb) < 0 ||
           PyList_SetItem (oreturn_value, 4, orxyzvertb) < 0 ||
           PyList_SetItem (oreturn_value, 5, orvalueb) < 0) {
          clearArrayList ();
          return ERRSS ("slice2: unable to assign to return list.");
          }
       else {
          /* don't want to decref returning values, just clear the list */
          array_list_length = 0;
          return oreturn_value;
          }
       }
    }

 if (list1_length != 0) {
    /* get dot products and keep list for the clipped polys */
    ndpd = (double *) malloc (list_length * sizeof (double));
    for (i = 0; i < list_length; i++)
       ndpd [i] = dpd [listd [i]] - _slice2_precision;
    freeArray (dp, 0);
    TRY (dp = arrayFromPointer (list_length, 'd', ndpd, 0),
       PyErr_NoMemory());
    dpd = ndpd;
    freeArray (keep, 0);
    TRY (keep = allocateArray (dp->size, 'b', 0), 
       PyErr_NoMemory());
    keepd = (Uchar *) (keep->data);
    TRY (prev = allocateArray (dp->size, 'i', 0), 
       PyErr_NoMemory());
    prevd = (int *) (prev->data);
    TRY (next = allocateArray (dp->size, 'i', 0), 
       PyErr_NoMemory());
    nextd = (int *) (next->data);
    for (i = 0; i < dp->size; i ++) {
       keepd [i] = (Uchar) (dpd [i] >= 0.0);
       prevd [i] = i - 1;
       nextd [i] = i - 1;
       }
    TRY (last = allocateArray (nvertc->size, 'i', 0), 
       PyErr_NoMemory());
    lastd = (int *) (last->data);
    for (i = 0, isum = 0; i < nvertc->size; i++) {
       isum += nvertcd [i];
       lastd [i] = isum;
       prevd [lastd [i] - nvertcd [i]] = lastd [i] - 1;
       }
    freeArray (nvertc, 0);
    for (i = 0; i < next->size; i++) {
       nextd [prevd [i]] = i;
       }
    if (avalues && node == 1) {
       TRY (_slice2_part (xyzc, keep, next, dp, prev, last, valuec, & rxyzc,
          & rnvertc, & rvaluec, FREE0, FREE0, atype), 
          PyErr_NoMemory());
       }
    else {
       TRY (_slice2_part (xyzc, keep, next, dp, prev, last,
          (ArrayObject *) NULL, & rxyzc, & rnvertc, & rvaluec, FREE0, SAVE, atype), 
          PyErr_NoMemory());
       if (valuec) {
          TRY (rvaluec = allocateArray (valuec->size, atype, 0), PyErr_NoMemory());
          if (atype == 'd') {
             rvaluecd = (double *) rvaluec->data;
             for ( i = 0; i < valuec->size; i ++)
                rvaluecd [i] = valuecd [i];
             }
          else {
             rvaluecc = (Uchar *) rvaluec->data;
             for ( i = 0; i < valuec->size; i ++)
                rvaluecc [i] = valuecc [i];
             }
          }
       }
    TRY (rnverts = concatenate (rnverts, rnvertc, FREE0, FREE0, 0), 
       PyErr_NoMemory());
    rnvertsd = (int *) rnverts->data;
    if (avalues) {
       if (rvaluec == (ArrayObject *) NULL)
          TRY (rvalues = concatenate (rvalues, valuec, FREE0, FREE0, 0),
             PyErr_NoMemory());
       else
          TRY (rvalues = concatenate (rvalues, rvaluec, FREE0, FREE0, 0), 
             PyErr_NoMemory());
       if (atype == 'd')
          rvaluesd = (double *) rvalues->data;
       else
          rvaluesc = (Uchar *) rvalues->data;
       }
    else {
       Py_INCREF (Py_None);
       rvalues = (ArrayObject *) Py_None;
       }
    TRY (rxyzverts = concatenate (rxyzverts, rxyzc, FREE0, FREE0, 0), 
       PyErr_NoMemory());
    rxyzvertsd = (double *) rxyzverts->data;
    freeArray (last, 0);
    freeArray (prev, 0);
    freeArray (next, 0);
    }

 freeArray (list, 0);

 if (_slice2x || list1_length == 0) {
    if (! _slice2_precision)
       for (i = 0; i < keep->size; i++)
          keepd [i] = ! keepd [i];
    else {
       TRY (ndp = allocateArray (listc_length, 'd', 0), 
          PyErr_NoMemory());
       ndpd = (double *) (ndp->data);
       freeArray (keep, 0);
       TRY (keep = allocateArray (listc_length, 'b', 0), 
          PyErr_NoMemory());
       keepd = (Uchar *) (keep->data);
       for (i = 0; i < listc_length; i ++) {
          ndpd [i] = dpd [listcd [i]] + _slice2_precision;
          keepd [i] = ndpd [i] >= 0.0;
          }
       freeArray (dp, 0);
       dp = ndp;
       }
    TRY (prev = allocateArray (keep->size, 'i', 0), 
       PyErr_NoMemory());
    prevd = (int *) (prev->data);
    TRY (next = allocateArray (keep->size, 'i', 0), 
       PyErr_NoMemory());
    nextd = (int *) (next->data);
    for (i = 0; i < keep->size; i ++) {
        prevd [i] = i - 1;
        nextd [i] = i - 1;
       }
    TRY (last = allocateArray (nvertc0->size, 'i', 0), 
       PyErr_NoMemory());
    lastd = (int *) (last->data);
    for (i = 0, isum = 0; i < nvertc0->size; i++) {
       isum += nvertc0d [i];
       lastd [i] = isum;
       prevd [lastd [i] - nvertc0d [i]] = lastd [i] - 1;
       }
    for (i = 0; i < next->size; i ++)
       nextd [prevd [i]] = i;
    if (avalues && node == 1) {
       TRY (_slice2_part (xyzc0, keep, next, dp, prev, last, valuec0,
          & rxyzc, & rnvertc, & rvaluec, FREE0, FREE0, atype), 
          PyErr_NoMemory());
       }
    else {
       TRY (_slice2_part (xyzc0, keep, next, dp, prev, last,
          (ArrayObject *) NULL, & rxyzc, & rnvertc, & rvaluec,
          FREE0, FREE0, atype), 
          PyErr_NoMemory());
       if (valuec0) {
          TRY (rvaluec = allocateArray (valuec0->size, atype, 0), PyErr_NoMemory());
          if (atype == 'd') {
             rvaluecd = (double *) rvaluec->data;
             for ( i = 0; i < valuec0->size; i ++)
                rvaluecd [i] = valuec0d [i];
             }
          else {
             rvaluecc = (Uchar *) rvaluec->data;
             for ( i = 0; i < valuec0->size; i ++)
                rvaluecc [i] = valuec0c [i];
             }
          }
       }
    if (rnvertb == (ArrayObject *) NULL) {
       rnvertb = rnvertc;
       rnvertbd = (int *) rnvertb->data ;
       if (avalues) {
          rvalueb = rvaluec;
          if (atype == 'd')
             rvaluebd = (double *) rvalueb->data;
          else
             rvaluebc = (Uchar *) rvalueb->data;
          }
       else {
          Py_INCREF (Py_None);
          rvalueb = (ArrayObject *) Py_None;
          }
       rxyzvertb = rxyzc;
       rxyzvertbd = (double *) rxyzvertb->data;
       }
    else if (rnvertc != (ArrayObject *) NULL) {
       TRY (rnvertb = concatenate (rnvertb, rnvertc, FREE0, FREE0, 0), 
          PyErr_NoMemory());
       rnvertbd = (int *) rnvertb->data ;
       if (avalues) {
          TRY (rvalueb = concatenate (rvalueb, rvaluec, FREE0, FREE0, 0), 
             PyErr_NoMemory());
          if (atype == 'd')
             rvaluebd = (double *) rvalueb->data;
          else
             rvaluebc = (Uchar *) rvalueb->data;
          }
       else {
          Py_INCREF (Py_None);
          rvalueb = (ArrayObject *) Py_None;
          }
       TRY (rxyzvertb = concatenate (rxyzvertb, rxyzc, FREE0, FREE0, 0), 
          PyErr_NoMemory());
       rxyzvertbd = (double *) rxyzvertb->data;
       }
    freeArray (prev, 0);
    freeArray (next, 0);
    freeArray (last, 0);
    }
  freeArray (listc, 0);
  freeArray (keep, 0);

 /* All done, set up return values. */
 RET_ARR (ornverts, 1, & (rnverts->size), PyArray_INT, (char *) rnvertsd,
    PyObject *);
 SET_OWN (ornverts);
 xdims [0] = rxyzverts->size / 3;
 xdims [1] = 3;
 RET_ARR (orxyzverts, 2, xdims, PyArray_DOUBLE, (char *) rxyzvertsd,
    PyObject *);
 SET_OWN (orxyzverts);
 if (rvalues != (ArrayObject *) Py_None) {
    i = rvalues->size;
    if (atype == 'd') {
       RET_ARR (orvalues, 1, & i, atype, (char *) rvaluesd,
          PyObject *);
       }
    else {
       RET_ARR (orvalues, 1, & i, atype, (char *) rvaluesc,
          PyObject *);
       }
    SET_OWN (orvalues);
    }
 else {
    Py_INCREF (Py_None);
    orvalues = Py_None;
    }
 if (rxyzvertb == NULL) {
    /* We remove the following three objects only because addresses */
    /* are in use.                                                  */
    removeArrayOnly (rnverts, 0);
    removeArrayOnly (rxyzverts, 0);
    removeArrayOnly (rvalues, 0);
    clearFreeList (0);
    TRY (addToArrayList (oreturn_value = PyList_New (3)), 
       PyErr_NoMemory());
    if (PyList_SetItem (oreturn_value, 0, ornverts) < 0 ||
        PyList_SetItem (oreturn_value, 1, orxyzverts) < 0 ||
        PyList_SetItem (oreturn_value, 2, orvalues) < 0) {
       clearArrayList ();
       return ERRSS ("slice2: unable to assign to return list.");
       }
    else {
       /* don't want to decref returning values, just clear the list */
       array_list_length = 0;
       return oreturn_value;
       }
    }
 else {
    /* Build the rest */
    RET_ARR (ornvertb, 1, & (rnvertb->size), PyArray_INT, (char *) rnvertbd,
       PyObject *);
    SET_OWN (ornvertb);
    xdims [0] = rxyzvertb->size / 3;
    xdims [1] = 3;
    RET_ARR (orxyzvertb, 2, xdims, PyArray_DOUBLE, (char *) rxyzvertbd,
       PyObject *);
    SET_OWN (orxyzvertb);
    if (rvalueb != (ArrayObject *) Py_None) {
       if (atype == 'd') {
          RET_ARR (orvalueb, 1, & (rvalueb->size), atype,
             (char *) rvaluebd, PyObject *);
          }
       else {
          RET_ARR (orvalueb, 1, & (rvalueb->size), atype,
             (char *) rvaluebc, PyObject *);
          }
       SET_OWN (orvalueb);
       }
    else {
       Py_INCREF (Py_None);
       orvalueb = Py_None;
       }
    /* We remove the following six objects only because addresses */
    /* are in use.                                                */
    removeArrayOnly (rnverts, 0);
    removeArrayOnly (rxyzverts, 0);
    removeArrayOnly (rvalues, 0);
    removeArrayOnly (rnvertb, 0);
    removeArrayOnly (rxyzvertb, 0);
    removeArrayOnly (rvalueb, 0);
    clearFreeList (0);
    TRY (addToArrayList (oreturn_value = PyList_New (6)), 
       PyErr_NoMemory());
    if (PyList_SetItem (oreturn_value, 0, ornverts) < 0 ||
        PyList_SetItem (oreturn_value, 1, orxyzverts) < 0 ||
        PyList_SetItem (oreturn_value, 2, orvalues) < 0 ||
        PyList_SetItem (oreturn_value, 3, ornvertb) < 0 ||
        PyList_SetItem (oreturn_value, 4, orxyzvertb) < 0 ||
        PyList_SetItem (oreturn_value, 5, orvalueb) < 0) {
       clearArrayList ();
       return ERRSS ("slice2: unable to assign to return list.");
       }
    else {
       /* don't want to decref returning values, just clear the list */
       array_list_length = 0;
       return oreturn_value;
       }
    }
} /* end of slice2 */

static PyObject *unzoom (PyObject * self, PyObject * args)
{
  SETJMP0;
  GdRevertLimits (1);
  Py_INCREF (Py_None);
  return Py_None;
}

/* Linear search for keyword in an array of allowable keywords */
static int verify_kw (char *keyword, char * kwlist[])
{
  int i;

  for (i = 0; kwlist[i]; i++){
    if (0 == strcmp (keyword, kwlist[i])){
      return 1; /* found it */
    }
  }
  return 0; /* didn't find it */
}

static PyObject *viewport (PyObject * self, PyObject * args)
{
  double xmin, xmax, ymin, ymax;
  xmin = gistD.trans.viewport.xmin;
  xmax = gistD.trans.viewport.xmax;
  ymin = gistD.trans.viewport.ymin;
  ymax = gistD.trans.viewport.ymax;
  return Py_BuildValue ("dddd", xmin, xmax, ymin, ymax);
}

static PyObject *window (PyObject * self, PyObject * args, PyObject * kd)
{
  int n, nGiven;
  Drawing *drawing;
  GpColorCell *palette;
  static char *windowKeys[]= {
    "display", "dpi", "private", "hcp", "legends", "dump", "style", "wait",
    "width", "height", 0 };
  PyObject * kwt[NELT(windowKeys) - 1];
  int nColors = 0;

  SETJMP0;			/* See Xerror_longjmp() */

  if (!PyArg_ParseTuple (args, "|i", &n))
    return ERRSS ("window takes zero or one non-keyword integer argument.");

  if(PyTuple_Size(args) == 1) { /* Window() was called with an argument. */
    if (n < 0 || n > 7)
      return ERRSS("graphics windows are numbered from 0 to 7");
    nGiven = (!ghDevices[n].display && !ghDevices[n].hcp);
  } else { /* No argument was given. */
    n = curPlotter;
    nGiven = (n < 0);
    if (nGiven)
      n = 0;
  }

  BUILD_KWT(kd, windowKeys, kwt);

  curElement = -1;

  /* get current palette for this graphics window */
  nColors = GhGetPalette (n, &palette);

    /* check for width and height specs */
#ifndef NO_XLIB
  if (kwt[8]) {
    extern int gx75width, gx100width;
    int width;
    SETKW(kwt[8], width, setkw_integer, windowKeys[8]);
    if (width>30) gx75width= gx100width= width;
    else { gx75width= 450; gx100width= 600; }
  }
  if (kwt[9]) {
    extern int gx75height, gx100height;
    int height;
    SETKW(kwt[9], height, setkw_integer, windowKeys[9]);
    if (height>30) gx75height= gx100height= height;
    else { gx75height= 450; gx100height= 600; }
  }
#endif

  if (nGiven || kwt[0] || kwt[1]) {
    /* display= and/or dpi= keywords */
    char *display = 0;
    int dpi = defaultDPI;
    Engine *engine = ghDevices[n].display;	/* current display engine */

    SETKW(kwt[0], display, setkw_string, windowKeys[0]);
    if (kwt[1]) {
      if (engine)
	return ERRSS ("cannot change dpi of an existing graphics window");
      SETKW(kwt[1], dpi, setkw_integer, windowKeys[1]);
      if (dpi != 100 && dpi != 75)
	return ERRSS ("dpi=100 or dpi=75 are only legal values");
    }
    if (engine) {
      ghDevices[n].display = 0;
      GpKillEngine (engine);
    }
    if (nGiven ? (!display || display[0]) : (display && display[0])) {
#ifndef NO_XLIB
      engine = DISPLAY_ENGINE (windowNames[n], 0, dpi, display);
      if (!engine)
	return ERRSS("failed to open X display or create X window");
      ghDevices[n].display = engine;
      if (palette)
	GhSetPalette (n, palette, nColors);
#else
      return ERRSS ("No interactive graphics in this Pygist -- hcp only");
#endif
    }
  }
  if (kwt[2]) {
    /* private= keyword -- turn on/off private X window colormap */
    int private;
    if (!ghDevices[n].display)
      return ERRSS ("private= keyword not legal without display engine");
    SETKW(kwt[2], private, setkw_boolean, windowKeys[2]);
    GhDumpColors (n, 0, private);
  }
  if (kwt[3]) {
    /* hcp= keyword -- make a new hcp file */
    Engine *engine = ghDevices[n].hcp;
    char *hcp = 0;

    SETKW(kwt[3], hcp, setkw_string, windowKeys[3]);
    if (engine) {
      ghDevices[n].hcp = 0;
      GpKillEngine (engine);
      SetHCPname (n, (char *) 0);
    }
    if (hcp && hcp[0]) {
      long len = strlen (hcp);
      if (len > 3 && strcmp (&hcp[len - 3], ".ps") == 0) {
	engine = GpPSEngine (windowNames[n], 0, hcpDump, SetHCPname (n, hcp));
	if (!engine)
	  return ERRSS ("failed to create PostScript file");
      } else {
	engine = GpCGMEngine (windowNames[n], 0, hcpDump, SetHCPname (n, hcp));
	if (!engine)
	  return ERRSS ("failed to create binary CGM file");
      }
      ghDevices[n].hcp = engine;
      if (palette)
	GhSetPalette (n, palette, nColors);
    }
  }
  if (kwt[4] || kwt[3] || nGiven || kwt[0] || kwt[1]) {
      /* legends= keyword -- turn on/off legend dumping to hcp file */
    int legends;
    SETKW(kwt[4], legends, setkw_boolean, windowKeys[4]);
    if (kwt[4])
      ghDevices[n].doLegends = legends;
    else
      ghDevices[n].doLegends = defaultLegends;
  }
  if (kwt[5]) {
    /* dump= keyword -- turn on/off colormap dumping to hcp file */
    int dump;
    if (!ghDevices[n].hcp)
      return ERRSS (
	"dump= keyword not legal without hcp engine -- use hcp_file");
    SETKW(kwt[5], dump, setkw_boolean, windowKeys[5]);
    GhDumpColors (n, 1, dump);
  }
  if (!ghDevices[n].display && !ghDevices[n].hcp) {
    /* shut down this graphics window completely */
    drawing = ghDevices[n].drawing;
    ghDevices[n].drawing = 0;
    if (drawing)
      GdKillDrawing (drawing);
    GhDeletePalette (n);
    paletteSize = 0;
    if (n == curPlotter) {
      /* highest numbered remaining window becomes current window */
      for (n = 7; n >= 0; n--)
	if (ghDevices[n].display || ghDevices[n].hcp)
	  break;
      curPlotter = n;
      GhSetPlotter (n);
      if (n >= 0) {
	Engine *engine = ghDevices[n].display;
	if (!engine)
	  engine = ghDevices[n].hcp;
	if (engine)
	  paletteSize = GpGetPalette (engine, &palette);
      }
    }
  } else {
    if (kwt[6]) {
      /* style= keyword -- make new drawing */
      char *style = 0;
      SETKW(kwt[6], style, setkw_string, windowKeys[6]);
      drawing = ghDevices[n].drawing;
      if (drawing) {
	ghDevices[n].drawing = 0;
	GdKillDrawing (drawing);
      }
      if (!style || !style[0])
	style = defaultStyle;
      ghDevices[n].drawing = drawing = GdNewDrawing (style ? style : "work.gs");

    } else if (!ghDevices[n].drawing) {
      /* supply default drawing */
      ghDevices[n].drawing = drawing =
	  GdNewDrawing (defaultStyle ? defaultStyle : "work.gs");

    } else {
      drawing = ghDevices[n].drawing;
    }

    if (!drawing) {
      ghDevices[n].drawing = drawing = GdNewDrawing ("work.gs");
      if (drawing)
	return ERRSS ("failed to create drawing -- bad style sheet name?");
      else
	return ERRSS (
	  "failed to create drawing -- Gist work.gs style sheet missing");
    }
    /* make this window current */
    curPlotter = n;
    GhSetPlotter (n);
    paletteSize = nColors;

    /* wait= keyword -- pause until X window is exposed */
    if (kwt[7]) {
      int wait;
      SETKW(kwt[7], wait, setkw_boolean, windowKeys[7]);
      if (1 == wait) GhWaitDisplay ();
    }
  }

  return Py_BuildValue ("i",n);
}

static PyObject *zoom_factor (PyObject * self, PyObject * args)
{
  if (!PyArg_ParseTuple (args, "d", &DISPLAY_ZOOM_FACTOR)) {
    return ERRSS ("Zoomfactor takes one floating point argument.");
  }
  /*
   * avert various disasters -- doesn't address DISPLAY_ZOOM_FACTOR==1.0,
   * which would be frustrating...
   */
  if (DISPLAY_ZOOM_FACTOR < 0.0)
    DISPLAY_ZOOM_FACTOR = -DISPLAY_ZOOM_FACTOR;
  if (DISPLAY_ZOOM_FACTOR < 0.05)
    DISPLAY_ZOOM_FACTOR = 0.05;
  else if (DISPLAY_ZOOM_FACTOR > 20.0)
    DISPLAY_ZOOM_FACTOR = 20.0;
  Py_INCREF (Py_None);
  return Py_None;
}

#ifdef __cplusplus
}
#endif
