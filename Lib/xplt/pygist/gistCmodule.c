/* 
 *  $Id$
 *  --------------------------------------------------------------------
 *  Copyright (c) 1996, 1997, The Regents of the University of California.
 *  All rights reserved.  See Legal.htm for full text and disclaimer. 
 *
 *  NAME:     gistCmodule.c
 *
 *  PURPOSE:  Gistmodule provides glue between Python and the Gist library 
 *            of graphics routines.  Much of the code is this module is 
 *            modeled after similar code in Dave Munro's Yorick interpreter.  
 *
 *  AUTHORS:  Lee Busby (originator)
 *            Zane Motteler (maintainer)
 *            Dave Munro (upgrade changes; event handling interface)
 *            Lila Chase (maintainer)
 *               Lawrence Livermore National Laboratory
 *            Michiel de Hoon (contributor; Windows/Cygwin ports) 
 *               University of Tokyo, Institute of Medical Science 
 *
 *  VERSION:  Original version was based on graph.c in Yorick 1.3.
 *            This version has been mostly upgraded to Yorick 1.5.
 *
 *  CHANGES:
 *  02/25/03 llc Fix mouse command so data is returned.
 *  02/13/03 llc Disable initial pygist printout upon import (except 
 *               when DEBUG flag is on).
 *  02/13/03 dpg Add plremove implementation. 
 *               Correct bug in pli (iMax anad JMax were reversed).
 *  01/13/03 mdh Add back on omitted bug fix for Windows (wait_for_expose).
 *  12/27/02 mdh Windows/Cygwin adjustments (e.g. pyg_on_idle calls for
 *               Windows); fix formats for pointers.
 *  12/04/02 dhm more fixes to get pygist working under windows, _tkinter
 *               Macros for proper printing to stdout and stderr.
 *               Add pyg_unhook, pyg_idler, pyg_pending, and pyg_register.
 *  11/26/02 llc Add documentation on plg markers.
 *  11/25/02 dhm Rework of play resulted in change to PyOS_InputHook and
 *               other improvements (multi-line strings).
 *  11/24/02 dhm fix event/exception handling interface
 *  11/15/02 mhd Add u_enterloop for third approach to interfacing Python
 *               to gist.  Save the last approach as USE_U_WAITER.
 *  11/13/02 llc Add back option to switch between two interfaces to Python.
 *               (using define USE_RL_GETC)
 *  11/11/02 llc Add Numeric to arrayobject.h include, so it is found
 *               automatically when the Python include directory is searched.
 *               (suggestion from MDH)
 *  11/08/02 llc Use the PyFPE_START_PROTECT and PyFPE_END_PROTECT 
 *               macros defined in the standard Python include pyfpe.h.
 *               Add a dummy argument to PyFPE_END_PROTECT.
 *  11/07/02 llc Clean up some compiler warnings:
 *               Get rid of SETERR macro; used just once.
 *               Correct format for writes of longs in debug_array, PrintColor.
 *               Add explict casts to ERRSS where required.
 *               Remove unused variables.
 *               Cast u_waiter properly, and note possible issue.
 *  11/01/02 mdh Replace rl_getc_function with PyOS_InputHook.
 *  09/04/02 llc Deep copy string returned from PyString_AsString
 *               in setkw_string.
 *  12/17/01 llc Add documentation lines to slice2.
 *  12/06/01 llc Remove fix for E_TEXT type entries in legend; these
 *               are legitimate.
 *  12/03/01 llc Discovered that new setting of zc in plf causes 
 *               unaligned accesses on DEC Alpha; go back to earlier
 *               form setting temporary zc1.
 *               With this change, test suite runs to completion.
 *  11/29/01 llc Set rgb for triples.
 *               bytscl: Change len to long.
 *  11/16/01 llc Update to 1.5's legend handling. 
 *  11/12/01 llc Unpack color triples.
 *               Colors changed from int to unsigned long (1.4 to 1.5).
 *               Check for negative colors.
 *  11/08/01 dhm Turn on idle function to do graphics tasks.
 *  11/05/01 llc Move all of the function documentation from gist.help
 *               to doc strings for each PyObject, so it can be accessed 
 *               with pydoc.
 *  11/03/01 llc Started from version with changes from Dave Munro:
 *               turn on rl_getc_function, and add readline.h include. 
 *               Must use with version of python built with readline.
 *  11/19/02 llc Use gist interface provided by Dave Munro for PyOS_InputHook.
 *               Remove old solutions now that play/unix has changed.
 *  01/17/03 llc Remove remaining warnings with gcc 2.96 by:
 *               - Reduce to one dummy for file
 *               - Remove dumpFreeList, Add1 cumsum, equal, greater, 
 *                 greater_equal, less, logical_and, not_equal, 
 *                 WeightedHist if not INCLUDE_EXTRA_GISTCODE.
 *  02/20/03 llc Correct mouse return; see MouseCallBack.
 *
 * Modified in December, 1997 to coexist with the readline library used
 * in Python's interpreter. William Magro, Cornell Theory Center.
 *
 *  --------------------------------------------------------------------
 */

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

#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <setjmp.h>
#include <string.h>

#include "Python.h"
#include "pyfpe.h"

#include "Numeric/arrayobject.h"
#include "hlevel.h"

#include "pstdlib.h"
#include "play.h"
#include "pmin.h"

/* primitive allowance for other non-X windows systems (untested) */
#ifndef DISPLAY_ENGINE
#  define DISPLAY_ENGINE GpFXEngine
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

static int dummy = 0;

/* We add a component to the default Gist search path, for style and
   palette files in our Python distribution.
*/
static char *gistpath = 0, *oldgistpath = 0;
/*  
 *  10/30/01 llc The new gist directory is no longer in graphics.
 *               Original OUR_SPECIAL_DIR "/graphics/gist" 
 */
#define OUR_SPECIAL_DIR "/gist"

/* Mouse() related stuff */
#ifndef NO_MOUSE
static int MouseCallBack (Engine * engine, int system,
			  int release, double x, double y,
			  int butmod, double xn, double yn);
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

static void clearFreeList (int n);
static void clearArrayList (void);
static void clearMemList (void);

/* PyErr_SetString returns a void; make an object and cast it */

/* 
 *  11/06/01 LLC: 
 *  On Teracluster, compiler warns of cast to int in some uses. 
 *
 *  ERRSS must evaluate to a non-null object connected to error
 *  The return value is used in SETKW and TRY.
 *
 *  EXPLANATION of ERRSS:
 *  ---------------------
 *  PyErr_SetString returns a void.  A trick is used here to 
 *  cast it into a NULL (0).  In the following:
 *     void c();
 *     char * y;
 *     y = (char *) (c(),0);
 *  y is 0 or NULL.
 */

#define ERRSS(s) ((PyObject *)(PyErr_SetString(GistError,s),0))
#define ERRMSG(s) (PyErr_SetString(GistError,s))
#define SETJMP0 if(setjmp(jmpbuf)){p_pending_events();return(0);}
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
#define TRYS(e) {char * errstr; \
                 if( (errstr=(e)) != NULL ) { \
                       clearArrayList(); \
                       clearFreeList(0); \
                       clearMemList(); \
                       return ERRSS(errstr); \
                 }} 
#define DECL_ZERO(type, var) type var = 0
#define DECREF_AND_ZERO(p) do{Py_XDECREF(p);p=0;}while(0)

/* (DHM) proper Python printing to stdout and stderr */
#define TO_STDOUT PySys_WriteStdout
#define TO_STDERR PySys_WriteStderr
static void flush_stdout(void);
static void flush_stderr(void);
static void flush_std(const char *s);
static int pyg_puts(const char *);

static void
flush_std(const char *s)
{
  /* code from Python/sysmodule.c:PySys_WriteStdout */
  PyObject *pstdout, *error_type, *error_value, *error_traceback;
  PyErr_Fetch(&error_type, &error_value, &error_traceback);
  pstdout = PySys_GetObject("stdout");
  fflush(pstdout? PyFile_AsFile(pstdout) : stdout);
  PyErr_Restore(error_type, error_value, error_traceback);
}
static void
flush_stdout(void)
{
  flush_std("stdout");
}
static void
flush_stderr(void)
{
  flush_std("stderr");
}
static int
pyg_puts(const char *s)
{
  if (s) {
    long len = strlen(s);
    if (len > 0) {
      TO_STDOUT(s);
      if (s[len-1] == '\n') return 0;
    }
    TO_STDOUT("\n");
  }
  return 0;
}

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
#ifdef INCLUDE_EXTRA_GISTCODE 
static void dumpFreeList (int n);
#endif

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

#ifdef INCLUDE_EXTRA_GISTCODE 
static void dumpFreeList (int n) {
   /* Useful for debugging ??? */
   int i;
   TO_STDOUT ("-----------start-%d-----------\n", n); flush_stdout();
   for (i =0; i < freeListLen [n]; i++) {
      TO_STDOUT ("entry %p points to %c data (%d) at %p.\n",
         freeList [n] [i], freeList [n] [i]->typecode, freeList [n] [i]->size,
         freeList [n] [i]->data); flush_stdout();
      }
   TO_STDOUT ("----------finish-------------\n"); flush_stdout();
   }
#endif

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
static PyObject *pyg_fma (PyObject * self, PyObject * args);
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
static PyObject *plremove (PyObject * self, PyObject * args);
static PyObject *plsys (PyObject * self, PyObject * args);
static PyObject *plt (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *plv (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *redraw (PyObject * self, PyObject * args);
static PyObject *set_slice2_precision (PyObject * self, PyObject * args);
static PyObject *slice2 (PyObject * self, PyObject * args);
static PyObject *unzoom (PyObject * self, PyObject * args);
static PyObject *viewport (PyObject * self, PyObject * args);
static PyObject *window (PyObject * self, PyObject * args, PyObject * kd);
static PyObject *zoom_factor (PyObject * self, PyObject * args);
static PyObject *pyg_unhook (PyObject * self, PyObject * args);
static PyObject *pyg_idler (PyObject * self, PyObject * args);
static PyObject *pyg_register (PyObject * self, PyObject * args);

/* Utility routines */
static GpColor *PushColors(double *z, long len, double zmin,
  double zmax, double scale, double offset);
static char *GetHCPname (int n);
static char *SetHCPname (int n, char *name);
static char *expand_pathname (const char *name);
static double Safe_dbl (double x);
static double *CopyLevels(double *levels, long nLevels);
static char *CheckDefaultWindow(void);
static int CheckPalette (void);
static int GetTypeface (char *s, int *f);
static int GrabByteScale ( PyObject **kwt, char **keywrds, double *scale,
  double *offset, double *zn, double *zx, double *z, int *reg, 
    int region, long iMax, long jMax, int zCompressed);
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
static int setkw_color (PyObject * v, unsigned long *t, char *kw);
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
static int unpack_color_tuple (PyObject * ob, unsigned long color_triple[3]);
static int unpack_limit_tuple (PyObject * ob, double limits[], int *flags);
static int verify_kw (char *keyword, char * kwlist[]);
static long FindMeshZone(double xx, double yy, double *x, double *y, 
  int *reg, long ix, long jx);
static long Safe_strlen(const char *s);
static void AllocTmpLegend(long len);
static void CheckDefaultPalette (void);
static void CleanUpGraphics (void);
static void ForceNewline (void);
static void FreeTmpLegend(void);
static long escape_count(char *arg);
static void escape_cat(char *leg, char *arg);
static void GetPCrange (double *zmn, double *zmx, double *z, 
  int *reg, int region, long iMax, long jMax);
static void GetZCrange(double *zmn, double *zmx, double *z, 
  int *reg, int region, long iMax, long jMax, int zCompressed);
static int LegendAndHide(char *func, char *arg1, char *arg2, char *arg3,
  char *arg4, PyObject *kwt[], char *keys[]);
static void PermitNewline (int nSpaces);
static void PrintColor (char *line, unsigned long color, int suffix);
static void PrintFunc (const char *s);
static void PrintHideLegend (char *line, int type);
static void PrintInit (int (*puts) (const char *));
static void PrintMarks (char *line, int suffix);
static void PrintRegion (char *line, int suffix);
static void PrintSuffix (int suffix);
static void PrintTypeWidth (char *line, int suffix);
static void clear_pyMsh(void);
static void get_mesh(GaQuadMesh *m);

static int pyg_on_idle(void);
static void pyg_on_connect(int dis, int fd);
static void pyg_on_keyline(char *msg);
static void pyg_got_expose(void);
static void pyg_got_alarm(void *context);

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

static int already_initialized = 0;
static int paletteSize = 0;
static int maxColors = 200; /* maximum number of colors for GpReadPalette */
static int hcpDump= 1;      /* whiners can't figure out how to dump colors */
static int hcpPSdefault = 1;  /* default now .ps instead of .cgm */
static int hcpOnFMA = 0;
static int defaultDPI = 75;
static int defaultLegends= 1;
static char *defaultStyle = 0;
static char *defaultPalette = 0;

static int curPlotter = -1;
static int curElement = -1;

static char *hcpNames[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static char *tmpLegend = 0;
static char *windowNames[8] = {
  "Pygist 0", "Pygist 1", "Pygist 2", "Pygist 3",
  "Pygist 4", "Pygist 5", "Pygist 6", "Pygist 7"
};

/* Next few variables are used by plq(), which does some fancy printing. */
/* static long prop3sizes[10] = {0, 8, 2, 5, 5, 3, 3, 7, 0, 0}; */
/* static long prop4sizes[10] = {0, 8, 1, 3, 1, 1, 3, 4, 4, 0}; */
/* static long prop5sizes[10] = {0, 3, 5, 2, 5, 6, 7, 9, 3, 5}; */
static int curIX = -1, curIXc = -1;
static char specialMarkers[5] = ".+*ox";
static int (*RawPrinter) (const char *s);
static int printLength = 79;	/* maximum number of characters on a line */
/* static int lenPrintBuf = 79; */
static char printBuf[80];
static long maxPrintLines = 5000;
static int printNow, permitNow;
static long printLines;
static double _slice2_precision = 0.0;

/* make p_abort work with SETJMP0
 * here's how it works:
 * p_abort calls pyg_abort_hook, which longjmps to SETJMP0 point
 * there, p_pending_events is called to process the
 * on_exception event that goes along with the p_abort
 */
static jmp_buf jmpbuf;
static jmp_buf pyg_jmpbuf;
static void pyg_abort_hook(void);
static void pyg_abort_hook(void)
{
  longjmp(pyg_jmpbuf, 1);
}
static int pyg_wait_flag = 0;
static void pyg_on_exception(int signal, char *errmsg);
static void pyg_on_exception(int signal, char *errmsg)
{
  pyg_wait_flag = 0;
  if (!errmsg) {
    PyErr_SetString (GistError, "unknown gist error");
  } else {
    PyErr_SetString (GistError, errmsg);
  }
}

static char gist_module_documentation[] =
"Gist Graphics Package, version1.5"
;

#define PYCF (PyCFunction)
/*#define PYCFWK (PyCFunctionWithKeywords)*/
#define PYCFWK (PyCFunction) /* Make the compiler shut up. */
#define KWFLG (METH_VARARGS | METH_KEYWORDS)


/*#######################################################################*/
/*           Auxiliary routines used by the slice2 suite                 */
/*#######################################################################*/

#ifdef INCLUDE_EXTRA_GISTCODE 
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
#endif
 
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

#ifdef INCLUDE_EXTRA_GISTCODE

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
   /* We may be using res uninitialized here -- MDH */
   for (; tar < (Uchar *) (res->data) + res->size; src1 ++, src2 ++, tar ++)
      * tar = * src1 && * src2;
   if (freea >= 0)
      freeArray (a, freea);
   if (freeb >= 0)
      freeArray (b, freea);
   return (res);
   }
#endif

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

#ifdef INCLUDE_EXTRA_GISTCODE
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
#endif

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

#ifdef INCLUDE_EXTRA_GISTCODE
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
#endif

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
   int * nextd=0,
       * nvertcd=0,
       * prevd=0,
       * list0d=0,
       * list1d=0,
       * listd=0,
       * xoldd=0;
   /*#######*/int * ndxsd;
   double * xyz0d=0,
          * valcd = (double *) NULL,
          * valc0d = (double *) NULL,
          * valc1d = (double *) NULL,
          * dpd=0,
          * xyzcd=0,
          * xyz1d=0,
          * xyzc_newd=0,
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
      if (valc != (ArrayObject *)NULL)  {
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
      if (valc != (ArrayObject *)NULL)  {
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
      if (valc != (ArrayObject *)NULL)  {
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

static char *CheckDefaultWindow(void)
{
  int i;
  for (i=0 ; i<8 ; i++) if (ghDevices[i].drawing) {
    if (!ghDevices[i].display && !ghDevices[i].hcp) {
      Drauing *drawing= ghDevices[i].drawing;
      ghDevices[i].drawing= 0;
      GdKillDrawing(drawing);
      curElement= -1;
    }
  }
  if (curPlotter<0) {
    for (i=0 ; i<8 ; i++) if (ghDevices[i].drawing)
      return ("graphics window killed -- use window command to re-select");
    ghDevices[0].drawing=
      GdNewDrawing(defaultStyle? defaultStyle : "work.gs");
    curElement= -1;
    if (!ghDevices[0].drawing)
      return ("failed to create drawing -- Gist work.gs style sheet missing");
    ghDevices[0].doLegends= defaultLegends;

#ifndef NO_XLIB
    gist_private_map = gist_rgb_hint = 0;
    ghDevices[0].display=
      DISPLAY_ENGINE(windowNames[0], 0, defaultDPI, (char *)0);
    if (!ghDevices[0].display)
      return ("failed to open X display or create X window");
#else
    ghDevices[0].display= 0;
    ghDevices[0].hcp= hcpDefault;
    hcpDefault= 0;
#endif

    curPlotter= 0;
    GhSetPlotter(0);
  }
  return NULL;
}

static void CheckDefaultPalette (void)
{
  GpColorCell *palette;
  GhGetPalette (curPlotter, &palette);
  if (!palette)
    paletteSize = GhReadPalette (curPlotter,
       defaultPalette ? defaultPalette : "earth.gp", &palette, maxColors);
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
  long i, j, k;
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

  if (!have_min_max)  {
     ERRMSG ( "Unable to find maximum and minimum of data??");
  }
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
    else  {
      return (int) ERRSS ( "illegal font keyword suffix -- B is bold, I is italic");
    }
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

/*  -------------------------------------------------------------------- */
/*  contour */

static char contour__doc__[] =
"[nc, yc, xc] = contour (level, z, y, x [, ireg] [, triangle = <vals>]\n"
"   [, region = num])\n"
"     returns the points on the contour curve that would have been\n"
"     plotted by plc.  Z, Y, X, and IREG are as for plc, and the\n"
"     triangle= and region= keywords are accepted and have the same\n"
"     meaning as for plc.  Unlike plc, the triangle array is an output\n"
"     as well as an input to contour; if supplied it may be modified\n"
"     to reflect any triangulations which were performed by contour.\n"
"\n"
"     either:\n"
"     LEVEL is a scalar z value to return the points at that contour\n"
"     level.  All such points lie on edges of the mesh.  If a contour\n"
"     curve closes, the final point is the same as the initial point\n"
"     (i.e.- that point is included twice in the returned list).\n"
"\n"
"     or:\n"
"     LEVEL is a pair of z values [z0,z1] to return the points of\n"
"     a set of polygons which outline the regions between the two\n"
"     contour levels.  These will include points on the mesh boundary\n"
"     which lie between the levels, in addition to the edge points\n"
"     for both levels.  The polygons are closed, simply connected,\n"
"     and will not contain more than about 4000 points (larger polygons\n"
"     are split into pieces with a few points repeated where the pieces\n"
"     join).\n"
"\n"
"     YC and XC are the output points on the curve(s), or None if there\n"
"     are no points. The return value NC is a list of the lengths of\n"
"     the polygons/polylines returned in (XC,YC), or None if there are\n"
"     none.  len(XC) == len(YC) == sum(NC).  For the level pair\n"
"     case, YC, XC, and NC are ready to be used as inputs to plfp.\n"
"\n"
"   KEYWORDS: triangle, region\n"
"\n"
"   SEE ALSO: plc, plfp\n";

#undef N_KEYWORDS
#define N_KEYWORDS 2
static char *cntrKeys[N_KEYWORDS+1]= { "triangle", "region", 0 };

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
  PyObject * olevels,
           * zop,
           * kwt [NELT (cntrKeys) - 1],
           * retval;
  PyArrayObject * zap,
                * alevels;
  PyObject      * anp,
                * axcp,
                * aycp;
  double levels [2],
         * lev,
         * xcp,
         * ycp;
  double * z;
  char *errstr =
    "contour requires 2D arguments (levels, z [region = num, triangle = <vals>] )";

  if (!pyMsh.y)  {
    return ERRSS ("contour: no current mesh - use plmesh(y, x) to initialize");
  }

  n = PyTuple_Size (args);
  /* contour (levels, z [, region = num]) */
  if (n != 2)  {
     return ERRSS ("contour requires 2 positional parameters (levels and z).");
  }
  BUILD_KWT (kd, cntrKeys, kwt);
  TRY ( PyArg_ParseTuple (args, "OO", &olevels, &zop),
     ERRSS ("contour: unable to parse arguments."));

  GET_ARR (zap, zop, PyArray_DOUBLE, 2, PyObject *);
  dims [0] = A_DIM (zap, 0);
  dims [1] = A_DIM (zap, 1);
  if (dims [0] != A_DIM (pyMsh.y, 0) || dims [1] != A_DIM (pyMsh.y, 1)) {
     return ERRSS ("z array must have same dimensions as mesh in contour.");
  }
  /* Provide a triangle if none supplied */
  if ( !pyMsh.triangle )
     TRY (pyMsh.triangle = (PyArrayObject *) PyArray_FromDims (2, dims, PyArray_SHORT),
        ERRSS ("contour: unable to create triangle array."));

  /* LLC:  
   *  1.5 has new first keyword "triangle" before "region", so
   *  change kwt and cntrKeys indices from 0 to 1.
   *  Add setz_mesh call. 
   */

  /* kwt[0] ("triangle=") is handled by setz_mesh. */
  /* Skip levels arg */

  {  PyObject * newargs;
     n = PyTuple_Size(args);
     TRY (newargs = PyTuple_GetSlice (args, 1, n), 0);
     TRY (setz_mesh (newargs, &zop, errstr, kwt[0]), (PyObject *) NULL);
  }
  if (!pyMsh.y)  {
    return ERRSS ("No current mesh - set (y, x) first");
  }

  gistD.region = 0;
  SETKW (kwt [1], gistD.region, setkw_integer, cntrKeys[1]);

  get_mesh (&mesh);

  /* Figure out the contour levels */
  if (isARRAY (olevels)) {
     GET_ARR (alevels, olevels, PyArray_DOUBLE, 1, PyObject *);
     lev = (double *) A_DATA (alevels);
     nlevels = A_SIZE (alevels);
     if (nlevels > 2) {
        clearArrayList ();
        return ERRSS ("contour: only 1 or 2 levels allowed."); 
     }
     for (i = 0; i < nlevels; i++)
        levels [i] = lev [i];
     removeFromArrayList ( (PyObject *) alevels);
     }
  /* levels argument can be scalar--allow Int or Float */
  else if (PyFloat_Check (olevels) || PyInt_Check (olevels)) {
     nlevels = 1;
     if (PyFloat_Check (olevels))
        levels [0] = (double) PyFloat_AsDouble (olevels);
     else
        levels [0] = (double) PyInt_AsLong (olevels);
     }
  else {
     clearArrayList ();
     return ERRSS ("contour: levels argument is wrong type."); 
  }

  z = (double *) A_DATA (zap);
  ntotal = (nlevels == 2) ?
     GcInit2 (&mesh, gistD.region, z, levels, 30L, &nparts):
     GcInit1 (&mesh, region, z, levels [0], &nparts);

  /* The following is necessary because I must own references to objects */
  /* that go in the list to be returned. Py_None will be INCREF'ed 3     */
  /* times when it is put on this list.                                  */
  if ( !(retval = Py_BuildValue ("[O,O,O]", Py_None, Py_None, Py_None))) {
     clearArrayList ();
     return ERRSS ("contour: unable to create return value list."); 
  }
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
  NEW_MEM (xcp, ntotal, double, PyObject *);
  RET_ARR ( axcp, 1, &ntotal, PyArray_DOUBLE, (char *) xcp, PyObject *);
  SET_OWN (axcp);
  NEW_MEM (ycp, ntotal, double, PyObject *);
  RET_ARR ( aycp, 1, &ntotal, PyArray_DOUBLE, (char *) ycp, PyObject *);
  SET_OWN (aycp);

  i = GcTrace (np, xcp, ycp);
  if ( i != ntotal) {
     clearArrayList ();
     clearMemList ();
     return ERRSS ("contour: GcTrace has failed.");
  }
  /* For some reason, if PyList_SetItem fails, it returns -1. */
  if (own_triangle) {
     Py_DECREF (kwt [0]);
  }
  if (PyList_SetItem (retval, 0, anp) < 0 ||
      PyList_SetItem (retval, 1, aycp) < 0 ||
      PyList_SetItem (retval, 2, axcp) < 0) {
     clearArrayList ();
     clearMemList ();
     return ERRSS ("contour was unable to build return list.");
  }
  removeFromArrayList ( (PyObject *) zap);
  mem_list_length = 0;
  array_list_length = 0;
  return retval;
}

#ifndef NO_MOUSE
static int MouseCallBack (Engine * engine, int system,
			  int release, double x, double y,
			  int butmod, double xn, double yn)
{
  int n = curPlotter;
  if (n < 0 || ghDevices[n].display != engine) {

    pyg_wait_flag = 0;
    /*  2/24/03 LLC:  
     *  Remove setting mouseError (Yorick does not use mouseError):
     *     mouseError = 1;
     *  Setting mouseError prevents data return for mouse.
     *  In Yorick, button press results in one call to MouseCallBack.
     *  On button release, one MouseCallBack call is made to capture
     *  end mouse coordinates, and another MouseCallBack is made with 
     *  engine = 0, at which time the results are pushed on the stack.
     */
    return 1;
  } else if (mouseError || release==-1) {
    pyg_wait_flag = 0;
    mouseError = 1;  /* Leave this one for now. LLC */
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
    pyg_wait_flag = 0;  /* LLC: Should this be here? */
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
/* 11/12/01 llc Change color from int to unsigned long */ 

static void PrintColor (char *line, unsigned long color, int suffix)
{
  if (color >= 0) {
    sprintf (line, "color= %ld,", color);
    PrintFunc (line);
  } else if (color == P_FG)
    PrintFunc ("color= \"fg\"");
  else if (color == P_BG)
    PrintFunc ("color= \"bg\"");
  else if (color == P_RED)
    PrintFunc ("color= \"red\"");
  else if (color == P_GREEN)
    PrintFunc ("color= \"green\"");
  else if (color == P_BLUE)
    PrintFunc ("color= \"blue\"");
  else if (color == P_CYAN)
    PrintFunc ("color= \"cyan\"");
  else if (color == P_MAGENTA)
    PrintFunc ("color= \"magenta\"");
  else if (color == P_YELLOW)
    PrintFunc ("color= \"yellow\"");
  else if (color == P_GREEN)
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
  sprintf (line, "marks= %d,  mcolor= 0x%02lx,  ",
           gistA.dl.marks, gistA.m.color);
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
  return (int) ERRSS ( "you appear to have Aa00 through Zz00 hcp files -- clean up");

got1:
  if (!hcpPSdefault)
    hcpDefault = GpCGMEngine ("Pygist default", 0, hcpDump,
			    SetHCPname (-1, hcpName));
  else
    hcpDefault = GpPSEngine ("Pygist default", 0, hcpDump,
			    SetHCPname (-1, hcpName));

  if (!hcpDefault)  {
    return (int) ERRSS ("failed to create default hcp file");
  }

  return 1;
}

static char *SetHCPname (int n, char *name)
{
  char *now;
  if (n < 0 || n > 7)
    n = 8;
  now = hcpNames[n];
#ifdef WINDOWS
#ifndef CYGWIN
  hcpNames[n] = name;
#else
  hcpNames[n] = expand_pathname (name);
#endif
#else
  hcpNames[n] = expand_pathname (name);
#endif
  if (now)
    free (now);
  return hcpNames[n];
}

/* Used only by mouse() */
static int YPrompt(const char *s)
{
  TO_STDOUT("%s", s);
  flush_stdout();
  return 0;
}

/*  -------------------------------------------------------------------- */

static char animate__doc__[] =
"animate()\n"
"or animate( 0/1 )\n"
"     Without any arguments, toggle animation mode; with argument 0,\n"
"     turn off animation mode; with argument 1 turn on animation mode.\n"
"     In animation mode, the X window associated with a graphics window\n"
"     is actually an offscreen pixmap which is bit-blitted onscreen\n"
"     when an fma() command is issued.  This is confusing unless you are\n"
"     actually trying to make a movie, but results in smoother animation\n"
"     if you are.  Generally, you should turn animation on, run your movie,\n"
"     then turn it off.\n"
"\n"
"   SEE ALSO: window, fma, plg\n";

static PyObject *animate (PyObject * self, PyObject * args)
{
  int i = 3;			/* default is to toggle */

  SETJMP0;
  if (!PyArg_ParseTuple (args, "|i", &i))  {
    return ERRSS ("Animate takes zero or one argument.");
  }
  PyFPE_START_PROTECT("animate", return 0)
  TRYS(CheckDefaultWindow())
  GhFMAMode (2, i);
  PyFPE_END_PROTECT(dummy)
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

  for (i = 0; (kw = kwlist[i]); i++) kwt[i] = 0;
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
      (int) ERRSS (errstr);
      return -1;
    }
  }
  Py_DECREF(keylist);

  /* Ok, all keywords were legal.  Now store pointers to their value.
   * Note that PyDict_GetItemString() puts 0 in kwt[i] if
     that key isn't found. */
  for (i = 0; (kw = kwlist[i]); i++)
    if((kwt[i] = PyDict_GetItemString (kd, kw)) != 0)
      ++nkw_set;
    /* I tried PyMapping_GetItemString() above, but kept getting
     * "KeyError: wait" messages back from Python.
     */

  return nkw_set;
}

/*  -------------------------------------------------------------------- */
/*  bytscl */

static char bytscl__doc__[] =
"bytscl(z)\n"
"or bytscl(z, top=max_byte, cmin=lower_cutoff, cmax=upper_cutoff)\n"
"     Returns a char array of the same shape as Z, with values linearly\n"
"     scaled to the range 0 to one less than the current palette size.\n"
"     If MAX_BYTE is specified, the scaled values will run from 0 to\n"
"     MAX_BYTE instead.\n"
"     If LOWER_CUTOFF and/or UPPER_CUTOFF are specified, Z values outside\n"
"     this range are mapped to the cutoff value; otherwise the linear\n"
"     scaling maps the extreme values of Z to 0 and MAX_BYTE.\n"
"\n"
"   SEE ALSO: plf, pli, histeq_scale\n";

#undef N_KEYWORDS
#define N_KEYWORDS 3
static char *bsKeys[N_KEYWORDS+1]= { "top", "cmin", "cmax", 0 };

static PyObject *bytscl (PyObject * self, PyObject * args, PyObject * kd)
{
  PyObject *zop, *kwt[NELT (bsKeys) - 1];
  PyArrayObject *zap, *zcap;
  double *z, zmin, zmax, scale, offset;
  GpColor *zc, *zc1;
  int i;
  long len;

  if (!PyArg_ParseTuple (args, "O", &zop))  {
    return ERRSS ("bytscl requires exactly one non-keyword argument");
  }

  TRY (addToArrayList((PyObject *)(zap = (PyArrayObject *)
      PyArray_ContiguousFromObject (zop, PyArray_DOUBLE, 1, 0))),
      (PyObject *)PyErr_NoMemory ());
  z = (double *) A_DATA (zap);
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

/*  -------------------------------------------------------------------- */
/*  current_window */

static char current_window__doc__[] =
"n = current_window()\n"
"     Return the number of the current graphics window, or -1 if none.\n";

static PyObject *current_window (PyObject * self, PyObject * args)
{
  return PyInt_FromLong (curPlotter);
}

/*  -------------------------------------------------------------------- */

/* The following routine has been added to check the integrity of what */
/* we think might be a NumPy array, including looking at addresses.    */

static char debug_array__doc__[] =
"None.";

static PyObject *debug_array (PyObject * self, PyObject * args)
{
 PyObject *oarray;
 PyArrayObject * aarray;
 int i;
 int max;
 long mmax;
 TRY (PyArg_ParseTuple (args, "O", &oarray),
    ERRSS ("debug_array: argument should be one PyObject*."));
 TO_STDOUT("Value of input pointer is %p.", oarray); flush_stdout();
 TO_STDOUT(" Reference count %d, size %d.\n", oarray->ob_refcnt,
           oarray->ob_type->ob_size);
 flush_stdout();
 if (! isARRAY (oarray)) {
    return ERRSS ("debug_array: argument should be a NumPy array.");
 }
 aarray = (PyArrayObject *) oarray;
 TO_STDOUT("Data pointer: %p; nd %d; dim1 %d; type %c.\n", aarray->data,
   aarray->nd, aarray->dimensions [0], aarray->descr->type); flush_stdout();
 if (aarray->descr->type == 'i') {
    TO_STDOUT ("%d ", ( (int *)(aarray->data)) [0]); flush_stdout();
    for (i = 1, max = ( (int *)(aarray->data)) [0]; i < aarray->dimensions [0]; i ++){
       if ( ( (int *)(aarray->data)) [i] > max) max = ( (int *)(aarray->data)) [i];
       TO_STDOUT ("%d ", ( (int *)(aarray->data)) [i]);
       if (i % 10 == 0) TO_STDOUT ("\n");
       flush_stdout();
       }
    TO_STDOUT ("maximum value is %d.\n", max); flush_stdout();
    }
 else if (aarray->descr->type == 'l') {
    TO_STDOUT ("%ld ", ( (long *)(aarray->data)) [0]); flush_stdout();
    for (i = 1, mmax = ( (long *)(aarray->data)) [0]; i < aarray->dimensions [0]; i ++){
       if ( ( (long *)(aarray->data)) [i] > mmax) mmax = ( (long *)(aarray->data)) [i];
       TO_STDOUT ("%ld ", ( (long *)(aarray->data)) [i]);
       if (i % 10 == 0) TO_STDOUT ("\n");
       flush_stdout();
       }
    TO_STDOUT ("maximum value is %ld.\n", mmax); flush_stdout();
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
  if(!PyErr_Occurred()) ERRSS (errstr ? errstr : "error in expand_path") ; 
  DECREF_AND_ZERO(p1);
  DECREF_AND_ZERO(p2);
  DECREF_AND_ZERO(p3);
  DECREF_AND_ZERO(p4);
  return 0;
}

/*  -------------------------------------------------------------------- */

static char fma__doc__[] =
"fma()\n"
"     Frame advance the current graphics window.  The current picture\n"
"     remains displayed in the associated X window until the next element\n"
"     is actually plotted.\n"
"\n"
"   SEE ALSO: window, hcp, animate, plg\n";

static PyObject *pyg_fma (PyObject * self, PyObject * args)
{
  SETJMP0;

  TRYS(CheckDefaultWindow())

  if (hcpOnFMA) {
    if (!CheckPalette ())
      return NULL;
  }
  curElement = -1;
  GhFMA ();
  Py_INCREF (Py_None);
  return Py_None;
}

/*  -------------------------------------------------------------------- */

/* Set pointers in the GaQuadMesh struct from values in the current
 * pyMsh struct.  Naturally, pyMsh must be fully initialized before
 * this is called.
 */
static void get_mesh(GaQuadMesh *m)
{
  m->iMax = A_DIM (pyMsh.y, 1);
  m->jMax = A_DIM (pyMsh.y, 0);
  m->y = (double *) A_DATA (pyMsh.y);
  m->x = (double *) A_DATA (pyMsh.x);
  m->reg = (int *) A_DATA (pyMsh.reg);
  if (isARRAY (pyMsh.triangle))
    m->triangle = (short *) A_DATA (pyMsh.triangle);
  else
    m->triangle = 0; /* Gist will provide a default in this case. */
}

/*  -------------------------------------------------------------------- */

static char get_slice2_precision__doc__[] =
"None.";

static PyObject* get_slice2_precision (PyObject * self, PyObject * args)
{
 if (PyTuple_Size (args) > 0)  {
    return ERRSS ("get_slice2_precision takes no arguments.") ;
 }
 return Py_BuildValue ( "d", _slice2_precision) ;
}

/*  -------------------------------------------------------------------- */
/*  gridxy */

static char gridxy__doc__[] =
"gridxy( flag )\n"
"or gridxy( xflag, yflag )\n"
"     Turns on or off grid lines according to FLAG.  In the first form, both\n"
"     the x and y axes are affected.  In the second form, XFLAG and YFLAG\n"
"     may differ to have different grid options for the two axes.  In either\n"
"     case, a FLAG value of 0 means no grid lines (the default), a value of\n"
"     1 means grid lines at all major ticks (the level of ticks which get\n"
"     grid lines can be set in the style sheet), and a FLAG value of 2 means\n"
"     that the coordinate origin only will get a grid line.  In styles with\n"
"     multiple coordinate systems, only the current coordinate system is\n"
"     affected.  The keywords can be used to affect the style of the grid\n"
"     lines.\n"
"\n"
"     You can also turn the ticks off entirely.  (You might want to do this\n"
"     to plot your own custom set of tick marks when the automatic tick\n"
"     generating machinery will never give the ticks you want.  For example\n"
"     a latitude axis in degrees might reasonably be labeled `0, 30, 60,\n"
"     90', but the automatic machinery considers 3 an `ugly' number - only\n"
"     1, 2, and 5 are `pretty' - and cannot make the required scale.  In\n"
"     this case, you can turn off the automatic ticks and labels, and use\n"
"     plsys, pldj, and plt to generate your own.)\n"
"     To fiddle with the tick flags in this general manner, set the\n"
"     0x200 bit of FLAG (or XFLAG or YFLAG), and `or-in' the 0x1ff bits\n"
"     however you wish.  The meaning of the various flags is described\n"
"     in the `work.gs' Gist style sheet.  Additionally, you can use the\n"
"     0x400 bit to turn on or off the frame drawn around the viewport.\n"
"     Here are some examples:\n"
"        gridxy(0x233)          work.gs default setting\n"
"        gridxy(0, 0x200)       like work.gs, but no y-axis ticks or labels\n"
"        gridxy(0, 0x231)       like work.gs, but no y-axis ticks on right\n"
"        gridxy(0x62b)          boxed.gs default setting\n"
"\n"
"   KEYWORDS: color, type, width\n"
"\n"
"   SEE ALSO: window, plsys, limits, ylimits, logxy\n";

#undef N_KEYWORDS
#define N_KEYWORDS 6
static char *gridKeys[N_KEYWORDS+1]= {
  "color", "type", "width", "base60", "degrees", "hhmm", 0 };

static PyObject *gridxy (PyObject * self, PyObject * args, PyObject * kd)
{
  int xgrid = 0, ygrid = 0, narg;
  PyObject * kwt[NELT(gridKeys) - 1];

  SETJMP0;

  if (!PyArg_ParseTuple (args, "|ii", &xgrid, &ygrid)) {
    return ERRSS ("gridxy takes zero, one or two non-keyword arguments.");
  }
  /* If one argument is given, use it for both x and y. */
  if((narg = PyTuple_Size(args)) == 1)
    ygrid = xgrid;

  TRYS(CheckDefaultWindow())

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
#ifdef WINDOWS
  pyg_on_idle();
#endif
  Py_INCREF (Py_None);
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  hcp */

static char hcp__doc__[] =
"hcp(), or hcpon(), or hcpoff()\n"
"     The hcp command sends the picture displayed in the current graphics\n"
"     window to the hardcopy file.  (The name of the default hardcopy file\n"
"     can be specified using hcp_file; each individual graphics window may\n"
"     have its own hardcopy file as specified by the window command.)\n"
"     The hcpon command causes every fma (frame advance) command to do\n"
"     and implicit hcp, so that every frame is sent to the hardcopy file.\n"
"     The hcpoff command reverts to the default \"demand only\" mode.\n"
"\n"
"   SEE ALSO: window, fma, plg\n";

static PyObject *hcp (PyObject * self, PyObject * args)
{
  SETJMP0;
  PyFPE_START_PROTECT("hcp", return 0)
  TRYS(CheckDefaultWindow())
  CheckPalette ();
  GhHCP ();
  PyFPE_END_PROTECT(dummy)
  Py_INCREF (Py_None);
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  hcp_file */

static char hcp_file__doc__[] =
"filename = hcp_finish( [n] )\n"
"     Close the current hardcopy file and return the filename.\n"
"     If N is specified, close the hcp file associated with window N\n"
"     and return its name; use hcp_finish(-1) to close the default\n"
"     hardcopy file.\n"
"\n"
"   SEE ALSO: window, fma, hcp, hcp_out, plg\n";

#undef N_KEYWORDS
#define N_KEYWORDS 2
static char *hcpKeys[N_KEYWORDS+1]= { "dump", "ps", 0 };
 
static PyObject *hcp_file (PyObject * self, PyObject * args, PyObject *kd)
{
  Engine *engine = hcpDefault;
  char *hcp = 0;
  int gotDump = 0;
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

/*  -------------------------------------------------------------------- */
/*  hcp_finish */

static char hcp_finish__doc__[] =
"filename = hcp_finish( [n] )\n"
"     Close the current hardcopy file and return the filename.\n"
"     If N is specified, close the hcp file associated with window N\n"
"     and return its name; use hcp_finish(-1) to close the default\n"
"     hardcopy file.\n"
"\n"
"   SEE ALSO: window, fma, hcp, hcp_out, plg\n";

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

/*  -------------------------------------------------------------------- */
/*  hcpoff */

static char hcpoff__doc__[] =
"hcp(), or hcpon(), or hcpoff()\n"
"     The hcp command sends the picture displayed in the current graphics\n"
"     window to the hardcopy file.  (The name of the default hardcopy file\n"
"     can be specified using hcp_file; each individual graphics window may\n"
"     have its own hardcopy file as specified by the window command.)\n"
"     The hcpon command causes every fma (frame advance) command to do\n"
"     and implicit hcp, so that every frame is sent to the hardcopy file.\n"
"     The hcpoff command reverts to the default `demand only' mode.\n"
"\n"
"   SEE ALSO: window, fma, plg\n";

static PyObject *hcpoff (PyObject * self, PyObject * args)
{
  SETJMP0;
  PyFPE_START_PROTECT("hcpoff", return 0)
  TRYS(CheckDefaultWindow())
  hcpOnFMA = 0;
  GhFMAMode (0, 2);
  PyFPE_END_PROTECT(dummy)
  Py_INCREF (Py_None);
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  hcpon */

static char hcpon__doc__[] =
"hcp(), or hcpon(), or hcpoff()\n"
"     The hcp command sends the picture displayed in the current graphics\n"
"     window to the hardcopy file.  (The name of the default hardcopy file\n"
"     can be specified using hcp_file; each individual graphics window may\n"
"     have its own hardcopy file as specified by the window command.)\n"
"     The hcpon command causes every fma (frame advance) command to do\n"
"     and implicit hcp, so that every frame is sent to the hardcopy file.\n"
"     The hcpoff command reverts to the default `demand only' mode.\n"
"\n"
"   SEE ALSO: window, fma, plg\n";

static PyObject *hcpon (PyObject * self, PyObject * args)
{
  SETJMP0;
  PyFPE_START_PROTECT("hcpon", return 0)
  TRYS(CheckDefaultWindow())
  hcpOnFMA = 1;
  GhFMAMode (1, 2);
  PyFPE_END_PROTECT(dummy)
  Py_INCREF (Py_None);
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  limits */

static char limits__doc__[] =
"old_limits = limits()\n"
"or old_limits = limits( xmin [, xmax, ymin, ymax,]\n"
"     [ square=0/1, nice=0/1, restrict=0/1 ] )\n"
"or limits( old_limits )\n"
"\n"
"     In the first form, restore all four plot limits to extreme values,\n"
"     and save the previous limits in the tuple old_limits.\n"
"\n"
"     In the second form, set the plot limits in the current coordinate\n"
"     system to XMIN, XMAX, YMIN, YMAX, which may each be a number to fix\n"
"     the corresponding limit to a specified value, or the string `e'\n"
"     to make the corresponding limit take on the extreme value of the\n"
"     currently displayed data. Arguments may be omitted from the right\n"
"     end only. (But see ``ylimits'' to set limits on the y-axis.)\n"
"\n"
"     If present, the square keyword determines whether limits marked as\n"
"     extreme values will be adjusted to force the x and y scales to be\n"
"     equal (square=1) or not (square=0, the default). If present, the\n"
"     nice keyword determines whether limits will be adjusted to nice\n"
"     values (nice=1) or not (nice=0, the default). There is a subtlety\n"
"     in the meaning of `extreme value' when one or both of the limits\n"
"     on the OPPOSITE axis have fixed values -- does the `extreme value'\n"
"     of the data include points which will not be plotted because their\n"
"     other coordinate lies outside the fixed limit on the opposite axis\n"
"     (restrict=0, the default), or not (restrict=1)?\n"
"\n"
"     Limits() always returns a tuple of 4 doubles and an integer;\n"
"     OLD_LIMITS[0:3] are the previous xmin, xmax, ymin, and ymax, and\n"
"     OLD_LIMITS[4] is a set of flags indicating extreme values and the\n"
"     square, nice, restrict, and log flags. This tuple can be saved and\n"
"     passed back to limits() in a future call to restore the limits to a\n"
"     previous state.\n"
"\n"
"     In an X window, the limits may also be adjusted interactively with\n"
"     the mouse. Drag left to zoom in and pan (click left to zoom in on a\n"
"     point without moving it), drag middle to pan, and click (and drag)\n"
"     right to zoom out (and pan). If you click just above or below the\n"
"     plot, these operations will be restricted to the x-axis; if you\n"
"     click just to the left or right, the operations are restricted to\n"
"     the y-axis. A shift-left click, drag, and release will expand the\n"
"     box you dragged over to fill the plot (other popular software zooms\n"
"     with this paradigm). If the rubber band box is not visible with\n"
"     shift-left zooming, try shift-middle or shift-right for alternate\n"
"     XOR masks. Such mouse-set limits are equivalent to a limits command\n"
"     specifying all four limits EXCEPT that the unzoom command can\n"
"     revert to the limits before a series of mouse zooms and pans.\n"
"\n"
"     The limits you set using the limits or ylimits functions carry over\n"
"     to the next plot -- that is, an fmaoperation does NOT reset the\n"
"     limits to extreme values.\n"
"\n"
"   SEE ALSO: plsys, ylimits, logxy, zoom_factor, unzoom, plg\n";

#undef N_KEYWORDS
#define N_KEYWORDS 3
static char *limKeys[N_KEYWORDS+1]= {
  "square", "nice", "restrict", 0 };

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
      return ( PyObject * ) NULL; 
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
    if (kwt[1])  {
      if(nice) gistD.flags |= D_NICE;
      else gistD.flags &= ~D_NICE;
    }
    if (kwt[2])  {
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
    if(0 == j)  { /* Error */
      return ERRSS ("bad xmin argument: Use float or 'e'");
    }
    else if(1 == j) /* Xmin changed or set to extreme. */
      ++changed;
  }
  if (xmax_ob) {
    j = set_limit (xmax_ob, &gistD.limits.xmax, &gistD.flags, D_XMAX);
    if(0 == j)  { /* Error */
      return ERRSS ("bad xmax argument: Use float or 'e'");
    }
    else if(1 == j) /* Xmax changed or set to extreme. */
      ++changed;
  }
  if (ymin_ob) {
    j = set_limit (ymin_ob, &gistD.limits.ymin, &gistD.flags, D_YMIN);
    if(0 == j) { /* Error */
      return ERRSS ("bad ymin argument: Use float or 'e'");
    }
    else if(1 == j) /* Ymin changed or set to extreme. */
      ++changed;
  }
  if (ymax_ob) {
    j = set_limit (ymax_ob, &gistD.limits.ymax, &gistD.flags, D_YMAX);
    if(0 == j)  { /* Error */
      return ERRSS ("bad ymax argument: Use float or 'e'");
    }
    else if(1 == j) /* Ymax changed or set to extreme. */
      ++changed;
  }

  if (changed) GdSetLimits ();

#ifdef WINDOWS
  pyg_on_idle();
#endif

  return Py_BuildValue ("ddddi",
     old_limits[0], old_limits[1], old_limits[2], old_limits[3], old_flags);
}

/*  -------------------------------------------------------------------- */
/*  logxy */

static char logxy__doc__[] =
"logxy( xflag, yflag )\n"
"     Sets the linear/log axis scaling flags for the current coordinate\n"
"     system. XFLAG and YFLAG may be 0 to select linear scaling, or 1 to\n"
"     select log scaling. YFLAG may be omitted (but not XFLAG).\n"
"\n"
"   SEE ALSO: plsys, limits, ylimits, plg, gridxy\n";

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
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  mesh_loc */

static char mesh_loc__doc__[] =
"mesh_loc(y0, x0)\n"
"or mesh_loc(y0, x0, y, x)\n"
"or mesh_loc(y0, x0, y, x, ireg)\n"
"     Returns the zone index (=i+imax*(j-1)) of the zone of the mesh\n"
"     (X,Y) (with optional region number array IREG) containing the\n"
"     point (X0,Y0).  If (X0,Y0) lies outside the mesh, returns 0.\n"
"     Thus, eg- ireg(mesh_loc(x0, y0, y, x, ireg)) is the region number of\n"
"     the region containing (x0,y0).  If no mesh specified, uses default.\n"
"     X0 and Y0 may be arrays as long as they are conformable.\n"
"\n"
"   SEE ALSO: plmesh, moush, mouse\n";

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

  if (PyTuple_Size (args) < 2)  {
    return ERRSS ("mesh_loc requires at least two arguments");
  }
  TRY (setvu_mesh (args, &y0op, &x0op, errstr), (PyObject *) NULL);
  if (!pyMsh.y)  {
    return ERRSS ("No current mesh - set (y, x) first");
  }
  get_mesh (&mesh);

  if (isARRAY (y0op)) {
    TRY (addToArrayList((PyObject *)(y0ap = (PyArrayObject *)
        PyArray_ContiguousFromObject (y0op, PyArray_DOUBLE, 1, 0))),
        (PyObject *)PyErr_NoMemory ());
    n = A_SIZE ( y0ap );
    TRY (addToArrayList((PyObject *)(x0ap = (PyArrayObject *)
        PyArray_ContiguousFromObject (x0op, PyArray_DOUBLE, 1, 0))),
        (PyObject *)PyErr_NoMemory ());
    if (n != A_SIZE ( x0ap )) {
      clearArrayList();
      return ERRSS ("(y0, x0) must be same size");
    }
    y0 = (double *) A_DATA (y0ap);
    x0 = (double *) A_DATA (x0ap);
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

/*  -------------------------------------------------------------------- */
/*  mfit */

static char mfit__doc__[] =
"Computes multiquadric fit to data; used for contour plotting\n"
"of random data. Calling sequence from Python:\n"
"   zcplot = mfit (alpha, x, xcplot, y, ycplot, rsqmqd)\n"
"where alpha are the interpolation coefficients, x and y\n"
"are the original randomly distributed coordinates\n"
"(alpha, x, and y are all the same length, say nzcplot).\n"
"xcplot (nxcplot long) and ycplot (nycplot long) specify\n"
"an overlying rectangular mesh. rsqmod is a scalar peculiar\n"
"to the problem.\n";

static PyObject *mfit (PyObject * self, PyObject * args)
{
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
          *oycplot;
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

/*  -------------------------------------------------------------------- */
/*  mouse */

static char mouse__doc__[] =
"result = mouse(system, style, prompt)\n"
"     Displays a PROMPT, then waits for a mouse button to be pressed,\n"
"     then released.  Returns tuple of length eleven:\n"
"       result= [x_pressed, y_pressed, x_released, y_released,\n"
"                xndc_pressed, yndc_pressed, xndc_released, yndc_released,\n"
"                system, button, modifiers]\n"
"\n"
"     If SYSTEM>=0, the first four coordinate values will be relative to\n"
"     that coordinate system.\n"
"     For SYSTEM<0, the first four coordinate values will be relative to\n"
"     the coordinate system under the mouse when the button was pressed.\n"
"     The second four coordinates are always normalized device coordinates,\n"
"     which start at (0,0) in the lower left corner of the 8.5x11 sheet of\n"
"     paper the picture will be printed on, with 0.0013 NDC unit being\n"
"     1/72.27 inch (1.0 point).  Look in the style sheet for the location\n"
"     of the viewport in NDC coordinates (see the style keyword).\n"
"\n"
"     If STYLE is 0, there will be no visual cues that the mouse\n"
"     command has been called; this is intended for a simple click.\n"
"     If STYLE is 1, a rubber band box will be drawn; if STYLE is 2,\n"
"     a rubber band line will be drawn.  These disappear when the\n"
"     button is released.\n"
"\n"
"     Clicking a second button before releasing the first cancels the\n"
"     mouse function, which will then return nil.\n"
"     Ordinary text input also cancels the mouse function, which again\n"
"     returns nil.\n"
"\n"
"     The left button reverses forground for background (by XOR) in\n"
"     order to draw the rubber band (if any).  The middle and right\n"
"     buttons use other masks, in case the rubber band is not visible\n"
"     with the left button.\n"
"\n"
"     result[8] is the coordinate system in which the first four\n"
"     coordinates are to be interpreted.\n"
"     result[9] is the button which was pressed, 1 for left, 2\n"
"     for middle, and 3 for right (4 and 5 are also possible).\n"
"     result[10] is a mask representing the modifier keys which\n"
"     were pressed during the operation: 1 for shift, 2 for shift lock,\n"
"     4 for control, 8 for mod1 (alt or meta), 16 for mod2, 32 for mod3,\n"
"     64 for mod4, and 128 for mod5.\n"
"\n"
"   SEE ALSO: moush\n";

static PyObject *mouse (PyObject * self, PyObject * args)
{
#ifdef DISPLAY_MOUSE
  char *prompt = 0;
  int system = -1, style = 0;
  int n = curPlotter;

  SETJMP0;
  if (n < 0 || !ghDevices[n].display)  {
    return ERRSS ("no current graphics window for mouse function");
  }

  if (!PyArg_ParseTuple (args, "|iis", &system, &style, &prompt))  {
    return ERRSS ("call with (system, style, prompt)");
  }

  /* GhWaitDisplay (); */   /* otherwise can lock up */
  GhBeforeWait ();          /* be sure display is current */
  if (!prompt)
    YPrompt (defaultPrompts[style != 0]);
  else if (prompt[0])
    YPrompt (prompt);
  mouseError = 0;
  mouseError |= DISPLAY_MOUSE (ghDevices[n].display, style, system,
			       &MouseCallBack);
  if (!prompt || prompt[0])
    YPrompt ("\n");

  if (!mouseError) {
    pyg_wait_flag = 1;
    p_wait_while(&pyg_wait_flag);
  }

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

/*  -------------------------------------------------------------------- */
/*  palette */

static char palette__doc__[] =
"palette( filename )\n"
"or palette( source_window_number )\n"
"or palette( red, green, blue, ntsc=1/0 )\n"
"or palette( red, green, blue, gray )\n"
"or palette( red, green, blue, query=1 )\n"
"or palette( red, green, blue, gray, query=1 )\n"
"     Set (or retrieve with query=1) the palette for the current\n"
"     graphics window.  The FILENAME is the name of a Gist palette file;\n"
"     the standard palettes are `earth.gp', `stern.gp', `rainbow.gp',\n"
"     `heat.gp', `gray.gp', and `yarg.gp'.  Use the maxcolors keyword\n"
"     in the pldefault command to put an upper limit on the number of\n"
"     colors which will be read from the palette in FILENAME.\n"
"\n"
"     In the second form, the palette for the current window is copied\n"
"     from the SOURCE_WINDOW_NUMBER.  If the X colormap for the window is\n"
"     private, there will still be two separate X colormaps for the two\n"
"     windows, but they will have the same color values.\n"
"\n"
"     In the third form, RED, GREEN, and BLUE are 1-D arrays of the same\n"
"     length specifying the palette you wish to install; the values\n"
"     should vary between 0 and 255, and your palette should have no\n"
"     more than 240 colors.  If ntsc=0, monochrome devices (such as most\n"
"     laser printers) will use the average brightness to translate your\n"
"     colors into gray; otherwise, the NTSC (television) averaging will\n"
"     be used (.30*RED+.59*GREEN+.11*BLUE).  Alternatively, you can specify\n"
"     GRAY explicitly.\n"
"\n"
"     Ordinarily, the palette is not dumped to a hardcopy file\n"
"     (color hardcopy is still rare and expensive), but you can\n"
"     force the palette to dump using the window() or hcp_file() commands.\n"
"\n"
"     See the dump= keyword for the hcp_file() and window() commands if you\n"
"     are having trouble getting color in your hardcopy files.\n"
"\n"
"   SEE ALSO: window, fma, hcp, pldefault, plg\n";

#undef N_KEYWORDS
#define N_KEYWORDS 2
static char *paletteKeys[N_KEYWORDS+1]= { "ntsc", "query", 0 };

static PyObject *palette (PyObject * self, PyObject * args, PyObject * kd)
{
  GpColorCell *palette = 0;
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
    if (query) {
       return ERRSS ("query requires (r,g,b) arrays as arguments");
    }
  
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

  TRYS(CheckDefaultWindow())
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
	red[i]   = P_R(palette[i]);
	green[i] = P_G(palette[i]);
	blue[i]  = P_B(palette[i]);
      }
      if (ngray)
	for (i = 0 ; i < nColors ; i++)
	  gray[i]  = (P_R(palette[i])+P_G(palette[i])+P_B(palette[i]))/3;

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
        palette[i] = P_RGB(red[i], green[i], blue[i]);
	/* if (gray) palette[i].gray = gray[i]; */
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

static void pyg_got_alarm(void *context)
{
  pyg_wait_flag = 0;
}

/*  -------------------------------------------------------------------- */

static char pause__doc__[] =
"pause( milliseconds )\n"
"     Pause for the specified number of milliseconds of wall clock\n"
"     time, or until input arrives from the keyboard.\n"
"     This is intended for use in creating animated sequences.\n";

static PyObject *pyg_pause (PyObject * self, PyObject * args)
{
  long timeout;

  if (!PyArg_ParseTuple (args, "i", &timeout)) {
    return ERRSS ("Pause requires one integer argument.");
  }
  if (timeout < 0)
    timeout = 0;

  p_set_alarm(0.001*timeout, pyg_got_alarm, 0);
  pyg_wait_flag = 1;
  p_wait_while(&pyg_wait_flag);

  Py_INCREF (Py_None);
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  plc */

static double *tmpLevels = 0; /* UPDATE */

static char plc__doc__[] =
"plc( z, y, x, levs=z_values )\n"
"or plc( z, y, x, ireg, levs=z_values )\n"
"or plc( z, levs=z_values )\n"
"     Plot contours of Z on the mesh Y versus X.  Y, X, and IREG are\n"
"     as for plm.  The Z array must have the same shape as Y and X.\n"
"     The function being contoured takes the value Z at each point\n"
"     (X,Y) -- that is, the Z array is presumed to be point-centered.\n"
"     The Y, X, and IREG arguments may all be omitted to default to the\n"
"     mesh set by the most recent plmesh call.\n"
"     The LEVS keyword is a list of the values of Z at which you want\n"
"     contour curves.  The default is eight contours spanning the\n"
"     range of Z.\n"
"     The following keywords are legal (each has a separate help entry):\n"
"\n"
"   KEYWORDS: legend, hide\n"
"             type, width, color, smooth\n"
"             marks, marker, mspace, mphase\n"
"             smooth, triangle, region\n"
"\n"
"   SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp, plmesh\n"
"             limits, logxy, ylimits, fma, hcp\n";

#undef N_KEYWORDS
#define N_KEYWORDS 15
static char *plcKeys[N_KEYWORDS+1]= {
  "legend", "hide", "region", "color", "type", "width",
  "marks", "mcolor", "marker", "msize", "mspace", "mphase",
  "smooth", "triangle", "levs", 0 };

static PyObject *plc (PyObject * self, PyObject * args, PyObject * kd)
{
  GaQuadMesh mesh;
  PyArrayObject *zap;
  PyObject *zop;
  int i;  
  char *z_name= 0, *y_name= 0, *x_name= 0, *r_name= 0;
  long iMax = 0, jMax = 0, nLevels = 0;
  double *z = 0, *levels = 0; /* UPDATE */
  PyObject * kwt[NELT(plcKeys) - 1];
  char *errstr =
    "plc requires 2D arguments (z [ , y, x, ireg, levs = levels ] )";

  SETJMP0;			/* See Xerror_longjmp() */

  if (PyTuple_Size (args) == 0)  {
    return ERRSS ("plc requires at least one argument");
  }
  BUILD_KWT(kd, plcKeys, kwt);
  TRY (setz_mesh (args, &zop, errstr, kwt[13]), (PyObject *) NULL);
  if (!pyMsh.y)  {
    return ERRSS ("No current mesh - set (y, x) first");
  }
  GET_ARR (zap, zop, PyArray_DOUBLE, 2, PyObject *);
  jMax = A_DIM(zap, 0);
  iMax = A_DIM(zap, 1);
  if (A_DIM (pyMsh.y, 0) != jMax || A_DIM (pyMsh.y, 1) != iMax) {
    clearArrayList ();
    return ERRSS ("Z array must match (y, x) mesh arrays in shape");
  }
  z = (double *) A_DATA (zap);
  get_mesh (&mesh);
  if (mesh.iMax!=iMax || mesh.jMax!=jMax)  {
     return ERRSS ("z array must have same dimensions as mesh in plc");
  }

  /* set legend and hide in gistD */
  TRYS(CheckDefaultWindow())
  if ( !LegendAndHide("\001: plc, ", z_name, y_name, x_name, r_name, kwt, plcKeys) )
     return ERRSS ( "Error in plc: LegendAndHide" );

  /* set properties, starting from defaults for decorated polylines */
  GhGetLines();
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

  /* set contour levels */
  if(kwt[14]) { /* levs= keyword */
    PyArrayObject *lap;
    double *lev;

    GET_ARR (lap, kwt[14], PyArray_DOUBLE, 1, PyObject *);
    lev = (double *) A_DATA (lap);
    nLevels = A_SIZE (lap);
    levels = p_malloc (sizeof(double) * nLevels);
    for(i = 0; i < nLevels; i++)
      levels[i] = lev[i];
    if (levels)  {
      levels= CopyLevels(levels, nLevels);
    }
    removeFromArrayList ( (PyObject *) lap);
  } 

  if (!levels) {
    /* create a default set of contour levels now */
    int i;
    double zmin, zmax, step;

    nLevels= 8;
    levels= CopyLevels((double *)0, nLevels);
    GetPCrange(&zmin, &zmax, z, mesh.reg, gistD.region, iMax, jMax);

    step= (zmax-zmin)/8.0;
    levels[0]= zmin+0.5*step;
    for (i=1 ; i<8 ; i++) levels[i]= levels[i-1]+step;
  }

  curElement = -1;
  PyFPE_START_PROTECT("plc", return 0)
  curElement =
    GdContours (NOCOPY_MESH, &mesh, gistD.region, z, levels, (int)nLevels);
  PyFPE_END_PROTECT(dummy)
  Py_DECREF (zap);
  SAFE_FREE (levels);
  array_list_length = 0;
  mem_list_length = 0;
  if (curElement < 0)  {
    return ERRSS ("Gist GdContour plotter failed");
  }
  tmpLevels = 0; /* Gist now owns this pointer */

  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

static double *CopyLevels(double *levels, long nLevels)
{
  long i;
  double *tmp= tmpLevels;
  tmpLevels= 0;
  if (tmp) p_free(tmp);
  tmpLevels= p_malloc(sizeof(double)*nLevels);
  if(!tmpLevels) return tmpLevels; 
  for (i=0 ; i<nLevels ; i++) tmpLevels[i]= levels? levels[i] : 0.0;
  return tmpLevels;
}

/*  -------------------------------------------------------------------- */
/*  pldefault */

static char pldefault__doc__[] =
"pldefault( key1=value1, key2=value2, ... )\n"
"     Set default values for the various properties of graphical elements.\n"
"\n"
"     The keywords can be most of the keywords that can be passed to the\n"
"     plotting commands:\n"
"       plg:  color, type, width,\n"
"             marks, mcolor, msize, mspace, mphase,\n"
"             rays, rspace, rphase, arrowl, arroww\n"
"       pldj: color, type, width\n"
"       plt:  color, font, height, path, justify, opaque\n"
"       plm:  color, type, width\n"
"       plv:  color, hollow, width, aspect\n"
"       plc:  color, type, width,\n"
"             marks, mcolor, marker, msize, mspace, mphase\n"
"       plf:  edges, ecolor, ewidth\n"
"\n"
"     The initial default values are:\n"
"       color=`fg', type=`solid', width=1.0 (1/2 point),\n"
"       marks=1, mcolor=`fg', msize=1.0 (10 points),\n"
"          mspace=0.16, mphase=0.14,\n"
"       rays=0, arrowl=1.0 (10 points), arroww=1.0 (4 points),\n"
"          rspace=0.13, rphase=0.11375,\n"
"       font=`helvetica', height=12.0, path=0, justify=`NN', opaque=0,\n"
"       hollow= 0, aspect=0.125,\n"
"       edges=0, ecolor=`fg', ewidth=1.0 (1/2 point)\n"
"\n"
"     Additional default keywords are:\n"
"       dpi, style, legends  (see window command)\n"
"       palette              (to set default filename as in palette command)\n"
"       maxcolors            (default 200)\n"
"\n"
"   SEE ALSO: window, plsys, plq, pledit, plg\n";

#undef N_KEYWORDS
#define N_KEYWORDS 29
static char *dfltKeys[N_KEYWORDS+1]= {
  "color", "type", "width",
  "marks", "mcolor", "marker", "msize", "mspace", "mphase",
  "rays", "arrowl", "arroww", "rspace", "rphase",
  "font", "height", "orient", "justify", "opaque",
  "hollow", "aspect", "dpi", "style", "legends", "palette", "maxcolors",
  "edges", "ecolor", "ewidth", 0 };

static PyObject *pldefault (PyObject * self, PyObject * args, PyObject * kd)
{
  PyObject * kwt[NELT(dfltKeys) - 1];
  char *errstr = "pldefault takes no non-keyword arguments";
  int dpi, type;

  if(PyTuple_Size(args) > 0) {
    return ERRSS (errstr);
  }
  
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
      return ERRSS ("orient= keyword must be 0, 1, 2, or 3");
    }
  }

  SETKW(kwt[17], dummy,             setkw_justify,  dfltKeys[17]);
  SETKW(kwt[18], gistA.t.opaque,    setkw_boolean,  dfltKeys[18]);
  SETKW(kwt[19], gistA.vect.hollow, setkw_boolean,  dfltKeys[19]);
  SETKW(kwt[20], gistA.vect.aspect, setkw_double,   dfltKeys[20]);

  if(kwt[21]) {
    SETKW(kwt[21], dpi,             setkw_integer,  dfltKeys[21]);
    if (dpi<25) dpi = 25;
    else if (dpi>300) dpi = 300;
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

/*  -------------------------------------------------------------------- */

#undef N_KEYWORDS
#define N_KEYWORDS 5
static char *pldjKeys[N_KEYWORDS+1]= {
  "legend", "hide", "color", "type", "width", 0 };

static char pldj__doc__[] =
"pldj( x0, y0, x1, y1 )\n"
"     Plot disjoint lines from (X0,Y0) to (X1,Y1).  X0, Y0, X1, and Y1\n"
"     may have any dimensionality, but all must have the same number of\n"
"     elements.\n"
"     The following keywords are legal (each has a separate help entry):\n"
"\n"
"   KEYWORDS: legend, hide\n"
"             type, width, color\n"
"\n"
"   SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp\n"
"             limits, logxy, ylimits, fma, hcp\n";

static PyObject *pldj (PyObject * self, PyObject * args, PyObject * kd)
{
  PyObject *op[4];
  PyArrayObject *ap[4];
  double *d[4];
  int i;
  char *x0_name= 0, *y0_name= 0, *x1_name= 0, *y1_name= 0;
  long n= 0;

  PyObject * kwt[NELT(pldjKeys) - 1];
  char *errstr = "pldj requires exactly four non-keyword arguments";

  SETJMP0;

  if (!PyArg_ParseTuple (args, "OOOO", &op[0], &op[1], &op[2], &op[3]))  {
    return ERRSS (errstr);
  }
  
  for (i=0; i<4; i++)
    TRY (addToArrayList ((PyObject *)(ap[i] = (PyArrayObject *)
        PyArray_ContiguousFromObject (op[i], PyArray_DOUBLE, 1, 0))),
        (PyObject *)PyErr_NoMemory ());

  n = A_SIZE ( ap[0] );
  for (i=1; i<4; i++)
    if ( A_SIZE (ap[i]) != n) {
      clearArrayList ();
      return ERRSS ("pldj arguments must all be the same size");
    }

  /* set legend and hide in gistD */
  TRYS(CheckDefaultWindow())

  /* set properties, starting from defaults for simple polylines */
  GhGetMesh();

  BUILD_KWT(kd, pldjKeys, kwt);
  if ( !LegendAndHide("pldj, ", x0_name, y0_name, x1_name, y1_name, kwt, pldjKeys) )
     return ERRSS ( "Error in pldj: LegendAndHide" );
  SETKW(kwt[0],  gistD.legend,    setkw_string,   pldjKeys[0]);
  SETKW(kwt[1],  gistD.hidden,    setkw_boolean,  pldjKeys[1]);
  SETKW(kwt[2],  gistA.l.color,   setkw_color,    pldjKeys[2]);
  SETKW(kwt[3],  gistA.l.type,    setkw_linetype, pldjKeys[3]);
  SETKW(kwt[4],  gistA.l.width,   setkw_double,   pldjKeys[4]);

  for (i=0; i<4; i++)
    d[i] = (double *) A_DATA (ap[i]);

  curElement = -1;
  PyFPE_START_PROTECT("pldj", return 0)
  curElement = GdDisjoint (n, d[0], d[1], d[2], d[3]);
  PyFPE_END_PROTECT(dummy)
  clearArrayList ();
  if (curElement < 0)  {
    return ERRSS ("Gist GdDisjoint plotter failed");
  }

  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

/*  -------------------------------------------------------------------- */

static char pledit__doc__[] =
"pledit( key1=value1, key2=value2, ... )\n"
"or pledit( n_element, key1=value1, key2=value2, ... )\n"
"or pledit( n_element, n_contour, key1=value1, key2=value2, ... )\n"
"     Changes some property of element number N_ELEMENT (and contour\n"
"     number N_CONTOUR of that element).  If N_ELEMENT and N_CONTOUR are\n"
"     omitted, the default is the most recently added element, or the\n"
"     element specified in the most recent plq query command.\n"
"\n"
"     The keywords can be any of the keywords that apply to the current\n"
"     element.  These are:\n"
"       plg:  color, type, width,\n"
"             marks, mcolor, marker, msize, mspace, mphase,\n"
"             rays, rspace, rphase, arrowl, arroww,\n"
"             closed, smooth\n"
"       pldj: color, type, width\n"
"       plt:  color, font, height, path, justify, opaque\n"
"       plm:  region, boundary, inhibit, color, type, width\n"
"       plf:  region\n"
"       plv:  region, color, hollow, width, aspect, scale\n"
"       plc:  region, color, type, width,\n"
"             marks, mcolor, marker, msize, mspace, mphase\n"
"             smooth, levs\n"
"     (For contours, if you aren't talking about a particular N_CONTOUR,\n"
"      any changes will affect ALL the contours.)\n"
"\n"
"     A plv (vector field) element can also take the scalem\n"
"     keyword to multiply all vector lengths by a specified factor.\n"
"\n"
"     A plt (text) element can also take the dx and/or dy\n"
"     keywords to adjust the text position by (dx,dy).\n"
"\n"
"   SEE ALSO: window, plsys, plq, pldefault, plg\n";

#undef N_KEYWORDS
#define N_KEYWORDS 36
static char *editKeys[N_KEYWORDS+1]= {
  "legend", "hide",
  "color", "type", "width",
  "marks", "mcolor", "marker", "msize", "mspace", "mphase",
  "rays", "arrowl", "arroww", "rspace", "rphase", "closed", "smooth",
  "font", "height", "orient", "justify", "opaque",
  "hollow", "aspect", "region", "boundary", "levs", "scale", "scalem",
  "dx", "dy", "edges", "ecolor", "ewidth", "inhibit", 0 };

static PyObject *pledit (PyObject * self, PyObject * args, PyObject * kd)
{
  int type = 0, n_element = 0, n_contour = 0;
  int changes = 0, resetLevs = 0; 
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
      if (type != E_CONTOURS)  {
	return ERRSS ("current graphical element is not contours in pledit");
      }
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
    if (type == 0)  {
      return ERRSS ("no such graphical element for pledit");
    }
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
      return ERRSS ("orient= keyword must be 0, 1, 2, or 3");
    }
  }

  SETKW(kwt[21], dummy,            setkw_justify,  editKeys[21]);
  SETKW(kwt[22], gistA.t.opaque,   setkw_boolean,  editKeys[22]);
  SETKW(kwt[23], gistA.vect.hollow, setkw_boolean,  editKeys[23]);
  SETKW(kwt[24], gistA.vect.aspect, setkw_double,   editKeys[24]);
  
  if (kwt[25]) {	/* region */
    if (type < 4 || type > 7)  {
      return ERRSS ("region = in pledit allowed only for plm, plf, plv, plc");
    }
    SETKW(kwt[25],  gistD.region,   setkw_integer,  editKeys[25]);
  }
  if (kwt[26]) {	/* boundary */
    if (type != 4)  {
      return ERRSS ("boundary = in pledit allowed only for plm");
    }
    SETKW(kwt[26],  gistD.boundary, setkw_boolean,  editKeys[26]);
  }

  if (kwt[27]) {	/* levs */
    double *levels;
    long nLevels = 0;
    PyArrayObject *lap;
    double *lev;
    int i;

    if (type != 7)  {
      return ERRSS ("levs = in pledit allowed only for plc");
    }

    GET_ARR (lap, kwt[27], PyArray_DOUBLE, 1, PyObject *);
    lev = (double *) A_DATA (lap);
    nLevels = A_SIZE (lap);
    if (0 == nLevels) {
      clearArrayList ();
      return ERRSS ("pledit cannot recompute default contour levels");
      }
    levels = malloc (sizeof(double) * nLevels);
    if(!levels) return PyErr_NoMemory();
    for(i = 0; i < nLevels; i++)
      levels[i] = lev[i];
    removeFromArrayList ( (PyObject *) lap);
    /* WARNING --
       this is a critical code section, since until GdEdit successfully
       completes, Gist owns a pointer to the freed levels -- no way to
       gracefully avoid this without "knowing" more about guts of Gist's
       data structures than seem reasonable here... */
    p_free (gistD.levels);
    gistD.levels = levels;
    gistD.nLevels = nLevels;
    changes |= CHANGE_Z;
    resetLevs = 1;
  }
  if (kwt[28]) {	/* scale */
    if (type != 6)  {
      return ERRSS ("scale = in pledit allowed only for plv");
    }
    SETKW(kwt[28],  gistD.scale,  setkw_double,  editKeys[28]);
  }
  if (kwt[29]) {	/* scalem */
    double scalem;
    if (type != 6)  {
      return ERRSS ("scalem = in pledit allowed only for plv");
    }
    SETKW(kwt[29],  scalem,       setkw_double,  editKeys[29]);
    PyFPE_START_PROTECT("pledit", return 0)
    gistD.scale *= scalem;
    PyFPE_END_PROTECT(dummy)
  }
  if (kwt[30]) {	/* dx */
    double x0;
    if (type != 3)  {
      return ERRSS ("dx = in pledit allowed only for plt");
    }
    SETKW(kwt[30],  x0,           setkw_double,  editKeys[30]);
    PyFPE_START_PROTECT("pledit", return 0)
    gistD.x0 += x0;
    PyFPE_END_PROTECT(dummy)
  }
  if (kwt[31]) {	/* dy */
    double y0;
    if (type != 3)  {
      return ERRSS ("dy = in pledit allowed only for plt");
    }
    SETKW(kwt[31],  y0,           setkw_double,  editKeys[31]);
    PyFPE_START_PROTECT("pledit", return 0)
    gistD.y0 += y0;
    PyFPE_END_PROTECT(dummy)
  }
  if (kwt[32]) {
    int edgetype = 0;
    SETKW(kwt[32],  edgetype,     setkw_boolean, editKeys[32]);
    gistA.e.type = edgetype ? L_SOLID : L_NONE;
  }
  SETKW(kwt[33],  gistA.e.color,  setkw_color,   editKeys[33]);
  SETKW(kwt[34],  gistA.e.width,  setkw_double,  editKeys[34]);

  if (kwt[35]) {	/* inhibit */
    if (type != 4)  {
      return ERRSS ("inhibit = in pledit allowed only for plm");
    }
    SETKW(kwt[35],  gistD.inhibit, setkw_integer, editKeys[35]);
  }
  if (legend) {
    /* Some jiggery-pokery necessary to get the old legend deleted properly,
       and the new legend allocated properly, so that Gist will delete it
       correctly when the graphical element is deleted.  */
    char *oldleg = gistD.legend;
    if (!(gistD.legend = p_malloc (strlen (legend) + 1)))
       return PyErr_NoMemory();
    strcpy (gistD.legend, legend);
    legend = oldleg;
  }
  GdEdit (changes);
  if (legend)
    p_free (legend);
  if ( resetLevs )  {
     tmpLevels = 0;
  }
  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  plf */

static char plf__doc__[] =
"plf( z, y, x )\n"
"or plf( z, y, x, ireg )\n"
"or plf( z )\n"
"     Plot a filled mesh Y versus X.  Y, X, and IREG are as for plm.\n"
"     The Z array must have the same shape as Y and X, or one smaller\n"
"     in both dimensions.  If Z is of type char, it is used `as is',\n"
"     otherwise it is linearly scaled to fill the current palette, as\n"
"     with the bytscl function.\n"
"     (See the bytscl function for explanation of top, cmin, cmax.)\n"
"     The mesh is drawn with each zone in the color derived from the Z\n"
"     function and the current palette; thus Z is interpreted as a\n"
"     zone-centered array.\n"
"     The Y, X, and IREG arguments may all be omitted to default to the\n"
"     mesh set by the most recent plmesh call.\n"
"     A solid edge can optionally be drawn around each zone by setting\n"
"     the EDGES keyword non-zero.  ECOLOR and EWIDTH determine the edge\n"
"     color and width.  The mesh is drawn zone by zone in order from\n"
"     IREG(2+imax) to IREG(jmax*imax) (the latter is IREG(imax,jmax)),\n"
"     so you can achieve 3D effects by arranging for this order to\n"
"     coincide with back-to-front order.  If Z is nil, the mesh zones\n"
"     are filled with the background color, which you can use to\n"
"     produce 3D wire frames.\n"
"     The following keywords are legal (each has a separate help entry):\n"
"\n"
"   KEYWORDS: legend, hide\n"
"             region, top, cmin, cmax, edges, ecolor, ewidth\n"
"\n"
"   SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp, plmesh,\n"
"             limits, logxy, ylimits, fma, hcp, palette, bytscl, histeq_scale\n";

#undef N_KEYWORDS
#define N_KEYWORDS 9
static char *plfKeys[N_KEYWORDS+1]= {
  "legend", "hide", "region", "top", "cmin", "cmax",
  "edges", "ecolor", "ewidth", 0 };

static PyObject *plf (PyObject * self, PyObject * args, PyObject * kd)
{
  PyArrayObject *zap;
  PyObject *zop = 0;
  char *z_name= 0, *y_name= 0, *x_name= 0, *r_name= 0;
  long iMax= 0, jMax= 0;
  double *z = 0;
  GpColor *zc = 0;
  GpColor *zc1 = 0;
  GaQuadMesh mesh;
  int convertedZ= 0;
  int rgb = 0;

  PyObject * kwt[NELT(plfKeys) - 1];
  char *errstr = "plf requires 2D arguments (z [ , y, x, ireg ] )";

  SETJMP0;			/* See Xerror_longjmp() */

  if (PyTuple_Size (args) == 0)  {
    return ERRSS ("plf requires at least one argument");
  }
  BUILD_KWT(kd, plfKeys, kwt);
  TRY (setz_mesh (args, &zop, errstr, 0), (PyObject *) NULL);
  if (!pyMsh.y)  {
    return ERRSS ("No current mesh - set (y, x) first");
  }

  get_mesh (&mesh);

  /*
   *  The first arg to plf, z, is of type unsigned char.
   *  The array size is (M-1,N-1) or
   *  (3,M-1,N-1), giving an (r,g,b) for each true color value.
   */

  if (isARRAY(zop) && (A_TYPE(zop) == Py_GpColor)) {

    if ( A_NDIM(zop) == 2 )  {
       /*  NXxNY */
       GET_ARR (zap, zop, Py_GpColor, 2, PyObject *);
       zc = (GpColor *) A_DATA (zap);
    }
    else if ( A_NDIM(zop) == 3 )  {
       /*  3xNXxNY */
       if ( A_DIM(zop,0) != 3 )  {
          return ERRSS ("expecting NXxNY or 3xNXxNY array as argument to plf");
       }
       GET_ARR (zap, zop, Py_GpColor, 3, PyObject *);
       zc = (GpColor *) A_DATA (zap);
       rgb = 1;
    }
    else  {
       return ERRSS ("expecting NXxNY or 3xNXxNY array as argument to plf"); 
    }

  } else {
    if (isARRAY(zop) && (A_TYPE(zop) == PyArray_DOUBLE)) {
      GET_ARR (zap, zop, PyArray_DOUBLE, 2, PyObject *);
      z = (double *) A_DATA (zap);
    } else {
      z = 0;
      zc = 0;
      zap = 0;
    }
  }

  if (zap) {
    jMax = A_DIM(zap, 0);
    iMax = A_DIM(zap, 1);
  } else {
    jMax = iMax = 0;
  }
  if ((z || zc) && ((mesh.iMax != iMax   || mesh.jMax != jMax) &&
		    (mesh.iMax != iMax+1 || mesh.jMax != jMax+1))) {
    removeFromArrayList ( (PyObject *) zap);
    return ERRSS (
      "z array must have same or 1 smaller dimensions as mesh in plf");
  }

  TRYS(CheckDefaultWindow())
  CheckDefaultPalette ();

  if ( !LegendAndHide("plf, ", z_name, y_name, x_name, r_name, kwt, plfKeys) )
     return ERRSS ( "Error in plf: LegendAndHide" );

  gistD.region = 0;
  SETKW(kwt[2],  gistD.region,    setkw_integer,  plfKeys[2]);

  if (!zc && z) {
    double zmin, zmax, scale, offset;

    TRY (GrabByteScale(&kwt[3], &plfKeys[3], &scale, &offset, &zmin, &zmax,
       z, mesh.reg, gistD.region, mesh.iMax, mesh.jMax,
       (int) (mesh.iMax != iMax) ), (PyObject *) NULL);
    TRY (zc = PushColors(z, iMax*jMax, zmin, zmax, scale, offset),
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
  gistA.rgb = rgb;

/*
 *  LLC:  For some reason, the following yields unaligned accesses on the Alpha.
  if (mesh.iMax==iMax) zc += rgb? 3*(iMax+1) : iMax+1;
 *  Use the following instead.
 */
  if (mesh.iMax==iMax) zc1 = zc + ( rgb? 3*(iMax+1) : iMax+1 );
  curElement = -1;
  PyFPE_START_PROTECT("plf", return 0)
  curElement = GdFillMesh(NOCOPY_MESH, &mesh, gistD.region, zc1, iMax);
  PyFPE_END_PROTECT(dummy)
  clearArrayList ();
  if (convertedZ && zc) free (zc);
  if (curElement < 0)  {
    return ERRSS ("Gist GdFillMesh plotter failed");
  }

  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  plfp */

static char plfp__doc__[] =
"plfp( z, y, x, n )\n"
"     Plot a list of filled polygons Y versus X, with colors Z.\n"
"     The N array is a 1D list of lengths (number of corners) of the\n"
"     polygons; the 1D colors array Z has the same length as N.  The\n"
"     X and Y arrays have length equal to the sum of all dimensions\n"
"     of N.\n"
"     The Z array must have the same shape as Y and X.  If Z is of\n"
"     type char, it is used `as is', otherwise it is linearly scaled\n"
"     to fill the current palette, as with the bytscl function.\n"
"     (See the bytscl function for explanation of top, cmin, cmax.)\n"
"     The following keywords are legal (each has a separate help entry):\n"
"\n"
"   KEYWORDS: legend, hide, top, cmin, cmax\n"
"\n"
"   SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj\n"
"             limits, logxy, ylimits, fma, hcp\n";

#undef N_KEYWORDS
#define N_KEYWORDS 8
static char *plfpKeys[N_KEYWORDS+1]= {
  "legend", "hide", "top", "cmin", "cmax", "edges", "ecolor", "ewidth", 0 };

static PyObject *plfp (PyObject * self, PyObject * args, PyObject * kd)
{
  PyArrayObject *zap = 0, *yap, *xap, *nap;
  PyObject *zop, *yop, *xop, *nop;
  int i;
  long nz, nx, nn, np;
  long ny =0, *pn= 0;
  double *z = 0, *x, *y;
  GpColor *zc = 0;
  int convertedZ= 0;
  int rgb = 0;

  PyObject * kwt[NELT(plfpKeys) - 1];
  char *errstr = "plfp requires arguments (z, y, x, n)";

  SETJMP0;			/* See Xerror_longjmp() */

  if (!PyArg_ParseTuple (args, "OOOO", &zop, &yop, &xop, &nop))  {
    return ERRSS (errstr);
  }

  /*
   *  The first arg to plfp, z, is of type unsigned char.
   *  The array size is N or 3xN, giving an (r,g,b) for each 
   *  true color value.
   */

  if (isARRAY(zop) && (A_TYPE(zop) == Py_GpColor)) {

    if ( A_NDIM(zop) == 1 )  {
       /*  N */
       GET_ARR (zap, zop, Py_GpColor, 1, PyObject *);
       zc = (GpColor *) A_DATA (zap);
    }
    else if ( A_NDIM(zop) == 2 )  {
       /*  3xN */
       if ( A_DIM(zop,0) != 3 )  {
          return ERRSS ("expecting N or 3xN array as argument to plfp");
       }
       GET_ARR (zap, zop, Py_GpColor, 2, PyObject *);
       zc = (GpColor *) A_DATA (zap);
       rgb = 1;
    }
    else  {
       return ERRSS ("expecting N or 3xN array as argument to plfp"); 
    }

  } else if (isARRAY(zop) && (A_TYPE(zop) == PyArray_DOUBLE)) {
    GET_ARR (zap, zop, PyArray_DOUBLE, 1, PyObject *);
    z = (double *) A_DATA (zap);
  }

  GET_ARR (yap, yop, PyArray_DOUBLE, 1, PyObject *);
  GET_ARR (xap, xop, PyArray_DOUBLE, 1, PyObject *);
  GET_ARR (nap, nop, PyArray_LONG, 1, PyObject *);
  nn = A_SIZE (nap);
  nx = A_SIZE (xap);
  ny = A_SIZE (yap);
  nz = (zap) ? A_SIZE (zap) : nn;
  y = (double *) A_DATA (yap);
  x = (double *) A_DATA (xap);
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

  TRYS(CheckDefaultWindow())
  CheckDefaultPalette ();
  /* would need to add plfp to quine list with YpQuine to get legend
     LegendAndHide("plfp, ", z_name, y_name, x_name, r_name, kwt, plfpKeys); */
  if ( !LegendAndHide((char *)0, (char *)0, (char *)0,
                (char *)0, (char *)0, kwt, plfpKeys) )
     return ERRSS ( "Error in plfp: LegendAndHide" );

  if (!zc && z) {
    double zmin, zmax, scale, offset;

    TRY (GrabByteScale(&kwt[2], &plfpKeys[2], &scale, &offset, &zmin, &zmax,
       z, (int *)0, 0, nz + 1, 2L, 1), (PyObject *) NULL);
    TRY (zc = PushColors(z, nz, zmin, zmax, scale, offset), (PyObject *) NULL);
    convertedZ= 1;
  }
  else {
    convertedZ= 0;
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
  gistA.rgb = rgb;

  curElement = -1;
  PyFPE_START_PROTECT("plfp", return 0)
  curElement = GdFill (nz, zc, x, y, pn);
  PyFPE_END_PROTECT(dummy)
  clearArrayList ();
  if (convertedZ) free (zc);
  if (curElement < 0)  {
    return ERRSS ("Gist GdFill plotter failed");
  }

  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  plg */

static char plg__doc__[] =
"plg( y [, x] )\n"
"     Plot a graph of Y versus X.  Y and X must be 1-D arrays of equal\n"
"     length; if X is omitted, it defaults to [1, 2, ..., 1+range(len(Y))].\n"
"     The following keywords are legal (each has a separate help entry):\n"
"\n"
"   KEYWORDS: legend, hide\n"
"             type, width, color, closed, smooth\n"
"             marks, marker, mspace, mphase\n"
"             rays, arrowl, arroww, rspace, rphase\n"
"\n"
"   Example:    plg ( y, x, type=0, marker=character )\n"
"\n"
"   If character is '\\1', '\\2', '\\3', '\\4', or '\\5', you get point, plus,\n"
"   asterisk, circle, and cross, respectively.  If you make the marker size\n"
"   small enough (the default is small enough), then '\\1' will plot points\n"
"   on the X display, which is the most usual request.  For 2-5 or for large\n"
"   marker size, you get the characters on the X display, but in the\n"
"   hardcopy files (postscript or cgm) those special markers will be rendered\n"
"   nicely.\n"
"\n"
"   SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp\n"
"             limits, logxy, ylimits, fma, hcp\n";

#undef N_KEYWORDS
#define N_KEYWORDS 19
static char *plgKeys[N_KEYWORDS+1]= {
  "legend", "hide", "color", "type", "width",
  "marks", "mcolor", "marker", "msize", "mspace", "mphase",
  "rays", "arrowl", "arroww", "rspace", "rphase",
  "closed", "smooth", "n", 0 };

static PyObject *plg (PyObject * self, PyObject * args, PyObject * kd)
{
  PyObject *xop = 0, *yop;
  PyArrayObject *xap, *yap;
  double *x= 0, *y= 0;

  int i;
  long length;
  PyObject * kwt[NELT(plgKeys) - 1];
  char *errstr =
    "plg requires one or two 1-D double arrays, of the same length";

  SETJMP0;			/* See Xerror_longjmp() */

  if (!PyArg_ParseTuple (args, "O|O", &yop, &xop)) {
    return ERRSS (errstr);
  }
  GET_ARR(yap, yop, PyArray_DOUBLE, 1, PyObject *);
  length = A_SIZE(yap);
  y = (double *) A_DATA(yap);

  TRYS(CheckDefaultWindow())

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
    GET_ARR(xap, xop, PyArray_DOUBLE, 1, PyObject *);
    if(A_SIZE(xap) != length) {
      clearArrayList ();
      return ERRSS (errstr);
    }
    x = (double *) A_DATA(xap);
  } else {
    NEW_MEM (x, length, double, PyObject *);
    for (i = 0; i < length; i++)
	x[i] = (double) (1+i);
  }

  curElement = -1;
  PyFPE_START_PROTECT("plg", return 0)
  curElement = GdLines (length, x, y);
  PyFPE_END_PROTECT(dummy)

  clearArrayList ();
  clearMemList ();

  if (curElement < 0)  {
    return ERRSS ("Gist GdLines plotter failed");
  }

  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  pli */

static char pli__doc__[] =
"pli( z )\n"
"or pli( z, x1, y1 )\n"
"or pli( z, x0, y0, x1, y1 )\n"
"     Plot the image Z as a cell array -- an array of equal rectangular\n"
"     cells colored according to the 2-D array Z.  The first dimension\n"
"     of Z is plotted along x, the second dimension is along y.\n"
"     If Z is of type char, it is used `as is', otherwise it is linearly\n"
"     scaled to fill the current palette, as with the bytscl function.\n"
"     (See the bytscl function for explanation of top, cmin, cmax.)\n"
"     \n"
"     As for plf and plfp, Z may also be a 3D array with 1st dimension 3\n"
"     of char giving the [r,g,b] components of each color.  See the\n"
"     color keyword for cautions about using this if you do not have\n"
"     a true color display.\n"
"     \n"
"     If X1 and Y1 are given, they represent the coordinates of the\n"
"     upper right corner of the image.  If X0, and Y0 are given, they\n"
"     represent the coordinates of the lower left corner, which is at\n"
"     (0,0) by default.  If only the Z array is given, each cell will be\n"
"     a 1x1 unit square, with the lower left corner of the image at (0,0).\n"
"     The following keywords are legal (each has a separate help entry):\n"
"\n"
"   KEYWORDS: legend, hide, top, cmin, cmax\n"
"\n"
"   SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp,\n"
"             limits, logxy, ylimits, fma, hcp, palette, bytscl, histeq_scale\n";

#undef N_KEYWORDS
#define N_KEYWORDS 5
static char *pliKeys[N_KEYWORDS+1]= {
  "legend", "hide", "top", "cmin", "cmax", 0 };

static PyObject *pli (PyObject * self, PyObject * args, PyObject * kd)
{
  PyArrayObject *zap;
  PyObject *zop = 0;
  char *z_name= 0;
  int nargs;
  double *z = 0, x0, y0, x1, y1;
  long iMax= 0, jMax= 0;
  GpColor *zc = 0;
  int convertedZ= 0;
  int rgb = 0;

  PyObject * kwt[NELT(pliKeys) - 1];
  char *errstr = "pli requires arguments (z [ , [ x0, y0, ] x1, y1 ] )";

  SETJMP0;			/* See Xerror_longjmp() */

  x0= y0= x1= y1= 0.0;

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

  /*  
   *  The first arg to pli, z, is of type unsigned char.
   *  The array size is (M-1,N-1) or 
   *  (3,M-1,N-1), giving an (r,g,b) for each true color value.
   */

  if (isARRAY(zop) && (A_TYPE(zop) == Py_GpColor)) {

    if ( A_NDIM(zop) == 2 )  {
       /*  NXxNY */
       GET_ARR (zap, zop, Py_GpColor, 2, PyObject *);
       zc = (GpColor *) A_DATA (zap);
    }
    else if ( A_NDIM(zop) == 3 )  {
       /*  3xNXxNY */
       if ( A_DIM(zop,0) != 3 )  {
          return ERRSS ("expecting NXxNY or 3xNXxNY array as argument to pli"); 
       }
       GET_ARR (zap, zop, Py_GpColor, 3, PyObject *);
       zc = (GpColor *) A_DATA (zap);
       rgb = 1;
    }
    else  {
       return ERRSS ("expecting NXxNY or 3xNXxNY array as argument to pli"); 
    }

  } else {
    GET_ARR (zap, zop, PyArray_DOUBLE, 2, PyObject *);
    z = (double *) A_DATA (zap);
  }

  iMax = A_DIM(zap, 1);
  jMax = A_DIM(zap, 0);

  if (1 == nargs) {
    x1 = (double) iMax;
    y1 = (double) jMax;
  }

  if (!z && !zc) return ERRSS ("pli needs at least one non-keyword argument");

  BUILD_KWT(kd, pliKeys, kwt);

  TRYS(CheckDefaultWindow())
  CheckDefaultPalette ();
  if ( !LegendAndHide("pli, ", z_name, (char *)0,(char *)0,(char *)0, kwt, pliKeys) )
     return ERRSS ( "Error in pli: LegendAndHide" );

  if (!zc) {
    /* need to generate colors array on stack now */
    double zmin, zmax, scale, offset;

    TRY (GrabByteScale(&kwt[2], &pliKeys[2], &scale, &offset, &zmin, &zmax,
       z, (int *)0, 0, iMax + 1, jMax + 1, 1), 
       (PyObject *) NULL);
    TRY (zc = PushColors(z, iMax*jMax, zmin, zmax, scale, offset), 
       (PyObject *) NULL);
    convertedZ= 1;
  }

  gistA.rgb = rgb;

  SETKW(kwt[0],  gistD.legend,    setkw_string,   pliKeys[0]);
  SETKW(kwt[1],  gistD.hidden,    setkw_boolean,  pliKeys[1]);

  curElement = -1;
  PyFPE_START_PROTECT("pli", return 0)
  curElement = GdCells (x0, y0, x1, y1, iMax, jMax, iMax, zc);
  PyFPE_END_PROTECT(dummy)
  removeFromArrayList ( (PyObject *) zap);
  if (convertedZ) free (zc);
  if (curElement < 0)  {
    return ERRSS ("Gist GdCells plotter failed");
  }

  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  plm */

static char plm__doc__[] =
"plm( y, x, boundary=0/1, inhibit=0/1/2 )\n"
"or plm( y, x, ireg, boundary=0/1, inhibit=0/1/2 )\n"
"or plm( boundary=0/1, inhibit=0/1/2 )\n"
"     Plot a mesh of Y versus X.  Y and X must be 2-D arrays with equal\n"
"     dimensions.  If present, IREG must be a 2-D region number array\n"
"     for the mesh, with the same dimensions as X and Y.  The values of\n"
"     IREG should be positive region numbers, and zero for zones which do\n"
"     not exist.  The first row and column of IREG never correspond to any\n"
"     zone, and should always be zero.  The default IREG is 1 everywhere\n"
"     else.  If present, the BOUNDARY keyword determines whether the\n"
"     entire mesh is to be plotted (boundary=0, the default), or just the\n"
"     boundary of the selected region (boundary=1).  If present, the\n"
"     INHIBIT keyword causes the (X(,j),Y(,j)) lines to not be plotted\n"
"     (inhibit=1), or the (X(i,),Y(i,)) lines to not be plotted (inhibit=2).\n"
"     By default (inhibit=0), mesh lines in both logical directions are\n"
"     plotted.\n"
"     The Y, X, and IREG arguments may all be omitted to default to the\n"
"     mesh set by the most recent plmesh call.\n"
"     The following keywords are legal (each has a separate help entry):\n"
"\n"
"   KEYWORDS: legend, hide\n"
"             type, width, color\n"
"             region\n"
"\n"
"   SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp, plmesh\n"
"             limits, logxy, ylimits, fma, hcp\n";

#undef N_KEYWORDS
#define N_KEYWORDS 8
static char *plmKeys[N_KEYWORDS+1]= {
  "legend", "hide", "color", "type", "width", "region", "boundary",
  "inhibit", 0 };

static PyObject *plm (PyObject * self, PyObject * args, PyObject * kd)
{
  char *y_name= 0, *x_name= 0, *r_name= 0;
  GaQuadMesh mesh;
  PyObject *kwt[NELT(plmKeys) - 1];
  char *errstr = "plm takes 1-3 non-keyword arguments: (y, x, ireg).";

  SETJMP0;
  if (PyTuple_Size (args) > 0)
    TRY (set_pyMsh (args, errstr, 0), 
       (PyObject *) NULL);

  get_mesh(&mesh);

  BUILD_KWT(kd, plmKeys, kwt);

  /* set legend and hide in gistD */
  TRYS(CheckDefaultWindow())
  if ( !LegendAndHide("plm, ", y_name, x_name, r_name, (char *)0, kwt, plmKeys) )
     return ERRSS ( "Error in plm: LegendAndHide" );

  /* set properties, starting from defaults for meshes */
  GhGetMesh();
  gistD.region = 0;
  gistD.boundary = 0;
  gistD.inhibit = 0;

  SETKW(kwt[0],  gistD.legend,    setkw_string,   plmKeys[0]);
  SETKW(kwt[1],  gistD.hidden,    setkw_boolean,  plmKeys[1]);
  SETKW(kwt[2],  gistA.l.color,   setkw_color,    plmKeys[2]);
  SETKW(kwt[3],  gistA.l.type,    setkw_linetype, plmKeys[3]);
  SETKW(kwt[4],  gistA.l.width,   setkw_double,   plmKeys[4]);
  SETKW(kwt[5],  gistD.region,    setkw_integer,  plmKeys[5]);
  SETKW(kwt[6],  gistD.boundary,  setkw_boolean,  plmKeys[6]);
  SETKW(kwt[7],  gistD.inhibit,   setkw_integer,  plmKeys[7]);

  if (!pyMsh.y) {
    return ERRSS ("no current mesh - use plmesh(y, x) to initialize");
  }

  TRYS(CheckDefaultWindow())
  curElement = -1;
  PyFPE_START_PROTECT("plm", return 0)
  curElement = GdMesh(NOCOPY_MESH, &mesh, gistD.region, gistD.boundary,
		     gistD.inhibit);
  PyFPE_END_PROTECT(dummy)

  if (curElement < 0)  {
    return ERRSS ("Gist GdMesh plotter failed");
  }

  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  plmesh */

static char plmesh__doc__[] =
"plmesh( y, x, ireg, triangle=tri_array )\n"
"or plmesh()\n"
"     Set the default mesh for subsequent plm, plc, plv, and plf calls.\n"
"     In the second form, deletes the default mesh (until you do this,\n"
"     or switch to a new default mesh, the default mesh arrays persist and\n"
"     take up space in memory).  The Y, X, and IREG arrays should all be\n"
"     the same shape; Y and X will be converted to double, and IREG will\n"
"     be converted to int.  If IREG is omitted, it defaults to IREG(1,)=\n"
"     IREG(,1)= 0, IREG(2:,2:)=1; that is, region number 1 is the whole\n"
"     mesh.  The triangulation array TRI_ARRAY is used by plc; the\n"
"     correspondence between TRI_ARRAY indices and zone indices is the\n"
"     same as for IREG, and its default value is all zero.\n"
"     The IREG or TRI_ARRAY arguments may be supplied without Y and X\n"
"     to change the region numbering or triangulation for a given set of\n"
"     mesh coordinates.  However, a default Y and X must already have been\n"
"     defined if you do this.\n"
"     If Y is supplied, X must be supplied, and vice-versa.\n"
"\n"
"   SEE ALSO: plm, plc, plv, plf, plfp\n";

#undef N_KEYWORDS
#define N_KEYWORDS 1
static char *meshKeys[N_KEYWORDS+1]= { "triangle", 0 };

static PyObject *plmesh (PyObject * self, PyObject * args, PyObject * kd)
{
  PyObject * kwt[NELT(meshKeys) - 1];
  char *errstr = "plmesh takes 0-3 non-keyword arguments: (y, x, ireg).";

  BUILD_KWT(kd, meshKeys, kwt);
  TRY (set_pyMsh (args, errstr, kwt[0]), 
     (PyObject *) NULL);

  Py_INCREF (Py_None);
  return Py_None;
}

/*  -------------------------------------------------------------------- */

static char plq__doc__[] =
"plq()\n"
"or plq( n_element )\n"
"or plq( n_element, n_contour )\n"
"or legend_list = plq() **** RETURN VALUE NOT YET IMPLEMENTED ****\n"
"or properties = plq(n_element, n_contour)\n"
"     Called as a subroutine, prints the list of legends for the current\n"
"     coordinate system (with an `(H)' to mark hidden elements), or prints\n"
"     a list of current properties of element N_ELEMENT (such as line type,\n"
"     width, font, etc.), or of contour number N_CONTOUR of element number\n"
"     N_ELEMENT (which must be contours generated using the plc command).\n"
"     Elements and contours are both numbered starting with one; hidden\n"
"     elements or contours are included in this numbering.\n"
"\n"
"     The plq function always operates on the current coordinate system\n"
"     in the current graphics window; use window and plsys to change these.\n"
"\n"
"   SEE ALSO: window, plsys, pledit, pldefault, plg\n";

static PyObject *plq (PyObject * self, PyObject * args)
{
  int type, n_element = 0, n_contour = 0;

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
      if (type != E_CONTOURS)  {
	return ERRSS ("current graphical element is not contours in pledit");
      }
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
      PrintInit (pyg_puts);

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
      PrintInit (pyg_puts);
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

/*  -------------------------------------------------------------------- */

static char plremove__doc__[] =
"plremove( n_element )\n"
"     Removes the element number N_ELEMENT from the display.\n"
"     If N_ELEMENT is omitted, the default is the most recently added\n"
"     element, or the element specified in the most recent plq query\n"
"     command.\n";

static PyObject *plremove(PyObject * self, PyObject * args)
{
  int type = 0, n_element = 0,result;

  switch (PyTuple_Size (args)) {
  case 1: /* (n_element) given */
    TRY (PyArg_ParseTuple (args, "i", &n_element), (PyObject *) NULL);
    break;
  case 0: /* () given */
    break;
  default:
    return ERRSS ("plremove function takes no more than one argument");
  }

  /* Pygist uses 1-origin element numbering, Gist uses 0-origin */
  n_element--;

  if (n_element < 0) {
    if (curElement >= 0) {
      n_element = GdFindIndex (curElement);
      if (n_element < 0) {
	curElement = -1;
	return ERRSS ("lost current graphical element for plremove (BUG?)");
      }
    } else if (curElement == -6666) {
      n_element = curIX;
    } else {
      return ERRSS ("no current graphical element for plremove");
    }
  }
  if (n_element >= 0) {
    /* retrieve specified element */
    type = GdSetElement (n_element);
  }
  curElement = -1;

  PyFPE_START_PROTECT("plremove", return 0)
  result = GdRemove();
  PyFPE_END_PROTECT(dummy)

  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

/*  -------------------------------------------------------------------- */

static char plsys__doc__[] =
"plsys( n )\n"
"     Set the current coordinate system to number N in the current\n"
"     graphics window.  If N equals 0, subsequent elements will be\n"
"     plotted in absolute NDC coordinates outside of any coordinate\n"
"     system.  The default style sheet `work.gs' defines only a single\n"
"     coordinate system, so the only other choice is N equal 1.  You\n"
"     can make up your own style sheet (using a text editor) which\n"
"     defines multiple coordinate systems.  You need to do this if\n"
"     you want to display four plots side by side on a single page,\n"
"     for example.  The standard style sheets `work2.gs' and `boxed2.gs'\n"
"     define two overlayed coordinate systems with the first labeled\n"
"     to the right of the plot and the second labeled to the left of\n"
"     the plot.  When using overlayed coordinate systems, it is your\n"
"     responsibility to ensure that the x-axis limits in the two\n"
"     systems are identical.\n"
"\n"
"   SEE ALSO: window, limits, plg\n";

static PyObject *plsys (PyObject * self, PyObject * args)
{
  int n = -9999, n0;
  char *errstr = "Error: plsys takes zero or one integer argument.";

  SETJMP0;
  if (!PyArg_ParseTuple (args, "|i", &n)) {
    return ERRSS ( errstr );
  }

  TRYS(CheckDefaultWindow())
  n0 = GdGetSystem();

  if (n != -9999){
    if (GdSetSystem (n) != E_SYSTEM && n != 0) {
      return ERRSS (
       "No such coordinate system exists in current graphics window.");
    }
  }
  return Py_BuildValue ("i",n0);
}

/*  -------------------------------------------------------------------- */
/*  plt */

static char plt__doc__[] =
"plt( text, x, y, tosys=0/1 )\n"
"     Plot TEXT (a string) at the point (X,Y).  The exact relationship\n"
"     between the point (X,Y) and the TEXT is determined by the\n"
"     justify keyword.  TEXT may contain newline (`\n') characters\n"
"     to output multiple lines of text with a single call.  The\n"
"     coordinates (X,Y) are NDC coordinates (outside of any coordinate\n"
"     system) unless the tosys keyword is present and non-zero, in\n"
"     which case the TEXT will be placed in the current coordinate\n"
"     system.  However, the character height is NEVER affected by the\n"
"     scale of the coordinate system to which the text belongs.\n"
"     Note that the pledit command takes dx and/or dy keywords to\n"
"     adjust the position of existing text elements.\n"
"     The following keywords are legal (each has a separate help entry):\n"
"\n"
"   KEYWORDS: legend, hide\n"
"             color, font, height, opaque, path, justify\n"
"\n"
"   SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp, pledit\n"
"             limits, ylimits, fma, hcp, pltitle\n";

#undef N_KEYWORDS
#define N_KEYWORDS 9
static char *pltKeys[N_KEYWORDS+1]= {
  "legend", "hide",
  "color", "font", "height", "orient", "justify", "opaque", "tosys", 0 };

static PyObject *plt (PyObject * self, PyObject * args, PyObject * kd)
{
  char *text = 0;
  double x = 0.0, y = 0.0;
  int toSys = 0;
  PyObject *kwt[NELT (pltKeys) - 1];

  SETJMP0;			/* See Xerror_longjmp() */

  if (!PyArg_ParseTuple (args, "sdd", &text, &x, &y)) {
    return ERRSS ("plt requires exactly three non-keyword arguments");
  }

  BUILD_KWT(kd, pltKeys, kwt);

  /* set legend and hide in gistD */
  TRYS(CheckDefaultWindow())
  if ( !LegendAndHide((char *)0, (char *)0, (char *)0,
                (char *)0, (char *)0, kwt, pltKeys) )
     return ERRSS ( "Error in plt: LegendAndHide" );

  /* set properties, starting from defaults for vectors */
  GhGetText ();
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
      return ERRSS ("orient= keyword must be 0, 1, 2, or 3");
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
  PyFPE_END_PROTECT(dummy)
  if (curElement < 0)  {
    return ERRSS ("Gist GdText plotter failed");
  }

  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

/*  -------------------------------------------------------------------- */
/*  plv */

static char plv__doc__[] =
"plv( vy, vx, y, x, scale=dt )\n"
"or plv( vy, vx, y, x, ireg, scale=dt )\n"
"or plv( vy, vx, scale=dt )\n"
"     Plot a vector field (VX,VY) on the mesh (X,Y).  Y, X, and IREG are\n"
"     as for plm.  The VY and VX arrays must have the same shape as Y and X.\n"
"     The Y, X, and IREG arguments may all be omitted to default to the\n"
"     mesh set by the most recent plmesh call.\n"
"     The SCALE keyword is the conversion factor from the units of\n"
"     (VX,VY) to the units of (X,Y) -- a time interval if (VX,VY) is a velocity\n"
"     and (X,Y) is a position -- which determines the length of the\n"
"     vector `darts' plotted at the (X,Y) points.  If omitted, SCALE is\n"
"     chosen so that the longest ray arrows have a length comparable\n"
"     to a `typical' zone size.\n"
"     You can use the scalem keyword in pledit to make adjustments to the\n"
"     SCALE factor computed by default.\n"
"     The following keywords are legal (each has a separate help entry):\n"
"\n"
"   KEYWORDS: legend, hide\n"
"             type, width, color, smooth\n"
"             marks, marker, mspace, mphase\n"
"             triangle, region\n"
"\n"
"   SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp, plmesh, pledit,\n"
"             limits, logxy, ylimits, fma, hcp\n";

#undef N_KEYWORDS
#define N_KEYWORDS 8
static char *plvKeys[N_KEYWORDS+1]= {
  "legend", "hide", "region",
  "color", "hollow", "width", "aspect", "scale", 0 };

static PyObject *plv (PyObject * self, PyObject * args, PyObject * kd)
{
  PyArrayObject *uap, *vap;
  PyObject *uop, *vop;
  char *v_name= 0, *u_name= 0, *y_name= 0, *x_name= 0;
  long iMax= 0, jMax= 0;
  double *u= 0, *v= 0, scale;
  GaQuadMesh mesh;

  PyObject * kwt[NELT(plvKeys) - 1];
  char *errstr =
    "plv requires 2D arguments (v, u [ , y, x, ireg, scale = dt ] )";

  SETJMP0;			/* See Xerror_longjmp() */

  BUILD_KWT(kd, plvKeys, kwt);

  TRYS(CheckDefaultWindow())
  if ( !LegendAndHide("plv, ", v_name, u_name, y_name, x_name, kwt, plvKeys) )
     return ERRSS ( "Error in plv: LegendAndHide" );

  if (PyTuple_Size (args) < 2)  {
    return ERRSS ("plv requires at least two arguments");
  }
  TRY (setvu_mesh (args, &vop, &uop, errstr), 
     (PyObject *) NULL);
  if (!pyMsh.y)  {
    return ERRSS ("No current mesh - set (y, x) first");
  }
  GET_ARR (vap, vop, PyArray_DOUBLE, 2, PyObject *);
  GET_ARR (uap, uop, PyArray_DOUBLE, 2, PyObject *);
  jMax = (A_DIM(vap, 0) == A_DIM(uap, 0)) ? A_DIM(vap, 0) : 0;
  iMax = (A_DIM(vap, 1) == A_DIM(uap, 1)) ? A_DIM(vap, 1) : 0;
  if (A_DIM (pyMsh.y, 0) != jMax || A_DIM (pyMsh.y, 1) != iMax ) {
    clearArrayList ();
    return ERRSS ("(v, u) arrays must match (y, x) mesh arrays in shape");
  }
  v = (double *) A_DATA (vap);
  u = (double *) A_DATA (uap);
  get_mesh (&mesh);
  if (mesh.iMax!=iMax || mesh.jMax!=jMax)
    return ERRSS ("v and u arrays must have same dimensions as mesh in plv");

  /* set legend and hide in gistD */
  TRYS(CheckDefaultWindow())
  if ( !LegendAndHide("plv, ", v_name, u_name, y_name, x_name, kwt, plvKeys) )
     return ERRSS ( "Error in plv: LegendAndHide" );

  /* set properties, starting from defaults for vectors */
  GhGetVectors();
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
    GetPCrange(&xmin, &xmax, mesh.x, mesh.reg, gistD.region, iMax, jMax);
    GetPCrange(&ymin, &ymax, mesh.y, mesh.reg, gistD.region, iMax, jMax);
    GetPCrange(&umin, &umax, u, mesh.reg, gistD.region, iMax, jMax);
    GetPCrange(&vmin, &vmax, v, mesh.reg, gistD.region, iMax, jMax);

    umax -= umin;
    vmax -= vmin;
    if (vmax > umax) umax = vmax;
    xmax = (xmax - xmin) + (ymax - ymin);
    xmax /= (iMax + jMax);

    if (umax > 0.0) scale = xmax / umax;
    else scale = 1.0;
    PyFPE_END_PROTECT(dummy)
  }

  curElement = -1;
  PyFPE_START_PROTECT("plv", return 0)
  curElement = GdVectors(NOCOPY_MESH, &mesh, gistD.region, u, v, scale);
  PyFPE_END_PROTECT(dummy)
  clearArrayList ();
  if (curElement < 0)  {
    return ERRSS ("Gist GdVectors plotter failed");
  }
  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
  return Py_None;
}

static long Safe_strlen(const char *s)
{
  if (s) return strlen(s);
  else return 0;
}

static void AllocTmpLegend(long len)
{
  if (tmpLegend) FreeTmpLegend();
  tmpLegend = p_malloc(len+1);
  tmpLegend[0] = '\0';
}

static void FreeTmpLegend(void)
{
  if (tmpLegend) {
    char *legend= tmpLegend;
    p_free(legend);
    tmpLegend= 0;
  }
}

static long escape_count(char *arg)
{
  long n= 0;
  if (arg) while (*arg) {
    if (*arg=='!' || *arg=='_' || *arg=='^') n++;
    arg++;
  }
  return n;
}

static void escape_cat(char *leg, char *arg)
{
  while (*arg) {
    if (*arg=='!' || *arg=='_' || *arg=='^') *(leg++)= '!';
    *(leg++)= *(arg++);
  }
  *leg= '\0';
}

static int LegendAndHide(char *func, char *arg1, char *arg2, char *arg3,
                          char *arg4, PyObject *kwt[], char *keys[])
{
  /* check for hide= keyword */

  gistD.hidden= 0;
  SETKW(kwt[1], gistD.hidden, setkw_boolean, keys[1]);

  if (tmpLegend) FreeTmpLegend();

  /* check for legend= keyword -- put legend into tmpLegend */
  /* legend=[] is same as legend=string() */

  SETKW(kwt[0], tmpLegend, setkw_string, keys[0]);
    
  else if (func) {
    /* construct default legend from up to 4 quined arguments */
    long len0= Safe_strlen(func);
    long len1= Safe_strlen(arg1)+escape_count(arg1);
    long len2= Safe_strlen(arg2)+escape_count(arg2);
    long len3= Safe_strlen(arg3)+escape_count(arg3);
    long len4= Safe_strlen(arg4)+escape_count(arg4);
    AllocTmpLegend(len0+len1+len2+len3+len4+6);
    if (func) strcat(tmpLegend, func);
    if (arg1) {
      escape_cat(tmpLegend+len0, arg1);
      len0+= len1;
      if (arg2) {
        strcat(tmpLegend+len0, ", ");
        escape_cat(tmpLegend+len0+2, arg2);
        len0+= 2+len2;
        if (arg3) {
          strcat(tmpLegend+len0, ", ");
          escape_cat(tmpLegend+len0+2, arg3);
          len0+= 2+len3;
          if (arg4) {
            strcat(tmpLegend+len0, ", ");
            escape_cat(tmpLegend+len0+2, arg4);
            len0+= 2+len4;
          }
        }
      }
    }
  }

  /* Put tmpLegend into gistD.legend -- it will be copied out when the
     element is created.  Only danger is pledit, since GdEdit just
     copies the pointer, not the string -- handle this case specially.  */
  gistD.legend= tmpLegend;
  return 1;
}

#if 0
static void print_array_stats(PyArrayObject *op)
{
  int i,ne;
  double *dp;
  TO_STDERR("Data pointer: %p Base pointer: %p\n", op->data, op->base);
  TO_STDERR("Num dims: %d Flags: %d\n", op->nd, op->flags);
  TO_STDERR("Dims & strides:");
  for(i=0; i<op->nd; i++)
    TO_STDERR(" i: %d dim: %d stride: %d",i,op->dimensions[i], op->strides[i]);
  TO_STDERR("\n");
  ne = op->dimensions[0];
  for(i=1; i<op->nd; i++)
    ne *= op->dimensions[i];
  TO_STDERR("Data: (ne = %d)", ne);
  for(i=0,dp = (double *)op->data; i < ne; i++, dp++)
    TO_STDERR(" %.1g", *dp);
  TO_STDERR("\n\n");
  flush_stderr();
}
#endif

/*  -------------------------------------------------------------------- */

static char redraw__doc__[] =
"redraw()\n"
"     Redraw the X window associated with the current graphics window.\n"
"\n"
"   SEE ALSO: window, fma, hcp, plg\n";

static PyObject *redraw (PyObject * self, PyObject * args)
{
  SETJMP0;
  PyFPE_START_PROTECT("redraw", return 0)
  TRYS(CheckDefaultWindow())
  GhRedraw ();
  PyFPE_END_PROTECT(dummy)
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

  if (!PyArg_ParseTuple (args, "|OOO", &op1, &op2, &op3))  {
    return (int) ERRSS (errstr);
  }

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
  if (!ok) {
     return (int) ERRSS ("(ireg) must be a 2-D int array");
  }

  if (!pyMsh.y)  {
    return (int) ERRSS ("No current mesh - ireg not set - set (y, x) first");
  }
  nr = A_DIM (op, 0);
  nc = A_DIM (op, 1);
  if (A_DIM (pyMsh.y, 0) != nr || A_DIM (pyMsh.y, 1) != nc)  {
    return (int) ERRSS ("(ireg) must match (y, x) in shape");
  }

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

/*  -------------------------------------------------------------------- */

static char set_slice2_position__doc__[] =
"None.";

static PyObject *set_slice2_precision (PyObject * self, PyObject * args)
{
 if ( ! PyArg_ParseTuple (args, "d", &_slice2_precision))  {
    return ERRSS ("set_slice2_precision: bad value.");
 }
 Py_INCREF (Py_None);
 return Py_None;
}

/* Create a triangulation (mesh) array. */
static int set_tri (PyObject *top)
{
  int nr, nc;

  if (!pyMsh.y)  {
    return (int) ERRSS ("No current mesh - triangle not set - set (y, x) first");
  }
  nr = A_DIM (pyMsh.y, 0);
  nc = A_DIM (pyMsh.y, 1);

  Py_XDECREF (pyMsh.triangle);
  GET_ARR (pyMsh.triangle, top, PyArray_SHORT, 2, int);

  if (A_DIM (pyMsh.triangle, 0) != nr || A_DIM (pyMsh.triangle, 1) != nc) {
    removeFromArrayList ((PyObject *)pyMsh.triangle);
    return (int) ERRSS ("triangle array must match shape of (y, x).");
  }
  array_list_length = 0;
  return 1;
}

static int set_yx (PyObject *yop, PyObject *xop)
{
  int nr, nc;

  clear_pyMsh();
  GET_ARR (pyMsh.y, yop, PyArray_DOUBLE, 2, int);
  nr = A_DIM (pyMsh.y, 0);
  nc = A_DIM (pyMsh.y, 1);
  if (nr < 2 || nc < 2) {
    clearArrayList ();
    return (int) ERRSS ("(y, x) arrays must be at least 2X2");
  }
  GET_ARR (pyMsh.x, xop, PyArray_DOUBLE, 2, int);
  if (A_DIM (pyMsh.x, 0) != nr || A_DIM (pyMsh.x, 1) != nc) {
    clearArrayList ();
    return (int) ERRSS ("x array must match shape of y");
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

/* Set value for "color=" keyword.  Value passed can be either a string,
 * an integer, or a triple.  All these setkw_*() functions return 0 on error,
 * non-zero otherwise. */

static int setkw_color (PyObject * v, unsigned long *t, char *kw) 
{
  unsigned long color = P_FG;
  unsigned long colors[3] = { 0, 0, 0 };

  if (PyString_Check (v)) {
    char *s = PyString_AsString (v);
    if (strcmp (s, "bg") == 0)
      color = P_BG;
    else if (strcmp (s, "fg") == 0)
      color = P_FG;
    else if (strcmp (s, "black") == 0)
      color = P_BLACK;
    else if (strcmp (s, "white") == 0)
      color = P_WHITE;
    else if (strcmp (s, "red") == 0)
      color = P_RED;
    else if (strcmp (s, "green") == 0)
      color = P_GREEN;
    else if (strcmp (s, "blue") == 0)
      color = P_BLUE;
    else if (strcmp (s, "cyan") == 0)
      color = P_CYAN;
    else if (strcmp (s, "magenta") == 0)
      color = P_MAGENTA;
    else if (strcmp (s, "yellow") == 0)
      color = P_YELLOW;
    else {
      char errstr[256];
      sprintf (errstr, "Unrecognized color keyword: %s: "
	       "Use fg, bg, or 8 primaries only", s);
      return (int) ERRSS (errstr);
    }
  } else if (PyInt_Check (v)) {
    int color1 = PyInt_AsLong (v);
    if ( color1 < 0 )  {
       color = (color1 & 0xff);  /* take right 8 bits */
    }
    else  {
       color = (unsigned long) color1;
    }
  } 

  /* Handle case of 3 element array for trucolor */

  else if ( v && PyTuple_Check (v) ) {  
     if (!unpack_color_tuple ( v, colors)) {
        return (int) NULL;
     }
     color = P_RGB(colors[0],colors[1],colors[2]);
  }
  else  {
    return (int) ERRSS ("Color keyword value must be string, integer, or a triple (r,g,b)");
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
    if (type < 0) type = 0;
    else if (type>5) type= 1 + (type-1)%5;
  } else {
      return (int) ERRSS (errstr);
  }

  *t = type;
  return 1;
}

static int setkw_string (PyObject * v, char **t, char *kw)
{
  char buf[256];
  char *format = "%s keyword requires string argument";
  char *tmpString;

  /*
   *  09/04/02 llc To address the comment below on the string returned
   *               from PyString_AsString: 
   *               Lee Taylor explained that the string returned is a 
   *               string in a string object (a PyObject), so
   *               the string object, rather than the string, was allocated. 
   *               To be on the same side, create space and 
   *               copy the string out of the Python string object.
   *               Error surfaced in FreeTmpLegend from LegendAndHide
   *               with seg fault on the LX cluster when trying to free 
   *               tmpLegend,
   *               but correct here instead of there to make sure 
   *               the problem does not recur.
   */

  if (PyString_Check (v)) {
    tmpString = PyString_AsString (v);
    if ( tmpString == NULL )  {
       *t = NULL;
    }
    else {
       *t = (char *) malloc ( strlen ( tmpString ) + 1 ); 
       strcpy ( *t, tmpString );
    }
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
    return (int) ERRSS (errstr);
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
    return (int) ERRSS (errstr);
  }
  return 1;
}

static int unpack_color_tuple (PyObject * ob, unsigned long color_triple[3])
{
  int i, size = PyTuple_Size (ob);
  PyObject *item;
  if ( size != 3 ) {
    return (int) ERRSS ("Color tuple must have 3 colors");
  }
  for (i = 0; i < 3; i++) {
    if ((item = PyTuple_GetItem (ob, i)) == 0) {
      return (int) ERRSS ("Error unpacking color tuple.");
    }
    if (PyInt_Check (item)) {
      color_triple[i] = PyInt_AsLong (item);
    } else {
      return (int) ERRSS ("Color tuple: expected integer value ");
    }
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

/*  -------------------------------------------------------------------- */

static char slice2__doc__[] =
"[nvf, xyzvf, colorf] = slice2 ((plane, nverts, xyzverts, values = None)\n"
"     Slice a polygon list, returning in nvf and xyzvf only those\n"
"     polygons or parts of polygons on the positive side of PLANE.\n"
"     If PLANE is a scalar real, then VALUES must be a function\n"
"     defined on the vertices of the mesh, and the mesh will\n"
"     be sliced where the function has that value.\n"
"     The NVERTS, XYZVERTS, and VALUES arrays have the meanings\n"
"     of the return values from the slice3 function. It is legal\n"
"     to omit the VALUES argument (e.g.- if there is no fcolor\n"
"     function).\n"
"     In order to plot two intersecting slices, one could\n"
"     slice (for example) the horizontal plane twice (slice2x) -\n"
"     first with the plane of the vertical slice, then with minus\n"
"     that same plane.  Then, plot first the back part of the\n"
"     slice, then the vertical slice, then the front part of the\n"
"     horizontal slice.  Of course, the vertical plane could\n"
"     be the one to be sliced, and `back' and `front' vary\n"
"     depending on the view point, but the general idea always\n"
"     works.\n"
"     \n"
"     slice2_precision= precision\n"
"     Controls how slice2 (or slice2x) handles points very close to\n"
"     the slicing plane.  PRECISION should be a positive number or zero.\n"
"     Zero PRECISION means to clip exactly to the plane, with points\n"
"     exactly on the plane acting as if they were slightly on the side\n"
"     the normal points toward.  Positive PRECISION means that edges\n"
"     are clipped to parallel planes a distance PRECISION on either\n"
"     side of the given plane.  (Polygons lying entirely between these\n"
"     planes are completely discarded.)\n"
"     \n"
"     Default value is 0.0.\n";

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
 ArrayObject * rnverts = (ArrayObject *) NULL,
             * rxyzverts = (ArrayObject *) NULL,
             * rvalues = (ArrayObject *) NULL,
             * rnvertb = (ArrayObject *) NULL,
             * rxyzvertb = (ArrayObject *) NULL,
             * rvalueb = (ArrayObject *) NULL,
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
             * mask2,
             * list = (ArrayObject *) NULL,
             * listc = (ArrayObject *) NULL;
 double * planed=0,
        * xyzvertsd=0,
        * valuesd = (double *) NULL,
        * rxyzvertsd=0,
        * rvaluecd = (double *) NULL,
        * rvaluesd = (double *) NULL,
        * rxyzvertbd=0,
        * rvaluebd = (double *) NULL,
        * dpd = (double *) NULL,
        * ndpd=0,
        * xyzcd=0,
        * xyzc0d=0,
        * valuecd = (double *) NULL,
        * valuec0d = (double *) NULL;
 int * nvertsd=0,
     * rnvertsd=0,
     * rnvertbd=0,
     * nvertcd=0,
     * nvertc0d=0,
     * nkeepd=0,
     * nkeep2d=0,
     * prevd=0,
     * nextd=0,
     * lastd=0,
     * listd=0,
     * listcd=0;
 Uchar * keepd=0,
       * mask2d=0,
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
 double dplane=0.;
 
 if (!PyArg_ParseTuple (args, "OOO|Oi", &oplane, &onverts, &oxyzverts,
    &ovalues, &_slice2x))  {
    return ERRSS ("slice2: unable to parse arguments.");
 }
 plane_is_scalar = PyFloat_Check (oplane);
 if (plane_is_scalar)
    dplane = PyFloat_AsDouble (oplane);
 else {
    GET_ARR (aplane, oplane, PyArray_DOUBLE, 1, PyObject *);
    }
 /* convert arguments to arrays */
 GET_ARR (anverts, onverts, PyArray_INT, 1, PyObject *);
 GET_ARR (axyzverts, oxyzverts, PyArray_DOUBLE, 2, PyObject *);
 if (isARRAY (ovalues)) {
    if (A_TYPE (ovalues) == PyArray_DOUBLE) {
       GET_ARR (avalues, ovalues, PyArray_DOUBLE, 1, PyObject *);
       valuesd = (double *) A_DATA (avalues);
       atype = 'd';
       }
    else if (A_TYPE (ovalues) == Py_GpColor) {
       GET_ARR (avalues, ovalues, Py_GpColor, 1, PyObject *);
       valuesc = (Uchar *) A_DATA (avalues);
       atype = 'b';
       }
    else {
       clearFreeList (0);
       return ERRSS ("Data type for values must be 'b' or 'd'.");
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
       clearFreeList (0);
       return ERRSS (
       "Number of data values must equal number of nodes or number of cells.");
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
    clearFreeList (0);
    return ERRSS ("Not sure what kind of slice you're asking for.");
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
   if (avalues)  {
      if (atype == 'd')
         valuecd = (double *) (valuec->data);
      else
         valuecc = (Uchar *) (valuec->data);
   }
   for (i = 0, k = 0, sumv = 0, sumt = 0; i < nkeep->size; i++) {
      if (nkeepd [i] != 0 && nkeepd [i] != nvertsd [i]) {
         nvertcd [k] = nvertsd [i];
         if (avalues && node == 0)  {
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
            if (avalues && node == 1)  {
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
                  if (avalues && node == 1)  {
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
               if (avalues && node != 0)  {
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
             if (avalues && node == 0)  {
                if (atype == 'd')
                   rvaluesd [k] = valuesd [i];
                else
                   rvaluesc [k] = valuesc [i];
             }
             for (j = 0; j < nvertsd [i]; j++) {
                rxyzvertsd [3 * (sumv + j)] = xyzvertsd [3 * (sumt + j)];
                rxyzvertsd [3 * (sumv + j) + 1] = xyzvertsd [3 * (sumt + j) + 1];
                rxyzvertsd [3 * (sumv + j) + 2] = xyzvertsd [3 * (sumt + j) + 2];
                if (avalues && node != 0)  {
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

/*  -------------------------------------------------------------------- */

static char unzoom__doc__[] =
"unzoom()\n"
"     Restore limits to their values before zoom and pan operations\n"
"     performed interactively using the mouse.\n"
"     Use    old_limits = limits()\n"
"            ...\n"
"            limits( old_limits )\n"
"     to save and restore plot limits generally.\n"
"\n"
"   SEE ALSO: limits, ylimits, zoom_factor, plg\n";

static PyObject *unzoom (PyObject * self, PyObject * args)
{
  GdRevertLimits (1);
  Py_INCREF (Py_None);
#ifdef WINDOWS
  pyg_on_idle();
#endif
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

/*  -------------------------------------------------------------------- */

static char viewport__doc__[] =
"Return viewport.";

static PyObject *viewport (PyObject * self, PyObject * args)
{
  double xmin, xmax, ymin, ymax;
  xmin = gistD.trans.viewport.xmin;
  xmax = gistD.trans.viewport.xmax;
  ymin = gistD.trans.viewport.ymin;
  ymax = gistD.trans.viewport.ymax;
  return Py_BuildValue ("dddd", xmin, xmax, ymin, ymax);
}

/*  -------------------------------------------------------------------- */
/*  window */

static char window__doc__[] =
"window( [n] [, display = `host:server.screen', dpi=100/75, wait=0/1,\n"
"                       private=0/1, hcp=`hcp_filename', dump=0/1,\n"
"                       legends=1/0, style=`style_sheet_filename' ] )\n"
"     select window N as the current graphics output window.  N may\n"
"     range from 0 to 7, inclusive.  Each graphics window corresponds to\n"
"     an X window, and optionally has its own associated hardcopy file.\n"
"     If N is omitted, it defaults to the current coordinate system.\n"
"\n"
"     The X window will appear on your default display at 75 dpi, unless\n"
"     you specify the display and/or dpi keywords.  A dpi=100 X window\n"
"     is larger than a dpi=75 X window; both represent the same thing\n"
"     on paper.  Use display=`' to create a graphics window which has\n"
"     no associated X window (you should do this if you want to make\n"
"     plots in a non-interactive batch mode).\n"
"\n"
"     By default, an X window will attempt to use shared colors, which\n"
"     permits several Pygist graphics windows (including windows from\n"
"     multiple instances of Python) to use a common palette.  You can\n"
"     force an X window to post its own colormap (set its colormap\n"
"     attribute) with the private=1 keyword.  You will most likely have\n"
"     to fiddle with your window manager to understand how it handles\n"
"     colormap focus if you do this.  Use private=0 to return to shared\n"
"     colors.\n"
"\n"
"     By default, Python will not wait for the X window to become visible;\n"
"     code which creates a new window, then plots a series of frames to\n"
"     that window should use wait=1 to assure that all frames are actually\n"
"     plotted.\n"
"\n"
"     By default, a graphics window does NOT have a hardcopy file\n"
"     of its own -- any request for hardcopy are directed to the\n"
"     default hardcopy file, so hardcopy output from any window goes\n"
"     to a single file.  By specifying the hcp keyword, however, a\n"
"     hardcopy file unique to this window will be created.  If the\n"
"     `hcp_filename' ends in `.ps', the hardcopy file will be a PostScript\n"
"     file; otherwise, hardcopy files are in binary CGM format.  Use\n"
"     hcp=`' to revert to the default hardcopy file (closing the window\n"
"     specific file, if any).  The legends keyword, if present, controls\n"
"     whether the curve legends are (legends=1, the default) or are not\n"
"     (legends=0) dumped to the hardcopy file.  The dump keyword, if\n"
"     present, controls whether all colors are converted to a gray scale\n"
"     (dump=0, the default), or the current palette is dumped at the\n"
"     beginning of each page of hardcopy output.  (The legends keyword\n"
"     applies to all pictures dumped to hardcopy from this graphics\n"
"     window.  The dump keyword applies only to the specific hardcopy\n"
"     file defined using the hcp keyword -- use the dump keyword in the\n"
"     hcp_file command to get the same effect in the default hardcopy\n"
"     file.)\n"
"\n"
"     If both display=`'; and hcp=`', the graphics window will be\n"
"     entirely eliminated.\n"
"\n"
"     The style keyword, if present, specifies the name of a Gist style\n"
"     sheet file; the default is `work.gs'.  The style sheet determines\n"
"     the number and location of coordinate systems, tick and label styles,\n"
"     and the like.  Other choices include `axes.gs', `boxed.gs',\n"
"     `work2.gs', and `boxed2.gs'.\n"
"\n"
"     Window(...) returns the current window number.\n"
"\n"
"   SEE ALSO: plsys, hcp_file, fma, hcp, redraw, palette, animate, plg,\n"
"             winkill, gridxy\n";

#undef N_KEYWORDS
#define N_KEYWORDS 11
static char *windowKeys[N_KEYWORDS+1]= {
  "display", "dpi", "private", "hcp", "legends", "dump", "style", "wait",
  "width", "height", "rgb", 0 };

static PyObject *window (PyObject * self, PyObject * args, PyObject * kd)
{
  int n, nGiven;
  Drauing *drawing;
  GpColorCell *palette;
  PyObject * kwt[NELT(windowKeys) - 1];
  int nColors = 0;
  int wait_for_expose = 0;
  int rgb = 0;

  SETJMP0;			/* See Xerror_longjmp() */

  if (!PyArg_ParseTuple (args, "|i", &n))  {
    return ERRSS ("window takes zero or one non-keyword integer argument."); 
  }

  if(PyTuple_Size(args) == 1) { /* Window() was called with an argument. */
    if (n < 0 || n > 7)  {
      return ERRSS ("graphics windows are numbered from 0 to 7");
    }
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

  if (nGiven || kwt[0] || kwt[1] || kwt[2]) {
    /* display= and/or dpi= keywords */
    char *display = 0;
    int dpi = defaultDPI;
    int privmap = 0;
    Engine *engine = ghDevices[n].display;	/* current display engine */

    SETKW(kwt[0], display, setkw_string, windowKeys[0]);
    if (kwt[1]) {
      if (engine)  {
	return ERRSS ("cannot change dpi of an existing graphics window");
      }
      SETKW(kwt[1], dpi, setkw_integer, windowKeys[1]);
      if (dpi<25) dpi = 25;
      else if (dpi>300) dpi = 300;
    }
    if (kwt[2]) {
      /* private= keyword -- turn on/off private X window colormap */
      if (engine)  {
        return ERRSS ("cannot give existing graphics window private colormap");
      }
      if (!(nGiven? (!display || display[0]) : (display && display[0])))  {
        return ERRSS ("private= keyword not legal without display engine");
      }
      SETKW(kwt[2], privmap, setkw_boolean, windowKeys[2]);
    }
    if (kwt[10]) {
      /* rgb= keyword -- maybe make this a true color window */
      if (engine)  {
        return ERRSS ("cannot use rgb= on existing graphics window");
      }
      if (!(nGiven? (!display || display[0]) : (display && display[0])))  {
        return ERRSS ("rgb= keyword not legal without display engine");
      }
      SETKW(kwt[2], rgb, setkw_boolean, windowKeys[2]);
    }

    if (engine) {
      ghDevices[n].display = 0;
      GpKillEngine (engine);
    }
    if (nGiven ? (!display || display[0]) : (display && display[0]))  {
#ifndef NO_XLIB
      gist_private_map = privmap;
      gist_rgb_hint = rgb;
      engine= DISPLAY_ENGINE(windowNames[n], 0, dpi, display);
      if (!engine)  {
	return ERRSS ("failed to open X display or create X window");
      } else {
        wait_for_expose = 1;
      }
      ghDevices[n].display = engine;
      if (palette)
	GhSetPalette (n, palette, nColors);
#else
      return ERRSS ("No interactive graphics in this Pygist -- hcp only");
#endif
    }
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
	if (!engine)  {
	  return ERRSS ("failed to create PostScript file");
        }
      } 
      else {
	engine = GpCGMEngine (windowNames[n], 0, hcpDump, SetHCPname (n, hcp));
	if (!engine)  {
	  return ERRSS ("failed to create binary CGM file");
        }
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
    if (!ghDevices[n].hcp)  {
      return ERRSS (
	"dump= keyword not legal without hcp engine -- use hcp_file");
    }
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
  } 
  else {
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
      if (drawing)  {
	return ERRSS ("failed to create drawing -- bad style sheet name?");
      }
      else  {
	return ERRSS (
	  "failed to create drawing -- Gist work.gs style sheet missing");
      }
    }
    /* make this window current */
    curPlotter = n;
    GhSetPlotter (n);
    paletteSize = nColors;

    /* wait= keyword -- pause until X window is exposed */
#ifdef WINDOWS
    /* The wait causes PythonWin to hang. Since it doesn't seem to be needed
     * anyway, just leave it out. Note that wait=1 is used in gistdemolow.
     */
    wait_for_expose = 0;
#else
    if (kwt[7] && wait_for_expose) {
      int wait;
      SETKW(kwt[7], wait, setkw_boolean, windowKeys[7]);


      wait_for_expose = wait_for_expose && (1 == wait);
    } else {
      wait_for_expose = 0;
    }
#endif
  }

  /* under MS Windows, the first expose event occurs synchronously inside
   * the GpFXEngine call -- hence the oops==2 logic */
  if (wait_for_expose) {
    int oops = gist_expose_wait(ghDevices[n].display, pyg_got_expose);
    if (oops==1 || pyg_wait_flag) {
      /* hopefully this is impossible */
      return ERRSS ("window,wait=1 while already waiting for a window");
    }
    /*
     *  Pause for window to pop up  
     */
    if (oops != 2) {
      pyg_wait_flag = 1;
      p_wait_while(&pyg_wait_flag);
    }
  }

  return Py_BuildValue ("i",n);
}

static void
pyg_got_expose(void)
{
  pyg_wait_flag = 0;
}

/*  -------------------------------------------------------------------- */

static char zoom_factor__doc__[] =
"zoom_factor( factor )\n"
"     Set the zoom factor for mouse-click zoom in and zoom out operations.\n"
"     The default FACTOR is 1.5; FACTOR should always be greater than 1.0.\n"
"\n"
"   SEE ALSO: limits, ylimits, unzoom, plg\n";

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

static char pyg_unhook__doc__[] =
"pyg_unhoook( )\n"
"     Remove pygist PyOS_InputHook if present (for _tkinter).\n";

static PyObject *pyg_unhook (PyObject * self, PyObject * args)
{
  /* _tkinter needs PyOs_InputHook */
#ifdef CYGWIN
  if (PyOS_InputHook == p_pending_events) PyOS_InputHook = 0;
#else
  if (PyOS_InputHook == p_wait_stdin) PyOS_InputHook = 0;
#endif
  Py_INCREF (Py_None);
  return Py_None;
}

/* pyg_idler
 * needs to be called by tcl's event loop
 * via tcl command "after idle pyg_idler" in order to make pygist
 * function together with _tkinter
 * -- this technique replaces PyOS_InputHook, which _tkinter needs
 *    (1) if PyOS_InputHook is already set when gistC loads,
 *        you *must* arrange to have pyg_idler called at idle time
 *    (2) if gistC loads before _tkinter, you *must* call pyg_unhook
 *        first in order for _tkinter to function properly
 */
static char pyg_idler__doc__[] =
"pyg_idler( )\n"
"     Do any deferred pygist window operations.\n";

static PyObject *pyg_idler (PyObject * self, PyObject * args)
{
  p_on_idle(0);
  Py_INCREF (Py_None);
  return Py_None;
}

static char pyg_pending__doc__[] =
"pyg_pending( )\n"
"     Handle any pending pygist window events.\n";

static PyObject *pyg_pending (PyObject * self, PyObject * args)
{
  p_pending_events();
  Py_INCREF (Py_None);
  return Py_None;
}

static PyObject *pyg_connector = 0;
static PyObject *pyg_keyhandler = 0;

static char pyg_register__doc__[] =
"pyg_register( connector, keyhandler )\n"
"     connector(dis,fd) will be called on gist connect/disconnect.\n"
"     keyhandler(line) will be called when line typed in gist window.\n";

/* wire up pyg_on_connect function */
static PyObject *pyg_register (PyObject * self, PyObject * args)
{
  PyObject *func1, *func2;
  if (!PyArg_ParseTuple(args, "|OO", &func1, &func2) ||
      (func1 && !PyCallable_Check(func1)) ||
      (func2 && !PyCallable_Check(func2)) ) {
    return ERRSS ("pyg_register takes two function arguments.");
  }
  if (pyg_connector) { Py_DECREF (pyg_connector); }
  pyg_connector = func1;
  Py_INCREF (func1);
  if (pyg_keyhandler) { Py_DECREF (pyg_keyhandler); }
  pyg_keyhandler = func2;
  if (func2) Py_INCREF (func2);
  Py_INCREF (Py_None);
  return Py_None;
}

/*
 *  10/30/01 llc Moved PyMethodDef to end, after doc strings are defined.
 *               Also move initgistC, which uses gist_methods.
 */

static struct PyMethodDef gist_methods[] =
{ 
  { "animate",        PYCF   animate,        1,     animate__doc__ },
  { "bytscl",         PYCFWK bytscl,         KWFLG, bytscl__doc__ },
  { "contour",        PYCFWK contour,        KWFLG, contour__doc__ },
  { "current_window", PYCF   current_window, 1,     current_window__doc__ },
  { "debug_array",    PYCF   debug_array,    1,     debug_array__doc__ },
  { "fma",            PYCF   pyg_fma,        1,     fma__doc__ },
  { "gridxy",         PYCFWK gridxy,         KWFLG, gridxy__doc__ },
  { "get_slice2_precision", PYCF get_slice2_precision, 1, get_slice2_precision__doc__ },
  { "hcp",            PYCF   hcp,            1,     hcp__doc__ },
  { "hcp_file",       PYCFWK hcp_file,       KWFLG, hcp_file__doc__ },
  { "hcp_finish",     PYCF   hcp_finish,     1,     hcp_finish__doc__ },
  { "hcpoff",         PYCF   hcpoff,         1,     hcpoff__doc__ },
  { "hcpon",          PYCF   hcpon,          1,     hcpon__doc__ },
  { "limits",         PYCFWK limits,         KWFLG, limits__doc__ },
  { "logxy",          PYCF   logxy,          1,     logxy__doc__ },
  { "mesh_loc",       PYCF   mesh_loc,       1,     mesh_loc__doc__ },
  { "mfit",           PYCF   mfit,           1,     mfit__doc__ },
  { "mouse",          PYCF   mouse,          1,     mouse__doc__ },
  { "palette",        PYCFWK palette,        KWFLG, palette__doc__ },
  { "pause",          PYCF   pyg_pause,      1,     pause__doc__ },
  { "plc",            PYCFWK plc,            KWFLG, plc__doc__ },
  { "pldefault",      PYCFWK pldefault,      KWFLG, pldefault__doc__ },
  { "pldj",           PYCFWK pldj,           KWFLG, pldj__doc__ },
  { "pledit",         PYCFWK pledit,         KWFLG, pledit__doc__ },
  { "plf",            PYCFWK plf,            KWFLG, plf__doc__ },
  { "plfp",           PYCFWK plfp,           KWFLG, plfp__doc__ },
  { "plg",            PYCFWK plg,            KWFLG, plg__doc__ },
  { "pli",            PYCFWK pli,            KWFLG, pli__doc__ },
  { "plm",            PYCFWK plm,            KWFLG, plm__doc__ },
  { "plmesh",         PYCFWK plmesh,         KWFLG, plmesh__doc__ },
  { "plq",            PYCF   plq,            1,     plq__doc__ },
  { "plremove",       PYCF   plremove,       1,     plremove__doc__ },
  { "plsys",          PYCF   plsys,          1,     plsys__doc__ },
  { "plt",            PYCFWK plt,            KWFLG, plt__doc__ },
  { "plv",            PYCFWK plv,            KWFLG, plv__doc__ },
  { "redraw",         PYCF   redraw,         1,     redraw__doc__ },
  { "set_slice2_precision", PYCF set_slice2_precision, 1, set_slice2_position__doc__ },
  { "slice2",         PYCF   slice2,         1,     slice2__doc__ },
  { "unzoom",         PYCF   unzoom,         1,     unzoom__doc__ },
  { "viewport",       PYCF   viewport,       1,     viewport__doc__ },
  { "window",         PYCFWK window,         KWFLG, window__doc__ },
  { "zoom_factor",    PYCF   zoom_factor,    1,     zoom_factor__doc__ },
  { "pyg_unhook",     PYCF   pyg_unhook,     1,     pyg_unhook__doc__ },
  { "pyg_idler",      PYCF   pyg_idler,      1,     pyg_idler__doc__ },
  { "pyg_pending",    PYCF   pyg_pending,    1,     pyg_pending__doc__ },
  { "pyg_register",   PYCF   pyg_register,   1,     pyg_register__doc__ },

  { 0, 0 }
};

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

#ifdef DEBUG
  TO_STDOUT ( "\nPyGist version 1.5.11 ($Id$)\n"
"    This version of gist uses pydoc for documentation.\n"
"       help(function_name)\n"
"    provides documentation for function_name\n"
"    Hit spacebar to page down, and 'q' to end documentation\n"
"    Single quotes and backquotes delimiting strings in documentation\n"
"    should be double quotes.\n\n" );
#endif

#ifdef import_array
  import_array();
#endif

  {
    /* DHM: in principal, gist might use argv[0] to try to figure out
     *      the PATH to this executable, from which it (in principle)
     *      might be able to find the g/ subdirectory for GISTPATH
     *      - this is done better below
     *      the other reason is that the argv[] might contain standard
     *      X resource switches, which, however, gist does not use */
    int argc = 0;
    char **argv = 0;
    g_initializer(&argc, argv);
  }

  if (0 != Py_AtExit (CleanUpGraphics)) {
    TO_STDERR("Gist: Warning: Exit procedure not registered\n");
    flush_stderr();
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

  {
    /* set up play p_abort to work with SETJMP0 macro used here */
    /* do not call p_handler -- would disturb python signal handling */
    p_xhandler(pyg_abort_hook, pyg_on_exception);

    /* note that g_on_keyline might be useful, especially for Windows */
    g_on_keyline = pyg_on_keyline;

    /*  Provide a way for gist to process its events
     * - the PyOS_InputHook is actually for exclusive use of the
     *   _tkinter module, so if it is already set, leave it alone
     */
    if (!PyOS_InputHook)
#ifdef CYGWIN
      PyOS_InputHook = p_pending_events;
#else
      PyOS_InputHook = p_wait_stdin;
#endif
    p_on_connect = pyg_on_connect;

    /* turn on idle function to do graphics tasks */
    /* Call p_idler to set up the idle callback   */
    p_idler(pyg_on_idle);
  }

  already_initialized = 1;

  if ( setjmp ( pyg_jmpbuf ) )  {
     p_pending_events();
     return;
  }
}

static int
pyg_on_idle(void)
{
  /* 
   *  Gist does all its drawing in GhBeforeWait, which should be called
   *  at idle time.
   */

  GhBeforeWait();
#ifdef CYGWIN
  p_pending_events();
#endif
  return 0;
}

/* p_on_connect(dis,fd)
 * called by gist/play when p_connect or p_disconnect makes or breaks
 * a display connection
 *   dis = 0 for connect, 1 for disconnect
 *   fd = file descriptor of X socket under UNIX/X11
 *        -1 under MS Windows
 * useful for wiring up pyg_pending to events arriving on fd
 *   under UNIX/X11, this requires createfilehandler _tkinter method
 *   under Windows, pyg_pending is unnecessary because each gist
 *     window has its own class and gets its messages delivered to
 *     the play/win/pscr.c w_winproc function automatically
 *
 * note that pyg_connector is set by calling pyg_register
 */
static void
pyg_on_connect(int dis, int fd)
{
  PyObject *args, *res;
  if (fd<0 || !pyg_connector) return;
  args = Py_BuildValue("(ii)", dis, fd);
  res = PyEval_CallObject(pyg_connector, args);
  Py_DECREF(args);
  /* nothing to be done if res==NULL?? */
  Py_XDECREF(res);
}

/* g_on_keyline(msg) is called when msg is typed in the gist
 * GpFXEngine window; the callback occurs when RET is pressed */
static void
pyg_on_keyline(char *msg)
{
  PyObject *args, *res;
  if (!pyg_keyhandler) return;
  args = Py_BuildValue("(s)", msg);
  res = PyEval_CallObject(pyg_keyhandler, args);
  Py_DECREF(args);
  /* nothing to be done if res==NULL?? */
  Py_XDECREF(res);
}

#ifdef __cplusplus
}
#endif
