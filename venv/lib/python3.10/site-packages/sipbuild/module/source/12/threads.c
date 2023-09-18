/*
 * Thread support for the SIP library.  This module provides the hooks for
 * C++ classes that provide a thread interface to interact properly with the
 * Python threading infrastructure.
 *
 * Copyright (c) 2020 Riverbank Computing Limited <info@riverbankcomputing.com>
 *
 * This file is part of SIP.
 *
 * This copy of SIP is licensed for use under the terms of the SIP License
 * Agreement.  See the file LICENSE for more details.
 *
 * This copy of SIP may also used under the terms of the GNU General Public
 * License v2 or v3 as published by the Free Software Foundation which can be
 * found in the files LICENSE-GPL2 and LICENSE-GPL3 included in this package.
 *
 * SIP is supplied WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */


#include "sipint.h"


/*
 * The data associated with pending request to wrap an object.
 */
typedef struct _pendingDef {
    void *cpp;                      /* The C/C++ object ot be wrapped. */
    sipWrapper *owner;              /* The owner of the object. */
    int flags;                      /* The flags. */
} pendingDef;


#ifdef WITH_THREAD

#include <pythread.h>


/*
 * The per thread data we need to maintain.
 */
typedef struct _threadDef {
    long thr_ident;                 /* The thread identifier. */
    pendingDef pending;             /* An object waiting to be wrapped. */
    struct _threadDef *next;        /* Next in the list. */
} threadDef;

static threadDef *threads = NULL;   /* Linked list of threads. */

static threadDef *currentThreadDef(int auto_alloc);

#endif


static pendingDef *get_pending(int auto_alloc);


/*
 * Get the address etc. of any C/C++ object waiting to be wrapped.
 */
int sipGetPending(void **pp, sipWrapper **op, int *fp)
{
    pendingDef *pd;

    if ((pd = get_pending(TRUE)) == NULL)
        return -1;

    *pp = pd->cpp;
    *op = pd->owner;
    *fp = pd->flags;

    /* Clear in case we execute Python code before finishing this wrapping. */
    pd->cpp = NULL;

    return 0;
}


/*
 * Return TRUE if anything is pending.
 */
int sipIsPending(void)
{
    pendingDef *pd;

    if ((pd = get_pending(FALSE)) == NULL)
        return FALSE;

    return (pd->cpp != NULL);
}


/*
 * Convert a new C/C++ pointer to a Python instance.
 */
PyObject *sipWrapInstance(void *cpp, PyTypeObject *py_type, PyObject *args,
        sipWrapper *owner, int flags)
{
    pendingDef old_pending, *pd;
    PyObject *self;

    if (cpp == NULL)
    {
        Py_INCREF(Py_None);
        return Py_None;
    }

    /*
     * Object creation can trigger the Python garbage collector which in turn
     * can execute arbitrary Python code which can then call this function
     * recursively.  Therefore we save any existing pending object before
     * setting the new one.
     */
    if ((pd = get_pending(TRUE)) == NULL)
        return NULL;

    old_pending = *pd;

    pd->cpp = cpp;
    pd->owner = owner;
    pd->flags = flags;

    self = PyObject_Call((PyObject *)py_type, args, NULL);

    *pd = old_pending;

    return self;
}


/*
 * Handle the termination of a thread.
 */
void sip_api_end_thread(void)
{
#ifdef WITH_THREAD
    threadDef *thread;
    PyGILState_STATE gil = PyGILState_Ensure();

    if ((thread = currentThreadDef(FALSE)) != NULL)
        thread->thr_ident = 0;

    PyGILState_Release(gil);
#endif
}


/*
 * Return the pending data for the current thread, allocating it if necessary,
 * or NULL if there was an error.
 */
static pendingDef *get_pending(int auto_alloc)
{
#ifdef WITH_THREAD
    threadDef *thread;

    if ((thread = currentThreadDef(auto_alloc)) == NULL)
        return NULL;

    return &thread->pending;
#else
    static pendingDef pending;

    return &pending;
#endif
}


#ifdef WITH_THREAD

/*
 * Return the thread data for the current thread, allocating it if necessary,
 * or NULL if there was an error.
 */
static threadDef *currentThreadDef(int auto_alloc)
{
    threadDef *thread, *empty = NULL;
    long ident = PyThread_get_thread_ident();

    /* See if we already know about the thread. */
    for (thread = threads; thread != NULL; thread = thread->next)
    {
        if (thread->thr_ident == ident)
            return thread;

        if (thread->thr_ident == 0)
            empty = thread;
    }

    if (!auto_alloc)
    {
        /* This is not an error. */
        return NULL;
    }

    if (empty != NULL)
    {
        /* Use an empty entry in the list. */
        thread = empty;
    }
    else if ((thread = sip_api_malloc(sizeof (threadDef))) == NULL)
    {
        return NULL;
    }
    else
    {
        thread->next = threads;
        threads = thread;
    }

    thread->thr_ident = ident;
    thread->pending.cpp = NULL;

    return thread;
}

#endif
