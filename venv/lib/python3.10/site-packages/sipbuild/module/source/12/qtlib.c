/*
 * The SIP library code that implements the interface to the optional module
 * supplied Qt support.
 *
 * Copyright (c) 2022 Riverbank Computing Limited <info@riverbankcomputing.com>
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


#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <assert.h>
#include <string.h>

#include "sipint.h"


/* This is how Qt "types" signals and slots. */
#define isQtSlot(s)     (*(s) == '1')
#define isQtSignal(s)   (*(s) == '2')


static PyObject *getWeakRef(PyObject *obj);
static char *sipStrdup(const char *);
static void *createUniversalSlot(sipWrapper *txSelf, const char *sig,
        PyObject *rxObj, const char *slot, const char **member, int flags);
static void *findSignal(void *txrx, const char **sig);
static void *newSignal(void *txrx, const char **sig);


/*
 * Find an existing signal.
 */
static void *findSignal(void *txrx, const char **sig)
{
    if (sipQtSupport->qt_find_universal_signal != NULL)
        txrx = sipQtSupport->qt_find_universal_signal(txrx, sig);

    return txrx;
}


/*
 * Return a usable signal, creating a new universal signal if needed.
 */
static void *newSignal(void *txrx, const char **sig)
{
    void *new_txrx = findSignal(txrx, sig);

    if (new_txrx == NULL && sipQtSupport->qt_create_universal_signal != NULL)
        new_txrx = sipQtSupport->qt_create_universal_signal(txrx, sig);

    return new_txrx;
}


/*
 * Create a universal slot.  Returns a pointer to it or 0 if there was an
 * error.
 */
static void *createUniversalSlot(sipWrapper *txSelf, const char *sig,
        PyObject *rxObj, const char *slot, const char **member, int flags)
{
    void *us;

    assert(sipQtSupport->qt_create_universal_slot);

    us = sipQtSupport->qt_create_universal_slot(txSelf, sig, rxObj, slot,
            member, flags);

    if (us && txSelf)
        sipSetPossibleProxy((sipSimpleWrapper *)txSelf);

    return us;
}


/*
 * Invoke a single slot (Qt or Python) and return the result.  Don't check if
 * any receiver C++ object still exists.
 */
PyObject *sip_api_invoke_slot(const sipSlot *slot, PyObject *sigargs)
{
    return sip_api_invoke_slot_ex(slot, sigargs, TRUE);
}


/*
 * Invoke a single slot (Qt or Python) and return the result.  Optionally check
 * that any receiver C++ object still exist.
 */
PyObject *sip_api_invoke_slot_ex(const sipSlot *slot, PyObject *sigargs,
        int no_receiver_check)
{
    PyObject *sa, *oxtype, *oxvalue, *oxtb, *sfunc, *sref;

    assert(sipQtSupport);

    /* Keep some compilers quiet. */
    oxtype = oxvalue = oxtb = NULL;

    /* Fan out Qt signals.  (Only PyQt3 would do this.) */
    if (slot->name != NULL && slot->name[0] != '\0')
    {
        assert(sipQtSupport->qt_emit_signal);

        if (sipQtSupport->qt_emit_signal(slot->pyobj, slot->name, sigargs) < 0)
            return NULL;

        Py_INCREF(Py_None);
        return Py_None;
    }

    /* Get the object to call, resolving any weak references. */
    if (slot->weakSlot == Py_True)
    {
        /*
         * The slot is guaranteed to be Ok because it has an extra reference or
         * is None.
         */
        sref = slot->pyobj;
        Py_INCREF(sref);
    }
    else if (slot -> weakSlot == NULL)
        sref = NULL;
    else if ((sref = PyWeakref_GetObject(slot -> weakSlot)) == NULL)
        return NULL;
    else
        Py_INCREF(sref);

    if (sref == Py_None)
    {
        /*
         * If the real object has gone then we pretend everything is Ok.  This
         * mimics the Qt behaviour of not caring if a receiving object has been
         * deleted.
         */
        Py_DECREF(sref);

        Py_INCREF(Py_None);
        return Py_None;
    }

    if (slot -> pyobj == NULL)
    {
        PyObject *self = (sref != NULL ? sref : slot->meth.mself);

        /*
         * If the receiver wraps a C++ object then ignore the call if it no
         * longer exists.
         */
        if (!no_receiver_check &&
            PyObject_TypeCheck(self, (PyTypeObject *)&sipSimpleWrapper_Type) &&
            sip_api_get_address((sipSimpleWrapper *)self) == NULL)
        {
            Py_XDECREF(sref);

            Py_INCREF(Py_None);
            return Py_None;
        }

        if ((sfunc = PyMethod_New(slot->meth.mfunc, self)) == NULL)
        {
            Py_XDECREF(sref);
            return NULL;
        }
    }
    else if (slot -> name != NULL)
    {
        char *mname = slot -> name + 1;
        PyObject *self = (sref != NULL ? sref : slot->pyobj);

        if ((sfunc = PyObject_GetAttrString(self, mname)) == NULL || !PyCFunction_Check(sfunc))
        {
            /*
             * Note that in earlier versions of SIP this error would be
             * detected when the slot was connected.
             */
            PyErr_Format(PyExc_NameError,"Invalid slot %s",mname);

            Py_XDECREF(sfunc);
            Py_XDECREF(sref);
            return NULL;
        }
    }
    else
    {
        sfunc = slot->pyobj;
        Py_INCREF(sfunc);
    }

    /*
     * We make repeated attempts to call a slot.  If we work out that it failed
     * because of an immediate type error we try again with one less argument.
     * We keep going until we run out of arguments to drop.  This emulates the
     * Qt ability of the slot to accept fewer arguments than a signal provides.
     */
    sa = sigargs;
    Py_INCREF(sa);

    for (;;)
    {
        PyObject *nsa, *xtype, *xvalue, *xtb, *resobj;

        if ((resobj = PyObject_CallObject(sfunc, sa)) != NULL)
        {
            Py_DECREF(sfunc);
            Py_XDECREF(sref);

            /* Remove any previous exception. */

            if (sa != sigargs)
            {
                Py_XDECREF(oxtype);
                Py_XDECREF(oxvalue);
                Py_XDECREF(oxtb);
                PyErr_Clear();
            }

            Py_DECREF(sa);

            return resobj;
        }

        /* Get the exception. */
        PyErr_Fetch(&xtype,&xvalue,&xtb);

        /*
         * See if it is unacceptable.  An acceptable failure is a type error
         * with no traceback - so long as we can still reduce the number of
         * arguments and try again.
         */
        if (!PyErr_GivenExceptionMatches(xtype,PyExc_TypeError) ||
            xtb != NULL ||
            PyTuple_GET_SIZE(sa) == 0)
        {
            /*
             * If there is a traceback then we must have called the slot and
             * the exception was later on - so report the exception as is.
             */
            if (xtb != NULL)
            {
                if (sa != sigargs)
                {
                    Py_XDECREF(oxtype);
                    Py_XDECREF(oxvalue);
                    Py_XDECREF(oxtb);
                }

                PyErr_Restore(xtype,xvalue,xtb);
            }
            else if (sa == sigargs)
                PyErr_Restore(xtype,xvalue,xtb);
            else
            {
                /*
                 * Discard the latest exception and restore the original one.
                 */
                Py_XDECREF(xtype);
                Py_XDECREF(xvalue);
                Py_XDECREF(xtb);

                PyErr_Restore(oxtype,oxvalue,oxtb);
            }

            break;
        }

        /* If this is the first attempt, save the exception. */
        if (sa == sigargs)
        {
            oxtype = xtype;
            oxvalue = xvalue;
            oxtb = xtb;
        }
        else
        {
            Py_XDECREF(xtype);
            Py_XDECREF(xvalue);
            Py_XDECREF(xtb);
        }

        /* Create the new argument tuple. */
        if ((nsa = PyTuple_GetSlice(sa,0,PyTuple_GET_SIZE(sa) - 1)) == NULL)
        {
            /* Tidy up. */
            Py_XDECREF(oxtype);
            Py_XDECREF(oxvalue);
            Py_XDECREF(oxtb);

            break;
        }

        Py_DECREF(sa);
        sa = nsa;
    }

    Py_DECREF(sfunc);
    Py_XDECREF(sref);

    Py_DECREF(sa);

    return NULL;
}


/*
 * Compare two slots to see if they are the same.
 */
int sip_api_same_slot(const sipSlot *sp, PyObject *rxObj, const char *slot)
{
    assert(sipQtSupport);
    assert(sipQtSupport->qt_same_name);

    /* See if they are signals or Qt slots, ie. they have a name. */
    if (slot != NULL)
    {
        if (sp->name == NULL || sp->name[0] == '\0')
            return 0;

        return (sipQtSupport->qt_same_name(sp->name, slot) && sp->pyobj == rxObj);
    }

    /* See if they are pure Python methods. */
    if (PyMethod_Check(rxObj))
    {
        if (sp->pyobj != NULL)
            return 0;

        return (sp->meth.mfunc == PyMethod_GET_FUNCTION(rxObj)
                && sp->meth.mself == PyMethod_GET_SELF(rxObj));
    }

    /* See if they are wrapped C++ methods. */
    if (PyCFunction_Check(rxObj))
    {
        if (sp->name == NULL || sp->name[0] != '\0')
            return 0;

        return (sp->pyobj == PyCFunction_GET_SELF(rxObj) &&
                strcmp(&sp->name[1], ((PyCFunctionObject *)rxObj)->m_ml->ml_name) == 0);
    }

    /* The objects must be the same. */
    return (sp->pyobj == rxObj);
}


/*
 * Convert a valid Python signal or slot to an existing universal slot.
 */
void *sipGetRx(sipSimpleWrapper *txSelf, const char *sigargs, PyObject *rxObj,
           const char *slot, const char **memberp)
{
    assert(sipQtSupport);
    assert(sipQtSupport->qt_find_slot);

    if (slot != NULL)
        if (isQtSlot(slot) || isQtSignal(slot))
        {
            void *rx;

            *memberp = slot;

            if ((rx = sip_api_get_cpp_ptr((sipSimpleWrapper *)rxObj, sipQObjectType)) == NULL)
                return NULL;

            if (isQtSignal(slot))
                rx = findSignal(rx, memberp);

            return rx;
        }

    /*
     * The slot was either a Python callable or PyQt3 Python signal so there
     * should be a universal slot.
     */
    return sipQtSupport->qt_find_slot(sip_api_get_address(txSelf), sigargs, rxObj, slot, memberp);
}


/*
 * Convert a Python receiver (either a Python signal or slot or a Qt signal or
 * slot) to a Qt receiver.  It is only ever called when the signal is a Qt
 * signal.  Return NULL is there was an error.
 */
void *sip_api_convert_rx(sipWrapper *txSelf, const char *sigargs,
        PyObject *rxObj, const char *slot, const char **memberp, int flags)
{
    assert(sipQtSupport);

    if (slot == NULL)
        return createUniversalSlot(txSelf, sigargs, rxObj, NULL, memberp, flags);

    if (isQtSlot(slot) || isQtSignal(slot))
    {
        void *rx;

        *memberp = slot;

        if ((rx = sip_api_get_cpp_ptr((sipSimpleWrapper *)rxObj, sipQObjectType)) == NULL)
            return NULL;

        if (isQtSignal(slot))
            rx = newSignal(rx, memberp);

        return rx;
    }

    /* The slot is a Python signal so we need a universal slot to catch it. */
    return createUniversalSlot(txSelf, sigargs, rxObj, slot, memberp, 0);
}


/*
 * Connect a Qt signal or a Python signal to a Qt slot, a Qt signal, a Python
 * slot or a Python signal.  This is all possible combinations.
 */
PyObject *sip_api_connect_rx(PyObject *txObj, const char *sig, PyObject *rxObj,
        const char *slot, int type)
{
    assert(sipQtSupport);
    assert(sipQtSupport->qt_connect);

    /* Handle Qt signals. */
    if (isQtSignal(sig))
    {
        void *tx, *rx;
        const char *member, *real_sig;
        int res;

        if ((tx = sip_api_get_cpp_ptr((sipSimpleWrapper *)txObj, sipQObjectType)) == NULL)
            return NULL;

        real_sig = sig;

        if ((tx = newSignal(tx, &real_sig)) == NULL)
            return NULL;

        if ((rx = sip_api_convert_rx((sipWrapper *)txObj, sig, rxObj, slot, &member, 0)) == NULL)
            return NULL;

        res = sipQtSupport->qt_connect(tx, real_sig, rx, member, type);

        return PyBool_FromLong(res);
    }

    /* Handle Python signals.  Only PyQt3 would get this far. */
    assert(sipQtSupport->qt_connect_py_signal);

    if (sipQtSupport->qt_connect_py_signal(txObj, sig, rxObj, slot) < 0)
        return NULL;

    Py_INCREF(Py_True);
    return Py_True;
}


/*
 * Disconnect a signal to a signal or a Qt slot.
 */
PyObject *sip_api_disconnect_rx(PyObject *txObj,const char *sig,
                PyObject *rxObj,const char *slot)
{
    assert(sipQtSupport);
    assert(sipQtSupport->qt_disconnect);
    assert(sipQtSupport->qt_destroy_universal_slot);

    /* Handle Qt signals. */
    if (isQtSignal(sig))
    {
        sipSimpleWrapper *txSelf = (sipSimpleWrapper *)txObj;
        void *tx, *rx;
        const char *member;
        int res;

        if ((tx = sip_api_get_cpp_ptr(txSelf, sipQObjectType)) == NULL)
            return NULL;

        if ((rx = sipGetRx(txSelf, sig, rxObj, slot, &member)) == NULL)
        {
            Py_INCREF(Py_False);
            return Py_False;
        }

        /* Handle Python signals. */
        tx = findSignal(tx, &sig);

        res = sipQtSupport->qt_disconnect(tx, sig, rx, member);

        /*
         * Delete it if it is a universal slot as this will be it's only
         * connection.  If the slot is actually a universal signal then it
         * should leave it in place.
         */
        sipQtSupport->qt_destroy_universal_slot(rx);

        return PyBool_FromLong(res);
    }

    /* Handle Python signals.  Only PyQt3 would get this far. */
    assert(sipQtSupport->qt_disconnect_py_signal);

    sipQtSupport->qt_disconnect_py_signal(txObj, sig, rxObj, slot);

    Py_INCREF(Py_True);
    return Py_True;
}


/*
 * Free the resources of a slot.
 */
void sip_api_free_sipslot(sipSlot *slot)
{
    assert(sipQtSupport);

    if (slot->name != NULL)
    {
        sip_api_free(slot->name);
    }
    else if (slot->weakSlot == Py_True)
    {
        Py_DECREF(slot->pyobj);
    }

    /* Remove any weak reference. */
    Py_XDECREF(slot->weakSlot);
}


/*
 * Implement strdup() using sip_api_malloc().
 */
static char *sipStrdup(const char *s)
{
    char *d;

    if ((d = (char *)sip_api_malloc(strlen(s) + 1)) != NULL)
        strcpy(d,s);

    return d;
}


/*
 * Initialise a slot, returning 0 if there was no error.  If the signal was a
 * Qt signal, then the slot may be a Python signal or a Python slot.  If the
 * signal was a Python signal, then the slot may be anything.
 */
int sip_api_save_slot(sipSlot *sp, PyObject *rxObj, const char *slot)
{
    assert(sipQtSupport);

    sp -> weakSlot = NULL;

    if (slot == NULL)
    {
        sp -> name = NULL;

        if (PyMethod_Check(rxObj))
        {
            /*
             * Python creates methods on the fly.  We could increment the
             * reference count to keep it alive, but that would keep "self"
             * alive as well and would probably be a circular reference.
             * Instead we remember the component parts and hope they are still
             * valid when we re-create the method when we need it.
             */
            sipSaveMethod(&sp -> meth,rxObj);

            /* Notice if the class instance disappears. */
            sp -> weakSlot = getWeakRef(sp -> meth.mself);

            /* This acts a flag to say that the slot is a method. */
            sp -> pyobj = NULL;
        }
        else
        {
            PyObject *self;

            /*
             * We know that it is another type of callable, ie. a
             * function/builtin.
             */

            if (PyCFunction_Check(rxObj) &&
                (self = PyCFunction_GET_SELF(rxObj)) != NULL &&
                PyObject_TypeCheck(self, (PyTypeObject *)&sipSimpleWrapper_Type))
            {
                /*
                 * It is a wrapped C++ class method.  We can't keep a copy
                 * because they are generated on the fly and we can't take a
                 * reference as that may keep the instance (ie. self) alive.
                 * We therefore treat it as if the user had specified the slot
                 * at "obj, SLOT('meth()')" rather than "obj.meth" (see below).
                 */

                const char *meth;

                /* Get the method name. */
                meth = ((PyCFunctionObject *)rxObj) -> m_ml -> ml_name;

                if ((sp -> name = (char *)sip_api_malloc(strlen(meth) + 2)) == NULL)
                    return -1;

                /*
                 * Copy the name and set the marker that it needs converting to
                 * a built-in method.
                 */
                sp -> name[0] = '\0';
                strcpy(&sp -> name[1],meth);

                sp -> pyobj = self;
                sp -> weakSlot = getWeakRef(self);
            }
            else
            {
                /*
                 * Give the slot an extra reference to keep it alive and
                 * remember we have done so by treating weakSlot specially.
                 */
                Py_INCREF(rxObj);
                sp->pyobj = rxObj;

                Py_INCREF(Py_True);
                sp->weakSlot = Py_True;
            }
        }
    }
    else if ((sp -> name = sipStrdup(slot)) == NULL)
        return -1;
    else if (isQtSlot(slot))
    {
        /*
         * The user has decided to connect a Python signal to a Qt slot and
         * specified the slot as "obj, SLOT('meth()')" rather than "obj.meth".
         */

        char *tail;

        /* Remove any arguments. */
        if ((tail = strchr(sp -> name,'(')) != NULL)
            *tail = '\0';

        /*
         * A bit of a hack to indicate that this needs converting to a built-in
         * method.
         */
        sp -> name[0] = '\0';

        /* Notice if the class instance disappears. */
        sp -> weakSlot = getWeakRef(rxObj);

        sp -> pyobj = rxObj;
    }
    else
        /* It's a Qt signal. */
        sp -> pyobj = rxObj;

    return 0;
}


/*
 * Return a weak reference to the given object.
 */
static PyObject *getWeakRef(PyObject *obj)
{
    PyObject *wr;

    if ((wr = PyWeakref_NewRef(obj,NULL)) == NULL)
        PyErr_Clear();

    return wr;
}
