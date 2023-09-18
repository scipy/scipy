/*
 * This module implements a hash table class for mapping C/C++ addresses to the
 * corresponding wrapped Python object.
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


#include <string.h>

#include "sipint.h"


#define hash_1(k,s) (((unsigned long)(k)) % (s))
#define hash_2(k,s) ((s) - 2 - (hash_1((k),(s)) % ((s) - 2)))


/* Prime numbers to use as hash table sizes. */
static unsigned long hash_primes[] = {
    521,        1031,       2053,       4099,
    8209,       16411,      32771,      65537,      131101,     262147,
    524309,     1048583,    2097169,    4194319,    8388617,    16777259,
    33554467,   67108879,   134217757,  268435459,  536870923,  1073741827,
    2147483659U,0
};


static sipHashEntry *newHashTable(unsigned long);
static sipHashEntry *findHashEntry(sipObjectMap *,void *);
static void add_object(sipObjectMap *om, void *addr, sipSimpleWrapper *val);
static void add_aliases(sipObjectMap *om, void *addr, sipSimpleWrapper *val,
        const sipClassTypeDef *base_ctd, const sipClassTypeDef *ctd);
static int remove_object(sipObjectMap *om, void *addr, sipSimpleWrapper *val);
static void remove_aliases(sipObjectMap *om, void *addr, sipSimpleWrapper *val,
        const sipClassTypeDef *base_ctd, const sipClassTypeDef *ctd);
static void reorganiseMap(sipObjectMap *om);
static void *getUnguardedPointer(sipSimpleWrapper *w);


/*
 * Initialise an object map.
 */
void sipOMInit(sipObjectMap *om)
{
    om -> primeIdx = 0;
    om -> unused = om -> size = hash_primes[om -> primeIdx];
    om -> stale = 0;
    om -> hash_array = newHashTable(om -> size);
}


/*
 * Finalise an object map.
 */
void sipOMFinalise(sipObjectMap *om)
{
    sip_api_free(om -> hash_array);
}


/*
 * Allocate and initialise a new hash table.
 */
static sipHashEntry *newHashTable(unsigned long size)
{
    size_t nbytes;
    sipHashEntry *hashtab;

    nbytes = sizeof (sipHashEntry) * size;

    if ((hashtab = (sipHashEntry *)sip_api_malloc(nbytes)) != NULL)
        memset(hashtab,0,nbytes);

    return hashtab;
}


/*
 * Return a pointer to the hash entry that is used, or should be used, for the
 * given C/C++ address.
 */
static sipHashEntry *findHashEntry(sipObjectMap *om,void *key)
{
    unsigned long hash, inc;
    void *hek;

    hash = hash_1(key,om -> size);
    inc = hash_2(key,om -> size);

    while ((hek = om -> hash_array[hash].key) != NULL && hek != key)
        hash = (hash + inc) % om -> size;

    return &om -> hash_array[hash];
}


/*
 * Return the wrapped Python object of a specific type for a C/C++ address or
 * NULL if it wasn't found.
 */
sipSimpleWrapper *sipOMFindObject(sipObjectMap *om, void *key,
        const sipTypeDef *td)
{
    sipHashEntry *he = findHashEntry(om, key);
    sipSimpleWrapper *sw;
    PyTypeObject *py_type = sipTypeAsPyTypeObject(td);

    /* Go through each wrapped object at this address. */
    for (sw = he->first; sw != NULL; sw = sw->next)
    {
        sipSimpleWrapper *unaliased;

        unaliased = (sipIsAlias(sw) ? (sipSimpleWrapper *)sw->data : sw);

        /*
         * If the reference count is 0 then it is in the process of being
         * deleted, so ignore it.  It's not completely clear how this can
         * happen (but it can) because it implies that the garbage collection
         * code is being re-entered (and there are guards in place to prevent
         * this).
         */
        if (Py_REFCNT(unaliased) == 0)
            continue;

        /* Ignore it if the C/C++ address is no longer valid. */
        if (sip_api_get_address(unaliased) == NULL)
            continue;

        /*
         * If this wrapped object is of the given type, or a sub-type of it,
         * then we assume it is the same C++ object.
         */
        if (PyObject_TypeCheck(unaliased, py_type))
            return unaliased;
    }

    return NULL;
}


/*
 * Add a C/C++ address and the corresponding wrapped Python object to the map.
 */
void sipOMAddObject(sipObjectMap *om, sipSimpleWrapper *val)
{
    void *addr = getUnguardedPointer(val);
    const sipClassTypeDef *base_ctd;

    /* Add the object. */
    add_object(om, addr, val);

    /* Add any aliases. */
    base_ctd = (const sipClassTypeDef *)((sipWrapperType *)Py_TYPE(val))->wt_td;
    add_aliases(om, addr, val, base_ctd, base_ctd);
}


/*
 * Add an alias for any address that is different when cast to a super-type.
 */
static void add_aliases(sipObjectMap *om, void *addr, sipSimpleWrapper *val,
        const sipClassTypeDef *base_ctd, const sipClassTypeDef *ctd)
{
    const sipEncodedTypeDef *sup;

    /* See if there are any super-classes. */
    if ((sup = ctd->ctd_supers) != NULL)
    {
        sipClassTypeDef *sup_ctd = sipGetGeneratedClassType(sup, ctd);

        /* Recurse up the hierachy for the first super-class. */
        add_aliases(om, addr, val, base_ctd, sup_ctd);

        /*
         * We only check for aliases for subsequent super-classes because the
         * first one can never need one.
         */
        while (!sup++->sc_flag)
        {
            void *sup_addr;

            sup_ctd = sipGetGeneratedClassType(sup, ctd);

            /* Recurse up the hierachy for the remaining super-classes. */
            add_aliases(om, addr, val, base_ctd, sup_ctd);

            sup_addr = (*base_ctd->ctd_cast)(addr, (sipTypeDef *)sup_ctd);

            if (sup_addr != addr)
            {
                sipSimpleWrapper *alias;

                /* Note that we silently ignore errors. */
                if ((alias = sip_api_malloc(sizeof (sipSimpleWrapper))) != NULL)
                {
                    /*
                     * An alias is basically a bit-wise copy of the Python
                     * object but only to ensure the fields we are subverting
                     * are in the right place.  An alias should never be passed
                     * to the Python API.
                     */
                    *alias = *val;

                    alias->sw_flags = (val->sw_flags & SIP_SHARE_MAP) | SIP_ALIAS;
                    alias->data = val;
                    alias->next = NULL;

                    add_object(om, sup_addr, alias);
                }
            }
        }
    }
}


/*
 * Add a wrapper (which may be an alias) to the map.
 */
static void add_object(sipObjectMap *om, void *addr, sipSimpleWrapper *val)
{
    sipHashEntry *he = findHashEntry(om, addr);

    /*
     * If the bucket is in use then we appear to have several objects at the
     * same address.
     */
    if (he->first != NULL)
    {
        /*
         * This can happen for three reasons.  A variable of one class can be
         * declared at the start of another class.  Therefore there are two
         * objects, of different classes, with the same address.  The second
         * reason is that the old C/C++ object has been deleted by C/C++ but we
         * didn't get to find out for some reason, and a new C/C++ instance has
         * been created at the same address.  The third reason is if we are in
         * the process of deleting a Python object but the C++ object gets
         * wrapped again because the C++ dtor called a method that has been
         * re-implemented in Python.  The absence of the SIP_SHARE_MAP flag
         * tells us that a new C++ instance has just been created and so we
         * know the second reason is the correct one so we mark the old
         * pointers as invalid and reuse the entry.  Otherwise we just add this
         * one to the existing list of objects at this address.
         */
        if (!(val->sw_flags & SIP_SHARE_MAP))
        {
            sipSimpleWrapper *sw = he->first;

            he->first = NULL;

            while (sw != NULL)
            {
                sipSimpleWrapper *next = sw->next;

                if (sipIsAlias(sw))
                {
                    sip_api_free(sw);
                }
                else
                {
                    /*
                     * We are removing it from the map here.  We first have to
                     * call the destructor as the destructor itself might end
                     * up trying to remove the wrapper and its aliases from the
                     * map.
                     */
                    sip_api_instance_destroyed(sw);
                }

                sw = next;
            }
        }

        val->next = he->first;
        he->first = val;

        return;
    }

    /* See if the bucket was unused or stale. */
    if (he->key == NULL)
    {
        he->key = addr;
        om->unused--;
    }
    else
    {
        om->stale--;
    }

    /* Add the rest of the new value. */
    he->first = val;
    val->next = NULL;

    reorganiseMap(om);
}


/*
 * Reorganise a map if it is running short of space.
 */
static void reorganiseMap(sipObjectMap *om)
{
    unsigned long old_size, i;
    sipHashEntry *ohe, *old_tab;

    /* Don't bother if it still has more than 12% available. */
    if (om -> unused > om -> size >> 3)
        return;

    /*
     * If reorganising (ie. making the stale buckets unused) using the same
     * sized table would make 25% available then do that.  Otherwise use a
     * bigger table (if possible).
     */
    if (om -> unused + om -> stale < om -> size >> 2 && hash_primes[om -> primeIdx + 1] != 0)
        om -> primeIdx++;

    old_size = om -> size;
    old_tab = om -> hash_array;

    om -> unused = om -> size = hash_primes[om -> primeIdx];
    om -> stale = 0;
    om -> hash_array = newHashTable(om -> size);

    /* Transfer the entries from the old table to the new one. */
    ohe = old_tab;

    for (i = 0; i < old_size; ++i)
    {
        if (ohe -> key != NULL && ohe -> first != NULL)
        {
            *findHashEntry(om,ohe -> key) = *ohe;
            om -> unused--;
        }

        ++ohe;
    }

    sip_api_free(old_tab);
}


/*
 * Remove a C/C++ object from the table.  Return 0 if it was removed
 * successfully.
 */
int sipOMRemoveObject(sipObjectMap *om, sipSimpleWrapper *val)
{
    void *addr;
    const sipClassTypeDef *base_ctd;

    /* Handle the trivial case. */
    if (sipNotInMap(val))
        return 0;

    if ((addr = getUnguardedPointer(val)) == NULL)
        return 0;

    /* Remove any aliases. */
    base_ctd = (const sipClassTypeDef *)((sipWrapperType *)Py_TYPE(val))->wt_td;
    remove_aliases(om, addr, val, base_ctd, base_ctd);

    /* Remove the object. */
    return remove_object(om, addr, val);
}


/*
 * Remove an alias for any address that is different when cast to a super-type.
 */
static void remove_aliases(sipObjectMap *om, void *addr, sipSimpleWrapper *val,
        const sipClassTypeDef *base_ctd, const sipClassTypeDef *ctd)
{
    const sipEncodedTypeDef *sup;

    /* See if there are any super-classes. */
    if ((sup = ctd->ctd_supers) != NULL)
    {
        sipClassTypeDef *sup_ctd = sipGetGeneratedClassType(sup, ctd);

        /* Recurse up the hierachy for the first super-class. */
        remove_aliases(om, addr, val, base_ctd, sup_ctd);

        /*
         * We only check for aliases for subsequent super-classes because the
         * first one can never need one.
         */
        while (!sup++->sc_flag)
        {
            void *sup_addr;

            sup_ctd = sipGetGeneratedClassType(sup, ctd);

            /* Recurse up the hierachy for the remaining super-classes. */
            remove_aliases(om, addr, val, base_ctd, sup_ctd);

            sup_addr = (*base_ctd->ctd_cast)(addr, (sipTypeDef *)sup_ctd);

            if (sup_addr != addr)
                remove_object(om, sup_addr, val);
        }
    }
}


/*
 * Remove a wrapper from the map.
 */
static int remove_object(sipObjectMap *om, void *addr, sipSimpleWrapper *val)
{
    sipHashEntry *he = findHashEntry(om, addr);
    sipSimpleWrapper **swp;

    for (swp = &he->first; *swp != NULL; swp = &(*swp)->next)
    {
        sipSimpleWrapper *sw, *next;
        int do_remove;

        sw = *swp;
        next = sw->next;

        if (sipIsAlias(sw))
        {
            if (sw->data == val)
            {
                sip_api_free(sw);
                do_remove = TRUE;
            }
            else
            {
                do_remove = FALSE;
            }
        }
        else
        {
            do_remove = (sw == val);
        }

        if (do_remove)
        {
            *swp = next;

            /*
             * If the bucket is now empty then count it as stale.  Note that we
             * do not NULL the key and count it as unused because that might
             * throw out the search for another entry that wanted to go here,
             * found it already occupied, and was put somewhere else.  In other
             * words, searches must be repeatable until we reorganise the
             * table.
             */
            if (he->first == NULL)
                om->stale++;

            return 0;
        }
    }

    return -1;
}


/*
 * Return the unguarded pointer to a C/C++ instance, ie. the pointer was valid
 * but may longer be.
 */
static void *getUnguardedPointer(sipSimpleWrapper *w)
{
    return (w->access_func != NULL) ? w->access_func(w, UnguardedPointer) : w->data;
}
