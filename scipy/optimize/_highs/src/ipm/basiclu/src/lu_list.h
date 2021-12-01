/*
 * lu_list.h
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Implementation of doubly linked lists (see [1] section 5.5)
 *
 * Maintain nelem elements in nlist doubly linked lists. Each element can belong
 * to zero or one list at a time.
 *
 * The implementation uses arrays
 *
 *  flink[0..nelem+nlist-1],
 *  blink[0..nelem+nlist-1].
 *
 * In each array, the leading nelem entries store links, the trailing nlist
 * entries store heads. That is, for 0 <= i < nelem and 0 <= j < nlist:
 *
 *  flink[i]        next element in the list containing element i
 *  blink[i]        previous element in the list containing element i
 *  flink[nelem+j]  first element in list j
 *  blink[nelem+j]  last element in list j
 *
 * The forward link of the last element in a list points to its flink-head. The
 * backward link of the first element in a list points to its blink-head. For
 * empty lists the heads point to themselves. When an element is not in any list
 * its links point to itself.
 *
 * Optionally the quantity min_list >= 1 can be updated such that lists
 * 1..min_list-1 are empty. Notice that list 0 is not covered by min_list.
 *
 * [1] Istvan Maros, Computational Techniques of the Simplex Method
 *
 * Methods:
 *
 *  lu_list_init
 *  lu_list_add
 *  lu_list_remove
 *  lu_list_move
 *  lu_list_swap
 *
 * The methods are defined in this header file as static inline. This header
 * file must be included after lu_def.h so that lu_int is defined.
 *
 */

/* ==========================================================================
 * lu_list_init
 *
 * Initialize all lists to empty.
 * ========================================================================== */

static inline void lu_list_init(
    lu_int *flink, lu_int *blink, lu_int nelem, lu_int nlist, lu_int *min_list)
{
    lu_int i;
    for (i = 0; i < nelem+nlist; i++) flink[i] = blink[i] = i;
    if (min_list)
        *min_list = MAX(1, nlist);
}


/* ==========================================================================
 * lu_list_add
 *
 * Add element @elem to list @list. @elem must not be in any list already.
 * If list > 0 and min_list != NULL, update *min_list = min(*min_list, list).
 * ========================================================================== */

static inline void lu_list_add(
    lu_int elem, lu_int list, lu_int *flink, lu_int *blink, lu_int nelem,
    lu_int *min_list)
{
    lu_int temp;
    assert(flink[elem] == elem);
    assert(blink[elem] == elem);
    /* append elem to the end of list */
    temp = blink[nelem+list];
    blink[nelem+list] = elem;
    blink[elem] = temp;
    flink[temp] = elem;
    flink[elem] = nelem+list;
    if (list > 0 && min_list && list < *min_list)
        *min_list = list;
}


/* ==========================================================================
 * lu_list_remove
 *
 * Remove element @elem from its list. If @elem was not in a list before,
 * then do nothing.
 * ========================================================================== */

static inline void lu_list_remove(
    lu_int *flink, lu_int *blink, lu_int elem)
{
    flink[blink[elem]] = flink[elem];
    blink[flink[elem]] = blink[elem];
    flink[elem] = elem;
    blink[elem] = elem;
}


/* ==========================================================================
 * lu_list_move
 *
 * Remove element @elem from its list (if in a list) and add it to list @list.
 * ========================================================================== */

static inline void lu_list_move(
    lu_int elem, lu_int list, lu_int *flink, lu_int *blink, lu_int nelem,
    lu_int *min_list)
{
    lu_list_remove(flink, blink, elem);
    lu_list_add(elem, list, flink, blink, nelem, min_list);
}


/* ==========================================================================
 * lu_list_swap
 *
 * Swap elements @e1 and @e2, which both must be in a list. If @e1 and @e2
 * are in the same list, then their positions are swapped. If they are in
 * different lists, then each is moved to the other's list.
 * ========================================================================== */

static inline void lu_list_swap(
    lu_int *flink, lu_int *blink, const lu_int e1, const lu_int e2)
{
    const lu_int e1next = flink[e1];
    const lu_int e2next = flink[e2];
    const lu_int e1prev = blink[e1];
    const lu_int e2prev = blink[e2];

    assert(e1next != e1);       /* must be in a list */
    assert(e2next != e2);

    if (e1next == e2)
    {
        flink[e2] = e1;
        blink[e1] = e2;
        flink[e1prev] = e2;
        blink[e2] = e1prev;
        flink[e1] = e2next;
        blink[e2next] = e1;
    }
    else if (e2next == e1)
    {
        flink[e1] = e2;
        blink[e2] = e1;
        flink[e2] = e1next;
        blink[e1next] = e2;
        flink[e2prev] = e1;
        blink[e1] = e2prev;
    }
    else
    {
        flink[e2] = e1next;
        blink[e1next] = e2;
        flink[e2prev] = e1;
        blink[e1] = e2prev;
        flink[e1prev] = e2;
        blink[e2] = e1prev;
        flink[e1] = e2next;
        blink[e2next] = e1;
    }
}
