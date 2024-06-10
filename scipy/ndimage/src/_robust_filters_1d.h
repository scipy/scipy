#include "Python.h"
#include "numpy/arrayobject.h"

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#ifndef HEADER
#define HEADER
struct Mediator // this is used for rank keeping
{
  int *pos;  // index into `heap` for each value
  int *heap; // max/rank/min heap holding indexes into `data`.
  int N;     // allocated size.
  int idx;   // position in circular queue
  int minCt; // count of items in min heap
  int maxCt; // count of items in max heap
};


typedef enum {
  NEAREST = 0,
  WRAP = 1,
  REFLECT = 2,
  MIRROR = 3,
  CONSTANT = 4,
} Mode;

inline void promoteIndex(Mediator* m)
{
   m->idx++;
   if (m->idx == m->N)
   {
        m->idx = 0;
   }
}

extern int mmexchange(Mediator *m, int i, int j);
template <typename T> inline int mmless(T *data, Mediator *m, int i, int j) {
  return (data[m->heap[i]] < data[m->heap[j]]);
}

// swaps items i & j if i < j;  returns true if swapped
template <typename T> inline int mmCmpExch(T *data, Mediator *m, int i, int j) {
  return (mmless(data, m, i, j) && mmexchange(m, i, j));
}

// maintains minheap property for all items below i.
template <typename T>
void minSortDown(T *data, Mediator *m, int i) {
  for (i *= 2; i <= m->minCt; i *= 2) {
    if (i < m->minCt && mmless(data, m, i + 1, i)) {
      ++i;
    }
    if (!mmCmpExch(data, m, i, i / 2)) {
      break;
    }
  }
}

// maintains minheap property for all items above i, including the rank
// returns true if rank changed
template <typename T> inline int minSortUp(T *data, Mediator *m, int i) {
  while (i > 0 && mmCmpExch(data, m, i, i / 2))
    i /= 2;
  return (i == 0);
}


// maintains maxheap property for all items below i. (negative indexes)
template <typename T>
void maxSortDown(T *data, Mediator *m, int i) {
  for (i *= 2; i >= -m->maxCt; i *= 2) {
    if (i > -m->maxCt && mmless(data, m, i, i - 1)) {
      --i;
    }
    if (!mmCmpExch(data, m, i / 2, i)) {
      break;
    }
  }
}

// maintains maxheap property for all items above i, including rank
// returns true if rank changed
template <typename T> inline int maxSortUp(T *data, Mediator *m, int i) {
  while (i < 0 && mmCmpExch(data, m, i / 2, i))
    i /= 2;
  return (i == 0);
}

//maintains rank in O(lg nItems)
template <typename T>
void sortHeap(T* data, Mediator* m, T v, int p, T old)
{
   if (p > 0) //new item is in minHeap
   {  if (v > old) { minSortDown(data, m, p); return; }
      if (minSortUp(data, m, p) && mmCmpExch(data, m, 0, -1)) { maxSortDown(data, m,-1); }
   }
   else if (p < 0) //new item is in maxheap
   {  if ( v < old) {maxSortDown(data, m, p); return; }
      if (maxSortUp(data, m, p) && mmCmpExch(data, m, 1, 0)) { minSortDown(data, m, 1); }
   }
   else //new item is at rank
   {  if (maxSortUp(data, m, -1)) { maxSortDown(data, m, -1); }
      if (minSortUp(data, m, 1)) { minSortDown(data, m, 1); }
   }
}

// Inserts item, maintains rank in O(lg nItems)
template <typename T>
void MediatorInsert(T *data, Mediator *m, T v) {
  int p = m->pos[m->idx];
  T old = data[m->idx];
  data[m->idx] = v;
  promoteIndex(m);
  sortHeap(data, m, v, p, old);
}

extern Mediator* MediatorNew(int nItems, int rank, bool buffer = false);
extern void MediatorDel(Mediator* m, int nItems, int rank, bool buffer = false);
#endif

