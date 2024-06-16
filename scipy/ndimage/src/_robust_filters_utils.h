#include <stdlib.h>
#include <stdio.h>

#ifndef ROBUST_FILTERS_UTILS_H
#define ROBUST_FILTERS_UTILS_H

struct Mediator//this is used for rank keeping
{
   int*  pos;   //index into `heap` for each value
   int*  heap;  //max/rank/min heap holding indexes into `data`.
   int   N;     //allocated size.
   int   idx;   //position in circular queue
   int   minCt; //count of items in min heap
   int   maxCt; //count of items in max heap
};

typedef enum {
   NEAREST = 0,
   WRAP = 1,
   REFLECT = 2,
   MIRROR = 3,
   CONSTANT = 4,
} Mode;

/*--- Helper Functions ---*/
//swaps items i&j in heap, maintains indexes
int mmexchange(Mediator* m, int i, int j)
{
   int t = m->heap[i];
   m->heap[i] = m->heap[j];
   m->heap[j] = t;
   m->pos[m->heap[i]] = i;
   m->pos[m->heap[j]] = j;
   return 1;
}

//returns 1 if heap[i] < heap[j]
template <typename T>
inline int mmless(T* data, Mediator* m, int i, int j)
{
   return (data[m->heap[i]] < data[m->heap[j]]);
}

//swaps items i & j if i < j;  returns true if swapped
template <typename T>
inline int mmCmpExch(T* data, Mediator* m, int i, int j)
{
   return (mmless(data, m,i,j) && mmexchange(m,i,j));
}

//maintains minheap property for all items below i.
template <typename T>
void minSortDown(T* data, Mediator* m, int i)
{
   for (i*=2; i <= m->minCt; i*=2)
   {  if (i < m->minCt && mmless(data, m, i+1, i)) { ++i; }
      if (!mmCmpExch(data, m, i, i/2)) { break; }
   }
}

//maintains maxheap property for all items below i. (negative indexes)
template <typename T>
void maxSortDown(T* data, Mediator* m, int i)
{
   for (i*=2; i >= -m->maxCt; i*=2)
   {  if (i > -m->maxCt && mmless(data, m, i, i-1)) { --i;}
      if (!mmCmpExch(data, m, i/2, i)) { break; }
   }
}

//maintains minheap property for all items above i, including the rank
//returns true if rank changed
template <typename T>
inline int minSortUp(T* data, Mediator* m, int i)
{
   while (i>0 && mmCmpExch(data, m, i, i/2)) {i/=2;};
   return (i==0);
}

//maintains maxheap property for all items above i, including rank
//returns true if rank changed
template <typename T>
inline int maxSortUp(T* data, Mediator* m, int i)
{
   while (i<0 && mmCmpExch(data, m, i/2, i)) {i/=2;};
   return (i==0);
}

//maintains rank in O(lg nItems)
template <typename T>
inline void sortHeap(T* data, Mediator* m, T v, int p, T old)
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

inline void promoteIndex(Mediator* m)
{
   m->idx++;
   if (m->idx == m->N)
   {
        m->idx = 0;
   }
}

//creates new Mediator: to calculate `nItems` running rank.
Mediator* MediatorNew(int nItems, int rank, bool buffer = false)
{
   Mediator* m =  new Mediator;
   int size = buffer ? 2 * nItems - 1 : nItems;
   m->pos = new int[size];
   m->heap = new int[size];
   if ((m == nullptr)||(m->pos == nullptr)||(m->heap == nullptr)){printf("out of memory\n"); exit(1);}
   m->heap += buffer ? rank + (nItems - 1) / 2 :rank; //points to rank
   m->N = nItems;
   m->idx = 0;
   m->minCt = nItems - rank - 1;
   m->maxCt = rank;
   while (nItems--)
   {
      m->pos[nItems]= nItems - rank;
      m->heap[m->pos[nItems]]=nItems;
   }
   return m;
}

template <typename T>
void MediatorInsert(T *data, Mediator *m, T v) {
  int p = m->pos[m->idx];
  T old = data[m->idx];
  data[m->idx] = v;
  promoteIndex(m);
  sortHeap(data, m, v, p, old);
}

// Deletes Mediator
void MediatorDel(Mediator* m, int nItems, bool buffer = false)
{
    m->heap -= buffer ? m->maxCt + (nItems - 1) / 2 :m->maxCt; //points to rank
    delete[] m->heap;
    m->heap = nullptr;
    delete[] m->pos;
    m->pos = nullptr;
    delete m;
    m = nullptr;
}
#endif