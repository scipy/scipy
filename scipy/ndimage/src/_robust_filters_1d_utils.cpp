#include "_robust_filters_1d.h"

/*--- Helper Functions ---*/
// swaps items i&j in heap, maintains indexes
int mmexchange(Mediator *m, int i, int j) {
  int t = m->heap[i];
  m->heap[i] = m->heap[j];
  m->heap[j] = t;
  m->pos[m->heap[i]] = i;
  m->pos[m->heap[j]] = j;
  return 1;
}

//creates new Mediator: to calculate `nItems` running rank.
Mediator* MediatorNew(int nItems, int rank, bool buffer)
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

// Deletes Mediator
void MediatorDel(Mediator* m, int nItems, int rank, bool buffer)
{
    m->heap -= buffer ? rank + (nItems - 1) / 2 :rank; //points to rank
    delete[] m->heap;
    m->heap = nullptr;
    delete[] m->pos;
    m->pos = nullptr;
    delete m;
    m = nullptr;
}

