#include "_robust_filters_utils.h"
#include <algorithm>

//Replaces item, maintains rank in O(lg nItems)
template <typename T>
void MediatorReplaceHampel(T* data, Mediator* m, Mediator* m_bottom, Mediator* m_top, T v)
{
   int p = m->pos[m->idx];
   int p_top = m_top->pos[m_top->idx];
   int p_bottom = m_bottom->pos[m_bottom->idx];
   
   T old = data[m->idx];
   data[m->idx] = v;

   // can use the same index for all mediators - later optimization
   promoteIndex(m);
   promoteIndex(m_top);
   promoteIndex(m_bottom);

   sortHeap(data, m, v, p, old);
   sortHeap(data, m_top, v, p_top, old);
   sortHeap(data, m_bottom, v, p_bottom, old);
}

//Replaces item, maintains rank in O(lg nItems)
template <typename T>
void MediatorReplaceRank(T* data, Mediator* m, T v)
{
   int p = m->pos[m->idx];
   T old = data[m->idx];
   data[m->idx] = v;
   promoteIndex(m);
   sortHeap(data, m, v, p, old);
}

// Mediator rank = rank - 1; in O(lg nItems)
template <typename T>
void rank_minus_1(Mediator* m, T* data)
{
   m->minCt++;
   m->maxCt--;
   m->heap[m->minCt] = m->heap[-m->maxCt - 1];
   m->pos[m->heap[m->minCt]] = m->pos[m->heap[-m->maxCt - 1]];
   if (minSortUp(data, m, m->minCt) && mmCmpExch(data, m, 0, -1)) { maxSortDown(data, m, -1); }
}

// Mediator rank = rank + 1; in O(lg nItems)
template <typename T>
void rank_plus_1(Mediator* m, T* data)
{
   m->minCt--;
   m->maxCt++;
   m->heap[-m->maxCt] = m->heap[m->minCt + 1];
   m->pos[m->heap[-m->maxCt]] = m->pos[m->heap[m->minCt + 1]];
   if (maxSortUp(data, m, -m->maxCt) && mmCmpExch(data, m, 1, 0)) { minSortDown(data, m, 1); }
}

template <typename T>
T get_mad(Mediator* m_top, Mediator* m_bottom, T median, T* data, const int win_len){
   T top_order_value = data[m_top->heap[0]];
   T bottom_order_value = data[m_bottom->heap[0]];
   T mad, top_diff, bottom_diff;
   while (top_order_value - median > median - data[m_bottom->heap[-1]])
   {
      if (m_bottom->maxCt == 1)
      {
         bottom_diff = median - data[m_bottom->heap[-1]];
         top_diff = data[m_top->heap[-1]] - median;
         mad = std::max(bottom_diff, top_diff);
         return mad;
      }
      rank_minus_1(m_top, data);
      rank_minus_1(m_bottom, data);
      top_order_value = data[m_top->heap[0]];
      bottom_order_value = data[m_bottom->heap[0]];
   }
   while (median - bottom_order_value > data[m_top->heap[1]] - median)
   {  
      if (m_top->minCt == 1)
         {
            bottom_diff = median - data[m_bottom->heap[1]];
            top_diff = data[m_top->heap[1]] - median;
            mad = std::max(bottom_diff, top_diff);
            return mad;
         }
         rank_plus_1(m_top, data);
         rank_plus_1(m_bottom, data);
         top_order_value = data[m_top->heap[0]];
         bottom_order_value = data[m_bottom->heap[0]];
   }
   top_diff = top_order_value - median;
   bottom_diff = median - bottom_order_value;
   mad = std::max(top_diff, bottom_diff);
   return mad;
}

template <typename T>
void _hampel_filter(T* in_arr, int arr_len, int win_len, T* median, T* mad, T* out_arr, T scale)
{
   if (win_len % 2 == 0)
   {
      printf("Window length must be odd\n");
      exit(1);
   }
   const int rank = (win_len - 1) / 2;
   const int arr_len_thresh = arr_len - 1;
   const int lim = (win_len - 1) / 2;
   int rank_top = win_len - 2;
   int rank_bottom = rank_top - rank;
   int i;
   int i_shift = 0;
   Mediator* m = MediatorNew(win_len, rank);
   Mediator* m_top = MediatorNew(win_len, rank_top, true);
   Mediator* m_bottom = MediatorNew(win_len, rank_bottom, true);
   T* data = (T*)malloc(sizeof(T) * win_len);
   T input_, median_, mad_;
   for (i=win_len - lim; i > 0; i--)
   {
      MediatorReplaceHampel(data, m, m_bottom, m_top, in_arr[0]);
      get_mad(m_top, m_bottom, data[m->heap[0]], data, win_len);
   }
   for (i=0; i < lim; i++)
   {
      MediatorReplaceHampel(data, m, m_bottom, m_top, in_arr[i]);
      get_mad(m_top, m_bottom, data[m->heap[0]], data, win_len);
   }
   for (i=lim; i < arr_len; i++)
   {
      MediatorReplaceHampel(data, m, m_bottom, m_top, in_arr[i]);
      median[i_shift] = median_ = data[m->heap[0]];
      mad[i_shift] = mad_ = get_mad(m_top, m_bottom, median_, data, win_len);
      input_ = in_arr[i_shift];
      out_arr[i_shift] = abs(input_ - median_) > scale * mad_ ? median_ : input_;
      i_shift++;
   }
   for (i=arr_len - lim; i < arr_len; i++)
   {
      MediatorReplaceHampel(data, m, m_bottom, m_top, in_arr[arr_len_thresh]);
      median[i] = median_ = data[m->heap[0]];
      mad[i] = mad_ = get_mad(m_top, m_bottom, median_, data, win_len);
      input_ = in_arr[i];
      out_arr[i] = abs(input_ - median_) > scale * mad_ ? median_ : input_;
   }
   delete[] data;
   data = nullptr;
   MediatorDel(m, win_len);
   MediatorDel(m_top, win_len, true);
   MediatorDel(m_bottom, win_len, true);
}