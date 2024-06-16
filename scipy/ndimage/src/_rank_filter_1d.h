// Copyright (c) 2011 ashelly.myopenid.com under
// <http://www.opensource.org/licenses/mit-license>
// Modified in 2024 by Gideon Genadi Kogan
#include "_robust_filters_utils.h"

template <typename T>
void _rank_filter(T *in_arr, int rank, int arr_len, int win_len, T *out_arr,
                  int mode, T cval, int origin) {
  int i, arr_len_thresh, lim = (win_len - 1) / 2 - origin;
  int lim2 = arr_len - lim;
  Mediator *m = MediatorNew(win_len, rank);
  T *data = new T[win_len];

  switch (mode) {
  case REFLECT:
    for (i = win_len - lim - 1; i > -1; i--) {
      MediatorInsert(data, m, in_arr[i]);
    }
    break;
  case CONSTANT:
    for (i = win_len - lim; i > 0; i--) {
      MediatorInsert(data, m, cval);
    }
    break;
  case NEAREST:
    for (i = win_len - lim; i > 0; i--) {
      MediatorInsert(data, m, in_arr[0]);
    }
    break;
  case MIRROR:
    for (i = win_len - lim; i > 0; i--) {
      MediatorInsert(data, m, in_arr[i]);
    }
    break;
  case WRAP:
    for (i = arr_len - lim - 1 - 2 * origin; i < arr_len; i++) {
      MediatorInsert(data, m, in_arr[i]);
    }
    break;
  }

  for (i = 0; i < lim; i++) {
    MediatorInsert(data, m, in_arr[i]);
  }
  for (i = lim; i < arr_len; i++) {
    MediatorInsert(data, m, in_arr[i]);
    out_arr[i - lim] = data[m->heap[0]];
  }
  switch (mode) {
  case REFLECT:
    arr_len_thresh = arr_len - 1;
    for (i = 0; i < lim; i++) {
      MediatorInsert(data, m, in_arr[arr_len_thresh - i]);
      out_arr[lim2 + i] = data[m->heap[0]];
    }
    break;
  case CONSTANT:
    for (i = 0; i < lim; i++) {
      MediatorInsert(data, m, cval);
      out_arr[lim2 + i] = data[m->heap[0]];
    }
    break;
  case NEAREST:
    arr_len_thresh = arr_len - 1;
    for (i = 0; i < lim; i++) {
      MediatorInsert(data, m, in_arr[arr_len_thresh]);
      out_arr[lim2 + i] = data[m->heap[0]];
    }
    break;
  case MIRROR:
    arr_len_thresh = arr_len - 2;
    for (i = 0; i < lim + 1; i++) {
      MediatorInsert(data, m, in_arr[arr_len_thresh - i]);
      out_arr[lim2 + i] = data[m->heap[0]];
    }
    break;
  case WRAP:
    for (i = 0; i < win_len; i++) {
      MediatorInsert(data, m, in_arr[i]);
      out_arr[lim2 + i] = data[m->heap[0]];
    }
    break;
  }

  MediatorDel(m, win_len);
  delete[] data;
  data = nullptr;
}

