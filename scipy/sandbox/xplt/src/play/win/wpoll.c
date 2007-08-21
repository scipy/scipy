/*
 * wpoll.c -- $Id$
 * MS Windows OnIdle methods for boss and worker threads
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"
#include "pstdlib.h"

DWORD w_id_worker = 0;
HANDLE w_worker = 0;

static void (*w_quit)(void)= 0;
static int (*w_quitter)(void)= 0;

static HANDLE *w_inputs = 0;
static struct w_input_item {
  void (*on_input)(void *);
  void *context;
} *w_in_items = 0;

int (*w_msg_hook)(MSG *msg)= 0;

static int w_nputs = 0;
static int w_nmsgs = 0;
static void (**w_msg_callbacks)(MSG *)= 0;

static int w_qclearing = 0;

static int w_poll(UINT timeout);

void
w_initialize(HINSTANCE i, HWND w, void (*wquit)(void),
             int (*wstdinit)(void(**)(char*,long), void(**)(char*,long)),
             HWND (*wparent)(int, int, char *, int))
{
  w_app_instance = i;
  w_main_window = w;
  w_quit = wquit;
  w_stdinit = wstdinit;
  w_parent = wparent;  
}

void
p_quit(void)
{
  w_quit();
}

int
w_on_quit(void)
{
  return w_quitter? w_quitter() : 0;
}

void
p_quitter(int (*on_quit)(void))
{
  w_quitter = on_quit;
}

void
p_qclear(void)
{
  /* causes w_app_msg (the PreTranslate method) to clear worker queue */
  w_qclearing = 1;
}

int
w_work_idle(void)
{
  if (w_poll(0) == WAIT_TIMEOUT) {
    double dtime = p_timeout();
    UINT timeout = (dtime<0. || dtime>2.e9)? INFINITE : (UINT)(1000.*dtime);
    if ((w_nputs+w_nwins) || timeout!=INFINITE) {
      if (w_poll(timeout) == WAIT_TIMEOUT)
        p_on_idle(0);
    } else {
      w_quit();
      return 0;
    }
  }
  return 1;
}

static int
w_poll(UINT timeout)
{
  int status = MsgWaitForMultipleObjects(w_nputs, w_inputs, 0, timeout,
                                         QS_ALLINPUT);
  if (status>=(int)WAIT_OBJECT_0 && status<(int)WAIT_OBJECT_0+w_nputs) {
    status -= WAIT_OBJECT_0;
    w_in_items[status].on_input(w_in_items[status].context);
    status += WAIT_OBJECT_0;
    p_on_idle(1);
  }
  return status;
}

void
w_pollinit(void)
{
  if (!w_inputs) {
    HANDLE heap = GetProcessHeap();
    /* force the system to create the message queue for worker thread */
    MSG msg;
    PeekMessage(&msg, NULL, WM_USER, WM_USER, PM_NOREMOVE);
    w_inputs = HeapAlloc(heap, HEAP_GENERATE_EXCEPTIONS,
                        sizeof(HANDLE *)*8);
    w_in_items = HeapAlloc(heap, HEAP_GENERATE_EXCEPTIONS,
                          sizeof(struct w_input_item *)*8);
    w_msg_callbacks = HeapAlloc(heap, HEAP_GENERATE_EXCEPTIONS,
                               sizeof(void (*)())*8);
    w_id_worker = GetCurrentThreadId();
    DuplicateHandle(GetCurrentProcess(), GetCurrentThread(),
                    GetCurrentProcess(), &w_worker, 0, FALSE,
                    DUPLICATE_SAME_ACCESS);
  }
}

int
w_add_input(HANDLE wait_obj, void (*on_input)(void *), void *context)
{
  int i;
  if (p_signalling) p_abort();
  for (i=0 ; i<w_nputs ; i++) if (w_inputs[i]==wait_obj) break;
  if (on_input) {
    if (i==w_nputs) {
      if (w_nputs && !(w_nputs&7)) {
        HANDLE heap = GetProcessHeap();
        w_inputs=
          HeapReAlloc(heap, HEAP_GENERATE_EXCEPTIONS,
                      w_inputs, sizeof(HANDLE *)*(w_nputs+8));
        w_in_items=
          HeapReAlloc(heap, HEAP_GENERATE_EXCEPTIONS,
                      w_in_items, sizeof(struct w_input_item *)*(w_nputs+8));
      }
      w_inputs[i] = wait_obj;
    }
    /* next three statements are a critical section */
    w_in_items[i].on_input = on_input;
    w_in_items[i].context = context;
    w_nputs++;
  } else {
    /* this whole else clause is a critical section */
    if (w_nputs>0) w_nputs--;
    for (; i<w_nputs ; i++) {
      w_inputs[i] = w_inputs[i+1];
      w_in_items[i].on_input = w_in_items[i+1].on_input;
      w_in_items[i].context = w_in_items[i+1].context;
    }
  }
  return 0;
}

UINT
w_add_msg(void (*on_msg)(MSG *))
{
  UINT msgid = 0;
  if (w_nmsgs && !(w_nmsgs&7)) {
    HANDLE heap = GetProcessHeap();
    w_msg_callbacks=
      HeapReAlloc(heap, HEAP_GENERATE_EXCEPTIONS,
                  w_msg_callbacks, sizeof(void (*)())*(w_nmsgs+8));
  }
  w_msg_callbacks[w_nmsgs] = on_msg;
  msgid = WM_APP + (w_nmsgs++);
  return msgid;
}

int
w_app_msg(MSG *msg)
{
  if (w_msg_hook) {
    int result = w_msg_hook(msg);
    if (result) return result;
  }
  if (p_signalling) w_caught();
  /* handle p_qclear if we are in worker thread */
  if (w_qclearing) {
    int got_quit = 0;
    w_qclearing = 0;
    while (PeekMessage(msg, NULL, 0,0, PM_REMOVE)) {
      if (msg->message==WM_PAINT) ValidateRect(msg->hwnd, 0);
      else if (msg->message==WM_DESTROY) DispatchMessage(msg);
      else if (msg->message==WM_QUIT) got_quit = 1;
    }
    if (got_quit) PostQuitMessage(0);
  }
  if (msg->hwnd==0 &&
      msg->message>=WM_APP && msg->message<WM_APP+(unsigned int)w_nmsgs) {
    w_msg_callbacks[msg->message-WM_APP](msg);
    p_on_idle(1);
    return 1;
  } else {
    return 0;
  }
}
