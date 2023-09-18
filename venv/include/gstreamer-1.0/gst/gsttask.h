/* GStreamer
 * Copyright (C) <1999> Erik Walthinsen <omega@cse.ogi.edu>
 *               <2005> Wim Taymans <wim@fluendo.com>
 *
 * gsttask.h: Streaming tasks
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __GST_TASK_H__
#define __GST_TASK_H__

#include <gst/gstobject.h>
#include <gst/gsttaskpool.h>

G_BEGIN_DECLS

/**
 * GstTaskFunction:
 * @user_data: user data passed to the function
 *
 * A function that will repeatedly be called in the thread created by
 * a #GstTask.
 */
typedef void         (*GstTaskFunction)          (gpointer user_data);

/* --- standard type macros --- */
#define GST_TYPE_TASK                   (gst_task_get_type ())
#define GST_TASK(task)                  (G_TYPE_CHECK_INSTANCE_CAST ((task), GST_TYPE_TASK, GstTask))
#define GST_IS_TASK(task)               (G_TYPE_CHECK_INSTANCE_TYPE ((task), GST_TYPE_TASK))
#define GST_TASK_CLASS(tclass)          (G_TYPE_CHECK_CLASS_CAST ((tclass), GST_TYPE_TASK, GstTaskClass))
#define GST_IS_TASK_CLASS(tclass)       (G_TYPE_CHECK_CLASS_TYPE ((tclass), GST_TYPE_TASK))
#define GST_TASK_GET_CLASS(task)        (G_TYPE_INSTANCE_GET_CLASS ((task), GST_TYPE_TASK, GstTaskClass))
#define GST_TASK_CAST(task)             ((GstTask*)(task))

typedef struct _GstTask GstTask;
typedef struct _GstTaskClass GstTaskClass;
typedef struct _GstTaskPrivate GstTaskPrivate;

/**
 * GstTaskState:
 * @GST_TASK_STARTED: the task is started and running
 * @GST_TASK_STOPPED: the task is stopped
 * @GST_TASK_PAUSED: the task is paused
 *
 * The different states a task can be in
 */
typedef enum {
  GST_TASK_STARTED,
  GST_TASK_STOPPED,
  GST_TASK_PAUSED
} GstTaskState;

/**
 * GST_TASK_STATE:
 * @task: Task to get the state of
 *
 * Get access to the state of the task.
 */
#define GST_TASK_STATE(task)            (GST_TASK_CAST(task)->state)

/**
 * GST_TASK_GET_COND:
 * @task: Task to get the cond of
 *
 * Get access to the cond of the task.
 */
#define GST_TASK_GET_COND(task)         (&GST_TASK_CAST(task)->cond)
/**
 * GST_TASK_WAIT:
 * @task: Task to wait for
 *
 * Wait for the task cond to be signalled
 */
#define GST_TASK_WAIT(task)             g_cond_wait(GST_TASK_GET_COND (task), GST_OBJECT_GET_LOCK (task))
/**
 * GST_TASK_SIGNAL:
 * @task: Task to signal
 *
 * Signal the task cond
 */
#define GST_TASK_SIGNAL(task)           g_cond_signal(GST_TASK_GET_COND (task))
/**
 * GST_TASK_BROADCAST:
 * @task: Task to broadcast
 *
 * Send a broadcast signal to all waiting task conds
 */
#define GST_TASK_BROADCAST(task)        g_cond_broadcast(GST_TASK_GET_COND (task))

/**
 * GST_TASK_GET_LOCK:
 * @task: Task to get the lock of
 *
 * Get access to the task lock.
 */
#define GST_TASK_GET_LOCK(task)         (GST_TASK_CAST(task)->lock)

/**
 * GstTaskThreadFunc:
 * @task: The #GstTask
 * @thread: The #GThread
 * @user_data: user data
 *
 * Custom GstTask thread callback functions that can be installed.
 */
typedef void (*GstTaskThreadFunc) (GstTask *task, GThread *thread, gpointer user_data);

/**
 * GstTask:
 * @state: the state of the task
 * @cond: used to pause/resume the task
 * @lock: The lock taken when iterating the task function
 * @func: the function executed by this task
 * @user_data: user_data passed to the task function
 * @notify: GDestroyNotify for @user_data
 * @running: a flag indicating that the task is running
 *
 * The #GstTask object.
 */
struct _GstTask {
  GstObject      object;

  /*< public >*/ /* with LOCK */
  GstTaskState     state;
  GCond            cond;

  GRecMutex       *lock;

  GstTaskFunction  func;
  gpointer         user_data;
  GDestroyNotify   notify;

  gboolean         running;

  /*< private >*/
  GThread         *thread;

  GstTaskPrivate  *priv;

  gpointer _gst_reserved[GST_PADDING];
};

struct _GstTaskClass {
  GstObjectClass parent_class;

  /*< private >*/
  GstTaskPool *pool;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
void            gst_task_cleanup_all    (void);

GST_API
GType           gst_task_get_type       (void);

GST_API
GstTask*        gst_task_new            (GstTaskFunction func,
                                         gpointer user_data, GDestroyNotify notify);
GST_API
void            gst_task_set_lock       (GstTask *task, GRecMutex *mutex);

GST_API
GstTaskPool *   gst_task_get_pool       (GstTask *task);

GST_API
void            gst_task_set_pool       (GstTask *task, GstTaskPool *pool);

GST_API
void            gst_task_set_enter_callback  (GstTask *task,
                                              GstTaskThreadFunc enter_func,
                                              gpointer user_data,
                                              GDestroyNotify notify);
GST_API
void            gst_task_set_leave_callback  (GstTask *task,
                                              GstTaskThreadFunc leave_func,
                                              gpointer user_data,
                                              GDestroyNotify notify);
GST_API
GstTaskState    gst_task_get_state      (GstTask *task);

GST_API
gboolean        gst_task_set_state      (GstTask *task, GstTaskState state);

GST_API
gboolean        gst_task_start          (GstTask *task);

GST_API
gboolean        gst_task_stop           (GstTask *task);

GST_API
gboolean        gst_task_pause          (GstTask *task);

GST_API
gboolean        gst_task_resume         (GstTask *task);

GST_API
gboolean        gst_task_join           (GstTask *task);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstTask, gst_object_unref)

G_END_DECLS

#endif /* __GST_TASK_H__ */
