/* GStreamer
 *
 * Common code for GStreamer unittests
 *
 * Copyright (C) <2004> Thomas Vander Stichele <thomas at apestaart dot org>
 * Copyright (C) <2008> Thijs Vermeir <thijsvermeir@gmail.com>
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

#ifndef __GST_CHECK_H__
#define __GST_CHECK_H__

#include <signal.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <gst/gst.h>
#include <gst/check/check-prelude.h>

#define CK_DLL_EXP GST_CHECK_API
#include <gst/check/internal-check.h>

G_BEGIN_DECLS

GST_CHECK_API GstDebugCategory *check_debug;
#define GST_CAT_DEFAULT check_debug

/* logging function for tests
 * a test uses g_message() to log a debug line
 * a gst unit test can be run with GST_TEST_DEBUG env var set to see the
 * messages
 */
GST_CHECK_API gboolean _gst_check_threads_running;
GST_CHECK_API gboolean _gst_check_raised_critical;
GST_CHECK_API gboolean _gst_check_raised_warning;
GST_CHECK_API gboolean _gst_check_expecting_log;
GST_CHECK_API gboolean _gst_check_list_tests;

/* global variables used in test methods */
GST_CHECK_API GList * buffers;

GST_CHECK_API GMutex check_mutex;
GST_CHECK_API GCond check_cond;

/**
 * GstCheckABIStruct:
 * @name: The name of the structure
 * @size: The current size of a structure
 * @abi_size: The reference size of the structure
 */
typedef struct
{
  const char *name;
  int size;
  int abi_size;
}
GstCheckABIStruct;

/**
 * GstCheckLogFilter:
 *
 * Opaque structure containing data about a log filter
 * function.
 */
typedef struct _GstCheckLogFilter GstCheckLogFilter;

/**
 * GstCheckLogFilterFunc:
 * @log_domain: the log domain of the message
 * @log_level: the log level of the message
 * @message: the message that has occurred
 * @user_data: user data
 *
 * A function that is called for messages matching the filter added by
 * @gst_check_add_log_filter.
 *
 * Returns: %TRUE if message should be discarded by GstCheck.
 *
 * Since: 1.12
 */
typedef gboolean (*GstCheckLogFilterFunc) (const gchar * log_domain,
    GLogLevelFlags log_level, const gchar * message, gpointer user_data);

GST_CHECK_API
void gst_check_init (int *argc, char **argv[]);

GST_CHECK_API
GstCheckLogFilter * gst_check_add_log_filter (const gchar * log_domain,
    GLogLevelFlags log_level, GRegex * regex, GstCheckLogFilterFunc func,
    gpointer user_data, GDestroyNotify destroy_data);

GST_CHECK_API
void gst_check_remove_log_filter (GstCheckLogFilter * filter);

GST_CHECK_API
void gst_check_clear_log_filter (void);

GST_CHECK_API
GstFlowReturn gst_check_chain_func (GstPad * pad, GstObject * parent, GstBuffer * buffer);

GST_CHECK_API
void gst_check_message_error (GstMessage * message, GstMessageType type,
    GQuark domain, gint code);

GST_CHECK_API
GstElement *gst_check_setup_element (const gchar * factory);

GST_CHECK_API
void gst_check_teardown_element (GstElement * element);

GST_CHECK_API
GstPad *gst_check_setup_src_pad (GstElement * element,
    GstStaticPadTemplate * tmpl);

GST_CHECK_API
GstPad *gst_check_setup_src_pad_from_template (GstElement * element,
    GstPadTemplate * tmpl);

GST_CHECK_API
GstPad * gst_check_setup_src_pad_by_name (GstElement * element,
          GstStaticPadTemplate * tmpl, const gchar *name);

GST_CHECK_API
GstPad * gst_check_setup_src_pad_by_name_from_template (GstElement * element,
          GstPadTemplate * tmpl, const gchar *name);

GST_CHECK_API
GstPad *gst_check_setup_sink_pad (GstElement * element,
    GstStaticPadTemplate * tmpl);

GST_CHECK_API
GstPad *gst_check_setup_sink_pad_from_template (GstElement * element,
    GstPadTemplate * tmpl);

GST_CHECK_API
GstPad * gst_check_setup_sink_pad_by_name (GstElement * element,
          GstStaticPadTemplate * tmpl, const gchar *name);

GST_CHECK_API
GstPad * gst_check_setup_sink_pad_by_name_from_template (GstElement * element,
          GstPadTemplate * tmpl, const gchar *name);

GST_CHECK_API
void gst_check_teardown_pad_by_name (GstElement * element, const gchar *name);

GST_CHECK_API
void gst_check_teardown_src_pad (GstElement * element);

GST_CHECK_API
void gst_check_drop_buffers (void);

GST_CHECK_API
void gst_check_caps_equal (GstCaps * caps1, GstCaps * caps2);

GST_CHECK_API
void gst_check_buffer_data (GstBuffer * buffer, gconstpointer data, gsize size);

GST_CHECK_API
void gst_check_element_push_buffer_list (const gchar * element_name,
    GList * buffer_in, GstCaps * caps_in, GList * buffer_out,
    GstCaps * caps_out, GstFlowReturn last_flow_return);

GST_CHECK_API
void gst_check_element_push_buffer (const gchar * element_name,
    GstBuffer * buffer_in, GstCaps * caps_in, GstBuffer * buffer_out,
    GstCaps *caps_out);

GST_CHECK_API
void gst_check_teardown_sink_pad (GstElement * element);

GST_CHECK_API
void gst_check_abi_list (GstCheckABIStruct list[], gboolean have_abi_sizes);

GST_CHECK_API
gint gst_check_run_suite (Suite * suite, const gchar * name,
    const gchar * fname);

GST_CHECK_API
void gst_check_setup_events (GstPad * srcpad, GstElement * element,
    GstCaps * caps, GstFormat format);

GST_CHECK_API
void gst_check_setup_events_with_stream_id (GstPad * srcpad,
    GstElement * element, GstCaps * caps, GstFormat format,
    const gchar * stream_id);

GST_CHECK_API
void gst_check_objects_destroyed_on_unref (gpointer object_to_unref, gpointer first_object, ...)
  G_GNUC_NULL_TERMINATED;

GST_CHECK_API
void gst_check_object_destroyed_on_unref (gpointer object_to_unref);

#ifndef __GI_SCANNER__

#define fail_unless_message_error(msg, domain, code)            \
gst_check_message_error (msg, GST_MESSAGE_ERROR,                \
  GST_ ## domain ## _ERROR, GST_ ## domain ## _ERROR_ ## code)
#define assert_message_error(m, d, c) fail_unless_message_error(m, d, c)

#ifdef GST_CHECK_TEST_ENVIRONMENT_BEACON
#define GST_DO_CHECK_TEST_ENVIRONMENT \
G_STMT_START {                        \
  if (g_getenv (GST_CHECK_TEST_ENVIRONMENT_BEACON) == NULL) \
    fail ("Test environment not set up correctly! Expected environment " \
       "variable '%s' to be set.", GST_CHECK_TEST_ENVIRONMENT_BEACON); \
} G_STMT_END

#else
#define GST_DO_CHECK_TEST_ENVIRONMENT /* nothing to check */
#endif

/**
 * GST_START_TEST:
 * @__testname: test function name
 *
 * wrapper for checks START_TEST
 */
/**
 * GST_END_TEST:
 *
 * wrapper for checks END_TEST
 */
#define GST_START_TEST(__testname) \
static void __testname (int G_GNUC_UNUSED __i__) \
{\
  GST_DEBUG ("test start"); \
  GST_DO_CHECK_TEST_ENVIRONMENT; \
  tcase_fn_start (""# __testname, __FILE__, __LINE__);

#define GST_END_TEST GST_LOG ("cleaning up tasks"); \
                     gst_task_cleanup_all (); \
                     END_TEST

/* additional fail macros */
/**
 * fail_unless_equals_int:
 * @a: a #gint value or expression
 * @b: a #gint value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to. This
 * macro is for use in unit tests.
 */
#define fail_unless_equals_int(a, b)                                    \
G_STMT_START {                                                          \
  int first = a;                                                        \
  int second = b;                                                       \
  fail_unless(first == second,                                          \
    "'" #a "' (%d) is not equal to '" #b"' (%d)", first, second);       \
} G_STMT_END;
/**
 * assert_equals_int:
 * @a: a #gint value or expression
 * @b: a #gint value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to. This
 * macro is for use in unit tests.
 */
#define assert_equals_int(a, b) fail_unless_equals_int(a, b)

/**
 * fail_unless_equals_int_hex:
 * @a: a #gint value or expression
 * @b: a #gint value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to in
 * hexadecimal format. This macro is for use in unit tests.
 *
 * Since: 1.2
 */
#define fail_unless_equals_int_hex(a, b)                                \
G_STMT_START {								\
  int first = a;							\
  int second = b;							\
  fail_unless(first == second,						\
    "'" #a "' (0x%08x) is not equal to '" #b"' (0x%08x)", first, second);\
} G_STMT_END;

/**
 * assert_equals_int_hex:
 * @a: a #gint value or expression
 * @b: a #gint value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to in
 * hexadecimal format. This macro is for use in unit tests.
 *
 * Since: 1.2
 */
#define assert_equals_int_hex(a, b) fail_unless_equals_int_hex(a, b)

/**
 * fail_unless_equals_int64:
 * @a: a #gint64 value or expression
 * @b: a #gint64 value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to. This
 * macro is for use in unit tests.
 */
#define fail_unless_equals_int64(a, b)                                  \
G_STMT_START {                                                          \
  gint64 first = a;                                                     \
  gint64 second = b;                                                    \
  fail_unless(first == second,                                          \
    "'" #a "' (%" G_GINT64_FORMAT") is not equal to '" #b"' (%"         \
    G_GINT64_FORMAT")", first, second);                                 \
} G_STMT_END;
/**
 * assert_equals_int64:
 * @a: a #gint64 value or expression
 * @b: a #gint64 value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to. This
 * macro is for use in unit tests.
 */
#define assert_equals_int64(a, b) fail_unless_equals_int64(a, b)

/**
 * fail_unless_equals_int64_hex:
 * @a: a #gint64 value or expression
 * @b: a #gint64 value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to in
 * hexadecimal format. This macro is for use in unit tests.
 *
 * Since: 1.2
 */
#define fail_unless_equals_int64_hex(a, b)                              \
G_STMT_START {								\
  gint64 first = a;							\
  gint64 second = b;							\
  fail_unless(first == second,						\
    "'" #a "' (0x%016x) is not equal to '" #b"' (0x%016x)", first, second);\
} G_STMT_END;
/**
 * assert_equals_int64_hex:
 * @a: a #gint64 value or expression
 * @b: a #gint64 value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to in
 * hexadecimal format. This macro is for use in unit tests.
 *
 * Since: 1.2
 */
#define assert_equals_int64_hex(a,b) fail_unless_equals_int64_hex(a,b)

/**
 * fail_unless_equals_uint64:
 * @a: a #guint64 value or expression
 * @b: a #guint64 value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to. This
 * macro is for use in unit tests.
 */
#define fail_unless_equals_uint64(a, b)                                 \
G_STMT_START {                                                          \
  guint64 first = a;                                                    \
  guint64 second = b;                                                   \
  fail_unless(first == second,                                          \
    "'" #a "' (%" G_GUINT64_FORMAT ") is not equal to '" #b"' (%"       \
    G_GUINT64_FORMAT ")", first, second);                               \
} G_STMT_END;
/**
 * assert_equals_uint64:
 * @a: a #guint64 value or expression
 * @b: a #guint64 value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to. This
 * macro is for use in unit tests.
 */
#define assert_equals_uint64(a, b) fail_unless_equals_uint64(a, b)

/**
 * fail_unless_equals_uint64_hex:
 * @a: a #gint64 value or expression
 * @b: a #gint64 value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to in
 * hexadecimal format. This macro is for use in unit tests.
 *
 * Since: 1.2
 */
#define fail_unless_equals_uint64_hex(a, b)                             \
G_STMT_START {								\
  guint64 first = a;							\
  guint64 second = b;							\
  fail_unless(first == second,						\
    "'" #a "' (0x%016x) is not equal to '" #b"' (0x%016x)", first, second);\
} G_STMT_END;
/**
 * assert_equals_uint64_hex:
 * @a: a #guint64 value or expression
 * @b: a #guint64 value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to in
 * hexadecimal format. This macro is for use in unit tests.
 *
 * Since: 1.2
 */
#define assert_equals_uint64_hex(a,b) fail_unless_equals_uint64_hex(a,b)

/**
 * fail_unless_equals_string:
 * @a: a string literal or expression
 * @b: a string literal or expression
 *
 * This macro checks that @a and @b are equal (as per g_strcmp0()) and aborts if
 * this is not the case, printing both expressions and the values they
 * evaluated to. This macro is for use in unit tests.
 */
#define fail_unless_equals_string(a, b)                             \
G_STMT_START {                                                      \
  const gchar * first = a;                                          \
  const gchar * second = b;                                         \
  fail_unless(g_strcmp0 (first, second) == 0,                          \
    "'" #a "' (%s) is not equal to '" #b"' (%s)", first, second);   \
} G_STMT_END;
/**
 * assert_equals_string:
 * @a: a string literal or expression
 * @b: a string literal or expression
 *
 * This macro checks that @a and @b are equal (as per g_strcmp0()) and aborts if
 * this is not the case, printing both expressions and the values they
 * evaluated to. This macro is for use in unit tests.
 */
#define assert_equals_string(a, b) fail_unless_equals_string(a, b)

/**
 * fail_unless_equals_float:
 * @a: a #gdouble or #gfloat value or expression
 * @b: a #gdouble or #gfloat value or expression
 *
 * This macro checks that @a and @b are (almost) equal and aborts if this
 * is not the case, printing both expressions and the values they evaluated
 * to. This macro is for use in unit tests.
 */
#define fail_unless_equals_float(a, b)                            \
G_STMT_START {                                                    \
  double first = a;                                               \
  double second = b;                                              \
  /* This will only work for 'normal' values and values around 0, \
   * which should be good enough for our purposes here */         \
  fail_unless(fabs (first - second) < 0.0000001,                  \
    "'" #a "' (%g) is not equal to '" #b "' (%g)", first, second);\
} G_STMT_END;

/**
 * assert_equals_float:
 * @a: a #gdouble or #gfloat value or expression
 * @b: a #gdouble or #gfloat value or expression
 *
 * This macro checks that @a and @b are (almost) equal and aborts if this
 * is not the case, printing both expressions and the values they evaluated
 * to. This macro is for use in unit tests.
 */
#define assert_equals_float(a, b) fail_unless_equals_float(a, b)

/**
 * fail_unless_equals_pointer:
 * @a: a pointer value or expression
 * @b: a pointer value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this
 * is not the case, printing both expressions and the values they
 * evaluated to. This macro is for use in unit tests.
 *
 * Since: 1.2
 */
#define fail_unless_equals_pointer(a, b)                          \
G_STMT_START {                                                    \
  gpointer first = a;                                             \
  gpointer second = b;                                            \
  fail_unless(first == second,                                    \
    "'" #a "' (%p) is not equal to '" #b "' (%p)", first, second);\
} G_STMT_END;

/**
 * assert_equals_pointer:
 * @a: a pointer value or expression
 * @b: a pointer value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this
 * is not the case, printing both expressions and the values they
 * evaluated to. This macro is for use in unit tests.
 *
 * Since: 1.2
 */
#define assert_equals_pointer(a, b) fail_unless_equals_pointer(a, b)

/**
 * fail_unless_equals_clocktime:
 * @a: a #GstClockTime value or expression
 * @b: a #GstClockTime value or expression
 *
 * This macro checks that @a and @b are equal and aborts if this is not the
 * case, printing both expressions and the values they evaluated to. This
 * macro is for use in unit tests.
 */
#define fail_unless_equals_clocktime(a, b)                              \
G_STMT_START {                                                          \
  GstClockTime first = a;                                                        \
  GstClockTime second = b;                                                       \
  fail_unless(first == second,                                          \
    "'" #a "' (%" GST_TIME_FORMAT") is not equal to '" #b"' (%" GST_TIME_FORMAT")", \
      GST_TIME_ARGS (first), GST_TIME_ARGS (second));       \
} G_STMT_END;

/***
 * thread test macros and variables
 */
GST_CHECK_API GList *thread_list;
GST_CHECK_API GMutex mutex;
GST_CHECK_API GCond start_cond;       /* used to notify main thread of thread startups */
GST_CHECK_API GCond sync_cond;        /* used to synchronize all threads and main thread */

#define MAIN_START_THREADS(count, function, data)               \
MAIN_INIT();                                                    \
MAIN_START_THREAD_FUNCTIONS(count, function, data);             \
MAIN_SYNCHRONIZE();

#define MAIN_INIT()                     \
G_STMT_START {                          \
  g_mutex_init (&mutex);                \
  g_cond_init (&start_cond);            \
  g_cond_init (&sync_cond);             \
  _gst_check_threads_running = TRUE;    \
} G_STMT_END;

#define MAIN_START_THREAD_FUNCTIONS(count, function, data)      \
G_STMT_START {                                                  \
  int i;                                                        \
  for (i = 0; i < count; ++i) {                                 \
    MAIN_START_THREAD_FUNCTION (i, function, data);             \
  }                                                             \
} G_STMT_END;

#define MAIN_START_THREAD_FUNCTION(i, function, data)           \
G_STMT_START {                                                  \
    GThread *thread = NULL;                                     \
    GST_DEBUG ("MAIN: creating thread %d", i);                  \
    g_mutex_lock (&mutex);                                      \
    thread = g_thread_try_new ("gst-check",                     \
        (GThreadFunc) function, data, NULL);                    \
    /* wait for thread to signal us that it's ready */          \
    GST_DEBUG ("MAIN: waiting for thread %d", i);               \
    g_cond_wait (&start_cond, &mutex);                          \
    g_mutex_unlock (&mutex);                                    \
                                                                \
    thread_list = g_list_append (thread_list, thread);          \
} G_STMT_END;


#define MAIN_SYNCHRONIZE()              \
G_STMT_START {                          \
  GST_DEBUG ("MAIN: synchronizing");    \
  g_cond_broadcast (&sync_cond);        \
  GST_DEBUG ("MAIN: synchronized");     \
} G_STMT_END;

#define MAIN_STOP_THREADS()                                     \
G_STMT_START {                                                  \
  _gst_check_threads_running = FALSE;                           \
                                                                \
  /* join all threads */                                        \
  GST_DEBUG ("MAIN: joining");                                  \
  g_list_foreach (thread_list, (GFunc) g_thread_join, NULL);    \
  g_list_free (thread_list);                                    \
  thread_list = NULL;                                           \
  g_mutex_clear (&mutex);                                       \
  g_cond_clear (&start_cond);                                   \
  g_cond_clear (&sync_cond);                                    \
  GST_DEBUG ("MAIN: joined");                                   \
} G_STMT_END;

#define THREAD_START()                                          \
THREAD_STARTED();                                               \
THREAD_SYNCHRONIZE();

#define THREAD_STARTED()                                        \
G_STMT_START {                                                  \
  /* signal main thread that we started */                      \
  GST_DEBUG ("THREAD %p: started", g_thread_self ());           \
  g_mutex_lock (&mutex);                                        \
  g_cond_signal (&start_cond);                                  \
} G_STMT_END;

#define THREAD_SYNCHRONIZE()                                    \
G_STMT_START {                                                  \
  /* synchronize everyone */                                    \
  GST_DEBUG ("THREAD %p: syncing", g_thread_self ());           \
  fail_if (g_mutex_trylock (&mutex),                            \
      "bug in unit test, mutex should be locked at this point");\
  g_cond_wait (&sync_cond, &mutex);                             \
  GST_DEBUG ("THREAD %p: synced", g_thread_self ());            \
  g_mutex_unlock (&mutex);                                      \
} G_STMT_END;

#define THREAD_SWITCH()                                         \
G_STMT_START {                                                  \
  g_thread_yield ();                                            \
} G_STMT_END;

#define THREAD_TEST_RUNNING()   (!!_gst_check_threads_running)

/* additional assertions */

#if GST_DISABLE_GLIB_CHECKS
#define ASSERT_CRITICAL(code)
#else
#define ASSERT_CRITICAL(code)                                   \
G_STMT_START {                                                  \
  _gst_check_expecting_log = TRUE;                              \
  _gst_check_raised_critical = FALSE;                           \
  code;                                                         \
  if (!_gst_check_raised_critical)                              \
    _ck_assert_failed (__FILE__, __LINE__,                      \
        "Expected g_critical, got nothing", NULL);              \
  _gst_check_expecting_log = FALSE;                             \
} G_STMT_END
#endif /* GST_DISABLE_GLIB_CHECKS */

#define ASSERT_WARNING(code)                                    \
G_STMT_START {                                                  \
  _gst_check_expecting_log = TRUE;                              \
  _gst_check_raised_warning = FALSE;                            \
  code;                                                         \
  if (!_gst_check_raised_warning)                               \
    _ck_assert_failed (__FILE__, __LINE__,                      \
        "Expected g_warning, got nothing", NULL);               \
  _gst_check_expecting_log = FALSE;                             \
} G_STMT_END


#define ASSERT_OBJECT_REFCOUNT(object, name, value)             \
G_STMT_START {                                                  \
  int rc;                                                       \
  rc = GST_OBJECT_REFCOUNT_VALUE (object);                      \
  fail_unless (rc == value,                                     \
      "%s (%p) refcount is %d instead of %d",                   \
      name, object, rc, value);                                 \
} G_STMT_END

#define ASSERT_OBJECT_REFCOUNT_BETWEEN(object, name, lower, upper)      \
G_STMT_START {                                                          \
  int rc = GST_OBJECT_REFCOUNT_VALUE (object);                          \
  int lo = lower;                                                       \
  int hi = upper;                                                       \
                                                                        \
  fail_unless (rc >= lo,                                                \
      "%s (%p) refcount %d is smaller than %d",                         \
      name, object, rc, lo);                                            \
  fail_unless (rc <= hi,                                                \
      "%s (%p) refcount %d is bigger than %d",                          \
      name, object, rc, hi);                                            \
} G_STMT_END


#define ASSERT_CAPS_REFCOUNT(caps, name, value)                 \
        ASSERT_MINI_OBJECT_REFCOUNT(caps, name, value)

#define ASSERT_BUFFER_REFCOUNT(buffer, name, value)             \
        ASSERT_MINI_OBJECT_REFCOUNT(buffer, name, value)

#define ASSERT_MINI_OBJECT_REFCOUNT(miniobj, name, value)       \
G_STMT_START {                                                  \
  int rc;                                                       \
  rc = GST_MINI_OBJECT_REFCOUNT_VALUE (miniobj);                \
  fail_unless (rc == value,                                     \
               name " (%p) refcount is %d instead of %d", miniobj, rc, value); \
} G_STMT_END

#define ASSERT_SET_STATE(element, state, ret)                   \
fail_unless (gst_element_set_state (GST_ELEMENT(element),       \
  state) == ret,                                                \
  "could not change state to " #state);

#define GST_CHECK_MAIN(name)                                    \
int main (int argc, char **argv)                                \
{                                                               \
  Suite *s;                                                     \
  gst_check_init (&argc, &argv);                                \
  s = name ## _suite ();                                        \
  return gst_check_run_suite (s, # name, __FILE__);             \
}

/* Hack to allow run-time selection of unit tests to run via the
 * GST_CHECKS environment variable (test function names globs, comma
 * separated), or GST_CHECKS_IGNORE with the same semantics */

GST_CHECK_API
gboolean _gst_check_run_test_func (const gchar * func_name);

static inline void
__gst_tcase_add_test (TCase * tc, TFun tf, const char * fname, int signal,
    int allowed_exit_value, int start, int end)
{
    if (_gst_check_list_tests) {
        g_print ("Test: %s\n", fname);
        return;
    }

    if (_gst_check_run_test_func (fname)) {
        _tcase_add_test (tc, tf, fname, signal, allowed_exit_value, start, end);
    }
}

#define _tcase_add_test __gst_tcase_add_test

/* A special variant to add broken tests. These are normally skipped, but can be
 * forced to run via GST_CHECKS */
#define tcase_skip_broken_test(chain,test_func) \
G_STMT_START {                                                  \
  const char *env = g_getenv ("GST_CHECKS");                    \
                                                                \
  if (env != NULL && g_pattern_match_simple (env, G_STRINGIFY (test_func))) {   \
    tcase_add_test(chain,test_func);                            \
  } else {                                                      \
    g_printerr ("FIXME: skipping test %s because it's broken\n", G_STRINGIFY (test_func)); \
  } \
} G_STMT_END

#define tcase_skip_broken_loop_test(chain,test_func,a,b)        \
  tcase_skip_broken_test (chain, test_func)

#endif /* !__GI_SCANNER__ */

G_END_DECLS

#endif /* __GST_CHECK_H__ */
