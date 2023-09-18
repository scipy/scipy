/* GStreamer
 * Copyright (C) 2004 Thomas Vander Stichele <thomas at apestaart dot org>
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

#ifndef __GST_ERROR_H__
#define __GST_ERROR_H__

#include <glib.h>
#include <glib-object.h>
#include <errno.h>

G_BEGIN_DECLS
/*
 * we define FIXME error domains:
 * GST_CORE_ERROR
 * GST_LIBRARY_ERROR
 * GST_RESOURCE_ERROR
 * GST_STREAM_ERROR
 *
 * Check GError API docs for rationale for naming.
 */
/**
 * GstCoreError:
 * @GST_CORE_ERROR_FAILED: a general error which doesn't fit in any other
 * category.  Make sure you add a custom message to the error call.
 * @GST_CORE_ERROR_TOO_LAZY: do not use this except as a placeholder for
 * deciding where to go while developing code.
 * @GST_CORE_ERROR_NOT_IMPLEMENTED: use this when you do not want to implement
 * this functionality yet.
 * @GST_CORE_ERROR_STATE_CHANGE: used for state change errors.
 * @GST_CORE_ERROR_PAD: used for pad-related errors.
 * @GST_CORE_ERROR_THREAD: used for thread-related errors.
 * @GST_CORE_ERROR_NEGOTIATION: used for negotiation-related errors.
 * @GST_CORE_ERROR_EVENT: used for event-related errors.
 * @GST_CORE_ERROR_SEEK: used for seek-related errors.
 * @GST_CORE_ERROR_CAPS: used for caps-related errors.
 * @GST_CORE_ERROR_TAG: used for negotiation-related errors.
 * @GST_CORE_ERROR_MISSING_PLUGIN: used if a plugin is missing.
 * @GST_CORE_ERROR_CLOCK: used for clock related errors.
 * @GST_CORE_ERROR_DISABLED: used if functionality has been disabled at
 *                           compile time.
 * @GST_CORE_ERROR_NUM_ERRORS: the number of core error types.
 *
 * Core errors are errors inside the core GStreamer library.
 */
/* FIXME: should we divide in numerical blocks so we can easily add
          for example PAD errors later ? */
typedef enum
{
  GST_CORE_ERROR_FAILED = 1,
  GST_CORE_ERROR_TOO_LAZY,
  GST_CORE_ERROR_NOT_IMPLEMENTED,
  GST_CORE_ERROR_STATE_CHANGE,
  GST_CORE_ERROR_PAD,
  GST_CORE_ERROR_THREAD,
  GST_CORE_ERROR_NEGOTIATION,
  GST_CORE_ERROR_EVENT,
  GST_CORE_ERROR_SEEK,
  GST_CORE_ERROR_CAPS,
  GST_CORE_ERROR_TAG,
  GST_CORE_ERROR_MISSING_PLUGIN,
  GST_CORE_ERROR_CLOCK,
  GST_CORE_ERROR_DISABLED,
  GST_CORE_ERROR_NUM_ERRORS
} GstCoreError;

/**
 * GstLibraryError:
 * @GST_LIBRARY_ERROR_FAILED: a general error which doesn't fit in any other
 * category.  Make sure you add a custom message to the error call.
 * @GST_LIBRARY_ERROR_TOO_LAZY: do not use this except as a placeholder for
 * deciding where to go while developing code.
 * @GST_LIBRARY_ERROR_INIT: used when the library could not be opened.
 * @GST_LIBRARY_ERROR_SHUTDOWN: used when the library could not be closed.
 * @GST_LIBRARY_ERROR_SETTINGS: used when the library doesn't accept settings.
 * @GST_LIBRARY_ERROR_ENCODE: used when the library generated an encoding error.
 * @GST_LIBRARY_ERROR_NUM_ERRORS: the number of library error types.
 *
 * Library errors are for errors from the library being used by elements
 * (initializing, finalizing, settings, ...)
 */
typedef enum
{
  GST_LIBRARY_ERROR_FAILED = 1,
  GST_LIBRARY_ERROR_TOO_LAZY,
  GST_LIBRARY_ERROR_INIT,
  GST_LIBRARY_ERROR_SHUTDOWN,
  GST_LIBRARY_ERROR_SETTINGS,
  GST_LIBRARY_ERROR_ENCODE,
  GST_LIBRARY_ERROR_NUM_ERRORS
} GstLibraryError;

/**
 * GstResourceError:
 * @GST_RESOURCE_ERROR_FAILED: a general error which doesn't fit in any other
 * category.  Make sure you add a custom message to the error call.
 * @GST_RESOURCE_ERROR_TOO_LAZY: do not use this except as a placeholder for
 * deciding where to go while developing code.
 * @GST_RESOURCE_ERROR_NOT_FOUND: used when the resource could not be found.
 * @GST_RESOURCE_ERROR_BUSY: used when resource is busy.
 * @GST_RESOURCE_ERROR_OPEN_READ: used when resource fails to open for reading.
 * @GST_RESOURCE_ERROR_OPEN_WRITE: used when resource fails to open for writing.
 * @GST_RESOURCE_ERROR_OPEN_READ_WRITE: used when resource cannot be opened for
 * both reading and writing, or either (but unspecified which).
 * @GST_RESOURCE_ERROR_CLOSE: used when the resource can't be closed.
 * @GST_RESOURCE_ERROR_READ: used when the resource can't be read from.
 * @GST_RESOURCE_ERROR_WRITE: used when the resource can't be written to.
 * @GST_RESOURCE_ERROR_SEEK: used when a seek on the resource fails.
 * @GST_RESOURCE_ERROR_SYNC: used when a synchronize on the resource fails.
 * @GST_RESOURCE_ERROR_SETTINGS: used when settings can't be manipulated on.
 * @GST_RESOURCE_ERROR_NO_SPACE_LEFT: used when the resource has no space left.
 * @GST_RESOURCE_ERROR_NOT_AUTHORIZED: used when the resource can't be opened
 *                                     due to missing authorization.
 *                                     (Since: 1.2.4)
 * @GST_RESOURCE_ERROR_NUM_ERRORS: the number of resource error types.
 *
 * Resource errors are for any resource used by an element:
 * memory, files, network connections, process space, ...
 * They're typically used by source and sink elements.
 */
typedef enum
{
  GST_RESOURCE_ERROR_FAILED = 1,
  GST_RESOURCE_ERROR_TOO_LAZY,
  GST_RESOURCE_ERROR_NOT_FOUND,
  GST_RESOURCE_ERROR_BUSY,
  GST_RESOURCE_ERROR_OPEN_READ,
  GST_RESOURCE_ERROR_OPEN_WRITE,
  GST_RESOURCE_ERROR_OPEN_READ_WRITE,
  GST_RESOURCE_ERROR_CLOSE,
  GST_RESOURCE_ERROR_READ,
  GST_RESOURCE_ERROR_WRITE,
  GST_RESOURCE_ERROR_SEEK,
  GST_RESOURCE_ERROR_SYNC,
  GST_RESOURCE_ERROR_SETTINGS,
  GST_RESOURCE_ERROR_NO_SPACE_LEFT,
  GST_RESOURCE_ERROR_NOT_AUTHORIZED,
  GST_RESOURCE_ERROR_NUM_ERRORS
} GstResourceError;

/**
 * GstStreamError:
 * @GST_STREAM_ERROR_FAILED: a general error which doesn't fit in any other
 * category.  Make sure you add a custom message to the error call.
 * @GST_STREAM_ERROR_TOO_LAZY: do not use this except as a placeholder for
 * deciding where to go while developing code.
 * @GST_STREAM_ERROR_NOT_IMPLEMENTED: use this when you do not want to implement
 * this functionality yet.
 * @GST_STREAM_ERROR_TYPE_NOT_FOUND: used when the element doesn't know the
 * stream's type.
 * @GST_STREAM_ERROR_WRONG_TYPE: used when the element doesn't handle this type
 * of stream.
 * @GST_STREAM_ERROR_CODEC_NOT_FOUND: used when there's no codec to handle the
 * stream's type.
 * @GST_STREAM_ERROR_DECODE: used when decoding fails.
 * @GST_STREAM_ERROR_ENCODE: used when encoding fails.
 * @GST_STREAM_ERROR_DEMUX: used when demuxing fails.
 * @GST_STREAM_ERROR_MUX: used when muxing fails.
 * @GST_STREAM_ERROR_FORMAT: used when the stream is of the wrong format
 * (for example, wrong caps).
 * @GST_STREAM_ERROR_DECRYPT: used when the stream is encrypted and can't be
 * decrypted because this is not supported by the element.
 * @GST_STREAM_ERROR_DECRYPT_NOKEY: used when the stream is encrypted and
 * can't be decrypted because no suitable key is available.
 * @GST_STREAM_ERROR_NUM_ERRORS: the number of stream error types.
 *
 * Stream errors are for anything related to the stream being processed:
 * format errors, media type errors, ...
 * They're typically used by decoders, demuxers, converters, ...
 */
typedef enum
{
  GST_STREAM_ERROR_FAILED = 1,
  GST_STREAM_ERROR_TOO_LAZY,
  GST_STREAM_ERROR_NOT_IMPLEMENTED,
  GST_STREAM_ERROR_TYPE_NOT_FOUND,
  GST_STREAM_ERROR_WRONG_TYPE,
  GST_STREAM_ERROR_CODEC_NOT_FOUND,
  GST_STREAM_ERROR_DECODE,
  GST_STREAM_ERROR_ENCODE,
  GST_STREAM_ERROR_DEMUX,
  GST_STREAM_ERROR_MUX,
  GST_STREAM_ERROR_FORMAT,
  GST_STREAM_ERROR_DECRYPT,
  GST_STREAM_ERROR_DECRYPT_NOKEY,
  GST_STREAM_ERROR_NUM_ERRORS
} GstStreamError;

/**
 * GST_LIBRARY_ERROR:
 *
 * Error domain for library loading. Errors in this domain will
 * be from the #GstLibraryError enumeration.
 * See #GError for information on error domains.
 */
#define GST_LIBRARY_ERROR   gst_library_error_quark ()
/**
 * GST_RESOURCE_ERROR:
 *
 * Error domain for resource handling. Errors in this domain will
 * be from the #GstResourceError enumeration.
 * See #GError for information on error domains.
 */
#define GST_RESOURCE_ERROR  gst_resource_error_quark ()
/**
 * GST_CORE_ERROR:
 *
 * Error domain for core system. Errors in this domain will
 * be from the #GstCoreError enumeration.
 * See #GError for information on error domains.
 */
#define GST_CORE_ERROR      gst_core_error_quark ()
/**
 * GST_STREAM_ERROR:
 *
 * Error domain for media stream processing. Errors in this domain will
 * be from the #GstStreamError enumeration.
 * See #GError for information on error domains.
 */
#define GST_STREAM_ERROR    gst_stream_error_quark ()

/**
 * GST_ERROR_SYSTEM: (skip):
 *
 * Builds a string using errno describing the previously failed system
 * call.  To be used as the debug argument in #GST_ELEMENT_ERROR.
 */
#define GST_ERROR_SYSTEM    ("system error: %s", g_strerror (errno))

GST_API
gchar *gst_error_get_message (GQuark domain, gint code);

/**
 * gst_stream_error_quark: (attributes doc.skip=true)
 */
GST_API
GQuark gst_stream_error_quark (void);

/**
 * gst_core_error_quark: (attributes doc.skip=true)
 */
GST_API
GQuark gst_core_error_quark (void);

/**
 * gst_resource_error_quark: (attributes doc.skip=true)
 */
GST_API
GQuark gst_resource_error_quark (void);

/**
 * gst_library_error_quark: (attributes doc.skip=true)
 */
GST_API
GQuark gst_library_error_quark (void);

G_END_DECLS
#endif /* __GST_ERROR_H__ */
