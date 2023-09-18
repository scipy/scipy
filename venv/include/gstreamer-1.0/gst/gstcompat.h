/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2004 Wim Taymans <wim@fluendo.com>
 *
 * gstcompat.h: backwards compatibility stuff
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

/* API compatibility stuff */
#ifndef __GSTCOMPAT_H__
#define __GSTCOMPAT_H__
/**
 * SECTION:gstcompat
 * @title: GstCompat
 * @short_description: Deprecated API entries
 *
 * Please do not use these in new code.
 * These symbols are only available by defining GST_DISABLE_DEPRECATED.
 * This can be done in CFLAGS for compiling old code.
 */

#include <glib.h>
#include <gst/gstpad.h>


G_BEGIN_DECLS

#define gst_buffer_new_and_alloc(s)            gst_buffer_new_allocate(NULL, s, NULL)

#define GST_BUFFER_TIMESTAMP          GST_BUFFER_PTS
#define GST_BUFFER_TIMESTAMP_IS_VALID GST_BUFFER_PTS_IS_VALID

static inline gboolean
gst_pad_set_caps (GstPad * pad, GstCaps * caps)
{
  GstEvent *event;
  gboolean res = TRUE;

  g_return_val_if_fail (GST_IS_PAD (pad), FALSE);
  g_return_val_if_fail (caps != NULL && gst_caps_is_fixed (caps), FALSE);

  event = gst_event_new_caps (caps);

  if (GST_PAD_IS_SRC (pad))
    res = gst_pad_push_event (pad, event);
  else
    res = gst_pad_send_event (pad, event);

  return res;
}


#ifndef GST_DISABLE_DEPRECATED

/* added to ease the transition to 0.11 */
#define gst_element_class_set_details_simple  gst_element_class_set_metadata

#define gst_element_factory_get_longname(f)    gst_element_factory_get_metadata(f, GST_ELEMENT_METADATA_LONGNAME)
#define gst_element_factory_get_klass(f)       gst_element_factory_get_metadata(f, GST_ELEMENT_METADATA_KLASS)
#define gst_element_factory_get_description(f) gst_element_factory_get_metadata(f, GST_ELEMENT_METADATA_DESCRIPTION)
#define gst_element_factory_get_author(f)      gst_element_factory_get_metadata(f, GST_ELEMENT_METADATA_AUTHOR)
#define gst_element_factory_get_documentation_uri(f)  gst_element_factory_get_metadata(f, GST_ELEMENT_METADATA_DOC_URI)
#define gst_element_factory_get_icon_name(f)   gst_element_factory_get_metadata(f, GST_ELEMENT_METADATA_ICON_NAME)

#define gst_pad_get_caps_reffed(p)             gst_pad_get_caps(p)
#define gst_pad_peer_get_caps_reffed(p)        gst_pad_peer_get_caps(p)

#define gst_adapter_prev_timestamp    gst_adapter_prev_pts

#define gst_tag_list_free(taglist)    gst_tag_list_unref(taglist)

#define GST_MESSAGE_DURATION GST_MESSAGE_DURATION_CHANGED
#define gst_message_new_duration(src,fmt,dur) \
    gst_message_new_duration_changed(src)
#define gst_message_parse_duration(msg,fmt,dur) \
G_STMT_START { \
  GstFormat *p_fmt = fmt; \
  gint64 *p_dur = dur; \
  if (p_fmt) \
    *p_fmt = GST_FORMAT_TIME; \
  if (p_dur) \
    *p_dur = GST_CLOCK_TIME_NONE; \
} G_STMT_END

#endif /* not GST_DISABLE_DEPRECATED */

G_END_DECLS

#endif /* __GSTCOMPAT_H__ */
