/*  GStreamer video sink base class
 *  Copyright (C) <2003> Julien Moutte <julien@moutte.net>
 *  Copyright (C) <2009> Tim-Philipp Müller <tim centricular net>
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

/* FIXME 0.11: turn this into a proper base class */

#ifndef __GST_VIDEO_SINK_H__
#define __GST_VIDEO_SINK_H__

#include <gst/gst.h>
#include <gst/base/gstbasesink.h>
#include <gst/video/video-prelude.h>
#include <gst/video/video-info.h>

G_BEGIN_DECLS

#define GST_TYPE_VIDEO_SINK (gst_video_sink_get_type())
#define GST_VIDEO_SINK(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_VIDEO_SINK, GstVideoSink))
#define GST_VIDEO_SINK_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_VIDEO_SINK, GstVideoSinkClass))
#define GST_IS_VIDEO_SINK(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_VIDEO_SINK))
#define GST_IS_VIDEO_SINK_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_VIDEO_SINK))
#define GST_VIDEO_SINK_GET_CLASS(klass) \
  (G_TYPE_INSTANCE_GET_CLASS ((klass), GST_TYPE_VIDEO_SINK, GstVideoSinkClass))

/**
 * GST_VIDEO_SINK_CAST:
 * @obj: a #GstVideoSink or derived object
 *
 * Cast @obj to a #GstVideoSink without runtime type check.
 */
#define GST_VIDEO_SINK_CAST(obj)  ((GstVideoSink *) (obj))

/**
 * GST_VIDEO_SINK_PAD:
 * @obj: a #GstVideoSink
 *
 * Get the sink #GstPad of @obj.
 */
#define GST_VIDEO_SINK_PAD(obj) GST_BASE_SINK_PAD(obj)

#define GST_VIDEO_SINK_WIDTH(obj) (GST_VIDEO_SINK_CAST (obj)->width)
#define GST_VIDEO_SINK_HEIGHT(obj) (GST_VIDEO_SINK_CAST (obj)->height)

typedef struct _GstVideoSink GstVideoSink;
typedef struct _GstVideoSinkClass GstVideoSinkClass;
typedef struct _GstVideoRectangle GstVideoRectangle;
typedef struct _GstVideoSinkPrivate GstVideoSinkPrivate;

/**
 * GstVideoRectangle:
 * @x: X coordinate of rectangle's top-left point
 * @y: Y coordinate of rectangle's top-left point
 * @w: width of the rectangle
 * @h: height of the rectangle
 *
 * Helper structure representing a rectangular area.
 */
struct _GstVideoRectangle {
  gint x;
  gint y;
  gint w;
  gint h;
};

/**
 * GstVideoSink:
 * @height: video height (derived class needs to set this)
 * @width: video width (derived class needs to set this)
 *
 * The video sink instance structure. Derived video sinks should set the
 * @height and @width members.
 */
struct _GstVideoSink {
  GstBaseSink element;    /* FIXME 0.11: this should not be called 'element' */

  /*< public >*/
  gint width, height;

  /*< private >*/
  GstVideoSinkPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstVideoSinkClass:
 * @parent_class: the parent class structure
 * @show_frame: render a video frame. Maps to #GstBaseSinkClass.render() and
 *     #GstBaseSinkClass.preroll() vfuncs. Rendering during preroll will be
 *     suppressed if the #GstVideoSink:show-preroll-frame property is set to
 *     %FALSE.
 *
 * The video sink class structure. Derived classes should override the
 * @show_frame virtual function.
 */
struct _GstVideoSinkClass {
  GstBaseSinkClass parent_class;

  GstFlowReturn  (*show_frame) (GstVideoSink *video_sink, GstBuffer *buf);

  /**
   * GstVideoSinkClass::set_info:
   * @caps: A #GstCaps.
   * @info: A #GstVideoInfo corresponding to @caps.
   *
   * Notifies the subclass of changed #GstVideoInfo.
   *
   * Since: 1.20
   */
  gboolean       (*set_info)   (GstVideoSink *video_sink, GstCaps *caps, const GstVideoInfo *info);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING-1];
};

GST_VIDEO_API
GType gst_video_sink_get_type (void);

GST_VIDEO_DEPRECATED_FOR(gst_video_center_rect)
void gst_video_sink_center_rect (GstVideoRectangle src, GstVideoRectangle dst,
                                 GstVideoRectangle *result, gboolean scaling);

GST_VIDEO_API
void gst_video_center_rect      (const GstVideoRectangle * src,
                                 const GstVideoRectangle * dst,
                                 GstVideoRectangle * result,
                                 gboolean scaling);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstVideoSink, gst_object_unref)

G_END_DECLS

#endif  /* __GST_VIDEO_SINK_H__ */
