/* GStreamer Base Class for Tag Demuxing
 * Copyright (C) 2005  Jan Schmidt <thaytan@mad.scientist.com>
 * Copyright (C) 2006  Tim-Philipp Müller <tim centricular net>
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

#ifndef __GST_TAG_DEMUX_H__
#define __GST_TAG_DEMUX_H__

#include <gst/gst.h>
#include <gst/tag/tag-enumtypes.h>

G_BEGIN_DECLS

#define GST_TYPE_TAG_DEMUX            (gst_tag_demux_get_type())
#define GST_TAG_DEMUX(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_TAG_DEMUX,GstTagDemux))
#define GST_TAG_DEMUX_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_TAG_DEMUX,GstTagDemuxClass))
#define GST_IS_TAG_DEMUX(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_TAG_DEMUX))
#define GST_IS_TAG_DEMUX_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_TAG_DEMUX))

typedef struct _GstTagDemux        GstTagDemux;
typedef struct _GstTagDemuxClass   GstTagDemuxClass;
typedef struct _GstTagDemuxPrivate GstTagDemuxPrivate;

/**
 * GstTagDemuxResult:
 * @GST_TAG_DEMUX_RESULT_BROKEN_TAG: cannot parse tag, just skip it
 * @GST_TAG_DEMUX_RESULT_AGAIN     : call again with less or more data
 * @GST_TAG_DEMUX_RESULT_OK	   : parsed tag successfully
 *
 * Result values from the parse_tag virtual function.
 */
typedef enum {
  GST_TAG_DEMUX_RESULT_BROKEN_TAG,
  GST_TAG_DEMUX_RESULT_AGAIN,
  GST_TAG_DEMUX_RESULT_OK
} GstTagDemuxResult;

/**
 * GstTagDemux:
 * @element: parent element
 *
 * Opaque #GstTagDemux structure.
 */
struct _GstTagDemux
{
  GstElement element;

  /*< private >*/
  GstTagDemuxPrivate  *priv;

  gpointer             reserved[GST_PADDING];
};

/**
 * GstTagDemuxClass:
 * @parent_class: the parent class.
 * @min_start_size: minimum size required to identify a tag at the start and
 * determine its total size. Set to 0 if not interested in start tags.
 * Subclasses should set this in their class_init function.
 * @min_end_size: minimum size required to identify a tag at the end and
 * determine its total size. Set to 0 if not interested in end tags.
 * Subclasses should set this in their class_init function.
 * @identify_tag: identify tag and determine the size required to parse the
 * tag. Buffer may be larger than the specified minimum size.
 * Subclassed MUST override this vfunc in their class_init function.
 * @parse_tag: parse the tag. Buffer will be exactly of the size determined by
 * the identify_tag vfunc before. The parse_tag vfunc may change the size
 * stored in *tag_size and return GST_TAG_DEMUX_RESULT_AGAIN to request a
 * larger or smaller buffer. It is also permitted to adjust the tag_size to a
 * smaller value and then return GST_TAG_DEMUX_RESULT_OK in one go.
 * Subclassed MUST override the parse_tag vfunc in their class_init function.
 * @merge_tags: merge start and end tags. Subclasses may want to override this
 * vfunc to allow prioritising of start or end tag according to user
 * preference.  Note that both start_tags and end_tags may be NULL. By default
 * start tags are preferred over end tags.
 *
 * The #GstTagDemuxClass structure.  See documentation at beginning of section
 * for details about what subclasses need to override and do.
 */
struct _GstTagDemuxClass
{
  GstElementClass parent_class;

  /* minimum size required to identify a tag at the start and determine
   * its total size */
  guint                  min_start_size;

  /* minimum size required to identify a tag at the end and determine
   * its total size */
  guint                  min_end_size;

  /* vtable */

  /* identify tag and determine the size required to parse the tag */
  gboolean               (*identify_tag)  (GstTagDemux * demux,
                                           GstBuffer   * buffer,
                                           gboolean      start_tag,
                                           guint       * tag_size);

  /* parse the tag once it is identified and its size is known */
  GstTagDemuxResult      (*parse_tag)     (GstTagDemux * demux,
                                           GstBuffer   * buffer,
                                           gboolean      start_tag,
                                           guint       * tag_size,
                                           GstTagList ** tags);

  /* merge start and end tags (optional) */
  GstTagList *           (*merge_tags)    (GstTagDemux      * demux,
                                           const GstTagList * start_tags,
                                           const GstTagList * end_tags);

  /*< private >*/
  gpointer               reserved[GST_PADDING];
};

GST_TAG_API
GType     gst_tag_demux_get_type (void);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstTagDemux, gst_object_unref)

G_END_DECLS

#endif /* __GST_TAG_DEMUX_H__ */

