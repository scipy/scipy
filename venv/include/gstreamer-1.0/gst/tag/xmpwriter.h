/* GStreamer XmpConfig
 * Copyright (C) 2011 Thiago Santos <thiago.sousa.santos@collabora.co.uk>
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

#ifndef __TAG_XMP_WRITER_H__
#define __TAG_XMP_WRITER_H__

#include <gst/gst.h>
#include <gst/tag/tag-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_TAG_XMP_WRITER \
  (gst_tag_xmp_writer_get_type ())
#define GST_TAG_XMP_WRITER(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_TAG_XMP_WRITER, GstTagXmpWriter))
#define GST_IS_TAG_XMP_WRITER(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_TAG_XMP_WRITER))
#define GST_TAG_XMP_WRITER_GET_INTERFACE(inst) \
  (G_TYPE_INSTANCE_GET_INTERFACE ((inst), GST_TYPE_TAG_XMP_WRITER, GstTagXmpWriterInterface))

typedef struct _GstTagXmpWriter GstTagXmpWriter;
typedef struct _GstTagXmpWriterInterface GstTagXmpWriterInterface;

struct _GstTagXmpWriterInterface {
  GTypeInterface parent;
};

GST_TAG_API
GType           gst_tag_xmp_writer_get_type		(void);

GST_TAG_API
void		gst_tag_xmp_writer_add_all_schemas	(GstTagXmpWriter * config);

GST_TAG_API
void		gst_tag_xmp_writer_add_schema	(GstTagXmpWriter * config,
						const gchar * schema);

GST_TAG_API
gboolean	gst_tag_xmp_writer_has_schema	(GstTagXmpWriter * config,
						const gchar * schema);

GST_TAG_API
void		gst_tag_xmp_writer_remove_schema	(GstTagXmpWriter * config,
						const gchar * schema);

GST_TAG_API
void		gst_tag_xmp_writer_remove_all_schemas (GstTagXmpWriter * config);

GST_TAG_API
GstBuffer*	gst_tag_xmp_writer_tag_list_to_xmp_buffer 	(GstTagXmpWriter * config,
							 const GstTagList * taglist,
							 gboolean read_only);

G_END_DECLS

#endif /* __TAG_XMP_WRITER_H__ */
