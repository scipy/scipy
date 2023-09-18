/* GStreamer
 * Copyright (C) 2010, 2012 Alexander Saprykin <xelfium@gmail.com>
 *
 * gsttocsetter.h: Interfaces for TOC
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

#ifndef __GST_TOC_SETTER_H__
#define __GST_TOC_SETTER_H__

#include <gst/gst.h>

G_BEGIN_DECLS

#define GST_TYPE_TOC_SETTER              (gst_toc_setter_get_type ())
#define GST_TOC_SETTER(obj)              (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_TOC_SETTER, GstTocSetter))
#define GST_IS_TOC_SETTER(obj)           (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_TOC_SETTER))
#define GST_TOC_SETTER_GET_IFACE(obj)    (G_TYPE_INSTANCE_GET_INTERFACE ((obj), GST_TYPE_TOC_SETTER, GstTocSetterInterface))
/**
 * GstTocSetter:
 *
 * Opaque #GstTocSetter data structure.
 */
typedef struct _GstTocSetter GstTocSetter;
typedef struct _GstTocSetterInterface GstTocSetterInterface;

/**
 * GstTocSetterInterface:
 * @g_iface: parent interface type.
 *
 * #GstTocSetterInterface interface.
 */

struct _GstTocSetterInterface
{
  GTypeInterface g_iface;

  /* signals */

  /* virtual table */
};

GST_API
GType         gst_toc_setter_get_type (void);

GST_API
void          gst_toc_setter_reset   (GstTocSetter *setter);

GST_API
GstToc *      gst_toc_setter_get_toc (GstTocSetter *setter);

GST_API
void          gst_toc_setter_set_toc (GstTocSetter *setter, GstToc *toc);

G_END_DECLS

#endif /* __GST_TOC_SETTER_H__ */

