/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *               2000,2004 Wim Taymans <wim@fluendo.com>
 *
 * gstelementmetadata.h: Metadata for GstElement classes
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

#ifndef __GST_ELEMENT_METADATA_H__
#define __GST_ELEMENT_METADATA_H__

#include <glib.h>

G_BEGIN_DECLS

/**
 * GST_ELEMENT_METADATA_LONGNAME:
 *
 * The long English name of the element. E.g. "File Sink"
 */
#define GST_ELEMENT_METADATA_LONGNAME      "long-name"
/**
 * GST_ELEMENT_METADATA_KLASS:
 *
 * String describing the type of element, as an unordered list
 * separated with slashes ('/'). See draft-klass.txt of the design docs
 * for more details and common types. E.g: "Sink/File"
 */
#define GST_ELEMENT_METADATA_KLASS         "klass"

/**
 * GST_ELEMENT_METADATA_DESCRIPTION:
 *
 * Sentence describing the purpose of the element.
 * E.g: "Write stream to a file"
 */
#define GST_ELEMENT_METADATA_DESCRIPTION   "description"
/**
 * GST_ELEMENT_METADATA_AUTHOR:
 *
 * Name and contact details of the author(s). Use \n to separate
 * multiple author details.
 * E.g: "Joe Bloggs &lt;joe.blogs at foo.com&gt;"
 */
#define GST_ELEMENT_METADATA_AUTHOR        "author"

/**
 * GST_ELEMENT_METADATA_DOC_URI:
 *
 * Set uri pointing to user documentation. Applications can use this to show
 * help for e.g. effects to users.
 */
#define GST_ELEMENT_METADATA_DOC_URI       "doc-uri"
/**
 * GST_ELEMENT_METADATA_ICON_NAME:
 *
 * Elements that bridge to certain other products can include an icon of that
 * used product. Application can show the icon in menus/selectors to help
 * identifying specific elements.
 */
#define GST_ELEMENT_METADATA_ICON_NAME     "icon-name"

G_END_DECLS

#endif /* __GST_ELEMENT_METADATA_H__ */
