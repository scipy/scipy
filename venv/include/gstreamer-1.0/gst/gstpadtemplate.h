/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2000 Wim Taymans <wim.taymans@chello.be>
 *
 * gstpadtemplate.h: Header for GstPadTemplate object
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


#ifndef __GST_PAD_TEMPLATE_H__
#define __GST_PAD_TEMPLATE_H__

#include <gst/gstconfig.h>

typedef struct _GstPadTemplate GstPadTemplate;
typedef struct _GstPadTemplateClass GstPadTemplateClass;
typedef struct _GstStaticPadTemplate GstStaticPadTemplate;

#include <gst/gstobject.h>
#include <gst/gstbuffer.h>
#include <gst/gstcaps.h>
#include <gst/gstevent.h>
#include <gst/gstquery.h>
#include <gst/gsttask.h>

G_BEGIN_DECLS

#define GST_TYPE_STATIC_PAD_TEMPLATE	(gst_static_pad_template_get_type ())

#define GST_TYPE_PAD_TEMPLATE		(gst_pad_template_get_type ())
#define GST_PAD_TEMPLATE(obj)		(G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_PAD_TEMPLATE,GstPadTemplate))
#define GST_PAD_TEMPLATE_CLASS(klass)	(G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_PAD_TEMPLATE,GstPadTemplateClass))
#define GST_IS_PAD_TEMPLATE(obj)	(G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_PAD_TEMPLATE))
#define GST_IS_PAD_TEMPLATE_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_PAD_TEMPLATE))

/**
 * GstPadPresence:
 * @GST_PAD_ALWAYS: the pad is always available
 * @GST_PAD_SOMETIMES: the pad will become available depending on the media stream
 * @GST_PAD_REQUEST: the pad is only available on request with
 *  gst_element_request_pad().
 *
 * Indicates when this pad will become available.
 */
typedef enum {
  GST_PAD_ALWAYS,
  GST_PAD_SOMETIMES,
  GST_PAD_REQUEST
} GstPadPresence;

/**
 * GST_PAD_TEMPLATE_NAME_TEMPLATE:
 * @templ: the template to query
 *
 * Get the nametemplate of the padtemplate.
 */
#define GST_PAD_TEMPLATE_NAME_TEMPLATE(templ)	(((GstPadTemplate *)(templ))->name_template)

/**
 * GST_PAD_TEMPLATE_DIRECTION:
 * @templ: the template to query
 *
 * Get the #GstPadDirection of the padtemplate.
 */
#define GST_PAD_TEMPLATE_DIRECTION(templ)	(((GstPadTemplate *)(templ))->direction)

/**
 * GST_PAD_TEMPLATE_PRESENCE:
 * @templ: the template to query
 *
 * Get the #GstPadPresence of the padtemplate.
 */
#define GST_PAD_TEMPLATE_PRESENCE(templ)	(((GstPadTemplate *)(templ))->presence)

/**
 * GST_PAD_TEMPLATE_CAPS:
 * @templ: the template to query
 *
 * Get a handle to the padtemplate #GstCaps
 */
#define GST_PAD_TEMPLATE_CAPS(templ)		(((GstPadTemplate *)(templ))->caps)

/**
 * GST_PAD_TEMPLATE_GTYPE:
 * @templ: the template to query
 *
 * Get the #GType of the padtemplate
 *
 * Since: 1.14
 */
#define GST_PAD_TEMPLATE_GTYPE(templ)		(((GstPadTemplate *)(templ))->ABI.abi.gtype)

/**
 * GstPadTemplateFlags:
 * @GST_PAD_TEMPLATE_FLAG_LAST: first flag that can be used by subclasses.
 *
 * Flags for the padtemplate
 */
typedef enum {
  /* padding */
  GST_PAD_TEMPLATE_FLAG_LAST    = (GST_OBJECT_FLAG_LAST << 4)
} GstPadTemplateFlags;

/**
 * GST_PAD_TEMPLATE_IS_FIXED:
 * @templ: the template to query
 *
 * Check if the properties of the padtemplate are fixed
 */
#define GST_PAD_TEMPLATE_IS_FIXED(templ)	(GST_OBJECT_FLAG_IS_SET(templ, GST_PAD_TEMPLATE_FIXED))

/**
 * GstPadTemplate:
 *
 * The padtemplate object.
 */
struct _GstPadTemplate {
  GstObject	   object;

  gchar           *name_template;
  GstPadDirection  direction;
  GstPadPresence   presence;
  GstCaps	  *caps;

  /*< private >*/
  union {
    gpointer _gst_reserved[GST_PADDING];
    struct {
      GType gtype;
      GstCaps *documentation_caps;
    } abi;
  } ABI;
};

struct _GstPadTemplateClass {
  GstObjectClass   parent_class;

  /* signal callbacks */
  void (*pad_created)	(GstPadTemplate *templ, GstPad *pad);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstStaticPadTemplate:
 * @name_template: the name of the template
 * @direction: the direction of the template
 * @presence: the presence of the template
 * @static_caps: the caps of the template.
 *
 * Structure describing the #GstStaticPadTemplate.
 */
struct _GstStaticPadTemplate {
  const gchar     *name_template;
  GstPadDirection  direction;
  GstPadPresence   presence;
  GstStaticCaps    static_caps;
};

/**
 * GST_STATIC_PAD_TEMPLATE:
 * @padname: the name template of the pad
 * @dir: the GstPadDirection of the pad
 * @pres: the GstPadPresence of the pad
 * @caps: the GstStaticCaps of the pad
 *
 * Convenience macro to fill the values of a #GstStaticPadTemplate
 * structure.
 * Example:
 * |[<!-- language="C" -->
 * static GstStaticPadTemplate my_src_template = \
 *   GST_STATIC_PAD_TEMPLATE("src", GST_PAD_SRC, GST_PAD_ALWAYS,
 *                           GST_STATIC_CAPS_ANY);
 * ]|
 */
#define GST_STATIC_PAD_TEMPLATE(padname, dir, pres, caps) \
{ \
  /* name_template */    padname, \
  /* direction */        dir, \
  /* presence */         pres, \
  /* caps */             caps \
}

/* templates and factories */

GST_API
GType			gst_pad_template_get_type		(void);

/**
 * gst_static_pad_template_get_type: (attributes doc.skip=true)
 */
GST_API
GType			gst_static_pad_template_get_type	(void);

GST_API
GstPadTemplate*		gst_pad_template_new			(const gchar *name_template,
								 GstPadDirection direction, GstPadPresence presence,
								 GstCaps *caps) G_GNUC_MALLOC;
GST_API
GstPadTemplate*		gst_pad_template_new_with_gtype		(const gchar *name_template,
								 GstPadDirection direction, GstPadPresence presence,
								 GstCaps *caps, GType pad_type) G_GNUC_MALLOC;
GST_API
GstPadTemplate *	gst_static_pad_template_get             (GstStaticPadTemplate *pad_template);

GST_API
GstPadTemplate * gst_pad_template_new_from_static_pad_template_with_gtype (
    GstStaticPadTemplate * pad_template,
    GType pad_type);

GST_API
GstCaps*		gst_static_pad_template_get_caps	(GstStaticPadTemplate *templ);

GST_API
GstCaps*		gst_pad_template_get_caps		(GstPadTemplate *templ);

GST_API
void        gst_pad_template_set_documentation_caps (GstPadTemplate *templ, GstCaps *caps);

GST_API
GstCaps*    gst_pad_template_get_documentation_caps (GstPadTemplate *templ);

GST_API
void                    gst_pad_template_pad_created            (GstPadTemplate * templ, GstPad * pad);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstPadTemplate, gst_object_unref)

G_END_DECLS

#endif /* __GST_PAD_TEMPLATE_H__ */

