/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2000 Wim Taymans <wtay@chello.be>
 *
 * gstbin.h: Header for GstBin container object
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


#ifndef __GST_BIN_H__
#define __GST_BIN_H__

#include <gst/gstelement.h>
#include <gst/gstiterator.h>
#include <gst/gstbus.h>

G_BEGIN_DECLS

#define GST_TYPE_BIN             (gst_bin_get_type ())
#define GST_IS_BIN(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_BIN))
#define GST_IS_BIN_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_BIN))
#define GST_BIN_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_BIN, GstBinClass))
#define GST_BIN(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_BIN, GstBin))
#define GST_BIN_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_BIN, GstBinClass))
#define GST_BIN_CAST(obj)        ((GstBin*)(obj))

/**
 * GstBinFlags:
 * @GST_BIN_FLAG_LAST: the last enum in the series of flags for bins.
 * Derived classes can use this as first value in a list of flags.
 *
 * GstBinFlags are a set of flags specific to bins. Most are set/used
 * internally. They can be checked using the GST_OBJECT_FLAG_IS_SET() macro,
 * and (un)set using GST_OBJECT_FLAG_SET() and GST_OBJECT_FLAG_UNSET().
 */
typedef enum {
  /**
   * GST_BIN_FLAG_NO_RESYNC:
   *
   * Don't resync a state change when elements are added or linked in the bin
   *
   * Since: 1.0.5
   */
  GST_BIN_FLAG_NO_RESYNC	= (GST_ELEMENT_FLAG_LAST << 0),

  /**
   * GST_BIN_FLAG_STREAMS_AWARE:
   *
   * Indicates whether the bin can handle elements that add/remove source pads
   * at any point in time without first posting a no-more-pads signal.
   *
   * Since: 1.10
   */
  GST_BIN_FLAG_STREAMS_AWARE	= (GST_ELEMENT_FLAG_LAST << 1),

  /* padding */

  /**
   * GST_BIN_FLAG_LAST:
   *
   * The last enum in the series of flags for bins. Derived classes can use this
   * as first value in a list of flags.
   */
  GST_BIN_FLAG_LAST		= (GST_ELEMENT_FLAG_LAST << 5)
} GstBinFlags;

/**
 * GST_BIN_IS_NO_RESYNC:
 * @bin: A #GstBin
 *
 * Check if @bin will resync its state change when elements are added and
 * removed.
 *
 * Since: 1.0.5
 */
#define GST_BIN_IS_NO_RESYNC(bin)        (GST_OBJECT_FLAG_IS_SET(bin,GST_BIN_FLAG_NO_RESYNC))

typedef struct _GstBin GstBin;
typedef struct _GstBinClass GstBinClass;
typedef struct _GstBinPrivate GstBinPrivate;

/**
 * GST_BIN_NUMCHILDREN:
 * @bin: a #GstBin
 *
 * Gets the number of children in a bin.
 */
#define GST_BIN_NUMCHILDREN(bin)	(GST_BIN_CAST(bin)->numchildren)
/**
 * GST_BIN_CHILDREN:
 * @bin: a #GstBin
 *
 * Gets the list of children in a bin.
 */
#define GST_BIN_CHILDREN(bin)		(GST_BIN_CAST(bin)->children)
/**
 * GST_BIN_CHILDREN_COOKIE:
 * @bin: a #GstBin
 *
 * Gets the children cookie that watches the children list.
 */
#define GST_BIN_CHILDREN_COOKIE(bin)	(GST_BIN_CAST(bin)->children_cookie)

/**
 * GstBin:
 * @numchildren: the number of children in this bin
 * @children: (element-type Gst.Element): the list of children in this bin
 * @children_cookie: updated whenever @children changes
 * @child_bus: internal bus for handling child messages
 * @messages: (element-type Gst.Message): queued and cached messages
 * @polling: the bin is currently calculating its state
 * @state_dirty: the bin needs to recalculate its state (deprecated)
 * @clock_dirty: the bin needs to select a new clock
 * @provided_clock: the last clock selected
 * @clock_provider: the element that provided @provided_clock
 *
 * The GstBin base class. Subclasses can access these fields provided
 * the LOCK is taken.
 */
struct _GstBin {
  GstElement	 element;

  /*< public >*/ /* with LOCK */
  /* our children, subclass are supposed to update these
   * fields to reflect their state with _iterate_*() */
  gint		 numchildren;
  GList		*children;
  guint32	 children_cookie;

  GstBus        *child_bus;
  GList         *messages;

  gboolean	 polling;
  gboolean       state_dirty;

  gboolean       clock_dirty;
  GstClock	*provided_clock;
  GstElement    *clock_provider;

  /*< private >*/
  GstBinPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstBinClass:
 * @parent_class: bin parent class
 *
 * Subclasses can override #GstBinClass::add_element and #GstBinClass::remove_element
 * to update the list of children in the bin.
 *
 * The #GstBinClass::handle_message method can be overridden to implement custom
 * message handling.
 *
 * #GstBinClass::deep_element_added will be called when a new element has been
 * added to any bin inside this bin, so it will also be called if a new child
 * was added to a sub-bin of this bin. #GstBin implementations that override
 * this message should chain up to the parent class implementation so the
 * #GstBin::deep-element-added signal is emitted on all parents.
 */
struct _GstBinClass {
  GstElementClass parent_class;

  /*< private >*/
  GThreadPool  *pool; /* deprecated */

  /* signals */

  /**
   * GstBinClass::element_added:
   * @bin: the #GstBin
   * @child: the element that was added
   *
   * Method called when an element was added to the bin.
   */
  void		(*element_added)	(GstBin *bin, GstElement *child);

  /**
   * GstBinClass::element_removed:
   * @bin: the #GstBin
   * @child: the element that was removed
   *
   * Method called when an element was removed from the bin.
   */
  void		(*element_removed)	(GstBin *bin, GstElement *child);

  /*< public >*/
  /* virtual methods for subclasses */

  /**
   * GstBinClass::add_element:
   * @bin: the #GstBin
   * @element: the element to be added
   *
   * Method to add an element to the bin.
   *
   * Returns: %TRUE if the @element was added
   */
  gboolean	(*add_element)		(GstBin *bin, GstElement *element);

  /**
   * GstBinClass::remove_element:
   * @bin: the #GstBin
   * @element: the element to be removed
   *
   * Method to remove an element from the bin.
   *
   * Returns: %TRUE if the @element was removed
   */
  gboolean	(*remove_element)	(GstBin *bin, GstElement *element);

  /**
   * GstBinClass::handle_message:
   * @bin: the #GstBin
   * @message: (transfer full): the message to be handled
   *
   * Method to handle a message from the children.
   */
  void		(*handle_message)	(GstBin *bin, GstMessage *message);

  /*< private >*/
  /* signal */
  gboolean	(*do_latency)           (GstBin *bin);

  /*< public >*/
  /* signal */

  /**
   * GstBinClass::deep_element_added:
   * @bin: the top level #GstBin
   * @sub_bin: the #GstBin to which the element was added
   * @child: the element that was added
   *
   * Method called when an element was added somewhere in the bin hierarchy.
   */
  void          (*deep_element_added)   (GstBin *bin, GstBin *sub_bin, GstElement *child);

  /**
   * GstBinClass::deep_element_removed:
   * @bin: the top level #GstBin
   * @sub_bin: the #GstBin from which the element was removed
   * @child: the element that was removed
   *
   * Method called when an element was removed somewhere in the bin hierarchy.
   */
  void          (*deep_element_removed) (GstBin *bin, GstBin *sub_bin, GstElement *child);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING-2];
};

GST_API
GType		gst_bin_get_type		(void);

GST_API
GstElement*	gst_bin_new			(const gchar *name);

/* add and remove elements from the bin */

GST_API
gboolean	gst_bin_add			(GstBin *bin, GstElement *element);

GST_API
gboolean	gst_bin_remove			(GstBin *bin, GstElement *element);

/* retrieve a single child */

GST_API
GstElement*	gst_bin_get_by_name		 (GstBin *bin, const gchar *name);

GST_API
GstElement*	gst_bin_get_by_name_recurse_up	 (GstBin *bin, const gchar *name);

GST_API
GstElement*	gst_bin_get_by_interface	 (GstBin *bin, GType iface);

/* retrieve multiple children */

GST_API
GstIterator*    gst_bin_iterate_elements	 (GstBin *bin);

GST_API
GstIterator*    gst_bin_iterate_sorted		 (GstBin *bin);

GST_API
GstIterator*    gst_bin_iterate_recurse		 (GstBin *bin);

GST_API
GstIterator*	gst_bin_iterate_sinks		 (GstBin *bin);

GST_API
GstIterator*	gst_bin_iterate_sources		 (GstBin *bin);

GST_API
GstIterator*	gst_bin_iterate_all_by_interface (GstBin *bin, GType iface);

GST_API
GstIterator*	gst_bin_iterate_all_by_element_factory_name (GstBin * bin, const gchar * factory_name);

/* latency */

GST_API
gboolean        gst_bin_recalculate_latency      (GstBin * bin);

/* set and get suppressed flags */

GST_API
void            gst_bin_set_suppressed_flags (GstBin * bin, GstElementFlags flags);

GST_API
GstElementFlags gst_bin_get_suppressed_flags (GstBin * bin);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstBin, gst_object_unref)

G_END_DECLS


#endif /* __GST_BIN_H__ */
