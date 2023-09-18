/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2000 Wim Taymans <wtay@chello.be>
 *
 * gst.h: Main header for GStreamer, apps should include this
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


#ifndef __GST_H__
#define __GST_H__

#include <glib.h>

#include <gst/glib-compat.h>

#include <gst/gstenumtypes.h>
#include <gst/gstversion.h>

#include <gst/gstatomicqueue.h>
#include <gst/gstbin.h>
#include <gst/gstbuffer.h>
#include <gst/gstbufferlist.h>
#include <gst/gstbufferpool.h>
#include <gst/gstcaps.h>
#include <gst/gstcapsfeatures.h>
#include <gst/gstchildproxy.h>
#include <gst/gstclock.h>
#include <gst/gstcontrolsource.h>
#include <gst/gstdatetime.h>
#include <gst/gstdebugutils.h>
#include <gst/gstdevice.h>
#include <gst/gstdevicemonitor.h>
#include <gst/gstdeviceprovider.h>
#include <gst/gstdynamictypefactory.h>
#include <gst/gstelement.h>
#include <gst/gstelementmetadata.h>
#include <gst/gsterror.h>
#include <gst/gstevent.h>
#include <gst/gstghostpad.h>
#include <gst/gstinfo.h>
#include <gst/gstiterator.h>
#include <gst/gstmessage.h>
#include <gst/gstmemory.h>
#include <gst/gstmeta.h>
#include <gst/gstminiobject.h>
#include <gst/gstobject.h>
#include <gst/gststreamcollection.h>
#include <gst/gstpad.h>
#include <gst/gstparamspecs.h>
#include <gst/gstpipeline.h>
#include <gst/gstplugin.h>
#include <gst/gstpoll.h>
#include <gst/gstpreset.h>
#include <gst/gstprotection.h>
#include <gst/gstquery.h>
#include <gst/gstregistry.h>
#include <gst/gstpromise.h>
#include <gst/gstsample.h>
#include <gst/gstsegment.h>
#include <gst/gststreams.h>
#include <gst/gststructure.h>
#include <gst/gstsystemclock.h>
#include <gst/gsttaglist.h>
#include <gst/gsttagsetter.h>
#include <gst/gsttask.h>
#include <gst/gsttaskpool.h>
#include <gst/gsttoc.h>
#include <gst/gsttocsetter.h>
#include <gst/gsttracer.h>
#include <gst/gsttracerfactory.h>
#include <gst/gsttracerrecord.h>
#include <gst/gsttypefind.h>
#include <gst/gsttypefindfactory.h>
#include <gst/gsturi.h>
#include <gst/gstutils.h>
#include <gst/gstvalue.h>

#include <gst/gstparse.h>

/* API compatibility stuff */
#include <gst/gstcompat.h>

#ifdef __APPLE__
#	include <TargetConditionals.h>
#	if TARGET_OS_MAC && !TARGET_OS_IPHONE
#	 include <gst/gstmacos.h>
#	endif
#endif

G_BEGIN_DECLS

GST_API
void		gst_init			(int *argc, char **argv[]);

GST_API
gboolean	gst_init_check			(int *argc, char **argv[],
						 GError ** error);
GST_API
gboolean        gst_is_initialized              (void);

GST_API
GOptionGroup *	gst_init_get_option_group	(void);

GST_API
void		gst_deinit			(void);

GST_API
void		gst_version			(guint *major, guint *minor,
						 guint *micro, guint *nano);
GST_API
gchar *		gst_version_string		(void);

GST_API
gboolean        gst_segtrap_is_enabled          (void);

GST_API
void            gst_segtrap_set_enabled         (gboolean enabled);

GST_API
gboolean        gst_registry_fork_is_enabled    (void);

GST_API
void            gst_registry_fork_set_enabled   (gboolean enabled);

GST_API
gboolean        gst_update_registry             (void);

GST_API
const gchar *   gst_get_main_executable_path    (void);

G_END_DECLS

#endif /* __GST_H__ */
