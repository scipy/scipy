/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *               2000,2004 Wim Taymans <wim@fluendo.com>
 *
 * gstelementfactory.h: Header for GstElementFactory
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


#ifndef __GST_ELEMENT_FACTORY_H__
#define __GST_ELEMENT_FACTORY_H__

/**
 * GstElementFactory:
 *
 * The opaque #GstElementFactory data structure.
 */
typedef struct _GstElementFactory GstElementFactory;
typedef struct _GstElementFactoryClass GstElementFactoryClass;

#include <gst/gstconfig.h>
#include <gst/gstelement.h>
#include <gst/gstpad.h>
#include <gst/gstplugin.h>
#include <gst/gstpluginfeature.h>
#include <gst/gsturi.h>

G_BEGIN_DECLS

#define GST_TYPE_ELEMENT_FACTORY                (gst_element_factory_get_type())
#define GST_ELEMENT_FACTORY(obj)                (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_ELEMENT_FACTORY,\
                                                 GstElementFactory))
#define GST_ELEMENT_FACTORY_CLASS(klass)        (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_ELEMENT_FACTORY,\
                                                 GstElementFactoryClass))
#define GST_IS_ELEMENT_FACTORY(obj)             (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_ELEMENT_FACTORY))
#define GST_IS_ELEMENT_FACTORY_CLASS(klass)     (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_ELEMENT_FACTORY))
#define GST_ELEMENT_FACTORY_CAST(obj)           ((GstElementFactory *)(obj))

GST_API
GType                   gst_element_factory_get_type            (void);

GST_API
GstElementFactory *     gst_element_factory_find                (const gchar *name);

GST_API
GType                   gst_element_factory_get_element_type    (GstElementFactory *factory);

GST_API
const gchar *           gst_element_factory_get_metadata        (GstElementFactory *factory, const gchar *key);

GST_API
gchar **                gst_element_factory_get_metadata_keys   (GstElementFactory *factory);

GST_API
guint                   gst_element_factory_get_num_pad_templates (GstElementFactory *factory);

GST_API
const GList *           gst_element_factory_get_static_pad_templates (GstElementFactory *factory);

GST_API
GstURIType              gst_element_factory_get_uri_type        (GstElementFactory *factory);

GST_API
const gchar * const *   gst_element_factory_get_uri_protocols   (GstElementFactory *factory);

GST_API
gboolean                gst_element_factory_has_interface       (GstElementFactory *factory,
                                                                 const gchar *interfacename);
GST_API
GstElement*             gst_element_factory_create              (GstElementFactory *factory,
                                                                 const gchar *name) G_GNUC_MALLOC;
GST_API
GstElement*             gst_element_factory_create_full         (GstElementFactory * factory,
                                                                 const gchar * first, ...) G_GNUC_MALLOC;
GST_API
GstElement *            gst_element_factory_create_valist       (GstElementFactory * factory,
                                                                 const gchar * first, va_list properties) G_GNUC_MALLOC;
GST_API
GstElement *            gst_element_factory_create_with_properties (GstElementFactory * factory,
                                                                 guint n, const gchar *names[], const GValue values[]) G_GNUC_MALLOC;
GST_API
GstElement*             gst_element_factory_make                (const gchar *factoryname, const gchar *name) G_GNUC_MALLOC;

GST_API
GstElement*             gst_element_factory_make_full           (const gchar *factoryname,
                                                                  const gchar *first, ...) G_GNUC_MALLOC;
GST_API
GstElement*             gst_element_factory_make_valist         (const gchar *factoryname,
                                                                 const gchar *first, va_list properties) G_GNUC_MALLOC;
GST_API
GstElement*             gst_element_factory_make_with_properties (const gchar *factoryname,
                                                                 guint n, const gchar *names[], const GValue values[]) G_GNUC_MALLOC;
GST_API
gboolean                gst_element_register                    (GstPlugin *plugin, const gchar *name,
                                                                 guint rank, GType type);

GST_API
void                    gst_element_type_set_skip_documentation (GType type);

GST_API
gboolean                gst_element_factory_get_skip_documentation (GstElementFactory * factory);

/* Factory list functions */

/**
 * GstFactoryListType:
 * @GST_ELEMENT_FACTORY_TYPE_DECODER: Decoder elements
 * @GST_ELEMENT_FACTORY_TYPE_ENCODER: Encoder elements
 * @GST_ELEMENT_FACTORY_TYPE_SINK: Sink elements
 * @GST_ELEMENT_FACTORY_TYPE_SRC: Source elements
 * @GST_ELEMENT_FACTORY_TYPE_MUXER: Muxer elements
 * @GST_ELEMENT_FACTORY_TYPE_DEMUXER: Demuxer elements
 * @GST_ELEMENT_FACTORY_TYPE_PARSER: Parser elements
 * @GST_ELEMENT_FACTORY_TYPE_PAYLOADER: Payloader elements
 * @GST_ELEMENT_FACTORY_TYPE_DEPAYLOADER: Depayloader elements
 * @GST_ELEMENT_FACTORY_TYPE_DECRYPTOR: Elements handling decryption (Since: 1.6)
 * @GST_ELEMENT_FACTORY_TYPE_ENCRYPTOR: Elements handling encryption (Since: 1.6)
 * @GST_ELEMENT_FACTORY_TYPE_HARDWARE: Hardware based elements (Since: 1.18)
 * @GST_ELEMENT_FACTORY_TYPE_MAX_ELEMENTS: Private, do not use
 * @GST_ELEMENT_FACTORY_TYPE_MEDIA_VIDEO: Elements handling video media types
 * @GST_ELEMENT_FACTORY_TYPE_MEDIA_AUDIO: Elements handling audio media types
 * @GST_ELEMENT_FACTORY_TYPE_MEDIA_IMAGE: Elements handling image media types
 * @GST_ELEMENT_FACTORY_TYPE_MEDIA_SUBTITLE: Elements handling subtitle media types
 * @GST_ELEMENT_FACTORY_TYPE_MEDIA_METADATA: Elements handling metadata media types
 *
 * The type of #GstElementFactory to filter.
 *
 * All @GstFactoryListType up to @GST_ELEMENT_FACTORY_TYPE_MAX_ELEMENTS are exclusive.
 *
 * If one or more of the MEDIA types are specified, then only elements
 * matching the specified media types will be selected.
 */

/**
 * GstElementFactoryListType:
 *
 * A type defining the type of an element factory.
 */
typedef guint64 GstElementFactoryListType;

#define  GST_ELEMENT_FACTORY_TYPE_DECODER        ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 0))
#define  GST_ELEMENT_FACTORY_TYPE_ENCODER        ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 1))
#define  GST_ELEMENT_FACTORY_TYPE_SINK           ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 2))
#define  GST_ELEMENT_FACTORY_TYPE_SRC            ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 3))
#define  GST_ELEMENT_FACTORY_TYPE_MUXER          ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 4))
#define  GST_ELEMENT_FACTORY_TYPE_DEMUXER        ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 5))
#define  GST_ELEMENT_FACTORY_TYPE_PARSER         ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 6))
#define  GST_ELEMENT_FACTORY_TYPE_PAYLOADER      ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 7))
#define  GST_ELEMENT_FACTORY_TYPE_DEPAYLOADER    ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 8))
#define  GST_ELEMENT_FACTORY_TYPE_FORMATTER      ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 9))
#define  GST_ELEMENT_FACTORY_TYPE_DECRYPTOR      ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 10))
#define  GST_ELEMENT_FACTORY_TYPE_ENCRYPTOR      ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 11))
#define  GST_ELEMENT_FACTORY_TYPE_HARDWARE      ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 12))

#define  GST_ELEMENT_FACTORY_TYPE_MAX_ELEMENTS   ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 48))

#define  GST_ELEMENT_FACTORY_TYPE_MEDIA_VIDEO    ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 49))
#define  GST_ELEMENT_FACTORY_TYPE_MEDIA_AUDIO    ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 50))
#define  GST_ELEMENT_FACTORY_TYPE_MEDIA_IMAGE    ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 51))
#define  GST_ELEMENT_FACTORY_TYPE_MEDIA_SUBTITLE ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 52))
#define  GST_ELEMENT_FACTORY_TYPE_MEDIA_METADATA ((GstElementFactoryListType)(G_GUINT64_CONSTANT (1) << 53))

/**
 * GST_ELEMENT_FACTORY_TYPE_ANY: (value 562949953421311) (type GstElementFactoryListType)
 *
 * Elements of any of the defined GST_ELEMENT_FACTORY_LIST types
 */
#define  GST_ELEMENT_FACTORY_TYPE_ANY ((GstElementFactoryListType)((G_GUINT64_CONSTANT (1) << 49) - 1))

/**
 * GST_ELEMENT_FACTORY_TYPE_MEDIA_ANY: (value 18446462598732840960) (type GstElementFactoryListType)
 *
 * Elements matching any of the defined GST_ELEMENT_FACTORY_TYPE_MEDIA types
 *
 * Note: Do not use this if you wish to not filter against any of the defined
 * media types. If you wish to do this, simply don't specify any
 * GST_ELEMENT_FACTORY_TYPE_MEDIA flag.
 */
#define GST_ELEMENT_FACTORY_TYPE_MEDIA_ANY ((GstElementFactoryListType)(~G_GUINT64_CONSTANT (0) << 48))

/**
 * GST_ELEMENT_FACTORY_TYPE_VIDEO_ENCODER: (value 2814749767106562) (type GstElementFactoryListType)
 *
 * All encoders handling video or image media types
 */
#define GST_ELEMENT_FACTORY_TYPE_VIDEO_ENCODER ((GstElementFactoryListType)(GST_ELEMENT_FACTORY_TYPE_ENCODER | GST_ELEMENT_FACTORY_TYPE_MEDIA_VIDEO | GST_ELEMENT_FACTORY_TYPE_MEDIA_IMAGE))

/**
 * GST_ELEMENT_FACTORY_TYPE_AUDIO_ENCODER: (value 1125899906842626) (type GstElementFactoryListType)
 *
 * All encoders handling audio media types
 */
#define GST_ELEMENT_FACTORY_TYPE_AUDIO_ENCODER ((GstElementFactoryListType)(GST_ELEMENT_FACTORY_TYPE_ENCODER | GST_ELEMENT_FACTORY_TYPE_MEDIA_AUDIO))

/**
 * GST_ELEMENT_FACTORY_TYPE_AUDIOVIDEO_SINKS: (value 3940649673949188) (type GstElementFactoryListType)
 *
 * All sinks handling audio, video or image media types
 */
#define GST_ELEMENT_FACTORY_TYPE_AUDIOVIDEO_SINKS ((GstElementFactoryListType)(GST_ELEMENT_FACTORY_TYPE_SINK | GST_ELEMENT_FACTORY_TYPE_MEDIA_AUDIO | GST_ELEMENT_FACTORY_TYPE_MEDIA_VIDEO | GST_ELEMENT_FACTORY_TYPE_MEDIA_IMAGE))

/**
 * GST_ELEMENT_FACTORY_TYPE_DECODABLE: (value 1377) (type GstElementFactoryListType)
 *
 * All elements used to 'decode' streams (decoders, demuxers, parsers, depayloaders)
 */
#define GST_ELEMENT_FACTORY_TYPE_DECODABLE \
  ((GstElementFactoryListType)(GST_ELEMENT_FACTORY_TYPE_DECODER | GST_ELEMENT_FACTORY_TYPE_DEMUXER | GST_ELEMENT_FACTORY_TYPE_DEPAYLOADER | GST_ELEMENT_FACTORY_TYPE_PARSER | GST_ELEMENT_FACTORY_TYPE_DECRYPTOR))

/* Element klass defines */
#define GST_ELEMENT_FACTORY_KLASS_DECODER               "Decoder"
#define GST_ELEMENT_FACTORY_KLASS_ENCODER               "Encoder"
#define GST_ELEMENT_FACTORY_KLASS_SINK                  "Sink"
#define GST_ELEMENT_FACTORY_KLASS_SRC                   "Source"
#define GST_ELEMENT_FACTORY_KLASS_MUXER                 "Muxer"
#define GST_ELEMENT_FACTORY_KLASS_DEMUXER               "Demuxer"
#define GST_ELEMENT_FACTORY_KLASS_PARSER                "Parser"
#define GST_ELEMENT_FACTORY_KLASS_PAYLOADER             "Payloader"
#define GST_ELEMENT_FACTORY_KLASS_DEPAYLOADER           "Depayloader"
#define GST_ELEMENT_FACTORY_KLASS_FORMATTER             "Formatter"
#define GST_ELEMENT_FACTORY_KLASS_DECRYPTOR             "Decryptor"
#define GST_ELEMENT_FACTORY_KLASS_ENCRYPTOR             "Encryptor"

#define GST_ELEMENT_FACTORY_KLASS_MEDIA_VIDEO           "Video"
#define GST_ELEMENT_FACTORY_KLASS_MEDIA_AUDIO           "Audio"
#define GST_ELEMENT_FACTORY_KLASS_MEDIA_IMAGE           "Image"
#define GST_ELEMENT_FACTORY_KLASS_MEDIA_SUBTITLE        "Subtitle"
#define GST_ELEMENT_FACTORY_KLASS_MEDIA_METADATA        "Metadata"

/**
 * GST_ELEMENT_FACTORY_KLASS_HARDWARE:
 *
 * Elements interacting with hardware devices should specify this classifier in
 * their metadata. You may need to put the element in "READY" state to test if
 * the hardware is present in the system.
 *
 * Since: 1.16
 */
#define GST_ELEMENT_FACTORY_KLASS_HARDWARE              "Hardware"

GST_API
gboolean      gst_element_factory_list_is_type      (GstElementFactory *factory,
                                                     GstElementFactoryListType type);

GST_API
GList *       gst_element_factory_list_get_elements (GstElementFactoryListType type,
                                                     GstRank minrank) G_GNUC_MALLOC;


GST_API
GList *       gst_element_factory_list_filter       (GList *list, const GstCaps *caps,
                                                     GstPadDirection direction,
                                                     gboolean subsetonly) G_GNUC_MALLOC;
G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstElementFactory, gst_object_unref)

G_END_DECLS

#endif /* __GST_ELEMENT_FACTORY_H__ */
