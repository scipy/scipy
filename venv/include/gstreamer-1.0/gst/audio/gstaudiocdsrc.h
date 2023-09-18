/* GStreamer Audio CD Source Base Class
 * Copyright (C) 2005 Tim-Philipp Müller <tim centricular net>
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

#ifndef __GST_AUDIO_AUDIO_H__
#include <gst/audio/audio.h>
#endif

#ifndef __GST_AUDIO_CD_SRC_H__
#define __GST_AUDIO_CD_SRC_H__

#include <gst/gst.h>
#include <gst/base/gstpushsrc.h>

G_BEGIN_DECLS

#define GST_TYPE_AUDIO_CD_SRC            (gst_audio_cd_src_get_type())
#define GST_AUDIO_CD_SRC(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj), GST_TYPE_AUDIO_CD_SRC, GstAudioCdSrc))
#define GST_AUDIO_CD_SRC_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass), GST_TYPE_AUDIO_CD_SRC, GstAudioCdSrcClass))
#define GST_IS_AUDIO_CD_SRC(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj), GST_TYPE_AUDIO_CD_SRC))
#define GST_IS_AUDIO_CD_SRC_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass), GST_TYPE_AUDIO_CD_SRC))
#define GST_AUDIO_CD_SRC_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_AUDIO_CD_SRC, GstAudioCdSrcClass))

typedef struct _GstAudioCdSrc GstAudioCdSrc;
typedef struct _GstAudioCdSrcClass GstAudioCdSrcClass;
typedef struct _GstAudioCdSrcTrack GstAudioCdSrcTrack;
typedef struct _GstAudioCdSrcPrivate GstAudioCdSrcPrivate;

/**
 * GstAudioCdSrcMode:
 * @GST_AUDIO_CD_SRC_MODE_NORMAL     : each single track is a stream
 * @GST_AUDIO_CD_SRC_MODE_CONTINUOUS : the entire disc is a single stream
 *
 * Mode in which the CD audio source operates. Influences timestamping,
 * EOS handling and seeking.
 */
typedef enum {
  GST_AUDIO_CD_SRC_MODE_NORMAL,          /* stream = one track  */
  GST_AUDIO_CD_SRC_MODE_CONTINUOUS       /* stream = whole disc */
} GstAudioCdSrcMode;

/**
 * GstAudioCdSrcTrack:
 * @is_audio: Whether this is an audio track
 * @num: Track number in TOC (usually starts from 1, but not always)
 * @start: The first sector of this track (LBA)
 * @end: The last sector of this track (LBA)
 * @tags: Track-specific tags (e.g. from cd-text information), or NULL
 *
 * CD track abstraction to communicate TOC entries to the base class.
 *
 * This structure is only for use by sub-classed in connection with
 * gst_audio_cd_src_add_track().
 *
 * Applications will be informed of the available tracks via a TOC message
 * on the pipeline's #GstBus instead.
 */
/* FIXME 2.0: remove this struct and pass values directly to _add_track() */
struct _GstAudioCdSrcTrack {
  gboolean     is_audio;      /* TRUE if this is an audio track             */
  guint        num;           /* real track number (usually starts from 1)  */
  guint        start;         /* first sector of track (LBA, not LSN!)      */
  guint        end;           /* last sector of track  (LBA, not LSN!)      */
  GstTagList  *tags;          /* NULL or tags for track (e.g. from cd-text) */

  /*< private >*/
  guint        _gst_reserved1[GST_PADDING/2];
  gpointer     _gst_reserved2[GST_PADDING/2];
};

struct _GstAudioCdSrc {
  GstPushSrc            pushsrc;

  /*< protected >*/ /* for use by sub-classes only */
  GstTagList           *tags;            /* tags that apply to all tracks   */

  /*< private >*/
  GstAudioCdSrcPrivate *priv;

  /*< private >*/
  guint                 _gst_reserved1[GST_PADDING/2];
  gpointer              _gst_reserved2[GST_PADDING/2];
};

/**
 * GstAudioCdSrcClass:
 * @pushsrc_class: the parent class
 * @open: opening the device
 * @close: closing the device
 * @read_sector: reading a sector
 * @get_default_device: getting the default device
 * @probe_devices: probing possible devices
 *
 * Audio CD source base class.
 */
struct _GstAudioCdSrcClass {
  GstPushSrcClass pushsrc_class;

  /* open/close the CD device */
  gboolean    (*open)               (GstAudioCdSrc *src, const gchar *device);
  void        (*close)              (GstAudioCdSrc *src);

  /* read one sector (LBA) */
  GstBuffer * (*read_sector)        (GstAudioCdSrc *src, gint sector);

#if 0
  /* return default device or NULL (optional) */
  gchar *     (*get_default_device) (GstAudioCdSrc *src);

  /* return NULL-terminated string array of CD devices, or NULL (optional) */
  /* FIXME 0.11: reconsider for new probing/device discovery API, remove if in doubt */
  gchar **    (*probe_devices)      (GstAudioCdSrc *src);
#endif

  /*< private >*/
  gpointer       _gst_reserved[GST_PADDING_LARGE];
};

GST_AUDIO_API
GType    gst_audio_cd_src_get_type (void);

GST_AUDIO_API
gboolean gst_audio_cd_src_add_track (GstAudioCdSrc      * src,
                                     GstAudioCdSrcTrack * track);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioCdSrc, gst_object_unref)

G_END_DECLS

#endif /* __GST_AUDIO_CD_SRC_H__ */
