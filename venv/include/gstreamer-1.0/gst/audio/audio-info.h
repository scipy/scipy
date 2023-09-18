/* GStreamer
 * Copyright (C) <1999> Erik Walthinsen <omega@cse.ogi.edu>
 * Library       <2001> Thomas Vander Stichele <thomas@apestaart.org>
 *               <2011> Wim Taymans <wim.taymans@gmail.com>
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

#ifndef __GST_AUDIO_INFO_H__
#define __GST_AUDIO_INFO_H__

G_BEGIN_DECLS

typedef struct _GstAudioInfo GstAudioInfo;

/**
 * GstAudioFlags:
 * @GST_AUDIO_FLAG_NONE: no valid flag
 * @GST_AUDIO_FLAG_UNPOSITIONED: the position array explicitly
 *     contains unpositioned channels.
 *
 * Extra audio flags
 */
typedef enum {
  GST_AUDIO_FLAG_NONE              = 0,
  GST_AUDIO_FLAG_UNPOSITIONED      = (1 << 0)
} GstAudioFlags;

/**
 * GstAudioInfo:
 * @finfo: the format info of the audio
 * @flags: additional audio flags
 * @layout: audio layout
 * @rate: the audio sample rate
 * @channels: the number of channels
 * @bpf: the number of bytes for one frame, this is the size of one
 *         sample * @channels
 * @position: the positions for each channel
 *
 * Information describing audio properties. This information can be filled
 * in from GstCaps with gst_audio_info_from_caps().
 *
 * Use the provided macros to access the info in this structure.
 */
struct _GstAudioInfo {
  const GstAudioFormatInfo *finfo;
  GstAudioFlags             flags;
  GstAudioLayout            layout;
  gint                      rate;
  gint                      channels;
  gint                      bpf;
  GstAudioChannelPosition   position[64];

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

#define GST_TYPE_AUDIO_INFO                  (gst_audio_info_get_type ())
GST_AUDIO_API
GType gst_audio_info_get_type                (void);

#define GST_AUDIO_INFO_IS_VALID(i)           ((i)->finfo != NULL && (i)->rate > 0 && (i)->channels > 0 && (i)->bpf > 0)

#define GST_AUDIO_INFO_FORMAT(i)             (GST_AUDIO_FORMAT_INFO_FORMAT((i)->finfo))
#define GST_AUDIO_INFO_NAME(i)               (GST_AUDIO_FORMAT_INFO_NAME((i)->finfo))
#define GST_AUDIO_INFO_WIDTH(i)              (GST_AUDIO_FORMAT_INFO_WIDTH((i)->finfo))
#define GST_AUDIO_INFO_DEPTH(i)              (GST_AUDIO_FORMAT_INFO_DEPTH((i)->finfo))
#define GST_AUDIO_INFO_BPS(info)             (GST_AUDIO_INFO_DEPTH(info) >> 3)

#define GST_AUDIO_INFO_IS_INTEGER(i)         (GST_AUDIO_FORMAT_INFO_IS_INTEGER((i)->finfo))
#define GST_AUDIO_INFO_IS_FLOAT(i)           (GST_AUDIO_FORMAT_INFO_IS_FLOAT((i)->finfo))
#define GST_AUDIO_INFO_IS_SIGNED(i)          (GST_AUDIO_FORMAT_INFO_IS_SIGNED((i)->finfo))

#define GST_AUDIO_INFO_ENDIANNESS(i)         (GST_AUDIO_FORMAT_INFO_ENDIANNESS((i)->finfo))
#define GST_AUDIO_INFO_IS_LITTLE_ENDIAN(i)   (GST_AUDIO_FORMAT_INFO_IS_LITTLE_ENDIAN((i)->finfo))
#define GST_AUDIO_INFO_IS_BIG_ENDIAN(i)      (GST_AUDIO_FORMAT_INFO_IS_BIG_ENDIAN((i)->finfo))

#define GST_AUDIO_INFO_FLAGS(info)           ((info)->flags)
#define GST_AUDIO_INFO_IS_UNPOSITIONED(info) (((info)->flags & GST_AUDIO_FLAG_UNPOSITIONED) != 0)
#define GST_AUDIO_INFO_LAYOUT(info)          ((info)->layout)

#define GST_AUDIO_INFO_RATE(info)            ((info)->rate)
#define GST_AUDIO_INFO_CHANNELS(info)        ((info)->channels)
#define GST_AUDIO_INFO_BPF(info)             ((info)->bpf)
#define GST_AUDIO_INFO_POSITION(info,c)      ((info)->position[c])

GST_AUDIO_API
GstAudioInfo * gst_audio_info_new         (void);

GST_AUDIO_API
GstAudioInfo * gst_audio_info_new_from_caps (const GstCaps * caps);

GST_AUDIO_API
void           gst_audio_info_init        (GstAudioInfo *info);

GST_AUDIO_API
GstAudioInfo * gst_audio_info_copy        (const GstAudioInfo *info);

GST_AUDIO_API
void           gst_audio_info_free        (GstAudioInfo *info);

GST_AUDIO_API
void           gst_audio_info_set_format  (GstAudioInfo *info, GstAudioFormat format,
                                           gint rate, gint channels,
                                           const GstAudioChannelPosition *position);

GST_AUDIO_API
gboolean       gst_audio_info_from_caps   (GstAudioInfo *info, const GstCaps *caps);

GST_AUDIO_API
GstCaps *      gst_audio_info_to_caps     (const GstAudioInfo *info);

GST_AUDIO_API
gboolean       gst_audio_info_convert     (const GstAudioInfo * info,
                                           GstFormat src_fmt, gint64 src_val,
                                           GstFormat dest_fmt, gint64 * dest_val);

GST_AUDIO_API
gboolean       gst_audio_info_is_equal    (const GstAudioInfo *info,
                                           const GstAudioInfo *other);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioInfo, gst_audio_info_free)

G_END_DECLS

#endif /* __GST_AUDIO_INFO_H__ */
