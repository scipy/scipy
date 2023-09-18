/* GStreamer
 * Copyright (C) 2004 Ronald Bultje <rbultje@ronald.bitfreak.net>
 *           (C) 2015 Wim Taymans <wim.taymans@gmail.com>
 *
 * audio-channel-mixer.h: setup of channel conversion matrices
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

#ifndef __GST_AUDIO_CHANNEL_MIXER_H__
#define __GST_AUDIO_CHANNEL_MIXER_H__

#include <gst/gst.h>
#include <gst/audio/audio.h>

typedef struct _GstAudioChannelMixer GstAudioChannelMixer;

/**
 * GstAudioChannelMixerFlags:
 * @GST_AUDIO_CHANNEL_MIXER_FLAGS_NONE: no flag
 * @GST_AUDIO_CHANNEL_MIXER_FLAGS_NON_INTERLEAVED_IN: input channels are not interleaved
 * @GST_AUDIO_CHANNEL_MIXER_FLAGS_NON_INTERLEAVED_OUT: output channels are not interleaved
 * @GST_AUDIO_CHANNEL_MIXER_FLAGS_UNPOSITIONED_IN: input channels are explicitly unpositioned
 * @GST_AUDIO_CHANNEL_MIXER_FLAGS_UNPOSITIONED_OUT: output channels are explicitly unpositioned
 *
 * Flags passed to gst_audio_channel_mixer_new()
 */
typedef enum {
  GST_AUDIO_CHANNEL_MIXER_FLAGS_NONE                = 0,
  GST_AUDIO_CHANNEL_MIXER_FLAGS_NON_INTERLEAVED_IN  = (1 << 0),
  GST_AUDIO_CHANNEL_MIXER_FLAGS_NON_INTERLEAVED_OUT = (1 << 1),
  GST_AUDIO_CHANNEL_MIXER_FLAGS_UNPOSITIONED_IN     = (1 << 2),
  GST_AUDIO_CHANNEL_MIXER_FLAGS_UNPOSITIONED_OUT    = (1 << 3)
} GstAudioChannelMixerFlags;

GST_AUDIO_API
GstAudioChannelMixer * gst_audio_channel_mixer_new   (GstAudioChannelMixerFlags flags,
                                                      GstAudioFormat format,
                                                      gint in_channels,
                                                      GstAudioChannelPosition *in_position,
                                                      gint out_channels,
                                                      GstAudioChannelPosition *out_position);

GST_AUDIO_API
GstAudioChannelMixer * gst_audio_channel_mixer_new_with_matrix (GstAudioChannelMixerFlags flags,
                                                                GstAudioFormat format,
                                                                gint in_channels,
                                                                gint out_channels,
                                                                gfloat **matrix);

GST_AUDIO_API
void                   gst_audio_channel_mixer_free  (GstAudioChannelMixer *mix);

/*
 * Checks for passthrough (= identity matrix).
 */

GST_AUDIO_API
gboolean        gst_audio_channel_mixer_is_passthrough  (GstAudioChannelMixer *mix);

/*
 * Do actual mixing.
 */

GST_AUDIO_API
void            gst_audio_channel_mixer_samples   (GstAudioChannelMixer * mix,
                                                   const gpointer         in[],
                                                   gpointer               out[],
                                                   gint                   samples);

#endif /* __GST_AUDIO_CHANNEL_MIXER_H__ */
