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

#ifndef __GST_AUDIO_CHANNELS_H__
#define __GST_AUDIO_CHANNELS_H__

G_BEGIN_DECLS

/**
 * GstAudioChannelPosition:
 * @GST_AUDIO_CHANNEL_POSITION_MONO: Mono without direction;
 *     can only be used with 1 channel
 * @GST_AUDIO_CHANNEL_POSITION_FRONT_LEFT: Front left
 * @GST_AUDIO_CHANNEL_POSITION_FRONT_RIGHT: Front right
 * @GST_AUDIO_CHANNEL_POSITION_FRONT_CENTER: Front center
 * @GST_AUDIO_CHANNEL_POSITION_LFE1: Low-frequency effects 1 (subwoofer)
 * @GST_AUDIO_CHANNEL_POSITION_REAR_LEFT: Rear left
 * @GST_AUDIO_CHANNEL_POSITION_REAR_RIGHT: Rear right
 * @GST_AUDIO_CHANNEL_POSITION_FRONT_LEFT_OF_CENTER: Front left of center
 * @GST_AUDIO_CHANNEL_POSITION_FRONT_RIGHT_OF_CENTER: Front right of center
 * @GST_AUDIO_CHANNEL_POSITION_REAR_CENTER: Rear center
 * @GST_AUDIO_CHANNEL_POSITION_LFE2: Low-frequency effects 2 (subwoofer)
 * @GST_AUDIO_CHANNEL_POSITION_SIDE_LEFT: Side left
 * @GST_AUDIO_CHANNEL_POSITION_SIDE_RIGHT: Side right
 * @GST_AUDIO_CHANNEL_POSITION_TOP_FRONT_LEFT: Top front left
 * @GST_AUDIO_CHANNEL_POSITION_TOP_FRONT_RIGHT: Top front right
 * @GST_AUDIO_CHANNEL_POSITION_TOP_FRONT_CENTER: Top front center
 * @GST_AUDIO_CHANNEL_POSITION_TOP_CENTER: Top center
 * @GST_AUDIO_CHANNEL_POSITION_TOP_REAR_LEFT: Top rear left
 * @GST_AUDIO_CHANNEL_POSITION_TOP_REAR_RIGHT: Top rear right
 * @GST_AUDIO_CHANNEL_POSITION_TOP_SIDE_LEFT: Top side right
 * @GST_AUDIO_CHANNEL_POSITION_TOP_SIDE_RIGHT: Top rear right
 * @GST_AUDIO_CHANNEL_POSITION_TOP_REAR_CENTER: Top rear center
 * @GST_AUDIO_CHANNEL_POSITION_BOTTOM_FRONT_CENTER: Bottom front center
 * @GST_AUDIO_CHANNEL_POSITION_BOTTOM_FRONT_LEFT: Bottom front left
 * @GST_AUDIO_CHANNEL_POSITION_BOTTOM_FRONT_RIGHT: Bottom front right
 * @GST_AUDIO_CHANNEL_POSITION_WIDE_LEFT: Wide left (between front left and side left)
 * @GST_AUDIO_CHANNEL_POSITION_WIDE_RIGHT: Wide right (between front right and side right)
 * @GST_AUDIO_CHANNEL_POSITION_SURROUND_LEFT: Surround left (between rear left and side left)
 * @GST_AUDIO_CHANNEL_POSITION_SURROUND_RIGHT: Surround right (between rear right and side right)
 * @GST_AUDIO_CHANNEL_POSITION_NONE: used for position-less channels, e.g.
 *     from a sound card that records 1024 channels; mutually exclusive with
 *     any other channel position
 * @GST_AUDIO_CHANNEL_POSITION_INVALID: invalid position
 *
 * Audio channel positions.
 *
 * These are the channels defined in SMPTE 2036-2-2008
 * Table 1 for 22.2 audio systems with the Surround and Wide channels from
 * DTS Coherent Acoustics (v.1.3.1) and 10.2 and 7.1 layouts. In the caps the
 * actual channel layout is expressed with a channel count and a channel mask,
 * which describes the existing channels. The positions in the bit mask correspond
 * to the enum values.
 * For negotiation it is allowed to have more bits set in the channel mask than
 * the number of channels to specify the allowed channel positions but this is
 * not allowed in negotiated caps. It is not allowed in any situation other
 * than the one mentioned below to have less bits set in the channel mask than
 * the number of channels.
 *
 * @GST_AUDIO_CHANNEL_POSITION_MONO can only be used with a single mono channel that
 * has no direction information and would be mixed into all directional channels.
 * This is expressed in caps by having a single channel and no channel mask.
 *
 * @GST_AUDIO_CHANNEL_POSITION_NONE can only be used if all channels have this position.
 * This is expressed in caps by having a channel mask with no bits set.
 *
 * As another special case it is allowed to have two channels without a channel mask.
 * This implicitly means that this is a stereo stream with a front left and front right
 * channel.
 */
typedef enum {
  /* These get negative indices to allow to use
   * the enum values of the normal cases for the
   * bit-mask position */
  GST_AUDIO_CHANNEL_POSITION_NONE = -3,
  GST_AUDIO_CHANNEL_POSITION_MONO = -2,
  GST_AUDIO_CHANNEL_POSITION_INVALID = -1,

  /* Normal cases */
  GST_AUDIO_CHANNEL_POSITION_FRONT_LEFT = 0,
  GST_AUDIO_CHANNEL_POSITION_FRONT_RIGHT,
  GST_AUDIO_CHANNEL_POSITION_FRONT_CENTER,
  GST_AUDIO_CHANNEL_POSITION_LFE1,
  GST_AUDIO_CHANNEL_POSITION_REAR_LEFT,
  GST_AUDIO_CHANNEL_POSITION_REAR_RIGHT,
  GST_AUDIO_CHANNEL_POSITION_FRONT_LEFT_OF_CENTER,
  GST_AUDIO_CHANNEL_POSITION_FRONT_RIGHT_OF_CENTER,
  GST_AUDIO_CHANNEL_POSITION_REAR_CENTER,
  GST_AUDIO_CHANNEL_POSITION_LFE2,
  GST_AUDIO_CHANNEL_POSITION_SIDE_LEFT,
  GST_AUDIO_CHANNEL_POSITION_SIDE_RIGHT,
  GST_AUDIO_CHANNEL_POSITION_TOP_FRONT_LEFT,
  GST_AUDIO_CHANNEL_POSITION_TOP_FRONT_RIGHT,
  GST_AUDIO_CHANNEL_POSITION_TOP_FRONT_CENTER,
  GST_AUDIO_CHANNEL_POSITION_TOP_CENTER,
  GST_AUDIO_CHANNEL_POSITION_TOP_REAR_LEFT,
  GST_AUDIO_CHANNEL_POSITION_TOP_REAR_RIGHT,
  GST_AUDIO_CHANNEL_POSITION_TOP_SIDE_LEFT,
  GST_AUDIO_CHANNEL_POSITION_TOP_SIDE_RIGHT,
  GST_AUDIO_CHANNEL_POSITION_TOP_REAR_CENTER,
  GST_AUDIO_CHANNEL_POSITION_BOTTOM_FRONT_CENTER,
  GST_AUDIO_CHANNEL_POSITION_BOTTOM_FRONT_LEFT,
  GST_AUDIO_CHANNEL_POSITION_BOTTOM_FRONT_RIGHT,
  GST_AUDIO_CHANNEL_POSITION_WIDE_LEFT,
  GST_AUDIO_CHANNEL_POSITION_WIDE_RIGHT,
  GST_AUDIO_CHANNEL_POSITION_SURROUND_LEFT,
  GST_AUDIO_CHANNEL_POSITION_SURROUND_RIGHT
} GstAudioChannelPosition;

#define GST_AUDIO_CHANNEL_POSITION_MASK(pos) (G_GUINT64_CONSTANT(1)<< GST_AUDIO_CHANNEL_POSITION_ ## pos)

GST_AUDIO_API
gboolean       gst_audio_buffer_reorder_channels (GstBuffer * buffer,
                                                  GstAudioFormat format,
                                                  gint channels,
                                                  const GstAudioChannelPosition * from,
                                                  const GstAudioChannelPosition * to);

GST_AUDIO_API
gboolean       gst_audio_reorder_channels        (gpointer data, gsize size,
                                                  GstAudioFormat format,
                                                  gint channels,
                                                  const GstAudioChannelPosition * from,
                                                  const GstAudioChannelPosition * to);

GST_AUDIO_API
gboolean       gst_audio_channel_positions_to_valid_order (GstAudioChannelPosition *position,
                                                           gint channels);

GST_AUDIO_API
gboolean       gst_audio_check_valid_channel_positions (const GstAudioChannelPosition *position,
                                                        gint channels, gboolean force_order);

GST_AUDIO_API
gboolean       gst_audio_channel_positions_to_mask  (const GstAudioChannelPosition *position,
                                                     gint channels, gboolean force_order,
                                                     guint64 *channel_mask);

GST_AUDIO_API
gboolean       gst_audio_channel_positions_from_mask (gint channels, guint64 channel_mask,
                                                      GstAudioChannelPosition * position);

GST_AUDIO_API
gboolean       gst_audio_get_channel_reorder_map (gint channels,
                                                  const GstAudioChannelPosition * from,
                                                  const GstAudioChannelPosition * to,
                                                  gint *reorder_map);

GST_AUDIO_API
guint64        gst_audio_channel_get_fallback_mask (gint channels);

GST_AUDIO_API
gchar*         gst_audio_channel_positions_to_string (const GstAudioChannelPosition * position,
                                                      gint channels);

G_END_DECLS

#endif /* __GST_AUDIO_CHANNELS_H__ */
