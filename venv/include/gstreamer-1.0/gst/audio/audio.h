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
#define __GST_AUDIO_AUDIO_H__

#include <gst/gst.h>
#include <gst/audio/audio-prelude.h>
#include <gst/audio/audio-enumtypes.h>
#include <gst/audio/audio-format.h>
#include <gst/audio/audio-channels.h>
#include <gst/audio/audio-channel-mixer.h>
#include <gst/audio/audio-info.h>
#include <gst/audio/audio-buffer.h>
#include <gst/audio/audio-quantize.h>
#include <gst/audio/audio-converter.h>
#include <gst/audio/audio-resampler.h>
#include <gst/audio/gstaudiostreamalign.h>
#include <gst/audio/gstaudioaggregator.h>

G_BEGIN_DECLS

/* conversion macros */
/**
 * GST_FRAMES_TO_CLOCK_TIME:
 * @frames: sample frames
 * @rate: sampling rate
 *
 * Calculate clocktime from sample @frames and @rate.
 */
#define GST_FRAMES_TO_CLOCK_TIME(frames, rate) \
  ((GstClockTime) gst_util_uint64_scale_round (frames, GST_SECOND, rate))

/**
 * GST_CLOCK_TIME_TO_FRAMES:
 * @clocktime: clock time
 * @rate: sampling rate
 *
 * Calculate frames from @clocktime and sample @rate.
 */
#define GST_CLOCK_TIME_TO_FRAMES(clocktime, rate) \
  gst_util_uint64_scale_round (clocktime, rate, GST_SECOND)

/* metadata macros */

/**
 * GST_META_TAG_AUDIO_STR:
 *
 * This metadata is relevant for audio streams.
 *
 * Since: 1.2
 */
#define GST_META_TAG_AUDIO_STR "audio"
/**
 * GST_META_TAG_AUDIO_CHANNELS_STR:
 *
 * This metadata stays relevant as long as channels are unchanged.
 *
 * Since: 1.2
 */
#define GST_META_TAG_AUDIO_CHANNELS_STR "channels"

/**
 * GST_META_TAG_AUDIO_RATE_STR:
 *
 * This metadata stays relevant as long as sample rate is unchanged.
 *
 * Since: 1.8
 */
#define GST_META_TAG_AUDIO_RATE_STR "rate"

/*
 * this library defines and implements some helper functions for audio
 * handling
 */

GST_AUDIO_API
GstBuffer *    gst_audio_buffer_clip     (GstBuffer *buffer,
                                          const GstSegment *segment,
                                          gint rate, gint bpf);

GST_AUDIO_API
GstBuffer *    gst_audio_buffer_truncate (GstBuffer *buffer,
                                          gint bpf, gsize trim, gsize samples);

G_END_DECLS

#include <gst/audio/gstaudioringbuffer.h>
#include <gst/audio/gstaudioclock.h>
#include <gst/audio/gstaudiofilter.h>
#include <gst/audio/gstaudiocdsrc.h>
#include <gst/audio/gstaudiodecoder.h>
#include <gst/audio/gstaudioencoder.h>
#include <gst/audio/gstaudiobasesink.h>
#include <gst/audio/gstaudiobasesrc.h>
#include <gst/audio/gstaudiometa.h>
#include <gst/audio/gstaudiosink.h>
#include <gst/audio/gstaudiosrc.h>
#include <gst/audio/streamvolume.h>
#include <gst/audio/gstaudioiec61937.h>

#endif /* __GST_AUDIO_AUDIO_H__ */
