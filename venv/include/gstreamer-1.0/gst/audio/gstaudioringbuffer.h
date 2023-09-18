/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2005 Wim Taymans <wim@fluendo.com>
 *
 * gstaudioringbuffer.h:
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

#ifndef __GST_AUDIO_RING_BUFFER_H__
#define __GST_AUDIO_RING_BUFFER_H__

G_BEGIN_DECLS

#define GST_TYPE_AUDIO_RING_BUFFER             (gst_audio_ring_buffer_get_type())
#define GST_AUDIO_RING_BUFFER(obj)             (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AUDIO_RING_BUFFER,GstAudioRingBuffer))
#define GST_AUDIO_RING_BUFFER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AUDIO_RING_BUFFER,GstAudioRingBufferClass))
#define GST_AUDIO_RING_BUFFER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_AUDIO_RING_BUFFER, GstAudioRingBufferClass))
#define GST_AUDIO_RING_BUFFER_CAST(obj)        ((GstAudioRingBuffer *)obj)
#define GST_IS_AUDIO_RING_BUFFER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AUDIO_RING_BUFFER))
#define GST_IS_AUDIO_RING_BUFFER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AUDIO_RING_BUFFER))

typedef struct _GstAudioRingBuffer GstAudioRingBuffer;
typedef struct _GstAudioRingBufferClass GstAudioRingBufferClass;
typedef struct _GstAudioRingBufferSpec GstAudioRingBufferSpec;

/**
 * GstAudioRingBufferCallback:
 * @rbuf: a #GstAudioRingBuffer
 * @data: (array length=len): target to fill
 * @len: amount to fill
 * @user_data: user data
 *
 * This function is set with gst_audio_ring_buffer_set_callback() and is
 * called to fill the memory at @data with @len bytes of samples.
 */
typedef void (*GstAudioRingBufferCallback) (GstAudioRingBuffer *rbuf, guint8* data, guint len, gpointer user_data);

/**
 * GstAudioRingBufferState:
 * @GST_AUDIO_RING_BUFFER_STATE_STOPPED: The ringbuffer is stopped
 * @GST_AUDIO_RING_BUFFER_STATE_PAUSED: The ringbuffer is paused
 * @GST_AUDIO_RING_BUFFER_STATE_STARTED: The ringbuffer is started
 * @GST_AUDIO_RING_BUFFER_STATE_ERROR: The ringbuffer has encountered an
 *     error after it has been started, e.g. because the device was
 *     disconnected (Since: 1.2)
 *
 * The state of the ringbuffer.
 */
typedef enum {
  GST_AUDIO_RING_BUFFER_STATE_STOPPED,
  GST_AUDIO_RING_BUFFER_STATE_PAUSED,
  GST_AUDIO_RING_BUFFER_STATE_STARTED,
  GST_AUDIO_RING_BUFFER_STATE_ERROR
} GstAudioRingBufferState;

/**
 * GstAudioRingBufferFormatType:
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_RAW: samples in linear or float
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MU_LAW: samples in mulaw
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_A_LAW: samples in alaw
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_IMA_ADPCM: samples in ima adpcm
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MPEG: samples in mpeg audio (but not AAC) format
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_GSM: samples in gsm format
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_IEC958: samples in IEC958 frames (e.g. AC3)
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_AC3: samples in AC3 format
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_EAC3: samples in EAC3 format
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_DTS: samples in DTS format
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MPEG2_AAC: samples in MPEG-2 AAC ADTS format
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MPEG4_AAC: samples in MPEG-4 AAC ADTS format
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MPEG2_AAC_RAW: samples in MPEG-2 AAC raw format (Since: 1.12)
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MPEG4_AAC_RAW: samples in MPEG-4 AAC raw format (Since: 1.12)
 * @GST_AUDIO_RING_BUFFER_FORMAT_TYPE_FLAC: samples in FLAC format (Since: 1.12)
 *
 * The format of the samples in the ringbuffer.
 */
typedef enum
{
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_RAW,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MU_LAW,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_A_LAW,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_IMA_ADPCM,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MPEG,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_GSM,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_IEC958,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_AC3,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_EAC3,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_DTS,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MPEG2_AAC,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MPEG4_AAC,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MPEG2_AAC_RAW,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_MPEG4_AAC_RAW,
  GST_AUDIO_RING_BUFFER_FORMAT_TYPE_FLAC
} GstAudioRingBufferFormatType;

/**
 * GstAudioRingBufferSpec:
 * @caps: The caps that generated the Spec.
 * @type: the sample type
 * @info: the #GstAudioInfo
 * @latency_time: the latency in microseconds
 * @buffer_time: the total buffer size in microseconds
 * @segsize: the size of one segment in bytes
 * @segtotal: the total number of segments
 * @seglatency: number of segments queued in the lower level device,
 *  defaults to segtotal
 *
 * The structure containing the format specification of the ringbuffer.
 */
struct _GstAudioRingBufferSpec
{
  /*< public >*/
  /* in */
  GstCaps  *caps;               /* the caps of the buffer */

  /* in/out */
  GstAudioRingBufferFormatType  type;
  GstAudioInfo                  info;


  guint64  latency_time;        /* the required/actual latency time, this is the
				 * actual the size of one segment and the
				 * minimum possible latency we can achieve. */
  guint64  buffer_time;         /* the required/actual time of the buffer, this is
				 * the total size of the buffer and maximum
				 * latency we can compensate for. */
  gint     segsize;             /* size of one buffer segment in bytes, this value
				 * should be chosen to match latency_time as
				 * well as possible. */
  gint     segtotal;            /* total number of segments, this value is the
				 * number of segments of @segsize and should be
				 * chosen so that it matches buffer_time as
				 * close as possible. */
  /* ABI added 0.10.20 */
  gint     seglatency;          /* number of segments queued in the lower
				 * level device, defaults to segtotal. */

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

#define GST_AUDIO_RING_BUFFER_GET_COND(buf) (&(((GstAudioRingBuffer *)buf)->cond))
#define GST_AUDIO_RING_BUFFER_WAIT(buf)     (g_cond_wait (GST_AUDIO_RING_BUFFER_GET_COND (buf), GST_OBJECT_GET_LOCK (buf)))
#define GST_AUDIO_RING_BUFFER_SIGNAL(buf)   (g_cond_signal (GST_AUDIO_RING_BUFFER_GET_COND (buf)))
#define GST_AUDIO_RING_BUFFER_BROADCAST(buf)(g_cond_broadcast (GST_AUDIO_RING_BUFFER_GET_COND (buf)))

/**
 * GstAudioRingBuffer:
 * @cond: used to signal start/stop/pause/resume actions
 * @open: boolean indicating that the ringbuffer is open
 * @acquired: boolean indicating that the ringbuffer is acquired
 * @memory: data in the ringbuffer
 * @size: size of data in the ringbuffer
 * @spec: format and layout of the ringbuffer data
 * @samples_per_seg: number of samples in one segment
 * @empty_seg: pointer to memory holding one segment of silence samples
 * @state: state of the buffer
 * @segdone: readpointer in the ringbuffer
 * @segbase: segment corresponding to segment 0 (unused)
 * @waiting: is a reader or writer waiting for a free segment
 *
 * The ringbuffer base class structure.
 */
struct _GstAudioRingBuffer {
  GstObject                   object;

  /*< public >*/ /* with LOCK */
  GCond                      cond;
  gboolean                    open;
  gboolean                    acquired;
  guint8                     *memory;
  gsize                       size;
  /*< private >*/
  GstClockTime               *timestamps;
  /*< public >*/ /* with LOCK */
  GstAudioRingBufferSpec      spec;
  gint                        samples_per_seg;
  guint8                     *empty_seg;

  /*< public >*/ /* ATOMIC */
  gint                        state;
  gint                        segdone;
  gint                        segbase;
  gint                        waiting;

  /*< private >*/
  GstAudioRingBufferCallback  callback;
  gpointer                    cb_data;

  gboolean                    need_reorder;
  /* gst[channel_reorder_map[i]] = device[i] */
  gint                        channel_reorder_map[64];

  gboolean                    flushing;
  /* ATOMIC */
  gint                        may_start;
  gboolean                    active;

  GDestroyNotify              cb_data_notify;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING - 1];
};

/**
 * GstAudioRingBufferClass:
 * @parent_class: parent class
 * @open_device:  open the device, don't set any params or allocate anything
 * @acquire: allocate the resources for the ringbuffer using the given spec
 * @release: free resources of the ringbuffer
 * @close_device: close the device
 * @start: start processing of samples
 * @pause: pause processing of samples
 * @resume: resume processing of samples after pause
 * @stop: stop processing of samples
 * @delay: get number of frames queued in device
 * @activate: activate the thread that starts pulling and monitoring the
 * consumed segments in the device.
 * @commit: write samples into the ringbuffer
 * @clear_all: Optional.
 *             Clear the entire ringbuffer.
 *             Subclasses should chain up to the parent implementation to
 *             invoke the default handler.
 *
 * The vmethods that subclasses can override to implement the ringbuffer.
 */
struct _GstAudioRingBufferClass {
  GstObjectClass parent_class;

  /*< public >*/
  gboolean     (*open_device)  (GstAudioRingBuffer *buf);
  gboolean     (*acquire)      (GstAudioRingBuffer *buf, GstAudioRingBufferSpec *spec);
  gboolean     (*release)      (GstAudioRingBuffer *buf);
  gboolean     (*close_device) (GstAudioRingBuffer *buf);

  gboolean     (*start)        (GstAudioRingBuffer *buf);
  gboolean     (*pause)        (GstAudioRingBuffer *buf);
  gboolean     (*resume)       (GstAudioRingBuffer *buf);
  gboolean     (*stop)         (GstAudioRingBuffer *buf);

  guint        (*delay)        (GstAudioRingBuffer *buf);

  /* ABI added */
  gboolean     (*activate)     (GstAudioRingBuffer *buf, gboolean active);

  guint        (*commit)       (GstAudioRingBuffer * buf, guint64 *sample,
                                guint8 * data, gint in_samples,
                                gint out_samples, gint * accum);

  void         (*clear_all)    (GstAudioRingBuffer * buf);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_AUDIO_API
GType gst_audio_ring_buffer_get_type(void);

/* callback stuff */

GST_AUDIO_API
void            gst_audio_ring_buffer_set_callback      (GstAudioRingBuffer *buf,
                                                         GstAudioRingBufferCallback cb,
                                                         gpointer user_data);

GST_AUDIO_API
void            gst_audio_ring_buffer_set_callback_full (GstAudioRingBuffer *buf,
                                                         GstAudioRingBufferCallback cb,
                                                         gpointer user_data,
                                                         GDestroyNotify notify);

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_parse_caps      (GstAudioRingBufferSpec *spec, GstCaps *caps);

GST_AUDIO_API
void            gst_audio_ring_buffer_debug_spec_caps (GstAudioRingBufferSpec *spec);

GST_AUDIO_API
void            gst_audio_ring_buffer_debug_spec_buff (GstAudioRingBufferSpec *spec);

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_convert         (GstAudioRingBuffer * buf, GstFormat src_fmt,
                                                       gint64 src_val, GstFormat dest_fmt,
                                                       gint64 * dest_val);

/* device state */

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_open_device     (GstAudioRingBuffer *buf);

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_close_device    (GstAudioRingBuffer *buf);

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_device_is_open  (GstAudioRingBuffer *buf);

/* allocate resources */

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_acquire         (GstAudioRingBuffer *buf, GstAudioRingBufferSpec *spec);

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_release         (GstAudioRingBuffer *buf);

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_is_acquired     (GstAudioRingBuffer *buf);

/* set the device channel positions */

GST_AUDIO_API
void            gst_audio_ring_buffer_set_channel_positions (GstAudioRingBuffer *buf, const GstAudioChannelPosition *position);

/* activating */

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_activate        (GstAudioRingBuffer *buf, gboolean active);

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_is_active       (GstAudioRingBuffer *buf);

/* flushing */

GST_AUDIO_API
void            gst_audio_ring_buffer_set_flushing    (GstAudioRingBuffer *buf, gboolean flushing);

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_is_flushing     (GstAudioRingBuffer *buf);

/* playback/pause */

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_start           (GstAudioRingBuffer *buf);

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_pause           (GstAudioRingBuffer *buf);

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_stop            (GstAudioRingBuffer *buf);

/* get status */

GST_AUDIO_API
guint           gst_audio_ring_buffer_delay           (GstAudioRingBuffer *buf);

GST_AUDIO_API
guint64         gst_audio_ring_buffer_samples_done    (GstAudioRingBuffer *buf);

GST_AUDIO_API
void            gst_audio_ring_buffer_set_sample      (GstAudioRingBuffer *buf, guint64 sample);

/* clear all segments */

GST_AUDIO_API
void            gst_audio_ring_buffer_clear_all       (GstAudioRingBuffer *buf);

/* commit samples */

GST_AUDIO_API
guint           gst_audio_ring_buffer_commit          (GstAudioRingBuffer * buf, guint64 *sample,
                                                       guint8 * data, gint in_samples,
                                                       gint out_samples, gint * accum);

/* read samples */

GST_AUDIO_API
guint           gst_audio_ring_buffer_read            (GstAudioRingBuffer *buf, guint64 sample,
                                                       guint8 *data, guint len, GstClockTime *timestamp);

/* Set timestamp on buffer */

GST_AUDIO_API
void            gst_audio_ring_buffer_set_timestamp   (GstAudioRingBuffer * buf, gint readseg, GstClockTime 
                                                       timestamp);

/* mostly protected */
/* not yet implemented
gboolean        gst_audio_ring_buffer_prepare_write   (GstAudioRingBuffer *buf, gint *segment, guint8 **writeptr, gint *len);
*/

GST_AUDIO_API
gboolean        gst_audio_ring_buffer_prepare_read    (GstAudioRingBuffer *buf, gint *segment,
                                                       guint8 **readptr, gint *len);

GST_AUDIO_API
void            gst_audio_ring_buffer_clear           (GstAudioRingBuffer *buf, gint segment);

GST_AUDIO_API
void            gst_audio_ring_buffer_advance         (GstAudioRingBuffer *buf, guint advance);

GST_AUDIO_API
void            gst_audio_ring_buffer_may_start       (GstAudioRingBuffer *buf, gboolean allowed);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioRingBuffer, gst_object_unref)

G_END_DECLS

#endif /* __GST_AUDIO_RING_BUFFER_H__ */
