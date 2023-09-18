/* GStreamer
 * Copyright (C) 2003 Benjamin Otte <in7y118@public.uni-hamburg.de>
 *
 * gsttaglist.h: Header for tag support
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


#ifndef __GST_TAGLIST_H__
#define __GST_TAGLIST_H__

#include <gst/gstdatetime.h>
#include <gst/gstsample.h>
#include <gst/gstbuffer.h>
#include <gst/glib-compat.h>

G_BEGIN_DECLS

/**
 * GstTagMergeMode:
 * @GST_TAG_MERGE_UNDEFINED: undefined merge mode
 * @GST_TAG_MERGE_REPLACE_ALL: replace all tags (clear list and append)
 * @GST_TAG_MERGE_REPLACE: replace tags
 * @GST_TAG_MERGE_APPEND: append tags
 * @GST_TAG_MERGE_PREPEND: prepend tags
 * @GST_TAG_MERGE_KEEP: keep existing tags
 * @GST_TAG_MERGE_KEEP_ALL: keep all existing tags
 * @GST_TAG_MERGE_COUNT: the number of merge modes
 *
 * The different tag merging modes are basically replace, overwrite and append,
 * but they can be seen from two directions. Given two taglists: (A) the tags
 * already in the element and (B) the ones that are supplied to the element (
 * e.g. via gst_tag_setter_merge_tags() / gst_tag_setter_add_tags() or a
 * %GST_EVENT_TAG), how are these tags merged?
 * In the table below this is shown for the cases that a tag exists in the list
 * (A) or does not exists (!A) and combinations thereof.
 *
 * | merge mode  | A + B | A + !B | !A + B | !A + !B |
 * | ----------- | ----- | ------ | ------ | ------- |
 * | REPLACE_ALL | B     | ø      | B      | ø       |
 * | REPLACE     | B     | A      | B      | ø       |
 * | APPEND      | A, B  | A      | B      | ø       |
 * | PREPEND     | B, A  | A      | B      | ø       |
 * | KEEP        | A     | A      | B      | ø       |
 * | KEEP_ALL    | A     | A      | ø      | ø       |
 */

typedef enum {
  GST_TAG_MERGE_UNDEFINED,
  GST_TAG_MERGE_REPLACE_ALL,
  GST_TAG_MERGE_REPLACE,
  GST_TAG_MERGE_APPEND,
  GST_TAG_MERGE_PREPEND,
  GST_TAG_MERGE_KEEP,
  GST_TAG_MERGE_KEEP_ALL,
  /* add more */
  GST_TAG_MERGE_COUNT
} GstTagMergeMode;

#define GST_TAG_MODE_IS_VALID(mode)     (((mode) > GST_TAG_MERGE_UNDEFINED) && ((mode) < GST_TAG_MERGE_COUNT))

/**
 * GstTagFlag:
 * @GST_TAG_FLAG_UNDEFINED: undefined flag
 * @GST_TAG_FLAG_META: tag is meta data
 * @GST_TAG_FLAG_ENCODED: tag is encoded
 * @GST_TAG_FLAG_DECODED: tag is decoded
 * @GST_TAG_FLAG_COUNT: number of tag flags
 *
 * Extra tag flags used when registering tags.
 */
/* FIXME: these are not really flags .. */
typedef enum {
  GST_TAG_FLAG_UNDEFINED,
  GST_TAG_FLAG_META,
  GST_TAG_FLAG_ENCODED,
  GST_TAG_FLAG_DECODED,
  GST_TAG_FLAG_COUNT
} GstTagFlag;

#define GST_TAG_FLAG_IS_VALID(flag)     (((flag) > GST_TAG_FLAG_UNDEFINED) && ((flag) < GST_TAG_FLAG_COUNT))

/**
 * GstTagList:
 * @mini_object: the parent type
 *
 * Object describing tags / metadata.
 */
typedef struct _GstTagList GstTagList;
struct _GstTagList {
  GstMiniObject mini_object;
};

GST_API GType _gst_tag_list_type;

#define GST_TAG_LIST(x)       ((GstTagList *) (x))
#define GST_TYPE_TAG_LIST     (_gst_tag_list_type)
#define GST_IS_TAG_LIST(obj)  (GST_IS_MINI_OBJECT_TYPE((obj), GST_TYPE_TAG_LIST))

/**
 * GstTagForeachFunc:
 * @list: the #GstTagList
 * @tag: a name of a tag in @list
 * @user_data: user data
 *
 * A function that will be called in gst_tag_list_foreach(). The function may
 * not modify the tag list.
 */
typedef void (*GstTagForeachFunc) (const GstTagList * list,
                                   const gchar      * tag,
                                   gpointer           user_data);

/**
 * GstTagMergeFunc:
 * @dest: the destination #GValue
 * @src: the source #GValue
 *
 * A function for merging multiple values of a tag used when registering
 * tags.
 */
typedef void (* GstTagMergeFunc) (GValue *dest, const GValue *src);

GST_API
GType        gst_tag_list_get_type (void);

/* tag registration */

GST_API
void         gst_tag_register      (const gchar     * name,
                                    GstTagFlag        flag,
                                    GType             type,
                                    const gchar     * nick,
                                    const gchar     * blurb,
                                    GstTagMergeFunc   func);
GST_API
void         gst_tag_register_static (const gchar   * name,
                                      GstTagFlag      flag,
                                      GType           type,
                                      const gchar   * nick,
                                      const gchar   * blurb,
                                      GstTagMergeFunc func);

/* some default merging functions */

GST_API
void      gst_tag_merge_use_first          (GValue * dest, const GValue * src);

GST_API
void      gst_tag_merge_strings_with_comma (GValue * dest, const GValue * src);

/* basic tag support */

GST_API
gboolean               gst_tag_exists          (const gchar * tag);

GST_API
GType                  gst_tag_get_type        (const gchar * tag);

GST_API
const gchar *          gst_tag_get_nick        (const gchar * tag);

GST_API
const gchar *          gst_tag_get_description (const gchar * tag);

GST_API
GstTagFlag             gst_tag_get_flag        (const gchar * tag);

GST_API
gboolean               gst_tag_is_fixed        (const gchar * tag);

/* tag lists */

/**
 * GstTagScope:
 * @GST_TAG_SCOPE_STREAM: tags specific to this single stream
 * @GST_TAG_SCOPE_GLOBAL: global tags for the complete medium
 *
 * GstTagScope specifies if a taglist applies to the complete
 * medium or only to one single stream.
 */
typedef enum {
  GST_TAG_SCOPE_STREAM,
  GST_TAG_SCOPE_GLOBAL
} GstTagScope;

GST_API
GstTagList * gst_tag_list_new_empty         (void) G_GNUC_MALLOC;

GST_API
GstTagList * gst_tag_list_new               (const gchar * tag, ...) G_GNUC_NULL_TERMINATED G_GNUC_MALLOC;

GST_API
GstTagList * gst_tag_list_new_valist        (va_list var_args) G_GNUC_MALLOC;

GST_API
void         gst_tag_list_set_scope         (GstTagList * list, GstTagScope scope);

GST_API
GstTagScope  gst_tag_list_get_scope         (const GstTagList * list);

GST_API
gchar      * gst_tag_list_to_string         (const GstTagList * list) G_GNUC_MALLOC;

GST_API
GstTagList * gst_tag_list_new_from_string   (const gchar      * str) G_GNUC_MALLOC;

GST_API
gint         gst_tag_list_n_tags            (const GstTagList * list);

GST_API
const gchar* gst_tag_list_nth_tag_name      (const GstTagList * list, guint index);

GST_API
gboolean     gst_tag_list_is_empty          (const GstTagList * list);

GST_API
gboolean     gst_tag_list_is_equal          (const GstTagList * list1,
                                             const GstTagList * list2);
GST_API
void         gst_tag_list_insert            (GstTagList       * into,
                                             const GstTagList * from,
                                             GstTagMergeMode    mode);
GST_API
GstTagList * gst_tag_list_merge             (const GstTagList * list1,
                                             const GstTagList * list2,
                                             GstTagMergeMode    mode) G_GNUC_MALLOC;
GST_API
guint        gst_tag_list_get_tag_size      (const GstTagList * list,
                                             const gchar      * tag);
GST_API
void         gst_tag_list_add               (GstTagList       * list,
                                             GstTagMergeMode    mode,
                                             const gchar      * tag,
                                             ...) G_GNUC_NULL_TERMINATED;
GST_API
void         gst_tag_list_add_values        (GstTagList       * list,
                                             GstTagMergeMode    mode,
                                             const gchar      * tag,
                                             ...) G_GNUC_NULL_TERMINATED;
GST_API
void         gst_tag_list_add_valist        (GstTagList       * list,
                                             GstTagMergeMode    mode,
                                             const gchar      * tag,
                                             va_list        var_args);
GST_API
void         gst_tag_list_add_valist_values (GstTagList       * list,
                                             GstTagMergeMode    mode,
                                             const gchar      * tag,
                                             va_list            var_args);
GST_API
void         gst_tag_list_add_value         (GstTagList       * list,
                                             GstTagMergeMode    mode,
                                             const gchar      * tag,
                                             const GValue     * value);
GST_API
void         gst_tag_list_remove_tag        (GstTagList       * list,
                                             const gchar      * tag);
GST_API
void         gst_tag_list_foreach           (const GstTagList * list,
                                             GstTagForeachFunc  func,
                                             gpointer           user_data);
GST_API
const GValue *
             gst_tag_list_get_value_index   (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index);
GST_API
gboolean     gst_tag_list_copy_value        (GValue           * dest,
                                             const GstTagList * list,
                                             const gchar      * tag);

/* simplifications (FIXME: do we want them?) */

GST_API
gboolean     gst_tag_list_get_boolean       (const GstTagList * list,
                                             const gchar      * tag,
                                             gboolean         * value);
GST_API
gboolean     gst_tag_list_get_boolean_index (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             gboolean         * value);
GST_API
gboolean     gst_tag_list_get_int           (const GstTagList * list,
                                             const gchar      * tag,
                                             gint             * value);
GST_API
gboolean     gst_tag_list_get_int_index     (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             gint             * value);
GST_API
gboolean     gst_tag_list_get_uint          (const GstTagList * list,
                                             const gchar      * tag,
                                             guint            * value);
GST_API
gboolean     gst_tag_list_get_uint_index    (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             guint            * value);
GST_API
gboolean     gst_tag_list_get_int64         (const GstTagList * list,
                                             const gchar      * tag,
                                             gint64           * value);
GST_API
gboolean     gst_tag_list_get_int64_index   (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             gint64           * value);
GST_API
gboolean     gst_tag_list_get_uint64        (const GstTagList * list,
                                             const gchar      * tag,
                                             guint64          * value);
GST_API
gboolean     gst_tag_list_get_uint64_index  (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             guint64          * value);
GST_API
gboolean     gst_tag_list_get_float         (const GstTagList * list,
                                             const gchar      * tag,
                                             gfloat           * value);
GST_API
gboolean     gst_tag_list_get_float_index   (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             gfloat           * value);
GST_API
gboolean     gst_tag_list_get_double        (const GstTagList * list,
                                             const gchar      * tag,
                                             gdouble          * value);
GST_API
gboolean     gst_tag_list_get_double_index  (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             gdouble          * value);
GST_API
gboolean     gst_tag_list_get_string        (const GstTagList * list,
                                             const gchar      * tag,
                                             gchar           ** value);
GST_API
gboolean     gst_tag_list_get_string_index  (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             gchar           ** value);
GST_API
gboolean     gst_tag_list_peek_string_index (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             const gchar     ** value);
GST_API
gboolean     gst_tag_list_get_pointer       (const GstTagList * list,
                                             const gchar      * tag,
                                             gpointer         * value);
GST_API
gboolean     gst_tag_list_get_pointer_index (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             gpointer         * value);
GST_API
gboolean     gst_tag_list_get_date          (const GstTagList * list,
                                             const gchar      * tag,
                                             GDate           ** value);
GST_API
gboolean     gst_tag_list_get_date_index    (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             GDate           ** value);
GST_API
gboolean     gst_tag_list_get_date_time     (const GstTagList * list,
                                             const gchar      * tag,
                                             GstDateTime     ** value);
GST_API
gboolean     gst_tag_list_get_date_time_index (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             GstDateTime     ** value);
GST_API
gboolean     gst_tag_list_get_sample        (const GstTagList * list,
                                             const gchar      * tag,
                                             GstSample       ** sample);
GST_API
gboolean     gst_tag_list_get_sample_index  (const GstTagList * list,
                                             const gchar      * tag,
                                             guint              index,
                                             GstSample       ** sample);

#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
/* refcounting */
static inline GstTagList *
gst_tag_list_ref (GstTagList * taglist)
{
  return (GstTagList *) gst_mini_object_ref (GST_MINI_OBJECT_CAST (taglist));
}

static inline void
gst_tag_list_unref (GstTagList * taglist)
{
  gst_mini_object_unref (GST_MINI_OBJECT_CAST (taglist));
}

static inline void
gst_clear_tag_list (GstTagList ** taglist_ptr)
{
  gst_clear_mini_object ((GstMiniObject **) taglist_ptr);
}
#else /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
GstTagList *  gst_tag_list_ref   (GstTagList * taglist);

GST_API
void          gst_tag_list_unref (GstTagList * taglist);

GST_API
void          gst_clear_tag_list (GstTagList ** taglist_ptr);
#endif /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */

GST_API
GstTagList* gst_tag_list_copy(const GstTagList* taglist);

#define gst_tag_list_copy(taglist) GST_TAG_LIST (gst_mini_object_copy (GST_MINI_OBJECT_CAST (taglist)))

#ifndef GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS
static inline gboolean
gst_tag_list_replace (GstTagList **old_taglist, GstTagList *new_taglist)
{
    return gst_mini_object_replace ((GstMiniObject **) old_taglist,
        (GstMiniObject *) new_taglist);
}

static inline gboolean
gst_tag_list_take (GstTagList **old_taglist, GstTagList *new_taglist)
{
  return gst_mini_object_take ((GstMiniObject **) old_taglist,
      (GstMiniObject *) new_taglist);
}
#else /* GST_DISABLE_MINIOBJECT_INLINE_FUNCTIONS */
GST_API
gboolean  gst_tag_list_replace (GstTagList ** old_taglist,
                                GstTagList * new_taglist);

GST_API
gboolean  gst_tag_list_take    (GstTagList ** old_taglist,
                                GstTagList * new_taglist);
#endif

/**
 * gst_tag_list_is_writable:
 * @taglist: a #GstTagList
 *
 * Tests if you can safely modify @taglist. It is only safe to modify taglist
 * when there is only one owner of the taglist - ie, the refcount is 1.
 */
#define gst_tag_list_is_writable(taglist)    gst_mini_object_is_writable (GST_MINI_OBJECT_CAST (taglist))

/**
 * gst_tag_list_make_writable:
 * @taglist: (transfer full): a #GstTagList
 *
 * Returns a writable copy of @taglist.
 *
 * If there is only one reference count on @taglist, the caller must be the
 * owner, and so this function will return the taglist object unchanged. If on
 * the other hand there is more than one reference on the object, a new taglist
 * object will be returned (which will be a copy of @taglist). The caller's
 * reference on @taglist will be removed, and instead the caller will own a
 * reference to the returned object.
 *
 * In short, this function unrefs the taglist in the argument and refs the
 * taglist that it returns. Don't access the argument after calling this
 * function. See also: gst_tag_list_ref().
 *
 * Returns: (transfer full): a writable taglist which may or may not be the
 *     same as @taglist
 */
#define gst_tag_list_make_writable(taglist)   GST_TAG_LIST (gst_mini_object_make_writable (GST_MINI_OBJECT_CAST (taglist)))

/* GStreamer core tags */
/**
 * GST_TAG_TITLE:
 *
 * commonly used title (string)
 *
 * The title as it should be displayed, e.g. 'The Doll House'
 */
#define GST_TAG_TITLE                  "title"
/**
 * GST_TAG_TITLE_SORTNAME:
 *
 * commonly used title, as used for sorting (string)
 *
 * The title as it should be sorted, e.g. 'Doll House, The'
 */
#define GST_TAG_TITLE_SORTNAME         "title-sortname"
/**
 * GST_TAG_ARTIST:
 *
 * person(s) responsible for the recording (string)
 *
 * The artist name as it should be displayed, e.g. 'Jimi Hendrix' or
 * 'The Guitar Heroes'
 */
#define GST_TAG_ARTIST                 "artist"
/**
 * GST_TAG_ARTIST_SORTNAME:
 *
 * person(s) responsible for the recording, as used for sorting (string)
 *
 * The artist name as it should be sorted, e.g. 'Hendrix, Jimi' or
 * 'Guitar Heroes, The'
 */
#define GST_TAG_ARTIST_SORTNAME        "artist-sortname"
/**
 * GST_TAG_ALBUM:
 *
 * album containing this data (string)
 *
 * The album name as it should be displayed, e.g. 'The Jazz Guitar'
 */
#define GST_TAG_ALBUM                  "album"
/**
 * GST_TAG_ALBUM_SORTNAME:
 *
 * album containing this data, as used for sorting (string)
 *
 * The album name as it should be sorted, e.g. 'Jazz Guitar, The'
 */
#define GST_TAG_ALBUM_SORTNAME         "album-sortname"
/**
 * GST_TAG_ALBUM_ARTIST:
 *
 * The artist of the entire album, as it should be displayed.
 */
#define GST_TAG_ALBUM_ARTIST           "album-artist"
/**
 * GST_TAG_ALBUM_ARTIST_SORTNAME:
 *
 * The artist of the entire album, as it should be sorted.
 */
#define GST_TAG_ALBUM_ARTIST_SORTNAME  "album-artist-sortname"
/**
 * GST_TAG_COMPOSER:
 *
 * person(s) who composed the recording (string)
 */
#define GST_TAG_COMPOSER               "composer"
/**
 * GST_TAG_CONDUCTOR:
 *
 * conductor/performer refinement (string)
 *
 * Since: 1.8
 */
#define GST_TAG_CONDUCTOR               "conductor"
/**
 * GST_TAG_DATE:
 *
 * date the data was created (#GDate structure)
 */
#define GST_TAG_DATE                   "date"
/**
 * GST_TAG_DATE_TIME:
 *
 * date and time the data was created (#GstDateTime structure)
 */
#define GST_TAG_DATE_TIME              "datetime"
/**
 * GST_TAG_GENRE:
 *
 * genre this data belongs to (string)
 */
#define GST_TAG_GENRE                  "genre"
/**
 * GST_TAG_COMMENT:
 *
 * free text commenting the data (string)
 */
#define GST_TAG_COMMENT                "comment"
/**
 * GST_TAG_EXTENDED_COMMENT:
 *
 * key/value text commenting the data (string)
 *
 * Must be in the form of 'key=comment' or
 * 'key[lc]=comment' where 'lc' is an ISO-639
 * language code.
 *
 * This tag is used for unknown Vorbis comment tags,
 * unknown APE tags and certain ID3v2 comment fields.
 */
#define GST_TAG_EXTENDED_COMMENT       "extended-comment"
/**
 * GST_TAG_TRACK_NUMBER:
 *
 * track number inside a collection (unsigned integer)
 */
#define GST_TAG_TRACK_NUMBER           "track-number"
/**
 * GST_TAG_TRACK_COUNT:
 *
 * count of tracks inside collection this track belongs to (unsigned integer)
 */
#define GST_TAG_TRACK_COUNT            "track-count"
/**
 * GST_TAG_ALBUM_VOLUME_NUMBER:
 *
 * disc number inside a collection (unsigned integer)
 */
#define GST_TAG_ALBUM_VOLUME_NUMBER    "album-disc-number"
/**
 * GST_TAG_ALBUM_VOLUME_COUNT:
 *
 * count of discs inside collection this disc belongs to (unsigned integer)
 */
#define GST_TAG_ALBUM_VOLUME_COUNT    "album-disc-count"
/**
 * GST_TAG_LOCATION:
 *
 * Origin of media as a URI (location, where the original of the file or stream
 * is hosted) (string)
 */
#define GST_TAG_LOCATION               "location"
/**
 * GST_TAG_HOMEPAGE:
 *
 * Homepage for this media (i.e. artist or movie homepage) (string)
 */
#define GST_TAG_HOMEPAGE               "homepage"
/**
 * GST_TAG_DESCRIPTION:
 *
 * short text describing the content of the data (string)
 */
#define GST_TAG_DESCRIPTION            "description"
/**
 * GST_TAG_VERSION:
 *
 * version of this data (string)
 */
#define GST_TAG_VERSION                "version"
/**
 * GST_TAG_ISRC:
 *
 * International Standard Recording Code - see http://www.ifpi.org/isrc/ (string)
 */
#define GST_TAG_ISRC                   "isrc"
/**
 * GST_TAG_ORGANIZATION:
 *
 * organization (string)
 */
#define GST_TAG_ORGANIZATION           "organization"
/**
 * GST_TAG_COPYRIGHT:
 *
 * copyright notice of the data (string)
 */
#define GST_TAG_COPYRIGHT              "copyright"
/**
 * GST_TAG_COPYRIGHT_URI:
 *
 * URI to location where copyright details can be found (string)
 */
#define GST_TAG_COPYRIGHT_URI          "copyright-uri"
/**
 * GST_TAG_ENCODED_BY:
 *
 * name of the person or organisation that encoded the file. May contain a
 * copyright message if the person or organisation also holds the copyright
 * (string)
 *
 * Note: do not use this field to describe the encoding application. Use
 * #GST_TAG_APPLICATION_NAME or #GST_TAG_COMMENT for that.
 */
#define GST_TAG_ENCODED_BY             "encoded-by"
/**
 * GST_TAG_CONTACT:
 *
 * contact information (string)
 */
#define GST_TAG_CONTACT                "contact"
/**
 * GST_TAG_LICENSE:
 *
 * license of data (string)
 */
#define GST_TAG_LICENSE                "license"
/**
 * GST_TAG_LICENSE_URI:
 *
 * URI to location where license details can be found (string)
 */
#define GST_TAG_LICENSE_URI            "license-uri"
/**
 * GST_TAG_PERFORMER:
 *
 * person(s) performing (string)
 */
#define GST_TAG_PERFORMER              "performer"
/**
 * GST_TAG_DURATION:
 *
 * length in GStreamer time units (nanoseconds) (unsigned 64-bit integer)
 */
#define GST_TAG_DURATION               "duration"
/**
 * GST_TAG_CODEC:
 *
 * codec the data is stored in (string)
 */
#define GST_TAG_CODEC                  "codec"
/**
 * GST_TAG_VIDEO_CODEC:
 *
 * codec the video data is stored in (string)
 */
#define GST_TAG_VIDEO_CODEC            "video-codec"
/**
 * GST_TAG_AUDIO_CODEC:
 *
 * codec the audio data is stored in (string)
 */
#define GST_TAG_AUDIO_CODEC            "audio-codec"
/**
 * GST_TAG_SUBTITLE_CODEC:
 *
 * codec/format the subtitle data is stored in (string)
 */
#define GST_TAG_SUBTITLE_CODEC         "subtitle-codec"
/**
 * GST_TAG_CONTAINER_FORMAT:
 *
 * container format the data is stored in (string)
 */
#define GST_TAG_CONTAINER_FORMAT       "container-format"
/**
 * GST_TAG_BITRATE:
 *
 * exact or average bitrate in bits/s (unsigned integer)
 */
#define GST_TAG_BITRATE                "bitrate"
/**
 * GST_TAG_NOMINAL_BITRATE:
 *
 * nominal bitrate in bits/s (unsigned integer). The actual bitrate might be
 * different from this target bitrate.
 */
#define GST_TAG_NOMINAL_BITRATE        "nominal-bitrate"
/**
 * GST_TAG_MINIMUM_BITRATE:
 *
 * minimum bitrate in bits/s (unsigned integer)
 */
#define GST_TAG_MINIMUM_BITRATE        "minimum-bitrate"
/**
 * GST_TAG_MAXIMUM_BITRATE:
 *
 * maximum bitrate in bits/s (unsigned integer)
 */
#define GST_TAG_MAXIMUM_BITRATE        "maximum-bitrate"
/**
 * GST_TAG_SERIAL:
 *
 * serial number of track (unsigned integer)
 */
#define GST_TAG_SERIAL                 "serial"
/**
 * GST_TAG_ENCODER:
 *
 * encoder used to encode this stream (string)
 */
#define GST_TAG_ENCODER                "encoder"
/**
 * GST_TAG_ENCODER_VERSION:
 *
 * version of the encoder used to encode this stream (unsigned integer)
 */
#define GST_TAG_ENCODER_VERSION        "encoder-version"
/**
 * GST_TAG_TRACK_GAIN:
 *
 * track gain in db (double)
 */
#define GST_TAG_TRACK_GAIN             "replaygain-track-gain"
/**
 * GST_TAG_TRACK_PEAK:
 *
 * peak of the track (double)
 */
#define GST_TAG_TRACK_PEAK             "replaygain-track-peak"
/**
 * GST_TAG_ALBUM_GAIN:
 *
 * album gain in db (double)
 */
#define GST_TAG_ALBUM_GAIN             "replaygain-album-gain"
/**
 * GST_TAG_ALBUM_PEAK:
 *
 * peak of the album (double)
 */
#define GST_TAG_ALBUM_PEAK             "replaygain-album-peak"
/**
 * GST_TAG_REFERENCE_LEVEL:
 *
 * reference level of track and album gain values (double)
 */
#define GST_TAG_REFERENCE_LEVEL        "replaygain-reference-level"
/**
 * GST_TAG_LANGUAGE_CODE:
 *
 * ISO-639-2 or ISO-639-1 code for the language the content is in (string)
 *
 * There is utility API in libgsttag in gst-plugins-base to obtain a translated
 * language name from the language code: `gst_tag_get_language_name()`
 */
#define GST_TAG_LANGUAGE_CODE          "language-code"
/**
 * GST_TAG_LANGUAGE_NAME:
 *
 * Name of the language the content is in (string)
 *
 * Free-form name of the language the content is in, if a language code
 * is not available. This tag should not be set in addition to a language
 * code. It is undefined what language or locale the language name is in.
 */
#define GST_TAG_LANGUAGE_NAME          "language-name"
/**
 * GST_TAG_IMAGE:
 *
 * image (sample) (sample taglist should specify the content type and preferably
 * also set "image-type" field as `GstTagImageType`)
 */
#define GST_TAG_IMAGE                  "image"
/**
 * GST_TAG_PREVIEW_IMAGE:
 *
 * image that is meant for preview purposes, e.g. small icon-sized version
 * (sample) (sample taglist should specify the content type)
 */
#define GST_TAG_PREVIEW_IMAGE          "preview-image"

/**
 * GST_TAG_ATTACHMENT:
 *
 * generic file attachment (sample) (sample taglist should specify the content
 * type and if possible set "filename" to the file name of the
 * attachment)
 */
#define GST_TAG_ATTACHMENT             "attachment"

/**
 * GST_TAG_BEATS_PER_MINUTE:
 *
 * number of beats per minute in audio (double)
 */
#define GST_TAG_BEATS_PER_MINUTE       "beats-per-minute"

/**
 * GST_TAG_KEYWORDS:
 *
 * comma separated keywords describing the content (string).
 */
#define GST_TAG_KEYWORDS               "keywords"

/**
 * GST_TAG_GEO_LOCATION_NAME:
 *
 * human readable descriptive location of where the media has been recorded or
 * produced. (string).
 */
#define GST_TAG_GEO_LOCATION_NAME               "geo-location-name"

/**
 * GST_TAG_GEO_LOCATION_LATITUDE:
 *
 * geo latitude location of where the media has been recorded or produced in
 * degrees according to WGS84 (zero at the equator, negative values for southern
 * latitudes) (double).
 */
#define GST_TAG_GEO_LOCATION_LATITUDE               "geo-location-latitude"

/**
 * GST_TAG_GEO_LOCATION_LONGITUDE:
 *
 * geo longitude location of where the media has been recorded or produced in
 * degrees according to WGS84 (zero at the prime meridian in Greenwich/UK,
 * negative values for western longitudes). (double).
 */
#define GST_TAG_GEO_LOCATION_LONGITUDE               "geo-location-longitude"

/**
 * GST_TAG_GEO_LOCATION_ELEVATION:
 *
 * geo elevation of where the media has been recorded or produced in meters
 * according to WGS84 (zero is average sea level) (double).
 */
#define GST_TAG_GEO_LOCATION_ELEVATION               "geo-location-elevation"
/**
 * GST_TAG_GEO_LOCATION_COUNTRY:
 *
 * The country (english name) where the media has been produced (string).
 */
#define GST_TAG_GEO_LOCATION_COUNTRY                 "geo-location-country"
/**
 * GST_TAG_GEO_LOCATION_CITY:
 *
 * The city (english name) where the media has been produced (string).
 */
#define GST_TAG_GEO_LOCATION_CITY                    "geo-location-city"
/**
 * GST_TAG_GEO_LOCATION_SUBLOCATION:
 *
 * A location 'smaller' than GST_TAG_GEO_LOCATION_CITY that specifies better
 * where the media has been produced. (e.g. the neighborhood) (string).
 *
 * This tag has been added as this is how it is handled/named in XMP's
 * Iptc4xmpcore schema.
 */
#define GST_TAG_GEO_LOCATION_SUBLOCATION             "geo-location-sublocation"
/**
 * GST_TAG_GEO_LOCATION_HORIZONTAL_ERROR:
 *
 * Represents the expected error on the horizontal positioning in
 * meters (double).
 */
#define GST_TAG_GEO_LOCATION_HORIZONTAL_ERROR   "geo-location-horizontal-error"
/**
 * GST_TAG_GEO_LOCATION_MOVEMENT_SPEED:
 *
 * Speed of the capturing device when performing the capture.
 * Represented in m/s. (double)
 *
 * See also #GST_TAG_GEO_LOCATION_MOVEMENT_DIRECTION
 */
#define GST_TAG_GEO_LOCATION_MOVEMENT_SPEED       "geo-location-movement-speed"
/**
 * GST_TAG_GEO_LOCATION_MOVEMENT_DIRECTION:
 *
 * Indicates the movement direction of the device performing the capture
 * of a media. It is represented as degrees in floating point representation,
 * 0 means the geographic north, and increases clockwise (double from 0 to 360)
 *
 * See also #GST_TAG_GEO_LOCATION_CAPTURE_DIRECTION
 */
#define GST_TAG_GEO_LOCATION_MOVEMENT_DIRECTION "geo-location-movement-direction"
/**
 * GST_TAG_GEO_LOCATION_CAPTURE_DIRECTION:
 *
 * Indicates the direction the device is pointing to when capturing
 * a media. It is represented as degrees in floating point representation,
 * 0 means the geographic north, and increases clockwise (double from 0 to 360)
 *
 * See also #GST_TAG_GEO_LOCATION_MOVEMENT_DIRECTION
 */
#define GST_TAG_GEO_LOCATION_CAPTURE_DIRECTION  "geo-location-capture-direction"
/**
 * GST_TAG_SHOW_NAME:
 *
 * Name of the show, used for displaying (string)
 */
#define GST_TAG_SHOW_NAME                         "show-name"
/**
 * GST_TAG_SHOW_SORTNAME:
 *
 * Name of the show, used for sorting (string)
 */
#define GST_TAG_SHOW_SORTNAME                     "show-sortname"
/**
 * GST_TAG_SHOW_EPISODE_NUMBER:
 *
 * Number of the episode within a season/show (unsigned integer)
 */
#define GST_TAG_SHOW_EPISODE_NUMBER               "show-episode-number"
/**
 * GST_TAG_SHOW_SEASON_NUMBER:
 *
 * Number of the season of a show/series (unsigned integer)
 */
#define GST_TAG_SHOW_SEASON_NUMBER                "show-season-number"
/**
 * GST_TAG_LYRICS:
 *
 * The lyrics of the media (string)
 */
#define GST_TAG_LYRICS                            "lyrics"
/**
 * GST_TAG_COMPOSER_SORTNAME:
 *
 * The composer's name, used for sorting (string)
 */
#define GST_TAG_COMPOSER_SORTNAME                 "composer-sortname"
/**
 * GST_TAG_GROUPING:
 *
 * Groups together media that are related and spans multiple tracks. An
 * example are multiple pieces of a concerto. (string)
 */
#define GST_TAG_GROUPING                          "grouping"
/**
 * GST_TAG_USER_RATING:
 *
 * Rating attributed by a person (likely the application user).
 * The higher the value, the more the user likes this media
 * (unsigned int from 0 to 100)
 */
#define GST_TAG_USER_RATING                       "user-rating"
/**
 * GST_TAG_DEVICE_MANUFACTURER:
 *
 * Manufacturer of the device used to create the media (string)
 */
#define GST_TAG_DEVICE_MANUFACTURER               "device-manufacturer"
/**
 * GST_TAG_DEVICE_MODEL:
 *
 * Model of the device used to create the media (string)
 */
#define GST_TAG_DEVICE_MODEL                      "device-model"
/**
 * GST_TAG_APPLICATION_NAME:
 *
 * Name of the application used to create the media (string)
 */
#define GST_TAG_APPLICATION_NAME                  "application-name"
/**
 * GST_TAG_APPLICATION_DATA:
 *
 * Arbitrary application data (sample)
 *
 * Some formats allow applications to add their own arbitrary data
 * into files. This data is application dependent.
 */
#define GST_TAG_APPLICATION_DATA          "application-data"
/**
 * GST_TAG_IMAGE_ORIENTATION:
 *
 * Represents the 'Orientation' tag from EXIF. Defines how the image
 * should be rotated and mirrored for display. (string)
 *
 * This tag has a predefined set of allowed values:
 *   "rotate-0"
 *   "rotate-90"
 *   "rotate-180"
 *   "rotate-270"
 *   "flip-rotate-0"
 *   "flip-rotate-90"
 *   "flip-rotate-180"
 *   "flip-rotate-270"
 *
 * The naming is adopted according to a possible transformation to perform
 * on the image to fix its orientation, obviously equivalent operations will
 * yield the same result.
 *
 * Rotations indicated by the values are in clockwise direction and
 * 'flip' means an horizontal mirroring.
 */
#define GST_TAG_IMAGE_ORIENTATION            "image-orientation"
/**
 * GST_TAG_PUBLISHER:
 *
 * Name of the label or publisher (string)
 *
 * Since: 1.2
 */
#define GST_TAG_PUBLISHER                         "publisher"
/**
 * GST_TAG_INTERPRETED_BY:
 *
 * Information about the people behind a remix and similar
 * interpretations of another existing piece (string)
 *
 * Since: 1.2
 */
#define GST_TAG_INTERPRETED_BY                    "interpreted-by"
/**
 * GST_TAG_MIDI_BASE_NOTE:
 *
 * [Midi note number](http://en.wikipedia.org/wiki/Note#Note_designation_in_accordance_with_octave_name)
 * of the audio track. This is useful for sample instruments and in particular
 * for multi-samples.
 *
 * Since: 1.4
 */
#define GST_TAG_MIDI_BASE_NOTE                    "midi-base-note"
/**
 * GST_TAG_PRIVATE_DATA:
 *
 * Any private data that may be contained in tags (sample).
 *
 * It is represented by #GstSample in which #GstBuffer contains the
 * binary data and the sample's info #GstStructure may contain any
 * extra information that identifies the origin or meaning of the data.
 *
 * Private frames in ID3v2 tags ('PRIV' frames) will be represented
 * using this tag, in which case the GstStructure will be named
 * "ID3PrivateFrame" and contain a field named "owner" of type string
 * which contains the owner-identification string from the tag.
 *
 * Since: 1.8
 */
#define GST_TAG_PRIVATE_DATA                         "private-data"

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstTagList, gst_tag_list_unref)

G_END_DECLS

#endif /* __GST_TAGLIST_H__ */
