/* GStreamer
 * Copyright (C) <2013> Wim Taymans <wim.taymans@gmail.com>
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

#ifndef __GST_VIDEO_TILE_H__
#define __GST_VIDEO_TILE_H__

#include <gst/gst.h>
#include <gst/video/video-prelude.h>

G_BEGIN_DECLS

/**
 * GstVideoTileType:
 * @GST_VIDEO_TILE_TYPE_INDEXED: Tiles are indexed. Use
 *   gst_video_tile_get_index () to retrieve the tile at the requested
 *   coordinates.
 *
 * Enum value describing the most common tiling types.
 */
typedef enum
{
  GST_VIDEO_TILE_TYPE_INDEXED = 0
} GstVideoTileType;

#define GST_VIDEO_TILE_TYPE_SHIFT     (16)

/**
 * GST_VIDEO_TILE_TYPE_MASK: (value 65535)
 */
#define GST_VIDEO_TILE_TYPE_MASK      ((1 << GST_VIDEO_TILE_TYPE_SHIFT) - 1)

/**
 * GST_VIDEO_TILE_MAKE_MODE:
 * @num: the mode number to create
 * @type: the tile mode type
 *
 * use this macro to create new tile modes.
 */
#define GST_VIDEO_TILE_MAKE_MODE(num, type) \
    (((num) << GST_VIDEO_TILE_TYPE_SHIFT) | (GST_VIDEO_TILE_TYPE_ ##type))

/**
 * GST_VIDEO_TILE_MODE_TYPE:
 * @mode: the tile mode
 *
 * Get the tile mode type of @mode
 */
#define GST_VIDEO_TILE_MODE_TYPE(mode)       ((mode) & GST_VIDEO_TILE_TYPE_MASK)

/**
 * GST_VIDEO_TILE_MODE_IS_INDEXED:
 * @mode: a tile mode
 *
 * Check if @mode is an indexed tile type
 */
#define GST_VIDEO_TILE_MODE_IS_INDEXED(mode) (GST_VIDEO_TILE_MODE_TYPE(mode) == GST_VIDEO_TILE_TYPE_INDEXED)


#define GST_VIDEO_TILE_Y_TILES_SHIFT     (16)

/**
 * GST_VIDEO_TILE_X_TILES_MASK: (value 65535)
 */
#define GST_VIDEO_TILE_X_TILES_MASK      ((1 << GST_VIDEO_TILE_Y_TILES_SHIFT) - 1)

/**
 * GST_VIDEO_TILE_MAKE_STRIDE:
 * @x_tiles: number of tiles in X
 * @y_tiles: number of tiles in Y
 *
 * Encode the number of tile in X and Y into the stride.
 */
#define GST_VIDEO_TILE_MAKE_STRIDE(x_tiles, y_tiles) \
    (((y_tiles) << GST_VIDEO_TILE_Y_TILES_SHIFT) | (x_tiles))

/**
 * GST_VIDEO_TILE_X_TILES:
 * @stride: plane stride
 *
 * Extract the number of tiles in X from the stride value.
 */
#define GST_VIDEO_TILE_X_TILES(stride) ((stride) & GST_VIDEO_TILE_X_TILES_MASK)

/**
 * GST_VIDEO_TILE_Y_TILES:
 * @stride: plane stride
 *
 * Extract the number of tiles in Y from the stride value.
 */
#define GST_VIDEO_TILE_Y_TILES(stride) ((stride) >> GST_VIDEO_TILE_Y_TILES_SHIFT)

typedef struct _GstVideoTileInfo GstVideoTileInfo;

/**
 * GstVideoTileMode:
 * @GST_VIDEO_TILE_MODE_UNKNOWN: Unknown or unset tile mode
 * @GST_VIDEO_TILE_MODE_ZFLIPZ_2X2: Every four adjacent blocks - two
 *    horizontally and two vertically are grouped together and are located
 *    in memory in Z or flipped Z order. In case of odd rows, the last row
 *    of blocks is arranged in linear order.
 * @GST_VIDEO_TILE_MODE_LINEAR: Tiles are in row order. (Since: 1.18)
 * @GST_VIDEO_TILE_MODE_LINEAR_SUBSAMPLED: Tiles are in row order, with
 *   variable tile size according to subsampling. (Since: 1.20)
 *
 * Enum value describing the available tiling modes.
 */
typedef enum
{
  GST_VIDEO_TILE_MODE_UNKNOWN = 0,
  GST_VIDEO_TILE_MODE_ZFLIPZ_2X2 = GST_VIDEO_TILE_MAKE_MODE (1, INDEXED),
  /**
   * GST_VIDEO_TILE_MODE_LINEAR:
   *
   * Tiles are in row order.
   *
   * Since: 1.18
   */
  GST_VIDEO_TILE_MODE_LINEAR = GST_VIDEO_TILE_MAKE_MODE (2, INDEXED),
} GstVideoTileMode;


/**
 * GstVideoTileInfo:
 *
 * Description of a tile. This structure allow to describe arbitrary tile
 * dimensions and sizes.
 *
 * Since: 1.22
 */
struct _GstVideoTileInfo
{
  /**
   * GstVideoTileInfo.width:
   *
   * The width in pixels of a tile. This value can be zero if the number of
   * pixels per line is not an integer value.
   *
   * Since: 1.22
   */
  guint width;

  /**
   * GstVideoTileInfo::height:
   *
   * The width in pixels of a tile. This value can be zero if the number of
   * pixels per line is not an integer value.
   *
   * Since: 1.22
   */
  guint height;

  /**
   * GstVideoTileInfo.stride:
   *
   * The stride (in bytes) of a tile line. Regardless if the tile have sub-tiles
   * this stride multiplied by the height should be equal to
   * #GstVideoTileInfo.size. This value is used to translate into linear stride
   * when older APIs are being used to expose this format.
   *
   * Since: 1.22
   */
  guint stride;

  /**
   * GstVideoTileInfo.size:
   *
   * The size in bytes of a tile. This value must be divisible by
   * #GstVideoTileInfo.stride.
   *
   * Since: 1.22
   */
  guint size;

  /* <private> */
  guint32 padding[GST_PADDING];
};

GST_VIDEO_API
guint           gst_video_tile_get_index                (GstVideoTileMode mode, gint x, gint y,
                                                         gint x_tiles, gint y_tiles);


G_END_DECLS

#endif /* __GST_VIDEO_TILE_H__ */
