/*
 *
 * Copyright © 2000 SuSE, Inc.
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of SuSE not be used in advertising or
 * publicity pertaining to distribution of the software without specific,
 * written prior permission.  SuSE makes no representations about the
 * suitability of this software for any purpose.  It is provided "as is"
 * without express or implied warranty.
 *
 * SuSE DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL SuSE
 * BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * Author:  Keith Packard, SuSE, Inc.
 */

/**
 * @file Xrender.h
 * @brief XRender library API.
 */

#ifndef _XRENDER_H_
#define _XRENDER_H_

#include <X11/Xfuncproto.h>
#include <X11/Xlib.h>
#include <X11/Xosdefs.h>
#include <X11/Xutil.h>

#include <X11/extensions/render.h>

/**
 * @mainpage libXrender API Documentation.
 *
 * Dummy text down here.
 */

/**
 * The direct component of a PictFormat.
 *
 * It contains a binary description of the color format used by the Picture.
 *
 * * A Zero bit alphaMask is declared to have an opaque alpha everywhere.
 * * A Zero bit redMask, greenMask and blueMask is declared to have red, green,
 *   blue == 0 everywhere.
 * * If any of redMask, greenMask or blueMask are zero, all other masks are
 *   zero.
 */
typedef struct {
    /** Red component binary displacement. */
    short red;
    /** Red component bit mask. */
    short redMask;
    /** Green component binary displacement. */
    short green;
    /** Green component bit mask. */
    short greenMask;
    /** Blue component binary displacement. */
    short blue;
    /** Blue component bit mask. */
    short blueMask;
    /** Alpha component binary displacement. */
    short alpha;
    /** Alpha component bit mask. */
    short alphaMask;
} XRenderDirectFormat;

/**
 * A Picture pixel format description.
 *
 * It describes the format used by the server to display colors.
 *
 * There are two types:
 * * Direct: Doesn't have a Colormap and the DirectFormat structure describes
 *   the pixel format.
 * * Indexed: Has a Colormap and it's DirectFormat structure is filled with
 *   zeros.
 */
typedef struct {
    /** XID of this structure server instance. */
    PictFormat id;
    /** Color management type. */
    int type;
    /** Pixel bit depth. */
    int depth;
    /** Color component description. */
    XRenderDirectFormat direct;
    /** XID of the map of indexed colors on the server. */
    Colormap colormap;
} XRenderPictFormat;

/*< XRenderPictFormat template field masks.
 * @{
 */
/** Include ID field. @hideinitializer */
#define PictFormatID (1 << 0)
/** Include Type field. @hideinitializer */
#define PictFormatType (1 << 1)
/** Include Depth field. @hideinitializer */
#define PictFormatDepth (1 << 2)

/*<--- XRenderPictFormat->direct fields. */
/** Include Direct->Red field. @hideinitializer */
#define PictFormatRed (1 << 3)
/** Include Direct->RedMask field. @hideinitializer */
#define PictFormatRedMask (1 << 4)
/** Include Direct->Green field. @hideinitializer */
#define PictFormatGreen (1 << 5)
/** Include Direct->GreenMask field. @hideinitializer */
#define PictFormatGreenMask (1 << 6)
/** Include Direct->Blue field. @hideinitializer */
#define PictFormatBlue (1 << 7)
/** Include Direct->BlueMask field. @hideinitializer */
#define PictFormatBlueMask (1 << 8)
/** Include Direct->Alpha field. @hideinitializer */
#define PictFormatAlpha (1 << 9)
/** Include Direct->AlphaMask field. @hideinitializer */
#define PictFormatAlphaMask (1 << 10)

/** Include Colormap field. @hideinitializer */
#define PictFormatColormap (1 << 11)
/** @} */

/**
 * Picture rendering attributes.
 */
typedef struct _XRenderPictureAttributes {
    /** How to repeat the picture. */
    int repeat;

    /** A replacement alpha-map. Must be a pixmap-containing Picture. */
    Picture alpha_map;
    /** Horizontal displacement of the replacement alpha-map. */
    int alpha_x_origin;
    /** Vertical displacement of the replacement alpha-map. */
    int alpha_y_origin;

    /** Horizontal displacement of the clip mask. */
    int clip_x_origin;
    /** Vertical displacement of the clip mask. */
    int clip_y_origin;
    /** A r/w restriction to the drawable. */
    Pixmap clip_mask;

    /** Whether to receive GraphicsExpose events. @note Ignored field. */
    Bool graphics_exposures;
    /** How to clip pixels on subwindow overlap. */
    int subwindow_mode;
    /** Alpha mask generation mode. */
    int poly_edge;
    /** Alpha value rasterization mode. */
    int poly_mode;
    /** Dithering mode. @note Ignored field. */
    Atom dither;
    /** Treat alpha channels independently. */
    Bool component_alpha;
} XRenderPictureAttributes;

/** An alpha-blended color with premultiplied components.
 *
 *  Values are in the range from 0 to 65535 inclusive, scaled down to the right
 *  hardware values by the server. Colors must be premultiplied by alpha by the
 *  client in all cases but gradient operations.
 */
typedef struct {
    /** Red color channel. */
    unsigned short red;
    /** Green color channel. */
    unsigned short green;
    /** Blue color channel. */
    unsigned short blue;
    /** Alpha color channel. */
    unsigned short alpha;
} XRenderColor;

/**
 * Glyph positioning and sizing information.
 *
 * A glyph is positioned by taking the requested position and substracting the
 * center offset.
 */
typedef struct _XGlyphInfo {
    /** Glyph width. */
    unsigned short width;
    /** Glyph height. */
    unsigned short height;

    /** Horizontal Glyph center offset relative to the upper-left corner. */
    short x;
    /** Vertical Glyph center offset relative to the upper-left corner. */
    short y;

    /** Horizontal margin to the next Glyph. */
    short xOff;
    /** Vertical margin to the next Glyph. */
    short yOff;
} XGlyphInfo;

/*< Glyph Elements.
 *  Group of glyphs to be rendered.
 *  While selecting the right element type, you should use as a reference the
 *  largest identifier in Elt->glyphset.
 */
/** @{ */

/**
 * 8-bit Glyph Element.
 */
typedef struct _XGlyphElt8 {
    /** Set of available glyphs. */
    GlyphSet glyphset;

    /** 8-bit glyph id array. */
    _Xconst char *chars;
    /** Glyph array size. */
    int nchars;

    /** Horizontal offset. */
    int xOff;
    /** Vertical offset. */
    int yOff;
} XGlyphElt8;

/**
 * 16-bit Glyph Element.
 */
typedef struct _XGlyphElt16 {
    /** Set of available glyphs. */
    GlyphSet glyphset;

    /** 16-bit glyph id array. */
    _Xconst unsigned short *chars;
    /** Glyph array size. */
    int nchars;

    /** Horizontal offset. */
    int xOff;
    /** Vertical offset. */
    int yOff;
} XGlyphElt16;

/**
 * 32-bit Glyph Element.
 */
typedef struct _XGlyphElt32 {
    /** Set of available glyphs. */
    GlyphSet glyphset;

    /** 32-bit glyph id array. */
    _Xconst unsigned int *chars;
    /** Glyph array size. */
    int nchars;

    /** Horizontal offset. */
    int xOff;
    /** Vertical offset. */
    int yOff;
} XGlyphElt32;
/**@} */

/*< Utility number types.
 *
 */
/**@{ */

/**
 * Floating-point number.
 */
typedef double XDouble;

/**
 * Fixed-point number.
 */
typedef int XFixed;

/** Turn XDouble into XFixed. @hideinitializer */
#define XDoubleToFixed(f) ((XFixed)((f)*65536))
/** Turn XFixed into XDouble. @hideinitializer */
#define XFixedToDouble(f) (((XDouble)(f)) / 65536)
/** @} */

/**
 * Point coordinates stored as floats.
 */
typedef struct _XPointDouble {
    XDouble x, y;
} XPointDouble;

/**
 * Point coordinates as integers.
 */
typedef struct _XPointFixed {
    XFixed x, y;
} XPointFixed;

/**
 * Line described by two points.
 */
typedef struct _XLineFixed {
    XPointFixed p1, p2;
} XLineFixed;

/**
 * Triangle described by it's vertices.
 * @see XTrap
 */
typedef struct _XTriangle {
    XPointFixed p1, p2, p3;
} XTriangle;

/**
 * Circle described by it's center point and a radius.
 */
typedef struct _XCircle {
    XFixed x;
    XFixed y;
    XFixed radius;
} XCircle;

/** A trapezoid.
 *
 *  @deprecated Use XTrap instead
 *  @see
 *    * XTriangle
 *    * XTrap
 */
typedef struct _XTrapezoid {
    XFixed top, bottom;
    XLineFixed left, right;
} XTrapezoid;

/**
 * A transform matrix.
 */
typedef struct _XTransform {
    XFixed matrix[3][3];
} XTransform;

/**
 * Group filters and filter aliases.
 */
typedef struct _XFilters {
    /** Filter names count. */
    int nfilter;
    /** Filter names array. */
    char **filter;
    /** Aliases array count. */
    int nalias;
    /** Array of «Index in .filter of the aliased filter or 0xffff». */
    short *alias;
} XFilters;

/**
 * The value of an indexed color.
 */
typedef struct _XIndexValue {
    /** Index ID. */
    unsigned long pixel;
    /** Color components. */
    unsigned short red, green, blue, alpha;
} XIndexValue;

/**
 * A single cursor frame.
 */
typedef struct _XAnimCursor {
    /** Existing cursor. */
    Cursor cursor;
    /** Animation delay. */
    unsigned long delay;
} XAnimCursor;

/**
 * An horizontal line.
 */
typedef struct _XSpanFix {
    XFixed left, right, y;
} XSpanFix;

/**
 * A trapezoid defined by two lines.
 * @see XTriangle
 */
typedef struct _XTrap {
    XSpanFix top, bottom;
} XTrap;

/**
 * Linear gradient shape.
 */
typedef struct _XLinearGradient {
    XPointFixed p1;
    XPointFixed p2;
} XLinearGradient;

/**
 * Radial gradient shape.
 */
typedef struct _XRadialGradient {
    XCircle inner;
    XCircle outer;
} XRadialGradient;

/**
 * Conical gradient shape.
 */
typedef struct _XConicalGradient {
    XPointFixed center;
    XFixed angle; /* in degrees */
} XConicalGradient;

_XFUNCPROTOBEGIN

/** @defgroup queries Early check queries.
 * @{
 */

/**
 * Ask for the Render extension presence and its base numbers.
 *
 * @param dpy Connection to the X server.
 * @param[out] event_basep first event number for the extension.
 * @param[out] error_basep first error number for the extension.
 * @return True if Render is present.
 */
Bool XRenderQueryExtension(Display *dpy, int *event_basep, int *error_basep);

/**
 * Ask for the extension version.
 *
 * @param dpy Connection to the X server.
 * @param[out] major_versionp Extension's major version.
 * @param[out] minor_versionp Extension's major version.
 * @return Status «1» on success.
 */
Status XRenderQueryVersion(Display *dpy, int *major_versionp,
                           int *minor_versionp);

/**
 * Check for and cache compatible picture formats.
 *
 * @param dpy Connection to the X server.
 * @return Status «1» on success.
 */
Status XRenderQueryFormats(Display *dpy);

/**
 * Ask for the current subpixel order of a screen.
 *
 * @param dpy Connection to the X server.
 * @param[in] screen Target screen number.
 * @return SubPixelUnknown on error, else a subpixel order.
 */
int XRenderQuerySubpixelOrder(Display *dpy, int screen);

/**
 * Change the subpixel order of a screen.
 *
 * @param dpy Connection to the X server
 * @param[in] screen Target screen number.
 * @param[in] subpixel Requested subpixel order.
 * @return True if the operation was successful.
 */
Bool XRenderSetSubpixelOrder(Display *dpy, int screen, int subpixel);
/** @} */

/**
 * Ask for the Picture format for a Visual.
 *
 * @param dpy Connection to the X server.
 * @param[in] visual Reference Visual object.
 * @return The requested Picture format.
 */
XRenderPictFormat *XRenderFindVisualFormat(Display *dpy,
                                           _Xconst Visual *visual);

/**
 * Ask for matching Picture formats from a template.
 *
 * @param dpy Connection to the X server.
 * @param[in] mask `templ` fields mask to use.
 * @param[in] templ Requested Picture format template.
 * @param[in] count Skip `count` formats.
 * @return NULL if no matching format found, else a Picture format.
 */
XRenderPictFormat *XRenderFindFormat(Display *dpy, unsigned long mask,
                                     _Xconst XRenderPictFormat *templ,
                                     int count);

/** Standard format specifiers.
 * @{
 */
/** 8-bit RGB with Alpha. @hideinitializer */
#define PictStandardARGB32 0
/** 8-bit RGB. @hideinitializer */
#define PictStandardRGB24 1
/** 8-bit Alpha map. @hideinitializer */
#define PictStandardA8 2
/** 4-bit Alpha map. @hideinitializer */
#define PictStandardA4 3
/** 1-bit Alpha map. @hideinitializer */
#define PictStandardA1 4
/** Supported standard formats count. @hideinitializer */
#define PictStandardNUM 5
/** @} */

/**
 * Ask for a predefined standard picture format.
 *
 * This is a shorthand to XRenderFindFormat for finding common formats.
 *
 * @param dpy Connection to the X server.
 * @param[in] format Desired format specifier.
 * @return NULL if no matching format found, else a Picture format.
 */
XRenderPictFormat *XRenderFindStandardFormat(Display *dpy, int format);

/**
 * Ask for the indexed colors of a Picture format.
 *
 * @param dpy Connection to the X server.
 * @param[in] format Queried picture format.
 * @param[out] num Size of the output array.
 * @return An array of XIndexValue.
 */
XIndexValue *XRenderQueryPictIndexValues(Display *dpy,
                                         _Xconst XRenderPictFormat *format,
                                         int *num);

/**
 * Creates a Picture for a drawable.
 *
 * @param dpy Connection to the X server.
 * @param[in] drawable Target Drawable.
 * @param[in] format Format for the Picture.
 * @param[in] valuemask `attributes` fields mask to use.
 * @param[in] attributes Desired attributes for the Picture.
 * @return A Picture tied to the drawable.
 */
Picture XRenderCreatePicture(Display *dpy, Drawable drawable,
                             _Xconst XRenderPictFormat *format,
                             unsigned long valuemask,
                             _Xconst XRenderPictureAttributes *attributes);

/**
 * Free allocated structures for a Picture.
 *
 * @warning A freed Picture shouldn't be used again.
 *
 * @param dpy Connection to the X server.
 * @param[in] picture Target Picture.
 */
void XRenderFreePicture(Display *dpy, Picture picture);

/**
 * Change a Picture's attributes structure.
 *
 * @param dpy Connection to the X server.
 * @param[in] picture Target Picture.
 * @param[in] valuemask `attributes` fields mask to use.
 * @param[in] attributes Desired attributes for the Picture.
 */
void XRenderChangePicture(Display *dpy, Picture picture,
                          unsigned long valuemask,
                          _Xconst XRenderPictureAttributes *attributes);

/**
 * Change a Picture's clip mask to the specified rectangles.
 *
 * @param dpy Connection to the X server.
 * @param[in] picture Target Picture.
 * @param[in] xOrigin Horizontal mask origin.
 * @param[in] yOrigin Vertical mask origin.
 * @param[in] rects Array of rectangles to clip with.
 * @param[in] n `rects` array size.
 */
void XRenderSetPictureClipRectangles(Display *dpy, Picture picture, int xOrigin,
                                     int yOrigin, _Xconst XRectangle *rects,
                                     int n);

/**
 * Change a Picture's clip mask to the specified Region.
 *
 * @param dpy Connection to the X server.
 * @param[in] picture Target Picture.
 * @param[in] r Region to clip with.
 */
void XRenderSetPictureClipRegion(Display *dpy, Picture picture, Region r);

/**
 * Change a Picture's Transform matrix.
 *
 * @param dpy Connection to the X server
 * @param[in] picture Target Picture.
 * @param[in] transform Transform matrix to use.
 */
void XRenderSetPictureTransform(Display *dpy, Picture picture,
                                XTransform *transform);

/**
 * Combines two Pictures with the specified compositing operation.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] mask Composition mask.
 * @param[in] dst Picture to combine into.
 * @param[in] src_x Horizontal `src` origin offset.
 * @param[in] src_y Vertical `src` origin offset
 * @param[in] mask_x Horizontal `mask` origin offset.
 * @param[in] mask_y Vertical `mask` origin offset.
 * @param[in] dst_x Horizontal `dst` origin offset.
 * @param[in] dst_y Vertical `dst` origin offset.
 * @param[in] width Maximum composition width.
 * @param[in] height Maximum composition height.
 */
void XRenderComposite(Display *dpy, int op, Picture src, Picture mask,
                      Picture dst, int src_x, int src_y, int mask_x, int mask_y,
                      int dst_x, int dst_y, unsigned int width,
                      unsigned int height);

/**
 * Create a Glyph Set.
 *
 * @param dpy Connection to the X server.
 * @param[in] format Desired format for the Glyphs Picture.
 * @return A GlyphSet.
 */
GlyphSet XRenderCreateGlyphSet(Display *dpy, _Xconst XRenderPictFormat *format);

/**
 * Generate a new reference for an existing Glyph Set.
 *
 * @param dpy Connection to the X server.
 * @param[in] existing Target Glyph Set.
 * @return A GlyphSet identical to `existing`.
 */
GlyphSet XRenderReferenceGlyphSet(Display *dpy, GlyphSet existing);

/**
 * Free allocated structures for a GlyphSet.
 *
 * If there's more references to the underlying GlyphSet structures, this will
 * remove only the specified GlyphSet reference.
 *
 * @warning A freed GlyphSet shouldn't be used again.
 *
 * @param dpy Connection to the X server.
 * @param[in] glyphset Target GlyphSet.
 */
void XRenderFreeGlyphSet(Display *dpy, GlyphSet glyphset);

/**
 * Add new Glyphs to a GlyphSet.
 *
 * @param dpy Connection to the X server.
 * @param[in] glyphset Glyph storage destination.
 * @param[in] gids Array of ids for the new Glyphs.
 * @param[in] glyphs Array of new Glyphs info.
 * @param[in] nglyphs Number of Glyphs to add.
 * @param[in] images Byte array containing the Glyphs graphics.
 * @param[in] nbyte_images Size of the `images` byte array.
 */
void XRenderAddGlyphs(Display *dpy, GlyphSet glyphset, _Xconst Glyph *gids,
                      _Xconst XGlyphInfo *glyphs, int nglyphs,
                      _Xconst char *images, int nbyte_images);

/**
 * Free allocated Glyphs.
 *
 * @param dpy Connection to the X server.
 * @param[in] glyphset GlyphSet storing the Glyphs.
 * @param[in] gids Identifier array of the Glyphs to dellocate.
 * @param[in] nglyphs Glyph count.
 */
void XRenderFreeGlyphs(Display *dpy, GlyphSet glyphset, _Xconst Glyph *gids,
                       int nglyphs);

/**
 * Draw a 8-bit character string into a Picture.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] dst Picture to combine into.
 * @param[in] maskFormat Picture format of the generated Picture mask.
 * @param[in] glyphset Glyph Source.
 * @param[in] xSrc Horizontal `src` origin offset.
 * @param[in] ySrc Vertical `src` origin offset.
 * @param[in] xDst Horizontal `dst` origin offset.
 * @param[in] yDst Vertical `dst` origin offset.
 * @param[in] string String to clip to.
 * @param[in] nchar String length.
 */
void XRenderCompositeString8(Display *dpy, int op, Picture src, Picture dst,
                             _Xconst XRenderPictFormat *maskFormat,
                             GlyphSet glyphset, int xSrc, int ySrc, int xDst,
                             int yDst, _Xconst char *string, int nchar);

/**
 * Draw a 16-bit character string into a Picture.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] dst Picture to combine into.
 * @param[in] maskFormat Picture format of the generated Picture mask.
 * @param[in] glyphset Glyph Source.
 * @param[in] xSrc Horizontal `src` origin offset.
 * @param[in] ySrc Vertical `src` origin offset.
 * @param[in] xDst Horizontal `dst` origin offset.
 * @param[in] yDst Vertical `dst` origin offset.
 * @param[in] string String to clip to.
 * @param[in] nchar String length.
 */
void XRenderCompositeString16(Display *dpy, int op, Picture src, Picture dst,
                              _Xconst XRenderPictFormat *maskFormat,
                              GlyphSet glyphset, int xSrc, int ySrc, int xDst,
                              int yDst, _Xconst unsigned short *string,
                              int nchar);

/**
 * Draw a 32-bit character string into a Picture.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] dst Picture to combine into.
 * @param[in] maskFormat Picture format of the generated Picture mask.
 * @param[in] glyphset Glyph Source.
 * @param[in] xSrc Horizontal `src` origin offset.
 * @param[in] ySrc Vertical `src` origin offset.
 * @param[in] xDst Horizontal `dst` origin offset.
 * @param[in] yDst Vertical `dst` origin offset.
 * @param[in] string String to clip to.
 * @param[in] nchar String length.
 */
void XRenderCompositeString32(Display *dpy, int op, Picture src, Picture dst,
                              _Xconst XRenderPictFormat *maskFormat,
                              GlyphSet glyphset, int xSrc, int ySrc, int xDst,
                              int yDst, _Xconst unsigned int *string,
                              int nchar);

/**
 * Draw several 8-bit Glyph Elements into a Picture.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] dst Picture to combine into.
 * @param[in] maskFormat Picture format of the generated Picture mask.
 * @param[in] xSrc Horizontal `src` origin offset.
 * @param[in] ySrc Vertical `src` origin offset.
 * @param[in] xDst Horizontal `dst` origin offset.
 * @param[in] yDst Vertical `dst` origin offset.
 * @param[in] elts Glyph Elements array to clip with.
 * @param[in] nelt Glyph Elements array size.
 */
void XRenderCompositeText8(Display *dpy, int op, Picture src, Picture dst,
                           _Xconst XRenderPictFormat *maskFormat, int xSrc,
                           int ySrc, int xDst, int yDst,
                           _Xconst XGlyphElt8 *elts, int nelt);

/**
 * Draw several 16-bit Glyph Elements into a Picture.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] dst Picture to combine into.
 * @param[in] maskFormat Picture format of the generated Picture mask.
 * @param[in] xSrc Horizontal `src` origin offset.
 * @param[in] ySrc Vertical `src` origin offset.
 * @param[in] xDst Horizontal `dst` origin offset.
 * @param[in] yDst Vertical `dst` origin offset.
 * @param[in] elts Glyph Elements array to clip with.
 * @param[in] nelt Glyph Elements array size.
 */
void XRenderCompositeText16(Display *dpy, int op, Picture src, Picture dst,
                            _Xconst XRenderPictFormat *maskFormat, int xSrc,
                            int ySrc, int xDst, int yDst,
                            _Xconst XGlyphElt16 *elts, int nelt);

/**
 * Draw several 32-bit Glyph Elements into a Picture.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] dst Picture to combine into.
 * @param[in] maskFormat Picture format of the generated Picture mask.
 * @param[in] xSrc Horizontal `src` origin offset.
 * @param[in] ySrc Vertical `src` origin offset.
 * @param[in] xDst Horizontal `dst` origin offset.
 * @param[in] yDst Vertical `dst` origin offset.
 * @param[in] elts Glyph Elements to clip with.
 * @param[in] nelt Glyph Elements array size.
 */
void XRenderCompositeText32(Display *dpy, int op, Picture src, Picture dst,
                            _Xconst XRenderPictFormat *maskFormat, int xSrc,
                            int ySrc, int xDst, int yDst,
                            _Xconst XGlyphElt32 *elts, int nelt);

/**
 * Fill a Rectangle with the given color.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] dst Picture to draw into.
 * @param[in] color Color to fill with.
 * @param[in] x Horizontal offset.
 * @param[in] y Vertical offset.
 * @param[in] width Rectangle width.
 * @param[in] height Rectangle height.
 */
void XRenderFillRectangle(Display *dpy, int op, Picture dst,
                          _Xconst XRenderColor *color, int x, int y,
                          unsigned int width, unsigned int height);

/**
 * Fill a bunch of Rectangle with the given color.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] dst Picture to draw into.
 * @param[in] color Color to fill with.
 * @param[in] rectangles Array of Rectangles to fill.
 * @param[in] n_rects `rectangles` array size.
 */
void XRenderFillRectangles(Display *dpy, int op, Picture dst,
                           _Xconst XRenderColor *color,
                           _Xconst XRectangle *rectangles, int n_rects);

/**
 * Combine two Pictures using a bunch of Trapezoids as the mask.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] dst Picture to combine into.
 * @param[in] maskFormat Picture format of the generated Picture mask.
 * @param[in] xSrc Horizontal `src` origin offset.
 * @param[in] ySrc Vertical `src` origin offset.
 * @param[in] traps Array of Trapezoids to clip with.
 * @param[in] ntrap `traps` Array size.
 */
void XRenderCompositeTrapezoids(Display *dpy, int op, Picture src, Picture dst,
                                _Xconst XRenderPictFormat *maskFormat, int xSrc,
                                int ySrc, _Xconst XTrapezoid *traps, int ntrap);

/**
 * Combine two Pictures using a bunch of Triangles as the mask.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] dst Picture to combine into.
 * @param[in] maskFormat Picture format of the generated Picture mask.
 * @param[in] xSrc Horizontal `src` origin offset.
 * @param[in] ySrc Vertical `src` origin offset.
 * @param[in] triangles Array of Triangles to clip with.
 * @param[in] ntriangle `triangles` array size.
 */
void XRenderCompositeTriangles(Display *dpy, int op, Picture src, Picture dst,
                               _Xconst XRenderPictFormat *maskFormat, int xSrc,
                               int ySrc, _Xconst XTriangle *triangles,
                               int ntriangle);

/**
 * Combine two Pictures using a Triangle Strip as the mask.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] dst Picture to combine into.
 * @param[in] maskFormat Picture format of the generated Picture mask.
 * @param[in] xSrc Horizontal `src` origin offset.
 * @param[in] ySrc Vertical `src` origin offset.
 * @param[in] points Array of Points to create Triangles with.
 * @param[in] npoint `points` array size.
 */
void XRenderCompositeTriStrip(Display *dpy, int op, Picture src, Picture dst,
                              _Xconst XRenderPictFormat *maskFormat, int xSrc,
                              int ySrc, _Xconst XPointFixed *points,
                              int npoint);

/**
 * Combine two Pictures using a Triangle Fan as the mask.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] dst Picture to combine into.
 * @param[in] maskFormat Picture format of the generated Picture mask.
 * @param[in] xSrc Horizontal `src` origin offset.
 * @param[in] ySrc Vertical `src` origin offset.
 * @param[in] points Array of Points to create Triangles with.
 * @param[in] npoint `points` array size.
 */
void XRenderCompositeTriFan(Display *dpy, int op, Picture src, Picture dst,
                            _Xconst XRenderPictFormat *maskFormat, int xSrc,
                            int ySrc, _Xconst XPointFixed *points, int npoint);

/**
 * Combine two Pictures using a Polygon as the mask.
 *
 * @param dpy Connection to the X server.
 * @param[in] op Compositing operation to perform.
 * @param[in] src Picture to combine with.
 * @param[in] dst Picture to combine into.
 * @param[in] maskFormat Picture format of the generated Picture mask.
 * @param[in] xSrc Horizontal `src` origin offset.
 * @param[in] ySrc Vertical `src` origin offset.
 * @param[in] xDst Horizontal `dst` origin offset.
 * @param[in] yDst Vertical `dst` origin offset.
 * @param[in] fpoints Array of DoublePoints to create a Polygon with.
 * @param[in] npoints `points` array size.
 * @param winding Unused.
 */
void XRenderCompositeDoublePoly(Display *dpy, int op, Picture src, Picture dst,
                                _Xconst XRenderPictFormat *maskFormat, int xSrc,
                                int ySrc, int xDst, int yDst,
                                _Xconst XPointDouble *fpoints, int npoints,
                                int winding);

/**
 * Parse a color string.
 *
 * @param dpy Connection to the X server.
 * @param[in] spec Null-terminated string.
 * @param[out] def Parsing result.
 * @return Status «1» on success.
 */
Status XRenderParseColor(Display *dpy, char *spec, XRenderColor *def);

/**
 * Creates a cursor looking like a Picture.
 *
 * @param dpy Connection to the X server.
 * @param[in] source Picture defining the cursor look.
 * @param[in] x Horizontal offset.
 * @param[in] y Vertical offset.
 * @return A Cursor.
 */
Cursor XRenderCreateCursor(Display *dpy, Picture source, unsigned int x,
                           unsigned int y);

/**
 * Ask for Filters applicable to some Drawable.
 *
 * @param dpy Connection to the X server.
 * @param[in] drawable Target Drawable.
 * @return Available Filters and Aliases.
 */
XFilters *XRenderQueryFilters(Display *dpy, Drawable drawable);

/**
 * Set the current filter of a Picture.
 *
 * @note On Picture creation, the «Nearest» filter is set by default.
 *
 * @param dpy Connection to the X server.
 * @param[in] picture Target.
 * @param[in] filter Filter name.
 * @param[in] params Filter parameters array.
 * @param[in] nparams `params` array size.
 */
void XRenderSetPictureFilter(Display *dpy, Picture picture, const char *filter,
                             XFixed *params, int nparams);

/**
 * Create an animated Cursor from the given Cursor frames.
 *
 * @param dpy Connection to the X server.
 * @param[in] ncursor Cursor frames count.
 * @param[in] cursors Cursor frames array.
 * @return An animated Cursor.
 */
Cursor XRenderCreateAnimCursor(Display *dpy, int ncursor, XAnimCursor *cursors);

/**
 * Add the given Trapezoids to a single-channel Picture.
 *
 * @param dpy Connection to the X server.
 * @param[in] picture An alpha-only Picture.
 * @param[in] xOff Horizontal offset.
 * @param[in] yOff Vertical offset.
 * @param[in] traps Array of trapezoids.
 * @param[in] ntrap `traps` array size.
 */
void XRenderAddTraps(Display *dpy, Picture picture, int xOff, int yOff,
                     _Xconst XTrap *traps, int ntrap);

/**
 * Create a Picture filled with a single Color.
 *
 * @param dpy Connection to the X server.
 * @param[in] color Desired filling.
 * @return A single Color Picture.
 */
Picture XRenderCreateSolidFill(Display *dpy, const XRenderColor *color);

/**
 * Create a Picture filled with a Linear Gradient.
 *
 * @param dpy Connection to the X server.
 * @param[in] gradient Gradient geometry.
 * @param[in] stops Stop sections.
 * @param[in] colors Stop colors.
 * @param[in] nstops Stops count.
 * @return A Picture filled with a Linear Gradient.
 */
Picture XRenderCreateLinearGradient(Display *dpy,
                                    const XLinearGradient *gradient,
                                    const XFixed *stops,
                                    const XRenderColor *colors, int nstops);

/**
 * Create a Picture filled with a Radial Gradient.
 *
 * @param dpy Connection to the X server.
 * @param[in] gradient Gradient geometry.
 * @param[in] stops Stop sections.
 * @param[in] colors Stop colors.
 * @param[in] nstops Stops count.
 * @return A Picture filled with a Radial Gradient.
 */
Picture XRenderCreateRadialGradient(Display *dpy,
                                    const XRadialGradient *gradient,
                                    const XFixed *stops,
                                    const XRenderColor *colors, int nstops);

/**
 * Create a Picture filled with a Conical Gradient.
 *
 * @param dpy Connection to the X server.
 * @param[in] gradient Gradient geometry.
 * @param[in] stops Stop sections.
 * @param[in] colors Stop colors.
 * @param[in] nstops Stops count.
 * @return A Picture filled with a Conical Gradient.
 */
Picture XRenderCreateConicalGradient(Display *dpy,
                                     const XConicalGradient *gradient,
                                     const XFixed *stops,
                                     const XRenderColor *colors, int nstops);

_XFUNCPROTOEND

#endif /* _XRENDER_H_ */
