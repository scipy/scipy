/*
 * $XFree86: xc/include/extensions/renderproto.h,v 1.12 2002/09/26 02:56:48 keithp Exp $
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

#ifndef _XRENDERP_H_
#define _XRENDERP_H_

#include <X11/Xmd.h>
#include <X11/extensions/render.h>

#define Window CARD32
#define Drawable CARD32
#define Font CARD32
#define Pixmap CARD32
#define Cursor CARD32
#define Colormap CARD32
#define GContext CARD32
#define Atom CARD32
#define VisualID CARD32
#define Time CARD32
#define KeyCode CARD8
#define KeySym CARD32

#define Picture	    CARD32
#define PictFormat  CARD32
#define Fixed	    INT32
#define Glyphset    CARD32

/*
 * data structures
 */

typedef struct {
    CARD16  red B16;
    CARD16  redMask B16;
    CARD16  green B16;
    CARD16  greenMask B16;
    CARD16  blue B16;
    CARD16  blueMask B16;
    CARD16  alpha B16;
    CARD16  alphaMask B16;
} xDirectFormat;

#define sz_xDirectFormat    16

typedef struct {
    PictFormat	id B32;
    CARD8	type;
    CARD8	depth;
    CARD16	pad1 B16;
    xDirectFormat   direct;
    Colormap	colormap;
} xPictFormInfo;

#define sz_xPictFormInfo    28

typedef struct {
    VisualID	visual;
    PictFormat	format;
} xPictVisual;

#define sz_xPictVisual	    8

typedef struct {
    CARD8	depth;
    CARD8	pad1;
    CARD16	nPictVisuals B16;
    CARD32	pad2 B32;
} xPictDepth;

#define sz_xPictDepth	8

typedef struct {
    CARD32	nDepth B32;
    PictFormat	fallback B32;
} xPictScreen;

#define sz_xPictScreen	8

typedef struct {
    CARD32	pixel B32;
    CARD16	red B16;
    CARD16	green B16;
    CARD16	blue B16;
    CARD16	alpha B16;
} xIndexValue;

#define sz_xIndexValue	12

typedef struct {
    CARD16	red B16;
    CARD16	green B16;
    CARD16	blue B16;
    CARD16	alpha B16;
} xRenderColor;

#define sz_xRenderColor	8

typedef struct {
    Fixed	x B32;
    Fixed	y B32;
} xPointFixed;

#define sz_xPointFixed	8

typedef struct {
    xPointFixed	p1;
    xPointFixed p2;
} xLineFixed;

#define sz_xLineFixed	16

typedef struct {
    xPointFixed	p1, p2, p3;
} xTriangle;

#define sz_xTriangle	24

typedef struct {
    Fixed	top B32;
    Fixed	bottom B32;
    xLineFixed	left;
    xLineFixed	right;
} xTrapezoid;

#define sz_xTrapezoid	40

typedef struct {
    CARD16  width B16;
    CARD16  height B16;
    INT16   x B16;
    INT16   y B16;
    INT16   xOff B16;
    INT16   yOff B16;
} xGlyphInfo;

#define sz_xGlyphInfo	12

typedef struct {
    CARD8   len;
    CARD8   pad1;
    CARD16  pad2;
    INT16   deltax;
    INT16   deltay;
} xGlyphElt;

#define sz_xGlyphElt	8

typedef struct {
    Fixed   l, r, y;
} xSpanFix;

#define sz_xSpanFix	12

typedef struct {
    xSpanFix	top, bot;
} xTrap;

#define sz_xTrap	24

/*
 * requests and replies
 */
typedef struct {
    CARD8   reqType;
    CARD8   renderReqType;
    CARD16  length B16;
    CARD32  majorVersion B32;
    CARD32  minorVersion B32;
} xRenderQueryVersionReq;

#define sz_xRenderQueryVersionReq   12

typedef struct {
    BYTE    type;   /* X_Reply */
    BYTE    pad1;
    CARD16  sequenceNumber B16;
    CARD32  length B32;
    CARD32  majorVersion B32;
    CARD32  minorVersion B32;
    CARD32  pad2 B32;
    CARD32  pad3 B32;
    CARD32  pad4 B32;
    CARD32  pad5 B32;
} xRenderQueryVersionReply;

#define sz_xRenderQueryVersionReply	32

typedef struct {
    CARD8   reqType;
    CARD8   renderReqType;
    CARD16  length B16;
} xRenderQueryPictFormatsReq;

#define sz_xRenderQueryPictFormatsReq	4

typedef struct {
    BYTE    type;   /* X_Reply */
    BYTE    pad1;
    CARD16  sequenceNumber B16;
    CARD32  length B32;
    CARD32  numFormats B32;
    CARD32  numScreens B32;
    CARD32  numDepths B32;
    CARD32  numVisuals B32;
    CARD32  numSubpixel B32;	    /* Version 0.6 */
    CARD32  pad5 B32;
} xRenderQueryPictFormatsReply;

#define sz_xRenderQueryPictFormatsReply	32

typedef struct {
    CARD8   reqType;
    CARD8   renderReqType;
    CARD16  length B16;
    PictFormat	format B32;
} xRenderQueryPictIndexValuesReq;

#define sz_xRenderQueryPictIndexValuesReq   8

typedef struct {
    BYTE    type;   /* X_Reply */
    BYTE    pad1;
    CARD16  sequenceNumber B16;
    CARD32  length B32;
    CARD32  numIndexValues;
    CARD32  pad2 B32;
    CARD32  pad3 B32;
    CARD32  pad4 B32;
    CARD32  pad5 B32;
    CARD32  pad6 B32;
} xRenderQueryPictIndexValuesReply;

#define sz_xRenderQueryPictIndexValuesReply 32

typedef struct {
    CARD8	reqType;
    CARD8	renderReqType;
    CARD16	length B16;
    Picture	pid B32;
    Drawable	drawable B32;
    PictFormat	format B32;
    CARD32	mask B32;
} xRenderCreatePictureReq;

#define sz_xRenderCreatePictureReq	    20

typedef struct {
    CARD8	reqType;
    CARD8	renderReqType;
    CARD16	length B16;
    Picture	picture B32;
    CARD32	mask B32;
} xRenderChangePictureReq;

#define sz_xRenderChangePictureReq	    12

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    Picture     picture B32;
    INT16	xOrigin B16;
    INT16	yOrigin B16;
} xRenderSetPictureClipRectanglesReq;

#define sz_xRenderSetPictureClipRectanglesReq	    12

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    Picture     picture B32;
} xRenderFreePictureReq;

#define sz_xRenderFreePictureReq	    8

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    CARD8	op;
    CARD8	pad1;
    CARD16	pad2 B16;
    Picture	src B32;
    Picture	mask B32;
    Picture	dst B32;
    INT16	xSrc B16;
    INT16	ySrc B16;
    INT16	xMask B16;
    INT16	yMask B16;
    INT16	xDst B16;
    INT16	yDst B16;
    CARD16	width B16;
    CARD16	height B16;
} xRenderCompositeReq;

#define sz_xRenderCompositeReq		    36

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    Picture	src B32;
    Picture	dst B32;
    CARD32	colorScale B32;
    CARD32	alphaScale B32;
    INT16	xSrc B16;
    INT16	ySrc B16;
    INT16	xDst B16;
    INT16	yDst B16;
    CARD16	width B16;
    CARD16	height B16;
} xRenderScaleReq;

#define sz_xRenderScaleReq			    32

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    CARD8	op;
    CARD8	pad1;
    CARD16	pad2 B16;
    Picture	src B32;
    Picture	dst B32;
    PictFormat	maskFormat B32;
    INT16	xSrc B16;
    INT16	ySrc B16;
} xRenderTrapezoidsReq;

#define sz_xRenderTrapezoidsReq			    24

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    CARD8	op;
    CARD8	pad1;
    CARD16	pad2 B16;
    Picture	src B32;
    Picture	dst B32;
    PictFormat	maskFormat B32;
    INT16	xSrc B16;
    INT16	ySrc B16;
} xRenderTrianglesReq;

#define sz_xRenderTrianglesReq			    24

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    CARD8	op;
    CARD8	pad1;
    CARD16	pad2 B16;
    Picture	src B32;
    Picture	dst B32;
    PictFormat	maskFormat B32;
    INT16	xSrc B16;
    INT16	ySrc B16;
} xRenderTriStripReq;

#define sz_xRenderTriStripReq			    24

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    CARD8	op;
    CARD8	pad1;
    CARD16	pad2 B16;
    Picture	src B32;
    Picture	dst B32;
    PictFormat	maskFormat B32;
    INT16	xSrc B16;
    INT16	ySrc B16;
} xRenderTriFanReq;

#define sz_xRenderTriFanReq			    24

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    Glyphset	gsid B32;
    PictFormat	format B32;
} xRenderCreateGlyphSetReq;

#define sz_xRenderCreateGlyphSetReq		    12

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    Glyphset    gsid B32;
    Glyphset    existing B32;
} xRenderReferenceGlyphSetReq;

#define sz_xRenderReferenceGlyphSetReq		    24

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    Glyphset    glyphset B32;
} xRenderFreeGlyphSetReq;

#define sz_xRenderFreeGlyphSetReq		    8

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    Glyphset    glyphset B32;
    CARD32	nglyphs;
} xRenderAddGlyphsReq;

#define sz_xRenderAddGlyphsReq			    12

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    Glyphset    glyphset B32;
} xRenderFreeGlyphsReq;

#define sz_xRenderFreeGlyphsReq			    8

typedef struct {
    CARD8       reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    CARD8	op;
    CARD8	pad1;
    CARD16	pad2 B16;
    Picture	src B32;
    Picture	dst B32;
    PictFormat	maskFormat B32;
    Glyphset    glyphset B32;
    INT16	xSrc B16;
    INT16	ySrc B16;
} xRenderCompositeGlyphsReq, xRenderCompositeGlyphs8Req,
xRenderCompositeGlyphs16Req, xRenderCompositeGlyphs32Req;

#define sz_xRenderCompositeGlyphs8Req		    28
#define sz_xRenderCompositeGlyphs16Req		    28
#define sz_xRenderCompositeGlyphs32Req		    28

/* 0.1 and higher */

typedef struct {
    CARD8	reqType;
    CARD8       renderReqType;
    CARD16      length B16;
    CARD8	op;
    CARD8	pad1;
    CARD16	pad2 B16;
    Picture	dst B32;
    xRenderColor    color;
} xRenderFillRectanglesReq;

#define sz_xRenderFillRectanglesReq		    20

/* 0.5 and higher */

typedef struct {
    CARD8	reqType;
    CARD8	renderReqType;
    CARD16	length B16;
    Cursor	cid B32;
    Picture	src B32;
    CARD16	x B16;
    CARD16	y B16;
} xRenderCreateCursorReq;

#define sz_xRenderCreateCursorReq		    16

/* 0.6 and higher */

/*
 * This can't use an array because 32-bit values may be in bitfields
 */
typedef struct {
    Fixed	matrix11 B32;
    Fixed	matrix12 B32;
    Fixed	matrix13 B32;
    Fixed	matrix21 B32;
    Fixed	matrix22 B32;
    Fixed	matrix23 B32;
    Fixed	matrix31 B32;
    Fixed	matrix32 B32;
    Fixed	matrix33 B32;
} xRenderTransform;

#define sz_xRenderTransform 36

typedef struct {
    CARD8		reqType;
    CARD8		renderReqType;
    CARD16		length B16;
    Picture		picture B32;
    xRenderTransform	transform;
} xRenderSetPictureTransformReq;

#define sz_xRenderSetPictureTransformReq	    44

typedef struct {
    CARD8		reqType;
    CARD8		renderReqType;
    CARD16		length B16;
    Drawable		drawable B32;
} xRenderQueryFiltersReq;

#define sz_xRenderQueryFiltersReq		    8

typedef struct {
    BYTE    type;   /* X_Reply */
    BYTE    pad1;
    CARD16  sequenceNumber B16;
    CARD32  length B32;
    CARD32  numAliases B32;	/* LISTofCARD16 */
    CARD32  numFilters B32;	/* LISTofSTRING8 */
    CARD32  pad2 B32;
    CARD32  pad3 B32;
    CARD32  pad4 B32;
    CARD32  pad5 B32;
} xRenderQueryFiltersReply;

#define sz_xRenderQueryFiltersReply		    32

typedef struct {
    CARD8		reqType;
    CARD8		renderReqType;
    CARD16		length B16;
    Picture		picture B32;
    CARD16		nbytes B16; /* number of bytes in name */
    CARD16		pad B16;
} xRenderSetPictureFilterReq;

#define sz_xRenderSetPictureFilterReq		    12

/* 0.8 and higher */

typedef struct {
    Cursor		cursor B32;
    CARD32		delay B32;
} xAnimCursorElt;

#define sz_xAnimCursorElt			    8

typedef struct {
    CARD8		reqType;
    CARD8		renderReqType;
    CARD16		length B16;
    Cursor		cid B32;
} xRenderCreateAnimCursorReq;

#define sz_xRenderCreateAnimCursorReq		    8

/* 0.9 and higher */

typedef struct {
    CARD8		reqType;
    CARD8		renderReqType;
    CARD16		length B16;
    Picture		picture;
    INT16		xOff B16;
    INT16		yOff B16;
} xRenderAddTrapsReq;

#define sz_xRenderAddTrapsReq			    12

/* 0.10 and higher */

typedef struct {
    CARD8	reqType;
    CARD8	renderReqType;
    CARD16	length B16;
    Picture	pid B32;
    xRenderColor color;
} xRenderCreateSolidFillReq;

#define sz_xRenderCreateSolidFillReq                 16

typedef struct {
    CARD8	reqType;
    CARD8	renderReqType;
    CARD16	length B16;
    Picture	pid B32;
    xPointFixed p1;
    xPointFixed p2;
    CARD32      nStops;
} xRenderCreateLinearGradientReq;

#define sz_xRenderCreateLinearGradientReq                 28

typedef struct {
    CARD8	reqType;
    CARD8	renderReqType;
    CARD16	length B16;
    Picture	pid B32;
    xPointFixed inner;
    xPointFixed outer;
    Fixed       inner_radius;
    Fixed       outer_radius;
    CARD32      nStops;
} xRenderCreateRadialGradientReq;

#define sz_xRenderCreateRadialGradientReq                 36

typedef struct {
    CARD8	reqType;
    CARD8	renderReqType;
    CARD16	length B16;
    Picture	pid B32;
    xPointFixed center;
    Fixed       angle; /* in degrees */
    CARD32      nStops;
} xRenderCreateConicalGradientReq;

#define sz_xRenderCreateConicalGradientReq                 24

#undef Window
#undef Drawable
#undef Font
#undef Pixmap
#undef Cursor
#undef Colormap
#undef GContext
#undef Atom
#undef VisualID
#undef Time
#undef KeyCode
#undef KeySym

#undef Picture
#undef PictFormat
#undef Fixed
#undef Glyphset

#endif /* _XRENDERP_H_ */
