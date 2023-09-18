/************************************************************
Copyright (c) 1997 by Silicon Graphics Computer Systems, Inc.
Permission to use, copy, modify, and distribute this
software and its documentation for any purpose and without
fee is hereby granted, provided that the above copyright
notice appear in all copies and that both that copyright
notice and this permission notice appear in supporting
documentation, and that the name of Silicon Graphics not be
used in advertising or publicity pertaining to distribution
of the software without specific prior written permission.
Silicon Graphics makes no representation about the suitability
of this software for any purpose. It is provided "as is"
without any express or implied warranty.
SILICON GRAPHICS DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE. IN NO EVENT SHALL SILICON
GRAPHICS BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE,
DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION  WITH
THE USE OR PERFORMANCE OF THIS SOFTWARE.
********************************************************/

#ifndef _EVIPROTO_H_
#define _EVIPROTO_H_

#include <X11/extensions/EVI.h>

#define X_EVIQueryVersion		0
#define X_EVIGetVisualInfo		1

#define VisualID CARD32

typedef CARD32 VisualID32;
#define sz_VisualID32 4

typedef struct _xExtendedVisualInfo {
    VisualID	core_visual_id B32;
    INT8	screen;
    INT8	level;
    CARD8	transparency_type;
    CARD8	pad0;
    CARD32	transparency_value B32;
    CARD8	min_hw_colormaps;
    CARD8	max_hw_colormaps;
    CARD16	num_colormap_conflicts B16;
} xExtendedVisualInfo;
#define sz_xExtendedVisualInfo 16

typedef struct _XEVIQueryVersion {
    CARD8	reqType;		/* always XEVIReqCode */
    CARD8	xeviReqType;		/* always X_EVIQueryVersion */
    CARD16	length B16;
} xEVIQueryVersionReq;
#define sz_xEVIQueryVersionReq	4

typedef struct {
    BYTE	type;			/* X_Reply */
    CARD8 	unused;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    CARD16	majorVersion B16;	/* major version of EVI protocol */
    CARD16	minorVersion B16;	/* minor version of EVI protocol */
    CARD32	pad0 B32;
    CARD32	pad1 B32;
    CARD32	pad2 B32;
    CARD32	pad3 B32;
    CARD32	pad4 B32;
} xEVIQueryVersionReply;
#define sz_xEVIQueryVersionReply	32

typedef struct _XEVIGetVisualInfoReq {
    CARD8	reqType;	/* always XEVIReqCode */
    CARD8	xeviReqType;	/* always X_EVIGetVisualInfo */
    CARD16      length B16;
    CARD32 	n_visual B32;
} xEVIGetVisualInfoReq;
#define sz_xEVIGetVisualInfoReq	8

typedef struct _XEVIGetVisualInfoReply {
    BYTE	type;  /* X_Reply */
    CARD8	unused;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    CARD32	n_info B32;
    CARD32	n_conflicts B32;
    CARD32	pad0 B32;
    CARD32	pad1 B32;
    CARD32	pad2 B32;
    CARD32	pad3 B32;
} xEVIGetVisualInfoReply;
#define sz_xEVIGetVisualInfoReply	32

#undef VisualID

#endif /* _EVIPROTO_H_ */
