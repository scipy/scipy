/************************************************************

Copyright 1989, 1998  The Open Group

Permission to use, copy, modify, distribute, and sell this software and its
documentation for any purpose is hereby granted without fee, provided that
the above copyright notice appear in all copies and that both that
copyright notice and this permission notice appear in supporting
documentation.

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
OPEN GROUP BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Except as contained in this notice, the name of The Open Group shall not be
used in advertising or otherwise to promote the sale, use or other dealings
in this Software without prior written authorization from The Open Group.

********************************************************/

/* THIS IS NOT AN X CONSORTIUM STANDARD OR AN X PROJECT TEAM SPECIFICATION */

#ifndef _SHMPROTO_H_
#define _SHMPROTO_H_

#include <X11/extensions/shm.h>

#define ShmSeg CARD32
#define Drawable CARD32
#define VisualID CARD32
#define GContext CARD32
#define Pixmap CARD32

#define X_ShmQueryVersion		0
#define X_ShmAttach			1
#define X_ShmDetach			2
#define X_ShmPutImage			3
#define X_ShmGetImage			4
#define X_ShmCreatePixmap		5
#define X_ShmAttachFd                   6
#define X_ShmCreateSegment              7

typedef struct _ShmQueryVersion {
    CARD8	reqType;		/* always ShmReqCode */
    CARD8	shmReqType;		/* always X_ShmQueryVersion */
    CARD16	length B16;
} xShmQueryVersionReq;
#define sz_xShmQueryVersionReq	4

typedef struct {
    BYTE	type;			/* X_Reply */
    BOOL	sharedPixmaps;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    CARD16	majorVersion B16;	/* major version of SHM protocol */
    CARD16	minorVersion B16;	/* minor version of SHM protocol */
    CARD16	uid B16;
    CARD16	gid B16;
    CARD8	pixmapFormat;
    CARD8	pad0;
    CARD16	pad1 B16;
    CARD32	pad2 B32;
    CARD32	pad3 B32;
    CARD32	pad4 B32;
} xShmQueryVersionReply;
#define sz_xShmQueryVersionReply	32

typedef struct _ShmAttach {
    CARD8	reqType;	/* always ShmReqCode */
    CARD8	shmReqType;	/* always X_ShmAttach */
    CARD16	length B16;
    ShmSeg	shmseg B32;
    CARD32	shmid B32;
    BOOL	readOnly;
    BYTE	pad0;
    CARD16	pad1 B16;
} xShmAttachReq;
#define sz_xShmAttachReq	16

typedef struct _ShmDetach {
    CARD8	reqType;	/* always ShmReqCode */
    CARD8	shmReqType;	/* always X_ShmDetach */
    CARD16	length B16;
    ShmSeg	shmseg B32;
} xShmDetachReq;
#define sz_xShmDetachReq	8

typedef struct _ShmPutImage {
    CARD8	reqType;	/* always ShmReqCode */
    CARD8	shmReqType;	/* always X_ShmPutImage */
    CARD16	length B16;
    Drawable	drawable B32;
    GContext	gc B32;
    CARD16	totalWidth B16;
    CARD16	totalHeight B16;
    CARD16	srcX B16;
    CARD16	srcY B16;
    CARD16	srcWidth B16;
    CARD16	srcHeight B16;
    INT16	dstX B16;
    INT16	dstY B16;
    CARD8	depth;
    CARD8	format;
    CARD8	sendEvent;
    CARD8	bpad;
    ShmSeg	shmseg B32;
    CARD32	offset B32;
} xShmPutImageReq;
#define sz_xShmPutImageReq	40

typedef struct _ShmGetImage {
    CARD8	reqType;	/* always ShmReqCode */
    CARD8	shmReqType;	/* always X_ShmGetImage */
    CARD16	length B16;
    Drawable	drawable B32;
    INT16	x B16;
    INT16	y B16;
    CARD16	width B16;
    CARD16	height B16;
    CARD32	planeMask B32;
    CARD8	format;
    CARD8	pad0;
    CARD8	pad1;
    CARD8	pad2;
    ShmSeg	shmseg B32;
    CARD32	offset B32;
} xShmGetImageReq;
#define sz_xShmGetImageReq	32

typedef struct _ShmGetImageReply {
    BYTE	type;  /* X_Reply */
    CARD8	depth;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    VisualID	visual B32;
    CARD32	size B32;
    CARD32	pad0 B32;
    CARD32	pad1 B32;
    CARD32	pad2 B32;
    CARD32	pad3 B32;
} xShmGetImageReply;
#define sz_xShmGetImageReply	32

typedef struct _ShmCreatePixmap {
    CARD8	reqType;	/* always ShmReqCode */
    CARD8	shmReqType;	/* always X_ShmCreatePixmap */
    CARD16	length B16;
    Pixmap	pid B32;
    Drawable	drawable B32;
    CARD16	width B16;
    CARD16	height B16;
    CARD8	depth;
    CARD8	pad0;
    CARD8	pad1;
    CARD8	pad2;
    ShmSeg	shmseg B32;
    CARD32	offset B32;
} xShmCreatePixmapReq;
#define sz_xShmCreatePixmapReq 28

typedef struct _ShmCompletion {
    BYTE	type;		/* always eventBase + ShmCompletion */
    BYTE	bpad0;
    CARD16	sequenceNumber B16;
    Drawable	drawable B32;
    CARD16	minorEvent B16;
    BYTE	majorEvent;
    BYTE	bpad1;
    ShmSeg	shmseg B32;
    CARD32	offset B32;
    CARD32	pad0 B32;
    CARD32	pad1 B32;
    CARD32	pad2 B32;
} xShmCompletionEvent;
#define sz_xShmCompletionEvent	32

/* Version 1.2 additions */
typedef struct _ShmAttachFd {
    CARD8	reqType;	/* always ShmReqCode */
    CARD8	shmReqType;	/* always X_ShmAttachFd */
    CARD16	length B16;
    ShmSeg	shmseg B32;
    BOOL	readOnly;
    BYTE	pad0;
    CARD16	pad1 B16;
} xShmAttachFdReq;
/* File descriptor is passed with this request */
#define sz_xShmAttachFdReq	12

typedef struct _ShmCreateSegment {
    CARD8	reqType;	/* always ShmReqCode */
    CARD8	shmReqType;	/* always X_ShmAttachFd */
    CARD16	length B16;
    ShmSeg	shmseg B32;
    CARD32      size B32;
    BOOL	readOnly;
    BYTE	pad0;
    CARD16	pad1 B16;
} xShmCreateSegmentReq;
#define sz_xShmCreateSegmentReq 16

typedef struct {
    CARD8	type;			/* must be X_Reply */
    CARD8	nfd;			/* must be 1	*/
    CARD16	sequenceNumber  B16;	/* last sequence number */
    CARD32	length  B32;		/* 0 */
    CARD32	pad2	B32;		/* unused */
    CARD32	pad3	B32;		/* unused */
    CARD32	pad4	B32;		/* unused */
    CARD32	pad5	B32;		/* unused */
    CARD32	pad6	B32;		/* unused */
    CARD32	pad7	B32;		/* unused */
} xShmCreateSegmentReply;
/* File descriptor is passed with this reply */
#define sz_xShmCreateSegmentReply	32

#undef ShmSeg
#undef Drawable
#undef VisualID
#undef GContext
#undef Pixmap

#endif /* _SHMPROTO_H_ */
