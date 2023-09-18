/*

Copyright 1991, 1993, 1994, 1998  The Open Group

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

*/

/***********************************************************
Copyright 1991,1993 by Digital Equipment Corporation, Maynard, Massachusetts,
and Olivetti Research Limited, Cambridge, England.

                        All Rights Reserved

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in
supporting documentation, and that the names of Digital or Olivetti
not be used in advertising or publicity pertaining to distribution of the
software without specific, written prior permission.

DIGITAL AND OLIVETTI DISCLAIM ALL WARRANTIES WITH REGARD TO THIS
SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS, IN NO EVENT SHALL THEY BE LIABLE FOR ANY SPECIAL, INDIRECT OR
CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

******************************************************************/

#ifndef _SYNCPROTO_H_
#define _SYNCPROTO_H_

#include <X11/extensions/syncconst.h>

#define X_SyncInitialize		0
#define X_SyncListSystemCounters	1
#define X_SyncCreateCounter		2
#define X_SyncSetCounter		3
#define X_SyncChangeCounter		4
#define X_SyncQueryCounter              5
#define X_SyncDestroyCounter		6
#define X_SyncAwait			7
#define X_SyncCreateAlarm               8
#define X_SyncChangeAlarm	        9
#define X_SyncQueryAlarm	       10
#define X_SyncDestroyAlarm	       11
#define X_SyncSetPriority   	       12
#define X_SyncGetPriority   	       13
#define X_SyncCreateFence	       14
#define X_SyncTriggerFence	       15
#define X_SyncResetFence	       16
#define X_SyncDestroyFence	       17
#define X_SyncQueryFence	       18
#define X_SyncAwaitFence	       19

/* cover up types from sync.h to make sure they're the right size for
 * protocol packaging.  These will be undef'ed after all the protocol
 * structures are defined.
 */
#define XSyncCounter CARD32
#define XSyncAlarm   CARD32
#define XSyncFence   CARD32
#define Drawable     CARD32

/*
 * Initialize
 */
typedef struct _xSyncInitialize {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    CARD8	majorVersion;
    CARD8	minorVersion;
    CARD16	pad B16;
} xSyncInitializeReq;
#define sz_xSyncInitializeReq		8

typedef struct {
    BYTE	type;
    CARD8	unused;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    CARD8	majorVersion;
    CARD8	minorVersion;
    CARD16	pad B16;
    CARD32	pad0 B32;
    CARD32	pad1 B32;
    CARD32	pad2 B32;
    CARD32	pad3 B32;
    CARD32	pad4 B32;
} xSyncInitializeReply;
#define sz_xSyncInitializeReply	32

/*
 * ListSystemCounters
 */
typedef struct _xSyncListSystemCounters
{
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
} xSyncListSystemCountersReq;
#define sz_xSyncListSystemCountersReq	4

typedef struct {
    BYTE	type;
    CARD8	unused;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    INT32	nCounters B32;
    CARD32	pad0 B32;
    CARD32	pad1 B32;
    CARD32	pad2 B32;
    CARD32	pad3 B32;
    CARD32	pad4 B32;
} xSyncListSystemCountersReply;
#define sz_xSyncListSystemCountersReply	32

typedef struct {
    XSyncCounter counter B32;
    INT32	resolution_hi B32;
    CARD32	resolution_lo B32;
    CARD16	name_length B16;
} xSyncSystemCounter;
#define sz_xSyncSystemCounter 14

/*
 * Create Counter
 */
typedef struct _xSyncCreateCounterReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncCounter cid B32;
    INT32       initial_value_hi B32;
    CARD32	initial_value_lo B32;
} xSyncCreateCounterReq;
#define sz_xSyncCreateCounterReq	16

/*
 * Change Counter
 */
typedef struct _xSyncChangeCounterReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncCounter cid B32;
    INT32       value_hi B32;
    CARD32	value_lo B32;
} xSyncChangeCounterReq;
#define sz_xSyncChangeCounterReq	16

/*
 * Set Counter
 */
typedef struct _xSyncSetCounterReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncCounter cid B32;
    INT32       value_hi B32;
    CARD32	value_lo B32;
} xSyncSetCounterReq;
#define sz_xSyncSetCounterReq	16

/*
 * Destroy Counter
 */
typedef struct _xSyncDestroyCounterReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncCounter counter B32;
} xSyncDestroyCounterReq;
#define sz_xSyncDestroyCounterReq	8

/*
 * Query Counter
 */
typedef struct _xSyncQueryCounterReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncCounter counter B32;
} xSyncQueryCounterReq;
#define sz_xSyncQueryCounterReq		8


typedef struct {
    BYTE	type;
    CARD8	unused;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    INT32	value_hi B32;
    CARD32	value_lo B32;
    CARD32	pad0 B32;
    CARD32	pad1 B32;
    CARD32	pad2 B32;
    CARD32	pad3 B32;
} xSyncQueryCounterReply;
#define sz_xSyncQueryCounterReply	32

/*
 * Await
 */
typedef struct _xSyncAwaitReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
} xSyncAwaitReq;
#define sz_xSyncAwaitReq		4

typedef struct _xSyncWaitCondition {
    XSyncCounter counter B32;
    CARD32	value_type B32;
    INT32       wait_value_hi B32;
    CARD32      wait_value_lo B32;
    CARD32	test_type B32;
    INT32	event_threshold_hi B32;
    CARD32	event_threshold_lo B32;
} xSyncWaitCondition;
#define sz_xSyncWaitCondition		28

/*
 * Create Alarm
 */
typedef struct _xSyncCreateAlarmReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncAlarm	id B32;
    CARD32      valueMask B32;
} xSyncCreateAlarmReq;
#define sz_xSyncCreateAlarmReq		12

/*
 * Destroy Alarm
 */
typedef struct _xSyncDestroyAlarmReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncAlarm	alarm B32;
} xSyncDestroyAlarmReq;
#define sz_xSyncDestroyAlarmReq		8

/*
 * Query Alarm
 */
typedef struct _xSyncQueryAlarmReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncAlarm	alarm B32;
} xSyncQueryAlarmReq;
#define sz_xSyncQueryAlarmReq		8

typedef struct {
    BYTE	type;
    CARD8	unused;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    XSyncCounter counter B32;
    CARD32	value_type B32;
    INT32	wait_value_hi B32;
    CARD32	wait_value_lo B32;
    CARD32	test_type      B32;
    INT32	delta_hi B32;
    CARD32	delta_lo B32;
    BOOL        events;
    BYTE        state;
    BYTE	pad0;
    BYTE	pad1;
} xSyncQueryAlarmReply;
#define sz_xSyncQueryAlarmReply		40

/*
 * Change Alarm
 */
typedef struct _xSyncChangeAlarmReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncAlarm	alarm B32;
    CARD32	valueMask B32;
} xSyncChangeAlarmReq;
#define sz_xSyncChangeAlarmReq		12

/*
 * SetPriority
 */
typedef struct _xSyncSetPriority{
    CARD8   	reqType;
    CARD8   	syncReqType;
    CARD16  	length B16;
    CARD32  	id B32;
    INT32  	priority B32;
} xSyncSetPriorityReq;
#define sz_xSyncSetPriorityReq	    	12

/*
 * Get Priority
 */
typedef struct _xSyncGetPriority{
    CARD8   	reqType;
    CARD8   	syncReqType;
    CARD16  	length B16;
    CARD32  	id B32; /*XXX XID? */
} xSyncGetPriorityReq;
#define sz_xSyncGetPriorityReq	    	 8

typedef struct {
    BYTE	type;
    CARD8	unused;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    INT32  	priority B32;
    CARD32  	pad0 B32;
    CARD32  	pad1 B32;
    CARD32  	pad2 B32;
    CARD32  	pad3 B32;
    CARD32  	pad4 B32;
} xSyncGetPriorityReply;
#define sz_xSyncGetPriorityReply	32

/*
 * Create Fence
 */
typedef struct _xSyncCreateFenceReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    Drawable	d B32;
    XSyncFence	fid B32;
    BOOL	initially_triggered;
    CARD8	pad0;
    CARD16	pad1;
} xSyncCreateFenceReq;
#define sz_xSyncCreateFenceReq		16

/*
 * Put a fence object in the "triggered" state.
 */
typedef struct _xSyncTriggerFenceReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncFence	fid B32;
} xSyncTriggerFenceReq;
#define sz_xSyncTriggerFenceReq		8

/*
 * Put a fence in the "untriggered" state.
 */
typedef struct _xSyncResetFenceReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncFence	fid B32;
} xSyncResetFenceReq;
#define sz_xSyncResetFenceReq		8

/*
 * Destroy a fence object
 */
typedef struct _xSyncDestroyFenceReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncFence	fid B32;
} xSyncDestroyFenceReq;
#define sz_xSyncDestroyFenceReq		8

/*
 * Query a fence object
 */
typedef struct _xSyncQueryFenceReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
    XSyncFence	fid B32;
} xSyncQueryFenceReq;
#define sz_xSyncQueryFenceReq		8

/*
 * Wait for any of a list of fence sync objects
 * to reach the "triggered" state.
 */
typedef struct _xSyncAwaitFenceReq {
    CARD8	reqType;
    CARD8	syncReqType;
    CARD16	length B16;
} xSyncAwaitFenceReq;
#define sz_xSyncAwaitFenceReq		4

typedef struct {
    BYTE	type;
    CARD8	unused;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    BOOL	triggered;
    BYTE	pad0;
    CARD16	pad1 B16;
    CARD32	pad2 B32;
    CARD32	pad3 B32;
    CARD32	pad4 B32;
    CARD32	pad5 B32;
    CARD32	pad6 B32;
} xSyncQueryFenceReply;
#define sz_xSyncQueryFenceReply		32

/*
 * Events
 */

typedef struct _xSyncCounterNotifyEvent {
    BYTE	type;
    BYTE	kind;
    CARD16	sequenceNumber B16;
    XSyncCounter counter B32;
    INT32	wait_value_hi B32;
    CARD32	wait_value_lo B32;
    INT32	counter_value_hi B32;
    CARD32	counter_value_lo B32;
    CARD32	time B32;
    CARD16	count B16;
    BOOL	destroyed;
    BYTE        pad0;
} xSyncCounterNotifyEvent;

typedef struct _xSyncAlarmNotifyEvent {
    BYTE	type;
    BYTE	kind;
    CARD16	sequenceNumber B16;
    XSyncAlarm	alarm B32;
    INT32	counter_value_hi B32;
    CARD32	counter_value_lo B32;
    INT32	alarm_value_hi B32;
    CARD32	alarm_value_lo B32;
    CARD32	time B32;
    CARD8       state;
    BYTE        pad0;
    BYTE        pad1;
    BYTE        pad2;
} xSyncAlarmNotifyEvent;

#undef XSyncCounter
#undef XSyncAlarm
#undef XSyncFence
#undef Drawable


#endif /* _SYNCPROTO_H_ */
