/*
 * DISPAT.H
 *
 * $Id$
 *
 * Declare generic dispatcher routines for event driven programs
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef DISPAT_H
#define DISPAT_H

/* A program may need to handle input events from several different
   devices-- stdin, one or more connections to X servers, UNIX pipelines,
   and the like.  A true event loop should:
     1) dispatch any events which have been read but not processed
        (pending events) BEFORE entering the loop, then loop:
     2) flush any pending output to the devices and
        pause until input arrives on any device,
     3) dispatch input to the appropriate event handler for each device
        on which input has arrived.

   Some devices, such as a connection to an X server, may require
   multiple additional handlers-- in the case of X, to distribute events
   to the handler for the toolkit which owns the window.

   A program creates a Dispatcher for each input connection for which it
   wants to receive events.  The program then calls DispatchEvents to
   enter a perpetual loop in which the Dispatch method for each Dispatcher
   is called as input becomes available on its file descriptor.
 */

/* Create a new dispatcher and add it to the list of active Dispatchers.
   The Pending or Flush methods may be 0 to skip those operations.
   Returns non-zero if the memory manager fails.  Any existing Dispatcher
   for the same file descriptor is replaced.
     Pending should handle any events which have already been read.  Its
      return values are the same as for Dispatch.
     Flush should flush output, if appropriate.
     Dispatch should read and process all available input, but must not
      block by attempting to read more than is available.  It should
      process all input that is read, then return 0.  If Dispatch
      returns >0, DispatchEvents returns immediately.  If Dispatch
      returns <0, DispatchEvents will return unless another Dispatcher
      exists which did not return -1 the last time it was invoked
      (this branch exists to prevent an infinite loop when the only
      input source is a terminal and the process is in the background).  */
extern int AddDispatcher(int fd, void *context, int (*Pending)(void *),
			 void (*Flush)(void *), int (*Dispatch)(void *));

/* Remove the dispatcher associated with the given file descriptor.
   RemoveDispatcher returns the context passed to AddDispatcher.  */
extern void *RemoveDispatcher(int fd);

/* Loop until some Dispatch returns non-0.  Note that Dispatchers may
   be created and removed within the DispatchEvents loop; if there are
   no active Dispatchers, DispatchEvents returns to its caller, as if
   a Dispatch returned non-0.  */
extern void DispatchEvents(void);

/* Loop until some Dispatch returns non-0, or until no more input events
   are available.  Same as DispatchEvents, but won't block.  */
extern void MaybeDispatch(void);

/* The DispatchWorker (if non-0) will be called if no input is available
   for any of the active Dispatchers.  It will be called repeatedly,
   checking for input between calls, until it returns 0.  */
extern int (*DispatchWorker)(void);

/* Test whether this file descriptor already has a Dispatcher.  */
extern int HasDispatcher(int fd);

/* The DispatchNext function is like DispatchEvents, but it only
   dispatches one event, returning the file descriptor for that event.
   Returns -1 if DispatchEvent would have exited without dispatching
   any events, -2 if a Dispatch returned >0, and -3 if a Pending
   returned >0.  */
extern int DispatchNext(void);

/* Replace DispatchError with routine of your choice to get non-default
   error behavior.  By default (DispatchError=0), a message is printed
   to stderr and processing continues.  If DispatchError returns non-zero,
   DispatchEvents will return to its caller, otherwise its loop
   continues.  */
extern int (*DispatchError)(void);

#endif
