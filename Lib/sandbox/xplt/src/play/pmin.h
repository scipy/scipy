/*
 * pmin.h -- $Id$
 * minimally intrusive play event handling (no stdio)
 */

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

/* call p_pending_events() to handle all pending
 *   events and expired alarms without blocking
 * call p_wait_while to set a flag and wait for play graphics events
 *   until flag is reset to zero
 * call p_xhandler to set abort_hook and on_exception routines
 *   a less obtrusive alternative to p_handler */
extern void p_pending_events(void);
extern void p_wait_while(int *flag);
extern void p_xhandler(void (*abort_hook)(void),
                       void (*on_exception)(int signal, char *errmsg));
extern int p_wait_stdin(void);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif
