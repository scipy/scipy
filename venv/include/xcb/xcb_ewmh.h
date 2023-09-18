#ifndef __XCB_EWMH_H__
#define __XCB_EWMH_H__

/*
 * Copyright (C) 2009-2011 Arnaud Fontaine <arnau@debian.org>
 *
 * Permission  is  hereby  granted,  free  of charge,  to  any  person
 * obtaining  a copy  of  this software  and associated  documentation
 * files   (the  "Software"),   to  deal   in  the   Software  without
 * restriction, including without limitation  the rights to use, copy,
 * modify, merge, publish,  distribute, sublicense, and/or sell copies
 * of  the Software, and  to permit  persons to  whom the  Software is
 * furnished to do so, subject to the following conditions:
 *
 * The  above copyright  notice and  this permission  notice  shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE  IS PROVIDED  "AS IS", WITHOUT  WARRANTY OF  ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT  NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY,   FITNESS    FOR   A   PARTICULAR    PURPOSE   AND
 * NONINFRINGEMENT. IN  NO EVENT SHALL  THE AUTHORS BE LIABLE  FOR ANY
 * CLAIM,  DAMAGES  OR  OTHER  LIABILITY,  WHETHER  IN  AN  ACTION  OF
 * CONTRACT, TORT OR OTHERWISE, ARISING  FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Except as  contained in  this notice, the  names of the  authors or
 * their institutions shall not be used in advertising or otherwise to
 * promote the  sale, use or  other dealings in this  Software without
 * prior written authorization from the authors.
 */

/**
 * @defgroup xcb__ewmh_t XCB EWMH Functions
 *
 * These functions  allow easy handling  of the protocol  described in
 * the Extended Window Manager  Hints specification. The list of Atoms
 * is stored as an M4 file  (atomlist.m4) where each Atom is stored as
 * a variable defined in the header.
 *
 * Replies of requests generating a  list of pointers (such as list of
 * windows, atoms and UTF-8 strings)  are simply stored as a structure
 * holding  the XCB  reply which  should (usually)  never  be accessed
 * directly and has  to be wipe afterwards. This  structure provides a
 * convenient access to the list given in the reply itself.
 *
 * @{
 */

#include <xcb/xcb.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Hold EWMH information specific to a screen
 */
typedef struct {
  /** The X connection */
  xcb_connection_t *connection;
  /** The screens on this connection */
  xcb_screen_t **screens;
  int nb_screens;
  /** _NET_WM_CM_Sn atoms depending on the number of screens */
  xcb_atom_t *_NET_WM_CM_Sn;
  /** The EWMH atoms of this connection */    
  xcb_atom_t _NET_SUPPORTED;
  xcb_atom_t _NET_CLIENT_LIST;
  xcb_atom_t _NET_CLIENT_LIST_STACKING;
  xcb_atom_t _NET_NUMBER_OF_DESKTOPS;
  xcb_atom_t _NET_DESKTOP_GEOMETRY;
  xcb_atom_t _NET_DESKTOP_VIEWPORT;
  xcb_atom_t _NET_CURRENT_DESKTOP;
  xcb_atom_t _NET_DESKTOP_NAMES;
  xcb_atom_t _NET_ACTIVE_WINDOW;
  xcb_atom_t _NET_WORKAREA;
  xcb_atom_t _NET_SUPPORTING_WM_CHECK;
  xcb_atom_t _NET_VIRTUAL_ROOTS;
  xcb_atom_t _NET_DESKTOP_LAYOUT;
  xcb_atom_t _NET_SHOWING_DESKTOP;
  xcb_atom_t _NET_CLOSE_WINDOW;
  xcb_atom_t _NET_MOVERESIZE_WINDOW;
  xcb_atom_t _NET_WM_MOVERESIZE;
  xcb_atom_t _NET_RESTACK_WINDOW;
  xcb_atom_t _NET_REQUEST_FRAME_EXTENTS;
  xcb_atom_t _NET_WM_NAME;
  xcb_atom_t _NET_WM_VISIBLE_NAME;
  xcb_atom_t _NET_WM_ICON_NAME;
  xcb_atom_t _NET_WM_VISIBLE_ICON_NAME;
  xcb_atom_t _NET_WM_DESKTOP;
  xcb_atom_t _NET_WM_WINDOW_TYPE;
  xcb_atom_t _NET_WM_STATE;
  xcb_atom_t _NET_WM_ALLOWED_ACTIONS;
  xcb_atom_t _NET_WM_STRUT;
  xcb_atom_t _NET_WM_STRUT_PARTIAL;
  xcb_atom_t _NET_WM_ICON_GEOMETRY;
  xcb_atom_t _NET_WM_ICON;
  xcb_atom_t _NET_WM_PID;
  xcb_atom_t _NET_WM_HANDLED_ICONS;
  xcb_atom_t _NET_WM_USER_TIME;
  xcb_atom_t _NET_WM_USER_TIME_WINDOW;
  xcb_atom_t _NET_FRAME_EXTENTS;
  xcb_atom_t _NET_WM_PING;
  xcb_atom_t _NET_WM_SYNC_REQUEST;
  xcb_atom_t _NET_WM_SYNC_REQUEST_COUNTER;
  xcb_atom_t _NET_WM_FULLSCREEN_MONITORS;
  xcb_atom_t _NET_WM_FULL_PLACEMENT;
  xcb_atom_t UTF8_STRING;
  xcb_atom_t WM_PROTOCOLS;
  xcb_atom_t MANAGER;
  xcb_atom_t _NET_WM_WINDOW_TYPE_DESKTOP;
  xcb_atom_t _NET_WM_WINDOW_TYPE_DOCK;
  xcb_atom_t _NET_WM_WINDOW_TYPE_TOOLBAR;
  xcb_atom_t _NET_WM_WINDOW_TYPE_MENU;
  xcb_atom_t _NET_WM_WINDOW_TYPE_UTILITY;
  xcb_atom_t _NET_WM_WINDOW_TYPE_SPLASH;
  xcb_atom_t _NET_WM_WINDOW_TYPE_DIALOG;
  xcb_atom_t _NET_WM_WINDOW_TYPE_DROPDOWN_MENU;
  xcb_atom_t _NET_WM_WINDOW_TYPE_POPUP_MENU;
  xcb_atom_t _NET_WM_WINDOW_TYPE_TOOLTIP;
  xcb_atom_t _NET_WM_WINDOW_TYPE_NOTIFICATION;
  xcb_atom_t _NET_WM_WINDOW_TYPE_COMBO;
  xcb_atom_t _NET_WM_WINDOW_TYPE_DND;
  xcb_atom_t _NET_WM_WINDOW_TYPE_NORMAL;
  xcb_atom_t _NET_WM_STATE_MODAL;
  xcb_atom_t _NET_WM_STATE_STICKY;
  xcb_atom_t _NET_WM_STATE_MAXIMIZED_VERT;
  xcb_atom_t _NET_WM_STATE_MAXIMIZED_HORZ;
  xcb_atom_t _NET_WM_STATE_SHADED;
  xcb_atom_t _NET_WM_STATE_SKIP_TASKBAR;
  xcb_atom_t _NET_WM_STATE_SKIP_PAGER;
  xcb_atom_t _NET_WM_STATE_HIDDEN;
  xcb_atom_t _NET_WM_STATE_FULLSCREEN;
  xcb_atom_t _NET_WM_STATE_ABOVE;
  xcb_atom_t _NET_WM_STATE_BELOW;
  xcb_atom_t _NET_WM_STATE_DEMANDS_ATTENTION;
  xcb_atom_t _NET_WM_ACTION_MOVE;
  xcb_atom_t _NET_WM_ACTION_RESIZE;
  xcb_atom_t _NET_WM_ACTION_MINIMIZE;
  xcb_atom_t _NET_WM_ACTION_SHADE;
  xcb_atom_t _NET_WM_ACTION_STICK;
  xcb_atom_t _NET_WM_ACTION_MAXIMIZE_HORZ;
  xcb_atom_t _NET_WM_ACTION_MAXIMIZE_VERT;
  xcb_atom_t _NET_WM_ACTION_FULLSCREEN;
  xcb_atom_t _NET_WM_ACTION_CHANGE_DESKTOP;
  xcb_atom_t _NET_WM_ACTION_CLOSE;
  xcb_atom_t _NET_WM_ACTION_ABOVE;
  xcb_atom_t _NET_WM_ACTION_BELOW;
} xcb_ewmh_connection_t;

/**
 * @brief Hold a GetProperty reply containing a list of Atoms
 */
typedef struct {
  /** The number of Atoms */
  uint32_t atoms_len;
  /** The list of Atoms */
  xcb_atom_t *atoms;
  /** The actual GetProperty reply */
  xcb_get_property_reply_t *_reply;
} xcb_ewmh_get_atoms_reply_t;

/**
 * @brief Hold a GetProperty reply containing a list of Windows
 */
typedef struct {
  /** The number of Windows */
  uint32_t windows_len;
  /** The list of Windows */
  xcb_window_t *windows;
  /** The actual GetProperty reply */
  xcb_get_property_reply_t *_reply;
} xcb_ewmh_get_windows_reply_t;

/**
 * @brief Hold a GetProperty reply containg a list of UTF-8 strings
 */
typedef struct {
  /** The number of UTF-8 strings */
  uint32_t strings_len;
  /** The list of UTF-8 strings */
  char *strings;
  /** The actual GetProperty reply */
  xcb_get_property_reply_t *_reply;
} xcb_ewmh_get_utf8_strings_reply_t;

/**
 * @brief Coordinates Property values
 */
typedef struct {
  /** The x coordinate */
  uint32_t x;
  /** The y coordinate */
  uint32_t y;
} xcb_ewmh_coordinates_t;

/**
 * @brief Hold reply of _NET_DESKTOP_VIEWPORT GetProperty
 */
typedef struct {
  /** The number of desktop viewports */
  uint32_t desktop_viewport_len;
  /** The desktop viewports */
  xcb_ewmh_coordinates_t *desktop_viewport;
  /** The actual GetProperty reply */
  xcb_get_property_reply_t *_reply;
} xcb_ewmh_get_desktop_viewport_reply_t;

/**
 * @brief Geometry Property values
 */
typedef struct {
  /** The x coordinate */
  uint32_t x;
  /** The y coordinate */
  uint32_t y;
  /** The width */
  uint32_t width;
  /** The height */
  uint32_t height;
} xcb_ewmh_geometry_t;

/**
 * @brief Hold reply of a _NET_WORKAREA GetProperty
 */
typedef struct {
  /** The number of desktop workarea */
  uint32_t workarea_len;
  /** The list of desktop workarea */
  xcb_ewmh_geometry_t *workarea;
  /** The actual GetProperty reply */
  xcb_get_property_reply_t *_reply;
} xcb_ewmh_get_workarea_reply_t;

/**
 * @brief Source indication in requests
 */
typedef enum {
  /** No source at all (for clients supporting an older version of
      EWMH specification) */
  XCB_EWMH_CLIENT_SOURCE_TYPE_NONE = 0,
  /** Normal application */
  XCB_EWMH_CLIENT_SOURCE_TYPE_NORMAL = 1,
  /** Pagers and other clients that represent direct user actions */
  XCB_EWMH_CLIENT_SOURCE_TYPE_OTHER = 2
} xcb_ewmh_client_source_type_t;

/**
 * @brief _NET_DESKTOP_LAYOUT orientation
 */
typedef enum {
  /** Horizontal orientation (desktops laid out in rows) */
  XCB_EWMH_WM_ORIENTATION_HORZ = 0,
  /** Vertical orientation (desktops laid out in columns) */
  XCB_EWMH_WM_ORIENTATION_VERT = 1
} xcb_ewmh_desktop_layout_orientation_t;

/**
 * @brief _NET_DESKTOP_LAYOUT starting corner
 */
typedef enum {
  /** Starting corner on the top left */
  XCB_EWMH_WM_TOPLEFT = 0,
  /** Starting corner on the top right */
  XCB_EWMH_WM_TOPRIGHT = 1,
  /** Starting corner on the bottom right */
  XCB_EWMH_WM_BOTTOMRIGHT = 2,
  /** Starting corner on the bottom left */
  XCB_EWMH_WM_BOTTOMLEFT = 3
} xcb_ewmh_desktop_layout_starting_corner_t;

/**
 * @brief Hold reply of a _NET_DESKTOP_LAYOUT GetProperty
 * @see xcb_ewmh_desktop_layout_orientation_t
 * @see xcb_ewmh_desktop_layout_starting_corner_t
 */
typedef struct {
  /** The desktops orientation */
  uint32_t orientation;
  /** The number of columns */
  uint32_t columns;
  /** The number of rows */
  uint32_t rows;
  /** The desktops starting corner */
  uint32_t starting_corner;
} xcb_ewmh_get_desktop_layout_reply_t;

/**
 * @brief _NET_WM_MOVERESIZE value when moving via keyboard
 * @see xcb_ewmh_moveresize_direction_t
 */
typedef enum {
  /** The window x coordinate */
  XCB_EWMH_MOVERESIZE_WINDOW_X = (1 << 8),
  /** The window y coordinate */
  XCB_EWMH_MOVERESIZE_WINDOW_Y = (1 << 9),
  /** The window width */
  XCB_EWMH_MOVERESIZE_WINDOW_WIDTH = (1 << 10),
  /** The window height */
  XCB_EWMH_MOVERESIZE_WINDOW_HEIGHT = (1 << 11)
} xcb_ewmh_moveresize_window_opt_flags_t;

/**
 * @brief _NET_WM_MOVERESIZE window movement or resizing
 */
typedef enum {
  /** Resizing applied on the top left edge */
  XCB_EWMH_WM_MOVERESIZE_SIZE_TOPLEFT = 0,
  /** Resizing applied on the top edge */
  XCB_EWMH_WM_MOVERESIZE_SIZE_TOP = 1,
  /** Resizing applied on the top right edge */
  XCB_EWMH_WM_MOVERESIZE_SIZE_TOPRIGHT = 2,
  /** Resizing applied on the right edge */
  XCB_EWMH_WM_MOVERESIZE_SIZE_RIGHT = 3,
  /** Resizing applied on the bottom right edge */
  XCB_EWMH_WM_MOVERESIZE_SIZE_BOTTOMRIGHT = 4,
  /** Resizing applied on the bottom edge */
  XCB_EWMH_WM_MOVERESIZE_SIZE_BOTTOM = 5,
  /** Resizing applied on the bottom left edge */
  XCB_EWMH_WM_MOVERESIZE_SIZE_BOTTOMLEFT = 6,
  /** Resizing applied on the left edge */
  XCB_EWMH_WM_MOVERESIZE_SIZE_LEFT = 7,
  /* Movement only */
  XCB_EWMH_WM_MOVERESIZE_MOVE = 8,
  /* Size via keyboard */
  XCB_EWMH_WM_MOVERESIZE_SIZE_KEYBOARD = 9,
  /* Move via keyboard */
  XCB_EWMH_WM_MOVERESIZE_MOVE_KEYBOARD = 10,
  /* Cancel operation */
  XCB_EWMH_WM_MOVERESIZE_CANCEL = 11
} xcb_ewmh_moveresize_direction_t;

/**
 * @brief Action on the _NET_WM_STATE property
 */
typedef enum {
  /* Remove/unset property */
  XCB_EWMH_WM_STATE_REMOVE = 0,
  /* Add/set property */
  XCB_EWMH_WM_STATE_ADD = 1,
  /* Toggle property  */
  XCB_EWMH_WM_STATE_TOGGLE = 2
} xcb_ewmh_wm_state_action_t;

/**
 * @brief Hold reply of _NET_WM_STRUT_PARTIAL GetProperty
 */
typedef struct {
  /** Reserved space on the left border of the screen */
  uint32_t left;
  /** Reserved space on the right border of the screen */
  uint32_t right;
  /** Reserved space on the top border of the screen */
  uint32_t top;
  /** Reserved space on the bottom border of the screen */
  uint32_t bottom;
  /** Beginning y coordinate of the left strut */
  uint32_t left_start_y;
  /** Ending y coordinate of the left strut */
  uint32_t left_end_y;
  /** Beginning y coordinate of the right strut */
  uint32_t right_start_y;
  /** Ending y coordinate of the right strut */
  uint32_t right_end_y;
  /** Beginning x coordinate of the top strut */
  uint32_t top_start_x;
  /** Ending x coordinate of the top strut */
  uint32_t top_end_x;
  /** Beginning x coordinate of the bottom strut */
  uint32_t bottom_start_x;
  /** Ending x coordinate of the bottom strut */
  uint32_t bottom_end_x;
} xcb_ewmh_wm_strut_partial_t;

/**
 * @brief Hold a single icon from reply of _NET_WM_ICON GetProperty
 */
typedef struct {
  /** Icon width */
  uint32_t width;
  /** Icon height */
  uint32_t height;
  /** Rows, left to right and top to bottom of the CARDINAL ARGB */
  uint32_t *data;
  /** Number of icons remaining */
  unsigned int rem;
  /** Index of the current icon in the array of icons */
  unsigned int index;
} xcb_ewmh_wm_icon_iterator_t;

/**
 * @brief Hold reply of _NET_WM_ICON GetProperty
 */
typedef struct {
  /** Number of icons */
  unsigned int num_icons;
  /** The actual GetProperty reply */
  xcb_get_property_reply_t *_reply;
} xcb_ewmh_get_wm_icon_reply_t;

/**
 * @brief Hold reply of _NET_REQUEST_FRAME_EXTENTS GetProperty
 */
typedef struct {
  /** Width of the left border */
  uint32_t left;
  /** Width of the right border */
  uint32_t right;
  /** Width of the top border */
  uint32_t top;
  /** Width of the bottom border */
  uint32_t bottom;
} xcb_ewmh_get_extents_reply_t;

/**
 * @brief Hold reply of _NET_WM_FULLSCREEN_MONITORS GetProperty
 */
typedef struct {
  /** Monitor whose top edge defines the top edge of the fullscreen
      window */
  uint32_t top;
  /** Monitor whose bottom edge defines the bottom edge of the
      fullscreen window */
  uint32_t bottom;
  /** Monitor whose left edge defines the left edge of the fullscreen
      window */
  uint32_t left;
  /** Monitor whose right edge defines the right edge of the
      fullscreen window */
  uint32_t right;
} xcb_ewmh_get_wm_fullscreen_monitors_reply_t;

/**
 * @brief Send InternAtom requests for the EWMH atoms and its required atoms
 *
 * @param c The connection to the X server
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The cookies corresponding to EWMH atoms
 */
xcb_intern_atom_cookie_t *xcb_ewmh_init_atoms(xcb_connection_t *c,
                                              xcb_ewmh_connection_t *ewmh);

/**
 * @brief Process  the replies  to the screen  initialisation requests
 * previously sent
 *
 * @param emwh The information relative to EWMH
 * @param ewmh_cookies The cookies corresponding to EWMH atoms
 * @param e Error if any
 * @return Return 1 on success, 0 otherwise
 */
uint8_t xcb_ewmh_init_atoms_replies(xcb_ewmh_connection_t *ewmh,
                                    xcb_intern_atom_cookie_t *ewmh_cookies,
                                    xcb_generic_error_t **e);

static inline void
xcb_ewmh_connection_wipe(xcb_ewmh_connection_t *ewmh)
{
  free(ewmh->screens);
  free(ewmh->_NET_WM_CM_Sn);
}

/**
 * @brief Send a SendEvent request containing a ClientMessage event
 *
 * This  function is  called  by all  the xcb_ewmh_request*  functions
 * whose should be used instead of calling directly this function
 *
 * @param c The X connection
 * @param window The window where the action will be applied
 * @param dest The destination window (usually the root window)
 * @param atom The type of the message
 * @param data_len The length of data to be sent
 * @param data The data to be sent
 * @return The cookie associated with the SendEvent request
 */
xcb_void_cookie_t xcb_ewmh_send_client_message(xcb_connection_t *c,
                                               xcb_window_t window,
                                               xcb_window_t dest,
                                               xcb_atom_t atom,
                                               uint32_t data_len,
                                               const uint32_t *data);

uint8_t xcb_ewmh_get_window_from_reply(xcb_window_t *window,
                                       xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_window_reply(xcb_ewmh_connection_t *ewmh,
                                  xcb_get_property_cookie_t cookie,
                                  xcb_window_t *window,
                                  xcb_generic_error_t **e);

uint8_t xcb_ewmh_get_cardinal_from_reply(uint32_t *cardinal,
                                         xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_cardinal_reply(xcb_ewmh_connection_t *ewmh,
                                    xcb_get_property_cookie_t cookie,
                                    uint32_t *cardinal,
                                    xcb_generic_error_t **e);

/**
 * @brief Get  a list  of atoms from  a given GetProperty  reply whose
 * type is ATOM
 *
 * This  function  is  called  by  all  the  xcb_ewmh_get_*_from_reply
 * functions  whose should  be used  instead of  calling  it directly.
 * Indeed,  The GetProperty request  has been  previously sent  by the
 * corresponding xcb_ewmh_get_*.
 *
 * @param atoms The atoms list
 * @param r The reply to get the atoms list from
 * @return Return 1 on success, 0 otherwise
 */
uint8_t xcb_ewmh_get_atoms_from_reply(xcb_ewmh_get_atoms_reply_t *atoms,
                                      xcb_get_property_reply_t *r);

/**
 * @brief Get a list of atoms  from the reply of a GetProperty request
 * whose type is ATOM
 *
 * This function  is called by all  the xcb_ewmh_get_*_reply functions
 * whose  should   be  used  instead  of  calling   it  directly.  The
 * GetProperty request  has been previously sent  by the corresponding
 * xcb_ewmh_get_*.
 *
 * @param ewmh The per-screen EWMH specific information
 * @param cookie The GetProperty cookie
 * @param atoms The atoms list
 * @param e The error if any
 * @return Return 1 on success, 0 otherwise
 */
uint8_t xcb_ewmh_get_atoms_reply(xcb_ewmh_connection_t *ewmh,
                                 xcb_get_property_cookie_t cookie,
                                 xcb_ewmh_get_atoms_reply_t *atoms,
                                 xcb_generic_error_t **e);

/**
 * @brief Wipe the Atoms list reply
 *
 * This function must be called to free the memory allocated for atoms
 * when the reply is requested in *_reply functions.
 *
 * @param data The X reply to be freed
 */
void xcb_ewmh_get_atoms_reply_wipe(xcb_ewmh_get_atoms_reply_t *data);

/**
 * @brief Get  a list  of atoms from  a given GetProperty  reply whose
 * type is WINDOW
 *
 * This  function  is  called  by  all  the  xcb_ewmh_get_*_from_reply
 * functions  whose should  be used  instead of  calling  it directly.
 * Indeed,  The GetProperty request  has been  previously sent  by the
 * corresponding xcb_ewmh_get_*.
 *
 * @param atoms The atoms list
 * @param r The reply to get the atoms list from
 * @return Return 1 on success, 0 otherwise
 */
uint8_t xcb_ewmh_get_windows_from_reply(xcb_ewmh_get_windows_reply_t *atoms,
                                        xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_utf8_strings_from_reply(xcb_ewmh_connection_t *ewmh,
                                             xcb_ewmh_get_utf8_strings_reply_t *data,
                                             xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_utf8_strings_reply(xcb_ewmh_connection_t *ewmh,
                                        xcb_get_property_cookie_t cookie,
                                        xcb_ewmh_get_utf8_strings_reply_t *data,
                                        xcb_generic_error_t **e);

/**
 * @brief Get a list of atoms  from the reply of a GetProperty request
 * whose type is WINDOW
 *
 * This function  is called by all  the xcb_ewmh_get_*_reply functions
 * whose  should   be  used  instead  of  calling   it  directly.  The
 * GetProperty request  has been previously sent  by the corresponding
 * xcb_ewmh_get_*.
 *
 * @param ewmh The per-screen EWMH specific information
 * @param cookie The GetProperty cookie
 * @param atoms The atoms list
 * @param e The error if any
 * @return Return 1 on success, 0 otherwise
 */
uint8_t xcb_ewmh_get_windows_reply(xcb_ewmh_connection_t *ewmh,
                                   xcb_get_property_cookie_t cookie,
                                   xcb_ewmh_get_windows_reply_t *atoms,
                                   xcb_generic_error_t **e);

/**
 * @brief Wipe the windows list reply
 *
 * This function must  be called to the free  the memory allocated for
 * windows when the reply is requested in '_reply' functions.
 *
 * @param data The X reply to be freed
 */
void xcb_ewmh_get_windows_reply_wipe(xcb_ewmh_get_windows_reply_t *data);

/**
 * @brief Wipe the UTF-8 strings list reply
 *
 * This function must  be called to the free  the memory allocated for
 * UTF-8 strings when the reply is requested in '_reply' functions.
 *
 * @param data The X reply to be freed
 */
void xcb_ewmh_get_utf8_strings_reply_wipe(xcb_ewmh_get_utf8_strings_reply_t *data);

/**
 * @brief Send a ChangeProperty request for _NET_SUPPORTED
 *
 * _NET_SUPPORTED, ATOM[]/32
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param list_len The number of Atoms supported by the WM
 * @param list The Atoms supported by the WM
 * @return Cookie associated with the ChangeProperty _NET_SUPPORTED request
 */
xcb_void_cookie_t xcb_ewmh_set_supported(xcb_ewmh_connection_t *ewmh,
                                         int screen_nbr,
                                         uint32_t list_len,
                                         xcb_atom_t *list);

/**
 * @see xcb_ewmh_set_supported
 */
xcb_void_cookie_t xcb_ewmh_set_supported_checked(xcb_ewmh_connection_t *ewmh,
                                                 int screen_nbr,
                                                 uint32_t list_len,
                                                 xcb_atom_t *list);

/**
 * @brief Send  GetProperty request to get  _NET_SUPPORTED root window
 *        property
 *
 * _NET_SUPPORTED, ATOM[]/32
 *
 * This property MUST  be set by the Window  Manager to indicate which
 * hints it supports. For example: considering _NET_WM_STATE both this
 * atom   and   all   supported  states   e.g.    _NET_WM_STATE_MODAL,
 * _NET_WM_STATE_STICKY, would be  listed. This assumes that backwards
 * incompatible changes will  not be made to the  hints (without being
 * renamed).
 *
 * This form can be used only if  the request will cause a reply to be
 * generated. Any returned error will be placed in the event queue.
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_SUPPORTED cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_supported_unchecked(xcb_ewmh_connection_t *ewmh,
                                                           int screen_nbr);

/**
 * @see xcb_ewmh_get_supported_unchecked
 */
xcb_get_property_cookie_t xcb_ewmh_get_supported(xcb_ewmh_connection_t *ewmh,
                                                 int screen_nbr);

/**
 * @brief Get the list of supported atoms
 *
 * @param supported The list of atoms contained in the reply
 * @param r GetProperty _NET_SUPPORTED reply
 */
static inline uint8_t
xcb_ewmh_get_supported_from_reply(xcb_ewmh_get_atoms_reply_t *supported,
                                  xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_atoms_from_reply(supported, r);
}

/**
 * @brief Get reply from the GetProperty _NET_SUPPORTED cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_window_supported_unchecked() is used.  Otherwise, it stores
 * the error if any.
 *
 * @param ewmh The information relative to EWMH
 * @param cookie The _NET_SUPPORTED GetProperty request cookie
 * @param supported The reply to be filled
 * @param The xcb_generic_error_t supplied
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_supported_reply(xcb_ewmh_connection_t *ewmh,
                             xcb_get_property_cookie_t cookie,
                             xcb_ewmh_get_atoms_reply_t *supported,
                             xcb_generic_error_t **e)
{
  return xcb_ewmh_get_atoms_reply(ewmh, cookie, supported, e);
}

/**
 * @brief Send a ChangeProperty request for _NET_CLIENT_LIST
 *
 * _NET_CLIENT_LIST, WINDOW[]/32
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param list_len The number of Atoms supported by the WM
 * @param list The Atoms supported by the WM
 * @return Cookie associated with the ChangeProperty _NET_CLIENT_LIST request
 */
xcb_void_cookie_t xcb_ewmh_set_client_list(xcb_ewmh_connection_t *ewmh,
                                           int screen_nbr,
                                           uint32_t list_len,
                                           xcb_window_t *list);

/**
 * @see xcb_ewmh_set_client_list
 */
xcb_void_cookie_t xcb_ewmh_set_client_list_checked(xcb_ewmh_connection_t *ewmh,
                                                   int screen_nbr,
                                                   uint32_t list_len,
                                                   xcb_window_t *list);

/**
 * @brief Send GetProperty request to get _NET_CLIENT_LIST root window
 *        property
 *
 * This  array   contains  all  X   Windows  managed  by   the  Window
 * Manager. _NET_CLIENT_LIST has  initial mapping order, starting with
 * the oldest window.  This property SHOULD be set  and updated by the
 * Window Manager.
 *
 * @param ewmh The information relative to EWMH.
 * @return The _NET_CLIENT_LIST cookie of the GetProperty request.
 */
xcb_get_property_cookie_t xcb_ewmh_get_client_list_unchecked(xcb_ewmh_connection_t *ewmh,
                                                             int screen_nbr);

/**
 * @brief Send GetProperty request to get _NET_CLIENT_LIST root window
 *        property
 *
 * @see xcb_ewmh_get_client_list_unchecked
 * @param ewmh The information relative to EWMH.
 * @return The _NET_CLIENT_LIST cookie of the GetProperty request.
 */
xcb_get_property_cookie_t xcb_ewmh_get_client_list(xcb_ewmh_connection_t *ewmh,
                                                   int screen_nbr);

/**
 * @brief  Get   the  list  of  client  windows   from  a  GetProperty
 * _NET_CLIENT_LIST reply
 *
 * @param clients The list of clients contained in the reply
 * @param r GetProperty _NET_CLIENT_LIST reply
 */
static inline uint8_t
xcb_ewmh_get_client_list_from_reply(xcb_ewmh_get_windows_reply_t *clients,
                                    xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_windows_from_reply(clients, r);
}

/**
 * @brief Get reply from the GetProperty _NET_CLIENT_LIST cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_window_client_list_unchecked()  is   used.   Otherwise,  it
 * stores the error if any.
 *
 * @param ewmh The information relative to EWMH
 * @param cookie The _NET_CLIENT_LIST GetProperty request cookie
 * @param clients The list of clients to be filled
 * @param The xcb_generic_error_t supplied
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_client_list_reply(xcb_ewmh_connection_t *ewmh,
                               xcb_get_property_cookie_t cookie,
                               xcb_ewmh_get_windows_reply_t *clients,
                               xcb_generic_error_t **e)
{
  return xcb_ewmh_get_windows_reply(ewmh, cookie, clients, e);
}

/**
 * @brief Send a ChangeProperty request for _NET_CLIENT_LIST_STACKING
 *
 * _NET_CLIENT_LIST_STACKING, WINDOW[]/32
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param list_len The number of Atoms supported by the WM
 * @param list The Atoms supported by the WM
 * @return Cookie associated with the ChangeProperty _NET_CLIENT_LIST_STACKING request
 */
xcb_void_cookie_t xcb_ewmh_set_client_list_stacking(xcb_ewmh_connection_t *ewmh,
                                                    int screen_nbr,
                                                    uint32_t list_len,
                                                    xcb_window_t *list);

/**
 * @see xcb_ewmh_set_client_list_stacking
 */
xcb_void_cookie_t xcb_ewmh_set_client_list_stacking_checked(xcb_ewmh_connection_t *ewmh,
                                                            int screen_nbr,
                                                            uint32_t list_len,
                                                            xcb_window_t *list);

/**
 * @brief  Send GetProperty  request to  get _NET_CLIENT_LIST_STACKING
 *        root window property
 *
 * This  array   contains  all  X   Windows  managed  by   the  Window
 * Manager.  _NET_CLIENT_LIST_STACKING   has  initial  mapping  order,
 * starting with the  oldest window.  This property SHOULD  be set and
 * updated by the Window Manager.
 *
 * @param ewmh The information relative to EWMH
 * @return The _NET_CLIENT_LIST_STACKING cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_client_list_stacking_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                      int screen_nbr);

/**
 * @see xcb_ewmh_get_client_list_unchecked
 */
xcb_get_property_cookie_t xcb_ewmh_get_client_list_stacking(xcb_ewmh_connection_t *ewmh,
                                                            int screen_nbr);

/**
 * @brief  Get   the  list  of  client  windows   from  a  GetProperty
 * _NET_CLIENT_LIST_STACKING reply
 *
 * @param clients The list of clients contained in the reply
 * @param r GetProperty _NET_CLIENT_LIST_STACKING reply
 */
static inline uint8_t
xcb_ewmh_get_client_list_stacking_from_reply(xcb_ewmh_get_windows_reply_t *clients,
                                             xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_windows_from_reply(clients, r);
}

/**
 * @brief  Get reply  from  the GetProperty  _NET_CLIENT_LIST_STACKING
 * cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_window_client_list_stacking_unchecked()       is      used.
 * Otherwise, it stores the error if any.
 *
 * @param ewmh The information relative to EWMH
 * @param cookie The _NET_CLIENT_LIST_STACKING GetProperty request cookie
 * @param clients The list of clients to be filled
 * @param The xcb_generic_error_t supplied
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_client_list_stacking_reply(xcb_ewmh_connection_t *ewmh,
                                        xcb_get_property_cookie_t cookie,
                                        xcb_ewmh_get_windows_reply_t *clients,
                                        xcb_generic_error_t **e)
{
  return xcb_ewmh_get_windows_reply(ewmh, cookie, clients, e);
}

/**
 * @brief Send a ChangeProperty request for _NET_NUMBER_OF_DESKTOPS
 *
 * _NET_NUMBER_OF_DESKTOPS? CARDINAL/32
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param number_of_desktops The number of desktops
 * @return Cookie associated with the ChangeProperty _NET_NUMBER_OF_DESKTOPS request
 */
xcb_void_cookie_t xcb_ewmh_set_number_of_desktops(xcb_ewmh_connection_t *ewmh,
                                                  int screen_nbr,
                                                  uint32_t number_of_desktops);

/**
 * @see xcb_ewmh_set_number_of_desktops
 */
xcb_void_cookie_t xcb_ewmh_set_number_of_desktops_checked(xcb_ewmh_connection_t *ewmh,
                                                          int screen_nbr,
                                                          uint32_t number_of_desktops);

/**
 * @brief Send GetProperty request to get _NET_NUMBER_OF_DESKTOPS root
 *        window property
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_NUMBER_OF_DESKTOPS cookie of the GetProperty request.
 */
xcb_get_property_cookie_t xcb_ewmh_get_number_of_desktops_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                    int screen_nbr);

/**
 * @see xcb_ewmh_get_number_of_desktops_unchecked
 */
xcb_get_property_cookie_t xcb_ewmh_get_number_of_desktops(xcb_ewmh_connection_t *ewmh,
                                                          int screen_nbr);

/**
 * @brief  Get   the  list  of  client  windows   from  a  GetProperty
 * _NET_NUMBER_OF_DESKTOPS reply
 *
 * @param clients The list of clients contained in the reply
 * @param r GetProperty _NET_NUMBER_OF_DESKTOPS reply
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_number_of_desktops_from_reply(uint32_t *number_of_desktops,
                                           xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_cardinal_from_reply(number_of_desktops, r);
}

/**
 * @brief Get reply from the GetProperty _NET_NUMBER_OF_DESKTOPS cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_window_number_of_desktops_unchecked()  is used.  Otherwise,
 * it stores the error if any.
 *
 * @param ewmh The information relative to EWMH
 * @param cookie The _NET_NUMBER_OF_DESKTOPS GetProperty request cookie
 * @param supported The reply to be filled
 * @param The xcb_generic_error_t supplied
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_number_of_desktops_reply(xcb_ewmh_connection_t *ewmh,
                                      xcb_get_property_cookie_t cookie,
                                      uint32_t *number_of_desktops,
                                      xcb_generic_error_t **e)
{
  return xcb_ewmh_get_cardinal_reply(ewmh, cookie, number_of_desktops, e);
}

static inline xcb_void_cookie_t
xcb_ewmh_request_change_number_of_desktops(xcb_ewmh_connection_t *ewmh,
                                           int screen_nbr,
                                           uint32_t new_number_of_desktops)
{
  return xcb_ewmh_send_client_message(ewmh->connection, XCB_NONE,
                                      ewmh->screens[screen_nbr]->root,
                                      ewmh->_NET_NUMBER_OF_DESKTOPS,
                                      sizeof(new_number_of_desktops),
                                      &new_number_of_desktops);
}

/**
 * @brief Send a ChangeProperty request for _NET_DESKTOP_GEOMETRY
 *
 * _NET_DESKTOP_GEOMETRY width, height, CARDINAL[2]/32
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param new_width The new desktop width
 * @param new_height The new desktop height
 * @return Cookie associated with the ChangeProperty _NET_DESKTOP_GEOMETRY request
 */
xcb_void_cookie_t xcb_ewmh_set_desktop_geometry(xcb_ewmh_connection_t *ewmh,
                                                int screen_nbr,
                                                uint32_t new_width,
                                                uint32_t new_height);

/**
 * @see xcb_ewmh_set_desktop_geometry
 */
xcb_void_cookie_t xcb_ewmh_set_desktop_geometry_checked(xcb_ewmh_connection_t *ewmh,
                                                        int screen_nbr,
                                                        uint32_t new_width,
                                                        uint32_t new_height);

/**
 * @brief Send  GetProperty request to  get _NET_DESKTOP_GEOMETRY root
 *        window property
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_DESKTOP_GEOMETRY cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_desktop_geometry_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                  int screen_nbr);

/**
 * @see xcb_ewmh_get_desktop_geometry_unchecked
 */
xcb_get_property_cookie_t xcb_ewmh_get_desktop_geometry(xcb_ewmh_connection_t *ewmh,
                                                        int screen_nbr);

/**
 * @brief Send ClientMessage requesting to change the _NET_DESKTOP_GEOMETRY
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param new_width The new desktop width
 * @param new_height The new desktop height
 * @return The SendEvent cookie
 */
xcb_void_cookie_t xcb_ewmh_request_change_desktop_geometry(xcb_ewmh_connection_t *ewmh,
                                                           int screen_nbr,
                                                           uint32_t new_width,
                                                           uint32_t new_height);

/**
 * @brief    Get   the   desktop    geometry   from    a   GetProperty
 * _NET_DESKTOP_GEOMETRY reply
 *
 * @param width The current desktop width
 * @param height The current desktop height
 * @param r GetProperty _NET_DESKTOP_GEOMETRY reply
 * @return Return 1 on success, 0 otherwise
 */
uint8_t xcb_ewmh_get_desktop_geometry_from_reply(uint32_t *width,
                                                 uint32_t *height,
                                                 xcb_get_property_reply_t *r);

/**
 * @brief Get reply from the GetProperty _NET_DESKTOP_GEOMETRY cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_desktop_geometry_unchecked() is used.  Otherwise, it stores
 * the error if any.
 *
 * @param ewmh The information relative to EWMH
 * @param cookie The _NET_DESKTOP_GEOMETRY GetProperty request cookie
 * @param width The current desktop width
 * @param width The current desktop height
 * @param The xcb_generic_error_t supplied
 * @return Return 1 on success, 0 otherwise
 */
uint8_t xcb_ewmh_get_desktop_geometry_reply(xcb_ewmh_connection_t *ewmh,
                                            xcb_get_property_cookie_t cookie,
                                            uint32_t *width, uint32_t *height,
                                            xcb_generic_error_t **e);

/**
 * @brief Send a ChangeProperty request for _NET_DESKTOP_VIEWPORT
 *
 * _NET_DESKTOP_VIEWPORT x, y, CARDINAL[][2]/32
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param list_len The number of desktop viewports
 * @param list The desktop viewports
 * @return Cookie associated with the ChangeProperty _NET_DESKTOP_VIEWPORT request
 */
xcb_void_cookie_t xcb_ewmh_set_desktop_viewport(xcb_ewmh_connection_t *ewmh,
                                                int screen_nbr,
                                                uint32_t list_len,
                                                xcb_ewmh_coordinates_t *list);

/**
 * @see xcb_ewmh_set_desktop_viewport
 */
xcb_void_cookie_t xcb_ewmh_set_desktop_viewport_checked(xcb_ewmh_connection_t *ewmh,
                                                        int screen_nbr,
                                                        uint32_t list_len,
                                                        xcb_ewmh_coordinates_t *list);

/**
 * @brief Send  GetProperty request to  get _NET_DESKTOP_VIEWPORT root
 *        window property
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_DESKTOP_VIEWPORT cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_desktop_viewport_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                  int screen_nbr);

/**
 * @see xcb_ewmh_get_desktop_viewport_unchecked
 */
xcb_get_property_cookie_t xcb_ewmh_get_desktop_viewport(xcb_ewmh_connection_t *ewmh,
                                                        int screen_nbr);

/**
 * @brief Send ClientMessage requesting to change the _NET_DESKTOP_VIEWPORT
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param new_x The new x coordinate
 * @param new_y The new y coordinate
 * @return The SendEvent cookie
 */
xcb_void_cookie_t xcb_ewmh_request_change_desktop_viewport(xcb_ewmh_connection_t *ewmh,
                                                           int screen_nbr,
                                                           uint32_t x, uint32_t y);

/**
 * @brief    Get   the   desktop    geometry   from    a   GetProperty
 * _NET_DESKTOP_VIEWPORT reply
 *
 * @param vp The  current desktop viewports
 * @param r GetProperty _NET_DESKTOP_VIEWPORT reply
 * @return Return 1 on success, 0 otherwise
 */
uint8_t xcb_ewmh_get_desktop_viewport_from_reply(xcb_ewmh_get_desktop_viewport_reply_t *vp,
                                                 xcb_get_property_reply_t *r);

/**
 * @brief Get reply from the GetProperty _NET_DESKTOP_VIEWPORT cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_desktop_viewport_unchecked() is used.  Otherwise, it stores
 * the error if any.
 *
 * @param ewmh The information relative to EWMH
 * @param cookie The _NET_DESKTOP_VIEWPORT GetProperty request cookie
 * @param vp The current desktop viewports
 * @param The xcb_generic_error_t supplied
 * @return Return 1 on success, 0 otherwise
 */
uint8_t xcb_ewmh_get_desktop_viewport_reply(xcb_ewmh_connection_t *ewmh,
                                            xcb_get_property_cookie_t cookie,
                                            xcb_ewmh_get_desktop_viewport_reply_t *vp,
                                            xcb_generic_error_t **e);

/**
 * @brief Wipe the desktop viewports list reply
 *
 * This function must be called to free the memory allocated for atoms
 * when the reply  is requested in xcb_ewmh_get_desktop_viewport_reply
 * function.
 *
 * @param r The X reply to be freed
 */
void xcb_ewmh_get_desktop_viewport_reply_wipe(xcb_ewmh_get_desktop_viewport_reply_t *r);

/**
 * @brief Send a ChangeProperty request for _NET_CURRENT_DESKTOP
 *
 * _NET_CURRENT_DESKTOP desktop, CARDINAL/32
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param new_current_desktop The new current desktop
 * @return Cookie associated with the ChangeProperty _NET_CURRENT_DESKTOP request
 */
xcb_void_cookie_t xcb_ewmh_set_current_desktop(xcb_ewmh_connection_t *ewmh,
                                               int screen_nbr,
                                               uint32_t new_current_desktop);

/**
 * @see xcb_ewmh_set_current_desktop
 */
xcb_void_cookie_t xcb_ewmh_set_current_desktop_checked(xcb_ewmh_connection_t *ewmh,
                                                       int screen_nbr,
                                                       uint32_t new_current_desktop);

/**
 * @brief  Send GetProperty request  to get  _NET_CURRENT_DESKTOP root
 *        window property
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_CURRENT_DESKTOP cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_current_desktop_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                 int screen_nbr);

/**
 * @see xcb_ewmh_get_current_desktop_unchecked
 */
xcb_get_property_cookie_t xcb_ewmh_get_current_desktop(xcb_ewmh_connection_t *ewmh,
                                                       int screen_nbr);

/**
 * @brief Send ClientMessage requesting to change the _NET_CURRENT_DESKTOP
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param new_desktop The new current desktop
 * @param timestamp The request timestamp
 * @return The SendEvent cookie
 */
xcb_void_cookie_t xcb_ewmh_request_change_current_desktop(xcb_ewmh_connection_t *ewmh,
                                                          int screen_nbr,
                                                          uint32_t new_desktop,
                                                          xcb_timestamp_t timestamp);

/**
 * @brief    Get   the   desktop    geometry   from    a   GetProperty
 * _NET_CURRENT_DESKTOP reply
 *
 * @param current_desktop The  current desktop
 * @param r GetProperty _NET_CURRENT_DESKTOP reply
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_current_desktop_from_reply(uint32_t *current_desktop,
                                        xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_cardinal_from_reply(current_desktop, r);
}

/**
 * @brief Get reply from the GetProperty _NET_CURRENT_DESKTOP cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_current_desktop_unchecked() is  used.  Otherwise, it stores
 * the error if any.
 *
 * @param ewmh The information relative to EWMH
 * @param cookie The _NET_CURRENT_DESKTOP GetProperty request cookie
 * @param vp The current desktop
 * @param The xcb_generic_error_t supplied
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_current_desktop_reply(xcb_ewmh_connection_t *ewmh,
                                   xcb_get_property_cookie_t cookie,
                                   uint32_t *current_desktop,
                                   xcb_generic_error_t **e)
{
  return xcb_ewmh_get_cardinal_reply(ewmh, cookie, current_desktop, e);
}

/**
 * @brief Send a ChangeProperty request for _NET_DESKTOP_NAMES
 *
 * _NET_DESKTOP_NAMES, UTF8_STRING[]
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param strings_len The number of desktop names
 * @param strings The desktop names
 * @return Cookie associated with the ChangeProperty _NET_DESKTOP_NAMES request
 */
xcb_void_cookie_t xcb_ewmh_set_desktop_names(xcb_ewmh_connection_t *ewmh,
                                             int screen_nbr,
                                             uint32_t strings_len,
                                             const char *strings);

/**
 * @see xcb_ewmh_set_desktop_names
 */
xcb_void_cookie_t xcb_ewmh_set_desktop_names_checked(xcb_ewmh_connection_t *ewmh,
                                                     int screen_nbr,
                                                     uint32_t strings_len,
                                                     const char *strings);

/**
 * @brief  Send  GetProperty request  to  get _NET_DESKTOP_NAMES  root
 *        window property
 *
 * @param ewmh The information relative to EWMH
 * @return The _NET_DESKTOP_NAMES cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_desktop_names_unchecked(xcb_ewmh_connection_t *ewmh,
                                                               int screen_nbr);

/**
 * @see xcb_ewmh_get_desktop_names_unchecked
 */
xcb_get_property_cookie_t xcb_ewmh_get_desktop_names(xcb_ewmh_connection_t *ewmh,
                                                     int screen_nbr);

/**
 * @brief    Get   the   desktop    geometry   from    a   GetProperty
 * _NET_DESKTOP_NAMES reply
 *
 * @param ewmh The information relative to EWMH
 * @param names The desktop names
 * @param r GetProperty _NET_DESKTOP_NAMES reply
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_desktop_names_from_reply(xcb_ewmh_connection_t *ewmh,
                                      xcb_ewmh_get_utf8_strings_reply_t *names,
                                      xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_utf8_strings_from_reply(ewmh, names, r);
}

/**
 * @brief Get reply from the GetProperty _NET_DESKTOP_NAMES cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_desktop_names_unchecked()  is used.   Otherwise,  it stores
 * the error if any.
 *
 * @param ewmh The information relative to EWMH
 * @param cookie The _NET_DESKTOP_NAMES GetProperty request cookie
 * @param names The desktop names
 * @param The xcb_generic_error_t supplied
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_desktop_names_reply(xcb_ewmh_connection_t *ewmh,
                                 xcb_get_property_cookie_t cookie,
                                 xcb_ewmh_get_utf8_strings_reply_t *names,
                                 xcb_generic_error_t **e)
{
  return xcb_ewmh_get_utf8_strings_reply(ewmh, cookie, names, e);
}

/**
 * @brief Send a ChangeProperty request for _NET_ACTIVE_WINDOW
 *
 * _NET_ACTIVE_WINDOW, WINDOW/32
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param new_active_window The window to make active
 * @return Cookie associated with the ChangeProperty _NET_ACTIVE_WINDOW request
 */
xcb_void_cookie_t xcb_ewmh_set_active_window(xcb_ewmh_connection_t *ewmh,
                                             int screen_nbr,
                                             xcb_window_t new_active_window);

/**
 * @see xcb_ewmh_set_active_window
 */
xcb_void_cookie_t xcb_ewmh_set_active_window_checked(xcb_ewmh_connection_t *ewmh,
                                                     int screen_nbr,
                                                     xcb_window_t new_active_window);

/**
 * @brief Send ClientMessage requesting to change the _NET_ACTIVE_WINDOW
 *
 * The window ID  of the currently active window or  None if no window
 * has  the focus.  This  is a  read-only property  set by  the Window
 * Manager. If a Client wants to activate another window, it MUST send
 * a  _NET_ACTIVE_WINDOW  client  message  to  the  root  window.  The
 * timestamp is Client's  last user activity timestamp at  the time of
 * the request, and the currently active window is the Client's active
 * toplevel window, if any (the Window Manager may be e.g. more likely
 * to obey  the request  if it will  mean transferring focus  from one
 * active window to another).
 *
 * @see xcb_ewmh_client_source_type_t
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @param window_to_active The window ID to activate
 * @param source_indication The source indication
 * @param timestamp The client's last user activity timestamp
 * @param current_active_window The currently active window or None
 */
xcb_void_cookie_t xcb_ewmh_request_change_active_window(xcb_ewmh_connection_t *ewmh,
                                                        int screen_nbr,
                                                        xcb_window_t window_to_activate,
                                                        xcb_ewmh_client_source_type_t source_indication,
                                                        xcb_timestamp_t timestamp,
                                                        xcb_window_t current_active_window);

/**
 * @brief  Send  GetProperty request  to  get _NET_ACTIVE_WINDOW  root
 *        window property
 *
 * The window ID  of the currently active window or  None if no window
 * has  the focus.  This is  a read-only  property set  by  the Window
 * Manager.  This property  SHOULD be  set and  updated by  the Window
 * Manager.
 *
 * This form can be used only if  the request will cause a reply to be
 * generated. Any returned error will be placed in the event queue.
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_ACTIVE_WINDOW cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_active_window_unchecked(xcb_ewmh_connection_t *ewmh,
                                                               int screen_nbr);

/**
 * @brief  Send  GetProperty request  to  get _NET_ACTIVE_WINDOW  root
 *        window property
 *
 * @see xcb_ewmh_get_active_window_unchecked
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_ACTIVE_WINDOW cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_active_window(xcb_ewmh_connection_t *ewmh,
                                                     int screen_nbr);

/**
 * @brief  Get   the  list  of  client  windows   from  a  GetProperty
 * _NET_ACTIVE_WINDOW reply
 *
 * @param active_window The current active window
 * @param r GetProperty _NET_ACTIVE_WINDOW_OF_DESKTOPS reply
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_active_window_from_reply(xcb_window_t *active_window,
                                      xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_window_from_reply(active_window, r);
}

/**
 * @brief Get reply from the GetProperty _NET_ACTIVE_WINDOW cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_active_window_unchecked()  is used.   Otherwise,  it stores
 * the error if any.
 *
 * @param ewmh The information relative to EWMH.
 * @param cookie The _NET_ACTIVE_WINDOW GetProperty request cookie.
 * @param active_window The reply to be filled.
 * @param The xcb_generic_error_t supplied.
 * @return Return 1 on success, 0 otherwise.
 */
static inline uint8_t
xcb_ewmh_get_active_window_reply(xcb_ewmh_connection_t *ewmh,
                                 xcb_get_property_cookie_t cookie,
                                 xcb_window_t *active_window,
                                 xcb_generic_error_t **e)
{
  return xcb_ewmh_get_window_reply(ewmh, cookie, active_window, e);
}

/**
 * @brief Send a ChangeProperty request for _NET_WORKAREA
 *
 * _NET_WORKAREA, x, y, width, height CARDINAL[][4]/32
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param list_len The number of desktops workareas
 * @param list The desktops workareas
 * @return Cookie associated with the ChangeProperty _NET_WORKAREA request
 */
xcb_void_cookie_t xcb_ewmh_set_workarea(xcb_ewmh_connection_t *ewmh,
                                        int screen_nbr,
                                        uint32_t list_len,
                                        xcb_ewmh_geometry_t *list);

/**
 * @see xcb_ewmh_set_workarea
 */
xcb_void_cookie_t xcb_ewmh_set_workarea_checked(xcb_ewmh_connection_t *ewmh,
                                                int screen_nbr,
                                                uint32_t list_len,
                                                xcb_ewmh_geometry_t *list);

/**
 * @brief  Send  GetProperty request  to  get _NET_WORKAREA  root
 *        window property
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_WORKAREA cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_workarea_unchecked(xcb_ewmh_connection_t *ewmh,
                                                          int screen_nbr);

/**
 * @see xcb_ewmh_get_virtual_roots_unchecked
 */
xcb_get_property_cookie_t xcb_ewmh_get_workarea(xcb_ewmh_connection_t *ewmh,
                                                int screen_nbr);

/**
 * @brief Get  the desktop  geometry from a  GetProperty _NET_WORKAREA
 * reply
 *
 * @param wa The  current workarea
 * @param r GetProperty _NET_WORKAREA reply
 * @return Return 1 on success, 0 otherwise
 */
uint8_t xcb_ewmh_get_workarea_from_reply(xcb_ewmh_get_workarea_reply_t *wa,
                                         xcb_get_property_reply_t *r);

/**
 * @brief Get reply from the GetProperty _NET_WORKAREA cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_workarea_unchecked()  is used.   Otherwise,  it stores  the
 * error if any.
 *
 * @param ewmh The information relative to EWMH
 * @param cookie The _NET_WORKAREA GetProperty request cookie
 * @param wa The current workareas of desktops
 * @param The xcb_generic_error_t supplied
 * @return Return 1 on success, 0 otherwise
 */
uint8_t xcb_ewmh_get_workarea_reply(xcb_ewmh_connection_t *ewmh,
                                    xcb_get_property_cookie_t cookie,
                                    xcb_ewmh_get_workarea_reply_t *wa,
                                    xcb_generic_error_t **e);

/**
 * @brief Wipe the workarea list reply
 *
 * This function must be called to free the memory allocated for atoms
 * when   the  reply   is  requested   in  xcb_ewmh_get_workarea_reply
 * function.
 *
 * @param r The X reply to be freed
 */
void xcb_ewmh_get_workarea_reply_wipe(xcb_ewmh_get_workarea_reply_t *r);

/**
 * @brief Send a ChangeProperty request for _NET_SUPPORTING_WM_CHECK
 *
 * _NET_SUPPORTING_WM_CHECK, WINDOW/32
 *
 * @param ewmh The per-screen EWMH information
 * @param parent_window The root window or child window created by the WM
 * @param child_window The child window created by the WM
 * @return Cookie associated with the ChangeProperty _NET_SUPPORTING_WM_CHECK request
 */
xcb_void_cookie_t xcb_ewmh_set_supporting_wm_check(xcb_ewmh_connection_t *ewmh,
                                                   xcb_window_t parent_window,
                                                   xcb_window_t child_window);

/**
 * @see xcb_ewmh_set_supporting_wm_check
 */
xcb_void_cookie_t xcb_ewmh_set_supporting_wm_check_checked(xcb_ewmh_connection_t *ewmh,
                                                           xcb_window_t parent_window,
                                                           xcb_window_t child_window);

/**
 * @brief  Send GetProperty  request  to get  _NET_SUPPORTING_WM_CHECK
 *        root window property
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_SUPPORTING_WM_CHECK cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_supporting_wm_check_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                     xcb_window_t window);

/**
 * @see xcb_ewmh_get_supporting_wm_check_unchecked
 */
xcb_get_property_cookie_t xcb_ewmh_get_supporting_wm_check(xcb_ewmh_connection_t *ewmh,
                                                           xcb_window_t window);

/**
 * @brief  Get   the  list  of  client  windows   from  a  GetProperty
 * _NET_SUPPORTING_WM_CHECK reply
 *
 * @param window The child window created by the WM
 * @param r GetProperty _NET_SUPPORTING_WM_CHECK reply
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_supporting_wm_check_from_reply(xcb_window_t *window,
                                            xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_window_from_reply(window, r);
}

/**
 * @brief  Get  reply  from the  GetProperty  _NET_SUPPORTING_WM_CHECK
 * cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_supporting_wm_check_unchecked()  is  used.   Otherwise,  it
 * stores the error if any.
 *
 * @param ewmh The information relative to EWMH
 * @param cookie The _NET_SUPPORTING_WM_CHECK GetProperty request cookie
 * @param window The reply to be filled
 * @param The xcb_generic_error_t supplied
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_supporting_wm_check_reply(xcb_ewmh_connection_t *ewmh,
                                       xcb_get_property_cookie_t cookie,
                                       xcb_window_t *window,
                                       xcb_generic_error_t **e)
{
  return xcb_ewmh_get_window_reply(ewmh, cookie, window, e);
}

/**
 * @brief Send a ChangeProperty request for _NET_VIRTUAL_ROOTS
 *
 * _NET_VIRTUAL_ROOTS, WINDOW[]/32
 *
 * @param ewmh The per-screen EWMH information
 * @param screen_nbr The screen number
 * @param list_len The number of virtual root windows
 * @param list The virtual root windows
 * @return Cookie associated with the ChangeProperty _NET_VIRTUAL_ROOTS request
 */
xcb_void_cookie_t xcb_ewmh_set_virtual_roots(xcb_ewmh_connection_t *ewmh,
                                             int screen_nbr,
                                             uint32_t list_len,
                                             xcb_window_t *list);

/**
 * @see xcb_ewmh_set_virtual_roots
 */
xcb_void_cookie_t xcb_ewmh_set_virtual_roots_checked(xcb_ewmh_connection_t *ewmh,
                                                     int screen_nbr,
                                                     uint32_t list_len,
                                                     xcb_window_t *list);

/**
 * @brief  Send  GetProperty request  to  get _NET_VIRTUAL_ROOTS  root
 *        window property
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_VIRTUAL_ROOTS cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_virtual_roots_unchecked(xcb_ewmh_connection_t *ewmh,
                                                               int screen_nbr);

/**
 * @see xcb_ewmh_get_virtual_roots_unchecked
 */
xcb_get_property_cookie_t xcb_ewmh_get_virtual_roots(xcb_ewmh_connection_t *ewmh,
                                                     int screen_nbr);

/**
 * @brief Get  the desktop  geometry from a  GetProperty _NET_WORKAREA
 * reply
 *
 * @param virtual_roots The current virtual root windows
 * @param r GetProperty _NET_VIRTUAL_ROOTS reply
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_virtual_roots_from_reply(xcb_ewmh_get_windows_reply_t *virtual_roots,
                                      xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_windows_from_reply(virtual_roots, r);
}

/**
 * @brief Get reply from the GetProperty _NET_VIRTUAL_ROOTS cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_virtual_roots_unchecked()  is used.   Otherwise,  it stores
 * the error if any.
 *
 * @param ewmh The information relative to EWMH
 * @param cookie The _NET_VIRTUAL_ROOTS GetProperty request cookie
 * @param virtual_roots The current virtual root windows
 * @param The xcb_generic_error_t supplied
 * @return Return 1 on success, 0 otherwise
 */
static inline uint8_t
xcb_ewmh_get_virtual_roots_reply(xcb_ewmh_connection_t *ewmh,
                                 xcb_get_property_cookie_t cookie,
                                 xcb_ewmh_get_windows_reply_t *virtual_roots,
                                 xcb_generic_error_t **e)
{
  return xcb_ewmh_get_windows_reply(ewmh, cookie, virtual_roots, e);
}

xcb_void_cookie_t xcb_ewmh_set_desktop_layout(xcb_ewmh_connection_t *ewmh,
                                              int screen_nbr,
                                              xcb_ewmh_desktop_layout_orientation_t orientation,
                                              uint32_t columns, uint32_t rows,
                                              xcb_ewmh_desktop_layout_starting_corner_t starting_corner);

xcb_void_cookie_t xcb_ewmh_set_desktop_layout_checked(xcb_ewmh_connection_t *ewmh,
                                                      int screen_nbr,
                                                      xcb_ewmh_desktop_layout_orientation_t orientation,
                                                      uint32_t columns, uint32_t rows,
                                                      xcb_ewmh_desktop_layout_starting_corner_t starting_corner);

/**
 * @brief  Send GetProperty  request to  get  _NET_DESKTOP_LAYOUT root
 *        window property
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_DESKTOP_LAYOUT cookie of the GetProperty request
 */
xcb_get_property_cookie_t xcb_ewmh_get_desktop_layout_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                int screen_nbr);

/**
 * @see xcb_ewmh_get_desktop_layout_unchecked
 */
xcb_get_property_cookie_t xcb_ewmh_get_desktop_layout(xcb_ewmh_connection_t *ewmh,
                                                      int screen_nbr);

uint8_t xcb_ewmh_get_desktop_layout_from_reply(xcb_ewmh_get_desktop_layout_reply_t *desktop_layouts,
                                               xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_desktop_layout_reply(xcb_ewmh_connection_t *ewmh,
                                          xcb_get_property_cookie_t cookie,
                                          xcb_ewmh_get_desktop_layout_reply_t *desktop_layouts,
                                          xcb_generic_error_t **e);

xcb_void_cookie_t xcb_ewmh_set_showing_desktop(xcb_ewmh_connection_t *ewmh,
                                               int screen_nbr,
                                               uint32_t desktop);

xcb_void_cookie_t xcb_ewmh_set_showing_desktop_checked(xcb_ewmh_connection_t *ewmh,
                                                       int screen_nbr,
                                                       uint32_t desktop);

xcb_get_property_cookie_t xcb_ewmh_get_showing_desktop_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                 int screen_nbr);

xcb_get_property_cookie_t xcb_ewmh_get_showing_desktop(xcb_ewmh_connection_t *ewmh,
                                                       int screen_nbr);

static inline uint8_t
xcb_ewmh_get_showing_desktop_from_reply(uint32_t *desktop,
                                        xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_cardinal_from_reply(desktop, r);
}

static inline uint8_t
xcb_ewmh_get_showing_desktop_reply(xcb_ewmh_connection_t *ewmh,
                                   xcb_get_property_cookie_t cookie,
                                   uint32_t *desktop,
                                   xcb_generic_error_t **e)
{
  return xcb_ewmh_get_cardinal_reply(ewmh, cookie, desktop, e);
}

static inline xcb_void_cookie_t
xcb_ewmh_request_change_showing_desktop(xcb_ewmh_connection_t *ewmh,
                                        int screen_nbr,
                                        uint32_t enter)
{
  return xcb_ewmh_send_client_message(ewmh->connection, XCB_NONE,
                                      ewmh->screens[screen_nbr]->root,
                                      ewmh->_NET_SHOWING_DESKTOP,
                                      sizeof(enter), &enter);
}

xcb_void_cookie_t xcb_ewmh_request_close_window(xcb_ewmh_connection_t *ewmh,
                                                int screen_nbr,
                                                xcb_window_t window_to_close,
                                                xcb_timestamp_t timestamp,
                                                xcb_ewmh_client_source_type_t source_indication);

xcb_void_cookie_t xcb_ewmh_request_moveresize_window(xcb_ewmh_connection_t *ewmh,
                                                     int screen_nbr,
                                                     xcb_window_t moveresize_window,
                                                     xcb_gravity_t gravity,
                                                     xcb_ewmh_client_source_type_t source_indication,
                                                     xcb_ewmh_moveresize_window_opt_flags_t flags,
                                                     uint32_t x, uint32_t y,
                                                     uint32_t width, uint32_t height);

xcb_void_cookie_t xcb_ewmh_request_wm_moveresize(xcb_ewmh_connection_t *ewmh,
                                                 int screen_nbr,
                                                 xcb_window_t moveresize_window,
                                                 uint32_t x_root, uint32_t y_root,
                                                 xcb_ewmh_moveresize_direction_t direction,
                                                 xcb_button_index_t button,
                                                 xcb_ewmh_client_source_type_t source_indication);

xcb_void_cookie_t xcb_ewmh_request_restack_window(xcb_ewmh_connection_t *ewmh,
                                                  int screen_nbr,
                                                  xcb_window_t window_to_restack,
                                                  xcb_window_t sibling_window,
                                                  xcb_stack_mode_t detail);

static inline xcb_void_cookie_t
xcb_ewmh_request_frame_extents(xcb_ewmh_connection_t *ewmh,
                               int screen_nbr,
                               xcb_window_t client_window)
{
  return xcb_ewmh_send_client_message(ewmh->connection, client_window,
                                      ewmh->screens[screen_nbr]->root,
                                      ewmh->_NET_REQUEST_FRAME_EXTENTS, 0, NULL);
}

xcb_void_cookie_t xcb_ewmh_set_wm_name(xcb_ewmh_connection_t *ewmh,
                                       xcb_window_t window,
                                       uint32_t strings_len,
                                       const char *strings);

xcb_void_cookie_t xcb_ewmh_set_wm_name_checked(xcb_ewmh_connection_t *ewmh,
                                               xcb_window_t window,
                                               uint32_t strings_len,
                                               const char *strings);

xcb_get_property_cookie_t xcb_ewmh_get_wm_name_unchecked(xcb_ewmh_connection_t *ewmh,
                                                         xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_name(xcb_ewmh_connection_t *ewmh,
                                               xcb_window_t window);

static inline uint8_t
xcb_ewmh_get_wm_name_from_reply(xcb_ewmh_connection_t *ewmh,
                                xcb_ewmh_get_utf8_strings_reply_t *data,
                                xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_utf8_strings_from_reply(ewmh, data, r);
}

static inline uint8_t
xcb_ewmh_get_wm_name_reply(xcb_ewmh_connection_t *ewmh,
                           xcb_get_property_cookie_t cookie,
                           xcb_ewmh_get_utf8_strings_reply_t *data,
                           xcb_generic_error_t **e)
{
  return xcb_ewmh_get_utf8_strings_reply(ewmh, cookie, data, e);
}

xcb_void_cookie_t xcb_ewmh_set_wm_visible_name(xcb_ewmh_connection_t *ewmh,
                                               xcb_window_t window,
                                               uint32_t strings_len,
                                               const char *strings);

xcb_void_cookie_t xcb_ewmh_set_wm_visible_name_checked(xcb_ewmh_connection_t *ewmh,
                                                       xcb_window_t window,
                                                       uint32_t strings_len,
                                                       const char *strings);

xcb_get_property_cookie_t xcb_ewmh_get_wm_visible_name_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                 xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_visible_name(xcb_ewmh_connection_t *ewmh,
                                                       xcb_window_t window);

static inline uint8_t
xcb_ewmh_get_wm_visible_name_from_reply(xcb_ewmh_connection_t *ewmh,
                                        xcb_ewmh_get_utf8_strings_reply_t *data,
                                        xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_utf8_strings_from_reply(ewmh, data, r);
}

static inline uint8_t
xcb_ewmh_get_wm_visible_name_reply(xcb_ewmh_connection_t *ewmh,
                                   xcb_get_property_cookie_t cookie,
                                   xcb_ewmh_get_utf8_strings_reply_t *data,
                                   xcb_generic_error_t **e)
{
  return xcb_ewmh_get_utf8_strings_reply(ewmh, cookie, data, e);
}

xcb_void_cookie_t xcb_ewmh_set_wm_icon_name(xcb_ewmh_connection_t *ewmh,
                                            xcb_window_t window,
                                            uint32_t strings_len,
                                            const char *strings);

xcb_void_cookie_t xcb_ewmh_set_wm_icon_name_checked(xcb_ewmh_connection_t *ewmh,
                                                    xcb_window_t window,
                                                    uint32_t strings_len,
                                                    const char *strings);

xcb_get_property_cookie_t xcb_ewmh_get_wm_icon_name_unchecked(xcb_ewmh_connection_t *ewmh,
                                                              xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_icon_name(xcb_ewmh_connection_t *ewmh,
                                                    xcb_window_t window);

static inline uint8_t
xcb_ewmh_get_wm_icon_name_from_reply(xcb_ewmh_connection_t *ewmh,
                                     xcb_ewmh_get_utf8_strings_reply_t *data,
                                     xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_utf8_strings_from_reply(ewmh, data, r);
}

static inline uint8_t
xcb_ewmh_get_wm_icon_name_reply(xcb_ewmh_connection_t *ewmh,
                                xcb_get_property_cookie_t cookie,
                                xcb_ewmh_get_utf8_strings_reply_t *data,
                                xcb_generic_error_t **e)
{
  return xcb_ewmh_get_utf8_strings_reply(ewmh, cookie, data, e);
}

xcb_void_cookie_t xcb_ewmh_set_wm_visible_icon_name(xcb_ewmh_connection_t *ewmh,
                                                    xcb_window_t window,
                                                    uint32_t strings_len,
                                                    const char *strings);

xcb_void_cookie_t xcb_ewmh_set_wm_visible_icon_name_checked(xcb_ewmh_connection_t *ewmh,
                                                            xcb_window_t window,
                                                            uint32_t strings_len,
                                                            const char *strings);

xcb_get_property_cookie_t xcb_ewmh_get_wm_visible_icon_name_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                      xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_visible_icon_name(xcb_ewmh_connection_t *ewmh,
                                                            xcb_window_t window);

static inline uint8_t
xcb_ewmh_get_wm_visible_icon_name_from_reply(xcb_ewmh_connection_t *ewmh,
                                             xcb_ewmh_get_utf8_strings_reply_t *data,
                                             xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_utf8_strings_from_reply(ewmh, data, r);
}

static inline uint8_t
xcb_ewmh_get_wm_visible_icon_name_reply(xcb_ewmh_connection_t *ewmh,
                                        xcb_get_property_cookie_t cookie,
                                        xcb_ewmh_get_utf8_strings_reply_t *data,
                                        xcb_generic_error_t **e)
{
  return xcb_ewmh_get_utf8_strings_reply(ewmh, cookie, data, e);
}

xcb_void_cookie_t xcb_ewmh_set_wm_desktop(xcb_ewmh_connection_t *ewmh,
                                          xcb_window_t window,
                                          uint32_t desktop);

xcb_void_cookie_t xcb_ewmh_set_wm_desktop_checked(xcb_ewmh_connection_t *ewmh,
                                                  xcb_window_t window,
                                                  uint32_t desktop);


xcb_get_property_cookie_t xcb_ewmh_get_wm_desktop_unchecked(xcb_ewmh_connection_t *ewmh,
                                                            xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_desktop(xcb_ewmh_connection_t *ewmh,
                                                  xcb_window_t window);

static inline uint8_t
xcb_ewmh_get_wm_desktop_from_reply(uint32_t *desktop,
                                   xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_cardinal_from_reply(desktop, r);
}

static inline uint8_t
xcb_ewmh_get_wm_desktop_reply(xcb_ewmh_connection_t *ewmh,
                              xcb_get_property_cookie_t cookie,
                              uint32_t *desktop,
                              xcb_generic_error_t **e)
{
  return xcb_ewmh_get_cardinal_reply(ewmh, cookie, desktop, e);
}

xcb_void_cookie_t xcb_ewmh_request_change_wm_desktop(xcb_ewmh_connection_t *ewmh,
                                                     int screen_nbr,
                                                     xcb_window_t client_window,
                                                     uint32_t new_desktop,
                                                     xcb_ewmh_client_source_type_t source_indication);

xcb_void_cookie_t xcb_ewmh_set_wm_window_type(xcb_ewmh_connection_t *ewmh,
                                              xcb_window_t window,
                                              uint32_t list_len,
                                              xcb_atom_t *list);

xcb_void_cookie_t xcb_ewmh_set_wm_window_type_checked(xcb_ewmh_connection_t *ewmh,
                                                      xcb_window_t window,
                                                      uint32_t list_len,
                                                      xcb_atom_t *list);

xcb_get_property_cookie_t xcb_ewmh_get_wm_window_type_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_window_type(xcb_ewmh_connection_t *ewmh,
                                                      xcb_window_t window);

uint8_t xcb_ewmh_get_wm_window_type_from_reply(xcb_ewmh_get_atoms_reply_t *wtypes,
                                               xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_wm_window_type_reply(xcb_ewmh_connection_t *ewmh,
                                          xcb_get_property_cookie_t cookie,
                                          xcb_ewmh_get_atoms_reply_t *name,
                                          xcb_generic_error_t **e);

xcb_void_cookie_t xcb_ewmh_set_wm_state(xcb_ewmh_connection_t *ewmh,
                                        xcb_window_t window,
                                        uint32_t list_len,
                                        xcb_atom_t *list);

xcb_void_cookie_t xcb_ewmh_set_wm_state_checked(xcb_ewmh_connection_t *ewmh,
                                                xcb_window_t window,
                                                uint32_t list_len,
                                                xcb_atom_t *list);

xcb_get_property_cookie_t xcb_ewmh_get_wm_state_unchecked(xcb_ewmh_connection_t *ewmh,
                                                          xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_state(xcb_ewmh_connection_t *ewmh,
                                                xcb_window_t window);

uint8_t xcb_ewmh_get_wm_state_from_reply(xcb_ewmh_get_atoms_reply_t *wtypes,
                                         xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_wm_state_reply(xcb_ewmh_connection_t *ewmh,
                                    xcb_get_property_cookie_t cookie,
                                    xcb_ewmh_get_atoms_reply_t *name,
                                    xcb_generic_error_t **e);

xcb_void_cookie_t xcb_ewmh_request_change_wm_state(xcb_ewmh_connection_t *ewmh,
                                                   int screen_nbr,
                                                   xcb_window_t client_window,
                                                   xcb_ewmh_wm_state_action_t action,
                                                   xcb_atom_t first_property,
                                                   xcb_atom_t second_property,
                                                   xcb_ewmh_client_source_type_t source_indication);

xcb_void_cookie_t xcb_ewmh_set_wm_allowed_actions(xcb_ewmh_connection_t *ewmh,
                                                  xcb_window_t window,
                                                  uint32_t list_len,
                                                  xcb_atom_t *list);

xcb_void_cookie_t xcb_ewmh_set_wm_allowed_actions_checked(xcb_ewmh_connection_t *ewmh,
                                                          xcb_window_t window,
                                                          uint32_t list_len,
                                                          xcb_atom_t *list);

xcb_get_property_cookie_t xcb_ewmh_get_wm_allowed_actions_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                    xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_allowed_actions(xcb_ewmh_connection_t *ewmh,
                                                          xcb_window_t window);

uint8_t xcb_ewmh_get_wm_allowed_actions_from_reply(xcb_ewmh_get_atoms_reply_t *wtypes,
                                                   xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_wm_allowed_actions_reply(xcb_ewmh_connection_t *ewmh,
                                              xcb_get_property_cookie_t cookie,
                                              xcb_ewmh_get_atoms_reply_t *name,
                                              xcb_generic_error_t **e);

xcb_void_cookie_t xcb_ewmh_set_wm_strut(xcb_ewmh_connection_t *ewmh,
                                        xcb_window_t window,
                                        uint32_t left, uint32_t right,
                                        uint32_t top, uint32_t bottom);

xcb_void_cookie_t xcb_ewmh_set_wm_strut_checked(xcb_ewmh_connection_t *ewmh,
                                                xcb_window_t window,
                                                uint32_t left, uint32_t right,
                                                uint32_t top, uint32_t bottom);

xcb_get_property_cookie_t xcb_ewmh_get_wm_strut_unchecked(xcb_ewmh_connection_t *ewmh,
                                                          xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_strut(xcb_ewmh_connection_t *ewmh,
                                                xcb_window_t window);

uint8_t xcb_ewmh_get_wm_strut_from_reply(xcb_ewmh_get_extents_reply_t *struts,
                                         xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_wm_strut_reply(xcb_ewmh_connection_t *ewmh,
                                    xcb_get_property_cookie_t cookie,
                                    xcb_ewmh_get_extents_reply_t *struts,
                                    xcb_generic_error_t **e);

xcb_void_cookie_t xcb_ewmh_set_wm_strut_partial(xcb_ewmh_connection_t *ewmh,
                                                xcb_window_t window,
                                                xcb_ewmh_wm_strut_partial_t wm_strut);

xcb_void_cookie_t xcb_ewmh_set_wm_strut_partial_checked(xcb_ewmh_connection_t *ewmh,
                                                        xcb_window_t window,
                                                        xcb_ewmh_wm_strut_partial_t wm_strut);

xcb_get_property_cookie_t xcb_ewmh_get_wm_strut_partial_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                  xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_strut_partial(xcb_ewmh_connection_t *ewmh,
                                                        xcb_window_t window);

uint8_t xcb_ewmh_get_wm_strut_partial_from_reply(xcb_ewmh_wm_strut_partial_t *struts,
                                                 xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_wm_strut_partial_reply(xcb_ewmh_connection_t *ewmh,
                                            xcb_get_property_cookie_t cookie,
                                            xcb_ewmh_wm_strut_partial_t *struts,
                                            xcb_generic_error_t **e);

xcb_void_cookie_t xcb_ewmh_set_wm_icon_geometry(xcb_ewmh_connection_t *ewmh,
                                                xcb_window_t window,
                                                uint32_t left, uint32_t right,
                                                uint32_t top, uint32_t bottom);

xcb_void_cookie_t xcb_ewmh_set_wm_icon_geometry_checked(xcb_ewmh_connection_t *ewmh,
                                                        xcb_window_t window,
                                                        uint32_t left, uint32_t right,
                                                        uint32_t top, uint32_t bottom);

xcb_get_property_cookie_t xcb_ewmh_get_wm_icon_geometry_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                  xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_icon_geometry(xcb_ewmh_connection_t *ewmh,
                                                        xcb_window_t window);

uint8_t xcb_ewmh_get_wm_icon_geometry_from_reply(xcb_ewmh_geometry_t *icons,
                                                 xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_wm_icon_geometry_reply(xcb_ewmh_connection_t *ewmh,
                                            xcb_get_property_cookie_t cookie,
                                            xcb_ewmh_geometry_t *icons,
                                            xcb_generic_error_t **e);

/**
 * @brief Send ChangeProperty request to set _NET_WM_ICON window
 *        property. The given data is considered to be already
 *        prepared, namely that it is an array such as: WIDTH1,
 *        HEIGHT1, IMG1, WIDTH2, HEIGHT2, IMG2.
 *
 *        If you only want to add or append a single icon, you may
 *        consider using xcb_ewmh_append_wm_icon_checked which is far
 *        easier to use.
 *
 * _NET_WM_ICON CARDINAL[][2+n]/32
 *
 * @param ewmh The information relative to EWMH
 * @param mode ChangeProperty mode (xcb_prop_mode_t)
 * @param window The window to set the property on
 * @param data_len Length of the data
 * @param data The data
 */
static inline xcb_void_cookie_t
xcb_ewmh_set_wm_icon_checked(xcb_ewmh_connection_t *ewmh,
                             uint8_t mode,
                             xcb_window_t window,
                             uint32_t data_len, uint32_t *data)
{
  return xcb_change_property_checked(ewmh->connection, mode,
                                     window, ewmh->_NET_WM_ICON,
                                     XCB_ATOM_CARDINAL, 32, data_len, data);
}

/**
 * @see xcb_ewmh_set_wm_icon_checked
 */
static inline xcb_void_cookie_t
xcb_ewmh_set_wm_icon(xcb_ewmh_connection_t *ewmh,
                     uint8_t mode,
                     xcb_window_t window,
                     uint32_t data_len, uint32_t *data)
{
  return xcb_change_property(ewmh->connection, mode, window,
                             ewmh->_NET_WM_ICON, XCB_ATOM_CARDINAL, 32,
                             data_len, data);
}

xcb_void_cookie_t xcb_ewmh_append_wm_icon_checked(xcb_ewmh_connection_t *ewmh,
                                                  xcb_window_t window,
                                                  uint32_t width, uint32_t height,
                                                  uint32_t img_len, uint32_t *img);

xcb_void_cookie_t xcb_ewmh_append_wm_icon(xcb_ewmh_connection_t *ewmh,
                                          xcb_window_t window,
                                          uint32_t width, uint32_t height,
                                          uint32_t img_len, uint32_t *img);

xcb_get_property_cookie_t xcb_ewmh_get_wm_icon_unchecked(xcb_ewmh_connection_t *ewmh,
                                                         xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_icon(xcb_ewmh_connection_t *ewmh,
                                               xcb_window_t window);

uint8_t xcb_ewmh_get_wm_icon_from_reply(xcb_ewmh_get_wm_icon_reply_t *wm_icon,
                                        xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_wm_icon_reply(xcb_ewmh_connection_t *ewmh,
                                   xcb_get_property_cookie_t cookie,
                                   xcb_ewmh_get_wm_icon_reply_t *wm_icon,
                                   xcb_generic_error_t **e);

xcb_ewmh_wm_icon_iterator_t xcb_ewmh_get_wm_icon_iterator(const xcb_ewmh_get_wm_icon_reply_t *wm_icon);

unsigned int xcb_ewmh_get_wm_icon_length(const xcb_ewmh_get_wm_icon_reply_t *wm_icon);

void xcb_ewmh_get_wm_icon_next(xcb_ewmh_wm_icon_iterator_t *iterator);

void xcb_ewmh_get_wm_icon_reply_wipe(xcb_ewmh_get_wm_icon_reply_t *wm_icon);

xcb_void_cookie_t xcb_ewmh_set_wm_pid(xcb_ewmh_connection_t *ewmh,
                                      xcb_window_t window,
                                      uint32_t pid);

xcb_void_cookie_t xcb_ewmh_set_wm_pid_checked(xcb_ewmh_connection_t *ewmh,
                                              xcb_window_t window,
                                              uint32_t pid);

xcb_get_property_cookie_t xcb_ewmh_get_wm_pid_unchecked(xcb_ewmh_connection_t *ewmh,
                                                        xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_pid(xcb_ewmh_connection_t *ewmh,
                                              xcb_window_t window);

static inline uint8_t
xcb_ewmh_get_wm_pid_from_reply(uint32_t *pid,
                               xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_cardinal_from_reply(pid, r);
}

static inline uint8_t
xcb_ewmh_get_wm_pid_reply(xcb_ewmh_connection_t *ewmh,
                          xcb_get_property_cookie_t cookie,
                          uint32_t *pid,
                          xcb_generic_error_t **e)
{
  return xcb_ewmh_get_cardinal_reply(ewmh, cookie, pid, e);
}

xcb_void_cookie_t xcb_ewmh_set_wm_handled_icons(xcb_ewmh_connection_t *ewmh,
                                                xcb_window_t window,
                                                uint32_t handled_icons);

xcb_void_cookie_t xcb_ewmh_set_wm_handled_icons_checked(xcb_ewmh_connection_t *ewmh,
                                                        xcb_window_t window,
                                                        uint32_t handled_icons);

xcb_get_property_cookie_t xcb_ewmh_get_wm_handled_icons_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                  xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_handled_icons(xcb_ewmh_connection_t *ewmh,
                                                        xcb_window_t window);

static inline uint8_t
xcb_ewmh_get_wm_handled_icons_from_reply(uint32_t *handled_icons,
                                         xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_cardinal_from_reply(handled_icons, r);
}

static inline uint8_t
xcb_ewmh_get_wm_handled_icons_reply(xcb_ewmh_connection_t *ewmh,
                                    xcb_get_property_cookie_t cookie,
                                    uint32_t *handled_icons,
                                    xcb_generic_error_t **e)
{
  return xcb_ewmh_get_cardinal_reply(ewmh, cookie, handled_icons, e);
}

xcb_void_cookie_t xcb_ewmh_set_wm_user_time(xcb_ewmh_connection_t *ewmh,
                                            xcb_window_t window,
                                            uint32_t xtime);

xcb_void_cookie_t xcb_ewmh_set_wm_user_time_checked(xcb_ewmh_connection_t *ewmh,
                                                    xcb_window_t window,
                                                    uint32_t pid);

xcb_get_property_cookie_t xcb_ewmh_get_wm_user_time_unchecked(xcb_ewmh_connection_t *ewmh,
                                                              xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_user_time(xcb_ewmh_connection_t *ewmh,
                                                    xcb_window_t window);

static inline uint8_t
xcb_ewmh_get_wm_user_time_from_reply(uint32_t *xtime,
                                     xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_cardinal_from_reply(xtime, r);
}

static inline uint8_t
xcb_ewmh_get_wm_user_time_reply(xcb_ewmh_connection_t *ewmh,
                                xcb_get_property_cookie_t cookie,
                                uint32_t *xtime,
                                xcb_generic_error_t **e)
{
  return xcb_ewmh_get_cardinal_reply(ewmh, cookie, xtime, e);
}

xcb_void_cookie_t xcb_ewmh_set_wm_user_time_window(xcb_ewmh_connection_t *ewmh,
                                                   xcb_window_t window,
                                                   uint32_t xtime);

xcb_void_cookie_t xcb_ewmh_set_wm_user_time_window_checked(xcb_ewmh_connection_t *ewmh,
                                                           xcb_window_t window,
                                                           uint32_t pid);

xcb_get_property_cookie_t xcb_ewmh_get_wm_user_time_window_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                     xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_user_time_window(xcb_ewmh_connection_t *ewmh,
                                                           xcb_window_t window);

static inline uint8_t
xcb_ewmh_get_wm_user_time_window_from_reply(uint32_t *xtime,
                                            xcb_get_property_reply_t *r)
{
  return xcb_ewmh_get_cardinal_from_reply(xtime, r);
}

static inline uint8_t
xcb_ewmh_get_wm_user_time_window_reply(xcb_ewmh_connection_t *ewmh,
                                       xcb_get_property_cookie_t cookie,
                                       uint32_t *xtime,
                                       xcb_generic_error_t **e)
{
  return xcb_ewmh_get_cardinal_reply(ewmh, cookie, xtime, e);
}

xcb_void_cookie_t xcb_ewmh_set_frame_extents(xcb_ewmh_connection_t *ewmh,
                                             xcb_window_t window,
                                             uint32_t left, uint32_t right,
                                             uint32_t top, uint32_t bottom);

xcb_void_cookie_t xcb_ewmh_set_frame_extents_checked(xcb_ewmh_connection_t *ewmh,
                                                     xcb_window_t window,
                                                     uint32_t left, uint32_t right,
                                                     uint32_t top, uint32_t bottom);

xcb_get_property_cookie_t xcb_ewmh_get_frame_extents_unchecked(xcb_ewmh_connection_t *ewmh,
                                                               xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_frame_extents(xcb_ewmh_connection_t *ewmh,
                                                     xcb_window_t window);

uint8_t xcb_ewmh_get_frame_extents_from_reply(xcb_ewmh_get_extents_reply_t *frame_extents,
                                              xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_frame_extents_reply(xcb_ewmh_connection_t *ewmh,
                                         xcb_get_property_cookie_t cookie,
                                         xcb_ewmh_get_extents_reply_t *frame_extents,
                                         xcb_generic_error_t **e);

xcb_void_cookie_t xcb_ewmh_send_wm_ping(xcb_ewmh_connection_t *ewmh,
                                        xcb_window_t window,
                                        xcb_timestamp_t timestamp);

xcb_void_cookie_t xcb_ewmh_set_wm_sync_request_counter(xcb_ewmh_connection_t *ewmh,
                                                       xcb_window_t window,
                                                       xcb_atom_t wm_sync_request_counter_atom,
                                                       uint32_t low, uint32_t high);

xcb_void_cookie_t xcb_ewmh_set_wm_sync_request_counter_checked(xcb_ewmh_connection_t *ewmh,
                                                               xcb_window_t window,
                                                               xcb_atom_t wm_sync_request_counter_atom,
                                                               uint32_t low, uint32_t high);

xcb_get_property_cookie_t xcb_ewmh_get_wm_sync_request_counter_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                         xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_sync_request_counter(xcb_ewmh_connection_t *ewmh,
                                                               xcb_window_t window);

uint8_t xcb_ewmh_get_wm_sync_request_counter_from_reply(uint64_t *counter,
                                                        xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_wm_sync_request_counter_reply(xcb_ewmh_connection_t *ewmh,
                                                   xcb_get_property_cookie_t cookie,
                                                   uint64_t *counter,
                                                   xcb_generic_error_t **e);

xcb_void_cookie_t xcb_ewmh_send_wm_sync_request(xcb_ewmh_connection_t *ewmh,
                                                xcb_window_t window,
                                                xcb_atom_t wm_protocols_atom,
                                                xcb_atom_t wm_sync_request_atom,
                                                xcb_timestamp_t timestamp,
                                                uint64_t counter);

xcb_void_cookie_t xcb_ewmh_set_wm_fullscreen_monitors(xcb_ewmh_connection_t *ewmh,
                                                      xcb_window_t window,
                                                      uint32_t top, uint32_t bottom,
                                                      uint32_t left, uint32_t right);

xcb_void_cookie_t xcb_ewmh_set_wm_fullscreen_monitors_checked(xcb_ewmh_connection_t *ewmh,
                                                              xcb_window_t window,
                                                              uint32_t top, uint32_t bottom,
                                                              uint32_t left, uint32_t right);

xcb_get_property_cookie_t xcb_ewmh_get_wm_fullscreen_monitors_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                        xcb_window_t window);

xcb_get_property_cookie_t xcb_ewmh_get_wm_fullscreen_monitors(xcb_ewmh_connection_t *ewmh,
                                                              xcb_window_t window);

uint8_t xcb_ewmh_get_wm_fullscreen_monitors_from_reply(xcb_ewmh_get_wm_fullscreen_monitors_reply_t *monitors,
                                                       xcb_get_property_reply_t *r);

uint8_t xcb_ewmh_get_wm_fullscreen_monitors_reply(xcb_ewmh_connection_t *ewmh,
                                                  xcb_get_property_cookie_t cookie,
                                                  xcb_ewmh_get_wm_fullscreen_monitors_reply_t *monitors,
                                                  xcb_generic_error_t **e);


xcb_void_cookie_t xcb_ewmh_request_change_wm_fullscreen_monitors(xcb_ewmh_connection_t *ewmh,
                                                                 int screen_nbr,
                                                                 xcb_window_t window,
                                                                 uint32_t top, uint32_t bottom,
                                                                 uint32_t left, uint32_t right,
                                                                 xcb_ewmh_client_source_type_t source_indication);

/**
 * @brief Set _NET_WM_CM_Sn ownership to the given window
 *
 * For  each  screen they  manage,  compositing  manager MUST  acquire
 * ownership of a selection named _NET_WM_CM_Sn, where n is the screen
 * number.
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @param owner The new owner of _NET_WM_CM_Sn selection
 * @param timestamp The client's last user activity timestamp
 * @param selection_data1 Optional data described by ICCCM
 * @param selection_data2 Optional data described by ICCCM
 */
xcb_void_cookie_t xcb_ewmh_set_wm_cm_owner(xcb_ewmh_connection_t *ewmh,
                                           int screen_nbr,
                                           xcb_window_t owner,
                                           xcb_timestamp_t timestamp,
                                           uint32_t selection_data1,
                                           uint32_t selection_data2);

/**
 * @see xcb_ewmh_set_wm_cm_owner
 */
xcb_void_cookie_t xcb_ewmh_set_wm_cm_owner_checked(xcb_ewmh_connection_t *ewmh,
                                                   int screen_nbr,
                                                   xcb_window_t owner,
                                                   xcb_timestamp_t timestamp,
                                                   uint32_t selection_data1,
                                                   uint32_t selection_data2);

/**
 * @brief   Send  GetSelectOwner   request   to  get   the  owner   of
 *        _NET_WM_CM_Sn root window property
 *
 * @param ewmh The information relative to EWMH
 * @param screen_nbr The screen number
 * @return The _NET_WM_CM_Sn cookie of the GetSelectionOwner request
 */
xcb_get_selection_owner_cookie_t xcb_ewmh_get_wm_cm_owner_unchecked(xcb_ewmh_connection_t *ewmh,
                                                                    int screen_nbr);

/**
 * @see xcb_ewmh_get_wm_cm_owner_unchecked
 */
xcb_get_selection_owner_cookie_t xcb_ewmh_get_wm_cm_owner(xcb_ewmh_connection_t *ewmh,
                                                          int screen_nbr);

uint8_t xcb_ewmh_get_wm_cm_owner_from_reply(xcb_window_t *owner,
                                            xcb_get_selection_owner_reply_t *r);

/**
 * @brief Get reply from the GetProperty _NET_CLIENT_LIST cookie
 *
 * The  parameter  e  supplied  to  this  function  must  be  NULL  if
 * xcb_get_window_client_list_unchecked()  is   used.   Otherwise,  it
 * stores the error if any.
 *
 * @param ewmh The information relative to EWMH.
 * @param cookie The _NET_WM_CM_Sn GetSelectionOwner request cookie.
 * @param owner The window ID which owns the selection or None.
 * @param The xcb_generic_error_t supplied.
 * @return Return 1 on success, 0 otherwise.
 */
uint8_t xcb_ewmh_get_wm_cm_owner_reply(xcb_ewmh_connection_t *ewmh,
                                       xcb_get_selection_owner_cookie_t cookie,
                                       xcb_window_t *owner,
                                       xcb_generic_error_t **e);

#ifdef __cplusplus
}
#endif

/**
 * @}
 */

#endif /* __XCB_EWMH_H__ */
