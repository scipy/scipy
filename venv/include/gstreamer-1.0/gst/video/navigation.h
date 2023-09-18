/* GStreamer Navigation
 * Copyright (C) 2003 Ronald Bultje <rbultje@ronald.bitfreak.net>
 * Copyright (C) 2003 David A. Schleef <ds@schleef.org>
 *
 * navigation.h: navigation interface design
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

#ifndef __GST_NAVIGATION_H__
#define __GST_NAVIGATION_H__

#include <gst/gst.h>
#include <gst/video/video-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_NAVIGATION \
  (gst_navigation_get_type ())
#define GST_NAVIGATION(obj) \
    (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_NAVIGATION, GstNavigation))
#define GST_IS_NAVIGATION(obj) \
      (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_NAVIGATION))
#define GST_NAVIGATION_GET_INTERFACE(obj) \
    (G_TYPE_INSTANCE_GET_INTERFACE ((obj), GST_TYPE_NAVIGATION, GstNavigationInterface))

typedef struct _GstNavigation GstNavigation;
typedef struct _GstNavigationInterface GstNavigationInterface;

/**
 * GstNavigationModifierType:
 * @GST_NAVIGATION_MODIFIER_SHIFT_MASK: the Shift key.
 * @GST_NAVIGATION_MODIFIER_CONTROL_MASK: the Control key.
 * @GST_NAVIGATION_MODIFIER_MOD1_MASK: the third modifier key
 * @GST_NAVIGATION_MODIFIER_MOD2_MASK: the fourth modifier key
 * @GST_NAVIGATION_MODIFIER_MOD3_MASK: the fifth modifier key
 * @GST_NAVIGATION_MODIFIER_MOD4_MASK: the sixth modifier key
 * @GST_NAVIGATION_MODIFIER_MOD5_MASK: the seventh modifier key
 * @GST_NAVIGATION_MODIFIER_BUTTON1_MASK: the first mouse button (usually the left button).
 * @GST_NAVIGATION_MODIFIER_BUTTON2_MASK: the second mouse button (usually the right button).
 * @GST_NAVIGATION_MODIFIER_BUTTON3_MASK: the third mouse button (usually the mouse wheel button or middle button).
 * @GST_NAVIGATION_MODIFIER_BUTTON4_MASK: the fourth mouse button (typically the "Back" button).
 * @GST_NAVIGATION_MODIFIER_BUTTON5_MASK: the fifth mouse button (typically the "forward" button).
 * @GST_NAVIGATION_MODIFIER_SUPER_MASK: the Super modifier
 * @GST_NAVIGATION_MODIFIER_HYPER_MASK: the Hyper modifier
 * @GST_NAVIGATION_MODIFIER_META_MASK: the Meta modifier
 * @GST_NAVIGATION_MODIFIER_MASK: A mask covering all entries in #GdkModifierType.
 *
 * Flags to indicate the state of modifier keys and mouse buttons
 * in events.
 *
 * Typical modifier keys are Shift, Control, Meta, Super, Hyper, Alt, Compose,
 * Apple, CapsLock or ShiftLock.
 *
 * Since: 1.22
 */
typedef enum
{
  GST_NAVIGATION_MODIFIER_NONE          = 0,
  GST_NAVIGATION_MODIFIER_SHIFT_MASK    = 1 << 0,
  GST_NAVIGATION_MODIFIER_LOCK_MASK     = 1 << 1,
  GST_NAVIGATION_MODIFIER_CONTROL_MASK  = 1 << 2,

  GST_NAVIGATION_MODIFIER_MOD1_MASK  = 1 << 3,
  GST_NAVIGATION_MODIFIER_MOD2_MASK  = 1 << 4,
  GST_NAVIGATION_MODIFIER_MOD3_MASK  = 1 << 5,
  GST_NAVIGATION_MODIFIER_MOD4_MASK  = 1 << 6,
  GST_NAVIGATION_MODIFIER_MOD5_MASK  = 1 << 7,

  GST_NAVIGATION_MODIFIER_BUTTON1_MASK  = 1 << 8,
  GST_NAVIGATION_MODIFIER_BUTTON2_MASK  = 1 << 9,
  GST_NAVIGATION_MODIFIER_BUTTON3_MASK  = 1 << 10,
  GST_NAVIGATION_MODIFIER_BUTTON4_MASK  = 1 << 11,
  GST_NAVIGATION_MODIFIER_BUTTON5_MASK  = 1 << 12,

  GST_NAVIGATION_MODIFIER_SUPER_MASK    = 1 << 26,
  GST_NAVIGATION_MODIFIER_HYPER_MASK    = 1 << 27,
  GST_NAVIGATION_MODIFIER_META_MASK     = 1 << 28,

  GST_NAVIGATION_MODIFIER_MASK = (
    GST_NAVIGATION_MODIFIER_NONE          |
    GST_NAVIGATION_MODIFIER_SHIFT_MASK    |
    GST_NAVIGATION_MODIFIER_LOCK_MASK     |
    GST_NAVIGATION_MODIFIER_CONTROL_MASK  |
    GST_NAVIGATION_MODIFIER_MOD1_MASK      |
    GST_NAVIGATION_MODIFIER_MOD2_MASK      |
    GST_NAVIGATION_MODIFIER_MOD3_MASK      |
    GST_NAVIGATION_MODIFIER_MOD4_MASK      |
    GST_NAVIGATION_MODIFIER_MOD5_MASK      |
    GST_NAVIGATION_MODIFIER_BUTTON1_MASK  |
    GST_NAVIGATION_MODIFIER_BUTTON2_MASK  |
    GST_NAVIGATION_MODIFIER_BUTTON3_MASK  |
    GST_NAVIGATION_MODIFIER_BUTTON4_MASK  |
    GST_NAVIGATION_MODIFIER_BUTTON5_MASK  |
    GST_NAVIGATION_MODIFIER_SUPER_MASK    |
    GST_NAVIGATION_MODIFIER_HYPER_MASK    |
    GST_NAVIGATION_MODIFIER_META_MASK
  )

} GstNavigationModifierType;

/**
 * GstNavigationInterface:
 * @iface: the parent interface
 * @send_event: sending a navigation event
 * @send_event_simple: sending a navigation event (Since: 1.22)
 *
 * Navigation interface.
 */
struct _GstNavigationInterface {
  GTypeInterface iface;

  /* virtual functions */

  /**
   * GstNavigationInterface::send_event:
   *
   * sending a navigation event.
   *
   * Deprecated: 1.22: Use #GstNavigationInterface.send_event_simple() instead.
   */
  void (*send_event) (GstNavigation *navigation, GstStructure *structure);

  /**
   * GstNavigationInterface::send_event_simple:
   * @navigation: The navigation interface instance
   * @event: (transfer full): The event to send
   *
   * sending a navigation event.
   *
   * Since: 1.22
   */
  void (*send_event_simple) (GstNavigation *navigation, GstEvent *event);
};

GST_VIDEO_API
GType           gst_navigation_get_type (void);

/* Navigation commands */

/**
 * GstNavigationCommand:
 * @GST_NAVIGATION_COMMAND_INVALID: An invalid command entry
 * @GST_NAVIGATION_COMMAND_MENU1: Execute navigation menu command 1. For DVD,
 * this enters the DVD root menu, or exits back to the title from the menu.
 * @GST_NAVIGATION_COMMAND_MENU2: Execute navigation menu command 2. For DVD,
 * this jumps to the DVD title menu.
 * @GST_NAVIGATION_COMMAND_MENU3: Execute navigation menu command 3. For DVD,
 * this jumps into the DVD root menu.
 * @GST_NAVIGATION_COMMAND_MENU4: Execute navigation menu command 4. For DVD,
 * this jumps to the Subpicture menu.
 * @GST_NAVIGATION_COMMAND_MENU5: Execute navigation menu command 5. For DVD,
 * the jumps to the audio menu.
 * @GST_NAVIGATION_COMMAND_MENU6: Execute navigation menu command 6. For DVD,
 * this jumps to the angles menu.
 * @GST_NAVIGATION_COMMAND_MENU7: Execute navigation menu command 7. For DVD,
 * this jumps to the chapter menu.
 * @GST_NAVIGATION_COMMAND_LEFT: Select the next button to the left in a menu,
 * if such a button exists.
 * @GST_NAVIGATION_COMMAND_RIGHT: Select the next button to the right in a menu,
 * if such a button exists.
 * @GST_NAVIGATION_COMMAND_UP: Select the button above the current one in a
 * menu, if such a button exists.
 * @GST_NAVIGATION_COMMAND_DOWN: Select the button below the current one in a
 * menu, if such a button exists.
 * @GST_NAVIGATION_COMMAND_ACTIVATE: Activate (click) the currently selected
 * button in a menu, if such a button exists.
 * @GST_NAVIGATION_COMMAND_PREV_ANGLE: Switch to the previous angle in a
 * multiangle feature.
 * @GST_NAVIGATION_COMMAND_NEXT_ANGLE: Switch to the next angle in a multiangle
 * feature.
 *
 * A set of commands that may be issued to an element providing the
 * #GstNavigation interface. The available commands can be queried via
 * the gst_navigation_query_new_commands() query.
 *
 * For convenience in handling DVD navigation, the MENU commands are aliased as:
 *    GST_NAVIGATION_COMMAND_DVD_MENU            = @GST_NAVIGATION_COMMAND_MENU1
 *    GST_NAVIGATION_COMMAND_DVD_TITLE_MENU      = @GST_NAVIGATION_COMMAND_MENU2
 *    GST_NAVIGATION_COMMAND_DVD_ROOT_MENU       = @GST_NAVIGATION_COMMAND_MENU3
 *    GST_NAVIGATION_COMMAND_DVD_SUBPICTURE_MENU = @GST_NAVIGATION_COMMAND_MENU4
 *    GST_NAVIGATION_COMMAND_DVD_AUDIO_MENU      = @GST_NAVIGATION_COMMAND_MENU5
 *    GST_NAVIGATION_COMMAND_DVD_ANGLE_MENU      = @GST_NAVIGATION_COMMAND_MENU6
 *    GST_NAVIGATION_COMMAND_DVD_CHAPTER_MENU    = @GST_NAVIGATION_COMMAND_MENU7
 */
typedef enum {
  GST_NAVIGATION_COMMAND_INVALID  = 0,

  GST_NAVIGATION_COMMAND_MENU1    = 1,
  GST_NAVIGATION_COMMAND_MENU2    = 2,
  GST_NAVIGATION_COMMAND_MENU3    = 3,
  GST_NAVIGATION_COMMAND_MENU4    = 4,
  GST_NAVIGATION_COMMAND_MENU5    = 5,
  GST_NAVIGATION_COMMAND_MENU6    = 6,
  GST_NAVIGATION_COMMAND_MENU7    = 7,

  GST_NAVIGATION_COMMAND_LEFT     = 20,
  GST_NAVIGATION_COMMAND_RIGHT    = 21,
  GST_NAVIGATION_COMMAND_UP       = 22,
  GST_NAVIGATION_COMMAND_DOWN     = 23,
  GST_NAVIGATION_COMMAND_ACTIVATE = 24,

  GST_NAVIGATION_COMMAND_PREV_ANGLE = 30,
  GST_NAVIGATION_COMMAND_NEXT_ANGLE = 31
} GstNavigationCommand;

/* Some aliases for the menu command types */
#define GST_NAVIGATION_COMMAND_DVD_MENU            GST_NAVIGATION_COMMAND_MENU1
#define GST_NAVIGATION_COMMAND_DVD_TITLE_MENU      GST_NAVIGATION_COMMAND_MENU2
#define GST_NAVIGATION_COMMAND_DVD_ROOT_MENU       GST_NAVIGATION_COMMAND_MENU3
#define GST_NAVIGATION_COMMAND_DVD_SUBPICTURE_MENU GST_NAVIGATION_COMMAND_MENU4
#define GST_NAVIGATION_COMMAND_DVD_AUDIO_MENU      GST_NAVIGATION_COMMAND_MENU5
#define GST_NAVIGATION_COMMAND_DVD_ANGLE_MENU      GST_NAVIGATION_COMMAND_MENU6
#define GST_NAVIGATION_COMMAND_DVD_CHAPTER_MENU    GST_NAVIGATION_COMMAND_MENU7

/* Queries */
/**
 * GstNavigationQueryType:
 * @GST_NAVIGATION_QUERY_INVALID: invalid query
 * @GST_NAVIGATION_QUERY_COMMANDS: command query
 * @GST_NAVIGATION_QUERY_ANGLES: viewing angle query
 *
 * Types of navigation interface queries.
 */
typedef enum
{
  GST_NAVIGATION_QUERY_INVALID     = 0,
  GST_NAVIGATION_QUERY_COMMANDS    = 1,
  GST_NAVIGATION_QUERY_ANGLES      = 2
} GstNavigationQueryType;

GST_VIDEO_API
GstNavigationQueryType gst_navigation_query_get_type (GstQuery *query);

GST_VIDEO_API
GstQuery *      gst_navigation_query_new_commands       (void);

GST_VIDEO_API
void            gst_navigation_query_set_commands       (GstQuery *query, gint n_cmds, ...);

GST_VIDEO_API
void            gst_navigation_query_set_commandsv      (GstQuery *query, gint n_cmds,
                                                         GstNavigationCommand *cmds);

GST_VIDEO_API
gboolean        gst_navigation_query_parse_commands_length     (GstQuery *query,
                                                                guint *n_cmds);

GST_VIDEO_API
gboolean        gst_navigation_query_parse_commands_nth        (GstQuery *query, guint nth,
                                                                GstNavigationCommand *cmd);

GST_VIDEO_API
GstQuery *      gst_navigation_query_new_angles         (void);

GST_VIDEO_API
void            gst_navigation_query_set_angles         (GstQuery *query, guint cur_angle,
                                                         guint n_angles);

GST_VIDEO_API
gboolean        gst_navigation_query_parse_angles       (GstQuery *query, guint *cur_angle,
                                                         guint *n_angles);

/* Element messages */
/**
 * GstNavigationMessageType:
 * @GST_NAVIGATION_MESSAGE_INVALID: Returned from
 * gst_navigation_message_get_type() when the passed message is not a
 * navigation message.
 * @GST_NAVIGATION_MESSAGE_MOUSE_OVER: Sent when the mouse moves over or leaves a
 * clickable region of the output, such as a DVD menu button.
 * @GST_NAVIGATION_MESSAGE_COMMANDS_CHANGED: Sent when the set of available commands
 * changes and should re-queried by interested applications.
 * @GST_NAVIGATION_MESSAGE_ANGLES_CHANGED: Sent when display angles in a multi-angle
 * feature (such as a multiangle DVD) change - either angles have appeared or
 * disappeared.
 * @GST_NAVIGATION_MESSAGE_EVENT: Sent when a navigation event was not handled
 * by any element in the pipeline (Since: 1.6)
 *
 * A set of notifications that may be received on the bus when navigation
 * related status changes.
 */
typedef enum {
  GST_NAVIGATION_MESSAGE_INVALID,
  GST_NAVIGATION_MESSAGE_MOUSE_OVER,
  GST_NAVIGATION_MESSAGE_COMMANDS_CHANGED,
  GST_NAVIGATION_MESSAGE_ANGLES_CHANGED,
  GST_NAVIGATION_MESSAGE_EVENT
} GstNavigationMessageType;

GST_VIDEO_API
GstNavigationMessageType gst_navigation_message_get_type (GstMessage *message);

GST_VIDEO_API
GstMessage *    gst_navigation_message_new_mouse_over       (GstObject *src,
                                                             gboolean active);

GST_VIDEO_API
gboolean        gst_navigation_message_parse_mouse_over     (GstMessage *message,
                                                             gboolean *active);

GST_VIDEO_API
GstMessage *    gst_navigation_message_new_commands_changed (GstObject *src);

GST_VIDEO_API
GstMessage *    gst_navigation_message_new_angles_changed   (GstObject *src,
                                                             guint cur_angle,
                                                             guint n_angles);

GST_VIDEO_API
gboolean        gst_navigation_message_parse_angles_changed (GstMessage *message,
                                                             guint *cur_angle,
                                                             guint *n_angles);

GST_VIDEO_API
GstMessage *    gst_navigation_message_new_event            (GstObject *src,
							     GstEvent *event);

GST_VIDEO_API
gboolean        gst_navigation_message_parse_event          (GstMessage *message,
							     GstEvent ** event);
/* event parsing functions */
/**
 * GstNavigationEventType:
 * @GST_NAVIGATION_EVENT_INVALID: Returned from
 * gst_navigation_event_get_type() when the passed event is not a navigation event.
 * @GST_NAVIGATION_EVENT_KEY_PRESS: A key press event. Use
 * gst_navigation_event_parse_key_event() to extract the details from the event.
 * @GST_NAVIGATION_EVENT_KEY_RELEASE: A key release event. Use
 * gst_navigation_event_parse_key_event() to extract the details from the event.
 * @GST_NAVIGATION_EVENT_MOUSE_BUTTON_PRESS: A mouse button press event. Use
 * gst_navigation_event_parse_mouse_button_event() to extract the details from the
 * event.
 * @GST_NAVIGATION_EVENT_MOUSE_BUTTON_RELEASE: A mouse button release event. Use
 * gst_navigation_event_parse_mouse_button_event() to extract the details from the
 * event.
 * @GST_NAVIGATION_EVENT_MOUSE_MOVE: A mouse movement event. Use
 * gst_navigation_event_parse_mouse_move_event() to extract the details from the
 * event.
 * @GST_NAVIGATION_EVENT_COMMAND: A navigation command event. Use
 * gst_navigation_event_parse_command() to extract the details from the event.
 * @GST_NAVIGATION_EVENT_MOUSE_SCROLL: A mouse scroll event. Use
 * gst_navigation_event_parse_mouse_scroll_event() to extract the details from
 * the event. (Since: 1.18)
 * @GST_NAVIGATION_EVENT_TOUCH_DOWN: An event describing a new touch point,
 * which will be assigned an identifier that is unique to it for the duration
 * of its movement on the screen. Use gst_navigation_event_parse_touch_event()
 * to extract the details from the event. (Since: 1.22)
 * @GST_NAVIGATION_EVENT_TOUCH_MOTION: An event describing the movement of an
 * active touch point across the screen. Use
 * gst_navigation_event_parse_touch_event() to extract the details from the
 * event. (Since: 1.22)
 * @GST_NAVIGATION_EVENT_TOUCH_UP: An event describing a removed touch point.
 * After this event, its identifier may be reused for any new touch points. Use
 * gst_navigation_event_parse_touch_up_event() to extract the details from the
 * event. (Since: 1.22)
 * @GST_NAVIGATION_EVENT_TOUCH_FRAME: An event signaling the end of a sequence
 * of simultaneous touch events. (Since: 1.22)
 * @GST_NAVIGATION_EVENT_TOUCH_CANCEL: An event cancelling all currently active
 * touch points. (Since: 1.22)
 *
 * Enum values for the various events that an element implementing the
 * GstNavigation interface might send up the pipeline. Touch events have been
 * inspired by the libinput API, and have the same meaning here.
 */
typedef enum {
  GST_NAVIGATION_EVENT_INVALID                    = 0,
  GST_NAVIGATION_EVENT_KEY_PRESS                  = 1,
  GST_NAVIGATION_EVENT_KEY_RELEASE                = 2,
  GST_NAVIGATION_EVENT_MOUSE_BUTTON_PRESS         = 3,
  GST_NAVIGATION_EVENT_MOUSE_BUTTON_RELEASE       = 4,
  GST_NAVIGATION_EVENT_MOUSE_MOVE                 = 5,
  GST_NAVIGATION_EVENT_COMMAND                    = 6,

  /**
   * GST_NAVIGATION_EVENT_MOUSE_SCROLL:
   *
   * A mouse scroll event. Use gst_navigation_event_parse_mouse_scroll_event()
   * to extract the details from the event.
   *
   * Since: 1.18
   */
  GST_NAVIGATION_EVENT_MOUSE_SCROLL               = 7,

  /**
   * GST_NAVIGATION_EVENT_TOUCH_DOWN:
   *
   * An event describing a new touch point, which will be assigned an identifier
   * that is unique to it for the duration of its movement on the screen.
   * Use gst_navigation_event_parse_touch_event() to extract the details
   * from the event.
   *
   * Since: 1.22
   */
  GST_NAVIGATION_EVENT_TOUCH_DOWN                 = 8,

  /**
   * GST_NAVIGATION_EVENT_TOUCH_MOTION:
   *
   * An event describing the movement of an active touch point across
   * the screen. Use gst_navigation_event_parse_touch_event() to extract
   * the details from the event.
   *
   * Since: 1.22
   */
  GST_NAVIGATION_EVENT_TOUCH_MOTION               = 9,

  /**
   * GST_NAVIGATION_EVENT_TOUCH_UP:
   *
   * An event describing a removed touch point. After this event,
   * its identifier may be reused for any new touch points.
   * Use gst_navigation_event_parse_touch_up_event() to extract the details
   * from the event.
   *
   * Since: 1.22
   */
  GST_NAVIGATION_EVENT_TOUCH_UP                   = 10,

  /**
   * GST_NAVIGATION_EVENT_TOUCH_FRAME:
   *
   * An event signaling the end of a sequence of simultaneous touch events.
   *
   * Since: 1.22
   */
  GST_NAVIGATION_EVENT_TOUCH_FRAME                = 11,

  /**
   * GST_NAVIGATION_EVENT_TOUCH_CANCEL:
   *
   * An event cancelling all currently active touch points.
   *
   * Since: 1.22
   */
  GST_NAVIGATION_EVENT_TOUCH_CANCEL               = 12,
} GstNavigationEventType;

GST_VIDEO_API
GstNavigationEventType gst_navigation_event_get_type          (GstEvent *event);

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_key_press            (const gchar * key,
                                                               GstNavigationModifierType state) G_GNUC_MALLOC;

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_key_release          (const gchar * key,
                                                               GstNavigationModifierType state) G_GNUC_MALLOC;

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_mouse_button_press   (gint button, gdouble x,
                                                               gdouble y,
                                                               GstNavigationModifierType state) G_GNUC_MALLOC;

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_mouse_button_release (gint button, gdouble x,
                                                               gdouble y,
                                                               GstNavigationModifierType state) G_GNUC_MALLOC;

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_mouse_move           (gdouble x,
                                                               gdouble y,
                                                               GstNavigationModifierType state) G_GNUC_MALLOC;

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_mouse_scroll         (gdouble x, gdouble y,
                                                               gdouble delta_x, gdouble delta_y,
                                                               GstNavigationModifierType state) G_GNUC_MALLOC;

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_command              (GstNavigationCommand command) G_GNUC_MALLOC;

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_touch_down           (guint identifier,
                                                               gdouble x, gdouble y,
                                                               gdouble pressure,
                                                               GstNavigationModifierType state) G_GNUC_MALLOC;

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_touch_motion         (guint identifier,
                                                               gdouble x, gdouble y,
                                                               gdouble pressure,
                                                               GstNavigationModifierType state) G_GNUC_MALLOC;

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_touch_up             (guint identifier,
                                                               gdouble x, gdouble y,
                                                               GstNavigationModifierType state) G_GNUC_MALLOC;

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_touch_frame          (GstNavigationModifierType state) G_GNUC_MALLOC;

GST_VIDEO_API
GstEvent*       gst_navigation_event_new_touch_cancel         (GstNavigationModifierType state) G_GNUC_MALLOC;

GST_VIDEO_API
gboolean        gst_navigation_event_parse_key_event          (GstEvent *event,
                                                               const gchar **key);

GST_VIDEO_API
gboolean        gst_navigation_event_parse_mouse_button_event (GstEvent *event,
                                                               gint *button, gdouble *x, gdouble *y);

GST_VIDEO_API
gboolean        gst_navigation_event_parse_mouse_move_event   (GstEvent *event,
                                                               gdouble *x, gdouble *y);

GST_VIDEO_API
gboolean        gst_navigation_event_parse_mouse_scroll_event (GstEvent *event,
                                                               gdouble *x, gdouble *y,
                                                               gdouble *delta_x, gdouble *delta_y);

GST_VIDEO_API
gboolean        gst_navigation_event_parse_command            (GstEvent *event,
                                                               GstNavigationCommand *command);

GST_VIDEO_API
gboolean        gst_navigation_event_parse_touch_event        (GstEvent * event,
                                                               guint * identifier,
                                                               gdouble * x, gdouble * y,
                                                               gdouble * pressure);

GST_VIDEO_API
gboolean        gst_navigation_event_parse_touch_up_event     (GstEvent * event,
                                                               guint * identifier,
                                                               gdouble * x, gdouble * y);

GST_VIDEO_API
gboolean  gst_navigation_event_get_coordinates (GstEvent * event,
                                                gdouble * x, gdouble * y);

GST_VIDEO_API
gboolean  gst_navigation_event_set_coordinates (GstEvent * event,
                                                gdouble x, gdouble y);

/* interface virtual function wrappers */

GST_VIDEO_DEPRECATED_FOR(gst_navigation_send_event_simple)
void    gst_navigation_send_event        (GstNavigation *navigation,
                                          GstStructure *structure);

GST_VIDEO_DEPRECATED_FOR(gst_navigation_send_event_simple)
void    gst_navigation_send_key_event    (GstNavigation *navigation,
                                          const char *event, const char *key);

GST_VIDEO_DEPRECATED_FOR(gst_navigation_send_event_simple)
void    gst_navigation_send_mouse_event  (GstNavigation *navigation,
                                          const char *event, int button, double x, double y);

GST_VIDEO_DEPRECATED_FOR(gst_navigation_send_event_simple)
void    gst_navigation_send_mouse_scroll_event (GstNavigation *navigation,
                                                double x, double y, double delta_x, double delta_y);

GST_VIDEO_DEPRECATED_FOR(gst_navigation_send_event_simple)
void    gst_navigation_send_command      (GstNavigation *navigation,
                                          GstNavigationCommand command);

GST_VIDEO_API
void    gst_navigation_send_event_simple (GstNavigation *navigation,
                                          GstEvent *event);

GST_VIDEO_API
gboolean        gst_navigation_event_parse_modifier_state (GstEvent *event,
                                                           GstNavigationModifierType *state);

G_END_DECLS

#endif /* __GST_NAVIGATION_H__ */
