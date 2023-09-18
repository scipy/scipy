/* GIO - GLib Input, Output and Streaming Library
 *
 * Copyright (C) 2006-2007 Red Hat, Inc.
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this library; if not, see <http://www.gnu.org/licenses/>.
 *
 * Author: Alexander Larsson <alexl@redhat.com>
 */

#ifndef __G_DESKTOP_APP_INFO_H__
#define __G_DESKTOP_APP_INFO_H__

#include <gio/gio.h>

G_BEGIN_DECLS

#define G_TYPE_DESKTOP_APP_INFO         (g_desktop_app_info_get_type ())
#define G_DESKTOP_APP_INFO(o)           (G_TYPE_CHECK_INSTANCE_CAST ((o), G_TYPE_DESKTOP_APP_INFO, GDesktopAppInfo))
#define G_DESKTOP_APP_INFO_CLASS(k)     (G_TYPE_CHECK_CLASS_CAST((k), G_TYPE_DESKTOP_APP_INFO, GDesktopAppInfoClass))
#define G_IS_DESKTOP_APP_INFO(o)        (G_TYPE_CHECK_INSTANCE_TYPE ((o), G_TYPE_DESKTOP_APP_INFO))
#define G_IS_DESKTOP_APP_INFO_CLASS(k)  (G_TYPE_CHECK_CLASS_TYPE ((k), G_TYPE_DESKTOP_APP_INFO))
#define G_DESKTOP_APP_INFO_GET_CLASS(o) (G_TYPE_INSTANCE_GET_CLASS ((o), G_TYPE_DESKTOP_APP_INFO, GDesktopAppInfoClass))

typedef struct _GDesktopAppInfo        GDesktopAppInfo;
typedef struct _GDesktopAppInfoClass   GDesktopAppInfoClass;

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GDesktopAppInfo, g_object_unref)

struct _GDesktopAppInfoClass
{
  GObjectClass parent_class;
};


GIO_AVAILABLE_IN_ALL
GType            g_desktop_app_info_get_type          (void) G_GNUC_CONST;

GIO_AVAILABLE_IN_ALL
GDesktopAppInfo *g_desktop_app_info_new_from_filename (const char      *filename);
GIO_AVAILABLE_IN_ALL
GDesktopAppInfo *g_desktop_app_info_new_from_keyfile  (GKeyFile        *key_file);

GIO_AVAILABLE_IN_ALL
const char *     g_desktop_app_info_get_filename      (GDesktopAppInfo *info);

GIO_AVAILABLE_IN_2_30
const char *     g_desktop_app_info_get_generic_name  (GDesktopAppInfo *info);
GIO_AVAILABLE_IN_2_30
const char *     g_desktop_app_info_get_categories    (GDesktopAppInfo *info);
GIO_AVAILABLE_IN_2_30
const char * const *g_desktop_app_info_get_keywords   (GDesktopAppInfo *info);
GIO_AVAILABLE_IN_2_30
gboolean         g_desktop_app_info_get_nodisplay     (GDesktopAppInfo *info);
GIO_AVAILABLE_IN_2_30
gboolean         g_desktop_app_info_get_show_in       (GDesktopAppInfo *info,
                                                       const gchar     *desktop_env);
GIO_AVAILABLE_IN_2_34
const char *     g_desktop_app_info_get_startup_wm_class (GDesktopAppInfo *info);

GIO_AVAILABLE_IN_ALL
GDesktopAppInfo *g_desktop_app_info_new               (const char      *desktop_id);
GIO_AVAILABLE_IN_ALL
gboolean         g_desktop_app_info_get_is_hidden     (GDesktopAppInfo *info);

GIO_DEPRECATED_IN_2_42
void             g_desktop_app_info_set_desktop_env   (const char      *desktop_env);

GIO_AVAILABLE_IN_2_36
gboolean         g_desktop_app_info_has_key           (GDesktopAppInfo *info,
                                                       const char      *key);
GIO_AVAILABLE_IN_2_36
char *           g_desktop_app_info_get_string        (GDesktopAppInfo *info,
                                                       const char      *key);
GIO_AVAILABLE_IN_2_56
char *           g_desktop_app_info_get_locale_string (GDesktopAppInfo *info,
                                                       const char      *key);
GIO_AVAILABLE_IN_2_36
gboolean         g_desktop_app_info_get_boolean       (GDesktopAppInfo *info,
                                                       const char      *key);

GIO_AVAILABLE_IN_2_60
gchar **         g_desktop_app_info_get_string_list (GDesktopAppInfo *info,
                                                     const char      *key,
                                                     gsize           *length);

GIO_AVAILABLE_IN_2_38
const gchar * const *   g_desktop_app_info_list_actions                 (GDesktopAppInfo   *info);

GIO_AVAILABLE_IN_2_38
void                    g_desktop_app_info_launch_action                (GDesktopAppInfo   *info,
                                                                         const gchar       *action_name,
                                                                         GAppLaunchContext *launch_context);

GIO_AVAILABLE_IN_2_38
gchar *                 g_desktop_app_info_get_action_name              (GDesktopAppInfo   *info,
                                                                         const gchar       *action_name);

#define G_TYPE_DESKTOP_APP_INFO_LOOKUP           (g_desktop_app_info_lookup_get_type ()) GIO_DEPRECATED_MACRO_IN_2_28
#define G_DESKTOP_APP_INFO_LOOKUP(obj)           (G_TYPE_CHECK_INSTANCE_CAST ((obj), G_TYPE_DESKTOP_APP_INFO_LOOKUP, GDesktopAppInfoLookup)) GIO_DEPRECATED_MACRO_IN_2_28
#define G_IS_DESKTOP_APP_INFO_LOOKUP(obj)	 (G_TYPE_CHECK_INSTANCE_TYPE ((obj), G_TYPE_DESKTOP_APP_INFO_LOOKUP)) GIO_DEPRECATED_MACRO_IN_2_28
#define G_DESKTOP_APP_INFO_LOOKUP_GET_IFACE(obj) (G_TYPE_INSTANCE_GET_INTERFACE ((obj), G_TYPE_DESKTOP_APP_INFO_LOOKUP, GDesktopAppInfoLookupIface)) GIO_DEPRECATED_MACRO_IN_2_28

/**
 * G_DESKTOP_APP_INFO_LOOKUP_EXTENSION_POINT_NAME:
 *
 * Extension point for default handler to URI association. See
 * [Extending GIO][extending-gio].
 *
 * Deprecated: 2.28: The #GDesktopAppInfoLookup interface is deprecated and
 *    unused by GIO.
 */
#define G_DESKTOP_APP_INFO_LOOKUP_EXTENSION_POINT_NAME "gio-desktop-app-info-lookup" GIO_DEPRECATED_MACRO_IN_2_28

/**
 * GDesktopAppInfoLookupIface:
 * @get_default_for_uri_scheme: Virtual method for
 *  g_desktop_app_info_lookup_get_default_for_uri_scheme().
 *
 * Interface that is used by backends to associate default
 * handlers with URI schemes.
 */
typedef struct _GDesktopAppInfoLookup GDesktopAppInfoLookup;
typedef struct _GDesktopAppInfoLookupIface GDesktopAppInfoLookupIface;

struct _GDesktopAppInfoLookupIface
{
  GTypeInterface g_iface;

  GAppInfo * (* get_default_for_uri_scheme) (GDesktopAppInfoLookup *lookup,
                                             const char            *uri_scheme);
};

GIO_DEPRECATED
GType     g_desktop_app_info_lookup_get_type                   (void) G_GNUC_CONST;

GIO_DEPRECATED
GAppInfo *g_desktop_app_info_lookup_get_default_for_uri_scheme (GDesktopAppInfoLookup *lookup,
                                                                const char            *uri_scheme);

/**
 * GDesktopAppLaunchCallback:
 * @appinfo: a #GDesktopAppInfo
 * @pid: Process identifier
 * @user_data: User data
 *
 * During invocation, g_desktop_app_info_launch_uris_as_manager() may
 * create one or more child processes.  This callback is invoked once
 * for each, providing the process ID.
 */
typedef void (*GDesktopAppLaunchCallback) (GDesktopAppInfo  *appinfo,
					   GPid              pid,
					   gpointer          user_data);

GIO_AVAILABLE_IN_2_28
gboolean    g_desktop_app_info_launch_uris_as_manager (GDesktopAppInfo            *appinfo,
						       GList                      *uris,
						       GAppLaunchContext          *launch_context,
						       GSpawnFlags                 spawn_flags,
						       GSpawnChildSetupFunc        user_setup,
						       gpointer                    user_setup_data,
						       GDesktopAppLaunchCallback   pid_callback,
						       gpointer                    pid_callback_data,
						       GError                    **error);

GIO_AVAILABLE_IN_2_58
gboolean    g_desktop_app_info_launch_uris_as_manager_with_fds (GDesktopAppInfo            *appinfo,
								GList                      *uris,
								GAppLaunchContext          *launch_context,
								GSpawnFlags                 spawn_flags,
								GSpawnChildSetupFunc        user_setup,
								gpointer                    user_setup_data,
								GDesktopAppLaunchCallback   pid_callback,
								gpointer                    pid_callback_data,
								gint                        stdin_fd,
								gint                        stdout_fd,
								gint                        stderr_fd,
								GError                    **error);

GIO_AVAILABLE_IN_2_40
gchar *** g_desktop_app_info_search (const gchar *search_string);

GIO_AVAILABLE_IN_2_42
GList *g_desktop_app_info_get_implementations (const gchar *interface);

G_END_DECLS

#endif /* __G_DESKTOP_APP_INFO_H__ */
