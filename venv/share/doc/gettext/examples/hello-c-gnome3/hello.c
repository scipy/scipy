/* Example for use of GNU gettext.
   This file is in the public domain.

   Source code of the C program.  */


/* Get GTK declarations.  */
#include <gtk/gtk.h>
#include <glib/gi18n.h>

/* Get getpid() declaration.  */
#if HAVE_UNISTD_H
# include <unistd.h>
#endif

#define UI_PATH "/org/gnu/gettext/examples/hello/hello.ui"
#define APPLICATION_ID "org.gnu.gettext.examples.hello"
#define GSETTINGS_SCHEMA "org.gnu.gettext.examples.hello"

/* Forward declaration of GObject types.  */

#define HELLO_TYPE_APPLICATION_WINDOW (hello_application_window_get_type ())
#define HELLO_APPLICATION_WINDOW(obj)                           \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj),                           \
                               HELLO_TYPE_APPLICATION_WINDOW,   \
                               HelloApplicationWindow))

typedef struct _HelloApplicationWindow HelloApplicationWindow;
typedef struct _HelloApplicationWindowClass HelloApplicationWindowClass;

#define HELLO_TYPE_APPLICATION (hello_application_get_type ())
#define HELLO_APPLICATION(obj)                          \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj),                   \
                               HELLO_TYPE_APPLICATION,  \
                               HelloApplication))

typedef struct _HelloApplication HelloApplication;
typedef struct _HelloApplicationClass HelloApplicationClass;

/* Custom application window implementation.  */

struct _HelloApplicationWindow
{
  GtkApplicationWindow parent;
  GtkWidget *label;
  GtkWidget *button;
  GSettings *settings;
  gsize label_id;
  gchar *labels[3];
};

struct _HelloApplicationWindowClass
{
  GtkApplicationWindowClass parent_class;
};

G_DEFINE_TYPE (HelloApplicationWindow, hello_application_window,
               GTK_TYPE_APPLICATION_WINDOW);

static void
update_content (HelloApplicationWindow *window)
{
  gtk_label_set_label (GTK_LABEL (window->label),
                       window->labels[window->label_id]);
  window->label_id = (window->label_id + 1) % G_N_ELEMENTS (window->labels);
}

static void
hello_application_window_init (HelloApplicationWindow *window)
{
  gtk_widget_init_template (GTK_WIDGET (window));

  window->settings = g_settings_new (GSETTINGS_SCHEMA);
  g_settings_bind (window->settings, "use-markup",
                   window->label, "use-markup",
                   G_SETTINGS_BIND_DEFAULT);

  window->labels[0]
    = g_strdup_printf (_("<big>Hello world!</big>\n"
                         "This program is running as "
                         "process number <b>%d</b>."),
                       getpid ());
  window->labels[1]
    = g_strdup (_("<big><u>This is another text</u></big>"));
  window->labels[2]
    = g_strdup (_("<big><i>This is yet another text</i></big>"));

  update_content (window);
}

static void
hello_application_window_dispose (GObject *object)
{
  HelloApplicationWindow *window = HELLO_APPLICATION_WINDOW (object);
  g_clear_object (&window->settings);
}

static void
hello_application_window_class_init (HelloApplicationWindowClass *klass)
{
  GObjectClass *gobject_class = G_OBJECT_CLASS (klass);

  gobject_class->dispose = hello_application_window_dispose;

  gtk_widget_class_set_template_from_resource (GTK_WIDGET_CLASS (klass),
                                               UI_PATH);
  gtk_widget_class_bind_template_child (GTK_WIDGET_CLASS (klass),
                                        HelloApplicationWindow, label);
  gtk_widget_class_bind_template_child (GTK_WIDGET_CLASS (klass),
                                        HelloApplicationWindow, button);
}

static HelloApplicationWindow *
hello_application_window_new (HelloApplication *application)
{
  return g_object_new (HELLO_TYPE_APPLICATION_WINDOW,
                       "application", application,
                       NULL);
}

/* Custom application implementation.  */

struct _HelloApplication
{
  GtkApplication parent;
};

struct _HelloApplicationClass
{
  GtkApplicationClass parent_class;
};

G_DEFINE_TYPE (HelloApplication, hello_application, GTK_TYPE_APPLICATION);

static void
hello_application_init (HelloApplication *application)
{
}

static void
clicked_callback (GtkWidget *widget, void *data)
{
  update_content (HELLO_APPLICATION_WINDOW (data));
}

static void
hello_application_activate (GApplication *application)
{
  HelloApplicationWindow *window;

  window = hello_application_window_new (HELLO_APPLICATION (application));
  g_signal_connect (window->button, "clicked",
                    G_CALLBACK (clicked_callback), window);
  gtk_window_present (GTK_WINDOW (window));
}

static void
hello_application_class_init (HelloApplicationClass *klass)
{
  G_APPLICATION_CLASS (klass)->activate = hello_application_activate;
}

static HelloApplication *
hello_application_new (void)
{
  return g_object_new (HELLO_TYPE_APPLICATION,
                       "application-id", APPLICATION_ID,
                       NULL);
}

int
main (int argc, char *argv[])
{
  GApplication *application;
  int status;

  /* Load the GSettings schema from the current directory.  */
  g_setenv ("GSETTINGS_SCHEMA_DIR", ".", FALSE);

  /* Initializations.  */
  textdomain ("hello-c-gnome3");
  bindtextdomain ("hello-c-gnome3", LOCALEDIR);

  /* Create application.  */
  application = G_APPLICATION (hello_application_new ());

  /* Start the application.  */
  status = g_application_run (application, argc, argv);
  g_object_unref (application);

  return status;
}
