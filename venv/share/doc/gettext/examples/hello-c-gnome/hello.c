/* Example for use of GNU gettext.
   This file is in the public domain.

   Source code of the C program.  */


/* Get GNOME declarations.  */
#include <gnome.h>

/* Get getpid() declaration.  */
#if HAVE_UNISTD_H
# include <unistd.h>
#endif

static void
quit_callback (GtkWidget *widget, void *data)
{
  gtk_main_quit ();
}

int
main (int argc, char *argv[])
{
  GtkWidget *window;
  GtkWidget *panel;
  GtkWidget *label1;
  GtkWidget *label1aligned;
  GtkWidget *label2;
  GtkWidget *label2aligned;
  GtkWidget *button;
  GtkWidget *buttonbar;

  /* Initializations.  */

  gnome_init (PACKAGE, VERSION, argc, argv);
  textdomain ("hello-c-gnome");
  bindtextdomain ("hello-c-gnome", LOCALEDIR);

  /* Create the GUI elements.  */

  window = gnome_app_new ("hello-c-gnome", "Hello example");
  gtk_widget_realize (window);
  gtk_signal_connect (GTK_OBJECT (window), "delete_event",
                      GTK_SIGNAL_FUNC (quit_callback), NULL);

  label1 = gtk_label_new (_("Hello, world!"));

  label1aligned = gtk_alignment_new (0.0, 0.5, 0, 0);
  gtk_container_add (GTK_CONTAINER (label1aligned), label1);

  label2 = gtk_label_new (g_strdup_printf (_("This program is running as process number %d."), getpid ()));

  label2aligned = gtk_alignment_new (0.0, 0.5, 0, 0);
  gtk_container_add (GTK_CONTAINER (label2aligned), label2);

  button = gtk_button_new_with_label ("OK");
  gtk_signal_connect (GTK_OBJECT (button), "clicked",
                      GTK_SIGNAL_FUNC (quit_callback), NULL);

  buttonbar = gtk_hbutton_box_new ();
  gtk_button_box_set_layout (GTK_BUTTON_BOX (buttonbar), GTK_BUTTONBOX_END);
  gtk_box_pack_start_defaults (GTK_BOX (buttonbar), button);

  panel = gtk_vbox_new (FALSE, GNOME_PAD_SMALL);
  gtk_box_pack_start_defaults (GTK_BOX (panel), label1aligned);
  gtk_box_pack_start_defaults (GTK_BOX (panel), label2aligned);
  gtk_box_pack_start_defaults (GTK_BOX (panel), buttonbar);

  gnome_app_set_contents (GNOME_APP (window), panel);

  /* Make the GUI elements visible.  */

  gtk_widget_show (label1);
  gtk_widget_show (label1aligned);
  gtk_widget_show (label2);
  gtk_widget_show (label2aligned);
  gtk_widget_show (button);
  gtk_widget_show (buttonbar);
  gtk_widget_show (panel);
  gtk_widget_show (window);

  /* Start the event loop.  */

  gtk_main ();

  return 0;
}
