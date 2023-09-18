/* Example for use of GNU gettext.
   This file is in the public domain.

   Source code of the Objective C program.  */


/* Get GNOME declarations.  */
#include <obgnome/obgnome.h>

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
  Gnome_App *application;
  Gtk_Window *window;
  Gtk_VBox *panel;
  Gtk_Label *label1;
  Gtk_Alignment *label1aligned;
  Gtk_Label *label2;
  Gtk_Alignment *label2aligned;
  Gtk_Button *button;
  Gtk_ButtonBox *buttonbar;

  /* Initializations.  */

  application = [[Gnome_App alloc] initApp: PACKAGE : VERSION : argc : argv];
  textdomain ("hello-objc-gnome");
  bindtextdomain ("hello-objc-gnome", LOCALEDIR);

  /* Create the GUI elements.  */

  window = [[Gtk_Window alloc] initWithWindowInfo: GTK_WINDOW_TOPLEVEL];
  [window set_title: "Hello example"];
  [window realize];
  [window signal_connect: "delete_event" signalFunc: quit_callback funcData: NULL];

  label1 = [[Gtk_Label alloc] initWithLabelInfo: _("Hello, world!")];

  label1aligned = [[Gtk_Alignment alloc] initWithAlignmentInfo: 0.0 : 0.5 : 0 : 0];
  [label1aligned add: label1];

  label2 = [[Gtk_Label alloc] initWithLabelInfo: g_strdup_printf (_("This program is running as process number %d."), getpid ())];

  label2aligned =  [[Gtk_Alignment alloc] initWithAlignmentInfo: 0.0 : 0.5 : 0 : 0];
  [label2aligned add: label2];

  button = [Gtk_Button alloc];
  [button initWithLabel: "OK"];
  [button signal_connect: "clicked" signalFunc: quit_callback funcData: NULL];

  buttonbar = [Gtk_HButtonBox new];
  [buttonbar set_layout: GTK_BUTTONBOX_END];
  [buttonbar pack_start_defaults: button];

  panel = [[Gtk_VBox alloc] initWithVBoxInfo: FALSE : GNOME_PAD_SMALL];
  [panel pack_start_defaults: label1aligned];
  [panel pack_start_defaults: label2aligned];
  [panel pack_start_defaults: buttonbar];

  [window add: panel];

  /* Make the GUI elements visible.  */

  [label1 show];
  [label1aligned show];
  [label2 show];
  [label2aligned show];
  [button show];
  [buttonbar show];
  [panel show];
  [window show];

  /* Start the event loop.  */

  gtk_main ();

  return 0;
}
