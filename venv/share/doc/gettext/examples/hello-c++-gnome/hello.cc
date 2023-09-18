/* Example for use of GNU gettext.
   This file is in the public domain.

   Source code of the C++ program.  */


/* Get GNOME declarations.  */
#include <gnome.h>
#include <gtk--.h>

/* Get getpid() declaration.  */
#if HAVE_UNISTD_H
# include <unistd.h>
#endif

static Gtk::Main *application;

static gint
quit_callback (GdkEventAny*)
{
  application->quit ();
}

int
main (int argc, char *argv[])
{
  Gtk::Window *window;
  Gtk::VBox *panel;
  Gtk::Label *label1;
  Gtk::Alignment *label1aligned;
  Gtk::Label *label2;
  Gtk::Alignment *label2aligned;
  Gtk::Button *button;
  Gtk::HButtonBox *buttonbar;

  /* Initializations.  */

  setlocale (LC_ALL, "");
  application = new Gtk::Main (argc, argv);
  textdomain ("hello-c++-gnome");
  bindtextdomain ("hello-c++-gnome", LOCALEDIR);

  /* Create the GUI elements.  */

  window = new Gtk::Window (GTK_WINDOW_TOPLEVEL);
  window->set_title ("Hello example");
  window->realize ();
  window->delete_event.connect (SigC::slot (quit_callback));

  label1 = new Gtk::Label (_("Hello, world!"));

  label1aligned = new Gtk::Alignment (0.0, 0.5, 0, 0);
  label1aligned->add (*label1);

  label2 = new Gtk::Label (g_strdup_printf (_("This program is running as process number %d."), getpid ()));

  label2aligned = new Gtk::Alignment (0.0, 0.5, 0, 0);
  label2aligned->add (*label2);

  button = new Gtk::Button ("OK");
  button->clicked.connect (Gtk::Main::quit.slot()); //slot (quit_callback));

  buttonbar = new Gtk::HButtonBox (GTK_BUTTONBOX_END);
  buttonbar->pack_start (*button);

  panel = new Gtk::VBox (false, GNOME_PAD_SMALL);
  panel->pack_start (*label1aligned);
  panel->pack_start (*label2aligned);
  panel->pack_start (*buttonbar);

  window->add (*panel);

  /* Make the GUI elements visible.  */

  label1->show ();
  label1aligned->show ();
  label2->show ();
  label2aligned->show ();
  button->show ();
  buttonbar->show ();
  panel->show ();
  window->show ();

  /* Start the event loop.  */

  application->run ();
}
