// Example for use of GNU gettext.
// This file is in the public domain.

// Source code of the C++ program.

#include <qapplication.h>
#include <qmainwindow.h>
#include <qlabel.h>
#include <qpushbutton.h>
#include <qstring.h>
#include <qvbox.h>
#include <qhbox.h>
#include <qtextcodec.h>

/* Get getpid() declaration.  */
#if HAVE_UNISTD_H
# include <unistd.h>
#endif

int
main (int argc, char *argv[])
{
  // Initializations.

  QApplication application (argc, argv);
#if 0
  GettextTranslator *translator =
    new GettextTranslator (&application, "hello-c++-qt", LOCALEDIR);
#else
  QTranslator *translator = new QTranslator (NULL);
  translator->load (QString ("hello-c++-qt") + "_" + QTextCodec::locale(),
                    PKGLOCALEDIR);
#endif
  application.installTranslator (translator);
#define _(string) application.translate ("", string)

  // Create the GUI elements.

  QMainWindow *window = new QMainWindow ();
  window->setCaption ("Hello example");

  QVBox *panel = new QVBox (window);
  panel->setSpacing (2);

  QLabel *label1 = new QLabel (_("Hello, world!"), panel);

  QString label2text;
  // NOT using QString::sprintf because it doesn't support reordering of
  // arguments.
  //label2text.sprintf (_("This program is running as process number %d"),
  //                    getpid ());
  label2text = _("This program is running as process number %1.").arg(getpid ());
  QLabel *label2 = new QLabel (label2text, panel);

  QHBox *buttonbar = new QHBox (panel);
  QWidget *filler = new QWidget (buttonbar); // makes the button right-aligned
  QPushButton *button = new QPushButton ("OK", buttonbar);
  button->setMaximumWidth (button->sizeHint().width() + 20);
  QObject::connect (button, SIGNAL (clicked ()), &application, SLOT (quit ()));

  panel->resize (panel->sizeHint ());
  window->resize (panel->frameSize ());

  application.setMainWidget (window);

  // Make the GUI elements visible.

  window->show ();

  // Start the event loop.

  return application.exec ();
}
