// Example for use of GNU gettext.
// Copyright (C) 2003, 2019 Free Software Foundation, Inc.
// This file is published under the GNU General Public License.

// Source code of the C++ program.

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
/* Declare KCmdLineArgs, KCmdLineOptions.  */
#include <kcmdlineargs.h>
/* Declare KApplication.  */
#include <kapplication.h>
/* Declare KAboutData.  */
#include <kaboutdata.h>
/* Declare main window widget.  */
#include "hellowindow.h"

// Comment line options.

static KCmdLineOptions options[] =
{
  { 0, 0, 0 } // End of options.
};

int
main (int argc, char *argv[])
{
  // Initializations.

  {
    // Add our installation directory to KDE's search list for message
    // catalogs.  By default it looks only in $KDEHOME/share/locale and
    // $KDEDIR/share/locale.
    QString kdedirs = getenv ("KDEDIRS");
    if (kdedirs.isEmpty ())
      kdedirs = PREFIX;
    else
      kdedirs = kdedirs + ":" + PREFIX;
    setenv ("KDEDIRS", (const char *) kdedirs.local8Bit(), true);
  }

  KAboutData aboutData ("hello-c++-kde",
                        I18N_NOOP ("Hello example"),
                        VERSION,
                        I18N_NOOP ("Hello world example"),
                        KAboutData::License_GPL,
                        "(C) 2003, 2019 Free Software Foundation",
                        NULL,
                        NULL,
                        "bug-gettext@gnu.org");
  KCmdLineArgs::init (argc, argv, &aboutData);
  KCmdLineArgs::addCmdLineOptions (options);
  KApplication application;

  // Create the GUI elements.

  HelloMainWindow *window = new HelloMainWindow ();
  QObject::connect (window->button, SIGNAL (clicked ()),
                    &application, SLOT (quit ()));

  application.setMainWidget (window);

  // Make the GUI elements visible.

  window->show ();

  // Start the event loop.

  return application.exec ();
}
