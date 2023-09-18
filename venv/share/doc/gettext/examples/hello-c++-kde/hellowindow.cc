// Example for use of GNU gettext.
// Copyright (C) 2003 Free Software Foundation, Inc.
// This file is published under the GNU General Public License.

#if HAVE_CONFIG_H
# include <config.h>
#endif

/* Specification.  */
#include "hellowindow.h"

/* Declare i18n.  */
#include <klocale.h>
/* Declare KMainWindow.  */
#include <kmainwindow.h>
/* Declare QLabel.  */
#include <qlabel.h>
/* Declare QPushButton.  */
#include <qpushbutton.h>
/* Declare QString.  */
#include <qstring.h>
/* Declare QVBox.  */
#include <qvbox.h>
/* Declare QHBox.  */
#include <qhbox.h>

/* Get getpid() declaration.  */
#if HAVE_UNISTD_H
# include <unistd.h>
#endif

// The main window widget.

HelloMainWindow::HelloMainWindow (QWidget * parent, const char * name)
  : KMainWindow (parent, name)
{
  setCaption ("Hello example");

  QVBox *panel = new QVBox (this);
  panel->setSpacing (2);

  QLabel *label1 = new QLabel (i18n ("Hello, world!"), panel);

  QString label2text;
  // NOT using QString::sprintf because it doesn't support reordering of
  // arguments.
  //label2text.sprintf (i18n ("This program is running as process number %d"),
  //                    getpid ());
  label2text = i18n ("This program is running as process number %1.").arg(getpid ());
  QLabel *label2 = new QLabel (label2text, panel);

  QHBox *buttonbar = new QHBox (panel);
  QWidget *filler = new QWidget (buttonbar); // makes the button right-aligned
  button = new QPushButton ("OK", buttonbar);
  button->setMaximumWidth (button->sizeHint().width() + 20);

  panel->resize (panel->sizeHint ());
  resize (panel->frameSize ());
}

HelloMainWindow::~HelloMainWindow ()
{
}
