// Example for use of GNU gettext.
// Copyright (C) 2003 Free Software Foundation, Inc.
// This file is published under the GNU General Public License.

/* Declare KMainWindow.  */
#include <kmainwindow.h>
/* Declare QPushButton.  */
#include <qpushbutton.h>

// The main window widget.

class HelloMainWindow : public KMainWindow
{
  Q_OBJECT
public:
  HelloMainWindow (QWidget * parent = NULL, const char * name = NULL);
  ~HelloMainWindow ();
public:
  QPushButton *button;
};
