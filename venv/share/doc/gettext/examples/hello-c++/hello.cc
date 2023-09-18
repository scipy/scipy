// Example for use of GNU gettext.
// This file is in the public domain.

// Source code of the C++ program.


// Avoid deprecation warnings from g++ 3.1 or newer.
#if defined __GNUG__ && defined __DEPRECATED
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif

// Get setlocale() declaration.
#include <locale.h>

// Get getpid() declaration.
#if HAVE_UNISTD_H
# include <unistd.h>
#endif

// Get gettext(), textdomain(), bindtextdomain() declaration.
#include "gettext.h"
// Define shortcut for gettext().
#define _(string) gettext (string)

// Get autosprintf class declaration.
#include "autosprintf.h"
using gnu::autosprintf;

int
main ()
{
  setlocale (LC_ALL, "");
  textdomain ("hello-c++");
  bindtextdomain ("hello-c++", LOCALEDIR);

  cout << _("Hello, world!") << endl;
  cout << autosprintf (_("This program is running as process number %d."),
                       getpid ())
       << endl;
}
