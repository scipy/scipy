/* Example for use of GNU gettext.
   This file is in the public domain.

   Source code of the C program.  */


/* Get setlocale() declaration.  */
#include <locale.h>

/* Get printf() declaration.  */
#include <stdio.h>

/* Get getpid() declaration.  */
#if HAVE_UNISTD_H
# include <unistd.h>
#endif

/* Get gettext(), textdomain(), bindtextdomain() declaration.  */
#include "gettext.h"
/* Define shortcut for gettext().  */
#define _(string) gettext (string)

int
main ()
{
  setlocale (LC_ALL, "");
  textdomain ("hello-c");
  bindtextdomain ("hello-c", LOCALEDIR);

  printf ("%s\n", _("Hello, world!"));
  printf (_("This program is running as process number %d."), getpid ());
  putchar ('\n');

  return 0;
}
