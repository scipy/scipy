#!@PHP@ -q
<?php
  // Example for use of GNU gettext.
  // This file is in the public domain.
  //
  // Source code of the PHP program.

  setlocale (LC_ALL, "");
  textdomain ("hello-php");
  bindtextdomain ("hello-php", "@localedir@");

  echo _("Hello, world!");
  echo "\n";
  printf (_("This program is running as process number %d."), posix_getpid());
  echo "\n";
?>
