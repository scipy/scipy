/* Reports run time, in seconds, for a command.
   The command argument can have multiple words, but then
   it has to be quoted, as for example

      time-it "command < file1 > file2"

   The time interval resolution is one whole second.  */


#include <time.h>
int system ();
int printf ();

int
main (argv, argc)
     int argv;
     char **argc;
{
  time_t t0, t1;

  if (argv < 2)
    {
      printf ("Usage: time-it name_of_program_to_be_timed\n");
      exit (1);
    }
  time (&t0);
  /* Wait til the clock changes before starting.  */
  do
    {
      time (&t1);
    }
  while (t1 == t0);
  system (argc[1]);
  t0 = t1;
  time (&t1);
  printf ("%ld seconds.\n", t1 - t0);
  exit (0);
}
