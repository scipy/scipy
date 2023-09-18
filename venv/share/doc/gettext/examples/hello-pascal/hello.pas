{ Example for use of GNU gettext.
  This file is in the public domain.

  Source code of the Pascal program.  }

program hello;
{$mode delphi}

uses gettext,  { translateresourcestrings }
     sysutils; { format }

resourcestring
  hello_world = 'Hello, world!';
  running_as = 'This program is running as process number %d.';

begin
  translateresourcestrings({$i %LOCALEDIR%}+'/%s/LC_MESSAGES/hello-pascal.mo');
  writeln(hello_world);
  writeln(format(running_as,[GetProcessID]));
end.
