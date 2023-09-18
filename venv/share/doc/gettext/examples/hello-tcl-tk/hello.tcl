#!@WISH@
# Example for use of GNU gettext.
# This file is in the public domain.
#
# Source code of the Tcl/Tk program.

package require msgcat
::msgcat::mcload [file join "@pkgdatadir@" "msgs"]
proc _ {s} {return [::msgcat::mc $s]}

frame .my
button .my.button \
  -text [_ "Hello, world!"] \
  -command exit
label .my.label \
  -text [format [_ "This program is running as process number %d."] [pid]]
pack .my.button -side top
pack .my.label -side bottom
pack .my
