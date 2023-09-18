#!@TCLSH@
# Example for use of GNU gettext.
# This file is in the public domain.
#
# Source code of the Tcl program.

package require msgcat
::msgcat::mcload [file join "@pkgdatadir@" "msgs"]
proc _ {s} {return [::msgcat::mc $s]}

puts [_ "Hello, world!"]
puts [format [_ "This program is running as process number %d."] [pid]]
