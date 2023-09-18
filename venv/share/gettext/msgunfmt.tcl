# Reading tcl/msgcat .msg files.
# Copyright (C) 2002 Free Software Foundation, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

namespace eval msgcat {
  namespace export mcset mcdump
  variable header ""
}

proc msgcat::puts_po_string {str} {
  # Replace " with \"
  regsub -all "\"" $str "\\\"" str
  # Replace \ with \\
  regsub -all "\\\\" $str "\\\\\\" str
  # Replace newline with \n
  regsub -all [subst "\n"] $str "\\n" str
  regsub -all [subst "\a"] $str "\\a" str
  regsub -all [subst "\b"] $str "\\b" str
  regsub -all [subst "\f"] $str "\\f" str
  regsub -all [subst "\r"] $str "\\r" str
  regsub -all [subst "\t"] $str "\\t" str
  regsub -all [subst "\v"] $str "\\v" str
  # Output it.
  puts -nonewline "\"$str\""
}

proc msgcat::write_po_message {msgid msgstr} {
  puts -nonewline "msgid "
  puts_po_string $msgid
  puts ""
  puts -nonewline "msgstr "
  puts_po_string $msgstr
  puts ""
  puts ""
}

# This gets called once for each message in the .msg catalog.
proc msgcat::mcset {locale src {dest ""}} {
  msgcat::write_po_message $src $dest
}

# Main function.
proc msgcat::mcdump {langfile} {
  if {[file exists $langfile]} {
    # msgunfmt expects the output in UTF-8 encoding.
    fconfigure stdout -encoding utf-8

    set msgcat::header ""

    set fd [open $langfile r]
    # In newer tcl versions, the .msg files are in UTF-8 encoding.
    fconfigure $fd -encoding utf-8
    eval [read $fd]
    close $fd

    if {$msgcat::header == ""} {
      # Provide a minimal header.
      set msgcat::header [subst "MIME-Version: 1.0\nContent-Type: text/plain; charset=UTF-8\nContent-Transfer-Encoding: 8bit\n"]
    }
    msgcat::write_po_message "" $msgcat::header
  } else {
    # Tell msgunfmt to emit an internationalized error message.
    exit 2
  }
}

# Main code: call the main function on the first and only argument.
msgcat::mcdump [lindex $argv 0]

exit 0
