# Copyright (C) 2014-2022 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import gdb

# Load the xmethods if GDB supports them.
def gdb_has_xmethods():
    try:
        import gdb.xmethod
        return True
    except ImportError:
        return False

def register_libstdcxx_printers(obj):
    # Load the pretty-printers.
    from .printers import register_libstdcxx_printers
    register_libstdcxx_printers(obj)

    if gdb_has_xmethods():
        from .xmethods import register_libstdcxx_xmethods
        register_libstdcxx_xmethods(obj)
