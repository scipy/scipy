# -*- Mode: Python -*-

# GDBus - GLib D-Bus Library
#
# Copyright (C) 2008-2011 Red Hat, Inc.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General
# Public License along with this library; if not, see <http://www.gnu.org/licenses/>.
#
# Author: David Zeuthen <davidz@redhat.com>

import distutils.version
import os
import sys


# pylint: disable=too-few-public-methods
class Color:
    """ANSI Terminal colors"""

    GREEN = "\033[1;32m"
    BLUE = "\033[1;34m"
    YELLOW = "\033[1;33m"
    RED = "\033[1;31m"
    END = "\033[0m"


def print_color(msg, color=Color.END, prefix="MESSAGE"):
    """Print a string with a color prefix"""
    if os.isatty(sys.stderr.fileno()):
        real_prefix = "{start}{prefix}{end}".format(
            start=color, prefix=prefix, end=Color.END
        )
    else:
        real_prefix = prefix
    sys.stderr.write("{prefix}: {msg}\n".format(prefix=real_prefix, msg=msg))


def print_error(msg):
    """Print an error, and terminate"""
    print_color(msg, color=Color.RED, prefix="ERROR")
    sys.exit(1)


def print_warning(msg, fatal=False):
    """Print a warning, and optionally terminate"""
    if fatal:
        color = Color.RED
        prefix = "ERROR"
    else:
        color = Color.YELLOW
        prefix = "WARNING"
    print_color(msg, color, prefix)
    if fatal:
        sys.exit(1)


def print_info(msg):
    """Print a message"""
    print_color(msg, color=Color.GREEN, prefix="INFO")


def strip_dots(s):
    ret = ""
    force_upper = False
    for c in s:
        if c == ".":
            force_upper = True
        else:
            if force_upper:
                ret += c.upper()
                force_upper = False
            else:
                ret += c
    return ret


def dots_to_hyphens(s):
    return s.replace(".", "-")


def camel_case_to_uscore(s):
    ret = ""
    insert_uscore = False
    prev_was_lower = False
    initial = True
    for c in s:
        # Keep initial underscores in camel case
        if initial and c == "_":
            ret += "_"
            continue
        initial = False

        if c.isupper():
            if prev_was_lower:
                insert_uscore = True
            prev_was_lower = False
        else:
            prev_was_lower = True
        if insert_uscore:
            ret += "_"
        ret += c.lower()
        insert_uscore = False
    return ret


def uscore_to_camel_case(s):
    return "".join([s[0].upper() + s[1:].lower() if s else "_" for s in s.split("_")])


def is_ugly_case(s):
    if s and s.find("_") > 0:
        return True
    return False


def lookup_annotation(annotations, key):
    if annotations:
        for a in annotations:
            if a.key == key:
                return a.value
    return None


def lookup_docs(annotations):
    s = lookup_annotation(annotations, "org.gtk.GDBus.DocString")
    if s is None:
        return ""
    else:
        return s


def lookup_since(annotations):
    s = lookup_annotation(annotations, "org.gtk.GDBus.Since")
    if s is None:
        return ""
    else:
        return s


def lookup_brief_docs(annotations):
    s = lookup_annotation(annotations, "org.gtk.GDBus.DocString.Short")
    if s is None:
        return ""
    else:
        return s


def version_cmp_key(key):
    # If the 'since' version is 'UNRELEASED', compare higher than anything else
    # If it is empty put a 0 in its place as this will
    # allow LooseVersion to work and will always compare lower.
    if key[0] == "UNRELEASED":
        v = "9999"
    elif key[0]:
        v = str(key[0])
    else:
        v = "0"
    return (distutils.version.LooseVersion(v), key[1])
