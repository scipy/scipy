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

from . import utils
from .utils import print_error


# See: variant_type_string_scan_internal()
def variant_type_string_scan(signature: str, depth_limit: int, i=0):
    beg_char = signature[i]
    i += 1
    if beg_char == "(":
        while signature[i] != ")":
            if depth_limit == 0:
                raise ValueError(
                    f'Bad signature "{signature}". Too much recursion beginning at {i}.'
                )
            i = variant_type_string_scan(signature, depth_limit - 1, i)
        i += 1
    elif beg_char == "{":
        if depth_limit == 0:
            raise ValueError(
                f'Bad signature "{signature}". Too much recursion beginning at {i}.'
            )
        elif signature[i] not in "bynqihuxtdsog?":
            raise ValueError(
                f'Bad signature "{signature}". "{signature[i]}" is not a valid type for dictionary keys at position {i}.'
            )
        i += 1
        i = variant_type_string_scan(signature, depth_limit - 1, i)
        if signature[i] != "}":
            raise ValueError(
                f'Bad signature "{signature}". Dict must end with "}}" at position {i}.'
            )
        i += 1
    elif beg_char == "a":
        if depth_limit == 0:
            raise ValueError(
                f'Bad signature "{signature}". Too much recursion beginning at {i}.'
            )
        i = variant_type_string_scan(signature, depth_limit - 1, i)
    elif beg_char not in "bynqiuxtdsogvr*?h":
        raise ValueError(
            f'Bad signature "{signature}". Unexpected value "{beg_char}" at position {i}.'
        )
    return i


# variant_check_signature() does not perform a strict validation check and
# should not be used in security sensitive contexts.
def variant_check_signature(signature: str):
    # See: gvariant-internal.h
    G_VARIANT_MAX_RECURSION_DEPTH = 128
    if len(signature) > 255:
        print_error("D-Bus maximum signature length of 255 exceeded.")
    for s in signature:
        if s not in "ybnqiuxthdvasog(){}":
            print_error(
                f'Bad signature "{signature}". "{s}" is not a valid D-Bus type.'
            )
    try:
        variant_type_string_scan(signature, G_VARIANT_MAX_RECURSION_DEPTH)
    except IndexError:
        print_error(
            f'Bad signature "{signature}". Error parsing string or brackets not closed.'
        )
    except ValueError as e:
        print_error(e.args[0])


class Annotation:
    def __init__(self, key, value):
        self.key = key
        self.value = value
        self.annotations = []
        self.since = ""

    def post_process(self, interface_prefix, cns, cns_upper, cns_lower, container):
        key = self.key
        overridden_key = utils.lookup_annotation(
            self.annotations, "org.gtk.GDBus.C.Name"
        )
        if utils.is_ugly_case(overridden_key):
            self.key_lower = overridden_key.lower()
        else:
            if overridden_key:
                key = overridden_key
            self.key_lower = (
                utils.camel_case_to_uscore(key)
                .lower()
                .replace("-", "_")
                .replace(".", "_")
            )

        if len(self.since) == 0:
            self.since = utils.lookup_since(self.annotations)
            if len(self.since) == 0:
                self.since = container.since

        for a in self.annotations:
            a.post_process(interface_prefix, cns, cns_upper, cns_lower, self)


class Arg:
    def __init__(self, name, signature):
        self.name = name
        self.signature = signature
        self.annotations = []
        self.doc_string = ""
        self.since = ""

    def post_process(self, interface_prefix, cns, cns_upper, cns_lower, arg_number):
        if len(self.doc_string) == 0:
            self.doc_string = utils.lookup_docs(self.annotations)
        if len(self.since) == 0:
            self.since = utils.lookup_since(self.annotations)

        if self.name is None:
            self.name = "unnamed_arg%d" % arg_number
        # default to GVariant
        self.ctype_in_g = "GVariant *"
        self.ctype_in = "GVariant *"
        self.ctype_in_dup = "GVariant *"
        self.ctype_in_default_value = "NULL"
        self.ctype_out = "GVariant **"
        self.gtype = "G_TYPE_VARIANT"
        self.free_func = "g_variant_unref"
        self.format_in = "@" + self.signature
        self.format_out = "@" + self.signature
        self.gvariant_get = "XXX"
        self.gvalue_type = "variant"
        self.gvalue_get = "g_marshal_value_peek_variant"
        self.gvalue_set = "g_value_take_variant"
        self.gclosure_marshaller = "g_cclosure_marshal_VOID__VARIANT"
        self.array_annotation = ""
        variant_check_signature(self.signature)

        if not utils.lookup_annotation(
            self.annotations, "org.gtk.GDBus.C.ForceGVariant"
        ):
            if self.signature == "b":
                self.ctype_in_g = "gboolean "
                self.ctype_in = "gboolean "
                self.ctype_in_default_value = "FALSE"
                self.ctype_out = "gboolean *"
                self.gtype = "G_TYPE_BOOLEAN"
                self.free_func = None
                self.format_in = "b"
                self.format_out = "b"
                self.gvariant_get = "g_variant_get_boolean"
                self.gvalue_type = "boolean"
                self.gvalue_get = "g_marshal_value_peek_boolean"
                self.gvalue_set = "g_value_set_boolean"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__BOOLEAN"
            elif self.signature == "y":
                self.ctype_in_g = "guchar "
                self.ctype_in = "guchar "
                self.ctype_in_default_value = "'\\0'"
                self.ctype_out = "guchar *"
                self.gtype = "G_TYPE_UCHAR"
                self.free_func = None
                self.format_in = "y"
                self.format_out = "y"
                self.gvariant_get = "g_variant_get_byte"
                self.gvalue_type = "uchar"
                self.gvalue_get = "g_marshal_value_peek_uchar"
                self.gvalue_set = "g_value_set_uchar"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__UCHAR"
            elif self.signature == "n":
                self.ctype_in_g = "gint "
                self.ctype_in = "gint16 "
                self.ctype_in_default_value = "0"
                self.ctype_out = "gint16 *"
                self.gtype = "G_TYPE_INT"
                self.free_func = None
                self.format_in = "n"
                self.format_out = "n"
                self.gvariant_get = "g_variant_get_int16"
                self.gvalue_type = "int"
                self.gvalue_get = "g_marshal_value_peek_int"
                self.gvalue_set = "g_value_set_int"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__INT"
            elif self.signature == "q":
                self.ctype_in_g = "guint "
                self.ctype_in = "guint16 "
                self.ctype_in_default_value = "0"
                self.ctype_out = "guint16 *"
                self.gtype = "G_TYPE_UINT"
                self.free_func = None
                self.format_in = "q"
                self.format_out = "q"
                self.gvariant_get = "g_variant_get_uint16"
                self.gvalue_type = "uint"
                self.gvalue_get = "g_marshal_value_peek_uint"
                self.gvalue_set = "g_value_set_uint"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__UINT"
            elif self.signature == "i":
                self.ctype_in_g = "gint "
                self.ctype_in = "gint "
                self.ctype_in_default_value = "0"
                self.ctype_out = "gint *"
                self.gtype = "G_TYPE_INT"
                self.free_func = None
                self.format_in = "i"
                self.format_out = "i"
                self.gvariant_get = "g_variant_get_int32"
                self.gvalue_type = "int"
                self.gvalue_get = "g_marshal_value_peek_int"
                self.gvalue_set = "g_value_set_int"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__INT"
            elif self.signature == "u":
                self.ctype_in_g = "guint "
                self.ctype_in = "guint "
                self.ctype_in_default_value = "0"
                self.ctype_out = "guint *"
                self.gtype = "G_TYPE_UINT"
                self.free_func = None
                self.format_in = "u"
                self.format_out = "u"
                self.gvariant_get = "g_variant_get_uint32"
                self.gvalue_type = "uint"
                self.gvalue_get = "g_marshal_value_peek_uint"
                self.gvalue_set = "g_value_set_uint"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__UINT"
            elif self.signature == "x":
                self.ctype_in_g = "gint64 "
                self.ctype_in = "gint64 "
                self.ctype_in_default_value = "0"
                self.ctype_out = "gint64 *"
                self.gtype = "G_TYPE_INT64"
                self.free_func = None
                self.format_in = "x"
                self.format_out = "x"
                self.gvariant_get = "g_variant_get_int64"
                self.gvalue_type = "int64"
                self.gvalue_get = "g_marshal_value_peek_int64"
                self.gvalue_set = "g_value_set_int64"
                self.gclosure_marshaller = None
            elif self.signature == "t":
                self.ctype_in_g = "guint64 "
                self.ctype_in = "guint64 "
                self.ctype_out = "guint64 *"
                self.ctype_in_default_value = "0"
                self.gtype = "G_TYPE_UINT64"
                self.free_func = None
                self.format_in = "t"
                self.format_out = "t"
                self.gvariant_get = "g_variant_get_uint64"
                self.gvalue_type = "uint64"
                self.gvalue_get = "g_marshal_value_peek_uint64"
                self.gvalue_set = "g_value_set_uint64"
                self.gclosure_marshaller = None
            elif self.signature == "d":
                self.ctype_in_g = "gdouble "
                self.ctype_in = "gdouble "
                self.ctype_in_default_value = "0.0"
                self.ctype_out = "gdouble *"
                self.gtype = "G_TYPE_DOUBLE"
                self.free_func = None
                self.format_in = "d"
                self.format_out = "d"
                self.gvariant_get = "g_variant_get_double"
                self.gvalue_type = "double"
                self.gvalue_get = "g_marshal_value_peek_double"
                self.gvalue_set = "g_value_set_double"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__DOUBLE"
            elif self.signature == "s":
                self.ctype_in_g = "const gchar *"
                self.ctype_in = "const gchar *"
                self.ctype_in_dup = "gchar *"
                self.ctype_in_default_value = "NULL"
                self.ctype_out = "gchar **"
                self.gtype = "G_TYPE_STRING"
                self.free_func = "g_free"
                self.format_in = "s"
                self.format_out = "s"
                self.gvariant_get = "g_variant_get_string"
                self.gvalue_type = "string"
                self.gvalue_get = "g_marshal_value_peek_string"
                self.gvalue_set = "g_value_set_string"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__STRING"
            elif self.signature == "o":
                self.ctype_in_g = "const gchar *"
                self.ctype_in = "const gchar *"
                self.ctype_in_dup = "gchar *"
                self.ctype_in_default_value = "NULL"
                self.ctype_out = "gchar **"
                self.gtype = "G_TYPE_STRING"
                self.free_func = "g_free"
                self.format_in = "o"
                self.format_out = "o"
                self.gvariant_get = "g_variant_get_string"
                self.gvalue_type = "string"
                self.gvalue_get = "g_marshal_value_peek_string"
                self.gvalue_set = "g_value_set_string"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__STRING"
            elif self.signature == "g":
                self.ctype_in_g = "const gchar *"
                self.ctype_in = "const gchar *"
                self.ctype_in_dup = "gchar *"
                self.ctype_in_default_value = "NULL"
                self.ctype_out = "gchar **"
                self.gtype = "G_TYPE_STRING"
                self.free_func = "g_free"
                self.format_in = "g"
                self.format_out = "g"
                self.gvariant_get = "g_variant_get_string"
                self.gvalue_type = "string"
                self.gvalue_get = "g_marshal_value_peek_string"
                self.gvalue_set = "g_value_set_string"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__STRING"
            elif self.signature == "ay":
                self.ctype_in_g = "const gchar *"
                self.ctype_in = "const gchar *"
                self.ctype_in_default_value = "NULL"
                self.ctype_in_dup = "gchar *"
                self.ctype_out = "gchar **"
                self.gtype = "G_TYPE_STRING"
                self.free_func = "g_free"
                self.format_in = "^ay"
                self.format_out = "^ay"
                self.gvariant_get = "g_variant_get_bytestring"
                self.gvalue_type = "string"
                self.gvalue_get = "g_marshal_value_peek_string"
                self.gvalue_set = "g_value_set_string"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__STRING"
            elif self.signature == "as":
                self.ctype_in_g = "const gchar *const *"
                self.ctype_in = "const gchar *const *"
                self.ctype_in_dup = "gchar **"
                self.ctype_in_default_value = "NULL"
                self.ctype_out = "gchar ***"
                self.gtype = "G_TYPE_STRV"
                self.free_func = "g_strfreev"
                self.format_in = "^as"
                self.format_out = "^as"
                self.gvariant_get = "g_variant_get_strv"
                self.gvalue_type = "boxed"
                self.gvalue_get = "g_marshal_value_peek_boxed"
                self.gvalue_set = "g_value_take_boxed"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__BOXED"
                self.array_annotation = "(array zero-terminated=1)"
            elif self.signature == "ao":
                self.ctype_in_g = "const gchar *const *"
                self.ctype_in = "const gchar *const *"
                self.ctype_in_dup = "gchar **"
                self.ctype_in_default_value = "NULL"
                self.ctype_out = "gchar ***"
                self.gtype = "G_TYPE_STRV"
                self.free_func = "g_strfreev"
                self.format_in = "^ao"
                self.format_out = "^ao"
                self.gvariant_get = "g_variant_get_objv"
                self.gvalue_type = "boxed"
                self.gvalue_get = "g_marshal_value_peek_boxed"
                self.gvalue_set = "g_value_take_boxed"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__BOXED"
                self.array_annotation = "(array zero-terminated=1)"
            elif self.signature == "aay":
                self.ctype_in_g = "const gchar *const *"
                self.ctype_in = "const gchar *const *"
                self.ctype_in_dup = "gchar **"
                self.ctype_in_default_value = "NULL"
                self.ctype_out = "gchar ***"
                self.gtype = "G_TYPE_STRV"
                self.free_func = "g_strfreev"
                self.format_in = "^aay"
                self.format_out = "^aay"
                self.gvariant_get = "g_variant_get_bytestring_array"
                self.gvalue_type = "boxed"
                self.gvalue_get = "g_marshal_value_peek_boxed"
                self.gvalue_set = "g_value_take_boxed"
                self.gclosure_marshaller = "g_cclosure_marshal_VOID__BOXED"
                self.array_annotation = "(array zero-terminated=1)"

        for a in self.annotations:
            a.post_process(interface_prefix, cns, cns_upper, cns_lower, self)


class Method:
    def __init__(self, name, h_type_implies_unix_fd=True):
        self.name = name
        self.h_type_implies_unix_fd = h_type_implies_unix_fd
        self.in_args = []
        self.out_args = []
        self.annotations = []
        self.doc_string = ""
        self.since = ""
        self.deprecated = False
        self.unix_fd = False

    def post_process(
        self, interface_prefix, cns, cns_upper, cns_lower, containing_iface
    ):
        if len(self.doc_string) == 0:
            self.doc_string = utils.lookup_docs(self.annotations)
        if len(self.since) == 0:
            self.since = utils.lookup_since(self.annotations)
            if len(self.since) == 0:
                self.since = containing_iface.since

        name = self.name
        overridden_name = utils.lookup_annotation(
            self.annotations, "org.gtk.GDBus.C.Name"
        )
        if utils.is_ugly_case(overridden_name):
            self.name_lower = overridden_name.lower()
        else:
            if overridden_name:
                name = overridden_name
            self.name_lower = utils.camel_case_to_uscore(name).lower().replace("-", "_")
        self.name_hyphen = self.name_lower.replace("_", "-")

        arg_count = 0
        for a in self.in_args:
            a.post_process(interface_prefix, cns, cns_upper, cns_lower, arg_count)
            arg_count += 1
            if self.h_type_implies_unix_fd and "h" in a.signature:
                self.unix_fd = True

        for a in self.out_args:
            a.post_process(interface_prefix, cns, cns_upper, cns_lower, arg_count)
            arg_count += 1
            if self.h_type_implies_unix_fd and "h" in a.signature:
                self.unix_fd = True

        if (
            utils.lookup_annotation(self.annotations, "org.freedesktop.DBus.Deprecated")
            == "true"
        ):
            self.deprecated = True

        if utils.lookup_annotation(self.annotations, "org.gtk.GDBus.C.UnixFD"):
            self.unix_fd = True

        self.marshaller_ret_arg = Arg("return", "b")
        self.marshaller_ret_arg.post_process(None, None, None, None, None)

        method_invocation_arg = Arg("method_invocation", None)
        method_invocation_arg.ctype_in = "GDBusMethodInvocation *"
        method_invocation_arg.gvalue_type = "object"
        method_invocation_arg.gvalue_get = "g_marshal_value_peek_object"
        method_invocation_arg.gclosure_marshaller = None
        self.marshaller_in_args = [method_invocation_arg] + self.in_args

        if self.unix_fd:
            fd_list_arg = Arg("fd_list", None)
            fd_list_arg.ctype_in = "GUnixFDList *"
            fd_list_arg.gvalue_type = "object"
            fd_list_arg.gvalue_get = "g_marshal_value_peek_object"
            fd_list_arg.gclosure_marshaller = None
            self.marshaller_in_args.insert(0, fd_list_arg)

        for a in self.annotations:
            a.post_process(interface_prefix, cns, cns_upper, cns_lower, self)


class Signal:
    def __init__(self, name):
        self.name = name
        self.args = []
        self.annotations = []
        self.doc_string = ""
        self.since = ""
        self.deprecated = False

    def post_process(
        self, interface_prefix, cns, cns_upper, cns_lower, containing_iface
    ):
        if len(self.doc_string) == 0:
            self.doc_string = utils.lookup_docs(self.annotations)
        if len(self.since) == 0:
            self.since = utils.lookup_since(self.annotations)
            if len(self.since) == 0:
                self.since = containing_iface.since

        name = self.name
        overridden_name = utils.lookup_annotation(
            self.annotations, "org.gtk.GDBus.C.Name"
        )
        if utils.is_ugly_case(overridden_name):
            self.name_lower = overridden_name.lower()
        else:
            if overridden_name:
                name = overridden_name
            self.name_lower = utils.camel_case_to_uscore(name).lower().replace("-", "_")
        self.name_upper = self.name_lower.upper()
        self.name_hyphen = self.name_lower.replace("_", "-")
        self.upper_id_name = "_".join(
            [cns_upper, containing_iface.name_upper, self.name_upper]
        )

        arg_count = 0
        for a in self.args:
            a.post_process(interface_prefix, cns, cns_upper, cns_lower, arg_count)
            arg_count += 1

        if (
            utils.lookup_annotation(self.annotations, "org.freedesktop.DBus.Deprecated")
            == "true"
        ):
            self.deprecated = True

        self.marshaller_ret_arg = None
        self.marshaller_in_args = self.args

        for a in self.annotations:
            a.post_process(interface_prefix, cns, cns_upper, cns_lower, self)


class Property:
    def __init__(self, name, signature, access):
        self.name = name
        self.signature = signature
        self.access = access
        self.annotations = []
        self.arg = Arg("value", self.signature)
        self.arg.annotations = self.annotations
        self.readable = False
        self.writable = False
        if self.access == "readwrite":
            self.readable = True
            self.writable = True
        elif self.access == "read":
            self.readable = True
        elif self.access == "write":
            self.writable = True
        else:
            print_error('Invalid access type "{}"'.format(self.access))
        self.doc_string = ""
        self.since = ""
        self.deprecated = False
        self.emits_changed_signal = True

    def post_process(
        self, interface_prefix, cns, cns_upper, cns_lower, containing_iface
    ):
        if len(self.doc_string) == 0:
            self.doc_string = utils.lookup_docs(self.annotations)
        if len(self.since) == 0:
            self.since = utils.lookup_since(self.annotations)
            if len(self.since) == 0:
                self.since = containing_iface.since

        name = self.name
        overridden_name = utils.lookup_annotation(
            self.annotations, "org.gtk.GDBus.C.Name"
        )
        if utils.is_ugly_case(overridden_name):
            self.name_lower = overridden_name.lower()
        else:
            if overridden_name:
                name = overridden_name
            self.name_lower = utils.camel_case_to_uscore(name).lower().replace("-", "_")
        self.name_hyphen = self.name_lower.replace("_", "-")
        # don't clash with the GType getter, e.g.:
        # GType foo_bar_get_type (void); G_GNUC_CONST
        if self.name_lower == "type":
            self.name_lower = "type_"

        # recalculate arg
        self.arg.annotations = self.annotations
        self.arg.post_process(interface_prefix, cns, cns_upper, cns_lower, 0)

        if (
            utils.lookup_annotation(self.annotations, "org.freedesktop.DBus.Deprecated")
            == "true"
        ):
            self.deprecated = True

        for a in self.annotations:
            a.post_process(interface_prefix, cns, cns_upper, cns_lower, self)

        # FIXME: for now we only support 'false' and 'const' on the signal itself,
        # see #674913 and
        # http://dbus.freedesktop.org/doc/dbus-specification.html#introspection-format
        # for details
        if utils.lookup_annotation(
            self.annotations, "org.freedesktop.DBus.Property.EmitsChangedSignal"
        ) in ("false", "const"):
            self.emits_changed_signal = False


class Interface:
    def __init__(self, name):
        self.name = name
        self.methods = []
        self.signals = []
        self.properties = []
        self.annotations = []
        self.doc_string = ""
        self.doc_string_brief = ""
        self.since = ""
        self.deprecated = False

    def post_process(self, interface_prefix, c_namespace):
        if len(self.doc_string) == 0:
            self.doc_string = utils.lookup_docs(self.annotations)
        if len(self.doc_string_brief) == 0:
            self.doc_string_brief = utils.lookup_brief_docs(self.annotations)
        if len(self.since) == 0:
            self.since = utils.lookup_since(self.annotations)

        if len(c_namespace) > 0:
            if utils.is_ugly_case(c_namespace):
                cns = c_namespace.replace("_", "")
                cns_upper = c_namespace.upper() + "_"
                cns_lower = c_namespace.lower() + "_"
            else:
                cns = c_namespace
                cns_upper = utils.camel_case_to_uscore(c_namespace).upper() + "_"
                cns_lower = utils.camel_case_to_uscore(c_namespace).lower() + "_"
        else:
            cns = ""
            cns_upper = ""
            cns_lower = ""

        overridden_name = utils.lookup_annotation(
            self.annotations, "org.gtk.GDBus.C.Name"
        )
        if utils.is_ugly_case(overridden_name):
            name = overridden_name.replace("_", "")
            name_with_ns = cns + name
            self.name_without_prefix = name
            self.camel_name = name_with_ns
            self.ns_upper = cns_upper
            self.name_lower = cns_lower + overridden_name.lower()
            self.name_upper = overridden_name.upper()

            # print_error('handle Ugly_Case "{}"'.format(overridden_name))
        else:
            if overridden_name:
                name = overridden_name
            else:
                name = self.name
                if name.startswith(interface_prefix):
                    name = name[len(interface_prefix) :]
            self.name_without_prefix = name
            name = utils.strip_dots(name)
            name_with_ns = utils.strip_dots(cns + "." + name)
            self.camel_name = name_with_ns
            self.ns_upper = cns_upper
            self.name_lower = cns_lower + utils.camel_case_to_uscore(name)
            self.name_upper = utils.camel_case_to_uscore(name).upper()

        self.name_hyphen = self.name_upper.lower().replace("_", "-")

        if (
            utils.lookup_annotation(self.annotations, "org.freedesktop.DBus.Deprecated")
            == "true"
        ):
            self.deprecated = True

        for m in self.methods:
            m.post_process(interface_prefix, cns, cns_upper, cns_lower, self)

        for s in self.signals:
            s.post_process(interface_prefix, cns, cns_upper, cns_lower, self)

        for p in self.properties:
            p.post_process(interface_prefix, cns, cns_upper, cns_lower, self)

        for a in self.annotations:
            a.post_process(interface_prefix, cns, cns_upper, cns_lower, self)

        if self.signals:
            self.signals_enum_name = "_".join([cns_upper, self.name_upper, "SIGNALS"])
