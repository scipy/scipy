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

import re
import textwrap
from os import path

from . import utils


# ----------------------------------------------------------------------------------------------------


class DocbookCodeGenerator:
    def __init__(self, ifaces):
        self.ifaces = ifaces
        self.generate_expand_dicts()

    def print_method_prototype(self, i, m, in_synopsis):
        max_method_len = 0
        if in_synopsis:
            for _m in i.methods:
                max_method_len = max(len(_m.name), max_method_len)
        else:
            max_method_len = max(len(m.name), max_method_len)

        max_signature_len = 0
        if in_synopsis:
            for _m in i.methods:
                for a in _m.in_args:
                    max_signature_len = max(len(a.signature), max_signature_len)
                for a in _m.out_args:
                    max_signature_len = max(len(a.signature), max_signature_len)
        else:
            for a in m.in_args:
                max_signature_len = max(len(a.signature), max_signature_len)
            for a in m.out_args:
                max_signature_len = max(len(a.signature), max_signature_len)

        if in_synopsis:
            self.out.write(
                '<link linkend="gdbus-method-%s.%s">%s</link>%*s ('
                % (
                    utils.dots_to_hyphens(i.name),
                    m.name,
                    m.name,
                    max_method_len - len(m.name),
                    "",
                )
            )
        else:
            self.out.write("%s%*s (" % (m.name, max_method_len - len(m.name), ""))
        count = 0
        for a in m.in_args:
            if count > 0:
                self.out.write(",\n%*s" % (max_method_len + 2, ""))
            self.out.write(
                "IN  %s%*s %s"
                % (a.signature, max_signature_len - len(a.signature), "", a.name)
            )
            count = count + 1
        for a in m.out_args:
            if count > 0:
                self.out.write(",\n%*s" % (max_method_len + 2, ""))
            self.out.write(
                "OUT %s%*s %s"
                % (a.signature, max_signature_len - len(a.signature), "", a.name)
            )
            count = count + 1
        self.out.write(");\n")

    def print_signal_prototype(self, i, s, in_synopsis):
        max_signal_len = 0
        if in_synopsis:
            for _s in i.signals:
                max_signal_len = max(len(_s.name), max_signal_len)
        else:
            max_signal_len = max(len(s.name), max_signal_len)

        max_signature_len = 0
        if in_synopsis:
            for _s in i.signals:
                for a in _s.args:
                    max_signature_len = max(len(a.signature), max_signature_len)
        else:
            for a in s.args:
                max_signature_len = max(len(a.signature), max_signature_len)

        if in_synopsis:
            self.out.write(
                '<link linkend="gdbus-signal-%s.%s">%s</link>%*s ('
                % (
                    utils.dots_to_hyphens(i.name),
                    s.name,
                    s.name,
                    max_signal_len - len(s.name),
                    "",
                )
            )
        else:
            self.out.write("%s%*s (" % (s.name, max_signal_len - len(s.name), ""))
        count = 0
        for a in s.args:
            if count > 0:
                self.out.write(",\n%*s" % (max_signal_len + 2, ""))
            self.out.write(
                "%s%*s %s"
                % (a.signature, max_signature_len - len(a.signature), "", a.name)
            )
            count = count + 1
        self.out.write(");\n")

    def print_property_prototype(self, i, p, in_synopsis):
        max_property_len = 0
        if in_synopsis:
            for _p in i.properties:
                max_property_len = max(len(_p.name), max_property_len)
        else:
            max_property_len = max(len(p.name), max_property_len)

        max_signature_len = 0
        if in_synopsis:
            for _p in i.properties:
                max_signature_len = max(len(_p.signature), max_signature_len)
        else:
            max_signature_len = max(len(p.signature), max_signature_len)

        if in_synopsis:
            self.out.write(
                '<link linkend="gdbus-property-%s.%s">%s</link>%*s'
                % (
                    utils.dots_to_hyphens(i.name),
                    p.name,
                    p.name,
                    max_property_len - len(p.name),
                    "",
                )
            )
        else:
            self.out.write("%s%*s" % (p.name, max_property_len - len(p.name), ""))
        if p.readable and p.writable:
            access = "readwrite"
        elif p.readable:
            access = "readable "
        else:
            access = "writable "
        self.out.write("  %s  %s\n" % (access, p.signature))

    def print_synopsis_methods(self, i):
        self.out.write('  <refsynopsisdiv role="synopsis">\n')
        self.out.write('    <title role="synopsis.title">Methods</title>\n')
        self.out.write("    <synopsis>\n")
        for m in i.methods:
            self.print_method_prototype(i, m, in_synopsis=True)
        self.out.write("</synopsis>\n")
        self.out.write("  </refsynopsisdiv>\n")

    def print_synopsis_signals(self, i):
        self.out.write('  <refsect1 role="signal_proto">\n')
        self.out.write('    <title role="signal_proto.title">Signals</title>\n')
        self.out.write("    <synopsis>\n")
        for s in i.signals:
            self.print_signal_prototype(i, s, in_synopsis=True)
        self.out.write("</synopsis>\n")
        self.out.write("  </refsect1>\n")

    def print_synopsis_properties(self, i):
        self.out.write('  <refsect1 role="properties">\n')
        self.out.write('    <title role="properties.title">Properties</title>\n')
        self.out.write("    <synopsis>\n")
        for p in i.properties:
            self.print_property_prototype(i, p, in_synopsis=True)
        self.out.write("</synopsis>\n")
        self.out.write("  </refsect1>\n")

    def print_method(self, i, m):
        self.out.write(
            '<refsect2 role="method" id="gdbus-method-%s.%s">\n'
            % (utils.dots_to_hyphens(i.name), m.name)
        )
        self.out.write("  <title>The %s() method</title>\n" % (m.name))
        self.out.write(
            '  <indexterm zone="gdbus-method-%s.%s"><primary sortas="%s.%s">%s.%s()</primary></indexterm>\n'
            % (
                utils.dots_to_hyphens(i.name),
                m.name,
                i.name_without_prefix,
                m.name,
                i.name,
                m.name,
            )
        )
        self.out.write("<programlisting>\n")
        self.print_method_prototype(i, m, in_synopsis=False)
        self.out.write("</programlisting>\n")
        self.out.write("%s\n" % (self.expand_paras(m.doc_string, True)))
        if m.in_args or m.out_args:
            self.out.write('<variablelist role="params">\n')
            for a in m.in_args:
                self.out.write("<varlistentry>\n")
                self.out.write(
                    "  <term><literal>IN %s <parameter>%s</parameter></literal>:</term>\n"
                    % (a.signature, a.name)
                )
                self.out.write(
                    "  <listitem>%s</listitem>\n"
                    % (self.expand_paras(a.doc_string, True))
                )
                self.out.write("</varlistentry>\n")
            for a in m.out_args:
                self.out.write("<varlistentry>\n")
                self.out.write(
                    "  <term><literal>OUT %s <parameter>%s</parameter></literal>:</term>\n"
                    % (a.signature, a.name)
                )
                self.out.write(
                    "  <listitem>%s</listitem>\n"
                    % (self.expand_paras(a.doc_string, True))
                )
                self.out.write("</varlistentry>\n")
            self.out.write("</variablelist>\n")
        if len(m.since) > 0:
            self.out.write('<para role="since">Since %s</para>\n' % (m.since))
        if m.deprecated:
            self.out.write(
                "<warning><para>The %s() method is deprecated.</para></warning>"
                % (m.name)
            )
        self.out.write("</refsect2>\n")

    def print_signal(self, i, s):
        self.out.write(
            '<refsect2 role="signal" id="gdbus-signal-%s.%s">\n'
            % (utils.dots_to_hyphens(i.name), s.name)
        )
        self.out.write('  <title>The "%s" signal</title>\n' % (s.name))
        self.out.write(
            '  <indexterm zone="gdbus-signal-%s.%s"><primary sortas="%s::%s">%s::%s</primary></indexterm>\n'
            % (
                utils.dots_to_hyphens(i.name),
                s.name,
                i.name_without_prefix,
                s.name,
                i.name,
                s.name,
            )
        )
        self.out.write("<programlisting>\n")
        self.print_signal_prototype(i, s, in_synopsis=False)
        self.out.write("</programlisting>\n")
        self.out.write("%s\n" % (self.expand_paras(s.doc_string, True)))
        if s.args:
            self.out.write('<variablelist role="params">\n')
            for a in s.args:
                self.out.write("<varlistentry>\n")
                self.out.write(
                    "  <term><literal>%s <parameter>%s</parameter></literal>:</term>\n"
                    % (a.signature, a.name)
                )
                self.out.write(
                    "  <listitem>%s</listitem>\n"
                    % (self.expand_paras(a.doc_string, True))
                )
                self.out.write("</varlistentry>\n")
            self.out.write("</variablelist>\n")
        if len(s.since) > 0:
            self.out.write('<para role="since">Since %s</para>\n' % (s.since))
        if s.deprecated:
            self.out.write(
                '<warning><para>The "%s" signal is deprecated.</para></warning>'
                % (s.name)
            )
        self.out.write("</refsect2>\n")

    def print_property(self, i, p):
        self.out.write(
            '<refsect2 role="property" id="gdbus-property-%s.%s">\n'
            % (utils.dots_to_hyphens(i.name), p.name)
        )
        self.out.write('  <title>The "%s" property</title>\n' % (p.name))
        self.out.write(
            '  <indexterm zone="gdbus-property-%s.%s"><primary sortas="%s:%s">%s:%s</primary></indexterm>\n'
            % (
                utils.dots_to_hyphens(i.name),
                p.name,
                i.name_without_prefix,
                p.name,
                i.name,
                p.name,
            )
        )
        self.out.write("<programlisting>\n")
        self.print_property_prototype(i, p, in_synopsis=False)
        self.out.write("</programlisting>\n")
        self.out.write("%s\n" % (self.expand_paras(p.doc_string, True)))
        if len(p.since) > 0:
            self.out.write('<para role="since">Since %s</para>\n' % (p.since))
        if p.deprecated:
            self.out.write(
                '<warning><para>The "%s" property is deprecated.</para></warning>'
                % (p.name)
            )
        self.out.write("</refsect2>\n")

    def expand(self, s, expandParamsAndConstants):
        for key in self.expand_member_dict_keys:
            s = s.replace(key, self.expand_member_dict[key])
        for key in self.expand_iface_dict_keys:
            s = s.replace(key, self.expand_iface_dict[key])
        if expandParamsAndConstants:
            # replace @foo with <parameter>foo</parameter>
            s = re.sub(
                "@[a-zA-Z0-9_]*",
                lambda m: "<parameter>" + m.group(0)[1:] + "</parameter>",
                s,
            )
            # replace e.g. %TRUE with <constant>TRUE</constant>
            s = re.sub(
                "%[a-zA-Z0-9_]*",
                lambda m: "<constant>" + m.group(0)[1:] + "</constant>",
                s,
            )
        return s

    def expand_paras(self, s, expandParamsAndConstants):
        s = textwrap.dedent(self.expand(s, expandParamsAndConstants)).rstrip()
        res = []
        if not s.startswith("<para>"):
            res.append("<para>")
        for line in s.split("\n"):
            line = line.rstrip()
            if not line:
                line = "</para><para>"
            res.append(line)
        if not s.endswith("</para>"):
            res.append("</para>")
        return "\n".join(res)

    def generate_expand_dicts(self):
        self.expand_member_dict = {}
        self.expand_iface_dict = {}
        for i in self.ifaces:
            key = "#%s" % (i.name)
            value = '<link linkend="gdbus-interface-%s.top_of_page">%s</link>' % (
                utils.dots_to_hyphens(i.name),
                i.name,
            )
            self.expand_iface_dict[key] = value
            for m in i.methods:
                key = "%s.%s()" % (i.name, m.name)
                value = '<link linkend="gdbus-method-%s.%s">%s()</link>' % (
                    utils.dots_to_hyphens(i.name),
                    m.name,
                    m.name,
                )
                self.expand_member_dict[key] = value
            for s in i.signals:
                key = "#%s::%s" % (i.name, s.name)
                value = '<link linkend="gdbus-signal-%s.%s">"%s"</link>' % (
                    utils.dots_to_hyphens(i.name),
                    s.name,
                    s.name,
                )
                self.expand_member_dict[key] = value
            for p in i.properties:
                key = "#%s:%s" % (i.name, p.name)
                value = '<link linkend="gdbus-property-%s.%s">"%s"</link>' % (
                    utils.dots_to_hyphens(i.name),
                    p.name,
                    p.name,
                )
                self.expand_member_dict[key] = value
        # Make sure to expand the keys in reverse order so e.g. #org.foo.Iface:MediaCompat
        # is evaluated before #org.foo.Iface:Media ...
        self.expand_member_dict_keys = sorted(
            self.expand_member_dict.keys(), reverse=True
        )
        self.expand_iface_dict_keys = sorted(
            self.expand_iface_dict.keys(), reverse=True
        )

    def generate(self, docbook, outdir):
        for i in self.ifaces:
            self.out = open(path.join(outdir, "%s-%s.xml" % (docbook, i.name)), "w")
            self.out.write("")
            self.out.write('<?xml version="1.0" encoding="utf-8"?>\n')
            self.out.write(
                '<!DOCTYPE refentry PUBLIC "-//OASIS//DTD DocBook XML V4.1.2//EN"\n'
            )
            self.out.write(
                '               "http://www.oasis-open.org/docbook/xml/4.1.2/docbookx.dtd" [\n'
            )
            self.out.write("]>\n")
            self.out.write('<refentry id="gdbus-%s">\n' % (i.name))
            self.out.write("  <refmeta>")
            self.out.write(
                '    <refentrytitle role="top_of_page" id="gdbus-interface-%s.top_of_page">%s</refentrytitle>\n'
                % (utils.dots_to_hyphens(i.name), i.name)
            )
            self.out.write(
                '  <indexterm zone="gdbus-interface-%s.top_of_page"><primary sortas="%s">%s</primary></indexterm>\n'
                % (utils.dots_to_hyphens(i.name), i.name_without_prefix, i.name)
            )
            self.out.write("  </refmeta>")

            self.out.write("  <refnamediv>")
            self.out.write("    <refname>%s</refname>" % (i.name))
            self.out.write("    <refpurpose>%s</refpurpose>" % (i.doc_string_brief))
            self.out.write("  </refnamediv>")

            if len(i.methods) > 0:
                self.print_synopsis_methods(i)
            if len(i.signals) > 0:
                self.print_synopsis_signals(i)
            if len(i.properties) > 0:
                self.print_synopsis_properties(i)

            self.out.write(
                '<refsect1 role="desc" id="gdbus-interface-%s">\n'
                % (utils.dots_to_hyphens(i.name))
            )
            self.out.write('  <title role="desc.title">Description</title>\n')
            self.out.write("  %s\n" % (self.expand_paras(i.doc_string, True)))
            if len(i.since) > 0:
                self.out.write('  <para role="since">Since %s</para>\n' % (i.since))
            if i.deprecated:
                self.out.write(
                    "<warning><para>The %s interface is deprecated.</para></warning>"
                    % (i.name)
                )
            self.out.write("</refsect1>\n")

            if len(i.methods) > 0:
                self.out.write(
                    '<refsect1 role="details" id="gdbus-methods-%s">\n' % (i.name)
                )
                self.out.write('  <title role="details.title">Method Details</title>\n')
                for m in i.methods:
                    self.print_method(i, m)
                self.out.write("</refsect1>\n")

            if len(i.signals) > 0:
                self.out.write(
                    '<refsect1 role="details" id="gdbus-signals-%s">\n' % (i.name)
                )
                self.out.write('  <title role="details.title">Signal Details</title>\n')
                for s in i.signals:
                    self.print_signal(i, s)
                self.out.write("</refsect1>\n")

            if len(i.properties) > 0:
                self.out.write(
                    '<refsect1 role="details" id="gdbus-properties-%s">\n' % (i.name)
                )
                self.out.write(
                    '  <title role="details.title">Property Details</title>\n'
                )
                for s in i.properties:
                    self.print_property(i, s)
                self.out.write("</refsect1>\n")

            self.out.write("</refentry>\n")
            self.out.write("\n")
