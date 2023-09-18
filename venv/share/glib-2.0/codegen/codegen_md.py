# SPDX-FileCopyrightText: 2023 Guido Günther
# base on # codegen_rst.py (C) 2022 Emmanuele Bassi
#
# SPDX-License-Identifier: LGPL-2.1-or-later

import os
import re

from . import utils
import textwrap

# Disable line length warnings as wrapping the templates would be hard
# flake8: noqa: E501


class MdCodeGenerator:
    """Generates documentation in Markdown format."""

    def __init__(self, ifaces):
        self.ifaces = ifaces
        self._generate_expand_dicts()

    def _expand(self, s, expandParamsAndConstants):
        """Expands parameters and constant literals."""
        res = []
        for line in textwrap.dedent(s).split("\n"):
            line = line.rstrip()
            if line == "":
                res.append("")
                continue
            for key in self._expand_member_dict_keys:
                line = line.replace(key, self._expand_member_dict[key])
            for key in self._expand_iface_dict_keys:
                line = line.replace(key, self._expand_iface_dict[key])
            if expandParamsAndConstants:
                # replace @foo with `foo`
                line = re.sub(
                    "@[a-zA-Z0-9_]*",
                    lambda m: "`" + m.group(0)[1:] + "`",
                    line,
                )
                # replace e.g. %TRUE with ``TRUE``
                line = re.sub(
                    "%[a-zA-Z0-9_]*",
                    lambda m: "`" + m.group(0)[1:] + "`",
                    line,
                )
            res.append(line)
        return "\n".join(res)

    def _generate_expand_dicts(self):
        """Generates the dictionaries used to expand gtk-doc sigils."""
        self._expand_member_dict = {}
        self._expand_iface_dict = {}
        for i in self.ifaces:
            key = f"#{i.name}"
            value = f"`{i.name}`_"
            self._expand_iface_dict[key] = value

            for m in i.methods:
                key = "%s.%s()" % (i.name, m.name)
                value = f"`{i.name}.{m.name}`_"
                self._expand_member_dict[key] = value

            for s in i.signals:
                key = "#%s::%s" % (i.name, s.name)
                value = f"`{i.name}::{s.name}`_"
                self._expand_member_dict[key] = value

            for p in i.properties:
                key = "#%s:%s" % (i.name, p.name)
                value = f"`{i.name}:{p.name}`_"
                self._expand_member_dict[key] = value

        # Make sure to expand the keys in reverse order so e.g. #org.foo.Iface:MediaCompat
        # is evaluated before #org.foo.Iface:Media ...
        self._expand_member_dict_keys = sorted(
            self._expand_member_dict.keys(), reverse=True
        )
        self._expand_iface_dict_keys = sorted(
            self._expand_iface_dict.keys(), reverse=True
        )

    def _generate_header(self, iface):
        """Generates the header and preamble of the document."""
        header_len = len(iface.name)
        res = [
            f"Title: {iface.name} D-Bus Interface",
            f"Slug: {iface.name}",
            "",
            "# " + iface.name,
            "",
            "## Description",
            "",
            iface.doc_string_brief.strip(),
            "",
            self._expand(iface.doc_string, True),
            "",
        ]
        if iface.since:
            res += [
                f"Interface available since: {iface.since}.",
                "",
            ]
        if iface.deprecated:
            res += [
                "*Warning*: This interface is deprecated.",
                "",
            ]
        res += [""]
        return "\n".join(res)

    def _generate_section(self, title, name):
        """Generates a section with the given title."""
        res = [
            "### " + title,
            "",
        ]
        return "\n".join(res)

    def _generate_properties(self, iface):
        """Generates the properties section."""
        res = []
        for p in iface.properties:
            title = f"{iface.name}:{p.name}"
            if p.readable and p.writable:
                access = "readwrite"
            elif p.writable:
                access = "writable"
            else:
                access = "readable"
            res += [
                "### " + title,
                "",
                "```",
                f"    {p.name} {access} {p.signature}",
                "```",
                "",
                self._expand(p.doc_string, True),
                "",
            ]
            if p.since:
                res += [
                    f"Property available since: {p.since}.",
                    "",
                ]
            if p.deprecated:
                res += [
                    "*Warning*: This property is deprecated.",
                    "",
                ]
            res += [""]
        return "\n".join(res)

    def _generate_method_signature(self, method):
        """Generates the method signature as a code block."""
        res = [
            "```",
        ]
        n_in_args = len(method.in_args)
        n_out_args = len(method.out_args)
        if n_in_args == 0 and n_out_args == 0:
            res += [
                f"    {method.name} ()",
            ]
        else:
            res += [
                f"    {method.name} (",
            ]
            for idx, arg in enumerate(method.in_args):
                if idx == n_in_args - 1 and n_out_args == 0:
                    res += [
                        f"      IN {arg.name} {arg.signature}",
                    ]
                else:
                    res += [
                        f"      IN {arg.name} {arg.signature},",
                    ]
            for idx, arg in enumerate(method.out_args):
                if idx == n_out_args - 1:
                    res += [
                        f"      OUT {arg.name} {arg.signature}",
                    ]
                else:
                    res += [
                        f"      OUT {arg.name} {arg.signature},",
                    ]
            res += [
                "    )",
            ]
        res += ["```"]
        return "\n".join(res)

    def _generate_methods(self, iface):
        """Generates the methods section."""
        res = []
        for m in iface.methods:
            title = f"{iface.name}.{m.name}"
            res += [
                "### " + title,
                "",
                self._generate_method_signature(m),
                "",
                self._expand(m.doc_string, True),
                "",
            ]
            for a in m.in_args:
                arg_desc = self._expand(a.doc_string, True)
                res += [
                    f"* {a.name}: {arg_desc}",
                    "",
                ]
            res += [""]
            if m.since:
                res += [
                    f"Method available since: {m.since}.",
                    "",
                ]
            if m.deprecated:
                res += [
                    "*Warning*: This method is deprecated.",
                    "",
                ]
            res += [""]
        return "\n".join(res)

    def _generate_signal_signature(self, signal):
        """Generates the signal signature."""
        res = [
            "```",
        ]
        n_args = len(signal.args)
        if n_args == 0:
            res += [
                f"    {signal.name} ()",
            ]
        else:
            res += [
                f"    {signal.name} (",
            ]
            for idx, arg in enumerate(signal.args):
                if idx == n_args - 1:
                    res += [
                        f"      {arg.name} {arg.signature}",
                    ]
                else:
                    res += [
                        f"      {arg.name} {arg.signature},",
                    ]
            res += [
                "    )",
            ]
        res += ["```"]
        return "\n".join(res)

    def _generate_signals(self, iface):
        """Generates the signals section."""
        res = []
        for s in iface.signals:
            title = f"{iface.name}::{s.name}"
            res += [
                "### " + title,
                "",
                self._generate_signal_signature(s),
                "",
                self._expand(s.doc_string, True),
                "",
            ]
            for a in s.args:
                arg_desc = self._expand(a.doc_string, True)
                res += [
                    f"{a.name}",
                    f"  {arg_desc}",
                    "",
                ]
            res += [""]
            if s.since:
                res += [
                    f"Signal available since: {s.since}.",
                    "",
                ]
            if s.deprecated:
                res += [
                    "*Warning*: This signal is deprecated.",
                    "",
                ]
            res += [""]
        return "\n".join(res)

    def generate(self, md, outdir):
        """Generates the Markdown file for each interface."""
        for i in self.ifaces:
            with open(os.path.join(outdir, f"{md}-{i.name}.md"), "w") as outfile:
                outfile.write(self._generate_header(i))
                if len(i.properties) > 0:
                    outfile.write(self._generate_section("Properties", i.name))
                    outfile.write(self._generate_properties(i))
                if len(i.methods) > 0:
                    outfile.write(self._generate_section("Methods", i.name))
                    outfile.write(self._generate_methods(i))
                if len(i.signals) > 0:
                    outfile.write(self._generate_section("Signals", i.name))
                    outfile.write(self._generate_signals(i))
