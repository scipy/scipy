# Copyright (c) 2023, Riverbank Computing Limited
# All rights reserved.
#
# This copy of SIP is licensed for use under the terms of the SIP License
# Agreement.  See the file LICENSE for more details.
#
# This copy of SIP may also used under the terms of the GNU General Public
# License v2 or v3 as published by the Free Software Foundation which can be
# found in the files LICENSE-GPL2 and LICENSE-GPL3 included in this package.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ('AS IS'
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


from ...specification import ValueType

from .argument import ArgumentFormatter
from .base_formatter import BaseFormatter
from .enum import EnumFormatter
from .variable import VariableFormatter


class ValueListFormatter(BaseFormatter):
    """ This creates various string representations of a list of values. """

    @property
    def cpp_expression(self):
        """ The C++ representation of the value list as an expression. """

        return self._expression()

    def py_expression(self, embedded=False, as_xml=False):
        """ The Python representation of the value list as an expression. """

        return self._expression(as_python=True, embedded=embedded,
                as_xml=as_xml)

    def as_rest_ref(self, as_xml=False):
        """ Return the Python representation of the value list as a reST
        reference.
        """

        value_list = self.object

        # The value must be a scoped name and we don't handle expressions.
        if len(value_list) != 1 or value_list[0].value_type is not ValueType.SCOPED:
            return None

        target = value_list[0].value

        # See if it is an attribute.
        for variable in self.spec.variables:
            if variable.fq_cpp_name == target:
                return VariableFormatter(self.spec, variable).as_rest_ref()

        # See if it is an enum member.
        target_scope = target.scope
        if target_scope is not None:
            target_scope.make_absolute()

        target_base_name = target.base_name

        for enum in self.spec.enums:
            # Look for the member name first before working out if it is the
            # correct enum.
            for member in enum.members:
                if member.cpp_name == target_base_name:
                    formatter = EnumFormatter(self.spec, enum)

                    if enum.is_scoped:
                        # It's a scoped enum so the fully qualified name of the
                        # enum must match the scope of the name.
                        if target_scope is not None and enum.fq_cpp_name == target_scope:
                            return formatter.member_as_rest_ref(member)
                    else:
                        # It's a traditional enum so the scope of the enum must
                        # match the scope of the name.
                        if (enum.scope is None and target_scope is None) or (enum.scope is not None and target_scope is not None and enum.scope.iface_file.fq_cpp_name == target_scope):
                            return formatter.member_as_rest_ref(member)

                    break

        return None

    def _expression(self, as_python=False, embedded=False, as_xml=False):
        """ The representation of the value list as an expression. """

        s = ''

        for value in self.object:
            if value.cast is not None and not as_python:
                s += '(' + value.cast.as_cpp + ')'

            if value.unary_operator is not None:
                s += value.unary_operator

            if value.value_type is ValueType.QCHAR:
                if value.value == '"' and embedded:
                    s += "'\\\"'"
                elif value.value.isprintable():
                    s += "'" + value.value + "'"
                else:
                    s += f"'\\{ord(value.value):03o}'"

            elif value.value_type is ValueType.STRING:
                quote = "\\\"" if embedded else "\""

                s += quote

                for ch in value.value:
                    escape = True

                    if ch in "\\\"":
                        pass
                    elif ch == '\n':
                        ch = 'n';
                    elif ch == '\r':
                        ch = 'r'
                    elif ch == '\t':
                        ch = 't'
                    else:
                        escape = False

                    if escape:
                        s += "\\"

                    s += ch

                s += quote

            elif value.value_type in (ValueType.NUMERIC, ValueType.REAL):
                s += str(value.value)

            elif value.value_type is ValueType.SCOPED:
                if as_python:
                    s += value.value.as_py
                else:
                    s += value.value.as_cpp

            elif value.value_type is ValueType.FCALL:
                arg_formatter = ArgumentFormatter(self.spec,
                        value.value.result)

                if as_python:
                    s += arg_formatter.as_py_type()
                else:
                    s += arg_formatter.cpp_type()

                args = [ValueListFormatter(self.spec, a)._expression(
                                as_python=as_python, embedded=embedded,
                                as_xml=as_xml)
                        for a in value.value.args]

                separator = ',' if as_xml else ', '
                s += '(' + separator.join(args) + ')'

            elif value.value_type is ValueType.EMPTY:
                s += '{}'

            if value.binary_operator is not None:
                s += value.binary_operator

        return s
