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


from ...scoped_name import STRIP_NONE
from ...specification import ArgumentType, ArrayArgument, ClassKey, ValueType

from ..type_hints import format_voidptr, TypeHintManager

from .base_formatter import BaseFormatter
from .utils import format_scoped_py_name


class ArgumentFormatter(BaseFormatter):
    """ This creates various string representations of an argument. """

    def cpp_type(self, *, name=None, scope=None, strip=STRIP_NONE,
            make_public=False, use_typename=True, as_xml=False):
        """ Return the argument as a C++ type. """

        arg = self.object

        original_typedef = arg.original_typedef
        nr_derefs = len(arg.derefs)
        is_reference = arg.is_reference

        s = ''

        if use_typename and original_typedef is not None and not original_typedef.no_type_name and arg.array is not ArrayArgument.ARRAY_SIZE:
            # Reverse the previous expansion of the typedef.
            if arg.is_const and not original_typedef.type.is_const:
                s += 'const '

            nr_derefs -= len(original_typedef.type.derefs)

            if original_typedef.type.is_reference:
                is_reference = False

            s += original_typedef.fq_cpp_name.cpp_stripped(strip)

        else:
            # A function type is handled differently because of the position of
            # the name.
            if arg.type is ArgumentType.FUNCTION:
                s += ArgumentFormatter(self.spec,
                        arg.definition.result).cpp_type(scope=scope,
                                strip=strip, as_xml=as_xml)

                s += ' (' + '*' * nr_derefs + name + ')('

                s += SignatureFormatter(self.spec,
                        arg.definition).cpp_arguments(scope=scope,
                                as_xml=as_xml)

                s += ')'

                return s

            if arg.is_const:
                s += 'const '

            if arg.type in (ArgumentType.SBYTE, ArgumentType.SSTRING):
                s += 'signed char'

            elif arg.type in (ArgumentType.UBYTE, ArgumentType.USTRING):
                s += 'unsigned char'

            elif arg.type is ArgumentType.WSTRING:
                s += 'wchar_t'

            elif arg.type in (ArgumentType.BYTE, ArgumentType.ASCII_STRING, ArgumentType.LATIN1_STRING, ArgumentType.UTF8_STRING, ArgumentType.STRING):
                s += 'char'

            elif arg.type is ArgumentType.USHORT:
                s += 'unsigned short'

            elif arg.type is ArgumentType.SHORT:
                s += 'short'

            elif arg.type is ArgumentType.UINT:
                # Qt4 moc uses "uint" in signal signatures.  We do all the time
                # and hope it is always defined.
                s += 'uint'

            elif arg.type in (ArgumentType.INT, ArgumentType.CINT):
                s += 'int'

            elif arg.type is ArgumentType.HASH:
                s += 'Py_hash_t'

            elif arg.type is ArgumentType.SSIZE:
                s += 'Py_ssize_t'

            elif arg.type is ArgumentType.SIZE:
                s += 'size_t'

            elif arg.type is ArgumentType.ULONG:
                s += 'unsigned long'

            elif arg.type is ArgumentType.LONG:
                s += 'long'

            elif arg.type is ArgumentType.ULONGLONG:
                s += 'unsigned long long'

            elif arg.type is ArgumentType.LONGLONG:
                s += 'long long'

            elif arg.type is ArgumentType.STRUCT:
                s += 'struct ' + arg.definition.as_cpp

            elif arg.type is ArgumentType.UNION:
                s += 'union ' + arg.definition.as_cpp

            elif arg.type is ArgumentType.CAPSULE:
                nr_derefs = 1
                s += 'void'

            elif arg.type in (ArgumentType.FAKE_VOID, ArgumentType.VOID):
                s += 'void'

            elif arg.type in (ArgumentType.BOOL, ArgumentType.CBOOL):
                s += 'bool'

            elif arg.type in (ArgumentType.FLOAT, ArgumentType.CFLOAT):
                s += 'float'

            elif arg.type in (ArgumentType.DOUBLE, ArgumentType.CDOUBLE):
                s += 'double'

            elif arg.type is ArgumentType.DEFINED:
                # The only defined types still remaining are arguments to
                # templates and default values.
                if as_xml:
                    s += arg.definition.as_py
                else:
                    if self.spec.c_bindings:
                        s += 'struct '

                    s += arg.definition.cpp_stripped(strip)

            elif arg.type is ArgumentType.MAPPED:
                s += ArgumentFormatter(self.spec,
                        arg.definition.type).cpp_type(scope=scope, strip=strip,
                                as_xml=as_xml)

            elif arg.type is ArgumentType.CLASS:
                from .klass import ClassFormatter

                if self.spec.c_bindings:
                    s += 'union ' if arg.definition.class_key is ClassKey.UNION else 'struct '

                s += ClassFormatter(self.spec, arg.definition).scoped_name(
                        scope=scope, strip=strip, make_public=make_public,
                        as_xml=as_xml)

            elif arg.type is ArgumentType.TEMPLATE:
                from .template import TemplateFormatter

                s += TemplateFormatter(self.spec, arg.definition, scope).cpp_type(
                        strip=strip, as_xml=as_xml)

            elif arg.type is ArgumentType.ENUM:
                enum = arg.definition

                if enum.fq_cpp_name is None or (enum.is_protected and not make_public):
                    s += 'int'
                else:
                    s += enum.fq_cpp_name.cpp_stripped(strip)

            elif arg.type in (ArgumentType.PYOBJECT, ArgumentType.PYTUPLE, ArgumentType.PYLIST, ArgumentType.PYDICT, ArgumentType.PYCALLABLE, ArgumentType.PYSLICE, ArgumentType.PYTYPE, ArgumentType.PYBUFFER, ArgumentType.PYENUM, ArgumentType.ELLIPSIS):
                s += 'PyObject *'

        space_before_name = True

        for i in range(nr_derefs):
            # Note that we don't put a space before the '*' so that Qt
            # normalised signal signatures are correct.
            s += '*'
            space_before_name = False

            if arg.derefs[i]:
                s += ' const'
                space_before_name = True

        if is_reference:
            s += '&'

        if name:
            if space_before_name:
                s += ' '

            s += name

        return s

    # The types that are implicitly pointers.
    _IMPLICIT_POINTERS = (ArgumentType.PYOBJECT, ArgumentType.PYTUPLE,
        ArgumentType.PYLIST, ArgumentType.PYDICT, ArgumentType.PYCALLABLE,
        ArgumentType.PYSLICE, ArgumentType.PYTYPE, ArgumentType.CAPSULE,
        ArgumentType.PYBUFFER, ArgumentType.PYENUM)

    def py_default_value(self, type_name, embedded=False, as_xml=False):
        """ Return the Python representation of the argument's default value.
        """

        from .value_list import ValueListFormatter

        arg = self.object

        # Use any explicitly provided documentation.
        if arg.type_hints is not None and arg.type_hints.default_value is not None:
            return arg.type_hints.default_value

        # Translate some special cases.
        if len(arg.default_value) == 1 and arg.default_value[0].value_type is ValueType.NUMERIC:
            value = arg.default_value[0].value

            if value == 0 and ('voidptr' in type_name or len(arg.derefs) > 0 or arg.type in self._IMPLICIT_POINTERS):
                return 'None'

            if arg.type in (ArgumentType.BOOL, ArgumentType.CBOOL):
                return 'True' if value else 'False'

        return ValueListFormatter(self.spec, arg.default_value).py_expression(
                embedded=embedded, as_xml=as_xml)

    def as_py_type(self, pep484=False, default_value=False, as_xml=False):
        """ Return the argument as a Python type. """

        arg = self.object

        scope, name = self._py_arg(pep484, as_xml)

        s = format_scoped_py_name(scope, name)

        if default_value and arg.default_value is not None:
            if arg.name is not None:
                s += ' ' + arg.name.name

            s += '=' + self.py_default_value(name)

        return s

    def as_rest_ref(self, out, as_xml=False):
        """ Return the argument as a reST reference. """

        arg = self.object

        s = ''

        hint = self._get_hint(out)

        if hint is None:
            if arg.type is ArgumentType.CLASS:
                from .klass import ClassFormatter

                s += ClassFormatter(self.spec, arg.definition).as_rest_ref()
            elif arg.type is ArgumentType.ENUM:
                if arg.definition.py_name is not None:
                    from .enum import EnumFormatter

                    s += EnumFormatter(self.spec, arg.definition).as_rest_ref()
                else:
                    s += 'int'
            elif arg.type is ArgumentType.MAPPED:
                # There would normally be a type hint.
                s += "unknown-type"
            else:
                s += self.as_py_type(as_xml=as_xml)
        else:
            s += TypeHintManager(self.spec).as_rest_ref(hint, out,
                    as_xml=as_xml)

        return s

    def as_type_hint(self, module, out, defined):
        """ Return the argument as a type hint. """

        arg = self.object

        s = ''

        hint = self._get_hint(out)

        if hint is None:
            if arg.type is ArgumentType.CLASS:
                from .klass import ClassFormatter

                s += ClassFormatter(self.spec, arg.definition).as_type_hint(
                        module, defined)
            elif arg.type is ArgumentType.ENUM:
                if arg.definition.py_name is not None:
                    from .enum import EnumFormatter

                    s += EnumFormatter(self.spec, arg.definition).as_type_hint(
                            module, defined)
                else:
                    s += 'int'
            elif arg.type is ArgumentType.MAPPED:
                # There would normally be a type hint.
                s += 'typing.Any'
            else:
                s += self.as_py_type(pep484=True)
        else:
            s += TypeHintManager(self.spec).as_type_hint(hint, module, out,
                    defined)

        return s

    def _get_hint(self, out):
        """ Return the raw type hint. """

        arg = self.object

        # Use any explicit type hint unless the argument is constrained.
        if arg.type_hints is None:
            hint = None
        elif out:
            hint = arg.type_hints.hint_out
        elif arg.is_constrained:
            hint = None
        else:
            hint = arg.type_hints.hint_in

        return hint

    def _py_arg(self, pep484, as_xml):
        """ Return an argument as a 2-tuple of scope and name. """

        type = self.object.type
        definition = self.object.definition

        sip_module_name = self.spec.sip_module + '.' if self.spec.sip_module else ''

        scope = None
        name = "unknown-type"

        if type is ArgumentType.CLASS:
            scope = definition.scope
            name = definition.py_name.name

        elif type is ArgumentType.MAPPED:
            if definition.py_name is not None:
                name = definition.py_name.name

        elif type is ArgumentType.ENUM:
            if definition.py_name is None:
                name = 'int'
            else:
                scope = definition.scope
                name = definition.py_name.name

        elif type is ArgumentType.DEFINED:
            name = definition.as_py

        elif type is ArgumentType.CAPSULE:
            name = definition.base_name

        elif type in (ArgumentType.STRUCT, ArgumentType.UNION, ArgumentType.VOID):
            name = format_voidptr(self.spec, pep484, as_xml)

        elif type in (ArgumentType.STRING, ArgumentType.SSTRING, ArgumentType.USTRING):
            name = 'bytes'

        elif type in (ArgumentType.WSTRING, ArgumentType.ASCII_STRING, ArgumentType.LATIN1_STRING, ArgumentType.UTF8_STRING):
            name = 'bytes' if self.object.array is ArrayArgument.ARRAY else 'str'

        elif type in (ArgumentType.BYTE, ArgumentType.SBYTE, ArgumentType.UBYTE, ArgumentType.USHORT, ArgumentType.UINT, ArgumentType.LONG, ArgumentType.LONGLONG, ArgumentType.ULONG, ArgumentType.ULONGLONG, ArgumentType.SHORT, ArgumentType.INT, ArgumentType.CINT, ArgumentType.SSIZE, ArgumentType.SIZE, ArgumentType.HASH):
            name = 'int'

        elif type in (ArgumentType.FLOAT, ArgumentType.CFLOAT, ArgumentType.DOUBLE, ArgumentType.CDOUBLE):
            name = 'float'

        elif type in (ArgumentType.BOOL, ArgumentType.CBOOL):
            name = 'bool'

        elif type in (ArgumentType.PYOBJECT, ArgumentType.ELLIPSIS):
            name = 'typing.Any' if pep484 else 'Any'

        elif type is ArgumentType.PYTUPLE:
            name = 'typing.Tuple' if pep484 else 'Tuple'

        elif type is ArgumentType.PYLIST:
            name = 'typing.List' if pep484 else 'List'

        elif type is ArgumentType.PYDICT:
            name = 'typing.Dict' if pep484 else 'Dict'

        elif type is ArgumentType.PYCALLABLE:
            name = 'typing.Callable[..., None]' if pep484 else 'Callable[..., None]'

        elif type is ArgumentType.PYSLICE:
            name = 'slice'

        elif type is ArgumentType.PYTYPE:
            name = 'type'

        elif type is ArgumentType.PYBUFFER:
            if pep484:
                name = sip_module_name + 'Buffer'
            else:
                # This replicates sip.pyi.
                name = f'Union[bytes, bytearray, memoryview, {sip_module_name}array, {sip_module_name}voidptr]'

        elif type is ArgumentType.PYENUM:
            name = 'enum.Enum'

        return scope, name
