# Copyright (c) 2022, Riverbank Computing Limited
# All rights reserved.
#
# This copy of SIP is licensed for use under the terms of the SIP License
# Agreement.  See the file LICENSE for more details.
#
# This copy of SIP may also used under the terms of the GNU General Public
# License v2 or v3 as published by the Free Software Foundation which can be
# found in the files LICENSE-GPL2 and LICENSE-GPL3 included in this package.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
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


from copy import copy

from .scoped_name import ScopedName
from .specification import ArgumentType, IfaceFileType
from .utils import append_iface_file, argument_as_str, same_base_type


def encoded_template_name(template):
    """ Return the encoded name of a template. """

    snd = ScopedName(template.cpp_name)

    for ad in template.types.args:
        flags = 0

        if ad.is_const:
            flags |= 1

        if ad.is_reference:
            flags |= 2

        # We use numbers so they don't conflict with names.
        encoding = '{:02d}{}{}'.format(ad.type.value, flags, len(ad.derefs))

        if ad.type is ArgumentType.DEFINED:
            arg_snd = ScopedName(ad.definition)
        elif ad.type is ArgumentType.TEMPLATE:
            arg_snd = encoded_template_name(ad.definition)
        elif ad.type is ArgumentType.STRUCT:
            arg_snd = ScopedName(ad.definition)
        else:
            arg_snd = None

        if arg_snd is None:
            snd.append(encoding)
        else:
            # Replace the first element of the argument name with a copy with
            # the encoding prepended.
            arg_snd[0] = encoding + arg_snd[0]

            for name in arg_snd:
                snd.append(name)

    return snd


def same_template_signature(sig1, sig2, deep=False):
    """ Return True if the template signatures are the same.  A deep comparison
    is used for mapped type templates where we want to recurse into any nested
    templates.
    """

    if len(sig1.args) != len(sig2.args):
        return False

    for type1, type2 in zip(sig1.args, sig2.args):
        # If we are doing a shallow comparision (ie. for class templates) then
        # a type name in the first signature matches anything in the second
        # signature.
        if type1.type is ArgumentType.DEFINED and not deep:
            continue

        # For type names only compare the references and pointers, and do the
        # same for any nested templates.
        if type1.type is ArgumentType.DEFINED and type2.type is ArgumentType.DEFINED:
            if type1.is_reference != type2.is_reference or len(type1.derefs) != len(type2.derefs):
                return False

        elif type1.type is ArgumentType.TEMPLATE and type2.type is ArgumentType.TEMPLATE:
            if not same_template_signature(type1.definition.types, type2.definition.types, deep=deep):
                return False

        elif not same_base_type(type1, type2):
            return False

    return True


def template_code(spec, used, proto_code, expansions):
    """ Return a copy of an optional CodeBlock object with sub-strings replaced
    by corresponding values.
    """

    # Handle the trivial case.
    if proto_code is None:
        return None

    return _template_code_block(spec, used, proto_code, expansions)


def template_code_blocks(spec, used, proto_code_blocks, expansions):
    """ Return a copy of a list of CodeBlock objects with sub-strings replaced
    by corresponding values.
    """

    return [_template_code_block(spec, used, pc, expansions)
            for pc in proto_code_blocks]


def template_expansions(template_names, instantiation_values,
        declared_names=None):
    """ Return a dict of expansions to be applied when instantiating mapped
    type of class templates (including handwritten code).  The key is the
    symbolic name of a template argument and the value is the replacement to be
    used in a particular instantiation.
    """

    # TODO: This is broken (or possibly just over-complicated).  The
    # declaration of template parameters is used for two purposes: firstly to
    # allow real types to be used to distinguish between overloaded templates;
    # secondly to define the names of the parameters that will be replaced by
    # values provided at instantiation.  It's possible that only the latter
    # will ever be expanded.

    expansions = {}

    for arg_nr, name_arg in enumerate(template_names.args):
        if name_arg.type is ArgumentType.DEFINED:
            # If the type names have been declared (as they are with a mapped
            # type template) check that this is one of them.
            if declared_names is not None:
                # Only consider unscoped names.
                if not name_arg.definition.is_simple:
                    continue

                for declared in declared_names.args:
                    # Skip anything but simple names.
                    if declared.type is not ArgumentType.DEFINED:
                        continue

                    if name_arg.definition.base_name == declared.definition.base_name:
                        name = name_arg.definition.base_name
                        break
                else:
                    continue
            else:
                name = name_arg.definition.base_name

            # Get the corresponding value.  For defined types we don't want any
            # indirection or references.
            value_arg = instantiation_values.args[arg_nr]

            if value_arg.type is ArgumentType.DEFINED:
                value = str(value_arg.definition)
            else:
                value = argument_as_str(value_arg)

            # We do want const.
            if value_arg.is_const:
                value = 'const ' + value;

            expansions[name] = value

        elif name_arg.type is ArgumentType.TEMPLATE:
            value_arg = instantiation_values.args[arg_nr]

            # These checks shouldn't be necessary, but...
            if value_arg.type is ArgumentType.TEMPLATE and len(name_arg.definition.types.args) == len(value_arg.definition.types.args):
                expansions.update(
                        template_expansions(name_arg.definition.types,
                                value_arg.definition.types, declared_names))

    return expansions


def template_string(proto_str, expansions, scope_replacement=None):
    """ Return a copy of a string with sub-strings replaced by corresponding
    values.
    """

    for name, value in expansions.items():
        value = _strip_const(value)

        # Translate any C++ scoping.
        if scope_replacement is not None:
            value = value.replace('::', scope_replacement)

        # Perform any replacement.
        proto_str = proto_str.replace(name, value)

    return proto_str


def _strip_const(s):
    """ Strip any leading 'const' from a string. """

    if s.startswith('const '):
        s = s[6:]

    return s


def _template_code_block(spec, used, proto_code, expansions):
    """ Return a copy of a CodeBlock object with sub-strings replaced by
    corresponding values or the original if there were no substitutions.
    """

    i_code = copy(proto_code)

    i_lines = []
    for proto_line in proto_code.text.split('\n'):
        i_line = proto_line

        # Don't do any substitution in lines that appear to be preprocessor
        # directives.  This prevents #include'd file names being broken.
        if not proto_line.lstrip().startswith('#'):
            # Go through each expansion.
            for name, value in expansions.items():
                # Look for the name at the current position in the current
                # line.
                pos = i_line.find(name)
                while pos >= 0:
                    # See if the name is referring to a generated type
                    # structure.
                    for gen_type in ('sipType_', 'sipException_'):
                        if i_line[:pos].endswith(gen_type):
                            value = _strip_const(value)

                            _add_used_from_code(spec, used, value)

                            # Convert the value to the rest of the name of the
                            # generated type structure.
                            if value.startswith('::'):
                                value = value[2:]

                            value = value.replace('::', '_')

                            break

                    # Perform the substitution and update the current position.
                    i_line = i_line[0:pos] + value + i_line[pos + len(name):]
                    pos = i_line.find(name, pos + len(value))

        i_lines.append(i_line)

    i_code.text = '\n'.join(i_lines)

    # Return the prototype itself if nothing changed.
    if proto_code.text == i_code.text:
        return proto_code

    return i_code


def _add_used_from_code(spec, used, name):
    """ Add any interface files to a used list that are defined for a name. """

    name = ScopedName.parse(name)
    name.make_absolute()

    for iface_file in spec.iface_files:
        if iface_file.type in (IfaceFileType.CLASS, IfaceFileType.EXCEPTION):
            if iface_file.fq_cpp_name == name:
                append_iface_file(used, iface_file)
                return

    for enum in spec.enums:
        if enum.scope is not None:
            if enum.fq_cpp_name == name:
                append_iface_file(used, enum.scope.iface_file)
                return
