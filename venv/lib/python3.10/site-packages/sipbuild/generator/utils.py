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


from .scoped_name import ScopedName
from .specification import ArgumentType, CachedName, IfaceFile, IfaceFileType


def append_iface_file(iface_file_list, iface_file):
    """ Append an IfaceFile object to a list of them. """

    # Make sure we don't try and add an interface file to its own list.
    if iface_file.used is iface_file_list:
        return

    # Don't bother if it is already there.
    if iface_file in iface_file_list:
        return

    iface_file_list.append(iface_file)


def argument_as_str(arg):
    """ Convert an Argument object to a string of valid C++. """

    if arg.original_typedef is None or arg.original_typedef.no_type_name:
        if arg.type is ArgumentType.TEMPLATE:
            s = str(arg.definition.cpp_name)
            s += '<'

            need_comma = False
            for sub_arg in arg.definition.types.args:
                if need_comma:
                    s += ','
                else:
                    need_comma = True

                s += argument_as_str(sub_arg)

            if s.endswith('>'):
                # For compilers earlier than C++11.
                s += ' '

            s += '>'

        elif arg.type in (ArgumentType.STRUCT, ArgumentType.DEFINED):
            s = str(arg.definition)

        elif arg.type in (ArgumentType.UBYTE, ArgumentType.USTRING):
            s = 'unsigned char'

        elif arg.type in (ArgumentType.BYTE, ArgumentType.STRING, ArgumentType.ASCII_STRING, ArgumentType.LATIN1_STRING, ArgumentType.UTF8_STRING):
            s = 'char'

        elif arg.type in (ArgumentType.SBYTE, ArgumentType.SSTRING):
            s = 'signed char'

        elif arg.type is ArgumentType.WSTRING:
            s = 'wchar_t'

        elif arg.type is ArgumentType.USHORT:
            s = 'unsigned short'

        elif arg.type is ArgumentType.SHORT:
            s = 'short'

        elif arg.type is ArgumentType.UINT:
            s = 'uint'

        elif arg.type in (ArgumentType.INT, ArgumentType.CINT):
            s = 'int'

        elif arg.type is ArgumentType.ULONG:
            s = 'unsigned long'

        elif arg.type is ArgumentType.LONG:
            s = 'long'

        elif arg.type is ArgumentType.ULONGLONG:
            s = 'unsigned long long'

        elif arg.type is ArgumentType.LONGLONG:
            s = 'long long'

        elif arg.type in (ArgumentType.FLOAT, ArgumentType.CFLOAT):
            s = 'float'

        elif arg.type in (ArgumentType.DOUBLE, ArgumentType.CDOUBLE):
            s = 'double'

        elif arg.type in (ArgumentType.BOOL, ArgumentType.CBOOL):
            s = 'bool'

        elif arg.type is ArgumentType.VOID:
            s = 'void'

        elif arg.type is ArgumentType.CAPSULE:
            s = 'void *'

        elif arg.type is ArgumentType.SSIZE:
            s = 'Py_ssize_t'

        elif arg.type is ArgumentType.SIZE:
            s = 'size_t'

        elif arg.type is ArgumentType.HASH:
            s = 'Py_hash_t'

    else:
        # Use the original typedef.
        s = str(arg.original_typedef.fq_cpp_name)

    # Remove any global scope specifier (simply to make generated less
    # cluttered).
    if s.startswith('::'):
        s = s[2:]

    for _ in arg.derefs:
        s += '*'

    if arg.is_reference:
        s += '&'

    return s


def cached_name(spec, name):
    """ Add a name to the cache if necessary and return the cached name. """

    # Get the line of the cache for the length of this name creating it if
    # necessary.
    line = spec.name_cache.setdefault(len(name), [])

    # See if the name has already been cached.
    for nd in line:
        if nd.name == name:
            return nd

    # Create a new entry.
    nd = CachedName(name)
    line.append(nd)

    return nd


def find_iface_file(spec, mod, fq_cpp_name, iface_file_type, error_logger,
        cpp_type=None, scope=None):
    """ Return an interface file for a fully qualified C/C++ name and type
    creating it if necessary.
    """

    # See if the name is already used.
    for iff in spec.iface_files:
        if not iff.fq_cpp_name.matches(fq_cpp_name, scope=scope):
            continue

        # They must be the same type except that we allow a class if we want an
        # exception.  This is because we allow classes to be used before they
        # are defined.
        if iff.type is not iface_file_type:
            if iface_file_type is not IfaceFileType.EXCEPTION or iff.type is not IfaceFileType.CLASS:
                error_logger(
                        "a class, exception, namespace or mapped type has already been defined with the same name")

        # Ignore an external class declared in another module.
        if iface_file_type is IfaceFileType.CLASS and iff.module is not None and iff.module is not mod:
            for cd in spec.classes:
                if cd.iface_file is iff:
                    external = cd.external
                    break
            else:
                external = False

            if external:
                continue

        # If this is a mapped type with the same name defined in a different
        # module, then check that this type isn't the same as any of the mapped
        # types defined in that module.
        if iface_file_type == IfaceFileType.MAPPED_TYPE and iff.module is not mod:
            for mtd in spec.mapped_types:
                if mtd.iface_file is not iff:
                    continue

                if cpp_type.type is not ArgumentType.TEMPLATE or \
                    mtd.type.type is not ArgumentType.TEMPLATE or \
                    same_base_type(cpp_type, mtd.type):
                    error_logger(
                            "the mapped type has already been defined in another module")

            # If we got here then we have a mapped type based on an existing
            # template, but with unique parameters.  We don't want to use
            # interface files from other modules, so skip this one.
            continue

        # Ignore a namespace defined in another module.
        if iface_file_type is IfaceFileType.NAMESPACE and iff.module is not mod:
            continue

        return iff

    # Create a new interface file.
    fq_cpp_name = normalised_scoped_name(fq_cpp_name, scope)

    iff = IfaceFile(iface_file_type,
                cpp_name=cached_name(spec, str(fq_cpp_name)),
                fq_cpp_name=fq_cpp_name)

    # Use the same ordering as the old parser.
    spec.iface_files.insert(0, iff)

    return iff


def find_method(klass, name):
    """ Return the Member object for a named member of a class or None if there
    was none.
    """

    for member in klass.members:
        if member.py_name.name == name:
            return member

    return None


def normalised_scoped_name(scoped_name, scope):
    """ Convert a scoped name to a fully qualified name. """

    # Clone the name.
    fq_scoped_name = ScopedName(scoped_name)

    if fq_scoped_name.is_absolute:
        pass
    elif scope is None:
        fq_scoped_name.make_absolute()
    elif fq_scoped_name.is_simple:
        fq_scoped_name.prepend(scope.iface_file.fq_cpp_name)
    else:
        # The relative name has a scope and appears within a scope so we need
        # lookup the name's scope within the current scope.
        names_scope = fq_scoped_name[0]
        scope_fq_cpp_name = scope.iface_file.fq_cpp_name

        while scope_fq_cpp_name is not None:
            if scope_fq_cpp_name.base_name == names_scope:
                del fq_scoped_name[0]
                fq_scoped_name.prepend(scope_fq_cpp_name)
                break

            scope_fq_cpp_name = scope_fq_cpp_name.scope
        else:
            # The lookup failed so just make the name absolute.
            fq_scoped_name.make_absolute()

    return fq_scoped_name


# Argument type qualifiers.
_PY_STRING = (ArgumentType.USTRING, ArgumentType.SSTRING, ArgumentType.STRING,
        ArgumentType.ASCII_STRING, ArgumentType.LATIN1_STRING,
        ArgumentType.UTF8_STRING)
_PY_FLOAT = (ArgumentType.CFLOAT, ArgumentType.FLOAT, ArgumentType.CDOUBLE,
        ArgumentType.DOUBLE)
_PY_INT = (ArgumentType.BOOL, ArgumentType.HASH, ArgumentType.SSIZE,
        ArgumentType.SIZE, ArgumentType.BYTE, ArgumentType.SBYTE,
        ArgumentType.UBYTE, ArgumentType.SHORT, ArgumentType.USHORT,
        ArgumentType.CINT, ArgumentType.INT, ArgumentType.UINT)
_PY_LONG = (ArgumentType.LONG, ArgumentType.LONGLONG)
_PY_ULONG = (ArgumentType.ULONG, ArgumentType.ULONGLONG)
_PY_AUTO = (ArgumentType.BOOL, ArgumentType.BYTE, ArgumentType.SBYTE,
        ArgumentType.UBYTE, ArgumentType.SHORT, ArgumentType.USHORT,
        ArgumentType.INT, ArgumentType.UINT, ArgumentType.FLOAT,
        ArgumentType.DOUBLE)
_PY_CONSTRAINED = (ArgumentType.CBOOL, ArgumentType.CINT, ArgumentType.CFLOAT,
        ArgumentType.CDOUBLE)

def same_argument_type(spec, arg1, arg2, strict=True):
    """ Compare two argument types and return True if they are the same.
    'strict' means as C++ would see it, rather than Python.
    """

    if arg1.is_reference != arg2.is_reference:
        return False

    if len(arg1.derefs) != len(arg2.derefs):
        return False

    if strict:
        # The const should be the same.
        if arg1.is_const != arg2.is_const:
            return False

        return same_base_type(arg1, arg2)

    # If both are constrained fundamental types then the types must match.
    if arg1.type in _PY_CONSTRAINED and arg2.type in _PY_CONSTRAINED:
        return arg1.type is arg2.type

    if spec.abi_version >= (13, 0):
        # Anonymous enums are ints.
        if arg1.type in _PY_INT and arg2.type is ArgumentType.ENUM and arg2.definition.fq_cpp_name is None:
            return True

        if arg1.type is ArgumentType.ENUM and arg1.definition.fq_cpp_name is None and arg2.type in _PY_INT:
            return True
    else:
        # An unconstrained enum also acts as a (very) constrained int.
        if arg1.type in _PY_INT and arg2.type is ArgumentType.ENUM and not arg2.is_constrained:
            return True

        if arg1.type is ArgumentType.ENUM and not arg1.is_constrained and arg2.type in _PY_INT:
            return True

    # Python will see all these as strings.
    if arg1.type in _PY_STRING and arg2.type in _PY_STRING:
        return True

    # Python will see all these as floats.
    if arg1.type in _PY_FLOAT and arg2.type in _PY_FLOAT:
        return True

    # Python will see all these as ints.
    if arg1.type in _PY_INT and arg2.type in _PY_INT:
        return True

    # Python will see all these as longs.
    if arg1.type in _PY_LONG and arg2.type in _PY_LONG:
        return True

    # Python will see all these as unsigned longs.
    if arg1.type in _PY_ULONG and arg2.type in _PY_ULONG:
        return True

    # Python will automatically convert between these.
    if arg1.type in _PY_AUTO and arg2.type in _PY_AUTO:
        return True

    # All the special cases have been handled.
    return same_base_type(arg1, arg2)


def same_base_type(type1, type2):
    """ Return True if two Argument objects refer to the same base type, ie.
    without taking into account const and pointers.
    """

    # The types must be the same.
    if type1.type is not type2.type:
        # If we are comparing a template with those that have already been used
        # to instantiate a class or mapped type then we need to compare with
        # the class or mapped type name.

        if type1.type is ArgumentType.CLASS and type2.type is ArgumentType.DEFINED:

            return type1.definition.iface_file.fq_cpp_name.matches(type2.definition)

        if type1.type is ArgumentType.DEFINED and type2.type is ArgumentType.CLASS:
            return type2.definition.iface_file.fq_cpp_name.matches(type1.definition)

        if type1.type is ArgumentType.MAPPED and type2.type is ArgumentType.DEFINED:
            return type1.definition.iface_file.fq_cpp_name.matches(type2.definition)

        if type1.type is ArgumentType.DEFINED and type2.type is ArgumentType.MAPPED:
            return type2.definition.iface_file.fq_cpp_name.matches(type1.definition)

        if type1.type is ArgumentType.ENUM and type2.type is ArgumentType.DEFINED:
            return type1.definition.fq_cpp_name.matches(type2.definition)

        if type1.type is ArgumentType.DEFINED and type2.type is ArgumentType.ENUM:
            return type2.definition.fq_cpp_name.matches(type1.definition)

        return False

    if type1.type is ArgumentType.CLASS:
        return type1.definition is type2.definition

    if type1.type is ArgumentType.ENUM:
        return type1.definition is type2.definition

    if type1.type is ArgumentType.TEMPLATE:
        td1 = type1.definition
        td2 = type2.definition

        if td1.cpp_name.absolute != td2.cpp_name.absolute:
            return False

        if len(td1.types.args) != len(td2.types.args):
            return False

        for ad1, ad2 in zip(td1.types.args, td2.types.args):
            if len(ad1.derefs) != len(ad2.derefs):
                return False

            if not same_base_type(ad1, ad2):
                return False

        return True

    if type1.type in (ArgumentType.STRUCT, ArgumentType.UNION):
        return type1.definition == type2.definition

    if type1.type is ArgumentType.DEFINED:
        return type1.definition == type2.definition

    if type1.type is ArgumentType.MAPPED:
        return type1.definition is type2.definition

    # They must be the same if we've got this far.
    return True


def same_signature(spec, sig1, sig2, strict=True):
    """ Compare two signatures and return True if they are the same. """

    if strict:
        # Count all the arguments.
        na1 = len(sig1.args)
        na2 = len(sig2.args)
    else:
        # Count only the compulsory arguments.
        na1 = 0

        for arg in sig1.args:
            if arg.default_value is not None:
                break

            na1 += 1

        na2 = 0

        for arg in sig2.args:
            if arg.default_value is not None:
                break;

            na2 += 1

    # The number of arguments must be the same.
    if na1 != na2:
        return False

    # The arguments must be the same.
    for a in range(len(sig1.args)):
        if not strict and sig1.args[a].default_value is not None:
            break

        if not same_argument_type(spec, sig1.args[a], sig2.args[a], strict=strict):
            return False

    # Must be the same if we've got this far.
    return True


def search_typedefs(spec, cpp_name, type):
    """ Search the typedefs and update the given type from any definition. """

    # Look for the name.
    fq_cpp_name = ScopedName(cpp_name)
    fq_cpp_name.make_absolute()

    for typedef in spec.typedefs:
        if typedef.fq_cpp_name == fq_cpp_name:
            break
    else:
        return

    # Update the type.
    type.type = typedef.type.type
    type.allow_none = type.allow_none or typedef.type.allow_none
    type.definition = typedef.type.definition
    type.derefs.extend(typedef.type.derefs)
    type.disallow_none = type.disallow_none or typedef.type.disallow_none
    type.is_const = type.is_const or typedef.type.is_const
    type.is_reference = type.is_reference or typedef.type.is_reference
    type.type_hints = type.type_hints or typedef.type.type_hints

    # Remember the original typedef.
    if type.original_typedef is None:
        type.original_typedef = typedef
