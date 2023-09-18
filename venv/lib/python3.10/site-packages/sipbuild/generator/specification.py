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


from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Union

from .scoped_name import ScopedName


class AccessSpecifier(Enum):
    """ The class access specifiers. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # Private access.
    PRIVATE = 0x04

    # Protected access.
    PROTECTED = 0x02

    # Public access.
    PUBLIC = 0x01


class ArgumentType(Enum):
    """ The types of either C/C++ or Python arguments.  The numerical values of
    these can occur in generated code so so types must always be appended and
    old (unused) types must never be removed.
    """

    # The type hasn't been specified.
    NONE = 0

    # A user defined type.
    DEFINED = 1

    # A class.
    CLASS = 2

    # A struct.
    STRUCT = 3

    # A void.
    VOID = 4

    # An enum.
    ENUM = 5

    # A template.
    TEMPLATE = 6

    # No longer used.
    SIGNAL_UNUSED = 7

    # No longer used.
    SLOT_UNUSED = 8

    # No longer used.
    RXCON_UNUSED = 9

    # No longer used.
    RXDIS_UNUSED = 10

    # No longer used.
    SLOTCON_UNUSED = 11

    # No longer used.
    SLOTDIS_UNUSED = 12

    # An unsigned char.
    USTRING = 13

    # A char.
    STRING = 14

    # A short.
    SHORT = 15

    # An unsigned short.
    USHORT = 16

    # A constrained int.
    CINT = 17

    # An int.
    INT = 18

    # An unsigned int.
    UINT = 19

    # A long.
    LONG = 20

    # An unsigned long.
    ULONG = 21

    # A float.
    FLOAT = 22

    # A constrained float.
    CFLOAT = 23

    # A double.
    DOUBLE = 24

    # A constrained double.
    CDOUBLE = 25

    # A bool.
    BOOL = 26

    # A mapped type.
    MAPPED = 27

    # A Python object.
    PYOBJECT = 28

    # A Python tuple.
    PYTUPLE = 29

    # A Python list.
    PYLIST = 30

    # A Python dict.
    PYDICT = 31

    # A Python callable.
    PYCALLABLE = 32

    # A Python slice.
    PYSLICE = 33

    # No longer used.
    QOBJECT_UNUSED = 34

    # A function.
    FUNCTION = 35

    # A Python type.
    PYTYPE = 36

    # An ellipsis.
    ELLIPSIS = 37

    # A long long.
    LONGLONG = 38

    # An unsigned long long.
    ULONGLONG = 39

    # No longer used.
    ANYSLOT_UNUSED = 40

    # A constrained bool.
    CBOOL = 41

    # A signed char.
    SSTRING = 42

    # A wchar_t.
    WSTRING = 43

    # A temporary void *.
    FAKE_VOID = 44

    # A Py_ssize_t.
    SSIZE = 45

    # An ASCII encoded string.
    ASCII_STRING = 46

    # A Latin-1 encoded string.
    LATIN1_STRING = 47

    # A UTF-8 encoded string.
    UTF8_STRING = 48

    # A char used as an int.
    BYTE = 49

    # A signed char used as an int.
    SBYTE = 50

    # An unsigned char used as an unsigned int.
    UBYTE = 51

    # A Python capsule.
    CAPSULE = 52

    # A Python object that implements the buffer protocol.
    PYBUFFER = 53

    # A size_t.
    SIZE = 54

    # A Python enum.
    PYENUM = 55

    # A union.
    UNION = 56

    # A Py_hash_t.
    HASH = 57


class ArrayArgument(Enum):
    """ The array support provided by an argument. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # /Array/ was specified.
    ARRAY = 0

    # /ArraySize/ was specified.
    ARRAY_SIZE = 1

    # The argument provides no array support.
    NONE = 2


class ClassKey(Enum):
    """ The key that identifies a particular type of class. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # A class.
    CLASS = 0

    # A struct.
    STRUCT = 1

    # A union.
    UNION = 2


class DocstringFormat(Enum):
    """ The formatting applied to the text of the docstring. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # The signature is appended to the docstring.
    # Any leading spaces common to all non-blank lines in the docstring are
    # removed.
    DEINDENTED = 1

    # The docstring is used as it is specified in the .sip file.
    RAW = 0


class DocstringSignature(Enum):
    """ The position of the automatically generated function or method
    signature relative to the docstring text.  In the context of a class's
    docstring then it applies to all the class's ctors.
    """

    # TODO: Change the values to auto() once the C code has been replaced.

    # The signature is appended to the docstring.
    APPENDED = 2

    # The signature is discard.
    DISCARDED = 0

    # The signature is prepended to the docstring.
    PREPENDED = 1


class EnumBaseType(Enum):
    """ The different base types fo an enum. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # enum.Enum
    ENUM = 0

    # enum.Flag
    FLAG = 1

    # enum.IntEnum
    INT_ENUM = 2

    # enum.IntFlag
    INT_FLAG = 3

    # enum.IntEnum with unsigned values.
    UNSIGNED_INT_ENUM = 4


class GILAction(Enum):
    """ The action to take with the GIL when calling C/C++ code. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # The default action.
    DEFAULT = 0

    # Hold the GIL.
    HOLD = 1

    # Release the GIL.
    RELEASE = 2


class IfaceFileType(Enum):
    """ The type of an interface file. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # An %Exception.
    EXCEPTION = 0

    # A %MappedType.
    MAPPED_TYPE = 1

    # A namespace.
    NAMESPACE = 2

    # A class.
    CLASS = 3


class KwArgs(Enum):
    """ The level of support for passing argument as keyword arguments. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # All named arguments can be passed as keyword arguments.
    ALL = 1

    # Keyword arguments are not supported.
    NONE = 0

    # All named optional arguments (ie. those with a default value) can be
    # passed as keyword arguments.
    OPTIONAL = 2


class PyQtMethodSpecifier(Enum):
    """ The PyQt-specific method specifier. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # A signal.
    SIGNAL = 0x10

    # A slot.
    SLOT = 0x08


class PySlot(Enum):
    """ The Python slots corresponding to entries in a type object. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # tp_str
    STR = 0

    # tp_as_number.nb_int
    INT = 1

    # tp_as_number.nb_float
    FLOAT = 2

    # tp.as_mapping.mp_length and tp.as_sequence.sq_length
    LEN = 3

    # tp.as_sequence.sq_contains
    CONTAINS = 4

    # tp_as_number.nb_add
    ADD = 5

    # tp.as_sequence.sq_concat
    CONCAT = 6

    # tp_as_number.nb_subtract
    SUB = 7

    # tp_as_number.nb_multiply
    MUL = 8

    # tp.as_sequence.sq_repeat
    REPEAT = 9

    # tp_as_number.nb_remainder
    MOD = 11

    # tp_as_number.nb_floor_divide
    FLOORDIV = 12

    # tp_as_number.nb_true_divide
    TRUEDIV = 13

    # tp_as_number.nb_and
    AND = 14

    # tp_as_number.nb_or
    OR = 15

    # tp_as_number.nb_xor
    XOR = 16

    # tp_as_number.nb_lshift
    LSHIFT = 17

    # tp_as_number.nb_rshift
    RSHIFT = 18

    # tp_as_number.nb_inplace_add
    IADD = 19

    # tp.as_sequence.sq_inplace_concat
    ICONCAT = 20

    # tp_as_number.nb_inplace_subtract
    ISUB = 21

    # tp_as_number.nb_inplace_multiply
    IMUL = 22

    # tp.as_sequence.sq_inplace_repeat
    IREPEAT = 23

    # tp_as_number.nb_inplace_remainder
    IMOD = 25

    # tp_as_number.nb_inplace_floor_divide
    IFLOORDIV = 26

    # tp_as_number.nb_inplace_true_divide
    ITRUEDIV = 27

    # tp_as_number.nb_inplace_and
    IAND = 28

    # tp_as_number.nb_inplace_or
    IOR = 29

    # tp_as_number.nb_inplace_xor
    IXOR = 30

    # tp_as_number.nb_inplace_lshift
    ILSHIFT = 31

    # tp_as_number.nb_inplace_rshift
    IRSHIFT = 32

    # tp_as_number.nb_invert
    INVERT = 33

    # tp_call
    CALL = 34

    # tp.as_mapping.mp_subscript and tp.as_sequence.sq_item
    GETITEM = 35

    # tp.as_mapping.mp_ass_subscript and tp.as_sequence.sq_ass_item
    SETITEM = 36

    # tp.as_mapping.mp_ass_subscript and tp.as_sequence.sq_ass_item
    DELITEM = 37

    # tp_richcompare (Py_LT)
    LT = 38

    # tp_richcompare (Py_LE)
    LE = 39

    # tp_richcompare (Py_EQ)
    EQ = 40

    # tp_richcompare (Py_NE)
    NE = 41

    # tp_richcompare (Py_GT)
    GT = 42

    # tp_richcompare (Py_GE)
    GE = 43

    # tp_as_number.nb_bool
    BOOL = 45

    # tp_as_number.nb_negative
    NEG = 46

    # tp_as_number.nb_positive
    POS = 47

    # tp_as_number.nb_absolute
    ABS = 48

    # tp_repr
    REPR = 49

    # tp_hash
    HASH = 50

    # tp_as_number.nb_index
    INDEX = 51

    # tp_iter
    ITER = 52

    # tp_iter_next
    NEXT = 53

    # tp_setattr
    SETATTR = 54

    # Internal to the parser (implemented as tp_setattr)
    DELATTR = 55

    # tp_as_number.nb_matrix_multiply
    MATMUL = 56

    # tp_as_number.nb_inplace_matrix_multiply
    IMATMUL = 57

    # tp_as_async.am_await
    AWAIT = 58

    # tp_as_async.am_aiter
    AITER = 59

    # tp_as_async.am_anext
    ANEXT = 60


class QualifierType(Enum):
    """ The type of a qualifier used in %If/%End directives. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # The qualifier is a feature.
    FEATURE = 2

    # The qualifier is a platform.
    PLATFORM = 1

    # The qualifier is part of a timeline.
    TIME = 0


class Transfer(Enum):
    """ The different types of ownership transfer. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # No transfer of ownership.
    NONE = 0

    # /Transfer/ was specified.
    TRANSFER = 1

    # /TransferBack/ was specified.
    TRANSFER_BACK = 2

    # /TransferThis/ was specified.
    TRANSFER_THIS = 3


class ValueType(Enum):
    """ The different types of a value in an expression. """

    # TODO: Change the values to auto() once the C code has been replaced.

    # A quoted character.
    QCHAR = 0

    # A string.
    STRING = 1

    # A number.
    NUMERIC = 2

    # A floating point number.
    REAL = 3

    # A scoped name.
    SCOPED = 4

    # A function call.
    FCALL = 5

    # A placeholder.
    EMPTY = 6


@dataclass
class Argument:
    """ Encapsulate a callable argument (or return value or variable type). """

    # The type.
    type: ArgumentType

    # Set if /AllowNone/ was specified.
    allow_none: bool = False

    # The support the argument provides for arrays.
    array: ArrayArgument = ArrayArgument.NONE

    # The optional default value.
    default_value: Optional[List['Value']] = None

    # The optional definition.  What this is depends on the type.
    definition: Any = None

    # The sequence of dereferences.  An element is True if the corresponding
    # dereference is const.
    derefs: List[bool] = field(default_factory=list)

    # Set if /DisallowNone/ was specified.
    disallow_none: bool = False

    # Set if /GetWrapper/ was specified.
    get_wrapper: bool = False

    # Set if the argument is const.
    is_const: bool = False

    # Set if /Constrained/ was specified.
    is_constrained: bool = False

    # Set if the argument is passing a value in.
    is_in: bool = False

    # Set if the argument is passing a value out.
    is_out: bool = False

    # Set if the argument is a reference.
    is_reference: bool = False

    # The key for a reference.  If None then /KeepReference/ wasn't specified.
    key: Optional[int] = None

    # The optional name.
    name: Optional['CachedName'] = None

    # Set if /NoCopy/ was specified.
    no_copy: bool = False

    # The original type if it was a typedef.
    original_typedef: Optional['WrappedTypedef'] = None

    # Set if /ResultSize/ was specified.
    result_size: bool = False

    # The value of /ScopesStripped/.
    scopes_stripped: int = 0

    # The source location.
    source_location: Optional['SourceLocation'] = None

    # Any transfer of ownership.
    transfer: Transfer = Transfer.NONE

    # The non-default type hints.
    type_hints: Optional['TypeHints'] = None


@dataclass
class CachedName:
    """ Encapsulate a name that may be needed as a string in the generated
    code.
    """

    # The name.
    name: str

    # Set if the name is a substring of another. (resolver)
    is_substring: bool = False

    # The offset of the name in the string pool. (resolver)
    offset: int = 0

    # Set if the name is used in the generated code.
    used: bool = False

    def __str__(self):
        """ Return the string representation. """

        return self.name


@dataclass
class CodeBlock:
    """ Encapsulate a code block, ie. the literal body of a directive. """

    # The name of the .sip file containing the code.
    sip_file: str

    # The line number in the .sip file that the code starts at.
    line_nr: int = 1

    # The text of the code block.
    text: str = ''


@dataclass
class Constructor:
    """ Encapsulate a constructor. """

    # The access specifier.
    access_specifier: AccessSpecifier

    # The Python signature.
    py_signature: 'Signature'

    # The C/C++ signature.  It will be none if /NoDerived/ was specified.
    cpp_signature: Optional['Signature'] = None

    # Set if /Deprecated/ was specified.
    deprecated: bool = False

    # The docstring.
    docstring: Optional['Docstring'] = None

    # The action required on the GIL.
    gil_action: GILAction = GILAction.DEFAULT

    # Set if the ctor is the implementation of a cast. (resolver)
    is_cast: bool = False

    # The keyword argument support.
    kw_args: KwArgs = KwArgs.NONE

    # The code specified by any %MethodCode directive.
    method_code: Optional[CodeBlock] = None

    # Set if the type hint should be suppressed.
    no_type_hint: bool = False

    # The /PostHook/ name.
    posthook: Optional[str] = None

    # The /PreHook/ name.
    prehook: Optional[str] = None

    # The code specified by any %PreMethodCode directive.
    premethod_code: Optional[CodeBlock] = None

    # Set if a Python exception is raised.
    raises_py_exception: bool = False

    # The optional throw arguments.
    throw_args: Optional['ThrowArguments'] = None

    # Any transfer of ownership.
    transfer: Transfer = Transfer.NONE


@dataclass
class Docstring:
    """ Encapsulate a docstring. """

    # The position of the automatically generated signature.
    signature: DocstringSignature

    # The text of the docstring.
    text: str


@dataclass
class Extract:
    """ Encapsulate a part of an extract. """

    # The ID of the extract.
    id: str

    # The order.  A negative value implies the part is appended to the extract.
    order: int

    # The text of the extract part.
    text: str


@dataclass
class FunctionCall:
    """ Encapsulate a call to a function in an expression. """

    # The type of the result.
    result: Argument

    # The list of arguments.
    args: List[List['Value']]


@dataclass
class IfaceFile:
    """ Encapsulate an interface file, ie. a generated source file. """

    # The type.
    type: IfaceFileType

    # The C/C++ name.  It will be None if the interface file relates to a
    # template.  Note that this is fully qualified and so overlaps with
    # 'fq_cpp_name'.
    cpp_name: Optional[CachedName] = None

    # The filename extension.
    file_extension: Optional[str] = None

    # The fully qualified C/C++ name.  It will be None if the interface file
    # relates to a template.
    fq_cpp_name: Optional[ScopedName] = None

    # The generated type number. (resolver)
    type_nr: int = -1

    # The defining module.  It will be None if the interface file relates to a
    # template.
    module: Optional['Module'] = None

    # Set if this interface file is needed by the module for which code is to
    # be generated.
    needed: bool = False

    # The %TypeHeaderCode.
    type_header_code: List[CodeBlock] = field(default_factory=list)

    # The interface files used by this one (either directly or indirectly).
    used: List['IfaceFile'] = field(default_factory=list)


@dataclass
class License:
    """ Encapsulate a license. """

    # The type of the license.
    type: str

    # The licensee.
    licensee: Optional[str] = None

    # The timestamp.
    timestamp: Optional[str] = None

    # The signature.
    signature: Optional[str] = None


@dataclass
class MappedType:
    """ Encapsulate a mapped type. """

    # The interface file.
    iface_file: IfaceFile

    # The type.
    type: Argument

    # The %ConvertFromTypeCode.
    convert_from_type_code: Optional[CodeBlock] = None

    # The %ConvertToTypeCode.
    convert_to_type_code: Optional[CodeBlock] = None

    # The C/C++ name.  It will be None for mapped type templates.
    cpp_name: Optional[CachedName] = None

    # Set if /AllowNone/ was specified.
    handles_none: bool = False

    # The %InstanceCode.
    instance_code: Optional[CodeBlock] = None

    # The member functions.
    members: List['Member'] = field(default_factory=list)

    # Set if the handwritten code requires user state information.
    needs_user_state: bool = False

    # Set if /NoAssignmentOperator/ was specified.
    no_assignment_operator: bool = False

    # Set if /NoCopyCtor/ was specified.
    no_copy_ctor: bool = False

    # Set if /NoDefaultCtor/ was specified.
    no_default_ctor: bool = False

    # Set if /NoRelease/ was specified.
    no_release: bool = False

    # The overloaded member functions.
    overloads: List['Overload'] = field(default_factory=list)

    # The Python name.  It will be None for mapped type templates.
    py_name: Optional[CachedName] = None

    # The /PyQtFlags/.
    pyqt_flags: int = 0

    # The %ReleaseCode.
    release_code: Optional[CodeBlock] = None

    # The %TypeCode.
    type_code: List[CodeBlock] = field(default_factory=list)

    # The %TypeHintCode.
    type_hint_code: Optional[CodeBlock] = None

    # The type hints.
    type_hints: Optional['TypeHints'] = None


@dataclass
class MappedTypeTemplate:
    """ Encapsulate a mapped type template. """

    # The prototype mapped type.
    mapped_type: MappedType

    # The template arguments.
    signature: 'Signature'


@dataclass
class Member:
    """ Encapsulate a member function. """

    # The defining module.
    module: 'Module'

    # The Python name.
    py_name: CachedName

    # Set if keyword arguments are allowed.
    allow_keyword_args: bool = False

    # Set if at least one of the overloads is protected.
    has_protected: bool = False

    # Set if /Numeric/ was specified.
    is_numeric: bool = False

    # Set if /Sequence/ was specified.
    is_sequence: bool = False

    # The original interface file if the function was defined in a namespace.
    namespace_iface_file: Optional[IfaceFile] = None

    # Set if /NoArgParser/ was specified.
    no_arg_parser: bool = False

    # The Python slot if it is not an ordinary member function.
    py_slot: Optional[PySlot] = None


@dataclass
class Module:
    """ Encapsulate a module. """

    # The list of all (explicit and implied) imports. (resolver)
    all_imports: List['Module'] = field(default_factory=list)

    # Set if wrapped ctors should support cooperative multi-inheritance.
    call_super_init: bool = False

    # The text specified by any %Copying directives.
    copying: List[CodeBlock] = field(default_factory=list)

    # The default docstring format.
    default_docstring_format: DocstringFormat = DocstringFormat.RAW

    # The default docstring signature position.
    default_docstring_signature: DocstringSignature = DocstringSignature.DISCARDED

    # The default exception.
    default_exception: Optional['WrappedException'] = None

    # The default meta-type.
    default_metatype: Optional[CachedName] = None

    # The default super-type.
    default_supertype: Optional[CachedName] = None

    # The handler called when a Python re-implementation of a virtual C++
    # function raises an exception.
    default_virtual_error_handler: Optional[str] = None

    # The module's docstring.
    docstring: Optional[Docstring] = None

    # The fully qualified name of the module.  It is only None until %Module is
    # specified.
    fq_py_name: Optional[CachedName] = None

    # The global functions.
    global_functions: List[Member] = field(default_factory=list)

    # Set if any class defined in the module has a delayed dtor.
    has_delayed_dtors: bool = False

    # The list of direct imports.
    imports: List['Module'] = field(default_factory=list)

    # The code specified by any %InitialisationCode directives.
    initialisation_code: List[CodeBlock] = field(default_factory=list)

    # The software license.
    license: Optional[License] = None

    # The %ModuleCode.
    module_code: List[CodeBlock] = field(default_factory=list)

    # The %ModuleHeaderCode.
    module_header_code: List[CodeBlock] = field(default_factory=list)

    # The next key to auto-allocate.
    next_key: int = -1

    # The number of exceptions defined in this module. (resolver)
    nr_exceptions: int = 0

    # The generated types needed by this module. (resolver)
    needed_types: List[Argument] = field(default_factory=list)

    # The number of typedefs defined in this module.
    nr_typedefs: int = 0

    # The number of virtual error handlers defined in this module. (resolver)
    nr_virtual_error_handlers: int = 0

    # The overloaded global functions.
    overloads: List['Overload'] = field(default_factory=list)

    # The code specified by any %PostInitialisationCode directives.
    postinitialisation_code: List[CodeBlock] = field(default_factory=list)

    # The proxy classes.
    proxies: List['WrappedClass'] = field(default_factory=list)

    # The code specified by any %PreInitialisationCode directives.
    preinitialisation_code: List[CodeBlock] = field(default_factory=list)

    # The name of the module. (resolver)
    py_name: Optional[str] = None

    # Set if the generated bindings are Py_ssize_t clean.
    py_ssize_t_clean: bool = False

    # The list of qualifiers.
    qualifiers: List['Qualifier'] = field(default_factory=list)

    # The %TypeHintCode.
    type_hint_code: List[CodeBlock] = field(default_factory=list)

    # The %UnitCode.
    unit_code: List[CodeBlock] = field(default_factory=list)

    # The %UnitPostIncludeCode.
    unit_postinclude_code: List[CodeBlock] = field(default_factory=list)

    # Set if the actual argument names to wrapped callables should be used in
    # the generated bindings rather than automatically generated ones.
    use_arg_names: bool = False

    # Set if the generated bindings should only use the limited Python API.
    use_limited_api: bool = False

    # The interface files used by the module.
    used: List[IfaceFile] = field(default_factory=list)


@dataclass
class Overload:
    """ Encapsulate an overloaded member function. """

    # The access specifier.
    access_specifier: Optional[AccessSpecifier]

    # The member common to all overloads.
    common: Member

    # The C/C++ name.
    cpp_name: str

    # The C/C++ signature.
    cpp_signature: 'Signature'

    # The Python signature.
    py_signature: 'Signature'

    # Set if /AbortOnException/ is specified.
    abort_on_exception: bool = False

    # Set if the overload is really protected.
    access_is_really_protected: bool = False

    # Set if /Deprecated/ was specified.
    deprecated: bool = False

    # The docstring.
    docstring: Optional[Docstring] = None

    # Set if /Factory/ was specified.
    factory: bool = False

    # The action required on the GIL.
    gil_action: GILAction = GILAction.DEFAULT

    # Set if the overload is abstract.
    is_abstract: bool = False

    # Set if /AutoGen/ was specified and the associated feature was enabled.
    is_auto_generated: bool = False

    # Set if the overload is a complementary slot. (resolver)
    is_complementary: bool = False

    # Set if the overload is const.
    is_const: bool = False

    # Set if the overload implements __delattr__ (as opposed to __setattr__).
    is_delattr: bool = False

    # Set if the overload is final.
    is_final: bool = False

    # Set if the C++ overload is global. (resolver)
    is_global: bool = False

    # Set if self should not be dereferenced. (resolver)
    dont_deref_self: bool = False

    # Set if the overload is a reflected slot. (resolver)
    is_reflected: bool = False

    # Set if the overload is static.
    is_static: bool = False

    # Set if the overload is virtual.
    is_virtual: bool = False

    # Set if the overload is a virtual reimplementation. (resolver)
    is_virtual_reimplementation: bool = False

    # The keyword argument support.
    kw_args: KwArgs = KwArgs.NONE

    # The code specified by any %MethodCode directive.
    method_code: Optional[CodeBlock] = None

    # Set if /NewThread/ was specified.
    new_thread: bool = False

    # Set if the type hint should be suppressed.
    no_type_hint: bool = False

    # Set if any virtual error handler should be ignored.
    no_virtual_error_handler: bool = False

    # The /PostHook/ name.
    posthook: Optional[str] = None

    # The /PreHook/ name.
    prehook: Optional[str] = None

    # The code specified by any %PreMethodCode directive.
    premethod_code: Optional[CodeBlock] = None

    # The PyQt method specifier.
    pyqt_method_specifier: Optional[PyQtMethodSpecifier] = None

    # Set if a Python exception is raised.
    raises_py_exception: bool = False

    # The source location.
    source_location: Optional['SourceLocation'] = None

    # The optional throw arguments.
    throw_args: Optional['ThrowArguments'] = None

    # Any transfer of ownership.
    transfer: Transfer = Transfer.NONE

    # The code specified by any %VirtualCallCode directive.
    virtual_call_code: Optional[CodeBlock] = None

    # The code specified by any %VirtualCatcherCode directive.
    virtual_catcher_code: Optional[CodeBlock] = None

    # The name of the virtual error handler to use.
    virtual_error_handler: Optional[str] = None


@dataclass
class Property:
    """ Encapsulate a property. """

    # The name of the getter.
    getter: str

    # The name.
    name: CachedName

    # The docstring.
    docstring: Optional[Docstring] = None

    # The name of the optional setter.
    setter: Optional[str] = None


@dataclass
class Qualifier:
    """ Encapsulate a qualifier used in %If/%End directives. """

    # The defining module.
    module: 'Module'

    # The name of the qualifier.
    name: str

    # The type of the qualifier
    type: QualifierType

    # Set if the qualifier is enabled by default.
    enabled_by_default: bool = False

    # The order if it is a TIME qualifier (ie. the position within the
    # timeline).
    order: int = 0

    # The timeline number within the defining module if it is a TIME qualifier.
    timeline: int = 0


@dataclass
class Signature:
    """ Encapsulate a signature (including the optional result). """

    # The list of arguments.
    args: List[Argument] = field(default_factory=list)

    # The type of the optional result.
    result: Optional[Argument] = None


@dataclass
class SourceLocation:
    """ Encapsulate a location in a .sip source file. """

    # The .sip file name.
    sip_file: str

    # The column number.
    column: int = 0

    # The line number.
    line: int = 0


@dataclass
class Specification:
    """ Encapsulate a parsed .sip file. """

    # The version of the ABI being targeted.
    abi_version: tuple

    # Set if the specification is strict.
    is_strict: bool

    # The fully qualified name of the sip module.  If it is None then there is
    # no shared sip module.
    sip_module: Optional[str]

    # Set if the bindings are for C rather than C++.
    c_bindings: bool = False

    # The list of classes.
    classes: List['WrappedClass'] = field(default_factory=list)

    # The list of enums.
    enums: List['WrappedEnum'] = field(default_factory=list)

    # The list of exceptions.
    exceptions: List['WrappedException'] = field(default_factory=list)

    # The %ExportedHeaderCode.
    exported_header_code: List[CodeBlock] = field(default_factory=list)

    # The %ExportedTypeHintCode.
    exported_type_hint_code: List[CodeBlock] = field(default_factory=list)

    # The extracts.
    extracts: List[Extract] = field(default_factory=list)

    # The interface files.
    iface_files: List[IfaceFile] = field(default_factory=list)

    # Set if the specification is for a composite module.
    is_composite: bool = False

    # The mapped type templates.
    mapped_type_templates: List[MappedTypeTemplate] = field(default_factory=list)

    # The mapped types.
    mapped_types: List[MappedType] = field(default_factory=list)

    # The module for which code is to be generated.
    module: Module = field(default_factory=Module)

    # The cache of names that may be required as strings in the generated code.
    name_cache: Dict[int, List[CachedName]] = field(default_factory=dict)

    # The number of virtual handlers. (resolver)
    nr_virtual_handlers: int = 0

    # The list of plugins.  Note that these are PyQt-specific and will be
    # removed in SIP v7.
    plugins: List[str] = field(default_factory=list)

    # The QObject class.
    pyqt_qobject: Optional['WrappedClass'] = None

    # The list of typedefs.
    typedefs: List['WrappedTypedef'] = field(default_factory=list)

    # The list of variables.
    variables: List['WrappedVariable'] = field(default_factory=list)

    # The list of virtual error handlers.
    virtual_error_handlers: List['VirtualErrorHandler'] = field(default_factory=list)

    # The list of virtual handlers. (resolver)
    virtual_handlers: List['VirtualHandler'] = field(default_factory=list)

    def __hash__(self):
        """ Reimplemented so a Specification object can be used as a dict key.
        """

        return id(self)


@dataclass
class Template:
    """ Encapsulate a template. """

    # The C++ name.
    cpp_name: ScopedName

    # The types.
    types: Signature


@dataclass
class ThrowArguments:
    """ Encapsulate the arguments to a C++ throw(). """

    # The list of the argument names. If it is None then 'noexcept' was
    # specified.
    arguments: Optional[List['WrappedException']] = None


@dataclass
class TypeHints:
    """ Encapsulate a set of PEP 484 type hints for a type. """

    # The type hint when used to pass a value into a callable.
    hint_in: Optional[str]

    # The type hint used to return a value from a callable.
    hint_out: Optional[str]

    # The representation of a default value in a type hint.
    default_value: Optional[str]


@dataclass
class Value:
    """ Encapsulate a literal value. """

    # The type of the value.
    value_type: ValueType

    # Any literal value.
    value: Optional[Union[str, int, float, FunctionCall, ScopedName]]

    # Any binary operator.
    binary_operator: Optional[str] = None

    # Any cast.
    cast: Optional[ScopedName] = None

    # Any unary operator.
    unary_operator: Optional[str] = None


@dataclass
class VirtualErrorHandler:
    """ Encapsulate a virtual error handler. """

    # The code implementing the handler.
    code: CodeBlock

    # The defining module.
    module: Module

    # The name of the handler.
    name: str

    # The number of the handler. (resolver)
    handler_nr: int = -1


@dataclass
class VirtualHandler:
    """ Encapsulate a virtual overload handler. (resolver) """

    # The C/C++ signature.
    cpp_signature: Signature

    # The Python signature.
    py_signature: Signature

    # The code specified by any %VirtualCatcherCode directive.
    virtual_catcher_code: Optional[CodeBlock]

    # The virtual error handler.
    virtual_error_handler: VirtualErrorHandler

    # Set if execution should abort if there is an exception.
    abort_on_exception: bool = False

    # The number of the handler.
    handler_nr: int = -1

    # Set if ownership of the result should be transferred.
    transfer_result: bool = False


@dataclass
class VirtualOverload:
    """ Encapsulate a virtual overloaded member function. (resolver) """

    # The overload
    overload: Overload

    # The handler for the overload.  It is only set for the module for which
    # code is being generated.
    handler: Optional[VirtualHandler]


@dataclass
class VisibleMember:
    """ Encapsulate a visible member function. (resolver) """

    # The member function.
    member: Member

    # The defining class.
    scope: 'WrappedClass'


@dataclass
class WrappedClass:
    """ Encapsulate a wrapped C/C++ namespace/class/struct/union. """

    # The interface file.
    iface_file: IfaceFile

    # The Python name.
    py_name: CachedName

    # The enclosing scope.
    scope: Optional['WrappedClass']

    # The %BIGetBufferCode.
    bi_get_buffer_code: Optional[CodeBlock] = None

    # The %BIReleaseBufferCode.
    bi_release_buffer_code: Optional[CodeBlock] = None

    # Set if the class has usable constructors.
    can_create: bool = False

    # Set if an instance of the class cannot be assigned.
    cannot_assign: bool = False

    # Set if an instance of the class cannot be copied. (resolver)
    cannot_copy: bool = False

    # The list of operator casts.
    casts: List[Argument] = field(default_factory=list)

    # The specific type of class.  It will be None for namespaces.
    class_key: Optional[ClassKey] = None

    # The %ConvertFromTypeCode.
    convert_from_type_code: Optional[CodeBlock] = None

    # The %ConvertToSubClassCode.
    convert_to_subclass_code: Optional[CodeBlock] = None

    # The %ConvertToTypeCode.
    convert_to_type_code: Optional[CodeBlock] = None

    # The constructors.
    ctors: List[Constructor] = field(default_factory=list)

    # The dtor's %PreMethodCode and %MethodCode.
    dealloc_code: List[CodeBlock] = field(default_factory=list)

    # The constructor that has /Default/ specified.
    default_ctor: Optional[Constructor] = None

    # Set if /DelayDtor/ was specified.
    delay_dtor: bool = False

    # Set if /Deprecated/ was specified.
    deprecated: bool = False

    # The docstring.
    docstring: Optional[Docstring] = None

    # The access specifier of any dtor.
    dtor: Optional[AccessSpecifier] = None

    # The action required on the GIL.
    dtor_gil_action: GILAction = GILAction.DEFAULT

    # The optional dtor throw arguments.
    dtor_throw_args: Optional[ThrowArguments] = None

    # The code specified by any dtor %VirtualCatcherCode directive.
    dtor_virtual_catcher_code: Optional[CodeBlock] = None

    # Set if /ExportDerived/ was specified.
    export_derived: bool = False

    # Set if /External/ was specified.
    external: bool = False

    # The %FinalisationCode.
    finalisation_code: Optional[CodeBlock] = None

    # The %GCClearCode.
    gc_clear_code: Optional[CodeBlock] = None

    # The %GCTraverseCode.
    gc_traverse_code: Optional[CodeBlock] = None

    # Set if /AllowNone/ was specified.
    handles_none: bool = False

    # Set if the class has a non-lazy method.
    has_nonlazy_method: bool = False

    # Set if the class actually has a shadow (ie. derived) class. (resolver)
    has_shadow: bool = False

    # Set if the class has variables that need handlers. (resolver)
    has_variable_handlers: bool = False

    # The %InstanceCode.
    instance_code: Optional[CodeBlock] = None

    # Set if the class is abstract.
    is_abstract: bool = False

    # Set if the class is a hidden namespace.
    is_hidden_namespace: bool = False

    # Set if the class is incomplete.
    is_incomplete: bool = False

    # Set if the class is opaque.
    is_opaque: bool = False

    # Set if the class is defined in a protected section.
    is_protected: bool = False

    # Set if the class is QObject or a sub-class. (resolver)
    is_qobject: bool = False

    # The methods.
    members: List[Member] = field(default_factory=list)

    # The value of /Metatype/ if specified.
    metatype: Optional[CachedName] = None

    # Set if /Mixin/ was specified.
    mixin: bool = False

    # The list of all classes in the class hierarchy starting with itself.
    # (resolver)
    mro: List['WrappedClass'] = field(default_factory=list)

    # Set if the class needs an array helper. (resolver)
    needs_array_helper: bool = False

    # Set if the class needs a copy helper. (resolver)
    needs_copy_helper: bool = False

    # Set if the class needs a shadow (ie. derived) class.
    needs_shadow: bool = False

    # Set if /NoDefaultCtors/ was specified.
    no_default_ctors: bool = False

    # Set if /NoTypeHint/ was specified.
    no_type_hint: bool = False

    # Set if the class name should not be used in the generated code (and the
    # instantiated template name should be used instead).
    no_type_name: bool = False

    # The overloaded methods.
    overloads: List[Overload] = field(default_factory=list)

    # The %PickleCode.
    pickle_code: Optional[CodeBlock] = None

    # The properties.
    properties: List[Property] = field(default_factory=list)

    # The /PyQtFlags/.
    pyqt_flags: int = 0

    # The /PyQtFlagsEnums/.
    pyqt_flags_enums: Optional[List[str]] = None

    # The /PyQtInterface/.
    pyqt_interface: Optional[str] = None

    # Set if /PyQtNoQMetaObject/ was specified.
    pyqt_no_qmetaobject: bool = False

    # The real class if this is a proxy or a namespace extender.
    real_class: Optional['WrappedClass'] = None

    # The sub-class base class. (resolver)
    subclass_base: Optional['WrappedClass'] = None

    # The super-classes.
    superclasses: List['WrappedClass'] = field(default_factory=list)

    # The value of /Supertype/ if specified.
    supertype: Optional[CachedName] = None

    # The template that was instantiated to create this class.
    template: Optional[Template] = None

    # The %TypeCode.
    type_code: List[CodeBlock] = field(default_factory=list)

    # The %TypeHintCode.
    type_hint_code: Optional[CodeBlock] = None

    # The type hints.
    type_hints: Optional[TypeHints] = None

    # The name of the virtual error handler to use.
    virtual_error_handler: Optional[str] = None

    # The virtual overloaded methods. (resolver)
    virtual_overloads: List[VirtualOverload] = field(default_factory=list)

    # The visible member functions. (resolver)
    visible_members: List[VisibleMember] = field(default_factory=list)


@dataclass
class WrappedEnum:
    """ Encapsulate a wrapped enum. """

    # The base type.
    base_type: EnumBaseType

    # The fully qualified C++ name.
    fq_cpp_name: Optional[ScopedName]

    # The defining module.
    module: Module

    # The cached fully qualified C++ name.
    cached_fq_cpp_name: Optional[CachedName] = None

    # Set if the enum is defined in the protected section.
    is_protected: bool = False

    # Set if the enum is a scoped enum.
    is_scoped: bool = False

    # The members.
    members: List['WrappedEnumMember'] = field(default_factory=list)

    # Set if this enum is needed by the module for which code is to be
    # generated. (resolver)
    needed: bool = False

    # Set if /NoScope/ was specified.
    no_scope: bool = False

    # Set if the type hint should be suppressed.
    no_type_hint: bool = False

    # The overloaded slot member functions. (resolver)
    overloads: List['Overload'] = field(default_factory=list)

    # The Python name.
    py_name: Optional[CachedName] = None

    # The enclosing scope.
    scope: Optional[Union[MappedType, WrappedClass]] = None

    # The slot member functions.  These can only be created by global operators
    # being moved. (resolver)
    slots: List[Member] = field(default_factory=list)

    # The generated type number. (resolver)
    type_nr: int = -1


@dataclass
class WrappedEnumMember:
    """ Encapsulate a member of a wrapped enum. """

    # The C++ name.
    cpp_name: str

    # The Python name.
    py_name: CachedName

    # The enclosing enum.
    scope: 'WrappedEnum'

    # Set if the type hint should be suppressed.
    no_type_hint: bool = False


@dataclass
class WrappedException:
    """ Encapsulate a wrapped exception. """

    # The interface file.
    iface_file: IfaceFile

    # The code specified by the %RaiseCode directive.
    raise_code: CodeBlock

    # The base exception if it is a builtin.
    builtin_base_exception: Optional[str] = None

    # The class that implements the exception (if the exception is not a Python
    # exception).
    class_exception: Optional[WrappedClass] = None

    # The base exception if it is defined in the specification.
    defined_base_exception: Optional['WrappedException'] = None

    # The number of the exception (only if a Python exception object is
    # required. (resolver)
    exception_nr: int = -1

    # Set if this exception is needed by the module for which code is to be
    # generated. (resolver)
    needed: bool = False

    # The Python name.
    py_name: Optional[str] = None


@dataclass
class WrappedTypedef:
    """ Encapsulate a wrapped typedef. """

    # The fully qualified C++ name.
    fq_cpp_name: ScopedName

    # The defining module.
    module: Module

    # The enclosing scope.
    scope: Optional[WrappedClass]

    # The type.
    type: Argument

    # Set if the typedef name should not be used in the generated code.
    no_type_name: bool = False


@dataclass
class WrappedVariable:
    """ Encapsulate a wrapped variable. """

    # The fully qualified C++ name.
    fq_cpp_name: ScopedName

    # The defining module.
    module: Module

    # The Python name.
    py_name: CachedName

    # The enclosing scope.
    scope: Optional[WrappedClass]

    # The type.
    type: Argument

    # The code specified by any %AccessCode directive.
    access_code: Optional[CodeBlock] = None

    # The code specified by any %GetCode directive.
    get_code: Optional[CodeBlock] = None

    # Set if the variable is static.
    is_static: bool = False

    # Set if the variable needs a handler. (resolver)
    needs_handler: bool = False

    # Set if the type hint should be suppressed.
    no_type_hint: bool = False

    # Set if the variable has no setter and will be read-only.
    no_setter: bool = False

    # The code specified by any %SetCode directive.
    set_code: Optional[CodeBlock] = None
