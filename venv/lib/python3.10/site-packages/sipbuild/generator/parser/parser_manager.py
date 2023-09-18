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


from functools import partial
import os

from ply import lex, yacc

from ...bindings_configuration import get_bindings_configuration
from ...exceptions import UserException

from ..error_log import ErrorLog
from ..instantiations import instantiate_class
from ..python_slots import invalid_global_slot, slot_name_detail_map
from ..scoped_name import ScopedName
from ..specification import (AccessSpecifier, Argument, ArgumentType,
        ArrayArgument, CachedName, ClassKey, CodeBlock, Constructor,
        DocstringFormat, DocstringSignature, EnumBaseType, GILAction,
        IfaceFile, IfaceFileType, KwArgs, MappedType, Member, Module, Overload,
        PyQtMethodSpecifier, PySlot, Qualifier, QualifierType, Signature,
        SourceLocation, Specification, Transfer, TypeHints, WrappedClass,
        WrappedException, WrappedEnum, WrappedEnumMember)
from ..templates import encoded_template_name, same_template_signature
from ..utils import (argument_as_str, cached_name, find_iface_file,
        normalised_scoped_name, same_base_type)

from . import rules
from . import tokens
from .annotations import InvalidAnnotation, validate_annotation_value


class ParserManager:
    """ This object manages the actual lexer and parser objects providing them
    with state and utility functions.
    """

    def __init__(self, hex_version, encoding, abi_version, tags,
            disabled_features, protected_is_public, include_dirs, sip_module,
            is_strict):
        """ Initialise the manager. """

        # Create the lexer.
        self._lexer = lex.lex(module=tokens)
        self._lexer.pm = self

        # Create the parser.
        self._parser = yacc.yacc(module=rules, debug=False)
        self._parser.pm = self

        # This is a hack to give p_error() access to the current parser object.
        rules.parser = self._parser

        # The list of class templates.  Each element is a 2-tuple of the
        # template arguments and the class itself.
        self.class_templates = []

        # Public state.
        self.tags = tags

        self.spec = Specification(
                tuple([int(v) for v in abi_version.split('.')]), is_strict,
                sip_module)

        # The module is initially unnamed.
        self.modules = [self.spec.module]

        self.c_bindings = None
        self.code_block = None
        self.module_state = None
        self.paren_depth = 0
        self.parsing_template = False
        self.parsing_virtual = False
        self.raw_sip_file = None
        self.skip_stack = [False]

        # Private state.
        self._hex_version = hex_version
        self._encoding = encoding
        self._disabled_features = disabled_features
        self._protected_is_public = protected_is_public
        self._include_dirs = include_dirs
        self._template_arg_classes = []

        self._scope_stack = []
        self._error_log = ErrorLog()
        self._file_stack = []
        self._pending_module_state = None
        self._input = None
        self._all_sip_files = []
        self._sip_file = None
        self._sip_files = []

    def complete_class(self, p, symbol, annotations, has_body):
        """ Complete the definition of the class that is the current scope, pop
        it and return it.
        """

        klass = self.scope

        if has_body:
            if klass.scope is not None:
                if klass.iface_file.fq_cpp_name.scope != klass.scope.iface_file.fq_cpp_name:
                    self.parser_error(p, symbol,
                            "a scoped name cannot be specified in the definition of a class/struct/union")
        elif len(klass.superclasses) != 0:
            self.parser_error(p, symbol,
                    "the class/struct has super-classes but no definition");
        else:
            klass.is_opaque = True

        # Get the Python name and see if it is different to the C++ name.
        py_name = self.get_py_name(klass.iface_file.fq_cpp_name.base_name,
                annotations)
        klass.py_name = cached_name(self.spec, py_name)

        klass.no_type_hint = annotations.get('NoTypeHint', False)

        metatype = annotations.get('Metatype')
        if metatype is not None:
            klass.metatype = cached_name(self.spec, metatype)

        supertype = annotations.get('Supertype')
        if supertype is not None:
            klass.supertype = cached_name(self.spec, supertype)

        klass.export_derived = annotations.get('ExportDerived', False)
        klass.mixin = annotations.get('Mixin', False)

        file_extension = annotations.get('FileExtension')
        if file_extension is not None:
            klass.iface_file.file_extension = file_extension

        pyqt_flags_enums = self._get_plugin_annotation(p, symbol, annotations,
                'PyQtFlagsEnums', 'PyQt5')
        if pyqt_flags_enums is not None:
            klass.pyqt_flags_enums = pyqt_flags_enums
            klass.pyqt_flags = 1

        pyqt_flags = self._get_plugin_annotation(p, symbol, annotations,
                'PyQtFlags', 'PyQt5')
        if pyqt_flags is not None:
            klass.pyqt_flags = pyqt_flags

        klass.pyqt_no_qmetaobject = annotations.get('PyQtNoQMetaObject', False)
        klass.pyqt_interface = annotations.get('PyQtInterface')

        if klass.is_opaque:
            klass.external = annotations.get('External', False)
        else:
            # A default dtor is public.
            if klass.dtor is None:
                klass.dtor = AccessSpecifier.PUBLIC

            klass.no_default_ctors = annotations.get('NoDefaultCtors', False)

            # Provide a default ctor if required.
            if len(klass.ctors) == 0 and not klass.no_default_ctors:
                signature = Signature(result=Argument(ArgumentType.VOID))
                ctor = Constructor(AccessSpecifier.PUBLIC,
                        py_signature=signature, cpp_signature=signature)

                klass.default_ctor = ctor
                klass.ctors.append(ctor)
                klass.can_create = True

            # Determine the default ctor if required.
            if klass.default_ctor is None:
                last_resort = None

                for ctor in klass.ctors:
                    if ctor.access_specifier is not AccessSpecifier.PUBLIC:
                        continue

                    if len(ctor.py_signature.args) == 0 or ctor.py_signature.args[0].default_value is not None:
                        klass.default_ctor = ctor
                        break

                    # The last resort is the first public ctor.
                    if last_resort is None:
                        last_resort = ctor
                else:
                    klass.default_ctor = last_resort

            klass.deprecated = annotations.get('Deprecated', False)

            if klass.convert_to_type_code is not None and annotations.get('AllowNone', False):
                klass.handles_none = True

            if annotations.get('Abstract', False):
                klass.is_abstract = True
                klass.is_incomplete = True
                klass.can_create = False

            klass.delay_dtor = annotations.get('DelayDtor', False)
            if klass.delay_dtor:
                self.module_state.module.has_delayed_dtors = True

            # There are subtle differences between the add and concat methods
            # and the multiply and repeat methods.  The number versions can
            # have their operands swapped and may return NotImplemented.  If
            # the user has used the /Numeric/ annotation or there are other
            # numeric operators then we use add/multiply.  Otherwise, if the
            # user has used the /Sequence/ annotation or there are indexing
            # operators then we use concat/repeat.
            seq_might = seq_not = False

            for md in klass.members:
                if md.py_slot in (PySlot.GETITEM, PySlot.SETITEM, PySlot.DELITEM):
                    # This might be a sequence.
                    seq_might = True

                if md.py_slot in (PySlot.SUB, PySlot.ISUB, PySlot.MOD, PySlot.IMOD, PySlot.FLOORDIV, PySlot.IFLOORDIV, PySlot.TRUEDIV, PySlot.ITRUEDIV, PySlot.POS, PySlot.NEG):
                    # This is definately not a sequence.
                    seq_not = True

            default_to_sequence = (not seq_not and seq_might)

            for md in klass.members:
                # Ignore if it is explicitly numeric.
                if md.is_numeric:
                    continue

                if md.is_sequence or default_to_sequence:
                    if md.py_slot is PySlot.ADD:
                        md.py_slot = PySlot.CONCAT
                    elif md.py_slot is PySlot.IADD:
                        md.py_slot = PySlot.ICONCAT
                    elif md.py_slot is PySlot.MUL:
                        md.py_slot = PySlot.REPEAT
                    elif md.py_slot is PySlot.IMUL:
                        md.py_slot = PySlot.IREPEAT

        if self.in_main_module:
            klass.iface_file.cpp_name.used = True
            klass.py_name.used = True

        self.pop_scope()

        # Check the name in the current scope (ie. the class's parent scope).
        self.check_attributes(p, symbol, py_name, ignore=klass)

        # Check that external classes have only been declared at the global
        # scope.
        if klass.external and self.scope is not None:
            self.parser_error(p, symbol,
                    "/External/ classes/structs/unions can only be declared in the global scope")

        return klass

    def define_class(self, p, symbol, class_key, scoped_name, annotations,
            superclasses=None):
        """ Create a new class and make it the current scope. """

        klass = self.new_class(p, symbol, IfaceFileType.CLASS,
                normalised_scoped_name(scoped_name, self.scope),
                virtual_error_handler=annotations.get('VirtualErrorHandler'),
                type_hints=self.get_type_hints(p, symbol, annotations))

        klass.class_key = class_key
        klass.superclasses = superclasses

        self.push_scope(klass,
                AccessSpecifier.PRIVATE if class_key is ClassKey.CLASS else AccessSpecifier.PUBLIC)

    def disambiguate_token(self, value, keywords):
        """ Disambiguate a token by inspecting its value. """

        # This seems to be needed because it's not possible to get lex() to do
        # it.  The problem seems to be that you can't control the order in
        # which lex() applies its regular expressions despite what the
        # documentation says.  It seems that tokens that are specific to a
        # state are always added after everything else no matter where they
        # appear in the file.

        if value in keywords:
            token_type = value
        elif value == '...':
            token_type = 'ELLIPSIS'
        elif value.startswith('.'):
            token_type = 'FILE_PATH'
        else:
            for marker in ('/', '..', '-'):
                if marker in value:
                    token_type = 'FILE_PATH'
                    break
            else:
                if '.' in value:
                    token_type = 'DOTTED_NAME'
                else:
                    token_type = 'NAME'

        return token_type

    def find_class(self, p, symbol, iface_file_type, fq_cpp_name,
            tmpl_arg=False):
        """ Return a WrappedClass object for a C++ name creating it if
        necessary.
        """

        return self._find_class_with_iface_file(
                self.find_iface_file(p, symbol, fq_cpp_name, iface_file_type),
                tmpl_arg=tmpl_arg)

    def find_exception(self, p, symbol, fq_cpp_name, raise_code=None):
        """ Find an exception, optionally creating a new one. """

        # Note that we don't normalise the name.
        iface_file = self.find_iface_file(p, symbol, fq_cpp_name,
                IfaceFileType.EXCEPTION)

        # See if it is an existing one.
        for w_exception in self.spec.exceptions:
            if w_exception.iface_file is iface_file:
                return w_exception

        # If it is an exception interface file then we have never seen this
        # name before.  We require that exceptions are defined before being
        # used, but don't make the same requirement of classes (for reasons of
        # backwards compatibility).  Therefore the name must be reinterpreted
        # as a (as yet undefined) class.
        if raise_code is not None:
            if iface_file.type is not IfaceFileType.EXCEPTION:
                self.parser_error(p, symbol,
                        "there is already a class with the same name or the exception has been used before being defined")

            class_exception = None
        else:
            # The C++ name must be the name of a class which implements an
            # exception.
            if iface_file.type is IfaceFileType.EXCEPTION:
                iface_file.type = IfaceFileType.CLASS

            class_exception = self._find_class_with_iface_file(iface_file)

        # Create a new one.
        w_exception = WrappedException(iface_file, raise_code,
                class_exception=class_exception)

        self.spec.exceptions.append(w_exception)

        return w_exception

    def new_class(self, p, symbol, iface_file_type, fq_cpp_name,
            virtual_error_handler=None, type_hints=None):
        """ Create a new, unannotated class and add it to the current scope.
        """

        if self.scope_access_specifier is AccessSpecifier.PRIVATE:
            self.parser_error(p, symbol,
                    "classes, structs, unions and namespaces must be in the public or protected sections")

        scope = self.scope

        if scope is not None:
            if self.scope_access_specifier is AccessSpecifier.PROTECTED and not self._protected_is_public:
                scope.is_protected = True

                if iface_file_type is IfaceFileType.CLASS:
                    scope.needs_shadow = True

            # Header code from outer scopes is also included.
            type_header_code = scope.iface_file.type_header_code
        else:
            type_header_code = None

        if bool(self.c_bindings):
            # C structs and unions are always global types no matter where they
            # are defined.
            fq_cpp_name = ScopedName(fq_cpp_name.base_name)
            scope = None

        klass = self.find_class(p, symbol, iface_file_type, fq_cpp_name)

        # Check it hasn't already been defined.
        if iface_file_type is not IfaceFileType.NAMESPACE and klass.iface_file.module is not None:
            self.parser_error(p, symbol,
                    "the class/struct/union has already been defined")

        # Complete the initialisation.
        klass.scope = scope
        klass.iface_file.module = self.module_state.module
        klass.virtual_error_handler = virtual_error_handler
        klass.type_hints = type_hints

        if type_header_code is not None:
            klass.iface_file.type_header_code.extend(type_header_code)

        # See if it is a namespace extender.
        if iface_file_type is IfaceFileType.NAMESPACE:
            for ns in self.spec.classes:
                if ns is klass:
                    continue

                if ns.iface_file.type is not IfaceFileType.NAMESPACE:
                    continue

                if ns.iface_file.fq_cpp_name != klass.iface_file.fq_cpp_name:
                    continue

                klass.real_class = ns

                if self.in_main_module:
                    ns.iface_file.needed = True

        return klass

    def _find_class_with_iface_file(self, iface_file, tmpl_arg=False):
        """ Return a WrappedClass object for an interface file creating it if
        necessary.
        """

        # See if it already exists.
        for klass in self.spec.classes:
            if klass.iface_file is iface_file:
                if not self.parsing_template:
                    try:
                        self._template_arg_classes.remove(klass)
                    except ValueError:
                        pass

                return klass

        # Create a new one.
        klass = WrappedClass(iface_file,
                cached_name(self.spec, iface_file.fq_cpp_name.base_name), None)

        if tmpl_arg:
            self._template_arg_classes.append(klass)

        # Use the same ordering as the old parser.
        self.spec.classes.insert(0, klass)

        return klass

    def add_ctor(self, p, symbol, arg_list, annotations, *, exceptions,
            cpp_signature, docstring, premethod_code, method_code):
        """ Create a Constructor and add it to the current scope. """

        scope = self.scope

        # Error checks.
        if p[symbol] != scope.iface_file.fq_cpp_name.base_name:
            self.parser_error(p, symbol,
                    "constructor doesn't have the same name as its class")

        if bool(self.c_bindings) and len(arg_list) != 0 and method_code is None:
            self.parser_error(p, symbol,
                    "constructors with arguments in C modules must specify %MethodCode")

        if self.scope_pyqt_method_specifier is not None:
            self.parser_error(p, symbol,
                    "constructors must be in the public, protected or private sections")

        # Handle the access specifier.
        access_specifier = self.scope_access_specifier

        if access_specifier is AccessSpecifier.PROTECTED and self._protected_is_public:
            access_specifier = AccessSpecifier.PUBLIC

        # Handle the signatures allowing it to be used like a function
        # signature.
        py_signature = Signature(args=arg_list,
                result=Argument(ArgumentType.VOID))
        self._check_ellipsis(p, symbol, py_signature)

        # Configure the constructor.
        ctor = Constructor(access_specifier, py_signature)

        if annotations.get("NoDerived", False):
            if cpp_signature is not None:
                self.parser_error(p, symbol,
                        "/NoDerived/ may not be specified with an explicit C++ signature")

            if method_code is None:
                self.parser_error(p, symbol,
                        "%MethodCode must also be specified if /NoDerived/ is specified")
        elif cpp_signature is None:
            ctor.cpp_signature = py_signature
        else:
            self._check_ellipsis(p, symbol, cpp_signature)
            ctor.cpp_signature = cpp_signature

        if annotations.get("Default", False):
            if scope.default_ctor is None:
                scope.default_ctor = ctor
            else:
                self.parser_error(p, symbol,
                        "/Default/ has already been specified for another constructor")

        ctor.docstring = docstring
        ctor.gil_action = self._get_gil_action(p, symbol, annotations)
        ctor.deprecated = annotations.get('Deprecated', False)

        if access_specifier is not AccessSpecifier.PRIVATE:
            ctor.kw_args = self._get_kw_args(p, symbol, annotations,
                    py_signature)
            scope.can_create = True

            if access_specifier is AccessSpecifier.PROTECTED:
                scope.needs_shadow = True

        ctor.method_code = method_code
        ctor.no_type_hint = annotations.get('NoTypeHint', False)
        ctor.posthook = annotations.get('PostHook')
        ctor.prehook = annotations.get('PreHook')
        ctor.premethod_code = premethod_code
        ctor.throw_args = exceptions

        if method_code is None and annotations.get('NoRaisesPyException') is None:
            if self.module_state.all_raise_py_exception or annotations.get('RaisesPyException', False):
                ctor.raises_py_exception = True

        ctor.transfer = self.get_transfer(p, symbol, annotations)

        scope.ctors.append(ctor)

    def add_dtor(self, p, symbol, name, annotations, *, exceptions,
            abstract, premethod_code, method_code, virtual_catcher_code):
        """ Add a dtor to the current scope. """

        scope = self.scope

        # Error checks.
        if name != scope.iface_file.fq_cpp_name.base_name:
            self.parser_error(p, symbol,
                    "destructor doesn't have the same name as its class")

        if scope.dtor is not None:
            self.parser_error(p, symbol, "destructor has already been defined")

        if bool(self.c_bindings) and method_code is None:
            self.parser_error(p, symbol,
                    "destructors in C modules must specify %MethodCode")

        if self.scope_pyqt_method_specifier is not None:
            self.parser_error(p, symbol,
                    "destructors must be in the public, protected or private sections")

        if self.parsing_virtual:
            if self.scope.class_key is ClassKey.UNION:
                self.parser_error(p, symbol,
                        "a union cannot have a virtual destructor")

            if virtual_catcher_code is not None:
                self.parser_error(p, symbol,
                        "destructors in C modules cannot specify %VirtualCatcherCode")
        elif abstract:
            self.parser_error(p, symbol,
                    "an abstract destructor must be virtual")

        # Configure the destructor.
        scope.dtor = self.scope_access_specifier
        scope.dtor_gil_action = self._get_gil_action(p, symbol, annotations)
        scope.dtor_throw_args = exceptions
        scope.dtor_virtual_catcher_code = virtual_catcher_code

        if premethod_code is not None:
            scope.dealloc_code.append(premethod_code)

        if method_code is not None:
            scope.dealloc_code.append(method_code)

        if abstract:
            scope.is_abstract = True

        if self.parsing_virtual or len(scope.dealloc_code) != 0:
            scope.needs_shadow = True

    def add_enum(self, p, symbol, cpp_name, is_scoped, annotations, members):
        """ Create a new enum and add it to the current scope. """

        if self.scope_access_specifier is AccessSpecifier.PRIVATE:
            self.parser_error(p, symbol, "class enums cannot be private")

        if is_scoped:
            self.cpp_only(p, symbol, "scoped enums")

        # Determine the base type.
        base_type_s = annotations.get('BaseType')
        base_type = EnumBaseType.ENUM

        if base_type_s is not None:
            if self.spec.abi_version < (13, 0):
                self.parser_error(p, symbol,
                        "/BaseType/ is only supported for ABI v13.0 and later")

            if base_type_s == 'Enum':
                base_type = EnumBaseType.ENUM
            elif base_type_s == 'Flag':
                base_type = EnumBaseType.FLAG
            elif base_type_s == 'IntEnum':
                base_type = EnumBaseType.INT_ENUM
            elif base_type_s == 'UIntEnum':
                base_type = EnumBaseType.UNSIGNED_INT_ENUM
            elif base_type_s == 'IntFlag':
                base_type = EnumBaseType.INT_FLAG
            else:
                self.parser_error(p, 1,
                        "'{0}' is an invalid value of /BaseType/".format(
                                base_type_s))

        if cpp_name is None:
            fq_cpp_name = None
            cached_fq_cpp_name = None
            py_name = None
        else:
            fq_cpp_name = normalised_scoped_name(cpp_name, self.scope)
            cached_fq_cpp_name = cached_name(self.spec, str(fq_cpp_name))
            py_name = cached_name(self.spec,
                    self.get_py_name(cpp_name, annotations))

            self.check_attributes(p, symbol, py_name)

            if self.in_main_module:
                cached_fq_cpp_name.used = True
                py_name.used = True

        w_enum = WrappedEnum(base_type, fq_cpp_name, self.module_state.module,
                cached_fq_cpp_name=cached_fq_cpp_name, is_scoped=is_scoped,
                py_name=py_name, scope=self.scope)

        if self.scope_access_specifier is AccessSpecifier.PROTECTED:
            if not self._protected_is_public:
                w_enum.is_protected = True
                self.scope.needs_shadow = True

        w_enum.no_scope = annotations.get('NoScope', False)
        w_enum.no_type_hint = annotations.get('NoTypeHint', False)

        # Create the members.
        for m_cpp_name, m_py_name, m_no_type_hint in members:
            if not is_scoped:
                self.check_attributes(p, symbol, m_py_name.name)

            w_enum.members.append(
                    WrappedEnumMember(m_cpp_name, m_py_name, w_enum,
                            no_type_hint=m_no_type_hint))

            if self.in_main_module:
                m_py_name.used = True

        self.spec.enums.insert(0, w_enum)

    def add_function(self, p, symbol, cpp_name, result, arg_list, annotations,
            *, const=False, final=False, exceptions=None, abstract=False,
            cpp_signature=None, docstring=None, premethod_code=None,
            method_code=None, virtual_catcher_code=None,
            virtual_call_code=None):
        """ Create and return an Overload and add it to the current scope. """

        # Get the Python name.
        py_name = self.get_py_name(cpp_name, annotations)

        # Handle the signatures.
        py_signature = Signature(args=arg_list, result=result)
        self._check_ellipsis(p, symbol, py_signature)

        if cpp_signature is None:
            cpp_signature = py_signature
        else:
            self._check_ellipsis(p, symbol, cpp_signature)

        # Find (or create) the member shared by overloads with the same Python
        # name.
        member = self._find_member(p, symbol, py_name, arg_list, annotations,
                method_code)

        # Create the overload.
        overload = Overload(self.scope_access_specifier, member, cpp_name,
                cpp_signature, py_signature)

        for m in self.module_state.module.global_functions:
            if m is member:
                self.module_state.module.overloads.append(overload)
                break
        else:
            self.scope.overloads.append(overload)

        overload.pyqt_method_specifier = self.scope_pyqt_method_specifier

        if overload.access_specifier is AccessSpecifier.PROTECTED and self._protected_is_public:
            overload.access_specifier = AccessSpecifier.PUBLIC
            overload.access_is_really_protected = True

        if overload.access_specifier is AccessSpecifier.PROTECTED:
            self.scope.needs_shadow = True
            member.has_protected = True

        if overload.access_specifier is AccessSpecifier.PUBLIC:
            if overload.pyqt_method_specifier is PyQtMethodSpecifier.SIGNAL:
                self.scope.needs_shadow = True

        overload.docstring = docstring
        overload.is_abstract = abstract
        overload.is_const = const
        overload.is_delattr = (py_name == '__delattr__')
        overload.is_final = final
        overload.premethod_code = premethod_code
        overload.throw_args = exceptions
        overload.virtual_call_code = virtual_call_code
        overload.virtual_catcher_code = virtual_catcher_code
        overload.source_location = self.get_source_location(p, symbol)

        # See if the function is a non-lazy method.  These are methods that
        # Python expects to see defined in the type before any instance of the
        # type is created.
        if self.scope is not None:
            NONLAZY_METHOD_NAMES = (
                '__getattribute__',
                '__getattr__',
                '__enter__',
                '__exit__',
                '__aenter__',
                '__aexit__',
            )

            if cpp_name in NONLAZY_METHOD_NAMES:
                self.scope.has_nonlazy_method = True

        # Handle any annotations.
        overload.abort_on_exception = annotations.get('AbortOnException',
                False)

        auto_gen = annotations.get('AutoGen')
        if auto_gen is not None:
            overload.is_auto_generated = self.evaluate_feature_or_platform(p,
                    symbol, name=auto_gen)

        overload.gil_action = self._get_gil_action(p, symbol, annotations)
        overload.factory = annotations.get('Factory', False)
        overload.deprecated = annotations.get('Deprecated', False)
        overload.new_thread = annotations.get('NewThread', False)
        overload.transfer = self.get_transfer(p, symbol, annotations)

        if overload.access_specifier is not AccessSpecifier.PRIVATE:
            if member.py_slot is None or member.py_slot is PySlot.CALL:
                overload.kw_args = self._get_kw_args(p, symbol, annotations,
                        py_signature, need_name=member.has_protected)

                if overload.kw_args is not KwArgs.NONE:
                    member.allow_keyword_args = True

                # If the overload is protected and defined in an imported
                # module then we need to make sure that any other overloads'
                # keyword argument names are marked as used.
                if overload.pyqt_method_specifier is not PyQtMethodSpecifier.SIGNAL and overload.access_specifier is AccessSpecifier.PROTECTED and not self.in_main_module:
                    for kwod in self.scope.overloads:
                        if kwod.common is not member:
                            continue

                        if kwod.kw_args is KwArgs.NONE:
                            continue

                        for arg in kwod.py_signature.args:
                            if kwod.kw_args is KwArgs.OPTIONAL and arg.default_value is None:
                                continue

                            if arg.name is not None:
                                arg.name.used = True

        overload.no_type_hint = annotations.get('NoTypeHint', False)
        overload.posthook = annotations.get('PostHook')
        overload.prehook = annotations.get('PreHook')

        if method_code is None and annotations.get('NoRaisesPyException') is None:
            if self.module_state.all_raise_py_exception or annotations.get('RaisesPyException', False):
                overload.raises_py_exception = True

        overload.virtual_error_handler = annotations.get('VirtualErrorHandler')
        overload.no_virtual_error_handler = annotations.get(
                'NoVirtualErrorHandler', False)

        if annotations.get('Numeric', False):
            if member.is_sequence:
                self.parser_error(p, symbol,
                        "an overload has already specified /Sequence/")
            else:
                member.is_numeric = True

        if annotations.get('Sequence', False):
            if member.is_numeric:
                self.parser_error(p, symbol,
                        "an overload has already specified /Numeric/")
            else:
                member.is_sequence = True

        self.apply_common_argument_annotations(p, symbol,
                overload.py_signature.result, annotations)

        overload.method_code = method_code

        # Add some auto-generated slots if required.
        if '__len__' in annotations:
            len_method_code = method_code
            if len_method_code is None:
                len_method_code = CodeBlock("Auto-generated",
                        text='            sipRes = (Py_ssize_t)sipCpp->{0}();\n'.format(cpp_name))

            len_py_signature = Signature(result=Argument(ArgumentType.SSIZE))

            self._add_auto_slot(p, symbol, annotations, '__len__',
                    len_py_signature, len_py_signature, len_method_code)

        if '__matmul__' in annotations:
            self._add_auto_slot(p, symbol, annotations, '__matmul__',
                    py_signature, cpp_signature, method_code)

        if '__imatmul__' in annotations:
            self._add_auto_slot(p, symbol, annotations, '__imatmul__',
                    py_signature, cpp_signature, method_code)

        return overload

    def add_mapped_type(self, p, symbol, cpp_type, annotations):
        """ Create a new mapped type and add it to the current scope. """

        # Check the type is one we want to map.
        if cpp_type.type is ArgumentType.DEFINED:
            fq_cpp_name = cpp_type.definition = normalised_scoped_name(
                    cpp_type.definition, self.scope)
            cpp_name = fq_cpp_name.base_name
        elif cpp_type.type is ArgumentType.TEMPLATE:
            cpp_type.definition.cpp_name = normalised_scoped_name(
                    cpp_type.definition.cpp_name, self.scope)
            fq_cpp_name = encoded_template_name(cpp_type.definition)
            cpp_name = None
        elif cpp_type.type is ArgumentType.STRUCT:
            fq_cpp_name = cpp_type.definition = normalised_scoped_name(
                    cpp_type.definition, self.scope)
            cpp_name = fq_cpp_name.base_name
        else:
            self.parser_error(p, symbol, "invalid type for %MappedType")
            fq_cpp_name = ''

        iface_file = self.find_iface_file(p, symbol, fq_cpp_name,
                IfaceFileType.MAPPED_TYPE, cpp_type=cpp_type)

        # Check it hasn't already been defined.
        for mtd in self.spec.mapped_types:
            if mtd.iface_file is iface_file:
                # We allow types based on the same template but with different
                # arguments.
                if cpp_type.type is not ArgumentType.TEMPLATE or same_base_type(cpp_type, mtd.type):
                    self.parser_error(p, symbol,
                            "the mapped type has already been defined in this module")

        # The module may not have been set yet.
        iface_file.module = self.module_state.module

        # Create a new mapped type.
        mapped_type = MappedType(iface_file, cpp_type)

        mapped_type.cpp_name = cached_name(self.spec,
                argument_as_str(cpp_type))

        if cpp_name is not None:
            mapped_type.py_name = cached_name(self.spec,
                    annotations.get('PyName', cpp_name))

        self.annotate_mapped_type(p, symbol, mapped_type, annotations)
        self.spec.mapped_types.insert(0, mapped_type)

        if self.in_main_module:
            mapped_type.cpp_name.used = True

            if mapped_type.py_name is not None:
                mapped_type.py_name.used = True

        self.push_scope(mapped_type)

    def add_qualifier(self, p, symbol, name, type, order=0, timeline=0):
        """ Create a Qualifier and add it to the current module. """

        module = self.module_state.module

        # See if it already exists.
        qualifier = self.find_qualifier(p, symbol, name, required=False)

        if qualifier is not None:
            # We allow versions to be defined more than once so long as they
            # are in different timelines.  It is sometimes necessary to define
            # the same timeline in multiple modules if a module that others
            # depend on is added during the timeline.
            if qualifier.type is not QualifierType.TIME or type is not QualifierType.TIME or (qualifier.module is module and qualifier.timeline == timeline):
                self.parser_error(p, symbol,
                        "'{0}' has already been defined as a qualifier".format(
                                name))

                return

        qualifier = Qualifier(module, name, type, order=order,
                timeline=timeline)

        if type is QualifierType.TIME or not self.skipping:
            qualifier.enabled_by_default = True

        module.qualifiers.insert(0, qualifier)

    def add_typedef(self, p, symbol, typedef):
        """ Add a typedef to the current scope. """

        if self.spec.is_strict:
            for td in self.spec.typedefs:
                if td.fq_cpp_name == typedef.fq_cpp_name:
                    self.parser_error(p, symbol,
                            "'{0}' has already been defined".format(
                                    typedef.fq_cpp_name))

        self.module_state.module.nr_typedefs += 1

        self.spec.typedefs.append(typedef)

    def annotate_mapped_type(self, p, symbol, mapped_type, annotations):
        """ Apply annotations to a mapped type. """

        mapped_type.handles_none = annotations.get('AllowNone', False)
        mapped_type.no_assignment_operator = annotations.get(
                'NoAssignmentOperator', False)
        mapped_type.no_copy_ctor = annotations.get('NoCopyCtor', False)
        mapped_type.no_default_ctor = annotations.get('NoDefaultCtor', False)
        mapped_type.no_release = annotations.get('NoRelease', False)
        mapped_type.type_hints = self.get_type_hints(p, symbol, annotations)

        if mapped_type.no_release:
            mapped_type.no_assignment_operator = True
            mapped_type.no_copy_ctor = True
            mapped_type.no_default_ctor = True

        pyqt_flags = self._get_plugin_annotation(p, symbol, annotations,
                'PyQtFlags', 'PyQt6')
        if pyqt_flags is not None:
            mapped_type.pyqt_flags = pyqt_flags

    def apply_common_argument_annotations(self, p, symbol, arg, annotations):
        """ Apply the annotations common to callable arguments and return type.
        """

        arg.allow_none = annotations.get('AllowNone', False)
        arg.disallow_none = annotations.get('DisallowNone', False)
        arg.no_copy = annotations.get('NoCopy', False)

        # We need to use an exception because we have to distinguish between
        # a missing annotation and one without a value specified.
        try:
            key = annotations['KeepReference']

            if key is None:
                key = self.module_state.module.next_key
                self.module_state.module.next_key -= 1
            elif key < 0:
                self.parser_error(p, symbol,
                        "a /KeepReference/ key cannot be negative")

            arg.key = key
        except KeyError:
            pass

    def apply_type_annotations(self, p, symbol, type, annotations):
        """ Apply the annotations for an argument type. """

        # The type annotations.
        type.type_hints = self.get_type_hints(p, symbol, annotations)

        # The PyInt annotation.
        if annotations.get('PyInt') is not None:
            if type.type is ArgumentType.STRING:
                type.type = ArgumentType.BYTE
            elif type.type is ArgumentType.SSTRING:
                type.type = ArgumentType.SBYTE
            elif type.type is ArgumentType.USTRING:
                type.type = ArgumentType.UBYTE

        # The Encoding annotation.
        can_be_encoded = type.type is ArgumentType.STRING and type.array is ArrayArgument.NONE and not type.is_reference

        if can_be_encoded:
            encoding = annotations.get('Encoding')
            if encoding is None:
                default_encoding = self.module_state.default_encoding
                if default_encoding is not None:
                    type.type = default_encoding
            else:
                type.type = self.convert_encoding(p, symbol, encoding)

    def check_annotations(self, p, symbol, context, annotations):
        """ Check that all the annotations provided as a dict of name/values
        are valid in a given context.
        """

        for name in p[symbol]:
            if name not in annotations:
                self.parser_error(p, symbol,
                        "{0} is not a valid {1} annotation".format(name,
                                context))

    def check_attributes(self, p, symbol, py_name, is_function=False,
            ignore=None):
        """ Check that a Python name will not clash with another object in the
        same Python scope.
        """

        # We don't do any check for a non-strict parse.
        if not self.spec.is_strict:
            return

        # Report a name clash with something.
        def clash(thing):
            self.parser_error(p, symbol,
                    "there is already {0} in scope called '{1}'".format(thing,
                            py_name))

        # Check the enums.
        for ed in self.spec.enums:
            if ed.py_name is None:
                continue

            if ed.scope is not self.scope:
                continue

            if ed.py_name.name == py_name:
                clash("an enum")

            if not ed.is_scoped:
                for emd in ed.members:
                    if emd.py_name.name == py_name:
                        clash("an enum member")

        # Only check the members if this attribute isn't a member because we
        # can handle members with the same name in the same scope.
        if not is_function:
            if self.scope is None:
                members = self.module_state.module.global_functions
                thing = "a function"
            else:
                members = self.scope.members
                thing = "a method"

            for md in members:
                if md.py_name.name == py_name:
                    clash(thing)
                    break

        # There is nothing more to check in mapped types.
        if isinstance(self.scope, MappedType):
            return

        # Check the variables.
        for vd in self.spec.variables:
            if vd.scope is not self.scope:
                continue

            if vd.py_name.name == py_name:
                clash("a variable")
                break

        # Check the classes.
        for cd in self.spec.classes:
            if cd.scope is not self.scope:
                continue

            # A class will have already been added to the scope and this will
            # tell us to ignore it.
            if cd is ignore:
                continue

            if cd.external:
                continue

            if cd.py_name.name == py_name:
                clash("a class or namespace")

        if self.scope is None:
            # Check the exceptions.
            for xd in self.spec.exceptions:
                if xd.py_name is not None and xd.py_name == py_name:
                    clash("an exception")
                    break
        else:
            # Check the properties.
            for pd in self.scope.properties:
                if pd.name.name == py_name:
                    clash("a property")
                    break

    def cpp_only(self, p, symbol, feature):
        """ Check that a C++ feature isn't being used in a C module. """

        if bool(self.c_bindings):
            self.parser_error(p, symbol,
                    feature + " are not allowed in a C module")

    def convert_docstring_format(self, p, symbol):
        """ Convert a string to the corresponding DocstringFormat member. """

        value = p[symbol]

        if value == 'deindented':
            return DocstringFormat.DEINDENTED

        if value == 'raw':
            return DocstringFormat.RAW

        self.parser_error(p, symbol,
                "'{0}' is not a valid value for a docstring format".format(
                        value))

        # Return any value of the right type.
        return DocstringFormat.RAW

    def convert_docstring_signature(self, p, symbol):
        """ Convert a string to the corresponding DocstringSignature member.
        """

        value = p[symbol]

        if value == 'appended':
            return DocstringSignature.APPENDED

        if value == 'discarded':
            return DocstringSignature.DISCARDED

        if value == 'prepended':
            return DocstringSignature.PREPENDED

        self.parser_error(p, symbol,
                "'{0}' is not a valid value for a docstring signature".format(
                        value))

        # Return any value of the right type.
        return DocstringSignature.PREPENDED

    def convert_encoding(self, p, symbol, value=None):
        """ Convert a string to the corresponding ArgumentType member. """

        if value is None:
            value = p[symbol]

        if value == 'ASCII':
            return ArgumentType.ASCII_STRING

        if value == 'Latin-1':
            return ArgumentType.LATIN1_STRING

        if value == 'None':
            return ArgumentType.STRING

        if value == 'UTF-8':
            return ArgumentType.UTF8_STRING

        self.parser_error(p, symbol,
                "'{0}' is not a valid encoding".format(value))

        # Return any value of the right type.
        return ArgumentType.NONE

    def convert_kw_args(self, p, symbol, value=None):
        """ Convert a string to the corresponding KwArgs member. """

        if value is None:
            value = p[symbol]

        if value == 'All':
            return KwArgs.ALL

        if value == 'None':
            return KwArgs.NONE

        if value == 'Optional':
            return KwArgs.OPTIONAL

        self.parser_error(p, symbol,
                "'{0}' is not a valid value for keyword arguments".format(
                        value))

        # Return any value of the right type.
        return KwArgs.OPTIONAL

    def ensure_import(self):
        """ We allow %Modules that are part of a %CompositeModule to be either
        %Imported or %Included.  In the case of the latter we need to adjust
        things so that it appears like the former.
        """

        sip_file, raw_sip_file, input, lineno, lexpos, module_state = self._file_stack.pop()

        if module_state is None:
            self._import_module(self._sip_file)
            module_state = self.module_state

        self._file_stack.append(
                (sip_file, raw_sip_file, input, lineno, lexpos, module_state))

    def evaluate_feature_or_platform(self, p, symbol, name=None,
            inverted=False):
        """ Evaluate a feature or platform qualifier. """

        if name is None:
            name = p[symbol]

        qual = self.find_qualifier(p, symbol, name)
        if qual is None:
            return False

        if qual.type is QualifierType.FEATURE:
            value = qual.name not in self._disabled_features
        elif qual.type is QualifierType.PLATFORM:
            # The platform is always ignored in non-strict mode.
            if not self.spec.is_strict:
                return True

            value = qual.name in self.tags
        else:
            self.parser_error(p, symbol,
                    "'{0}' is a %Timeline qualifier which can only be used in a range".format(name))
            return False

        if inverted:
            value = not value

        return value

    def evaluate_timeline(self, p, symbol_lower, symbol_upper):
        """ Evaluate a timeline qualifier. """

        # Get the lower and upper qualifiers if specified.
        lower_name = p[symbol_lower]
        if lower_name:
            lower_qual = self._find_timeline_qualifier(p, symbol_lower)
            if lower_qual is None:
                return False
        else:
            lower_qual = None

        upper_name = p[symbol_upper]
        if upper_name:
            upper_qual = self._find_timeline_qualifier(p, symbol_upper)
            if upper_qual is None:
                return False
        else:
            upper_qual = None

        if lower_qual is None and upper_qual is None:
            self.parser_error(p, symbol_lower,
                    "the lower and upper bounds cannot both be omitted")
            return False

        if lower_qual is not None and upper_qual is not None:
            if lower_qual.module is not upper_qual.module or lower_qual.timeline != upper_qual.timeline:
                self.parser_error(p, symbol_lower,
                        "'{0}' and '{1}' are defined in different %Timelines".format(lower_name, upper_name))
                return False

            if lower_qual is upper_qual:
                self.parser_error(p, symbol_lower,
                        "the lower and upper bounds must be different")
                return False

            if lower_qual.order > upper_qual.order:
                self.parser_error(p, symbol_lower,
                        "'{0}' is later in the %Timeline than '{1}'".format(lower_name, upper_name))
                return False

        # Get the module and timeline.
        if lower_qual is not None:
            module = lower_qual.module
            timeline = lower_qual.timeline
        else:
            module = upper_qual.module
            timeline = upper_qual.timeline

        # Handle the SIP version number pseudo-timeline.
        if timeline < 0:
            if lower_qual is not None and self._hex_version < lower_qual.order:
                return False

            if upper_qual is not None and self._hex_version >= upper_qual.order:
                return False

            return True

        # See if there is a selected qualifier withing range.
        for qual in module.qualifiers:
            if qual.type is QualifierType.TIME and qual.timeline == timeline and qual.name in self.tags:
                if lower_qual is not None and qual.order < lower_qual.order:
                    return False
                if upper_qual is not None and qual.order >= upper_qual.order:

                    return False

                return True

        return upper_qual is None

    def find_qualifier(self, p, symbol, name, required=True):
        """ Return a Qualifier or None if one doesn't exist. """

        for module in self.modules:
            for qual in module.qualifiers:
                if qual.name == name:
                    return qual

        # Qualifiers corresponding to the SIP version are created on the fly.
        if name.startswith('SIP_'):
            parts = name.split('_')[1:]
            if len(parts) > 3:
                order = -1
            else:
                while len(parts) < 3:
                    parts.append('0')

                order = 0

                for part in parts:
                    try:
                        order = (order << 8) + int(part)
                    except ValueError:
                        order = -1
                        break

            if order >= 0:
                module = self.module_state.module

                qualifier = Qualifier(module, name, QualifierType.TIME,
                        order=order)
                module.qualifiers.append(qualifier)

                return qualifier

        if required:
            self.parser_error(p, symbol,
                    "'{0}' is not a known qualifier".format(name))

        return None

    def get_py_name(self, cpp_name, annotations):
        """ Return a valid Python name given a C/C++ name. """

        # Use any name specified by annotation.
        try:
            return annotations['PyName']
        except KeyError:
            pass

        # Use the C/C++ name.
        py_name = cpp_name

        # Apply any automatic naming rules.
        for rule in self.module_state.auto_py_name_rules:
            if rule[0] == 'REMOVE_LEADING':
                leading = rule[1]
                if py_name.startswith(leading):
                    py_name = py_name[len(leading):]

        # Fix any Python keywords.
        if py_name in self._PYTHON_KEYWORDS:
            py_name += '_'

        return py_name

    def get_transfer(self, p, symbol, annotations):
        """ Return the a Transfer value from a dict of annotations. """

        has_transfer = annotations.get('Transfer', False)
        has_transfer_back = annotations.get('TransferBack', False)
        has_transfer_this = annotations.get('TransferThis', False)

        transfer = None

        if has_transfer:
            if has_transfer_back or has_transfer_this:
                pass
            else:
                transfer = Transfer.TRANSFER
        elif has_transfer_back:
            if has_transfer or has_transfer_this:
                pass
            else:
                transfer = Transfer.TRANSFER_BACK
        elif has_transfer_this:
            if has_transfer or has_transfer_back:
                pass
            else:
                transfer = Transfer.TRANSFER_THIS
        else:
            transfer = Transfer.NONE

        if transfer is None:
            self.parser_error(p, symbol,
                    "Only one transfer annotation may be specified")
            transfer = Transfer.NONE

        return transfer

    def get_type_hints(self, p, symbol, annotations):
        """ Return a TypeHints object constructed from a dict of annotations or
        None if none were specified.
        """

        th = annotations.get('TypeHint')
        th_in = annotations.get('TypeHintIn')
        th_out = annotations.get('TypeHintOut')
        th_value = annotations.get('TypeHintValue')

        if th_in is None:
            th_in = th
        elif th is not None:
            self.parser_error(p, symbol,
                    "'TypeHint' and 'TypeHintIn' cannot both be specified")

            return None

        if th_out is None:
            th_out = th
        elif th is not None:
            self.parser_error(p, symbol,
                    "'TypeHint' and 'TypeHintOut' cannot both be specified")

            return None

        if th_in is not None or th_out is not None or th_value is not None:
            # Check that type hints haven't been suppressed.
            if annotations.get('NoTypeHint') is not None:
                self.parser_error(p, symbol,
                        "'NoTypeHint' cannot be specified with a type hint")

            return TypeHints(th_in, th_out, th_value)

        return None

    @property
    def in_main_module(self):
        """ Set if the current module is the main one, ie. the one for which
        code will be generated for.
        """

        return self.module_state.module is self.spec.module

    def instantiate_class_template(self, p, symbol, fq_cpp_name, template,
            py_name, no_type_name, docstring):
        """ Try and instantiate a class template and return True if one was
        found.
        """

        # Look for an appropriate class template.
        for tmpl_names, proto_class in self.class_templates:
            if proto_class.iface_file.fq_cpp_name.matches(template.cpp_name, scope=self.scope) and same_template_signature(tmpl_names, template.types):
                break
        else:
            # There was no class template to instantiate.
            return False

        instantiate_class(p, symbol, fq_cpp_name, tmpl_names, proto_class,
                template, py_name, no_type_name, docstring, self)

        return True

    def lexer_error(self, t, text):
        """ Record an error caused by a token. """

        self._error_log.log(text,
                SourceLocation(self._sip_file, line=t.lineno,
                        column=self._get_column_from_lexpos(t.lexpos)))

    def parse(self, sip_file):
        """ Parse a .sip file and return a 3-tuple of a Specification object, a
        list of Module objects and a list of the .sip files that specify the
        module to be generated.  A UserException is raised if there was an
        error.
        """

        # Note that the retention of the 'raw' filename, ie. that which was
        # specified by the user is only done so that generated '#line'
        # directives match those from older versions of SIP.

        raw_sip_file = sip_file
        sip_file = os.path.abspath(sip_file)

        self.module_state = ModuleState(self.spec.module, sip_file)

        try:
            self._parser.parse(self._read(sip_file, raw_sip_file),
                    lexer=self._lexer, tracking=True)
        except UnexpectedEOF:
            self._unexpected_eof_error()

        self._handle_eom()
        self._error_log.as_exception()

        self.spec.c_bindings = bool(self.c_bindings)

        self.spec.typedefs.sort(key=lambda k: k.fq_cpp_name)
        self.spec.variables.sort(key=lambda k: k.py_name.name)

        # Remove all template classes and anything they contain.
        template_classes = [k for _, k in self.class_templates]

        for enum in list(self.spec.enums):
            if enum.scope in template_classes:
                spec.spec.enums.remove(enum)

        for typedef in list(self.spec.typedefs):
            if typedef.scope in template_classes:
                spec.spec.typedefs.remove(typedef)

        for variable in list(self.spec.variables):
            if variable.scope in template_classes:
                spec.spec.variables.remove(variable)

        for klass in template_classes:
            self.spec.classes.remove(klass)
            self.spec.iface_files.remove(klass.iface_file)

        # Remove all classes that are only template arguments.
        for klass in self._template_arg_classes:
            self.spec.classes.remove(klass)

        return self.spec, self.modules, self._sip_files

    def parser_error(self, p, symbol, text):
        """ Record an error caused by a symbol in a production. """

        self._error_log.log(text, self.get_source_location(p, symbol))

    def pop_file(self):
        """ Restore the current .sip file from the stack and make it current.
        An IndexError is raised if the stack is empty.
        """

        # Restore the state of the previous .sip file.  Note that we don't
        # restore the module state until after the EOF has been seen.
        self._sip_file, self.raw_sip_file, self._input, self._lexer.lineno, lexpos, self._pending_module_state = self._file_stack.pop()

        self._lexer.input(self._input)
        self._lexer.lexpos = lexpos

    def pop_module_state(self):
        """ Restore the current module state. """

        if self._pending_module_state is None:
            return

        self._handle_eom()

        # Inherit any default encoding.
        if self._pending_module_state.default_encoding is None:
            self._pending_module_state.default_encoding = self.module_state.default_encoding

        # Inherit any call_super_init.
        if self._pending_module_state.call_super_init is None:
            self._pending_module_state.call_super_init = self.module_state.call_super_init

        self.module_state = self._pending_module_state
        self._pending_module_state = None

    def pop_scope(self):
        """ Pop the current scope. """

        self._scope_stack.pop()

    def push_file(self, p, symbol, sip_file=None, new_module=False,
            optional=False):
        """ Push the current .sip file onto the stack and make the new one
        current.  The new .sip file may be part of a new module (ie. %Import
        rather than %Include).
        """

        if sip_file is None:
            sip_file = p[symbol]

        raw_sip_file = sip_file

        # Make the name platform-native.
        sip_file = sip_file.replace('/', os.sep)

        # See if the file can be found.
        if os.path.isfile(sip_file):
            pass
        else:
            found = None

            # If the name is relative then check the directory containing the
            # current file and any include directories.
            if not os.path.isabs(sip_file):
                inc_dirs = [os.path.dirname(self._sip_file)]
                inc_dirs.extend(self._include_dirs)

                for inc_dir in inc_dirs:
                    fn = os.path.join(inc_dir, sip_file)
                    if os.path.isfile(fn):
                        found = fn
                        break

            if found is None:
                if not optional:
                    self.parser_error(p, symbol,
                            "'{0}' could not be found".format(raw_sip_file))

                return

            # For historic reasons we keep the absolute name for the raw name.
            raw_sip_file = sip_file = found

        sip_file = os.path.abspath(sip_file)

        # Check we aren't reading the file recursively.
        for detail in self._file_stack:
            if detail[0] == sip_file:
                self.parser_error(p, symbol,
                        "'{0}' is being read recursively".format(sip_file))
                return

        # Ignore the file if we have already read it.
        if sip_file in self._all_sip_files:
            return

        if new_module:
            old_module_state = self.module_state
            self._import_module(sip_file)
        else:
            # This means that the file was %Included rather than %Imported.
            old_module_state = None

        # Save the state of the current .sip file.
        self._file_stack.append(
                (self._sip_file, self.raw_sip_file, self._input,
                        self._lexer.lineno, self._lexer.lexpos,
                        old_module_state))

        # Make the new one current and give it's content to the lexer.
        self._lexer.lineno = 1
        self._lexer.input(self._read(sip_file, raw_sip_file))

    def push_scope(self, scope, access_specifier=None):
        """ Push a new scope. """

        self._scope_stack.append(ScopeState(scope, access_specifier))

    def set_lexer_state(self, state='INITIAL'):
        """ Set the lexer state. """

        self._lexer.begin(state)

    @property
    def scope(self):
        """ The current scope if any. """

        try:
            return self._scope_stack[-1].scope
        except IndexError:
            pass

        return None

    @property
    def scope_access_specifier(self):
        """ The current access specifier. """

        return None if len(self._scope_stack) == 0 else self._scope_stack[-1].access_specifier

    @scope_access_specifier.setter
    def scope_access_specifier(self, access_specifier):
        """ Set the current access specifier. """

        self._scope_stack[-1].access_specifier = access_specifier

    @property
    def scope_pyqt_method_specifier(self):
        """ The current method specifier. """

        return None if len(self._scope_stack) == 0 else self._scope_stack[-1].pyqt_method_specifier

    @scope_pyqt_method_specifier.setter
    def scope_pyqt_method_specifier(self, pyqt_method_specifier):
        """ Set the current method specifier. """

        self._scope_stack[-1].pyqt_method_specifier = pyqt_method_specifier

    @property
    def skipping(self):
        """ True if symbols are currently being skipped. """

        return self.skip_stack[-1]

    def validate_annotation(self, p, symbol, value):
        """ Validate an annotation and its value and return a valid version of
        the value.
        """

        try:
            value = validate_annotation_value(self, p, symbol, p[symbol],
                    value)
        except InvalidAnnotation as e:
            self.parser_error(p, symbol, str(e))
            value = e.use

        return value

    def validate_function(self, p, symbol, overload):
        """ Validate a completed function. """

        # Shortcuts.
        cpp_only = partial(self.cpp_only, p, symbol)
        error = partial(self.parser_error, p, symbol)

        if overload.access_specifier is None:
            if overload.is_abstract:
                error("abstract non-member functions cannot be specified")
        else:
            cpp_only("struct/union member functions")

        if overload.new_thread:
            # This is an arbitary limitation to make the code generator
            # slightly easier - laziness on my part.
            result = overload.cpp_signature.result

            if result.type is not ArgumentType.VOID or len(result.derefs) != 0:
                error("/NewThread/ may only be specified for void functions")

        if overload.is_static:
            cpp_only("static struct/union data members")

            if overload.pyqt_method_specifier is PyQtMethodSpecifier.SIGNAL:
                error("signals cannot be static")

        if overload.throw_args is not None:
            cpp_only("exceptions")

        if overload.is_virtual:
            if overload.virtual_error_handler is not None and overload.no_virtual_error_handler:
                error("/VirtualErrorHandler/ and /NoVirtualErrorHandler/ cannot both be specified")
        else:
            if overload.abort_on_exception:
                error("/AbortOnException/ cannot be specified for non-virtual member functions")

            if overload.new_thread:
                error("/NewThread/ cannot be specified for non-virtual member functions")

            if overload.virtual_call_code is not None:
                error("%VirtualCallCode cannot be specified for non-virtual member functions")

            if overload.virtual_catcher_code is not None:
                error("%VirtualCatcherCode cannot be specified for non-virtual member functions")

            if overload.virtual_error_handler is not None:
                error("/VirtualErrorHandler/ cannot be specified for non-virtual member functions")

            if overload.no_virtual_error_handler:
                error("/NoVirtualErrorHandler/ cannot be specified for non-virtual member functions")

        if overload.factory:
            if overload.transfer is Transfer.TRANSFER_BACK:
                error("/TransferBack/ and /Factory/ cannot both be specified")
        elif self.scope_access_specifier is None or overload.is_static:
            for arg in overload.py_signature.args:
                if arg.transfer is Transfer.TRANSFER_THIS:
                    error("/TransferThis/ may only be specified for constructors and member functions")
                    break

        if overload.common.no_arg_parser and overload.method_code is None:
            error("%MethodCode must be specified when /NoArgParser/ is specified")

    def validate_mapped_type(self, p, symbol, mapped_type):
        """ Validate a completed mapped type. """

        if self.spec.abi_version >= (13, 0):
            convert_to_us = mapped_type.convert_to_type_code is not None and 'sipUserState' in mapped_type.convert_to_type_code.text

            release_us = mapped_type.release_code is not None and 'sipUserState' in mapped_type.release_code.text

            if convert_to_us != release_us:
                self.parser_error(p, symbol,
                        "both %ConvertToTypeCode and %ReleaseCode must use user state or neither must")

            mapped_type.needs_user_state = convert_to_us or release_us
        else:
            if mapped_type.convert_to_type_code is None:
                self.parser_error(p, symbol,
                        "%MappedType must have a %ConvertToTypeCode directive")

            if mapped_type.convert_from_type_code is None:
                self.parser_error(p, symbol,
                        "%MappedType must have a %ConvertFromTypeCode directive")

    def validate_variable(self, p, symbol, variable):
        """ Validate a completed variable. """

        if variable.type.type is ArgumentType.CAPSULE:
            self.parser_error(p, symbol,
                    "capsule variables are not yet supported")

        access_specifier = self.scope_access_specifier
        if access_specifier is not None:
            if access_specifier is not AccessSpecifier.PUBLIC:
                self.parser_error(p, symbol,
                        "class variables must be in the public section")

            if variable.is_static:
                self.cpp_only(p, symbol, "static members in C structures")
            elif variable.access_code is not None:
                self.parser_error(p, symbol,
                        "%AccessCode cannot be specified for non-static class variables")

        if variable.get_code is not None or variable.set_code is not None:
            if variable.access_code is not None:
                self.parser_error(p, symbol,
                        "%AccessCode cannot be specified with %GetCode or %SetCode")

            if self.scope is None:
                # TODO: this can be supported for versions of Python that
                # support module descriptors.
                self.parser_error(p, symbol,
                        "%GetCode or %SetCode cannot be specified for global variables")

        if self.scope is not None and self.scope.iface_file.type is IfaceFileType.NAMESPACE:
            variable.is_static = True

        self.check_attributes(p, symbol, variable.py_name)

    def _add_auto_slot(self, p, symbol, annotations, py_name, py_signature,
            cpp_signature, method_code):
        """ Add an automatically generated slot. """

        member = self._find_member(p, symbol, py_name, py_signature.args,
                annotations, method_code)

        overload = Overload(AccessSpecifier.PUBLIC, member, py_name,
                cpp_signature, py_signature, method_code=method_code)

        if self.scope is None:
            self.module_state.module.overloads.append(overload)
        else:
            self.scope.overloads.append(overload)

    def _check_ellipsis(self, p, symbol, signature):
        """ Check any ellipsis in a signature. """

        seen_ellipsis = False

        for arg in signature.args:
            if arg.type is ArgumentType.ELLIPSIS:
                # Give the argument the standard name so that it appears in
                # docstrings and type hints.  It is never used in generated
                # code so there is no need to add it to the module cache.
                arg.name = CachedName('*args')

                if seen_ellipsis:
                    self.parser_error(p, symbol,
                            "'...' may only be specified once")
                    break

                seen_ellipsis = True
            elif seen_ellipsis and arg.default_value is None:
                self.parser_error(p, symbol,
                        "'...' must be at the end of the positional arguments")
                break

    def find_iface_file(self, p, symbol, fq_cpp_name, iface_file_type,
            cpp_type=None):
        """ Return an interface file for a fully qualified C/C++ name and type
        creating it if necessary.
        """

        def error_logger(text):
            self.parser_error(p, symbol, text)

        return find_iface_file(self.spec, self.module_state.module,
                fq_cpp_name, iface_file_type, error_logger, cpp_type=cpp_type,
                scope=self.scope)

    def _find_member(self, p, symbol, py_name, args, annotations, method_code):
        """ Return (creating if necessary) the member for a Python name. """

        # See if it is a Python slot rather than an ordinary member function
        # and check its requirements.
        slot_detail = slot_name_detail_map.get(py_name)

        if slot_detail is None:
            py_slot = None
        else:
            py_slot, needs_method_code, nr_args_needed = slot_detail

            if needs_method_code and method_code is None:
                self.parser_error("'{0}' requires %MethodCode".format(py_name))

            if nr_args_needed >= 0:
                # Global operators need an extra argument.
                if self.scope_access_specifier is None:
                    nr_args_needed += 1

                    # Global operators can only be numeric or comparision
                    # slots.
                    if invalid_global_slot(py_slot):
                        self.parser_error(p, symbol,
                                "'{0}' cannot be specified at the module level".format(py_name))

                nr_args = len(args)

                if nr_args_needed != nr_args:
                    self.parser_error(p, symbol,
                            "{0} arguments are needed but {1} are provided".format(
                                    nr_args_needed, nr_args))

            # __delattr__ is implemented as __setattr__.
            if py_slot is PySlot.DELATTR:
                if self.in_main_module:
                    cached_name(self.spec, py_name).used = True

                py_slot = PySlot.SETATTR
                py_name = '__setattr__'

        # Check for name clashes.
        self.check_attributes(p, symbol, py_name, is_function=True)

        # Create a new member if necessary.
        no_arg_parser = annotations.get('NoArgParser', False)

        if self.scope is None:
            members = self.module_state.module.global_functions
            namespace_iface_file = None
        elif self.scope.iface_file.type is IfaceFileType.NAMESPACE and py_slot is not None:
            # The scope is a namespace and the function is an operator so
            # handle it as a global, but remember it's C++ scope.
            members = self.module_state.module.global_functions
            namespace_iface_file = self.scope.iface_file
        else:
            members = self.scope.members
            namespace_iface_file = None

        for member in members:
            if member.py_name.name == py_name:
                # If /NoArgParser/ has been specified then there should only be
                # one overload.
                if member.no_arg_parser:
                    self.parser_error(p, symbol,
                            "an overloaded member function has already been specified with /NoArgParser/")

                break
        else:
            member = Member(self.module_state.module,
                    cached_name(self.spec, py_name))
            member.namespace_iface_file = namespace_iface_file
            member.no_arg_parser = no_arg_parser
            member.py_slot = py_slot

            if self.in_main_module:
                member.py_name.used = True

            members.insert(0, member)

        return member

    def _find_timeline_qualifier(self, p, symbol):
        """ Return an optional timeline qualifier. """

        name = p[symbol]

        qual = self.find_qualifier(p, symbol, name)
        if qual is None:
            return None

        if qual.type is not QualifierType.TIME:
            self.parser_error(p, symbol,
                    "'{0}' is not a %Timeline qualifier".format(name))
            return None

        return qual

    def _get_column_from_lexpos(self, lexpos):
        """ Return the column within the current line corresponding to a lexer
        position.
        """

        line_start = self._input.rfind('\n', 0, lexpos) + 1

        return lexpos - line_start + 1

    def _get_gil_action(self, p, symbol, annotations):
        """ Return an appropriate GILAction according to the annotations. """

        hold = annotations.get('HoldGIL', False)
        release = annotations.get('ReleaseGIL', False)

        if hold:
            if release:
                self.parser_error(p, symbol,
                        "/HoldGIL/ and /ReleaseGIL' cannot both be specified")

            return GILAction.HOLD

        if release:
            return GILAction.RELEASE

        return GILAction.DEFAULT

    # The Python keywords.
    _PYTHON_KEYWORDS = (
        'False', 'None', 'True', 'and', 'as', 'assert', 'async', 'await',
        'break', 'class', 'continue', 'def', 'del', 'elif', 'else', 'except',
        'finally', 'for', 'from', 'global', 'if', 'import', 'in', 'is',
        'lambda', 'nonlocal', 'not', 'or', 'pass', 'raise', 'return', 'try',
        'while', 'with', 'yield',
    )

    def _get_kw_args(self, p, symbol, annotations, signature, need_name=False):
        """ Return the keyword argument support. """

        kw_args = annotations.get('KeywordArgs')
        if kw_args is not None:
            kw_args = self.convert_kw_args(p, symbol, kw_args)
        else:
            kw_args = self.module_state.kw_args

        # An ellipsis cannot be used with keyword arguments.
        if len(signature.args) > 0 and signature.args[-1].type is ArgumentType.ELLIPSIS:
            kw_args = KwArgs.NONE

        # Check that there is at least one optional argument.
        if kw_args is not KwArgs.NONE:
            no_name = True

            for arg in signature.args:
                if kw_args is KwArgs.OPTIONAL and arg.default_value is None:
                    continue

                if arg.name is not None:
                    if need_name or self.in_main_module:
                        arg.name.used = True

                    no_name = False

            if no_name:
                kw_args = KwArgs.NONE

        return kw_args

    def _get_plugin_annotation(self, p, symbol, annotations, name, plugin):
        """ Return an annotation that is only supported by a plugin. """

        anno = annotations.get(name)

        if anno is not None and plugin not in self.spec.plugins:
            self.parser_error(p, symbol,
                    "/{0}/ is only supported for {1}".format(name, plugin))

        return anno

    def get_source_location(self, p, symbol):
        """ Return a SourceLocation object for a symbol. """

        return SourceLocation(self._sip_file, line=p.lineno(symbol),
                column=self._get_column_from_lexpos(p.lexpos(symbol)))

    def _handle_eom(self):
        """ Check that the current module is complete. """

        module_state = self.module_state
        module = module_state.module

        if module.fq_py_name is None:
            self._error_log.log("%Module has not been specified",
                    SourceLocation(self._sip_file))

        # call_super_init defaults to False if it wasn't specified.
        module.call_super_init = bool(module_state.call_super_init)

    def _import_module(self, sip_file):
        """ Create a new Module object and corresponding ModuleState object for
        a .sip file and make it current.
        """

        importing_from = self.module_state.module

        module = Module()
        self.modules.append(module)

        module.default_exception = self.module_state.module.default_exception
        self.module_state = ModuleState(module, sip_file)

        # Get the configuration of the new module.
        mod_tags, mod_disabled = get_bindings_configuration(
                self.spec.abi_version[0], sip_file, self._include_dirs)

        for tag in mod_tags:
            if tag not in self.tags:
                self.tags.append(tag)

        for feature in mod_disabled:
            if feature not in self._disabled_features:
                self._disabled_features.append(feature)

        # Add the new import.
        importing_from.imports.append(module)

    def _read(self, sip_file, raw_sip_file):
        """ Return the contents of the current .sip file. """

        try:
            with open(sip_file, encoding=self._encoding) as f:
                self._input = f.read()
        except FileNotFoundError:
            raise UserException("unable to read '{0}'".format(sip_file))
        except UnicodeDecodeError as e:
            raise UserException(
                    "'{0}' doesn't appear to use the '{1}' encoding".format(
                            sip_file, self._encoding),
                    detail=str(e))

        self.raw_sip_file = raw_sip_file
        self._sip_file = sip_file
        self._all_sip_files.append(sip_file)

        if self.in_main_module:
            self._sip_files.append(sip_file)

        return self._input

    def _unexpected_eof_error(self):
        """ Record an error caused by an unexpected EOF. """

        self._error_log.log("unexpected end of file",
                SourceLocation(self._sip_file))


class ModuleState:
    """ Encapsulate the parser-related state for a module. """

    def __init__(self, module, sip_file):
        """ Initialise the state. """

        self.module = module
        self.sip_file = sip_file

        self.all_raise_py_exception = False
        self.auto_py_name_rules = []
        self.call_super_init = None
        self.default_encoding = None
        self.kw_args = KwArgs.NONE
        self.nr_timelines = 0


class ScopeState:
    """ Encapsulate the parser-related state for a scope. """

    def __init__(self, scope, access_specifier):
        """ Initialise the state. """

        self.scope = scope
        self.access_specifier = access_specifier
        self.pyqt_method_specifier = None


class UnexpectedEOF(Exception):
    """ This is raised by p_error() when an unexpected EOF is seen. """
