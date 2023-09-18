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



from copy import copy

from ..error_log import ErrorLog
from ..instantiations import instantiate_type_hints
from ..python_slots import (is_hash_return_slot, is_int_return_slot,
        is_inplace_number_slot, is_rich_compare_slot, is_ssize_return_slot,
        is_void_return_slot, is_zero_arg_slot)
from ..scoped_name import ScopedName
from ..specification import (AccessSpecifier, Argument, ArgumentType, ClassKey,
        Constructor, IfaceFileType, MappedType, Member, PyQtMethodSpecifier,
        PySlot, Signature, Transfer, ValueType, VirtualHandler,
        VirtualOverload, VisibleMember, WrappedClass)
from ..templates import (encoded_template_name, same_template_signature,
        template_code, template_code_blocks, template_expansions)
from ..utils import (append_iface_file, argument_as_str, cached_name,
        find_iface_file, find_method, same_argument_type, same_base_type,
        same_signature, search_typedefs)


def resolve(spec, modules):
    """ Resolve all types of a parsed specification and create additional views
    so that code can be generated.
    """

    error_log = ErrorLog()

    # Build the list of all imports for each module.
    for mod in modules:
        _set_all_imports(mod, error_log)

        # Set the base name of the module.  This is done for efficiency.
        mod.py_name = mod.fq_py_name.name.split('.')[-1]

    # Set the default meta-type for the main module if it doesn't have one
    # explicitly set.
    if spec.module.default_metatype is None:
        for mod in spec.module.all_imports:
            if mod.default_metatype is None:
                continue

            if spec.module.default_metatype is None:
                spec.module.default_metatype = mod.default_metatype
            elif spec.module.default_metatype is not mod.default_metatype:
                error_log.log(
                        "'{0}' module has imported different default meta-types '{1}' and '{2}'".format(
                                spec.module.fq_py_name,
                                spec.module.default_metatype,
                                mod.default_metatype))

    # Check each class has been defined.
    for klass in spec.classes:
        if klass.iface_file.module is None:
            error_log.log(
                    "class '{0}' has not been defined".format(
                            klass.iface_file.fq_cpp_name))

        # Mark any QObject class.  This flag will ripple through all derived
        # classes when we set the hierarchy.
        if klass.iface_file.fq_cpp_name.base_name == 'QObject':
            klass.is_qobject = True

    # The class list has the main module's classes at the front and the ones
    # from the module at the most nested %Import at the end.  Set the MRO for
    # each class and re-order the list of classes so that no class appears
    # before a super class or an enclosing scope class.
    reversed_classes = reversed(spec.classes)
    spec.classes = []

    for klass in reversed_classes:
        # Ignore undefined classes.
        if klass.iface_file.module is None:
            continue

        _set_mro(spec, klass, error_log)

    # Resolve the various types in the modules.
    _resolve_module(spec, spec.module, error_log)

    # Handle default ctors now that the argument types are resolved.
    for klass in spec.classes:
        if klass.no_default_ctors or klass.is_opaque or klass.iface_file.type is IfaceFileType.NAMESPACE:
            continue

        _add_default_copy_ctor(klass)

    # Add any automatically generated methods.
    for klass in spec.classes:
        for overload in klass.overloads:
            if overload.is_auto_generated:
                _add_auto_overload(spec, klass, overload)

    # Move casts and slots around to their correct classes (if in the same
    # module) or create proxies for them (if cross-module).
    _move_main_module_casts_slots(spec, error_log)

    # Automatically generate missing complementary slots.
    for klass in spec.classes:
        _add_complementary_slots(spec, klass)

    for klass in spec.module.proxies:
        _add_complementary_slots(spec, klass)

    # Generate the different class views.
    for klass in spec.classes:
        if klass.iface_file.type is IfaceFileType.CLASS:
            if klass.needs_shadow and not klass.is_incomplete and klass.dtor is not AccessSpecifier.PRIVATE and klass.can_create:
                klass.has_shadow = True

            # Get the list of visible Python member functions.
            _get_visible_py_members(spec, klass)

            # Get the virtual members.
            if klass.needs_shadow:
                _get_virtuals(spec, klass, error_log)

        elif klass.iface_file.type is IfaceFileType.NAMESPACE:
            for overload in klass.overloads:
                _iface_files_are_used_by_overload(spec, klass.iface_file.used,
                        overload)

    for mod in modules:
        # Create the list of numbered types sorted by type name.
        _create_sorted_numbered_types(spec, mod, error_log)

        for overload in mod.overloads:
            _iface_files_are_used_by_overload(spec, mod.used, overload)

        # Update proxies with some information from the real classes.
        for klass in mod.proxies:
            klass.iface_file.type_nr = klass.real_class.iface_file.type_nr

    # Additional class specific checks.
    for klass in spec.classes:
        # Ignore undefined classes.
        if klass.iface_file.module is None:
            continue

        _check_helpers(spec, klass)
        _check_properties(klass, error_log)

    # Number the exceptions as they will be seen by the main module.
    for exception in spec.exceptions:
        exception_mod = exception.iface_file.module

        # Include the %TypeHeaderCode for exceptions defined in the main
        # module.
        if spec.abi_version >= (13, 1) or (spec.abi_version >= (12, 9) and spec.abi_version < (13, 0)):
            if exception_mod is spec.module:
                append_iface_file(spec.module.used, exception.iface_file)

        # Skip those that don't require a Python exception object to be
        # created.
        if exception.iface_file.type is not IfaceFileType.EXCEPTION:
            continue

        if exception.builtin_base_exception is None and exception.defined_base_exception is None:
            continue

        if exception_mod is spec.module or exception.needed:
            exception.exception_nr = exception_mod.nr_exceptions
            exception_mod.nr_exceptions += 1

    # For PyQt6 mark all enum interface files as being used.
    if 'PyQt6' in spec.plugins:
        for enum in spec.enums:
            if enum.module is spec.module:
                _enum_iface_file_is_used(enum, spec.module)

    # Create the name cache as seen by the legacy code.
    name_cache = spec.name_cache
    spec.name_cache = []

    for k in sorted(name_cache.keys(), reverse=True):
        spec.name_cache.extend(sorted(name_cache[k], key=lambda k: k.name))

    _set_string_pool_offsets(spec)

    # Raise an exception for any errors.
    error_log.as_exception()


def _resolve_module(spec, mod, error_log, seen=None):
    """ Resolve a module and the modules it imports. """

    if seen is None:
        seen = []
    elif mod in seen:
        # Nothing left to do.
        return

    # The modules on which this one depends must be done first because they
    # might generate new template-based types and they must be defined in the
    # right module.
    for imported_mod in mod.imports:
        _resolve_module(spec, imported_mod, error_log, seen=seen)

    # Resolve typedefs, variables and global functions.
    _resolve_typedefs(spec, mod, error_log)
    _resolve_variables(spec, mod, error_log)
    _resolve_scope_overloads(spec, mod.overloads, error_log)

    # Resolve class ctors, functions and casts.
    for klass in spec.classes:
        if klass.iface_file.module is mod:
            _resolve_ctors(spec, klass, error_log)

            # Handle any dtor exceptions.
            _set_needed_exceptions(spec, mod, klass.dtor_throw_args)

            _resolve_scope_overloads(spec, klass.overloads, error_log,
                    scope=klass)
            _transform_casts(spec, klass, error_log)

    # Resolve mapped types based on templates.
    _resolve_mapped_types(spec, mod, error_log)

    seen.append(mod)


def _set_string_pool_offsets(spec):
    """ Set the offset into the string pool for every used name. """

    offset = 0

    for name in spec.name_cache:
        if not name.used:
            continue

        name_len = len(name.name)

        # See if the tail of a previous used name could be used instead.
        for prev_name in spec.name_cache:
            prev_name_len = len(prev_name.name)

            if prev_name_len <= name_len:
                break

            if not prev_name.used or prev_name.is_substring:
                continue

            if prev_name.name.endswith(name.name):
                name.is_substring = True
                name.offset = prev_name.offset + prev_name_len - name_len;
                break

        if not name.is_substring:
            name.offset = offset
            offset += name_len + 1


# The map of slots and their compliments.
_SLOT_MAP = {
    PySlot.LT:  (PySlot.GE, '__ge__'),
    PySlot.LE:  (PySlot.GT, '__gt__'),
    PySlot.GT:  (PySlot.LE, '__le__'),
    PySlot.GE:  (PySlot.LT, '__lt__'),
    PySlot.EQ:  (PySlot.NE, '__ne__'),
    PySlot.NE:  (PySlot.EQ, '__eq__'),
}

def _add_complementary_slots(spec, klass):
    """ Add any missing complementary slots to a class.  This emulates the C++
    behaviour of automatically interpreting (for example) >= as !<.
    """

    for member in list(klass.members):
        try:
            compl, compl_name = _SLOT_MAP[member.py_slot]
        except KeyError:
            continue

        _add_complementary_slot(spec, klass, member, compl, compl_name)


def _add_complementary_slot(spec, klass, member, compl, compl_name):
    """ Add a complementary slot if it is missing. """

    member2 = None

    for overload in klass.overloads:
        if overload.common is not member or overload.is_complementary or overload.method_code is not None:
            continue

        # Try and find an existing complementary slot.
        for overload2 in klass.overloads:
            if overload2.common.py_slot is compl and same_signature(spec, overload.py_signature, overload2.py_signature):
                break
        else:
            # There is no explicit complementary slot so create a new member if
            # needed.
            if member2 is None:
                for member2 in klass.members:
                    if member2.py_slot is compl:
                        break
                else:
                    member2 = copy(member)

                    member2.py_name = cached_name(spec, compl_name)
                    member2.py_slot = compl

                    if member.py_name.used:
                        member2.py_name.used = True

                    klass.members.insert(0, member2)

            # Create the complementary slot.
            overload2 = copy(overload)

            overload2.is_virtual = False
            overload2.is_complementary = True
            overload2.common = member2
            overload2.cpp_name = compl_name

            klass.overloads.insert(0, overload2)
            break


def _check_helpers(spec, klass):
    """ See if a class supports array and copy helpers. """

    # We disregard classes that are abstract, or have a private assignment
    # operator or don't have a public dtor.
    if klass.is_abstract:
        return

    if klass.cannot_assign:
        return

    if klass.dtor is not AccessSpecifier.PUBLIC:
        return

    # See if the class has a default ctor and a public copy ctor.
    public_default_ctor = public_copy_ctor = False

    for ctor in klass.ctors:
        if ctor.cpp_signature is None or ctor.access_specifier is not AccessSpecifier.PUBLIC:
            continue

        if len(ctor.cpp_signature.args) == 0 or ctor.cpp_signature.args[0].default_value is not None:
            # The ctor either has no arguments or all arguments have defaults.
            public_default_ctor = True
        elif len(ctor.cpp_signature.args) == 1:
            arg = ctor.cpp_signature.args[0]

            if arg.type is not ArgumentType.CLASS:
                continue

            if arg.definition is klass and arg.is_reference and arg.is_const and len(arg.derefs) == 0 and arg.default_value is None:
                public_copy_ctor = True

    if public_default_ctor:
        klass.needs_array_helper = True
        append_iface_file(klass.iface_file.module.used, klass.iface_file)

    if public_copy_ctor:
        klass.needs_copy_helper = True
        append_iface_file(klass.iface_file.module.used, klass.iface_file)


def _set_all_imports(mod, error_log, seen=None):
    """ Set the list of all imports for a module.  The list is ordered so that
    a module appears before any module that imports it.
    """

    # Check for recursive imports.
    if seen is None:
        seen = []
    elif mod in seen:
        error_log.log(
                "module '{0}' is imported recursively".format(mod.fq_py_name))

    seen.append(mod)

    # Make sure all the direct imports are done first.
    for direct_import in mod.imports:
        _set_all_imports(direct_import, error_log, seen=seen)

    # Now build the list from our direct imports' lists but ignoring
    # duplicates.
    def _append_unique_import(imported):
        if imported not in mod.all_imports:
            mod.all_imports.append(imported)

    for direct_import in mod.imports:
        for imported in direct_import.all_imports:
            _append_unique_import(imported)

        _append_unique_import(direct_import)

    seen.remove(mod)


def _move_main_module_casts_slots(spec, error_log):
    """ Move the casts and slots to the correct place for a main module (ie.
    the one we are generating code for).
    """

    for klass in spec.classes:
        if klass.iface_file.module is spec.module:
            _move_class_casts(spec, klass, error_log)

    for member in spec.module.global_functions:
        if member.py_slot is not None and member.module is spec.module:
            _move_global_slot(spec, member, error_log)


def _move_class_casts(spec, klass, error_log):
    """ Move any class casts to its correct class, or publish as a ctor
    extender.
    """

    for cast in klass.casts:
        dst_klass = cast.definition

        # Create the new ctor.
        arg = Argument(ArgumentType.CLASS, definition=klass,
                derefs=list(cast.derefs), is_reference=cast.is_reference,
                is_const=cast.is_const, is_in=True,
                source_location=cast.source_location)
        signature = Signature(args=[arg])
        ctor = Constructor(AccessSpecifier.PUBLIC, py_signature=signature,
                cpp_signature=signature, is_cast=True)

        # If the destination class is in a different module then use a proxy.
        if dst_klass.iface_file.module is not spec.module:
            _iface_file_is_used(spec.module.used, arg)
            dst_klass = _get_proxy(spec.module, dst_klass)
            ctor.no_typehint = True

        _iface_file_is_used(dst_klass.iface_file.used, arg)

        # Check it hasn't already been defined.
        for dst_ctor in dst_klass.ctors:
            if same_signature(spec, dst_ctor.py_signature, ctor.py_signature, strict=False):
                error_log.log(
                        " operator '{0}::{0}({1})' already defined".format(
                                dst_klass.iface_file.fq_cpp_name,
                                klass.iface_file.fq_cpp_name))

        dst_klass.ctors.append(ctor)


def _move_global_slot(spec, global_slot, error_log):
    """ If possible, move a global slot to its correct class. """

    for overload in list(spec.module.overloads):
        if overload.common is not global_slot:
            continue

        # We know that the slot has the right number of arguments, but the
        # first or second one needs to be a class or enum defined in the same
        # module.  Otherwise we leave it as it is and publish it as a slot
        # extender.
        arg0 = overload.py_signature.args[0]

        try:
            arg1 = overload.py_signature.args[1]
        except IndexError:
            arg1 = None

        arg_members = None
        arg_overloads = None
        arg_module = None
        arg_enum = None
        is_second = False

        if arg0.type is ArgumentType.CLASS:
            arg_members = arg0.definition.members
            arg_overloads = arg0.definition.overloads
            arg_module = arg0.definition.iface_file.module

        elif arg0.type is ArgumentType.ENUM:
            arg_members = arg0.definition.slots
            arg_overloads = arg0.definition.overloads;
            arg_module = arg0.definition.module
            arg_enum = arg0.definition

        elif arg1 is None:
            if arg0.type is ArgumentType.NONE:
                # The error has already been logged.
                pass
            else:
                _log_overload_error(error_log,
                        "the argument must be a class or enum", overload)
            continue

        elif arg1.type is ArgumentType.CLASS:
            arg_members = arg1.definition.members
            arg_overloads = arg1.definition.overloads
            arg_module = arg1.definition.iface_file.module
            is_second = True

        elif arg1.type is ArgumentType.ENUM:
            arg_members = arg1.definition.slots
            arg_overloads = arg1.definition.overloads
            arg_module = arg1.definition.module
            arg_enum = arg1.definition
            is_second = True

        else:
            if arg0.type is ArgumentType.NONE or arg1.type is ArgumentType.NONE:
                # The error has already been logged.
                pass
            else:
                _log_overload_error(error_log,
                        "one of the arguments must be a class or enum",
                        overload)
            continue

        # For rich comparisons the first argument must be a class or an enum.
        # For cross-module slots then it may only be a class.  (This latter
        # limitation is artificial, but is unlikely to be a problem in
        # practice.)
        if is_rich_compare_slot(global_slot.py_slot):
            if is_second:
                _log_overload_error(error_log,
                        "argument 1 must be a class or enum", overload)
                continue

            if arg_module is not global_slot.module and arg0.type is ArgumentType.ENUM:
                _log_overload_error(error_log, "argument 1 must be a class",
                        overload)
                continue

        if arg_module is not global_slot.module:
            if is_rich_compare_slot(global_slot.py_slot):
                proxy = _get_proxy(arg_module, arg0.definition)

                # Create a new proxy member if needed.
                for proxy_member in proxy.members:
                    if proxy_member.py_slot is global_slot.py_slot:
                        break
                else:
                    proxy_member = Member(arg_module, global_slot.py_name,
                            py_slot=global_slot.py_slot,
                            namespace_iface_file=global_slot.namespace_iface_file)

                    proxy.members.insert(0, proxy_member)

                # Remove the overload from the list.
                spec.module.overloads.remove(overload)

                # Add the overload to the proxy.
                overload.common = proxy_member
                overload.no_typehint = True

                proxy.overloads.insert(0, overload)

                # Remove the overload's first argument.
                del overload.py_signature.args[0]

            continue

        # Remove from the list.
        spec.module.overloads.remove(overload)

        if arg_enum is not None:
            _enum_iface_file_is_used(arg_enum, arg_module)
            arg_enum.py_name.used = True

        # See if there is already a member or create a new one.
        inject_equality_slot = False

        for arg_member in arg_members:
            if arg_member.py_slot is global_slot.py_slot:
                break
        else:
            arg_member = copy(global_slot)

            arg_member.module = arg_module
            arg_members.insert(0, arg_member)

            # Legacy enum members, when accessed as scoped values, are created
            # on the fly.  By default these members compare for equality
            # correctly (ie. 'E.M == E.M' works as expected).  However if there
            # is another equality operator defined then it will fail so we have
            # to explicitly inject the comparison.
            if spec.abi_version < (13, 0) and arg0.type is ArgumentType.ENUM and arg_member.py_slot is PySlot.EQ and not is_second:
                inject_equality_slot = True

        # Move the overload to the end of the destination list.
        if is_second:
            overload.is_reflected = True

        overload.access_specifier = AccessSpecifier.PUBLIC
        overload.common = arg_member
        overload.is_global = True

        arg_overloads.append(overload)

        # Inject an additional equality slot if necessary.
        if inject_equality_slot:
            eq_overload = copy(overload)
            eq_overload.py_signature = copy(eq_overload.py_signature)
            eq_overload.py_signature.args = copy(eq_overload.py_signature.args)
            eq_overload.cpp_signature = eq_overload.py_signature

            eq_overload.py_signature.args[0].derefs = [False]
            del eq_overload.py_signature.args[1]

            arg_overloads.append(eq_overload)

        # Remove the first argument of inplace numeric operators and comparison
        # operators.
        if is_inplace_number_slot(arg_member.py_slot) or is_rich_compare_slot(arg_member.py_slot):
            # Remember if the argument was a pointer.
            if len(arg0.derefs) > 0:
                overload.dont_deref_self = True

            del overload.py_signature.args[0]

        # Remove the only argument of unary operators.
        if is_zero_arg_slot(arg_member.py_slot):
            del overload.py_signature.args[0]


def _get_proxy(mod, klass):
    """ Create a proxy for a class if it doesn't already exist.  Proxies are
    used as containers for cross-module extenders.
    """

    for proxy in mod.proxies:
        if proxy.iface_file is klass.iface_file:
            return proxy

    proxy = WrappedClass(klass.iface_file, klass.py_name, scope=klass.scope,
            mro=klass.mro, real_class=klass, superclasses=klass.superclasses)

    mod.proxies.insert(0, proxy)

    return proxy


def _add_auto_overload(spec, auto_klass, auto_overload):
    """ Add an overload that is automatically generated. """

    # Find every class that has this one in its hierarchy.
    for klass in spec.classes:
        if klass is auto_klass:
            continue

        for mro_klass in klass.mro:
            if mro_klass is auto_klass:
                # Another overload may already exist.
                for member in klass.members:
                    if member.py_name is auto_overload.common.py_name:
                        break
                else:
                    member = copy(auto_overload.common)
                    member.module = klass.iface_file.module
                    klass.members.insert(0, member)

                overload = copy(auto_overload)
                overload.common = member
                overload.is_autogenerated = False
                klass.overloads.insert(0, overload)

                if klass.iface_file.module is spec.module:
                    member.py_name.used = True

                break


def _set_mro(spec, klass, error_log, seen=None):
    """ Set the MRO for a class and add it to the list of classes so that it is
    after any classes it depends on.
    """

    # See if it has already been done.
    if klass in spec.classes:
        return

    # Initialise the detection of recursive hierarchies.
    if seen is None:
        seen = []

    # Handle any enclosing scope.
    if klass.scope is not None:
        _set_mro(spec, klass.scope, error_log, seen=seen)

        if klass.scope.deprecated:
            klass.deprecated = True

    if klass.iface_file.type is IfaceFileType.CLASS:
        # The first thing is itself.
        klass.mro.append(klass)

        if klass.convert_to_subclass_code is not None:
            klass.subclass_base = klass

        # Now do it's super-classes.
        seen.append(klass)

        for superklass in klass.superclasses:
            if superklass in seen:
                error_log.log(
                        "recursive class hierarchy detected: '{0}' and '{1}'".format(
                                klass.iface_file.fq_cpp_name,
                                superklass.iface_file.fq_cpp_name))
                continue

            # Unions cannot be super-classes.
            if superklass.class_key is ClassKey.UNION:
                error_log.log(
                        "union '{0}' cannot be a super-class".format(
                                superklass.iface_file.fq_cpp_name))

            # Make sure the super-class's hierarchy has been done. */
            _set_mro(spec, superklass, error_log, seen=seen)

            # Append the super-class's MRO.
            for superklass_mro in superklass.mro:
                if superklass_mro not in klass.mro:
                    klass.mro.append(superklass_mro)

                if klass.iface_file.module is spec.module:
                    superklass_mro.iface_file.needed = True

                if superklass_mro.deprecated:
                    klass.deprecated = True

                # If the super-class is a QObject sub-class then this one is as
                # well.
                if superklass_mro.is_qobject:
                    klass.is_qobject = True

                # If the super-class can't be assigned to then this one cannot
                # either.
                if superklass_mro.cannot_assign:
                    klass.cannot_assign = True

                # If the super-class needs a shadow then this one should have
                # one as well.
                if superklass_mro.needs_shadow:
                    klass.needs_shadow = True

                # Ensure that the sub-class base class is the furthest up the
                # hierarchy.
                if superklass_mro.subclass_base is not None:
                    klass.subclass_base = superklass_mro.subclass_base;

        seen.remove(klass)

        # If the class doesn't have an explicit meta-type then inherit from the
        # module's default.
        if klass.metatype is None and len(klass.superclasses) == 0:
            # The class may not have been defined.
            if klass.iface_file.module is not None:
                klass.metatype = klass.iface_file.module.default_metatype

        if klass.metatype is not None and klass.iface_file.module is spec.module:
            klass.metatype.used = True

        # If the class doesn't have an explicit super-type then inherit from
        # the module's default.
        if klass.supertype is None and len(klass.superclasses) == 0:
            # The class may not have been defined.
            if klass.iface_file.module is not None:
                klass.supertype = klass.iface_file.module.default_supertype

        if klass.supertype is not None:
            # If the super-type ends with 'sip.wrapper' then assume it is the
            # default.
            if klass.supertype.name.endswith('sip.wrapper'):
                klass.supertype = None

        if klass.supertype is not None and klass.iface_file.module is spec.module:
            klass.supertype.used = True

    # Make sure that the module in which a sub-class convertor will be created
    # knows about the base class.
    if klass.subclass_base is not None:
        append_iface_file(klass.iface_file.module.used,
                klass.subclass_base.iface_file)

    # We can't have a shadow if the specification is incomplete, there is a
    # private dtor, there are no non-private ctors or there are private
    # abstract methods.
    if klass.is_incomplete or klass.dtor is AccessSpecifier.PRIVATE or not klass.can_create:
        klass.has_shadow = False
    else:
        # Note that we should be able to provide better support for abstract
        # private methods than we do at the moment.
        for overload in klass.overloads:
            if overload.is_abstract and overload.access_specifier is AccessSpecifier.PRIVATE:
                klass.has_shadow = False

                # It also means we cannot create an instance from Python.
                klass.can_create = False

                break

    # Add it to the new list of classes.
    spec.classes.append(klass)


def _resolve_typedefs(spec, mod, error_log):
    """ Resolve the base types for all typedefs of a module. """

    for typedef in spec.typedefs:
        if typedef.module is mod:
            _resolve_type(spec, typedef.module, typedef.scope, typedef.type,
                    error_log)


def _resolve_mapped_types(spec, mod, error_log):
    """ Resolve the data types for mapped types based on a template. """

    for mapped_type in spec.mapped_types:
        if mapped_type.iface_file.module is mod:
            if mapped_type.type.type is ArgumentType.TEMPLATE:
                _resolve_mapped_type_types(spec, mapped_type, error_log)
            else:
                _resolve_scope_overloads(spec, mapped_type.overloads,
                        error_log, scope=mapped_type)


def _resolve_ctors(spec, klass, error_log):
    """ Resolve the data types for a class's ctors. """

    for ctor in klass.ctors:
        _resolve_ctor_types(spec, klass, ctor, error_log)

        # Now check that the Python signature doesn't conflict with an earlier
        # one.  If there is %MethodCode then assume that it will handle any
        # potential conflicts.
        if ctor.method_code is None:
            for previous_ctor in klass.ctors:
                if previous_ctor is ctor:
                    break

                if previous_ctor.method_code is not None:
                    continue

                sig_state = _same_python_signature(spec,
                        previous_ctor.py_signature, ctor.py_signature)

                if sig_state is None:
                    # The error has already been logged.
                    break

                if sig_state:
                    error_log.log(
                            "class '{0}' has ctors with the same Python signature".format(
                                    klass.iface_file.fq_cpp_name))
                    break

        if klass.deprecated:
            ctor.deprecated = True


def _transform_casts(spec, klass, error_log):
    """ Resolve the data type for a list of casts. """

    for cast in klass.casts:
        _resolve_type(spec, klass.iface_file.module, klass, cast, error_log)

        if cast.type is ArgumentType.NONE:
            # The error has already been logged.
            pass
        elif cast.type is not ArgumentType.CLASS:
            error_log.log(
                    "operator cast '{0}' must be to a class".format(
                            klass.iface_file.fq_cpp_name))


def _add_default_copy_ctor(klass):
    """ Add a default copy ctor if required. """

    # See if there is a private copy ctor in the hierarchy.
    for mro_klass in klass.mro:
        for ctor in mro_klass.ctors:
            # See if is a copy ctor.
            if len(ctor.py_signature.args) != 1:
                continue

            arg = ctor.py_signature.args[0]
 
            if len(arg.derefs) == 0 and arg.is_reference:
                if arg.type is ArgumentType.CLASS and arg.definition.iface_file is mro_klass.iface_file:
                    break
        else:
            continue

        # If the copy ctor is private then the class can't be copied.
        if ctor.access_specifier is AccessSpecifier.PRIVATE:
            klass.cannot_copy = True
            return

        # If the ctor is in the class itself then there is nothing to do.
        if mro_klass is klass.mro[0]:
            return
 
        # Otherwise we need to create a default.
        break
 
    # Create a default public copy ctor.
    arg = Argument(ArgumentType.CLASS, definition=klass, is_reference=True,
            is_const=True, is_in=True)
    result = Argument(ArgumentType.VOID)
    signature = Signature(args=[arg], result=result)
    ctor = Constructor(AccessSpecifier.PUBLIC, py_signature=signature,
            cpp_signature=signature)
 
    if klass.deprecated:
        ctor.deprecated = True

    if not klass.is_abstract:
        klass.can_create = True

    # Append it to the list.
    klass.ctors.append(ctor)


def _resolve_scope_overloads(spec, overloads, error_log, scope=None):
    """ Resolve the data types for a scope's overloads. """

    for overload in overloads:
        _resolve_func_types(spec, overload.common.module, scope, overload,
                error_log)

        # Now check that the Python signature doesn't conflict with an earlier
        # one.  If there is %MethodCode then assume that it will handle any
        # potential conflicts.
        if overload.method_code is None and spec.is_strict:
            for previous_overload in overloads:
                if previous_overload is overload:
                    break

                if previous_overload.common is not overload.common:
                    continue

                if previous_overload.method_code is not None:
                    continue

                sig_state = _same_python_signature(spec,
                        previous_overload.py_signature, overload.py_signature)

                if sig_state is None:
                    # The error has already been logged.
                    break

                if sig_state:
                    _log_overload_error(error_log,
                            "has overloaded functions with the same Python signature",
                            overload, scope=scope)
                    break

        if isinstance(scope, WrappedClass):
            if scope.deprecated:
                overload.deprecated = True

            if overload.is_abstract:
                scope.is_abstract = True


def _resolve_variables(spec, mod, error_log):
    """ Resolve the data types for the variables of a module. """

    for variable in spec.variables:
        if variable.module is mod:
            _resolve_variable_type(spec, variable, error_log)


def _get_visible_py_members(spec, klass):
    """ Set the list of visible Python member functions for a class. """

    for mro_klass in klass.mro:
        for member in mro_klass.members:
            # See if it is already in the list.  This has the desired side
            # effect of eliminating any functions that have an implementation
            # closer to this class in the hierarchy.  This is the only reason
            # to define private functions.
            for visible_member in klass.visible_members:
                if visible_member.member.py_name is member.py_name:
                    break
            else:
                visible_member = VisibleMember(member, mro_klass)
                klass.visible_members.insert(0, visible_member)

                for overload in mro_klass.overloads:
                    if overload.common is member:
                        need_types = False

                        # If the visible overload is abstract then it hasn't
                        # had a concrete implementation so this class must also
                        # be abstract.
                        if overload.is_abstract:
                            klass.is_abstract = True

                        if klass.iface_file.module is spec.module and (klass is mro_klass or (overload.access_specifier is AccessSpecifier.PROTECTED and klass.has_shadow)):
                            need_types = True
                            member.py_name.used = True

                        _iface_files_are_used_by_overload(spec,
                                klass.iface_file.used, overload,
                                need_types=need_types);


def _get_virtuals(spec, klass, error_log):
    """ Get all the virtuals for a particular class. """

    # Copy the collected virtuals of each super-class updating from what we
    # find in this class.
    for superklass in klass.superclasses:
        for virtual_overload in superklass.virtual_overloads:
            implicit = True

            for overload in klass.overloads:
                if virtual_overload.overload.cpp_name != overload.cpp_name:
                    continue

                implicit = False

                if overload.is_final:
                    break

                # See if it re-implements rather than hides.
                if _same_cpp_overload(spec, virtual_overload.overload, overload):
                    # What if is is private?
                    overload.is_virtual = True
                    overload.is_virtual_reimplementation = True

                    # Use the base implementation's virtual handler code if
                    # there is any.  We cannot just use its virtual handler
                    # because this re-implementation may have different
                    # annotations which means the complete handler would be
                    # different.  In practice there is no reason why it would
                    # be different (and maybe this should be detected as an
                    # error) but if they are the same then the same handler
                    # will eventually be chosen.
                    if overload.virtual_catcher_code is None:
                        overload.virtual_catcher_code = virtual_overload.overload.virtual_catcher_code

                    # Use the base implementation's virtual error handler if
                    # one isn't explicitly specified.
                    if overload.virtual_error_handler is None:
                        overload.virtual_error_handler = virtual_overload.overload.virtual_error_handler

                    _add_virtual_overload(spec, overload, klass, error_log)

            # Add it if it wasn't explicitly mentioned in the class.
            if implicit:
                _add_virtual_overload(spec, virtual_overload.overload, klass,
                        error_log)

    # Handle any new virtuals.
    for overload in klass.overloads:
        if overload.is_virtual and not overload.is_virtual_reimplementation and not overload.is_final:
            _add_virtual_overload(spec, overload, klass, error_log)


def _add_virtual_overload(spec, overload, klass, error_log):
    """ Add an overload to the list of virtuals for a class. """

    # If this class is defined in the main module then make sure the virtuals
    # have a handler.
    if klass.iface_file.module is spec.module:
        virtual_handler = _get_virtual_handler(spec, overload, klass,
                error_log)

        # Make sure we get the name.
        overload.common.py_name.used = True

        # Make sure we have the interface files and type definitions for the
        # virtual handler.
        _iface_files_are_used_by_overload(spec, spec.module.used, overload,
                need_types=True)
    else:
        virtual_handler = None

    # Add it to the class.
    virtual_overload = VirtualOverload(overload, virtual_handler)
    klass.virtual_overloads.insert(0, virtual_overload)


def _get_virtual_error_handler(spec, overload, klass, error_log):
    """ Get the virtual error handler for a function. """

    # Handle the trivial case.
    if overload.no_virtual_error_handler:
        return None

    klass_mod = klass.iface_file.module

    # Check the function itself.
    handler_name = overload.virtual_error_handler

    if handler_name is None:
        # Check the class hierarchy.
        for mro_klass in klass.mro:
            handler_name = mro_klass.virtual_error_handler

            if handler_name is not None:
                break
        else:
            # Check the class's module.
            handler_name = klass_mod.default_virtual_error_handler

            if handler_name is None:
                # Check the module hierarchy.
                for mod in klass_mod.all_imports:
                    handler_name = mod.default_virtual_error_handler

                    if handler_name is not None:
                        break
                else:
                    return None

    # Find the handler with the name.
    for handler in spec.virtual_error_handlers:
        if handler.name == handler_name:
            break
    else:
        error_log.log(
                "unknown virtual error handler '{0}'".format(handler_name))
        return None

    # Assign it an index if we need to import the handler.
    if klass_mod is not handler.module and handler.handler_nr < 0:
        handler.handler_nr = handler.module.nr_virtual_error_handlers
        handler.module.nr_virtual_error_handlers += 1

    return handler


def _get_virtual_handler(spec, overload, klass, error_log):
    """ Get the virtual handler for an overload. """

    # See if there is an existing handler that is suitable.
    for handler in spec.virtual_handlers:
        if _check_virtual_handler(spec, overload, handler):
            return handler

    # Create a new one.
    handler = VirtualHandler(overload.cpp_signature, overload.py_signature,
            overload.virtual_catcher_code,
            _get_virtual_error_handler(spec, overload, klass, error_log))

    handler.handler_nr = spec.nr_virtual_handlers
    spec.nr_virtual_handlers += 1

    if overload.factory or overload.transfer is Transfer.TRANSFER_BACK:
        handler.transfer_result = True

    if overload.abort_on_exception:
        handler.abort_on_exception = True

    spec.virtual_handlers.insert(0, handler)

    return handler


def _check_virtual_handler(spec, overload, virtual_handler):
    """ Return True if a virtual handler is appropriate for an overload. """

    if overload.virtual_catcher_code is not virtual_handler.virtual_catcher_code:
        return False

    # If the overload has an explicit error handler then it must be the same as
    # the candidate.
    if overload.virtual_error_handler is not None:
        if virtual_handler.virtual_error_handler is None or overload.virtual_error_handler != virtual_handler.virtual_error_handler.name:
            return False

    if (overload.factory or overload.transfer is Transfer.TRANSFER_BACK) and  not virtual_handler.transfer_result:
        return False

    if overload.abort_on_exception is not virtual_handler.abort_on_exception:
        return False

    if not same_argument_type(spec, overload.py_signature.result, virtual_handler.py_signature.result):
        return False

    if overload.py_signature.result.allow_none is not virtual_handler.py_signature.result.allow_none:
        return False

    if overload.py_signature.result.disallow_none is not virtual_handler.py_signature.result.disallow_none:
        return False

    if not same_signature(spec, overload.py_signature, virtual_handler.py_signature):
        return False

    # Take into account the argument directions in the Python signatures.
    for arg1, arg2 in zip(overload.py_signature.args, virtual_handler.py_signature.args):
        if arg1.is_in is not arg2.is_in:
            return False

        if arg1.is_out is not arg2.is_out:
            return False

    if overload.py_signature is overload.cpp_signature and virtual_handler.py_signature is virtual_handler.cpp_signature:
        return True

    if not same_argument_type(spec, overload.cpp_signature.result, virtual_handler.cpp_signature.result):
        return False

    return same_signature(spec, overload.cpp_signature,
            virtual_handler.cpp_signature)


def _resolve_mapped_type_types(spec, mapped_type, error_log):
    """ Resolve the types of a mapped type based on a template. """

    template_types = mapped_type.type.definition.types

    for template_type in template_types.args:
        # Leave templates as they are.
        if template_type.type is not ArgumentType.TEMPLATE:
            _resolve_type(spec, mapped_type.iface_file.module, None,
                    template_type, error_log, allow_defined=True)

    # Make sure that the signature result won't cause problems.
    template_types.result = None

    _iface_files_are_used_by_signature(mapped_type.iface_file.used,
            template_types)


def _resolve_ctor_types(spec, scope, ctor, error_log):
    """ Resolve the types of a ctor. """

    # Handle any exceptions.
    _set_needed_exceptions(spec, scope.iface_file.module, ctor.throw_args)

    # Handle any C++ signature.
    if ctor.cpp_signature is not None and ctor.cpp_signature is not ctor.py_signature:
        for arg in ctor.cpp_signature.args:
            _resolve_type(spec, scope.iface_file.module, scope, arg, error_log,
                    allow_defined=True)
 
    # Handle the Python signature.
    for arg_nr, arg in enumerate(ctor.py_signature.args):
        _resolve_type(spec, scope.iface_file.module, scope, arg, error_log)

        if arg.type is ArgumentType.NONE:
            # The error has already been logged.
            continue

        if not _supported_type(scope, None, arg, error_log):
            error_log.log(
                    "argument {0} of ctor '{1}' has an unsupported type for a Python signature - provide a valid type, %MethodCode and a C++ signature".format(
                            arg_nr + 1, scope.iface_file.fq_cpp_name))
            continue

        _iface_file_is_used(scope.iface_file.used, arg)
        _scope_default_value(spec, scope, arg)


def _resolve_func_types(spec, mod, scope, overload, error_log):
    """ Resolve the types of a function. """

    # Handle any exceptions.
    _set_needed_exceptions(spec, mod, overload.throw_args)

    # Handle any C++ signature.
    if overload.cpp_signature is not overload.py_signature:
        result = overload.cpp_signature.result

        _resolve_type(spec, mod, scope, result, error_log, allow_defined=True)

        if (result.type is not ArgumentType.VOID or len(result.derefs) != 0) and overload.is_virtual and not _supported_type(scope, overload, result, error_log) and overload.virtual_catcher_code is None:
            _log_overload_error(error_log,
                    "has an unsupported virtual function return type - provide %VirtualCatcherCode",
                    overload, scope=scope)

        for arg in overload.cpp_signature.args:
            _resolve_type(spec, mod, scope, arg, error_log, allow_defined=True)
 
    # Handle the Python signature.
    _resolve_py_signature_types(spec, mod, scope, overload, error_log)

    result = overload.py_signature.result

    # These slots must return Py_ssize_t.
    if is_ssize_return_slot(overload.common.py_slot):
        if spec.abi_version >= (13, 0):
            required_types = (ArgumentType.SSIZE, )
        else:
            required_types = (ArgumentType.SSIZE, ArgumentType.INT)

        if result.type not in required_types or len(result.derefs) != 0 or result.is_reference or result.is_const:
            _log_overload_error(error_log, "must return a Py_ssize_t",
                    overload, scope=scope)

    # These slots must return int.
    if is_int_return_slot(overload.common.py_slot):
        if result.type is not ArgumentType.INT or len(result.derefs) != 0 or result.is_reference or result.is_const:
            _log_overload_error(error_log, "must return an int", overload,
                    scope=scope)

    # These slots must return void.
    if is_void_return_slot(overload.common.py_slot):
        if result.type is not ArgumentType.VOID or len(result.derefs) != 0 or result.is_reference or result.is_const:
            _log_overload_error(error_log, "must return void", overload,
                    scope=scope)

    # These slots must return Py_hash_t.
    if is_hash_return_slot(overload.common.py_slot):
        if spec.abi_version >= (13, 0):
            required_type = ArgumentType.HASH
            required_type_name = 'Py_hash_t'
        else:
            required_type = ArgumentType.LONG
            required_type_name = 'long'

        if result.type is not required_type or len(result.derefs) != 0 or result.is_reference or result.is_const:
            _log_overload_error(error_log,
                    "must return a {0}".format(required_type_name), overload,
                    scope=scope)


def _resolve_py_signature_types(spec, mod, scope, overload, error_log):
    """ Resolve the types of a Python signature. """

    result = overload.py_signature.result

    if result.type is not ArgumentType.VOID or len(result.derefs) != 0:
        if overload.pyqt_method_specifier is PyQtMethodSpecifier.SIGNAL:
            _log_overload_error(error_log, "is a signal and must return void",
                    overload, scope=scope)

        _resolve_type(spec, mod, scope, result, error_log)

        # Results must be simple.
        if result.type is ArgumentType.NONE:
            # The error has already been logged.
            pass
        elif not _supported_type(scope, overload, result, error_log):
            if overload.cpp_signature is overload.py_signature or overload.method_code is None:
                _log_overload_error(error_log,
                        "has an unsupported return type - provide %MethodCode and a {0} signature".format(
                                'C' if spec.c_bindings else 'C++'),
                        overload, scope=scope)

    for arg_nr, arg in enumerate(overload.py_signature.args):
        _resolve_type(spec, mod, scope, arg, error_log)

        if arg.type is ArgumentType.NONE:
            continue

        # Note signal arguments are restricted in their types because we don't
        # (yet) support handwritten code for them.
        if overload.pyqt_method_specifier is PyQtMethodSpecifier.SIGNAL:
            if not _supported_type(scope, overload, arg, error_log):
                _log_overload_error(error_log,
                        "argument {0} has an unsupported type for a Python signature".format(
                                arg_nr + 1),
                        overload, scope=scope)

        elif not _supported_type(scope, overload, arg, error_log, outputs=True):
            if overload.is_virtual:
                _log_overload_error(error_log,
                        "argument {0} has an unsupported type for a Python signature - provide a valid type, %MethodCode, %VirtualCatcherCode and a C++ signature".format(
                                arg_nr + 1),
                        overload, scope=scope)

            _log_overload_error(error_log,
                    "argument {0} has an unsupported type for a Python signature - provide a valid type, %MethodCode and a C++ signature".format(
                            arg_nr + 1),
                    overload, scope=scope)

        if scope is not None:
            _scope_default_value(spec, scope, arg)


# Various type classifications.
_CLASS_TYPES = (ArgumentType.MAPPED, ArgumentType.CLASS)
_SIMPLE_TYPES = (ArgumentType.CFLOAT, ArgumentType.FLOAT, ArgumentType.CDOUBLE,
        ArgumentType.DOUBLE, ArgumentType.ENUM, ArgumentType.BOOL,
        ArgumentType.CBOOL, ArgumentType.BYTE, ArgumentType.SBYTE,
        ArgumentType.UBYTE, ArgumentType.USHORT, ArgumentType.SHORT,
        ArgumentType.UINT, ArgumentType.CINT, ArgumentType.INT,
        ArgumentType.ULONG, ArgumentType.LONG, ArgumentType.ULONGLONG,
        ArgumentType.LONGLONG, ArgumentType.HASH, ArgumentType.SSIZE,
        ArgumentType.SIZE, ArgumentType.PYOBJECT, ArgumentType.PYTUPLE,
        ArgumentType.PYLIST, ArgumentType.PYDICT, ArgumentType.PYCALLABLE,
        ArgumentType.PYSLICE, ArgumentType.PYTYPE, ArgumentType.PYBUFFER,
        ArgumentType.PYENUM, ArgumentType.CAPSULE)
_POINTER_TYPES = (ArgumentType.STRUCT, ArgumentType.UNION, ArgumentType.VOID)
_STRING_TYPES = (ArgumentType.ASCII_STRING, ArgumentType.LATIN1_STRING,
        ArgumentType.UTF8_STRING, ArgumentType.SSTRING, ArgumentType.USTRING,
        ArgumentType.STRING, ArgumentType.WSTRING)


def _resolve_variable_type(spec, variable, error_log):
    """ Resolve the type of a variable. """

    bad_type = True
    variable_type = variable.type

    _resolve_type(spec, variable.module, variable.scope, variable_type,
            error_log)

    if variable_type is ArgumentType.NONE:
        return

    if variable_type.type in _CLASS_TYPES:
        # Class, Class & and Class * are supported.
        if len(variable_type.derefs) <= 1:
            bad_type = False

    elif variable_type.type in _STRING_TYPES:
        # (signed/unsigned) char, (signed/unsigned) char *, wchar_t, wchar_t *
        # are supported.
        if not variable_type.is_reference and len(variable_type.derefs) <= 1:
            bad_type = False

    elif variable_type.type in _SIMPLE_TYPES:
        # These are supported without pointers or references.
        if not variable_type.is_reference and len(variable_type.derefs) == 0:
            bad_type = False

    elif variable_type.type in _POINTER_TYPES:
        # A simple pointer is supported.
        if not variable_type.is_reference and len(variable_type.derefs) == 1:
            bad_type = False

    if bad_type and (variable.get_code is None or (not variable.no_setter and variable.set_code is None)):
        if variable.no_setter:
            set_s = ''
        else:
            set_s = " and %SetCode"

        error_log.log(
                "'{0}' has an unsupported type - provide %GetCode{1}".format(
                    variable.fq_cpp_name, set_s))
 
    if variable_type.type is not ArgumentType.CLASS and variable.access_code is not None:
        error_log.log(
                "'{0}' has %AccessCode but isn't a class instance".format(
                    variable.fq_cpp_name))

    if variable.scope is not None:
        _iface_file_is_used(variable.scope.iface_file.used, variable_type)
    else:
        _iface_file_is_used(variable.module.used, variable_type)

    # Scoped variables need a handler unless they have %AccessCode.
    if variable.access_code is None:
        if variable.scope is not None and not variable.scope.is_hidden_namespace:
            variable.needs_handler = True
            variable.scope.has_variable_handlers = True


def _supported_type(klass, overload, arg, error_log, outputs=False):
    """ See if a type is supported by the generated code. """

    if arg.type in _CLASS_TYPES:
        if arg.is_reference:
            if len(arg.derefs) == 0:
                _default_input(arg)
                return True

            if len(arg.derefs) == 1 and outputs:
                _default_output(arg)
                return True

        elif len(arg.derefs) == 0:
            _ensure_input(klass, overload, arg, error_log)
            return True

        elif len(arg.derefs) == 1:
            if outputs:
                _default_input(arg)
            else:
                _ensure_input(klass, overload, arg, error_log)

            return True

        elif len(arg.derefs) == 2 and outputs:
            _default_output(arg)
            return True

    elif arg.type in _STRING_TYPES:
        if arg.is_reference:
            if len(arg.derefs) <= 1 and outputs:
                _default_output(arg)
                return True

        elif len(arg.derefs) == 0:
            _ensure_input(klass, overload, arg, error_log)
            return True

        elif len(arg.derefs) == 1:
            if outputs:
                _default_input(arg)
            else:
                _ensure_input(klass, overload, arg, error_log)

            return True

        elif len(arg.derefs) == 2 and outputs:
            _default_output(arg)
            return True

    elif arg.type in _SIMPLE_TYPES:
        if arg.is_reference:
            if arg.is_const:
                _ensure_input(klass, overload, arg, error_log)
                return True

            if len(arg.derefs) == 0 and outputs:
                _default_output(arg)
                return True

        elif len(arg.derefs) == 0:
            _ensure_input(klass, overload, arg, error_log)
            return True

        elif len(arg.derefs) == 1 and outputs:
            _default_output(arg)
            return True

    elif arg.type in _POINTER_TYPES:
        if arg.is_reference:
            if len(arg.derefs) == 1 and outputs:
                _default_output(arg)
                return True

        elif len(arg.derefs) == 1:
            _ensure_input(klass, overload, arg, error_log)
            return True

        elif len(arg.derefs) == 2 and outputs:
            _default_output(arg)
            return True

    elif arg.type is ArgumentType.ELLIPSIS:
        # This can only appear in argument lists without * or &.
        _ensure_input(klass, overload, arg, error_log)
        return True

    # Unsupported if we got this far.
    return False


def _ensure_input(klass, overload, arg, error_log):
    """ Ensure the direction of an argument is an input. """

    if arg.is_out:
        _log_overload_error(error_log,
                "has an invalid argument type for /Out/", overload,
                scope=klass)

    arg.is_in = True


def _default_input(arg):
    """ Default the direction of an argument to an input. """

    if not arg.is_out:
        arg.is_in = True


def _default_output(arg):
    """ Default the direction of an argument to an output unless the argument
    is const.
    """

    if not arg.is_out and not arg.is_in:
        if arg.is_const:
            arg.is_in = True
        else:
            arg.is_out = True


def _same_cpp_overload(spec, overload1, overload2):
    """ Compare two overloads and return True if they are the same. """

    # They must both be const, or both not.
    if overload1.is_const is not overload2.is_const:
        return False

    return same_signature(spec, overload1.cpp_signature,
            overload2.cpp_signature)


def _same_python_signature(spec, signature1, signature2):
    """ See if two Python signatures are the same as far as Python is
    concerned.  Return None if any argument's type is unknown.
    """

    a1 = a2 = -1

    while True:
        a1, arg1 = _next_significant_arg(signature1, a1)
        a2, arg2 = _next_significant_arg(signature2, a2)

        if arg1 is not None and arg1.type is ArgumentType.NONE:
            return None

        if arg2 is not None and arg2.type is ArgumentType.NONE:
            return None

        if a1 < 0 or a2 < 0:
            break

        if not same_argument_type(spec, arg1, arg2, strict=False):
            return False

    return a1 < 0 and a2 < 0


def _next_significant_arg(signature, a):
    """ Return the next significant argument from a Python signature (ie. one
    that is not optional or an output only argument.  Return -1 if there isn't
    one.
    """

    a += 1

    while a < len(signature.args):
        arg = signature.args[a]

        if arg.default_value is not None:
            break

        if arg.type is ArgumentType.NONE or arg.is_in:
            return a, arg

        a += 1

    return -1, None


def _scope_default_value(spec, klass, arg):
    """ Add an explicit scope to the default value of an argument if possible.
    """

    # We do a quick check to see if we need to do anything.  This means we can
    # limit the times we need to copy the default value.  It needs to be copied
    # because it will be shared by class versions that have been created on the
    # fly and it may need to be scoped differently for each of those versions.
    # Note that, as we no longer support API versions, this comment may no
    # longer apply.
    if arg.default_value is None:
        return

    for value in arg.default_value:
        if value.value_type is ValueType.SCOPED and value.value.is_simple:
            break
    else:
        return

    # It's not certain that we will do anything, but we assume we will and
    # start copying.
    new_default_value = []

    for value in arg.default_value:
        # Make the copy.
        new_value = copy(value)
        new_default_value.append(new_value)

        # Skip this part of the expression if it isn't a named value or it
        # already has a scope.
        if value.value_type is not ValueType.SCOPED or not value.value.is_simple:
            continue

        # Search the class hierarchy for an enum member with the same name.  If
        # we don't find one, leave it as it is (the compiler will find out if
        # this is a problem).
        original_name = value.value.base_name

        for mro_klass in klass.mro:
            for enum in spec.enums:
                if enum.scope is not mro_klass:
                    continue

                for enum_member in enum.members:
                    if enum_member.cpp_name == original_name:
                        # Take the scope from the class that the enum was
                        # defined in.
                        new_name = ScopedName(mro_klass.iface_file.fq_cpp_name)
                        new_name.append(original_name)
                        new_value.value = new_name

                        # Nothing more to do.
                        break
                else:
                    continue

                break
            else:
                continue

            break

    arg.default_value = new_default_value


def _resolve_type(spec, mod, scope, type, error_log, allow_defined=False):
    """ Resolve a type if possible. """

    # Loop until we've got to a base type.
    while type.type is ArgumentType.DEFINED:
        scoped_name = type.definition
        type.type = ArgumentType.NONE

        # Search the local scopes unless what we are looking for has an
        # explicit global scope.
        if not scoped_name.is_absolute:
            klass_scope = scope
            while klass_scope is not None:
                if klass_scope.iface_file.type is IfaceFileType.CLASS:
                    _search_class_scope(spec, klass_scope, scoped_name, type)
                else:
                    _search_scope(spec, klass_scope, scoped_name, type)

                if type.type is not ArgumentType.NONE:
                    break

                klass_scope = klass_scope.scope

            if type.type is not ArgumentType.NONE:
                break

            # We now need an absolute name.
            scoped_name = scoped_name.absolute

        _name_lookup(spec, mod, scoped_name, type)

        if type.type is ArgumentType.NONE:
            if allow_defined:
                type.type = ArgumentType.DEFINED
            else:
                error_log.log("'{0}' is undefined".format(scoped_name),
                        source_location=type.source_location)

            return

    # See if the type refers to an instantiated template.
    _resolve_instantiated_class_template(spec, type)

    # Replace the base type if it has been mapped.
    if type.type in (ArgumentType.STRUCT, ArgumentType.UNION, ArgumentType.TEMPLATE):
        _search_mapped_types(spec, mod, type)

        # If we still have a template then see if we need to automatically
        # instantiate it.
        if type.type is ArgumentType.TEMPLATE:
            for mapped_type_template in spec.mapped_type_templates:
                if mapped_type_template.mapped_type.type.definition.cpp_name == type.definition.cpp_name and same_template_signature(mapped_type_template.mapped_type.type.definition.types, type.definition.types, deep=True):
                    _instantiate_mapped_type_template(spec, mod,
                            mapped_type_template, type, error_log)
                    break

    # If we are in the main module then mark any generated types as being
    # needed.
    if mod is spec.module:
        _set_needed_type(type)


def _set_needed_type(arg):
    """ Specify that a generated type is needed. """

    if arg.type in (ArgumentType.CLASS, ArgumentType.MAPPED):
        arg.definition.iface_file.needed = True
    elif arg.type is ArgumentType.ENUM:
        arg.definition.needed = True


def _set_needed_exceptions(spec, mod, throw_args):
    """ Specify that a set of thrown arguments are needed. """

    if mod is spec.module and throw_args is not None and throw_args.arguments is not None:
        for exception in throw_args.arguments:
            _set_needs_exception(exception)


def _set_needs_exception(exception):
    """ Specify that an exception is needed. """

    if exception.class_exception is not None:
        exception.class_exception.iface_file.needed = True
    else:
        exception.needed = True


def _resolve_instantiated_class_template(spec, type):
    """ If the type corresponds to a previously instantiated class template
    then replace it with the class that was created.
    """

    if type.type is not ArgumentType.TEMPLATE:
        return

    template = type.definition
    template_signature = template.types

    for arg in template_signature.args:
        _resolve_instantiated_class_template(spec, arg)

    for klass in spec.classes:
        if klass.template is not None and klass.template.cpp_name == template.cpp_name and same_signature(spec, klass.template.types, template_signature):
            type.type = ArgumentType.CLASS
            type.definition = klass
            break


def _instantiate_mapped_type_template(spec, mod, mapped_type_template, type,
        error_log):
    """ Instantiate a mapped type template. """

    expansions = template_expansions(
            mapped_type_template.mapped_type.type.definition.types,
            type.definition.types,
            declared_names=mapped_type_template.signature)

    iface_file = find_iface_file(spec, mod,
            encoded_template_name(type.definition), IfaceFileType.MAPPED_TYPE,
            error_log.log, cpp_type=type)
    iface_file.module = mod

    mapped_type = MappedType(iface_file, copy(type))
    mapped_type.type.derefs = []
    mapped_type.type.is_const = False
    mapped_type.type.is_reference = False
    mapped_type.cpp_name = cached_name(spec, argument_as_str(mapped_type.type))

    if mod is spec.module:
        mapped_type.cpp_name.used = True

    proto_mapped_type = mapped_type_template.mapped_type

    mapped_type.handles_none = proto_mapped_type.handles_none
    mapped_type.needs_user_state = proto_mapped_type.needs_user_state
    mapped_type.no_assignment_operator = proto_mapped_type.no_assignment_operator
    mapped_type.no_copy_ctor = proto_mapped_type.no_copy_ctor
    mapped_type.no_default_ctor = proto_mapped_type.no_default_ctor
    mapped_type.no_release = proto_mapped_type.no_release
    mapped_type.pyqt_flags = proto_mapped_type.pyqt_flags

    if proto_mapped_type.type_hints is not None:
        mapped_type.type_hints = instantiate_type_hints(spec,
                proto_mapped_type.type_hints, expansions)

    used = mapped_type.iface_file.used

    mapped_type.iface_file.type_header_code = template_code_blocks(spec, used,
            proto_mapped_type.iface_file.type_header_code, expansions)

    if proto_mapped_type.convert_from_type_code is not None:
        mapped_type.convert_from_type_code = template_code(spec, used,
                proto_mapped_type.convert_from_type_code, expansions)

    if proto_mapped_type.convert_to_type_code is not None:
        mapped_type.convert_to_type_code = template_code(spec, used,
                proto_mapped_type.convert_to_type_code, expansions)

    if proto_mapped_type.release_code is not None:
        mapped_type.release_code = template_code(spec, used,
                proto_mapped_type.release_code, expansions)

    spec.mapped_types.insert(0, mapped_type)

    _replace_template_type(mapped_type, type)


def _search_class_scope(spec, scope, scoped_name, type):
    """ Search for a name in a class scope and resolve the corresponding type.
    """

    for mro_klass in scope.mro:
        _search_scope(spec, mro_klass, scoped_name, type)

        if type.type is not ArgumentType.NONE:
            break


def _search_scope(spec, scope, scoped_name, type):
    """ Search for a name in a scope and resolve the corresponding type. """

    # Prepend the scope and see if it exists.
    fq_name = ScopedName(scoped_name)
    fq_name.prepend(scope.iface_file.fq_cpp_name)

    _name_lookup(spec, scope.iface_file.module, fq_name, type)


def _name_lookup(spec, mod, scoped_name, type):
    """ Look up a name and resole the corresponding type. """

    _search_mapped_types(spec, mod, type, scoped_name)
    if type.type is not ArgumentType.NONE:
        return

    search_typedefs(spec, scoped_name, type)
    if type.type is not ArgumentType.NONE:
        return

    _search_enums(spec, scoped_name, type)
    if type.type is not ArgumentType.NONE:
        return

    _search_classes(spec, mod, scoped_name, type)


def _search_mapped_types(spec, mod, type, scoped_name=None):
    """ Search the mapped types for a name and resolve the type. """

    # Patch back to defined types so we can use same_base_type().
    if scoped_name is not None:
        orig_name = type.definition
        type.definition = scoped_name
        type.type = ArgumentType.DEFINED

    for mapped_type in spec.mapped_types:
        if same_base_type(mapped_type.type, type):
            break
    else:
        # Restore because we didn't find anything.
        if scoped_name is not None:
            type.definition = orig_name
            type.type = ArgumentType.NONE

        return

    _replace_template_type(mapped_type, type)


def _replace_template_type(mapped_type, type):
    """ If a mapped type is based on a template then update the type with a
    copy that keeps the original types of the template arguments.
    """

    if mapped_type.type.type is ArgumentType.TEMPLATE:
        dst_types = None

        for arg_nr, arg in enumerate(type.definition.types.args):
            if arg.original_typedef is None:
                continue

            # Create an appropriately deep copy now that we know it is needed
            # and if it hasn't already been done.
            if dst_types is None:
                mapped_type = copy(mapped_type)
                mapped_type.type = copy(mapped_type.type)
                mapped_type.type.definition = copy(mapped_type.type.definition)
                mapped_type.type.definition.types = copy(mapped_type.type.definition.types)
                mapped_type.type.definition.types.args = [copy(arg) for arg in
                        mapped_type.type.definition.types.args]
                dst_types = mapped_type.type.definition.types

            dst_types.args[arg_nr].original_typedef = arg.original_typedef

    # Replace the template with the mapped type.
    type.type = ArgumentType.MAPPED
    type.definition = mapped_type
    type.type_hints = mapped_type.type_hints


def _merged_type_hints(type_hints, defaults):
    """ Return a TypeHints object that is the merge of another and some
    defaults.
    """

    if type_hints is None:
        return defaults

    if defaults is None:
        return type_hints

    if type_hints.hint_in is None:
        type_hints.hint_in = defaults.hint_in

    if type_hints.hint_out is None:
        type_hints.hint_out = defaults.hint_out

    if type_hints.default_value is None:
        type_hints.default_value = defaults.default_value

    return type_hints


def _search_enums(spec, scoped_name, type):
    """ Search the enums for a name and resolve the type. """

    for enum in spec.enums:
        if enum.fq_cpp_name is None:
            continue

        if enum.fq_cpp_name == scoped_name:
            type.type = ArgumentType.ENUM
            type.definition = enum
            break


def _search_classes(spec, mod, scoped_name, type):
    """ Search the classes for one with a particular name and resolve it as a
    type.
    """

    for klass in spec.classes:
        # Ignore an external class unless it was declared in the same module as
        # the name is being used.
        if klass.external and klass.iface_file.module is not mod:
            continue

        if klass.iface_file.fq_cpp_name == scoped_name:
            type.type = ArgumentType.CLASS
            type.definition = klass
            type.type_hints = _merged_type_hints(type.type_hints,
                    klass.type_hints)

            break


def _iface_files_are_used_by_signature(used, signature, need_types=False):
    """ Make sure all interface files for a signature are used. """

    if signature.result is not None:
        _iface_file_is_used(used, signature.result, need_types=need_types)

    for arg in signature.args:
        _iface_file_is_used(used, arg, need_types=need_types)


def _iface_files_are_used_by_overload(spec, used, overload, need_types=False):
    """ Make sure all interface files for an overload are used. """

    _iface_files_are_used_by_signature(used, overload.py_signature,
            need_types=need_types)

    if overload.cpp_signature is not overload.py_signature:
        _iface_files_are_used_by_signature(used, overload.cpp_signature,
                need_types=need_types)

    # Don't bother with %TypeHeaderCode from %Exception for later ABI versions.
    if spec.abi_version >= (13, 1) or (spec.abi_version >= (12, 9) and spec.abi_version < (13, 0)):
        return

    throw_args = overload.throw_args
    if throw_args is not None and throw_args.arguments is not None:
        for exception in throw_args.arguments:
            append_iface_file(used, exception.iface_file)

            if need_types:
                _set_needs_exception(exception)


def _iface_file_is_used(used, arg, need_types=False):
    """ If a type has an interface file then add it to the the given list of
    used interface files so that the header file is #included in the generated
    code.
    """

    if arg.type in (ArgumentType.CLASS, ArgumentType.MAPPED):
        iface_file = arg.definition.iface_file
    elif arg.type is ArgumentType.ENUM:
        iface_file = _get_iface_file_for_enum(arg.definition)
    else:
        iface_file = None

    if iface_file is not None:
        append_iface_file(used, iface_file)

        # For mapped type templates we also need the template arguments.  These
        # will be in the mapped type's used list (which itself will be empty
        # for non-template mapped types).
        if arg.type is ArgumentType.MAPPED:
            for used_iface_file in iface_file.used:
                append_iface_file(used, used_iface_file)

    if need_types:
        _set_needed_type(arg)


def _enum_iface_file_is_used(enum, mod):
    """ Add an enum's interface file to that used by a module. """

    enum_iface_file = _get_iface_file_for_enum(enum)

    if enum_iface_file is not None:
        append_iface_file(mod.used, enum_iface_file)


def _get_iface_file_for_enum(enum):
    """ Return the interface file for an enum, or None if it doesn't have one.
    """

    if enum.fq_cpp_name is not None:
        if enum.scope is not None:
            return enum.scope.iface_file

    return None


def _create_sorted_numbered_types(spec, mod, error_log):
    """ Create the sorted list of numbered types for a module.  For the main
    module this will be every type defined in the module.  For other modules
    this will be every type needed by the main module.
    """

    # Collect the needed types.
    for klass in spec.classes:
        if klass.iface_file.module is not mod:
            continue

        if mod is spec.module or klass.iface_file.needed:
            if not klass.is_hidden_namespace:
                mod.needed_types.append(Argument(ArgumentType.CLASS,
                        definition=klass, name=klass.iface_file.cpp_name))

    for mapped_type in spec.mapped_types:
        if mapped_type.iface_file.module is not mod:
            continue

        if mod is spec.module or mapped_type.iface_file.needed:
            mod.needed_types.append(Argument(ArgumentType.MAPPED,
                    definition=mapped_type, name=mapped_type.cpp_name))

    for enum in spec.enums:
        if enum.module is not mod:
            continue

        if enum.fq_cpp_name is None:
            continue

        if mod is spec.module or enum.needed:
            mod.needed_types.append(Argument(ArgumentType.ENUM,
                    definition=enum, name=enum.cached_fq_cpp_name))

    # Sort the list and assign type numbers.
    mod.needed_types.sort(key=lambda t: t.name.name)

    needed_type_nr = 0

    for needed_type in mod.needed_types:
        if needed_type.type is ArgumentType.CLASS:
            needed_type.definition.iface_file.type_nr = needed_type_nr

            # If we find a class called QObject, assume it's Qt.
            if needed_type.name.name == 'QObject':
                if spec.pyqt_qobject is not None:
                    error_log.log(
                            "class 'QObject' has been defined more than once")

                spec.pyqt_qobject = needed_type.definition

        elif needed_type.type is ArgumentType.MAPPED:
            needed_type.definition.iface_file.type_nr = needed_type_nr

        elif needed_type.type is ArgumentType.ENUM:
            needed_type.definition.type_nr = needed_type_nr

        needed_type_nr += 1


def _check_properties(klass, error_log):
    """ Check that any properties are valid. """

    for prop in klass.properties:
        if find_method(klass, prop.getter) is None:
            error_log.log(
                    "property '{0}.{1}' has no getter '{3}'".format(
                            klass.py_name.name, prop.name.name, prop.getter))

        if prop.setter is not None and find_method(klass, prop.setter) is None:
            error_log.log(
                    "property '{0}.{1}' has no setter '{3}'".format(
                            klass.py_name.name, prop.name.name, prop.setter))


def _log_overload_error(error_log, text, overload, scope=None):
    """ Log an error about an overload. """

    if scope is None:
        fq_cpp_name = overload.cpp_name
    else:
        fq_cpp_name = f'{scope.iface_file.fq_cpp_name}::{overload.cpp_name}'

    error_log.log(f"'{fq_cpp_name}' {text}",
            source_location=overload.source_location)
