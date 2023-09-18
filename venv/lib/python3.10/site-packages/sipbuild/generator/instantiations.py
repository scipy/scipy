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
from .specification import (Argument, ArgumentType, FunctionCall,
        IfaceFileType, KwArgs, Signature, TypeHints, Value, ValueType)
from .templates import (template_code, template_code_blocks,
        template_expansions, template_string)
from .utils import append_iface_file, cached_name, normalised_scoped_name


def instantiate_class(p, symbol, fq_cpp_name, tmpl_names, proto_class,
            template, py_name, no_type_name, docstring, pm):
    """ Instantiate a class template. """

    # Create the expansions.
    expansions = template_expansions(tmpl_names, template.types)

    # Add a mapping from the template name to the instantiated name.
    expansions[str(proto_class.iface_file.fq_cpp_name)] = str(fq_cpp_name)

    # Create the new class starting with a shallow copy.
    i_class = copy(proto_class)

    if docstring is not None:
        i_class.docstring = docstring

    i_class.mro = []
    i_class.virtual_overloads = []
    i_class.visible_members = []
    i_class.py_name = cached_name(pm.spec, py_name)
    i_class.template = template
    i_class.no_type_name = no_type_name

    # Handle the interface file.
    i_class.iface_file = pm.find_iface_file(p, symbol, fq_cpp_name,
            IfaceFileType.CLASS)
    i_class.iface_file.module = pm.module_state.module

    used = i_class.iface_file.used

    i_class.iface_file.type_header_code = template_code_blocks(pm.spec, used,
            proto_class.iface_file.type_header_code, expansions)
    i_class.iface_file.type_header_code.extend(
            proto_class.iface_file.type_header_code)

    # Make a copy of the used list and add the enclosing scope.
    for iface_file in proto_class.iface_file.used:
        append_iface_file(i_class.iface_file.used, iface_file)

    # Include any scope header code.
    i_class.scope = pm.scope

    if i_class.scope is not None:
        i_class.iface_file.type_header_code.extend(
                i_class.scope.iface_file.type_header_code)

    if pm.in_main_module:
        i_class.iface_file.cpp_name.used = True
        i_class.py_name.used = True

    # Handle any type hints.
    if proto_class.type_hints is not None:
        i_class.type_hints = instantiate_type_hints(pm.spec,
                proto_class.type_hints, expansions)

    # Handle any flagged enums.
    if proto_class.pyqt_flags_enums is not None:
        i_class.pyqt_flags_enums = [template_string(s, expansions)
                for s in proto_class.pyqt_flags_enums]

    # Handle the super-classes.
    i_class.superclasses = []

    for superclass in proto_class.superclasses:
        superclass_name = superclass.iface_file.fq_cpp_name

        # Don't try and expand defined or scoped classes.
        if superclass.iface_file.module is not None or not superclass_name.is_simple:
            i_class.superclasses.append(superclass)
            continue

        for a, arg in enumerate(tmpl_names.args):
            if superclass_name.base_name == arg.definition.base_name:
                tmpl_arg = template.types.args[a]

                if tmpl_arg.type is ArgumentType.DEFINED:
                    i_superclass = pm.find_class(p, symbol,
                            IfaceFileType.CLASS, tmpl_arg.definition)
                elif tmpl_arg.type is ArgumentType.CLASS:
                    i_superclass = tmpl_arg.definition
                else:
                    pm.parser_error(p, symbol,
                            "template argument '{0}' must expand to a class".format(superclass_name))
                    i_superclass = superclass

                i_class.superclasses.append(i_superclass)

    # Handle the enums.
    _instantiate_enums(tmpl_names, proto_class, template, i_class, expansions,
            pm)

    # Handle the variables.
    _instantiate_vars(tmpl_names, proto_class, template, i_class, expansions,
            pm)

    # Handle the typedefs.
    _instantiate_typedefs(p, symbol, tmpl_names, proto_class, template,
            i_class, expansions, pm)

    # Handle the ctors.
    i_class.ctors = _instantiate_ctors(tmpl_names, proto_class, template,
            i_class, expansions, pm)

    # Handle the methods.
    i_class.members = _instantiate_methods(proto_class.members,
            i_class.iface_file.module, pm)
    i_class.overloads = _instantiate_overloads(proto_class.overloads,
            proto_class.members, i_class.members, tmpl_names, proto_class,
            template, i_class, expansions, pm)

    # Handle the remaining handwritten code.
    i_class.bi_get_buffer_code = template_code(pm.spec, used,
            proto_class.bi_get_buffer_code, expansions)
    i_class.bi_release_buffer_code = template_code(pm.spec, used,
            proto_class.bi_release_buffer_code, expansions)
    i_class.convert_from_type_code = template_code(pm.spec, used,
            proto_class.convert_from_type_code, expansions)
    i_class.convert_to_subclass_code = template_code(pm.spec, used,
            proto_class.convert_to_subclass_code, expansions)
    i_class.convert_to_type_code = template_code(pm.spec, used,
            proto_class.convert_to_type_code, expansions)
    i_class.dealloc_code = template_code_blocks(pm.spec, used,
            proto_class.dealloc_code, expansions)
    i_class.dtor_virtual_catcher_code = template_code(pm.spec, used,
            proto_class.dtor_virtual_catcher_code, expansions)
    i_class.finalisation_code = template_code(pm.spec, used,
            proto_class.finalisation_code, expansions)
    i_class.gc_clear_code = template_code(pm.spec, used,
            proto_class.gc_clear_code, expansions)
    i_class.gc_traverse_code = template_code(pm.spec, used,
            proto_class.gc_traverse_code, expansions)
    i_class.instance_code = template_code(pm.spec, used,
            proto_class.instance_code, expansions)
    i_class.pickle_code = template_code(pm.spec, used, proto_class.pickle_code,
            expansions)
    i_class.type_code = template_code_blocks(pm.spec, used,
            proto_class.type_code, expansions)
    i_class.type_hint_code = template_code(pm.spec, used,
            proto_class.type_hint_code, expansions)

    pm.spec.classes.insert(0, i_class)


def _instantiate_argument(proto_arg, proto_class, tmpl_names, template,
        i_class, expansions, pm):
    """ Return an instantiated Argument object. """

    # Start with a shallow copy.
    i_arg = copy(proto_arg)

    # Descend into any sub-templates.
    if proto_arg.type is ArgumentType.TEMPLATE:
        proto_template = proto_arg.definition
        i_template = copy(proto_template)
        i_template.types = _instantiate_signature(proto_template.types,
                proto_class, tmpl_names, template, i_class, expansions, pm)
        i_arg.definition = i_template

    # Handle any default value.
    if proto_arg.default_value is not None:
        i_arg.default_value = [_instantiate_value(v, expansions)
                for v in proto_arg.default_value]

    # Handle any type hints.
    if proto_arg.type_hints is not None:
        i_arg.type_hints = instantiate_type_hints(pm.spec,
                proto_arg.type_hints, expansions)

    # Handle arguments that are unscoped names.
    if proto_arg.type is ArgumentType.DEFINED and proto_arg.definition.is_simple:
        name = proto_arg.definition.base_name

        for a, arg in enumerate(tmpl_names.args):
            if name == arg.definition.base_name:
                tad = template.types.args[a]

                i_arg.type = tad.type
                i_arg.definition = tad.definition

                # We take the constrained flag from the real type.
                i_arg.is_constrained = tad.is_constrained

                break
        else:
            # Handle the class name itself.
            if name == proto_class.iface_file.fq_cpp_name.base_name:
                i_arg.type = ArgumentType.CLASS
                i_arg.definition = i_class
                i_arg.original_type = None

    return i_arg


def _instantiate_ctors(tmpl_names, proto_class, template, i_class, expansions,
        pm):
    """ Return a list of the instantiated ctors of a template class. """

    i_ctors = []
    used = i_class.iface_file.used

    for proto_ctor in proto_class.ctors:
        # Start with a shallow copy.
        i_ctor = copy(proto_ctor)

        i_ctor.py_signature = _instantiate_signature(proto_ctor.py_signature,
                proto_class, tmpl_names, template, i_class, expansions, pm,
                kw_args=proto_ctor.kw_args)

        if proto_ctor.cpp_signature is proto_ctor.py_signature:
            i_ctor.cpp_signature = i_ctor.py_signature
        else:
            i_ctor.cpp_signature = _instantiate_signature(
                    proto_ctor.cpp_signature.cpp_signature, proto_class,
                    tmpl_names, template, i_class, expansions, pm)

        i_ctor.method_code = template_code(pm.spec, used,
                proto_ctor.method_code, expansions)
        i_ctor.premethod_code = template_code(pm.spec, used,
                proto_ctor.premethod_code, expansions)

        # Handle the default ctor.
        if proto_class.default_ctor is proto_ctor:
            i_class.default_ctor = i_ctor

        i_ctors.append(i_ctor)

    return i_ctors


def _instantiate_enums(tmpl_names, proto_class, template, i_class, expansions,
        pm):
    """ Instantiate the enums for a template class. """

    for proto_enum in list(pm.spec.enums):
        if proto_enum.scope is not proto_class:
            continue

        # Start with a shallow copy.
        i_enum = copy(proto_enum)

        if proto_enum.fq_cpp_name is not None:
            i_enum.fq_cpp_name = normalised_scoped_name(proto_enum.fq_cpp_name,
                    i_class)
            i_enum.cached_fq_cpp_name = cached_name(pm.spec,
                    str(i_enum.fq_cpp_name))

        if pm.in_main_module:
            if i_enum.py_name is not None:
                i_enum.py_name = True

            if i_enum.cached_fq_cpp_name is not None:
                i_enum.cached_fq_cpp_name.used = True

        i_enum.scope = i_class
        i_enum.module = i_class.iface_file.module
        i_enum.members = []

        for proto_member in proto_enum.members:
            # Start with a shallow copy.
            w_member = copy(proto_member)

            w_member.scope = i_enum

            i_enum.members.append(w_member)

        pm.spec.enums.insert(0, i_enum)


def _instantiate_methods(proto_methods, target_module, pm):
    """ Return a list of the instantiated methods of a template class or enum.
    """

    i_methods = []

    for proto_method in proto_methods:
        # Start with a shallow copy.
        i_method = copy(proto_method)

        i_method.module = target_module

        if pm.in_main_module:
            i_method.py_name.used = True

        i_methods.append(i_method)

    return i_methods


def _instantiate_overloads(proto_overloads, proto_methods, i_methods,
        tmpl_names, proto_class, template, i_class, expansions, pm):
    """ Return a list of the instantiated overloads of a template class or
    enum.
    """

    i_overloads = []
    used = i_class.iface_file.used

    for proto_overload in proto_overloads:
        # Start with a shallow copy.
        i_overload = copy(proto_overload)

        for i_method, proto_method in zip(i_methods, proto_methods):
            if proto_overload.common is proto_method:
                i_overload.common = i_method
                break

        i_overload.py_signature = _instantiate_signature(
                proto_overload.py_signature, proto_class, tmpl_names, template,
                i_class, expansions, pm, kw_args=proto_overload.kw_args)

        if proto_overload.cpp_signature is proto_overload.py_signature:
            i_overload.cpp_signature = i_overload.py_signature
        else:
            i_overload.cpp_signature = _instantiate_signature(
                    proto_overload.cpp_signature.cpp_signature, proto_class,
                    tmpl_names, template, i_class, expansions, pm)

        i_overload.method_code = template_code(pm.spec, used,
                proto_overload.method_code, expansions)
        i_overload.premethod_code = template_code(pm.spec, used,
                proto_overload.premethod_code, expansions)
        i_overload.virtual_call_code = template_code(pm.spec, used,
                proto_overload.virtual_call_code, expansions)
        i_overload.virtual_catcher_code = template_code(pm.spec, used,
                proto_overload.virtual_catcher_code, expansions)

        i_overloads.append(i_overload)

    return i_overloads


def _instantiate_signature(proto_signature, proto_class, tmpl_names, template,
        i_class, expansions, pm, kw_args=KwArgs.NONE):
    """ Return an instantiated Signature object. """

    i_signature = Signature()

    for proto_arg in proto_signature.args:
        i_arg = _instantiate_argument(proto_arg, proto_class, tmpl_names,
                template, i_class, expansions, pm)

        i_signature.args.append(i_arg)

        # Make sure we have the name of any keyword argument.
        if pm.in_main_module and i_arg.name is not None:
            if kw_args is KwArgs.ALL or (kw_args is KwArgs.OPTIONAL and i_arg.default_value is not None):
                i_arg.name.used = True

    if proto_signature.result is not None:
        i_signature.result = _instantiate_argument(proto_signature.result,
                proto_class, tmpl_names, template, i_class, expansions, pm)

    return i_signature


def instantiate_type_hints(spec, proto_type_hints, expansions):
    """ Return an instantiated TypeHints object. """

    if proto_type_hints.hint_in is not None:
        hint_in = template_string(proto_type_hints.hint_in, expansions,
                        scope_replacement='.')
    else:
        hint_in = None

    if proto_type_hints.hint_out is not None:
        hint_out = template_string(proto_type_hints.hint_out, expansions,
                        scope_replacement='.')
    else:
        hint_out = None

    return TypeHints(hint_in=hint_in, hint_out=hint_out,
            default_value=proto_type_hints.default_value)


def _instantiate_typedefs(p, symbol, tmpl_names, proto_class, template,
        i_class, expansions, pm):
    """ Instantiate the typedefs of a template class. """

    for proto_typedef in pm.spec.typedefs:
        if proto_typedef.scope is not proto_class:
            continue

        # Start with a shallow copy.
        i_typedef = copy(proto_typedef)

        i_typedef.fq_cpp_name = normalised_scoped_name(
                proto_typedef.fq_cpp_name, i_class)
        i_typedef.scope = i_class
        i_typedef.module = i_class.iface_file.module

        i_typedef.type = _instantiate_argument(proto_typedef.type, proto_class,
                tmpl_names, template, i_class, expansions, pm)

        pm.add_typedef(p, symbol, i_typedef)


def _instantiate_value(proto_value, expansions):
    """ Return an instantiated Value object. """

    # We only handle the subset where the value is an function call, ie. a
    # template ctor.
    i_value = proto_value

    if proto_value.value_type is ValueType.FCALL and proto_value.value.result.type is ArgumentType.DEFINED:
        proto_name = proto_value.value.result.definition

        if proto_name.is_simple:
            i_name = ScopedName.parse(
                    template_string(proto_name.base_name, expansions))
            i_result = Argument(type=ArgumentType.DEFINED, definition=i_name)
            i_fcall = FunctionCall(result=i_result,
                    args=proto_value.value.args)
            i_value = Value(ValueType.FCALL, value=i_fcall)

    return i_value


def _instantiate_vars(tmpl_names, proto_class, template, i_class, expansions,
        pm):
    """ Instantiate the enums for a template class. """

    used = i_class.iface_file.used

    for proto_var in pm.spec.variables:
        if proto_var.scope is not proto_class:
            continue

        # Start with a shallow copy.
        i_var = copy(proto_var)

        if pm.in_main_module:
            i_var.py_name.used = True

        i_var.fq_cpp_name = normalised_scoped_name(proto_var.fq_cpp_name,
                i_class)
        i_var.scope = i_class
        i_var.module = i_class.iface_file.module

        i_var.type = _instantiate_argument(proto_var.type, proto_class,
                tmpl_names, template, i_class, expansions, pm)

        i_var.access_code = template_code(pm.spec, used, proto_var.access_code,
                expansions)
        i_var.get_code = template_code(pm.spec, used, proto_var.get_code,
                expansions)
        i_var.set_code = template_code(pm.spec, used, proto_var.set_code,
                expansions)

        pm.spec.variables.append(i_var)
