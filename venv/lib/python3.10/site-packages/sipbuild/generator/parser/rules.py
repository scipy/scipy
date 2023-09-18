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


from ..scoped_name import ScopedName
from ..specification import (AccessSpecifier, Argument, ArgumentType,
        ArrayArgument, ClassKey, Docstring, DocstringFormat, Extract,
        FunctionCall, IfaceFile, IfaceFileType, KwArgs, License, MappedType,
        MappedTypeTemplate, Overload, Property, PyQtMethodSpecifier,
        QualifierType, Signature, Template, ThrowArguments, Value, ValueType,
        VirtualErrorHandler, WrappedTypedef, WrappedVariable)
from ..templates import same_template_signature
from ..utils import cached_name, normalised_scoped_name, search_typedefs

from .annotations import DottedName
from .tokens import tokens


# This is set externally by the parser manager.
parser = None

# The start symbol as we prefer to keep productions in alphabetical order.
start = 'specification'


def p_error(t):
    """ Invoked when a syntax error occurs. """

    if t is None:
        from .parser_manager import UnexpectedEOF

        raise UnexpectedEOF()

    t.lexer.pm.lexer_error(t, "'{0}' is unexpected".format(t.value))

    # Skip past the next '}', ';' or '%End' and restart.
    t.lexer.pm.set_lexer_state()

    t = parser.token()
    while t is not None and t.type not in ('}', ';', 'End'):
        t = parser.token()

    parser.restart()


# Top-level productions. ######################################################

def p_specification(p):
    """specification : statement
        | specification statement"""


def p_statement(p):
    """statement : eof
        | namespace_statement
        | composite_module
        | copying
        | defdocstringfmt
        | defdocstringsig
        | defencoding
        | defmetatype
        | defsupertype
        | exported_header_code
        | exported_type_hint_code
        | extract
        | feature
        | hidden_ns
        | import
        | include
        | init_code
        | license
        | mapped_type
        | mapped_type_template
        | module
        | module_code
        | module_header_code
        | platforms
        | plugin
        | preinit_code
        | postinit_code
        | timeline
        | type_hint_code
        | unit_code
        | unit_postinclude_code
        | virtual_error_handler"""


def p_namespace_statement(p):
    """namespace_statement : if_start
        | if_end
        | class_decl
        | class_template
        | enum_decl
        | exception
        | function
        | namespace_decl
        | struct_decl
        | typedef_decl
        | union_decl
        | variable
        | type_header_code"""


def p_eof(p):
    "eof : EOF"

    p.parser.pm.pop_module_state()


# State changing productions. #################################################

def p_begin_args(p):
    "begin_args :"

    p.parser.pm.set_lexer_state('directive')


def p_end_args(p):
    "end_args :"

    p.parser.pm.set_lexer_state()


def p_need_eol(p):
    "need_eol :"

    p.parser.pm.set_lexer_state('needeol')


# %AutoPyName #################################################################

def p_autopyname(p):
    "autopyname : AutoPyName begin_args '(' remove_leading '=' STRING end_args ')'"

    if p.parser.pm.skipping:
        return

    p[0] = ('REMOVE_LEADING', p[6])


# %BIGetBufferCode ############################################################

def p_get_buffer_code(p):
    "get_buffer_code : BIGetBufferCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.bi_get_buffer_code is not None:
        pm.parser_error(p, 1, "%BIGetBufferCode can only be specified once")

    pm.scope.bi_get_buffer_code = p[2]


# %BIReleaseBufferCode ########################################################

def p_release_buffer_code(p):
    "release_buffer_code : BIReleaseBufferCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.bi_release_buffer_code is not None:
        pm.parser_error(p, 1,
                "%BIReleaseBufferCode can only be specified once")

    pm.scope.bi_release_buffer_code = p[2]


# %CompositeModule ############################################################

def p_composite_module(p):
    """composite_module : CompositeModule dotted_name c_module_body
        | CompositeModule begin_args '(' c_module_args end_args ')' c_module_body"""

    pm = p.parser.pm

    if pm.skipping:
        return

    # A composite module must be the first one in the specification.
    if not pm.in_main_module:
        pm.parser_error(p, 1, "a %CompositeModule cannot be %Imported")

    pm.spec.is_composite = True

    if len(p) == 4:
        name = p[2]
        body = p[3]
    else:
        name = p[4]['name']
        body = p[7]

    module = pm.module_state.module

    module.fq_py_name = cached_name(pm.spec, str(name))

    for directive in body:
        if isinstance(directive, Docstring):
            module.docstring = directive


def p_c_module_args(p):
    """c_module_args : c_module_arg
        | c_module_args ',' c_module_arg"""

    if len(p) == 4:
        p[1].update(p[3])

    p[0] = p[1]


def p_c_module_arg(p):
    "c_module_arg : name '=' dotted_name"

    p[0] = {p[1]: p[3]}


def p_c_module_body(p):
    """c_module_body : '{' c_module_body_directives '}' ';'
        | empty"""

    p[0] = p[2] if len(p) == 4 else []


def p_c_module_body_directives(p):
    """c_module_body_directives : c_module_body_directive
        | c_module_body_directives c_module_body_directive"""

    if len(p) == 2:
        value = []
        directive = p[1]
    else:
        value = p[1]
        directive = p[2]

    if directive is not None:
        value.append(directive)

    p[0] = value


def p_c_module_body_directive(p):
    """c_module_body_directive : if_start
        | if_end
        | docstring"""

    p[0] = p[1]


# %ConvertFromTypeCode ########################################################

def p_convert_from_type_code(p):
    "convert_from_type_code : ConvertFromTypeCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.convert_from_type_code is not None:
        pm.parser_error(p, 1,
                "%ConvertFromTypeCode can only be specified once")

    pm.scope.convert_from_type_code = p[2]


# %ConvertToSubClassCode ######################################################

def p_convert_to_subclass_code(p):
    "convert_to_subclass_code : ConvertToSubClassCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.convert_to_subclass_code is not None:
        pm.parser_error(p, 1,
                "%ConvertToSubClassCode can only be specified once")

    pm.scope.convert_to_subclass_code = p[2]


# %ConvertToTypeCode ##########################################################

def p_convert_to_type_code(p):
    "convert_to_type_code : ConvertToTypeCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.convert_to_type_code is not None:
        pm.parser_error(p, 1, "%ConvertToTypeCode can only be specified once")

    pm.scope.convert_to_type_code = p[2]


# %Copying ####################################################################

def p_copying(p):
    "copying : Copying CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.module_state.module.copying.append(p[2])


# %DefaultDocstringFormat #####################################################

def p_defdocstringfmt(p):
    """defdocstringfmt : DefaultDocstringFormat STRING
        | DefaultDocstringFormat begin_args '(' name '=' STRING end_args ')'"""

    pm = p.parser.pm

    if pm.skipping:
        return

    symbol = 2 if len(p) == 3 else 6

    pm.module_state.module.default_docstring_format = pm.convert_docstring_format(p, symbol)


# %DefaultDocstringSignature ##################################################

def p_defdocstringsig(p):
    """defdocstringsig : DefaultDocstringSignature STRING
        | DefaultDocstringSignature begin_args '(' name '=' STRING end_args ')'"""

    pm = p.parser.pm

    if pm.skipping:
        return

    symbol = 2 if len(p) == 3 else 6

    pm.module_state.module.default_docstring_signature = pm.convert_docstring_signature(p, symbol)


# %DefaultEncoding ############################################################

def p_defencoding(p):
    """defencoding : DefaultEncoding STRING
        |  DefaultEncoding begin_args '(' name '=' STRING end_args ')'"""

    pm = p.parser.pm

    if pm.skipping:
        return

    symbol = 2 if len(p) == 3 else 6

    pm.module_state.default_encoding = pm.convert_encoding(p, symbol)


# %DefaultMetatype ############################################################

def p_defmetatype(p):
    """defmetatype : DefaultMetatype dotted_name
        | DefaultMetatype begin_args '(' name '=' dotted_name end_args ')'"""

    pm = p.parser.pm

    if pm.skipping:
        return

    module = pm.module_state.module

    if module.default_metatype:
        pm.parser_error(p, 1,
                "%DefaultMetatype has already been defined for this module")

    symbol = 2 if len(p) == 3 else 6

    module.default_metatype = cached_name(pm.spec, str(p[symbol]))


# %DefaultSupertype ###########################################################

def p_defsupertype(p):
    """defsupertype : DefaultSupertype dotted_name
        | DefaultSupertype begin_args '(' name '=' dotted_name end_args ')'"""

    pm = p.parser.pm

    if pm.skipping:
        return

    module = pm.module_state.module

    if module.default_supertype:
        pm.parser_error(p, 1,
                "%DefaultSupertype has already been defined for this module")

    symbol = 2 if len(p) == 3 else 6

    module.default_supertype = cached_name(pm.spec, str(p[symbol]))


# %Docstring ##################################################################

def p_docstring(p):
    "docstring : Docstring docstring_args CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    module = pm.module_state.module

    ds_format = p[2].get('format', module.default_docstring_format)

    if ds_format == DocstringFormat.DEINDENTED:
        from textwrap import dedent

        text = dedent(p[3].text)
    else:
        text = p[3].text

    p[0] = Docstring(p[2].get('signature', module.default_docstring_signature),
            text)


def p_docstring_args(p):
    """docstring_args : empty
    | STRING
    | begin_args '(' docstring_arg_list end_args ')'"""

    if len(p) == 2:
        pm = p.parser.pm

        if p[1] is None:
            value = {}
        else:
            value = {'format': pm.convert_docstring_format(p, 1)}
    else:
        value = p[3]

    p[0] = value


def p_docstring_arg_list(p):
    """docstring_arg_list : docstring_arg
        | docstring_arg_list ',' docstring_arg"""

    if len(p) == 4:
        p[1].update(p[3])

    p[0] = p[1]


def p_docstring_arg(p):
    """docstring_arg : format '=' STRING
    | signature '=' STRING"""

    pm = p.parser.pm

    if p[1] == 'format':
        value = pm.convert_docstring_format(p, 3)
    else:
        value = pm.convert_docstring_signature(p, 3)

    p[0] = {p[1]: value}


# %ExportedHeaderCode #########################################################

def p_exported_header_code(p):
    "exported_header_code : ExportedHeaderCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.spec.exported_header_code.append(p[2])


# %ExportedTypeHintCode #######################################################

def p_exported_type_hint_code(p):
    "exported_type_hint_code : ExportedTypeHintCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if not pm.in_main_module:
        pm.spec.exported_type_hint_code.append(p[2])


# %Extract ####################################################################

def p_extract(p):
    """extract : Extract NAME CODE_BLOCK
        | Extract begin_args '(' extract_args end_args ')' CODE_BLOCK"""

    pm = p.parser.pm

    if pm.skipping:
        return

    if len(p) == 4:
        id = p[2]
        order = -1
        part = p[3]
    else:
        args = p[4]

        try:
            id = args['id']
        except KeyError:
            pm.parser_error(p, 1, "%Extract must specify 'id'")
            id = ''

        order = args.get('order', -1)
        part = p[7]

    pm.spec.extracts.append(Extract(id, order, part.text))


def p_extract_args(p):
    """extract_args : extract_arg
        | extract_args ',' extract_arg"""

    if len(p) == 4:
        p[1].update(p[3])

    p[0] = p[1]


def p_extract_arg(p):
    """extract_arg : id '=' NAME
        | order '=' NUMBER"""

    p[0] = {p[1]: p[3]}


# %Feature ####################################################################

def p_feature(p):
    """feature : Feature NAME
        | Feature begin_args '(' name '=' NAME end_args ')'"""

    symbol = 2 if len(p) == 3 else 6

    # Note that qualifiers are never skipped.
    p.parser.pm.add_qualifier(p, symbol, p[symbol], QualifierType.FEATURE)


# %FinalisationCode ###########################################################

def p_finalisation_code(p):
    "finalisation_code : FinalisationCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.finalisation_code is not None:
        pm.parser_error(p, 1, "%FinalisationCode can only be specified once")

    pm.scope.finalisation_code = p[2]


# %GCClearCode ################################################################

def p_gc_clear_code(p):
    "gc_clear_code : GCClearCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.gc_clear_code is not None:
        pm.parser_error(p, 1, "%GCClearCode can only be specified once")

    pm.scope.gc_clear_code = p[2]


# %GCTraverseCode #############################################################

def p_gc_traverse_code(p):
    "gc_traverse_code : GCTraverseCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.gc_traverse_code is not None:
        pm.parser_error(p, 1, "%GCTraverseCode can only be specified once")

    pm.scope.gc_traverse_code = p[2]


# %HiddenNamespace ############################################################

def p_hidden_ns(p):
    """hidden_ns : HideNamespace scoped_name
        | HideNamespace begin_args '(' hidden_ns_args end_args ')'"""

    pm = p.parser.pm

    if pm.skipping:
        return

    if len(p) == 3:
        name = p[2]
    else:
        try:
            name = p[4]['name']
        except KeyError:
            pm.parser_error(p, 1,
                    "the name of the namespace has not been specified")
            name = ''

    namespace = pm.new_class(p, 1, IfaceFileType.NAMESPACE,
            normalised_scoped_name(name, pm.scope))
    namespace.is_hidden_namespace = True


def p_hidden_ns_args(p):
    """hidden_ns_args : hidden_ns_arg
        | hidden_ns_args ',' hidden_ns_arg"""

    if len(p) == 4:
        p[1].update(p[3])

    p[0] = p[1]


def p_hidden_ns_arg(p):
    "hidden_ns_arg : name '=' scoped_name"

    p[0] = {p[1]: p[3]}


# %If/%End ####################################################################

def p_if_start(p):
    "if_start : If '(' qualifiers ')'"

    pm = p.parser.pm

    if pm.skipping:
        skipping = True
    else:
        skipping = not p[3]

    pm.skip_stack.append(skipping)


def p_if_end(p):
    "if_end : End"

    pm = p.parser.pm

    if len(pm.skip_stack) >= 2:
        pm.skip_stack.pop()
    else:
        pm.parser_error(p, 1, "too many '%End' directives")


# %Import #####################################################################

def p_import(p):
    """import : Import need_eol import_simple EOL
        | Import begin_args '(' import_compound end_args ')'"""


def p_import_simple(p):
    "import_simple : file_path"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.push_file(p, 1, new_module=True)


def p_import_compound(p):
    "import_compound : import_args"

    pm = p.parser.pm

    if pm.skipping:
        return

    try:
        sip_file = p[1]['name']
    except KeyError:
        pm.parser_error(p, 1, "the name of the file has not been specified")

    pm.push_file(p, 1, sip_file, new_module=True)


def p_import_args(p):
    """import_args : import_arg
        | import_args ',' import_arg"""

    if len(p) == 4:
        p[1].update(p[3])

    p[0] = p[1]


def p_import_arg(p):
    "import_arg : name '=' file_path"

    p[0] = {p[1]: p[3]}


# %Include ####################################################################

def p_include(p):
    """include : Include need_eol include_simple EOL
        | Include begin_args '(' include_compound end_args ')'"""


def p_include_simple(p):
    "include_simple : file_path"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.push_file(p, 1)


def p_include_compund(p):
    "include_compound : include_args"

    pm = p.parser.pm

    if pm.skipping:
        return

    try:
        sip_file = p[1]['name']
    except KeyError:
        pm.parser_error(p, 1, "the name of the file has not been specified")

    pm.push_file(p, 1, sip_file, optional=p[1].get('optional', False))


def p_include_args(p):
    """include_args : include_arg
        | include_args ',' include_arg"""

    if len(p) == 4:
        p[1].update(p[3])

    p[0] = p[1]


def p_include_arg(p):
    """include_arg : name '=' file_path
        | optional '=' bool_value"""

    p[0] = {p[1]: p[3]}


# %InitialisationCode #########################################################

def p_init_code(p):
    "init_code : InitialisationCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.module_state.module.initialisation_code.append(p[2])


# %InstanceCode ###############################################################

def p_instance_code(p):
    "instance_code : InstanceCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.instance_code is not None:
        pm.parser_error(p, 1, "%InstanceCode can only be specified once")

    pm.scope.instance_code = p[2]


# %License ####################################################################

def p_license(p):
    """license : License STRING
        | License begin_args '(' license_args end_args ')'"""

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.module_state.module.license is not None:
        pm.parser_error(p, 1, "%License can only be specified once")

    if len(p) == 3:
        license = License(p[2])
    else:
        args = p[4]

        try:
            license_type = args['type']
        except KeyError:
            pm.parser_error(p, 1, "%License must specify 'type'")
            license_type = ''

        license = License(license_type)

        license.licensee = args.get('licensee')
        license.signature = args.get('signature')
        license.timestamp = args.get('timestamp')

    pm.module_state.module.license = license


def p_license_args(p):
    """license_args : license_arg
        | license_args ',' license_arg"""

    if len(p) == 4:
        p[1].update(p[3])

    p[0] = p[1]


def p_license_arg(p):
    """license_arg : licensee '=' STRING
        | signature '=' STRING
        | timestamp '=' STRING
        | type '=' STRING"""

    p[0] = {p[1]: p[3]}


# %MappedType #################################################################

# The mapped type annotations.
_MAPPED_TYPE_ANNOTATIONS = (
    'AllowNone',
    'NoAssignmentOperator',
    'NoCopyCtor',
    'NoDefaultCtor',
    'NoRelease',
    'PyName',
    'PyQtFlags',
    'TypeHint',
    'TypeHintIn',
    'TypeHintOut',
    'TypeHintValue',
)


def p_mapped_type(p):
    "mapped_type : mapped_type_head '{' mapped_type_body '}' ';'"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.validate_mapped_type(p, 1, pm.scope)

    pm.pop_scope()


def p_mapped_type_template(p):
    "mapped_type_template : mapped_type_template_head '{' mapped_type_body '}' ';'"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.validate_mapped_type(p, 1, pm.scope)

    pm.parsing_template = False

    pm.pop_scope()


def p_mapped_type_head(p):
    "mapped_type_head : MappedType base_type opt_annos"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.check_annotations(p, 3, "mapped type", _MAPPED_TYPE_ANNOTATIONS)

    pm.add_mapped_type(p, 1, p[2], p[3])


def p_mapped_type_template_head(p):
    "mapped_type_template_head : template_decl MappedType base_type opt_annos"

    # Note that we only use the template arguments to confirm that any simple
    # (ie. unscoped) names in the base type are to be substituted when the
    # template is instantiated.  Anything else is quietly ignored.

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.cpp_only(p, 1, "%MappedType templates")
    pm.check_annotations(p, 4, "mapped type", _MAPPED_TYPE_ANNOTATIONS)

    if p[3].type is not ArgumentType.TEMPLATE:
        pm.parser_error(p, 3, "%MappedType template must map a template type")

    # Check a template hasn't already been provided.
    for mtt in pm.spec.mapped_type_templates:
        if mtt.mapped_type.type.definition.cpp_name == p[3].definition.cpp_name and same_template_signature(mtt.mapped_type.type.definition.types, p[3].definition.types, deep=True):
            pm.parser_error(p, 1,
                    "a %MappedType template for this type has already been defined")

    # Use a dummy interface file.
    mapped_type = MappedType(IfaceFile(IfaceFileType.MAPPED_TYPE), p[3])
    pm.annotate_mapped_type(p, 4, mapped_type, p[4])

    mtt = MappedTypeTemplate(mapped_type, p[1])
    pm.spec.mapped_type_templates.insert(0, mtt)

    pm.push_scope(mapped_type)


def p_mapped_type_body(p):
    """mapped_type_body : mapped_type_line
        | mapped_type_body mapped_type_line"""


def p_mapped_type_line(p):
    """mapped_type_line : if_start
        | if_end
        | convert_from_type_code
        | convert_to_type_code
        | enum_decl
        | instance_code
        | mapped_type_function
        | release_code
        | type_code
        | type_header_code"""


def p_mapped_type_function(p):
    "mapped_type_function : static cpp_type NAME '(' opt_arg_list ')' opt_const opt_exceptions opt_annos opt_signature ';' opt_docstring premethod_code method_code"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.check_annotations(p, 9, "function", _FUNCTION_ANNOTATIONS)
    pm.apply_type_annotations(p, 9, p[2], p[9])

    overload = pm.add_function(p, 1, p[3], p[2], p[5], p[9], const=p[7],
            exceptions=p[8], cpp_signature=p[10], docstring=p[12],
            premethod_code=p[13], method_code=p[14])

    overload.is_static = True

    pm.validate_function(p, 1, overload)


# %ModuleHeaderCode ###########################################################

def p_module_header_code(p):
    "module_header_code : ModuleHeaderCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.module_state.module.module_header_code.append(p[2])


# %Module #####################################################################

def p_module(p):
    """module : Module dotted_name module_body
        | Module begin_args '(' module_args end_args ')' module_body"""

    pm = p.parser.pm

    if pm.skipping:
        return

    # See if this %Module is part of a %CompositeModule.
    if pm.spec.is_composite:
        # Historically we %Include modules although conceptually we actually
        # %Import them.  Ensure that the scopes etc. are correct in either
        # case.
        pm.ensure_import()

    module_state = pm.module_state
    module = module_state.module

    if module.fq_py_name is not None:
        pm.parser_error(p, 1, "%Module has already been specified")
        return

    if len(p) == 4:
        args = {'name': p[2]}
        body = p[3]
    else:
        args = p[4]
        body = p[7]

    try:
        module.fq_py_name = cached_name(pm.spec, str(args['name']))

        if pm.in_main_module:
            module.fq_py_name.used = True
    except KeyError:
        pm.parser_error(p, 1, "the name of the module has not been specified")

    # Configure the current module.
    module_state.all_raise_py_exception = args.get('all_raise_py_exception',
            False)
    module.default_virtual_error_handler = args.get(
            'default_VirtualErrorHandler')
    module_state.call_super_init = args.get('call_super_init')
    module_state.kw_args = args.get('keyword_arguments', KwArgs.NONE)

    c_bindings = args.get('language')
    if c_bindings is not None:
        if pm.c_bindings is None:
            pm.c_bindings = c_bindings
        elif pm.c_bindings != c_bindings:
            pm.parser_error(p, 1, "cannot mix 'C' and 'C++' modules")

    if 'py_ssize_t_clean' in args:
        module.py_ssize_t_clean = args['py_ssize_t_clean']

    if 'use_argument_names' in args:
        module.use_arg_names = args['use_argument_names']

    if 'use_limited_api' in args:
        module.use_limited_api = args['use_limited_api']

    for directive in body:
        if isinstance(directive, tuple):
            module_state.auto_py_name_rules.append(directive)
        elif isinstance(directive, Docstring):
            module.docstring = directive


def p_module_args(p):
    """module_args : module_arg
        | module_args ',' module_arg"""

    if len(p) == 4:
        p[1].update(p[3])

    p[0] = p[1]


def p_module_arg(p):
    """module_arg : all_raise_py_exception '=' bool_value
        | call_super_init '=' bool_value
        | default_VirtualErrorHandler '=' NAME
        | keyword_arguments '=' STRING
        | language '=' STRING
        | name '=' dotted_name
        | py_ssize_t_clean '=' bool_value
        | use_argument_names '=' bool_value
        | use_limited_api '=' bool_value"""

    pm = p.parser.pm

    if p[1] == 'keyword_arguments':
        value = pm.convert_kw_args(p, 3)
    elif p[1] == 'language':
        if p[3] == 'C':
            value = True
        elif p[3] == 'C++':
            value = False
        else:
            pm.parser_error(p, 3, "unsupported language '{0}'".format(p[3]))
            value = None
    else:
        value = p[3]

    p[0] = {p[1]: value}


def p_module_body(p):
    """module_body : '{' module_body_directives '}' ';'
        | empty"""

    p[0] = p[2] if len(p) == 5 else []


def p_module_body_directives(p):
    """module_body_directives : module_body_directive
        | module_body_directives module_body_directive"""

    if len(p) == 2:
        body = []
        value = p[1]
    else:
        body = p[1]
        value = p[2]

    if value is not None:
        body.append(value)

    p[0] = body


def p_module_body_directive(p):
    """module_body_directive : if_start
        | if_end
        | autopyname
        | docstring"""

    p[0] = p[1]


# %ModuleCode #################################################################

def p_module_code(p):
    "module_code : ModuleCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.module_state.module.module_code.append(p[2])


# %PickleCode #################################################################

def p_pickle_code(p):
    "pickle_code : PickleCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.pickle_code is not None:
        pm.parser_error(p, 1, "%PickleCode can only be specified once")

    pm.scope.pickle_code = p[2]


# %Platforms ##################################################################

def p_platforms(p):
    "platforms : Platforms '{' qualifier_list '}'"

    pm = p.parser.pm

    # Remember which platforms were selected.  Note that qualifiers are never
    # skipped.
    selected = []

    for qual_name in p[3]:
        if qual_name in pm.tags:
            selected.append(qual_name)

        pm.add_qualifier(p, 1, qual_name, QualifierType.PLATFORM)

    # Check that no more than one platform was selected.
    if len(selected) > 1:
        pm.parser_error(p, 1,
                "only one of {0} can be selected from %Platforms".format(
                        ','.join(["'" + s + "'" for s in selected])))


# %Plugin #####################################################################

def p_plugin(p):
    "plugin : Plugin NAME"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.spec.plugins.append(p[2])


# %PostInitialisationCode #####################################################

def p_postinit_code(p):
    "postinit_code : PostInitialisationCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.module_state.module.postinitialisation_code.append(p[2])


# %PreInitialisationCode ######################################################

def p_preinit_code(p):
    "preinit_code : PreInitialisationCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.module_state.module.preinitialisation_code.append(p[2])


# %Property ###################################################################

def p_property(p):
    "property : Property begin_args '(' property_args end_args ')' opt_property_body"

    pm = p.parser.pm

    if pm.skipping:
        return

    name = p[4].get('name')
    if name is None:
        p.parser.pm.parser_error(p, 1,
                "a name must be specified for %Property")
        return

    pm.check_attributes(p, 4, name)

    name = cached_name(pm.spec, name)
    if pm.in_main_module:
        name.used = True

    getter = p[4].get('get')
    if getter is None:
        p.parser.pm.parser_error(p, 1,
                "a getter must be specified for %Property")
        return

    prop = Property(name=name, getter=getter, setter=p[4].get('set'))

    for directive in p[7]:
        if isinstance(directive, Docstring):
            prop.docstring = directive

    pm.scope.properties.insert(0, prop)


def p_property_args(p):
    """property_args : property_arg
        | property_args ',' property_arg"""

    if len(p) == 4:
        p[1].update(p[3])

    p[0] = p[1]


def p_property_arg(p):
    """property_arg : get '=' NAME
        | name '=' NAME
        | set '=' NAME"""

    p[0] = {p[1]: p[3]}


def p_opt_property_body(p):
    """opt_property_body : empty
        | '{' property_body '}' ';'"""

    p[0] = [] if len(p) == 2 else p[2]


def p_property_body(p):
    """property_body : property_line
        | property_body property_line"""

    if len(p) == 2:
        value = []
        directive = p[1]
    else:
        value = p[1]
        directive = p[2]

    if directive is not None:
        value.append(directive)

    p[0] = value


def p_property_line(p):
    """property_line : if_start
        | if_end
        | docstring"""

    p[0] = p[1]


# %ReleaseCode ################################################################

def p_release_code(p):
    "release_code : ReleaseCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.release_code is not None:
        pm.parser_error(p, 1, "%ReleaseCode can only be specified once")

    pm.scope.release_code = p[2]


# %Timeline ###################################################################

def p_timeline(p):
    "timeline : Timeline '{' qualifier_list '}'"

    pm = p.parser.pm

    # Allocate the timeline number.
    module_state = pm.module_state
    timeline = module_state.nr_timelines
    module_state.nr_timelines += 1

    # Remember which tags were selected.  Note that qualifiers are never
    # skipped.
    selected = []

    for order, qual_name in enumerate(p[3]):
        if qual_name in pm.tags:
            selected.append(qual_name)

        pm.add_qualifier(p, 1, qual_name, QualifierType.TIME, order=order,
                timeline=timeline)

    # Check that no more than one qualifier was selected.
    if len(selected) > 1:
        pm.parser_error(p, 1,
                "only one of {0} can be selected from this %Timeline".format(
                        ','.join(["'" + s + "'" for s in selected])))


# %TypeCode ###################################################################

def p_type_code(p):
    "type_code : TypeCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.scope.type_code.append(p[2])


# %TypeHeaderCode #############################################################

def p_type_header_code(p):
    "type_header_code : TypeHeaderCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.scope.iface_file.type_header_code.append(p[2])


# %TypeHintCode ###############################################################

def p_type_hint_code(p):
    "type_hint_code : TypeHintCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    # This can be specified multiple times for a module but only once for a
    # class.
    scope = pm.scope

    if scope is None:
        pm.module_state.module.type_hint_code.append(p[2])
    else:
        if scope.type_hint_code is not None:
            pm.parser_error(p, 1, "%TypeHintCode can only be specified once")

        scope.type_hint_code = p[2]


# %UnitCode ###################################################################

def p_unit_code(p):
    "unit_code : UnitCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.module_state.module.unit_code.append(p[2])


# %UnitPostIncludeCode ########################################################

def p_unit_postinclude_code(p):
    "unit_postinclude_code : UnitPostIncludeCode CODE_BLOCK"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.module_state.module.unit_postinclude_code.append(p[2])


# %VirtualErrorHandler ########################################################

def p_virtual_error_handler(p):
    """virtual_error_handler : VirtualErrorHandler NAME CODE_BLOCK
        | VirtualErrorHandler begin_args '(' veh_args end_args ')' CODE_BLOCK"""

    pm = p.parser.pm

    if pm.skipping:
        return

    if len(p) == 4:
        name = p[2]
        code = p[3]
    else:
        name = p[4].get('name')
        if name is None:
            pm.parser_error(p, 1,
                    "a name must be specified for %VirtualErrorHandler")
            return

        code = p[7]

    for handler in pm.spec.virtual_error_handlers:
        if handler.name == name:
            pm.parser_error(p, 1,
                    "a virtual error handler called '{0}' has already been defined".format(name))
            break
    else:
        pm.spec.virtual_error_handlers.append(
                VirtualErrorHandler(code, pm.module_state.module, name))


def p_veh_args(p):
    """veh_args : veh_arg
        | veh_args ',' veh_arg"""

    if len(p) == 4:
        p[1].update(p[3])

    p[0] = p[1]


def p_veh_arg(p):
    "veh_arg : name '=' NAME"

    p[0] = {p[1]: p[3]}


# A C/C++ type. ###############################################################

def p_cpp_type(p):
    """cpp_type : const base_type derefs opt_ref
        | base_type derefs opt_ref"""

    if len(p) == 5:
        value = p[2]
        value.is_const = True
        value.derefs.extend(p[3])
        value.is_reference = p[4]
    else:
        value = p[1]
        value.derefs.extend(p[2])
        value.is_reference = p[3]

        # PyObject * is a synonym for SIP_PYOBJECT.
        if value.type is ArgumentType.DEFINED and len(value.definition) == 1 and value.definition.base_name == 'PyObject' and len(value.derefs) == 1 and not value.is_reference:
            value.type = ArgumentType.PYOBJECT
            value.definition = None
            value.derefs = []

    p[0] = value


def p_base_type(p):
    """base_type : pod_type
        | scoped_name
        | scoped_name '<' cpp_types '>'
        | struct scoped_name
        | union scoped_name"""

    pm = p.parser.pm

    if isinstance(p[1], ArgumentType):
        p[0] = Argument(p[1])
    elif len(p) == 2:
        # Resolve it if it is the name of a typedef.  This is done early as a
        # workaround for allowing /PyInt/ to be applied to typedef'ed types.
        ad = Argument(ArgumentType.DEFINED, definition=p[1])

        while ad.type is ArgumentType.DEFINED:
            ad.type = ArgumentType.NONE
            search_typedefs(pm.spec, ad.definition, ad)

            # Don't resolve to a template type as it may be superceded later on
            # by a more specific mapped type.
            if ad.type in (ArgumentType.NONE, ArgumentType.TEMPLATE):
                ad = Argument(ArgumentType.DEFINED, definition=p[1])
                break

        p[0] = ad
    elif len(p) == 3:
        # In a C module all structs and unions must be defined.
        if bool(p.parser.pm.c_bindings):
            type = ArgumentType.DEFINED
        elif p[1] == 'struct':
            type = ArgumentType.STRUCT
        else:
            type = ArgumentType.UNION

        p[0] = Argument(type, definition=p[2])
    else:
        p[0] = Argument(ArgumentType.TEMPLATE,
                definition=Template(p[1], Signature(args=p[3])))

    p[0].source_location = pm.get_source_location(p, 1)


# Map unsigned two-word POD types.
_UNSIGNED_MAP = {
    'char':     ArgumentType.USTRING,
    'int':      ArgumentType.UINT,
    'long':     ArgumentType.ULONG,
    'short':    ArgumentType.USHORT,
}

# Map one-word POD types.
_ONE_WORD_MAP = {
    '...':              ArgumentType.ELLIPSIS,
    'void':             ArgumentType.VOID,
    'bool':             ArgumentType.BOOL,
    'char':             ArgumentType.STRING,
    'double':           ArgumentType.DOUBLE,
    'float':            ArgumentType.FLOAT,
    'short':            ArgumentType.SHORT,
    'int':              ArgumentType.INT,
    'unsigned':         ArgumentType.UINT,
    'long':             ArgumentType.LONG,
    'size_t':           ArgumentType.SIZE,
    'wchar_t':          ArgumentType.WSTRING,
    'Py_hash_t':        ArgumentType.HASH,
    'Py_ssize_t':       ArgumentType.SSIZE,
    'SIP_PYBUFFER':     ArgumentType.PYBUFFER,
    'SIP_PYCALLABLE':   ArgumentType.PYCALLABLE,
    'SIP_PYDICT':       ArgumentType.PYDICT,
    'SIP_PYENUM':       ArgumentType.PYENUM,
    'SIP_PYLIST':       ArgumentType.PYLIST,
    'SIP_PYOBJECT':     ArgumentType.PYOBJECT,
    'SIP_PYSLICE':      ArgumentType.PYSLICE,
    'SIP_PYTUPLE':      ArgumentType.PYTUPLE,
    'SIP_PYTYPE':       ArgumentType.PYTYPE,
    'SIP_SSIZE_T':      ArgumentType.SSIZE,
}


def p_pod_type(p):
    """pod_type : unsigned long long
        | signed char
        | long long
        | unsigned char
        | unsigned short
        | unsigned int
        | unsigned long
        | unsigned
        | short
        | int
        | long
        | float
        | double
        | bool
        | char
        | wchar_t
        | void
        | SIP_PYOBJECT
        | SIP_PYTUPLE
        | SIP_PYLIST
        | SIP_PYDICT
        | SIP_PYCALLABLE
        | SIP_PYSLICE
        | SIP_PYTYPE
        | SIP_PYBUFFER
        | SIP_PYENUM
        | SIP_SSIZE_T
        | Py_hash_t
        | Py_ssize_t
        | size_t
        | ELLIPSIS"""

    if len(p) == 4:
        p[0] = ArgumentType.ULONGLONG
    elif len(p) == 3:
        if p[1] == 'signed':
            p[0] = ArgumentType.SSTRING
        elif p[1] == 'long':
            p[0] = ArgumentType.LONGLONG
        else:
            p[0] = _UNSIGNED_MAP[p[2]]
    else:
        p[0] = _ONE_WORD_MAP[p[1]]


def p_cpp_types(p):
    """cpp_types : cpp_type
        | cpp_types ',' cpp_type"""

    if len(p) == 2:
        p[0] = [p[1]]
    else:
        if p[1][-1].type == ArgumentType.ELLIPSIS:
            p.parser.pm.parser_error(p, 3,
                    "an ellipsis must be at the end of the argument list")

        p[0] = p[1]
        p[0].append(p[3])


def p_derefs(p):
    """derefs : empty
        | derefs '*'
        | derefs '*' const"""

    if len(p) == 2:
        p[0] = []
    else:
        p[0] = p[1]
        p[0].append(len(p) == 4)


def p_opt_ref(p):
    """opt_ref : '&'
        | empty"""

    if p[1] is None:
        p[0] = False
    else:
        p.parser.pm.cpp_only(p, 1, "references")

        p[0] = True


# C++ classes and structs. ####################################################

# The class annotations.
_CLASS_ANNOTATIONS = (
    'Abstract',
    'AllowNone',
    'DelayDtor',
    'Deprecated',
    'ExportDerived',
    'External',
    'FileExtension',
    'Metatype',
    'Mixin',
    'NoDefaultCtors',
    'NoTypeHint',
    'PyName',
    'PyQtFlags',
    'PyQtFlagsEnums',
    'PyQtInterface',
    'PyQtNoQMetaObject',
    'Supertype',
    'TypeHint',
    'TypeHintIn',
    'TypeHintOut',
    'TypeHintValue',
    'VirtualErrorHandler',
)


def p_class_template(p):
    "class_template : template_decl class_decl"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.cpp_only(p, 1, "class templates")

    pm.class_templates.append((p[1], p[2]))

    pm.parsing_template = False


def p_class_docstring(p):
    "class_docstring : docstring"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope.docstring is None:
        pm.scope.docstring = p[1]
    else:
        pm.parser_error(p, 1,
                "%Docstring has already been defined for this class")


def p_class_decl(p):
    "class_decl : class class_head opt_class_definition ';'"

    pm = p.parser.pm

    if pm.skipping:
        return

    # Return the class for any class template.
    p[0] = pm.complete_class(p, 2, p[2], p[3])


def p_class_head(p):
    "class_head : scoped_name superclasses opt_annos"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.check_annotations(p, 3, "class", _CLASS_ANNOTATIONS)

    if p[2] is not None:
        pm.cpp_only(p, 2, "super-classes")

    pm.define_class(p, 1, ClassKey.CLASS, p[1], p[3], superclasses=p[2])

    # Return the annotations.
    p[0] = p[3]


def p_struct_decl(p):
    "struct_decl : struct struct_head opt_class_definition ';'"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.complete_class(p, 2, p[2], p[3])


def p_struct_head(p):
    "struct_head : scoped_name superclasses opt_annos"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.check_annotations(p, 3, "class", _CLASS_ANNOTATIONS)

    pm.define_class(p, 1, ClassKey.STRUCT, p[1], p[3], superclasses=p[2])

    # Return the annotations.
    p[0] = p[3]


def p_superclasses(p):
    """superclasses : ':' superclass_list
        | empty"""

    p[0] = [] if len(p) == 2 else p[2]


def p_superclass_list(p):
    """superclass_list : superclass
        | superclass_list ',' superclass"""

    if len(p) == 2:
        l = []
        value = p[1]
    else:
        l = p[1]
        value = p[3]

    if value is not None:
        l.append(value)

    p[0] = l


def p_superclass(p):
    "superclass : class_access scoped_name"

    pm = p.parser.pm

    if pm.skipping:
        return

    if p[1] != 'public':
        p[0] = None
        return

    # This is a hack to allow typedef'ed classes to be used before we have
    # resolved the typedef definitions.  Unlike elsewhere, we require that the
    # typedef is defined before being used.
    ad = Argument(ArgumentType.DEFINED, definition=p[2])

    while ad.type is ArgumentType.DEFINED:
        ad.type = ArgumentType.NONE
        search_typedefs(pm.spec, ad.definition, ad)

    if ad.type is not ArgumentType.NONE or len(ad.derefs) != 0 or ad.is_const or ad.is_reference:
        pm.parser_error(p, 2, "super-class list contains an invalid type")

    # Find the actual class.
    p[0] = pm.find_class(p, 2, IfaceFileType.CLASS, ad.definition,
            tmpl_arg=pm.parsing_template)


def p_class_access(p):
    """class_access : empty
        | public
        | protected
        | private"""

    if p[1] is None:
        p[1] = 'public'

    p[0] = p[1]


def p_opt_class_definition(p):
    """opt_class_definition : '{' opt_class_body '}'
        | empty"""

    pm = p.parser.pm

    if pm.skipping:
        return

    p[0] = (len(p) == 4)


def p_opt_class_body(p):
    """opt_class_body : class_body
        | empty"""


def p_class_body(p):
    """class_body : class_line
        | class_body class_line"""


def p_class_line(p):
    """class_line : if_start
        | if_end
        | class_decl
        | class_docstring
        | class_template
        | ctor
        | dtor
        | enum_decl
        | exception
        | typedef_decl
        | method_variable
        | namespace_decl
        | struct_decl
        | union_decl
        | public_specifier
        | protected_specifier
        | private_specifier
        | signals_specifier
        | convert_from_type_code
        | convert_to_subclass_code
        | convert_to_type_code
        | finalisation_code
        | gc_clear_code
        | gc_traverse_code
        | get_buffer_code
        | instance_code
        | pickle_code
        | property
        | release_buffer_code
        | type_code
        | type_header_code
        | type_hint_code
        | BIGetReadBufferCode CODE_BLOCK
        | BIGetWriteBufferCode CODE_BLOCK
        | BIGetSegCountCode CODE_BLOCK
        | BIGetCharBufferCode CODE_BLOCK"""


# The ctor annotations.
_CTOR_ANNOTATIONS = (
    'Default',
    'Deprecated',
    'HoldGIL',
    'KeywordArgs',
    'NoDerived',
    'NoRaisesPyException',
    'NoTypeHint',
    'PostHook',
    'PreHook',
    'RaisesPyException',
    'ReleaseGIL',
    'Transfer',
)


def p_ctor(p):
    """ctor : explicit ctor_decl
        | ctor_decl"""

    pm = p.parser.pm

    if pm.skipping:
        return

    if len(p) == 3:
        pm.cpp_only(p, 1, "explicit constructors")


def p_ctor_decl(p):
    "ctor_decl : NAME '(' opt_arg_list ')' opt_exceptions opt_annos opt_ctor_signature ';' opt_docstring premethod_code method_code"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.check_annotations(p, 6, "constructor", _CTOR_ANNOTATIONS)

    pm.add_ctor(p, 1, p[3], p[6], exceptions=p[5], cpp_signature=p[7],
            docstring=p[9], premethod_code=p[10], method_code=p[11])


def p_opt_ctor_signature(p):
    """opt_ctor_signature : '[' '(' opt_arg_list ')' ']'
        | empty"""

    if len(p) == 6:
        p[0] = Signature(args=p[3])
    else:
        p[0] = None


# The dtor annotations.
_DTOR_ANNOTATIONS = (
    'HoldGIL',
    'ReleaseGIL',
)


def p_dtor(p):
    "dtor : opt_virtual '~' NAME '(' ')' opt_exceptions opt_abstract opt_annos ';' premethod_code method_code virtual_catcher_code"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.check_annotations(p, 8, "destructor", _DTOR_ANNOTATIONS)

    pm.add_dtor(p, 1, p[3], p[8], exceptions=p[6], abstract=p[7],
            premethod_code=p[10], method_code=p[11],
            virtual_catcher_code=p[12])

    pm.parsing_virtual = False


def p_method_variable(p):
    """method_variable : Q_SIGNAL simple_method_variable
        | Q_SLOT simple_method_variable
        | simple_method_variable"""

    pm = p.parser.pm

    if len(p) == 3:
        item = p[2]

        if isinstance(item, Overload):
            item.pyqt_method_specifier = PyQtMethodSpecifier.SIGNAL if p[1] == 'Q_SIGNAL' else PyQtMethodSpecifier.SLOT
        else:
            pm.parser_error(p, 1,
                    "a PyQt method specifier can only be applied to member functions")
    else:
        item = p[1]

    if isinstance(item, Overload):
        pm.validate_function(p, 1, item)


def p_simple_method_variable(p):
    """simple_method_variable : virtual function
        | static plain_method_variable
        | plain_method_variable"""

    if len(p) == 3:
        item = p[2]

        if item is not None:
            if p[1] == 'static':
                item.is_static = True
            elif not item.is_final:
                item.is_virtual = True
                p.parser.pm.scope.needs_shadow = True
    else:
        item = p[1]

    p[0] = item


def p_plain_method_variable(p):
    """plain_method_variable : function
        | variable"""

    p[0] = p[1]


def p_public_specifier(p):
    "public_specifier : public opt_slots ':'"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.scope_access_specifier = AccessSpecifier.PUBLIC
    pm.scope_pyqt_method_specifier = p[2]


def p_protected_specifier(p):
    "protected_specifier : protected opt_slots ':'"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.scope_access_specifier = AccessSpecifier.PROTECTED
    pm.scope_pyqt_method_specifier = p[2]


def p_private_specifier(p):
    "private_specifier : private opt_slots ':'"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.scope_access_specifier = AccessSpecifier.PRIVATE
    pm.scope_pyqt_method_specifier = p[2]


def p_signals_specifier(p):
    """signals_specifier : signals ':'
        | Q_SIGNALS ':'"""

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.scope_access_specifier = AccessSpecifier.PUBLIC
    pm.scope_pyqt_method_specifier = PyQtMethodSpecifier.SIGNAL


def p_opt_slots(p):
    """opt_slots : slots
        | Q_SLOTS
        | empty"""

    p[0] = None if p[1] is None else PyQtMethodSpecifier.SLOT


# C/C++ enums. ################################################################

# The enum annotations.
_ENUM_ANNOTATIONS = (
    'BaseType',
    'NoScope',
    'NoTypeHint',
    'PyName',
)

# The enum member annotations.
_ENUM_MEMBER_ANNOTATIONS = (
    'NoTypeHint',
    'PyName',
)


def p_enum_decl(p):
    "enum_decl : enum opt_enum_key opt_name opt_annos '{' opt_enum_body '}' ';'"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.check_annotations(p, 4, "enum", _ENUM_ANNOTATIONS)

    pm.add_enum(p, 1, p[3], p[2], p[4], p[6])


def p_opt_enum_key(p):
    """opt_enum_key : class
        | struct
        | union
        | empty"""

    # Return True if the enum is scoped.
    p[0] = p[1] is not None


def p_opt_enum_body(p):
    """opt_enum_body : enum_body
        | empty"""

    p[0] = [] if p[1] is None else p[1]


def p_enum_body(p):
    """enum_body : enum_line
        | enum_body enum_line"""

    if len(p) == 2:
        p[0] = []
        line = p[1]
    else:
        p[0] = p[1]
        line = p[2]

    if line is not None:
        p[0].append(line)


def p_enum_line(p):
    """enum_line : if_start
        | if_end
        | NAME opt_enum_assign opt_annos opt_comma"""

    pm = p.parser.pm

    if pm.skipping:
        p[0] = None
        return

    if len(p) == 5:
        pm.check_annotations(p, 3, "enum member", _ENUM_MEMBER_ANNOTATIONS)

        p[0] = (p[1], cached_name(pm.spec, pm.get_py_name(p[1], p[3])),
                p[3].get('NoTypeHint', False))
    else:
        p[0] = None


def p_opt_enum_assign(p):
    """opt_enum_assign : '=' value
        | empty"""


def p_opt_comma(p):
    """opt_comma : empty
        | ','"""


# C++ exceptions. #############################################################

# The exception annotations.
_EXCEPTION_ANNOTATIONS = (
    'Default',
    'PyName',
)


def p_exception(p):
    "exception : Exception scoped_name opt_base_exception opt_annos '{' exception_body '}' ';'"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.cpp_only(p, 1, "%Exception")
    pm.check_annotations(p, 4, "exception", _EXCEPTION_ANNOTATIONS)

    cpp_name = p[2]
    py_name = pm.get_py_name(cpp_name.base_name, p[4])
    pm.check_attributes(p, 2, py_name)

    builtin_base, defined_base = p[3]
    raise_code = p[6].get('%RaiseCode')
    type_header_code = p[6].get('%TypeHeaderCode')

    if raise_code is None:
        pm.parser_error(p, 1,
                "%Exception must have a %RaiseCode sub-directive")

    xd = pm.find_exception(p, 1, cpp_name, raise_code=raise_code)

    module = pm.module_state.module

    if xd.iface_file.module is not None:
        pm.parser_error(p, 2,
                "an %Exception with this name has already been defined")

    xd.iface_file.module = module

    if type_header_code is not None:
        xd.iface_file.type_header_code.append(type_header_code)

    xd.builtin_base_exception = builtin_base
    xd.defined_base_exception = defined_base
    xd.py_name = py_name

    if p[4].get('Default'):
        if module.default_exception is None:
            module.default_exception = xd
        else:
            pm.parser_error(p, 1,
                    "another %Exception has already been annotated with '/Default/'")


def p_opt_base_exception(p):
    """opt_base_exception : '(' scoped_name ')'
        | empty"""

    builtin_base = None
    defined_base = None

    if len(p) == 4:
        base = p[2]
        base.make_absolute()

        # See if it is a project-defined exception.
        for xd in p.parser.pm.spec.exceptions:
            if xd.iface_file.fq_cpp_name == base:
                defined_base = xd
                break
        else:
            # See if it is a builtin Python exception.
            if len(base) == 1 and base.base_name.startswith('SIP_'):
                py_name = base.base_name[4:]

                from .python_exceptions import PYTHON_EXCEPTIONS

                if py_name in PYTHON_EXCEPTIONS:
                    builtin_base = py_name

    p[0] = (builtin_base, defined_base)


def p_exception_body(p):
    """exception_body : exception_line
        | exception_body exception_line"""

    if len(p) == 3:
        p[1].update(p[2])

    p[0] = p[1]


def p_exception_line(p):
    """exception_line : if_start
        | if_end
        | RaiseCode CODE_BLOCK
        | TypeHeaderCode CODE_BLOCK"""

    p[0] = {} if len(p) == 2 or p.parser.pm.skipping else {p[1]: p[2]}


# C/C++ functions. ############################################################

# The function annotations.
_FUNCTION_ANNOTATIONS = (
    '__len__',
    '__imatmul__',
    '__matmul__',
    'AbortOnException',
    'AllowNone',
    'AutoGen',
    'Deprecated',
    'DisallowNone',
    'Encoding',
    'Factory',
    'HoldGIL',
    'KeywordArgs',
    'KeepReference',
    'NewThread',
    'NoArgParser',
    'NoCopy',
    'NoRaisesPyException',
    'NoTypeHint',
    'NoVirtualErrorHandler',
    'Numeric',
    'PostHook',
    'PreHook',
    'PyInt',
    'PyName',
    'PyQtSignalHack',
    'RaisesPyException',
    'ReleaseGIL',
    'Sequence',
    'VirtualErrorHandler',
    'Transfer',
    'TransferBack',
    'TransferThis',
    'TypeHint',
)


def p_function(p):
    """function : function_decl
        | assignment_operator_decl
        | operator_decl
        | operator_cast_decl"""

    p[0] = p[1]


def p_function_decl(p):
    "function_decl : cpp_type NAME '(' opt_arg_list ')' opt_const opt_final opt_exceptions opt_abstract opt_annos opt_signature ';' opt_docstring premethod_code method_code virtual_catcher_code virtual_call_code"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.check_annotations(p, 10, "function", _FUNCTION_ANNOTATIONS)
    pm.apply_type_annotations(p, 10, p[1], p[10])

    p[0] = pm.add_function(p, 1, p[2], p[1], p[4], p[10], const=p[6],
            final=p[7], exceptions=p[8], abstract=p[9], cpp_signature=p[11],
            docstring=p[13], premethod_code=p[14], method_code=p[15],
            virtual_catcher_code=p[16], virtual_call_code=p[17])


def p_assignment_operator_decl(p):
    "assignment_operator_decl : cpp_type operator '=' '(' cpp_type ')' ';'"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.scope_access_specifier is AccessSpecifier.PRIVATE:
        pm.scope.cannot_assign = True
    else:
        pm.parser_error(p, 2,
                "an assignment operator may only be specified in the private section of a class")


def p_operator_decl(p):
    "operator_decl : cpp_type operator operator_name '(' opt_arg_list ')' opt_const opt_final opt_exceptions opt_abstract opt_annos opt_signature ';' premethod_code method_code virtual_catcher_code virtual_call_code"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.check_annotations(p, 11, "function", _FUNCTION_ANNOTATIONS)
    pm.apply_type_annotations(p, 11, p[1], p[11])

    scope = pm.scope

    # Handle the unary '+' and '-' operators.
    if (scope is not None and len(p[5]) == 0) or (scope is None and len(p[5]) == 1):
        if p[3] == '__add__':
            p[3] = '__pos__'
        elif p[3] == '__sub__':
            p[3] = '__neg__'

    p[0] = pm.add_function(p, 1, p[3], p[1], p[5], p[11], const=p[7],
            final=p[8], exceptions=p[9], abstract=p[10], cpp_signature=p[12],
            premethod_code=p[14], method_code=p[15],
            virtual_catcher_code=p[16], virtual_call_code=p[17])


# Types that can be cast to a Python int.
_INT_TYPES = (
    ArgumentType.BOOL,
    ArgumentType.CBOOL,
    ArgumentType.BYTE,
    ArgumentType.SBYTE,
    ArgumentType.UBYTE,
    ArgumentType.SHORT,
    ArgumentType.USHORT,
    ArgumentType.INT,
    ArgumentType.CINT,
    ArgumentType.UINT,
    ArgumentType.LONG,
    ArgumentType.ULONG,
    ArgumentType.LONGLONG,
    ArgumentType.ULONGLONG,
)

# Types that can be cast to a Python float.
_FLOAT_TYPES = (
    ArgumentType.FLOAT,
    ArgumentType.CFLOAT,
    ArgumentType.DOUBLE,
    ArgumentType.CDOUBLE,
)

def p_operator_cast_decl(p):
    "operator_cast_decl : operator cpp_type '(' opt_arg_list ')' opt_const opt_final opt_exceptions opt_abstract opt_annos opt_signature ';' premethod_code method_code virtual_catcher_code virtual_call_code"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.check_annotations(p, 10, "function", _FUNCTION_ANNOTATIONS)
    pm.apply_type_annotations(p, 10, p[2], p[10])

    if pm.scope is None:
        pm.parser_error(p, 1,
                "operator casts can only be specified as part of a class")

    if len(p[4]) != 0:
        pm.parser_error(p, 1, "operator casts cannot have any arguments")

    # Get the slot name if it is implemented as a Python cast.
    slot_name = None
    cpp_type = p[2]

    if cpp_type.type in _INT_TYPES:
        slot_name = '__int__'
    elif cpp_type.type in _FLOAT_TYPES:
        slot_name = '__float__'
    elif cpp_type.type is not ArgumentType.DEFINED:
        pm.parser_error(p, 2, "unsupported operator cast")

    if slot_name is None:
        # Check the cast hasn't already been specified.
        for cast in pm.scope.casts:
            if cast.definition == cpp_type.definition:
                pm.parser_error(p, 2,
                        "operator cast has already been specified for this class")
                break
        else:
            pm.scope.casts.insert(0, cpp_type)
    else:
        p[0] = pm.add_function(p, 1, slot_name, p[2], p[4], p[10], const=p[6],
                final=p[7], exceptions=p[8], abstract=p[9],
                cpp_signature=p[11], premethod_code=p[13], method_code=p[14],
                virtual_catcher_code=p[15], virtual_call_code=p[16])


def p_opt_arg_list(p):
    """opt_arg_list : arg_list
        | empty"""

    if p[1] is None:
        arg_list = []
    else:
        arg_list = p[1]

        # Handle C void prototypes.
        if bool(p.parser.pm.c_bindings) and len(arg_list) == 1:
            arg = arg_list[0]

            if arg.type is ArgumentType.VOID and len(arg.derefs) == 0:
                arg_list = []

    p[0] = arg_list


def p_arg_list(p):
    """arg_list : arg_value
        | arg_list ',' arg_value"""

    if len(p) == 2:
        arg_list = [p[1]]
    else:
        arg_list = p[1]
        arg_list.append(p[3])

    p[0] = arg_list


def p_arg_value(p):
    "arg_value : arg_type opt_assign"

    p[0] = p[1]
    p[0].default_value = p[2]


# The argument annotations.
_ARGUMENT_ANNOTATIONS = (
    'AllowNone',
    'Array',
    'ArraySize',
    'Constrained',
    'DisallowNone',
    'Encoding',
    'GetWrapper',
    'In',
    'KeepReference',
    'NoCopy',
    'Out',
    'PyInt',
    'ResultSize',
    'ScopesStripped',
    'Transfer',
    'TransferBack',
    'TransferThis',
    'TypeHint',
    'TypeHintIn',
    'TypeHintOut',
    'TypeHintValue',
)


def p_arg_type(p):
    "arg_type : cpp_type opt_name opt_annos"

    pm = p.parser.pm

    pm.check_annotations(p, 3, "argument", _ARGUMENT_ANNOTATIONS)

    arg = p[1]
    annotations = p[3]

    if p[2] is not None:
        arg.name = cached_name(pm.spec, p[2])

    pm.apply_common_argument_annotations(p, 3, arg, annotations)

    has_array = 'Array' in annotations
    has_array_size = 'ArraySize' in annotations

    if has_array:
        if has_array_size:
            pm.parser_error(p, 3,
                    "/Array/ and /ArraySize/ cannot both be specified")

        arg.array = ArrayArgument.ARRAY
    elif has_array_size:
        arg.array = ArrayArgument.ARRAY_SIZE

    pm.apply_type_annotations(p, 3, arg, annotations)

    if 'Constrained' in annotations:
        arg.is_constrained = True

        if arg.type is ArgumentType.BOOL:
            arg.type = ArgumentType.CBOOL
        elif arg.type is ArgumentType.INT:
            arg.type = ArgumentType.CINT
        elif arg.type is ArgumentType.FLOAT:
            arg.type = ArgumentType.CFLOAT
        elif arg.type is ArgumentType.DOUBLE:
            arg.type = ArgumentType.CDOUBLE

    arg.get_wrapper = annotations.get('GetWrapper', False)
    arg.is_in = annotations.get('In', False)
    arg.is_out = annotations.get('Out', False)
    arg.result_size = annotations.get('ResultSize', False)

    try:
        arg.scopes_stripped = annotations['ScopesStripped']

        if arg.scopes_stripped <= 0:
            pm.parser_error(p, 3, "/ScopesStripped/ must be greater than 0")
    except KeyError:
        arg.scopes_stripped = 0

    arg.transfer = pm.get_transfer(p, 3, annotations)

    p[0] = arg


def p_opt_assign(p):
    """opt_assign : '=' expr
        | empty"""

    p[0] = p[2] if len(p) == 3 else None


def p_expr(p):
    """expr : value
        | expr binop value"""

    if len(p) == 2:
        expr = [p[1]]
    else:
        expr = p[1]
        expr[-1].binary_operator = p[2]
        expr.append(p[3])

    p[0] = expr


def p_value(p):
    "value : opt_cast opt_unop simple_value"

    value = p[3]
    value.cast = p[1]
    value.unary_operator = p[2]

    p[0] = value


def p_simple_value(p):
    """simple_value : empty_value
        | function_call_value
        | null_value
        | number_value
        | quoted_char_value
        | real_value
        | scoped_name_value
        | string_value"""

    p[0] = p[1]


def p_empty_value(p):
    "empty_value : '{' '}'"

    p[0] = Value(ValueType.EMPTY, None)


def p_function_call_value(p):
    "function_call_value : base_type '(' opt_expr_list ')'"

    p[0] = Value(ValueType.FCALL, FunctionCall(p[1], p[3]))


def p_null_value(p):
    "null_value : NULL"

    p[0] = Value(ValueType.NUMERIC, 0)


def p_number_value(p):
    """number_value : NUMBER
        | bool_value"""

    p[0] = Value(ValueType.NUMERIC, p[1])


def p_quoted_char_value(p):
    "quoted_char_value : QUOTED_CHAR"

    p[0] = Value(ValueType.QCHAR, p[1])


def p_real_value(p):
    "real_value : REAL"

    p[0] = Value(ValueType.REAL, p[1])


def p_scoped_name_value(p):
    "scoped_name_value : scoped_name"

    p[0] = Value(ValueType.SCOPED, p[1])


def p_string_value(p):
    "string_value : STRING"

    p[0] = Value(ValueType.STRING, p[1])


def p_opt_expr_list(p):
    """opt_expr_list : expr_list
        | empty"""

    p[0] = [] if p[1] is None else p[1]


def p_expr_list(p):
    """expr_list : expr
        | expr_list ',' expr"""

    if len(p) == 2:
        value = [p[1]]
    else:
        value = p[1]
        value.append(p[3])

    p[0] = value


def p_opt_cast(p):
    """opt_cast : '(' scoped_name ')'
        | empty"""

    p[0] = p[2] if len(p) == 4 else None


def p_binop(p):
    """binop : '-'
        | '+'
        | '*'
        | '/'
        | '&'
        | '|'"""

    p[0] = p[1]


def p_opt_unop(p):
    """opt_unop : empty
        | '!'
        | '~'
        | '-'
        | '+'
        | '*'
        | '&'"""

    p[0] = p[1]


def p_opt_exceptions(p):
    """opt_exceptions : empty
        | noexcept
        | throw '(' opt_exception_list ')'"""

    if p[1] == 'throw':
        p[0] = p[3]
    elif p[1] == 'noexcept':
        p[0] = ThrowArguments()
    else:
        p[0] = None


def p_opt_exception_list(p):
    """opt_exception_list : exception_list
        | empty"""

    p[0] = ThrowArguments(arguments=[]) if p[1] is None else p[1]


def p_exception_list(p):
    """exception_list : scoped_name
        | exception_list ',' scoped_name"""

    pm = p.parser.pm

    if len(p) == 2:
        throw_args = ThrowArguments(arguments=[pm.find_exception(p, 1, p[1])])
    else:
        throw_args = p[1]
        throw_args.arguments.append(pm.find_exception(p, 3, p[3]))

    p[0] = throw_args


def p_opt_abstract(p):
    """opt_abstract : '=' NUMBER
        | empty"""

    if len(p) == 2:
        p[0] = False
    else:
        if p[2] != 0:
            p.parser.pm.parser_error(p, 2, "'0' expected")

        p[0] = True


def p_opt_signature(p):
    """opt_signature : '[' cpp_type '(' opt_arg_list ')' ']'
        | empty"""

    if len(p) == 7:
        p[0] = Signature(args=p[4], result=p[2])
    else:
        p[0] = None


# The map of operators to dunder method names.
_OPERATOR_NAME_MAP = {
    '+':    '__add__',
    '-':    '__sub__',
    '*':    '__mul__',
    '/':    '__truediv__',
    '%':    '__mod__',
    '&':    '__and__',
    '|':    '__or__',
    '^':    '__xor__',
    '<<':   '__lshift__',
    '>>':   '__rshift__',
    '+=':   '__iadd__',
    '-=':   '__isub__',
    '*=':   '__imul__',
    '/=':   '__itruediv__',
    '%=':   '__imod__',
    '&=':   '__iand__',
    '|=':   '__ior__',
    '^=':   '__ixor__',
    '<<=':  '__ilshift__',
    '>>=':  '__irshift__',
    '~':    '__invert__',
    '()':   '__call__',
    '[]':   '__getitem__',
    '<':    '__lt__',
    '<=':   '__le__',
    '==':   '__eq__',
    '!=':   '__ne__',
    '>':    '__gt__',
    '>=':   '__ge__',
}

def p_operator_name(p):
    """operator_name : '+'
        | '-'
        | '*'
        | '/'
        | '%'
        | '&'
        | '|'
        | '^'
        | '<' '<'
        | '>' '>'
        | '+' '='
        | '-' '='
        | '*' '='
        | '/' '='
        | '%' '='
        | '&' '='
        | '|' '='
        | '^' '='
        | '<' '<' '='
        | '>' '>' '='
        | '~'
        | '(' ')'
        | '[' ']'
        | '<'
        | '<' '='
        | '=' '='
        | '!' '='
        | '>'
        | '>' '='"""

    p[0] = _OPERATOR_NAME_MAP[''.join(p[1:])]


def p_method_code(p):
    """method_code : MethodCode CODE_BLOCK
        | empty"""

    p[0] = p[2] if len(p) == 3 else None


def p_premethod_code(p):
    """premethod_code : PreMethodCode CODE_BLOCK
        | empty"""

    p[0] = p[2] if len(p) == 3 else None


def p_virtual_call_code(p):
    """virtual_call_code : VirtualCallCode CODE_BLOCK
        | empty"""

    p[0] = p[2] if len(p) == 3 else None


def p_virtual_catcher_code(p):
    """virtual_catcher_code : VirtualCatcherCode CODE_BLOCK
        | empty"""

    p[0] = p[2] if len(p) == 3 else None


# C++ namespaces. #############################################################

# The namespace annotations.
_NAMESPACE_ANNOTATIONS = (
    'PyQtNoQMetaObject',
)


def p_namespace_decl(p):
    "namespace_decl : namespace namespace_head opt_namespace_body ';'"

    pm = p.parser.pm

    if pm.skipping:
        return

    if pm.in_main_module:
        pm.scope.iface_file.cpp_name.used = True
        pm.scope.py_name.used = True

    pm.pop_scope()


def p_namespace_head(p):
    "namespace_head : scoped_name opt_annos"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.cpp_only(p, 1, "namespaces")
    pm.check_annotations(p, 2, "namespace", _NAMESPACE_ANNOTATIONS)

    namespace = pm.new_class(p, 1, IfaceFileType.NAMESPACE,
            normalised_scoped_name(p[1], pm.scope))

    namespace.pyqt_no_qmetaobject = p[2].get('PyQtNoQMetaObject', False)

    pm.push_scope(namespace)

def p_opt_namespace_body(p):
    """opt_namespace_body : '{' namespace_body '}'
        | empty"""


def p_namespace_body(p):
    """namespace_body : namespace_statement
        | namespace_body namespace_statement"""


# C/C++ typedefs. #############################################################

# The typedef annotations.
_TYPEDEF_ANNOTATIONS = (
    'Capsule',
    'Encoding',
    'NoTypeName',
    'PyInt',
    'PyName',
    'TypeHint',
    'TypeHintIn',
    'TypeHintOut',
)


def p_typedef_decl(p):
    """typedef_decl : typedef cpp_type NAME opt_annos ';' opt_docstring
        | typedef cpp_type '(' '*' NAME ')' '(' cpp_types ')' opt_annos ';' opt_docstring"""

    pm = p.parser.pm

    if pm.skipping:
        return

    if len(p) == 7:
        type = p[2]
        name_symbol = 3
        annos_symbol = 4
        docstring = p[6]
    else:
        type = p[2]
        name_symbol = 5
        annos_symbol = 10
        docstring = p[12]

    cpp_name = p[name_symbol]
    fq_cpp_name = normalised_scoped_name(ScopedName(cpp_name), pm.scope)
    annotations = p[annos_symbol]

    pm.check_annotations(p, annos_symbol, "typedef", _TYPEDEF_ANNOTATIONS)
    pm.apply_type_annotations(p, annos_symbol, type, annotations)

    no_type_name = annotations.get('NoTypeName', False)

    # See if we are instantiating a class template.
    if type.type is ArgumentType.TEMPLATE:
        instantiated = pm.instantiate_class_template(p, name_symbol,
                fq_cpp_name, type.definition,
                pm.get_py_name(cpp_name, annotations), no_type_name, docstring)

        if instantiated:
            return

    # Handle any 'Capsule' annotation.
    if 'Capsule' in annotations:
        # Make sure the type is 'void *'.
        if type.type is ArgumentType.VOID and len(type.derefs) == 1 and not type.is_const and not type.is_reference:
            type.type = ArgumentType.CAPSULE
            type.derefs = []
        else:
            pm.parser_error(p, annos_symbol,
                    "/Capsule/ can only be specified for a void* type")

    pm.add_typedef(p, name_symbol,
            WrappedTypedef(fq_cpp_name, pm.module_state.module, pm.scope, type,
                    no_type_name=no_type_name))


# C/C++ unions. ###############################################################

# The union annotations.
_UNION_ANNOTATIONS = (
    'AllowNone',
    'DelayDtor',
    'Deprecated',
    'ExportDerived',
    'External',
    'FileExtension',
    'Metatype',
    'Mixin',
    'NoDefaultCtors',
    'NoTypeHint',
    'PyName',
    'PyQtFlags',
    'PyQtFlagsEnums',
    'PyQtInterface',
    'PyQtNoQMetaObject',
    'Supertype',
    'TypeHint',
    'TypeHintIn',
    'TypeHintOut',
    'TypeHintValue',
)


def p_union_decl(p):
    "union_decl : union union_head opt_class_definition ';'"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.complete_class(p, 2, p[2], p[3])


def p_union_head(p):
    "union_head : scoped_name opt_annos"

    pm = p.parser.pm

    if pm.skipping:
        return

    pm.check_annotations(p, 2, "union", _UNION_ANNOTATIONS)

    pm.define_class(p, 1, ClassKey.UNION, p[1], p[2])

    # Return the annotations.
    p[0] = p[2]



# C/C++ variables. ############################################################

# The variable annotations.
_VARIABLE_ANNOTATIONS = (
    'Encoding',
    'NoSetter',
    'NoTypeHint',
    'PyInt',
    'PyName',
    'TypeHint',
)


def p_variable(p):
    "variable : cpp_type NAME opt_annos variable_body ';'"

    pm = p.parser.pm

    if pm.skipping:
        return

    type = p[1]
    cpp_name = p[2]
    annos_symbol = 3
    body = p[4]

    annotations = p[annos_symbol]

    pm.check_annotations(p, annos_symbol, "variable", _VARIABLE_ANNOTATIONS)
    pm.apply_type_annotations(p, annos_symbol, type, annotations)

    py_name = cached_name(pm.spec, pm.get_py_name(cpp_name, p[3]))

    if pm.in_main_module:
        py_name.used = True

    fq_cpp_name = normalised_scoped_name(cpp_name, pm.scope)

    variable = WrappedVariable(fq_cpp_name, pm.module_state.module, py_name,
            pm.scope, type)

    variable.no_setter = annotations.get('NoSetter', False)
    variable.no_type_hint = annotations.get('NoTypeHint', False)

    variable.access_code = body.get('%AccessCode')
    variable.get_code = body.get('%GetCode')
    variable.set_code = body.get('%SetCode')

    pm.validate_variable(p, 1, variable)

    pm.spec.variables.append(variable)

    p[0] = variable


def p_variable_body(p):
    """variable_body : '{' variable_body_directives '}'
        | empty"""

    p[0] = p[2] if len(p) == 4 else {}


def p_variable_body_directives(p):
    """variable_body_directives : variable_body_directive
        | variable_body_directives variable_body_directive"""

    if len(p) == 3:
        p[1].update(p[2])

    p[0] = p[1]


def p_variable_body_directive(p):
    """variable_body_directive : if_start
        | if_end
        | AccessCode CODE_BLOCK
        | GetCode CODE_BLOCK
        | SetCode CODE_BLOCK"""

    p[0] = {} if len(p) == 2 or p.parser.pm.skipping else {p[1]: p[2]}


# Annotations. ################################################################

def p_opt_annos(p):
    """opt_annos : '/' annotations '/'
        | empty"""

    p[0] = p[2] if len(p) == 4 else {}


def p_annotations(p):
    """annotations : annotation
        | annotations ',' annotation"""

    if len(p) == 4:
        p[1].update(p[3])

    p[0] = p[1]


def p_annotation(p):
    """annotation : NAME
        | NAME '=' annotation_value"""

    value = None if len(p) == 2 else p[3]
    value = p.parser.pm.validate_annotation(p, 1, value)

    p[0] = {p[1]: value}


def p_annotation_value(p):
    """annotation_value : dotted_name
        | STRING
        | NUMBER"""

    p[0] = p[1]


# A (possibly) scoped name. ###################################################

# Note that support for the global scope (ie. a leading '::') requires the
# following to eliminate a shift/reduce conflict.  However that still leaves
# one other shift/reduce conflict which I can't get rid of.
precedence = (
    ('left', 'SCOPE'),
)


def p_scoped_name(p):
    """scoped_name : SCOPE relative_scoped_name
        | relative_scoped_name"""

    if len(p) == 3:
        p.parser.pm.cpp_only(p, 1, "scoped names")
        p[2].make_absolute()
        p[0] = p[2]
    else:
        p[0] = p[1]


def p_relative_scoped_name(p):
    """relative_scoped_name : NAME
        | relative_scoped_name SCOPE NAME"""

    if len(p) == 2:
        p[0] = ScopedName(p[1])
    else:
        p.parser.pm.cpp_only(p, 2, "scoped names")
        p[1].append(p[3])
        p[0] = p[1]


# The remaining value productions. ############################################

def p_template_decl(p):
    "template_decl : template '<' cpp_types '>'"

    pm = p.parser.pm

    if pm.skipping:
        return

    p[0] = Signature(args=p[3])

    pm.parsing_template = True


def p_bool_value(p):
    """bool_value : true
        | True
        | false
        | False"""

    p[0] = (p[1].lower() == 'true')


def p_dotted_name(p):
    """dotted_name : NAME
        | DOTTED_NAME"""

    p[0] = DottedName(p[1])


def p_file_path(p):
    """file_path : NAME
        | DOTTED_NAME
        | FILE_PATH"""

    p[0] = p[1]


def p_empty(p):
    "empty :"

    p[0] = None


def p_opt_const(p):
    """opt_const : const
        | empty"""

    p[0] = p[1] is not None


def p_opt_docstring(p):
    """opt_docstring : docstring
        | empty"""

    p[0] = p[1]


def p_opt_final(p):
    """opt_final : final
        | empty"""

    p[0] = p[1] is not None


def p_opt_name(p):
    """opt_name : NAME
        | empty"""

    p[0] = p[1]


def p_opt_virtual(p):
    """opt_virtual : virtual
        | empty"""

    pm = p.parser.pm

    if pm.skipping:
        return

    if p[1] is not None:
        pm.parsing_virtual = True


def p_ored_qualifiers(p):
    """ored_qualifiers : NAME
        | '!' NAME
        | ored_qualifiers LOGICAL_OR NAME
        | ored_qualifiers LOGICAL_OR '!' NAME"""

    pm = p.parser.pm

    nr_symbols = len(p)

    if nr_symbols == 2:
        p[0] = pm.evaluate_feature_or_platform(p, 1)
    elif nr_symbols == 3:
        p[0] = pm.evaluate_feature_or_platform(p, 2, inverted=True)
    elif nr_symbols == 4:
        p[0] = p[1] or pm.evaluate_feature_or_platform(p, 3)
    else:
        p[0] = p[1] or pm.evaluate_feature_or_platform(p, 4, inverted=True)


def p_qualifier_list(p):
    """qualifier_list : NAME
        | qualifier_list NAME"""

    if len(p) == 2:
        value = [p[1]]
    else:
        value = p[1]
        value.append(p[2])

    p[0] = value


def p_qualifiers(p):
    """qualifiers : ored_qualifiers
        | opt_name '-' opt_name"""

    if len(p) == 2:
        p[0] = p[1]
    else:
        p[0] = p.parser.pm.evaluate_timeline(p, 1, 3)
