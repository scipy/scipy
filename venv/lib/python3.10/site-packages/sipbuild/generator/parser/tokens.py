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


from ply.lex import TOKEN

from ..specification import CodeBlock


# The lexer states.
states = (
    ('code',            'exclusive'),
    ('ccomment',        'exclusive'),
    ('directive',       'inclusive'),
    ('needeol',         'inclusive'),
)


# The single character tokens.
literals = '(){}.,;:=!-+*/&|~<>[]%^'


# The non-code directives.
directives = {
    'AutoPyName', 'CompositeModule', 'DefaultDocstringFormat',
    'DefaultDocstringSignature', 'DefaultEncoding', 'DefaultMetatype',
    'DefaultSupertype', 'End', 'Exception', 'Feature', 'HideNamespace', 'If',
    'Import', 'Include', 'License', 'MappedType', 'Module', 'Platforms',
    'Property', 'Timeline',

    # Remove in SIP v7.
    'Plugin',
}


# The code directives.
code_directives = {
    'AccessCode', 'BIGetBufferCode', 'BIReleaseBufferCode',
    'ConvertFromTypeCode', 'ConvertToSubClassCode', 'ConvertToTypeCode',
    'Copying', 'Docstring', 'ExportedHeaderCode', 'ExportedTypeHintCode',
    'Extract', 'FinalisationCode', 'GCClearCode', 'GCTraverseCode', 'GetCode',
    'InitialisationCode', 'InstanceCode', 'MethodCode', 'ModuleCode',
    'ModuleHeaderCode', 'PickleCode', 'PostInitialisationCode',
    'PreInitialisationCode', 'PreMethodCode', 'RaiseCode', 'ReleaseCode',
    'SetCode', 'TypeCode', 'TypeHeaderCode', 'TypeHintCode', 'UnitCode',
    'UnitPostIncludeCode', 'VirtualCallCode', 'VirtualCatcherCode',
    'VirtualErrorHandler',

    # Remove in SIP v7.
    'BIGetCharBufferCode', 'BIGetReadBufferCode', 'BIGetSegCountCode',
    'BIGetWriteBufferCode',
}


# The plain keywords.
keywords = {
    'bool', 'char', 'class', 'const', 'double', 'enum', 'explicit', 'false',
    'final', 'float', 'int', 'long', 'namespace', 'noexcept', 'NULL',
    'operator', 'private', 'protected', 'public', 'Py_hash_t', 'Py_ssize_t',
    'Q_SIGNAL', 'Q_SIGNALS', 'Q_SLOT', 'Q_SLOTS', 'short', 'signals', 'signed',
    'SIP_PYBUFFER', 'SIP_PYCALLABLE', 'SIP_PYDICT', 'SIP_PYENUM', 'SIP_PYLIST',
    'SIP_PYOBJECT', 'SIP_PYSLICE', 'SIP_PYTUPLE', 'SIP_PYTYPE', 'size_t',
    'slots', 'static', 'struct', 'template', 'throw', 'true', 'typedef',
    'union', 'unsigned', 'virtual', 'void', 'wchar_t',

    # Remove in SIP v7.
    'SIP_SSIZE_T',
}


# The directive keywords.
directive_keywords = {
    'all_raise_py_exception', 'call_super_init', 'default_VirtualErrorHandler',
    'False', 'format', 'get', 'id', 'keyword_arguments', 'language',
    'licensee', 'name', 'optional', 'order', 'remove_leading', 'set',
    'signature', 'timestamp', 'True', 'type', 'py_ssize_t_clean',
    'use_argument_names', 'use_limited_api',
}


# The lexer tokens.
tokens = [
    'CODE_BLOCK', 'DOTTED_NAME', 'ELLIPSIS', 'EOF', 'EOL', 'FILE_PATH',
    'LOGICAL_OR', 'NAME', 'NUMBER', 'QUOTED_CHAR', 'REAL', 'SCOPE', 'STRING',
]

tokens.extend(directives)
tokens.extend(code_directives)
tokens.extend(keywords)
tokens.extend(directive_keywords)


# Handle EOF.
def t_eof(t):

    try:
        t.lexer.pm.pop_file()
    except IndexError:
        return None

    # Return an explicit EOF token.  This stops the parser looking too far into
    # the popped file.
    t.type = 'EOF'

    return t


# Handle errors.
def t_ANY_error(t):

    t.lexer.pm.lexer_error(t, "'{0}' is unexpected".format(t.value[0]))
    t.lexer.skip(1)


# Ignore whitespace except when reading code blocks.
t_ANY_ignore = ' \t\r'
t_code_ignore = ''


# Handle newlines outside of code blocks and comments.
def t_newline(t):
    r'\n'

    lexer = t.lexer
    pm = lexer.pm

    # Maintain the line number.
    lexer.lineno += 1

    # Enter the 'code' state if we are at the end of a code directive name and
    # arguments.
    if pm.code_block is not None and pm.paren_depth == 0:
        pm.code_block.line_nr = lexer.lineno
        pm.set_lexer_state('code')


# Maintain the parenthesis depth.
def t_LPAREN(t):
    r'\('

    t.lexer.pm.paren_depth += 1
    t.type = '('

    return t


def t_RPAREN(t):
    r'\)'

    t.lexer.pm.paren_depth -= 1
    t.type = ')'

    return t


# Handle directives.
def t_DIRECTIVE(t):
    r'%[a-zA-Z][a-zA-Z]*'

    # The name of the directive is used as its type.
    name = t.value[t.value.index('%') + 1:]

    if name in code_directives:
        t.lexer.pm.code_block = CodeBlock(t.lexer.pm.raw_sip_file)
        t.type = name
    elif name in directives:
        t.type = name

    return t


# Handle the %End of a code directive.
def t_code_END(t):
    r'%End'

    t.type = 'CODE_BLOCK'
    t.value = t.lexer.pm.code_block
    t.lexer.pm.code_block = None
    t.lexer.begin('INITIAL')

    return t


# Handle a newline when an end-of-line needs to be reported to the parser.
def t_needeol_newline(t):
    r'\n'

    # Maintain the line number.
    t.lexer.lineno += 1

    t.lexer.pm.set_lexer_state()
    t.type = 'EOL'

    return t


# Handle a newline in a code directive.
def t_code_newline(t):
    r'\n'

    # Maintain the line number.
    t.lexer.lineno += 1

    t.lexer.pm.code_block.text += t.value

    # Discard the token.
    return None


# Handle a character in a code directive.
def t_code_CH(t):
    r'.'

    t.lexer.pm.code_block.text += t.value

    # Discard the token.
    return None


# Handle keywords, ellipsis, names, dotted name and file paths.
ambiguous = r'[._A-Za-z][._/A-Za-z\d\-]*[._A-Za-z\d]'

@TOKEN(ambiguous)
def t_AMBIGUOUS(t):

    t.type = t.lexer.pm.disambiguate_token(t.value, keywords)

    return t


# Handle directive keywords (ie. keywords that are only recognised in the
# context of a directive), ellipsis, names, dotted name and file paths.
@TOKEN(ambiguous)
def t_directive_AMBIGUOUS(t):

    t.type = t.lexer.pm.disambiguate_token(t.value, directive_keywords)

    return t


# Handle a C++-style comment.
def t_CPPCOMMENT(t):
    r'//.*'

    # Discard the token.
    return None


# Handle the start of a C-style comment.
def t_COMMENTSTART(t):
    r'/\*'

    t.lexer.push_state('ccomment')

    # Discard the token.
    return None


# Handle the end of a C-style comment.
def t_ccomment_COMMENTEND(t):
    r'\*/'

    t.lexer.pop_state()

    # Discard the token.
    return None


# Handle a newline in a C-style comment.
def t_ccomment_newline(t):
    r'\n'

    # Maintain the line number.
    t.lexer.lineno += 1

    # Discard the token.
    return None


# Handle the content of a C-style comment.
def t_ccomment_CH(t):
    r'.'

    # Maintain the line number.
    if t.value == '\n':
        t.lexer.lineno += 1

    # Discard the token.
    return None


# Handle an unsigned hexadecimal number.
def t_HEXNUMBER(t):
    r'0x[\da-fA-F]+'

    t.type = 'NUMBER'
    t.value = int(t.value, base=16)

    return t


# Handle a number.
def t_NUMBER(t):
    r'-?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?[fFlLuU]?'

    # Remove any suffix character.
    value = t.value
    if not value[-1].isdigit():
        value = value[:-1]

    try:
        t.type = 'NUMBER'
        t.value = int(value)
    except ValueError:
        t.type = 'REAL'
        t.value = float(value)

    return t


# Handle a double-quoted string.
def t_STRING(t):
    r'"(\\.|[^\\"])*"'

    # Strip the quotes and handle any standard escape characters.
    value = t.value.strip('"')
    value = value.replace(r'\n', '\n')
    value = value.replace(r'\r', '\r')
    value = value.replace(r'\t', '\t')
    value = value.replace(r'\\', '\\')

    t.type = 'STRING'
    t.value = value

    return t


# Handle a single-quoted hex encoded character.
def t_QHEXCH(t):
    r"'\\x[\da-fA-F]+'"

    t.type = 'QUOTED_CHAR'
    t.value = chr(int(t.value.strip("'")[2:], base=16))

    return t


# Handle a single-quoted character.
def t_QCH(t):
    r"'[^'\n]*['\n]"

    # Make sure these is only one quoted character.  If not then report the
    # error and carry on with a fudged value.
    n_ch = len(t.value)

    if n_ch != 3:
        t.lexer.pm.lexer_error(t,
                "exactly one character expected between single quotes")

        if n_ch == 0:
            t.value = '?'

    t.type = 'QUOTED_CHAR'
    t.value = t.value[1]

    return t


# The remaining trivial token definitions.
t_LOGICAL_OR = r'\|\|'
t_SCOPE = r'::'

# We only deal with a single character as everything else is handled by
# AMBIGUOUS.
t_NAME = r'[_A-Za-z]'
