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


from enum import IntEnum

from ..specification import AccessSpecifier

from .formatters import (ClassFormatter, EnumFormatter, OverloadFormatter,
        SignatureFormatter, VariableFormatter)


class IconNumber(IntEnum):
    """ The numbers of the different icons.  The values are those used by the
    eric IDE.
    """

    CLASS = 1
    METHOD = 4
    VARIABLE = 7
    ENUM = 10


def output_api(spec, api_filename):
    """ Output a QScintilla API file. """

    with open(api_filename, 'w') as af:
        module = spec.module

        _enums(af, spec, module)
        _variables(af, spec, module)

        for overload in module.overloads:
            if overload.common.module is module and overload.common.py_slot is None:
                _overload(af, spec, module, overload)

        for klass in spec.classes:
            if klass.iface_file.module is module and not klass.external:
                _enums(af, spec, module, scope=klass)
                _variables(af, spec, module, scope=klass)

                for ctor in klass.ctors:
                    if ctor.access_specifier is not AccessSpecifier.PRIVATE:
                        _ctor(af, spec, module, ctor, klass)

                for overload in klass.overloads:
                    if overload.access_specifier is not AccessSpecifier.PRIVATE and overload.common.py_slot is None:
                        _overload(af, spec, module, overload, scope=klass)


def _ctor(af, spec, module, ctor, scope):
    """ Generate an API ctor. """

    py_class = module.py_name + '.' + ClassFormatter(spec, scope).fq_py_name
    py_arguments = SignatureFormatter(spec, ctor.py_signature).py_arguments

    # Do the callable type form.
    af.write(f'{py_class}?{IconNumber.CLASS}({py_arguments})\n')

    # Do the call __init__ form.
    if py_arguments:
        py_arguments = 'self, ' + py_arguments
    else:
        py_arguments = 'self'

    af.write(f'{py_class}.__init__?{IconNumber.CLASS}({py_arguments})\n')


def _enums(af, spec, module, scope=None):
    """ Generate the APIs for all the enums in a scope. """

    for enum in spec.enums:
        if enum.module is module and enum.scope is scope:
            formatter = EnumFormatter(spec, enum)

            if enum.py_name is not None:
                af.write(f'{module.py_name}.{formatter.fq_py_name}?{IconNumber.ENUM}\n')

            for member_s in formatter.fq_py_member_names:
                af.write(f'{member_s}?{IconNumber.ENUM}\n')


def _variables(af, spec, module, scope=None):
    """ Generate the APIs for all the variables in a scope. """

    for variable in spec.variables:
        if variable.module is module and variable.scope is scope:
            formatter = VariableFormatter(spec, variable)

            af.write(f'{module.py_name}.{formatter.fq_py_name}?{IconNumber.VARIABLE}\n')


def _overload(af, spec, module, overload, scope=None):
    """ Generate a single API overload. """

    sig_formatter = SignatureFormatter(spec, overload.py_signature)

    s = module.py_name + '.' + OverloadFormatter(spec, overload, scope).fq_py_name

    s += f'?{IconNumber.METHOD}({sig_formatter.py_arguments})'

    results = sig_formatter.py_results
    if results:
        s += ' -> '

        if ', ' in results:
            s += '(' + results + ')'
        else:
            s += results

    s += '\n'

    af.write(s)
