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


from ...specification import IfaceFileType

from .scoped import EmbeddedScopeFormatter
from .utils import format_scoped_py_name, iface_is_defined


class EnumFormatter(EmbeddedScopeFormatter):
    """ This creates various string representations of an enum. """

    @property
    def fq_py_member_names(self):
        """ An iterator over the fully qualified Python names of the members of
        the enum.
        """

        enum = self.object

        enum_name = enum.module.py_name + '.'

        if enum.py_name is not None:
            enum_name += format_scoped_py_name(self.scope, enum.py_name.name)
            enum_name += '.'

        for member in enum.members:
            yield enum_name + member.py_name.name

    def member_as_rest_ref(self, member):
        """ Return the fully qualified Python name of a member as a reST
        reference.
        """

        enum = self.object
        module_name = enum.module.fq_py_name.name

        if enum.py_name is None:
            member_name = format_scoped_py_name(self.scope,
                    member.py_name.name)

            return f':sip:ref:`~{module_name}.{member_name}`'

        enum_name = format_scoped_py_name(self.scope, enum.py_name.name)
        member_name = member.py_name.name

        return f':sip:ref:`~{module_name}.{enum_name}.{member_name}`'

    def as_rest_ref(self):
        """ Return the fully qualified Python name as a reST reference. """

        enum = self.object
        module_name = enum.module.fq_py_name.name
        enum_name = format_scoped_py_name(self.scope, enum.py_name.name)

        return f':sip:ref:`~{module_name}.{enum_name}`'

    def as_type_hint(self, module, defined):
        """ Return the type hint. """

        enum = self.object

        if self.scope is None:
            # Global enums are defined early on.
            is_defined = True
        else:
            scope_iface = self.scope.iface_file
            outer_scope = self.scope.scope if scope_iface.type is IfaceFileType.CLASS else None

            is_defined = iface_is_defined(scope_iface, outer_scope, module,
                    defined)

        quote = '' if is_defined else "'"

        # Include the module name if it is not the current one.
        module_name = enum.module.py_name + '.' if enum.module is not module else ''

        return f'{quote}{module_name}{self.fq_py_name}{quote}'
