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


from .base_formatter import BaseFormatter
from .utils import format_scoped_py_name


class ScopedFormatter(BaseFormatter):
    """ A base class for formatters of objects that can be contained by a
    scope.
    """

    def __init__(self, spec, object, scope):
        """ Initialise the object. """

        super().__init__(spec, object)

        self.scope = scope

    @property
    def cpp_scope(self):
        """ The C++ scope as a string. """

        if self.scope is None:
            return ''

        return str(self.scope.iface_file.fq_cpp_name) + '::'

    @property
    def fq_cpp_name(self):
        """ The fully qualified C++ name. """

        return self.cpp_scope + self.object.cpp_name

    @property
    def fq_py_name(self):
        """ The fully qualified Python name. """

        return format_scoped_py_name(self.scope, self.object.py_name.name)


class EmbeddedScopeFormatter(ScopedFormatter):
    """ A base class for formatters of objects that have a reference to their
    scope.
    """

    def __init__(self, spec, object):
        """ Initialise the object. """

        super().__init__(spec, object, object.scope)
