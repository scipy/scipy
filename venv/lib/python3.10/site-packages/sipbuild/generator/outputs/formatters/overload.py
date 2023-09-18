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


from ...specification import PySlot

from .scoped import ScopedFormatter
from .utils import format_scoped_py_name


# Map the Python slot types to C++.
_PYSLOT_CPP_MAP = {
    PySlot.ADD: '+',
    PySlot.SUB: '-',
    PySlot.MUL: '*',
    PySlot.TRUEDIV: '/',
    PySlot.MOD: '%',
    PySlot.AND: '&',
    PySlot.OR: '|',
    PySlot.XOR: '^',
    PySlot.LSHIFT: '<<',
    PySlot.RSHIFT: '>>',
    PySlot.IADD: '+=',
    PySlot.ISUB: '-=',
    PySlot.IMUL: '*=',
    PySlot.ITRUEDIV: '/=',
    PySlot.IMOD: '%=',
    PySlot.IAND: '&=',
    PySlot.IOR: '|=',
    PySlot.IXOR: '^=',
    PySlot.ILSHIFT: '<<=',
    PySlot.IRSHIFT: '>>=',
    PySlot.INVERT: '~',
    PySlot.CALL: '()',
    PySlot.GETITEM: '[]',
    PySlot.LT: '<',
    PySlot.LE: '<=',
    PySlot.EQ: '==',
    PySlot.NE: '!=',
    PySlot.GT: '>',
    PySlot.GE: '>=',
}


class OverloadFormatter(ScopedFormatter):
    """ This creates various string representations of an overload. """

    @property
    def fq_cpp_name(self):
        """ The fully qualified C++ name. """

        try:
            cpp_op = _PYSLOT_CPP_MAP[self.object.common.py_slot]
        except KeyError:
            return super().fq_cpp_name

        return self.cpp_scope + 'operator' + cpp_op

    @property
    def fq_py_name(self):
        """ The fully qualified Python name. """

        return format_scoped_py_name(self.scope,
                self.object.common.py_name.name)
