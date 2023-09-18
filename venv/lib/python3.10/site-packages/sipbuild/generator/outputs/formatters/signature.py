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


from ...scoped_name import STRIP_NONE
from ...specification import ArgumentType, ArrayArgument

from .argument import ArgumentFormatter
from .base_formatter import BaseFormatter


class SignatureFormatter(BaseFormatter):
    """ This creates various string representations of a signature. """

    def cpp_arguments(self, *, scope=None, strip=STRIP_NONE, make_public=False,
            as_xml=False):
        """ Return the C++ representation of the signature arguments. """

        args = [ArgumentFormatter(self.spec, arg).cpp_type(scope=scope,
                strip=strip, make_public=make_public, as_xml=as_xml)
                for arg in self.object.args]

        # Note the lack of separating space (although this may be for XML only
        # to be consistent with previous implementations).
        return ','.join(args)

    @property
    def py_arguments(self):
        """ The Python representation of the signature arguments. """

        args = []

        for arg in self.object.args:
            if arg.array is not ArrayArgument.ARRAY_SIZE and arg.is_in:
                args.append(
                        ArgumentFormatter(self.spec,
                                arg).as_py_type(default_value=True))

        return ', '.join(args)

    @property
    def py_results(self):
        """ The Python representation of the signature results. """

        sig = self.object

        results = []

        if sig.result is not None:
            if sig.result.type is not ArgumentType.VOID or len(sig.result.derefs) != 0:
                results.append(ArgumentFormatter(self.spec,
                        sig.result).as_py_type())

        for arg in sig.args:
            if arg.is_out:
                results.append(ArgumentFormatter(self.spec, arg).as_py_type())

        return ', '.join(results)
