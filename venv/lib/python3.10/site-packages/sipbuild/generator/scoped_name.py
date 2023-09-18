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



# Specify how scopes should be stripped.  Any other value is the number of
# scopes to strip.
STRIP_NONE = 0      # Don't strip any scopes.
STRIP_GLOBAL = -1   # Strip the global scope.


class ScopedName:
    """ Encapsulate a scoped name. """

    def __init__(self, name):
        """ Initialise the scoped name. """

        if isinstance(name, list):
            # This use is internal to this class.
            self._name = name
        elif isinstance(name, ScopedName):
            # Make a copy of another ScopedName object.
            self._name = list(name._name)
        else:
            # Create a simple name.
            self._name = [name]

    def __eq__(self, other):
        """ Compare with another scoped name for equality. """

        if isinstance(other, ScopedName):
            return self._name == other._name

        return NotImplemented

    def __delitem__(self, index):
        """ Remove the requested name. """

        del self._name[self._normalised_index(index)]

    def __getitem__(self, index):
        """ Get the requested name. """

        return self._name[self._normalised_index(index)]

    def __len__(self):
        """ Return the length of the name, ie. the number of indivual names.
        """

        nr_names = len(self._name)

        if self.is_absolute:
            nr_names -= 1

        return nr_names

    def __lt__(self, other):
        """ Compare with another scoped name to allow lists to be sorted. """

        if isinstance(other, ScopedName):
            return self._name < other._name

        return NotImplemented

    def __setitem__(self, index, name):
        """ Set the requested name. """

        self._name[self._normalised_index(index)] = name

    def __str__(self):
        """ Return the C++ string representation. """

        # The special treatment is a hack to simplify how the result is used.
        # We ignore any leading '::' and truncate the conversion if we think we
        # are converting an encoded template name (ie. we find a word that
        # starts with a digit).
        s_l = []
        for s in self:
            if s[0].isdigit():
                break

            s_l.append(s)

        return '::'.join(s_l)

    @property
    def absolute(self):
        """ The absolute version of the name. """

        if self.is_absolute:
            return self

        copy = type(self)(self)
        copy.make_absolute()

        return copy

    def append(self, name):
        """ Append a simple name. """

        self._name.append(name)

    @property
    def as_cpp(self):
        """ The C++ representation of the name. """

        cpp = '::'.join(self._name)

        # Precede an explicit global scope with a space to avoid '<::'.
        if cpp.startswith('::'):
            cpp = ' ' + cpp

        return cpp

    @property
    def as_py(self):
        """ The Python representation of the name. """

        start = 1 if self.is_absolute else 0

        return '.'.join(self._name[start:])

    @property
    def as_word(self):
        """ The word representation of the name. """

        start = 1 if self.is_absolute else 0

        return '_'.join(self._name[start:])

    @property
    def base_name(self):
        """ The base name of the scoped name. """

        return self._name[-1]

    def cpp_stripped(self, strip):
        """ Return the C++ representation of the name with leading scopes
        stripped.
        """

        if strip == STRIP_NONE:
            start = 0
        else:
            start = 1 if self.is_absolute else 0

            if strip != STRIP_GLOBAL:
                start += strip

        return '::'.join(self._name[start:])

    @property
    def is_absolute(self):
        """ True if the scoped name is absolute. """

        return self._name[0] == ''

    @property
    def is_simple(self):
        """ Return True if the name is simple, ie. relative and unscoped. """

        return len(self) == 1

    def make_absolute(self):
        """ Make sure the scoped name is absolute. """

        if self._name[0] != '':
            self._name.insert(0, '')

    def matches(self, scoped_name, scope=None):
        """ Return True if a scoped name matches this taking account of an
        optional scope if the scoped name is relative.
        """

        assert self.is_absolute

        # Check for the simple case.
        if scoped_name.is_absolute:
            return self == scoped_name

        # Try each scope, inner to outer.
        while scope is not None:
            fq_name = ScopedName(scope.iface_file.fq_cpp_name)
            fq_name._name.extend(scoped_name._name)

            if self == fq_name:
                return True

            scope = scope.scope

        # Treat the name as if it were absolute.
        return self._name[1:] == scoped_name._name

    @classmethod
    def parse(cls, raw):
        """ Return a ScopedName object by parsing a raw string. """

        return cls(raw.split('::'))

    def prepend(self, scoped_name):
        """ Prepend a scoped name. """

        new_name = list(scoped_name._name)
        new_name.extend(self._name)
        self._name = new_name

    @property
    def scope(self):
        """ The scoped name that is the enclosing scope of this one.  It will
        be None if there isn't one.
        """

        if len(self) == 1:
            return None

        return type(self)(self._name[:-1])

    def _normalised_index(self, index):
        """ Return a normalised index. """

        if not isinstance(index, int):
            raise TypeError("'{}' is an invalid index".format(type(index)))

        if index >= 0 and self.is_absolute:
            index += 1

        return index
