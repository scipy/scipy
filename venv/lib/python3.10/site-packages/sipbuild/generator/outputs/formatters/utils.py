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


def format_copying(copying, comment):
    """ Return a formatted %Copying text. """

    s = comment + ' ' + ''.join([b.text for b in copying]).rstrip().replace('\n', '\n' + comment + ' ')

    if s:
        s = comment + '\n' + s

    return s


def format_scoped_py_name(scope, py_name):
    """ Return a formatted scoped Python name. """

    if scope is None or scope.is_hidden_namespace:
        scope_s = ''
    else:
        scope_s = format_scoped_py_name(scope.scope, None) + scope.py_name.name + '.'

    if py_name is None:
        py_name_s = ''
    else:
        py_name_s = py_name

    return scope_s + py_name_s


def iface_is_defined(iface_file, scope, module, defined):
    """ Return True if a type corresponding to an interface file has been
    defined in the context of a module.
    """

    # A type in another module would have been imported.
    if iface_file.module is not module:
        return True

    if iface_file not in defined:
        return False

    # Check all enclosing scopes have been defined as well.
    while scope is not None:
        if scope.iface_file not in defined:
            return False

        scope = scope.scope

    return True
