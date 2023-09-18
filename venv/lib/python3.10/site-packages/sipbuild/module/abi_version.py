# Copyright (c) 2021, Riverbank Computing Limited
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


import os

from ..exceptions import UserException


# The directory containing the different module implementations.
_module_source_dir = os.path.join(os.path.dirname(__file__), 'source')


def get_module_source_dir(abi_major_version):
    """ Return the name of the directory containing the latest source of the
    sip module that implements the given ABI version.
    """

    return os.path.join(_module_source_dir, abi_major_version)


def get_sip_module_version(abi_major_version):
    """ Return the version number of the latest implementation of the sip
    module with the given ABI as a string.
    """

    abi_minor_version = patch_version = None

    # Read the version from the header file shared with the code generator.
    with open(os.path.join(get_module_source_dir(abi_major_version), 'sip.h.in')) as vf:
        for line in vf:
            parts = line.strip().split()
            if len(parts) == 3 and parts[0] == '#define':
                name = parts[1]
                value = parts[2]

                if name == 'SIP_ABI_MINOR_VERSION':
                    abi_minor_version = value
                elif name == 'SIP_MODULE_PATCH_VERSION':
                    patch_version = value

    # These are internal errors and should never happen.
    if abi_minor_version is None:
        raise ValueError(
                f"'SIP_ABI_MINOR_VERSION' not found for ABI {abi_major_version}")

    if patch_version is None:
        raise ValueError(
                f"'SIP_MODULE_PATCH_VERSION' not found for ABI {abi_major_version}")

    return f'{abi_major_version}.{abi_minor_version}.{patch_version}'


def resolve_abi_version(abi_version, module=True):
    """ Return a valid ABI version or the latest if none was given. """

    # Get the major and minimum minor version of what we need.
    if abi_version:
        parts = abi_version.split('.')
        abi_major_version = parts[0]

        if not os.path.isdir(get_module_source_dir(abi_major_version)):
            raise UserException(
                    f"'{abi_version}' is not a supported ABI version")

        if len(parts) == 1:
            minimum_minor_version = 0
        else:
            try:
                minimum_minor_version = int(parts[1])
            except ValueError:
                minimum_minor_version = None

            if len(parts) > 2 or minimum_minor_version is None:
                raise UserException(
                        f"'{abi_version}' is not a valid ABI version")
    else:
        abi_major_version = sorted(os.listdir(_module_source_dir), key=int)[-1]
        minimum_minor_version = 0

    # Get the minor version of what we actually have.
    module_version = get_sip_module_version(abi_major_version)
    _, abi_minor_version, _ = module_version.split('.')

    # Check we meet the minimum requirement.
    if int(abi_minor_version) < minimum_minor_version:
        raise UserException(f"'{abi_version}' is not a supported ABI version")

    if module:
        # Return the module's version.
        return f'{abi_major_version}.{abi_minor_version}'

    # Return the required version.
    return f'{abi_major_version}.{minimum_minor_version}'
