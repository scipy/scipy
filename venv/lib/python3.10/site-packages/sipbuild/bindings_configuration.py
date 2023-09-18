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


import os

from .exceptions import UserFileException, UserParseException
from .module import resolve_abi_version
from .toml import toml_load


def get_bindings_configuration(abi_major, sip_file, sip_include_dirs):
    """ Get the configuration of a set of bindings. """

    # We make no assumption about the name of the .sip file but we assume that
    # the directory it is in is the name of the bindings.
    bindings_name = os.path.basename(os.path.dirname(sip_file))

    # See if there is a .toml file.
    for sip_dir in sip_include_dirs:
        toml_file = os.path.join(sip_dir, bindings_name,
                bindings_name + '.toml')

        if os.path.isfile(toml_file):
            break
    else:
        return [], []

    # Read the configuration.
    try:
        cfg = toml_load(toml_file)
    except Exception as e:
        raise UserParseException(toml_file, detail=str(e))

    # Check the ABI version is compatible.
    cfg_abi_version = cfg.get('sip-abi-version')
    if cfg_abi_version is None or not isinstance(cfg_abi_version, str):
        raise UserFileException(toml_file,
                "'sip-abi-version' must be specified as a string")

    cfg_abi_major = int(resolve_abi_version(cfg_abi_version).split('.')[0])

    if cfg_abi_major != abi_major:
        raise UserFileException(toml_file,
                "'{0}' was built against ABI v{1} but this module is being "
                        "built against ABI v{2}".format(bindings_name,
                                cfg_abi_major, abi_major))

    # Return the tags and disabled features.
    return (_get_string_list(toml_file, cfg, 'module-tags'),
            _get_string_list(toml_file, cfg, 'module-disabled-features'))


def _get_string_list(toml_file, cfg, name):
    """ Get an option from the configuration and check it is a list of strings.
    """

    option_list = cfg.get(name)
    if option_list is None:
        option_list = list()
    elif not isinstance(option_list, list):
        raise UserFileException(toml_file, "'{0}' must be a list".format(name))

    for option in option_list:
        if not isinstance(option, str):
            raise UserFileException(toml_file,
                    "elements of '{0}' must be strings".format(name))

    return option_list
