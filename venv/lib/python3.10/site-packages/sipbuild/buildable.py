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


import importlib
import os

from .exceptions import UserException
from .installable import Installable
from .py_versions import OLDEST_SUPPORTED_MINOR
from .version import SIP_VERSION_STR


class Buildable:
    """ Encapsulate the components used to build something that can be
    installed.
    """

    def __init__(self, project, name):
        """ Initialise the buildable. """

        self.project = project
        self.name = name
        self.builder_settings = []
        self.installables = []

        self.build_dir = os.path.join(project.build_dir, name)
        os.makedirs(self.build_dir, exist_ok=True)


class BuildableFromSources(Buildable):
    """ Encapsulate the sources used to build an extension module, executable
    etc.
    """

    def __init__(self, project, name, target, *, uses_limited_api=False):
        """ Initialise the buildable. """

        super().__init__(project, name)

        if project.py_debug:
            uses_limited_api = False
        elif uses_limited_api and not project.sip_module:
            raise UserException(
                    "{0} cannot use the limited API without using a shared "
                    "'sip' module".format(name))

        self.target = target
        self.uses_limited_api = uses_limited_api

        self.define_macros = []
        self.sources = []
        self.headers = []
        self.include_dirs = []
        self.libraries = []
        self.library_dirs = []
        self.extra_compile_args = []
        self.extra_link_args = []
        self.extra_objects = []
        self.debug = False

        if self.uses_limited_api:
            self.define_macros.append(
                    'Py_LIMITED_API=0x03{0:02x}0000'.format(
                            OLDEST_SUPPORTED_MINOR))

    def make_names_relative(self):
        """ Make all file and directory names relative to the build directory.
        """

        # Make the file names relative to the build directory.
        self.include_dirs = self._relative_names(self.include_dirs)
        self.headers = self._relative_names(self.headers)
        self.sources = self._relative_names(self.sources)
        self.library_dirs = self._relative_names(self.library_dirs)

    def _relative_names(self, names):
        """ Return a list of times made relative to build directory.  Note that
        we only really do this for cosmetic reasons to simplify what the user
        might see.
        """

        rel_names = []

        for fn in names:
            try:
                common = os.path.commonpath([fn, self.build_dir])
                _, common = os.path.splitdrive(common)

                if len(common) > 1:
                    # Only convert to a relative name if there is at least one
                    # parent directory in common.
                    fn = os.path.relpath(fn, self.build_dir)
            except ValueError:
                # This is most likely to happen if the build directory is on a
                # different Windows drive.
                pass

            rel_names.append(fn)

        return rel_names


class BuildableExecutable(BuildableFromSources):
    """ Encapsulate the sources used to build an executable. """


class BuildableModule(BuildableFromSources):
    """ Encapsulate the sources used to build an extension module. """

    def __init__(self, project, name, fq_name, *, uses_limited_api=False):
        """ Initialise the sources. """

        super().__init__(project, name, fq_name.split('.')[-1],
                uses_limited_api=uses_limited_api)

        self.fq_name = fq_name

        self.exceptions = False
        self.static = False

    def get_install_subdir(self):
        """ Return the sub-directory the extension module should be installed
        in relative to the eventual target directory.
        """

        return os.sep.join(self.fq_name.split('.')[:-1])

    def get_module_extension(self):
        """ Return the filename extension that a module should have. """

        if self.project.py_platform == 'win32':
            return '.pyd'

        suffixes = importlib.machinery.EXTENSION_SUFFIXES

        if self.uses_limited_api:
            for s in suffixes:
                if '.abi3' in s:
                    return s

        return suffixes[0]


class BuildableBindings(BuildableModule):
    """ Encapsulate the sources used to build the extension module for a set of
    bindings.
    """

    def __init__(self, bindings, fq_name, *, uses_limited_api=False):
        """ Initialise the sources. """

        super().__init__(bindings.project, fq_name.split('.')[-1], fq_name,
                uses_limited_api=uses_limited_api)

        self.bindings = bindings

    def get_bindings_installable(self, name):
        """ Return an installable for the buildable's bindings directory. """

        target_subdir = os.path.join(self.project.get_bindings_dir(),
                self.name)

        return Installable(name, target_subdir=target_subdir)

    def write_configuration(self, bindings_dir):
        """ Write the configuration of the bindings and add it as an
        installable.
        """

        # Make sure the bindings directory exists.
        bindings_dir = os.path.join(bindings_dir, self.name)
        os.makedirs(bindings_dir, exist_ok=True)

        config_path = os.path.join(bindings_dir, self.name + '.toml')

        # Create an installable for the configuration file.
        installable = self.get_bindings_installable('config')
        installable.files.append(config_path)
        self.installables.append(installable)

        # Write the configuration file.
        bindings = self.bindings

        with open(config_path, 'w') as cf:
            sip_version_str = SIP_VERSION_STR if self.project.version_info else ''
            tags = ', '.join(['"{}"'.format(t) for t in bindings.tags])
            disabled = ', '.join(
                    ['"{}"'.format(f) for f in bindings.disabled_features])

            cf.write("# Automatically generated configuration for {0}.\n".format(self.fq_name))
            cf.write('''
sip-version = "{}"
sip-abi-version = "{}"
module-tags = [{}]
module-disabled-features = [{}]
'''.format(sip_version_str, self.project.abi_version, tags, disabled))
