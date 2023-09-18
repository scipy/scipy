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


import collections
import os
import packaging
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import warnings

from .abstract_builder import AbstractBuilder
from .abstract_project import AbstractProject
from .bindings import Bindings
from .configurable import Configurable, Option
from .exceptions import UserException
from .module import resolve_abi_version
from .py_versions import OLDEST_SUPPORTED_MINOR
from .pyproject import PyProjectException, PyProjectOptionException


class Project(AbstractProject, Configurable):
    """ Encapsulate a project containing one or more sets of bindings. """

    # The configurable options.
    _options = (
        # The minimum required ABI version of the sip module.
        Option('abi_version'),

        # The callable that will return an Bindings instance.  This is used for
        # bindings implicitly defined in the .toml file.
        Option('bindings_factory'),

        # The callable that will return an AbstractBuilder instance.
        Option('builder_factory'),

        # The list of console script entry points.
        Option('console_scripts', option_type=list),

        # Set if an __init__.py should be installed.
        Option('dunder_init', option_type=bool, default=False),

        # The list of GUI script entry points.
        Option('gui_scripts', option_type=list),

        # The minimum macOS version required by the project.  This is used to
        # determine the correct platform tag to use for macOS wheels.
        Option('minimum_macos_version'),

        # Set if building for a debug version of Python.
        Option('py_debug', option_type=bool),

        # The name of the directory containing Python.h.
        Option('py_include_dir', default=sysconfig.get_path('include')),

        # The name of the target Python platform.
        Option('py_platform'),

        # The major version number of the target Python installation.
        Option('py_major_version', option_type=int),

        # The minor version number of the target Python installation.
        Option('py_minor_version', option_type=int),

        # The name of the directory containing the .sip files.  If the sip
        # module is shared then each set of bindings is in its own
        # sub-directory.
        Option('sip_files_dir', default='.'),

        # The list of files and directories, specified as glob patterns
        # relative to the project directory, that should be excluded from an
        # sdist.
        Option('sdist_excludes', option_type=list),

        # The list of additional directories to search for .sip files.
        Option('sip_include_dirs', option_type=list),

        # The fully qualified name of the sip module.
        Option('sip_module'),

        # The list of files and directories, specified as glob patterns, that
        # should be included in a wheel.  If a pattern is relative then it is
        # taken as being relative to the project directory.  If an element of
        # the list is a string then it is a pattern and files and directories
        # are installed in the target directory.  If an element is a 2-tuple
        # then the first part is the pattern and the second part is the name of
        # a sub-directory relative to the target directory where the files and
        # directories are installed.
        Option('wheel_includes', option_type=list),

        # The user-configurable options.
        Option('quiet', option_type=bool,
                help="disable all progress messages"),
        Option('verbose', option_type=bool,
                help="enable verbose progress messages"),
        Option('name', help="the name used in sdist and wheel file names",
                metavar="NAME", tools=['sdist', 'wheel']),
        Option('build_dir', help="the build directory", metavar="DIR"),
        Option('build_tag', help="the build tag to be used in the wheel name",
                metavar="TAG", tools=['wheel']),
        Option('distinfo', option_type=bool, inverted=True,
                help="don't create a .dist-info directory", tools=['install']),
        Option('manylinux', option_type=bool, inverted=True,
                help="disable the use of manylinux in the platform tag used "
                        "in the wheel name",
                tools=['wheel']),
        Option('minimum_glibc_version',
                help="the minimum GLIBC version to be used in the platform "
                        "tag of Linux wheels",
                metavar="M.N", tools=['wheel']),
        Option('scripts_dir', default=os.path.dirname(sys.executable),
                help="the scripts installation directory", metavar="DIR",
                tools=['build', 'install']),
        Option('target_dir', default=sysconfig.get_path('platlib'),
                help="the target installation directory", metavar="DIR",
                tools=['build', 'install']),
        Option('api_dir', help="generate a QScintilla .api file in DIR",
                metavar="DIR"),
        Option('compile', option_type=bool, inverted=True,
                help="disable the compilation of the generated code",
                tools=['build']),
        Option('version_info', option_type=bool, inverted=True,
                help="disable any reference to the SIP version number in "
                        "generated code",
                tools=['build']),
    )

    # The configurable options for multiple bindings.
    _multibindings_options = (
        Option('disable', option_type=list, help="disable the NAME bindings",
                metavar="NAME"),
        Option('enable', option_type=list, help="enable the NAME bindings",
                metavar="NAME"),
    )

    def __init__(self, **kwargs):
        """ Initialise the project. """

        super().__init__()

        # The current directory should contain the .toml file.
        self.root_dir = os.getcwd()
        self.arguments = None
        self.bindings = collections.OrderedDict()
        self.bindings_factories = []
        self.builder = None
        self.buildables = []
        self.installables = []

        self._metadata_overrides = None
        self._temp_build_dir = None

        self.initialise_options(kwargs)

    def apply_nonuser_defaults(self, tool):
        """ Set default values for non-user options that haven't been set yet.
        """

        if self.bindings_factory is None:
            self.bindings_factory = Bindings
        elif isinstance(self.bindings_factory, str):
            # Convert the name to a callable.
            self.bindings_factory = self.import_callable(self.bindings_factory,
                    Bindings)

        if self.py_major_version is None or self.py_minor_version is None:
            self.py_major_version = sys.hexversion >> 24
            self.py_minor_version = (sys.hexversion >> 16) & 0x0ff

        if self.builder_factory is None:
            if (self.py_major_version, self.py_minor_version) >= (3, 10):
                from .setuptools_builder import SetuptoolsBuilder
                self.builder_factory = SetuptoolsBuilder
            else:
                from .distutils_builder import DistutilsBuilder
                self.builder_factory = DistutilsBuilder
        elif isinstance(self.builder_factory, str):
            # Convert the name to a callable.
            self.builder_factory = self.import_callable(self.builder_factory,
                    AbstractBuilder)

        if self.py_platform is None:
            self.py_platform = sys.platform

        if self.py_debug is None:
            self.py_debug = hasattr(sys, 'gettotalrefcount')

        super().apply_nonuser_defaults(tool)

    def apply_user_defaults(self, tool):
        """ Set default values for user options that haven't been set yet. """

        # This is only used when creating sdist and wheel files.
        if self.name is None:
            self.name = self.metadata['name']

        # For the build tool we want build_dir to default to a local 'build'
        # directory (which we won't remove).  However, for other tools (and for
        # PEP 517 frontends) we want to use a temporary directory in case the
        # current directory is read-only.
        if self.build_dir is None:
            if tool == 'build':
                self.build_dir = 'build'
            else:
                self._temp_build_dir = tempfile.TemporaryDirectory()
                self.build_dir = self._temp_build_dir.name

        super().apply_user_defaults(tool)

        # Adjust the list of bindings according to what has been explicitly
        # enabled and disabled.
        self._enable_disable_bindings()

        # Set the user defaults for the builder and bindings.
        self.builder.apply_user_defaults(tool)

        for bindings in self.bindings.values():
            bindings.apply_user_defaults(tool)

    def build(self):
        """ Build the project in-situ. """

        self.builder.build()

    def build_sdist(self, sdist_directory):
        """ Build an sdist for the project and return the name of the sdist
        file.
        """

        sdist_file = self.builder.build_sdist(sdist_directory)
        self._remove_build_dir()

        return sdist_file

    def build_wheel(self, wheel_directory):
        """ Build a wheel for the project and return the name of the wheel
        file.
        """

        wheel_file = self.builder.build_wheel(wheel_directory)
        self._remove_build_dir()

        return wheel_file

    def get_bindings_dir(self):
        """ Return the name of the 'bindings' directory relative to the
        eventual target directory.
        """

        return os.path.join(self.get_package_dir(), 'bindings')

    def get_distinfo_dir(self, target_dir):
        """ Return the name of the .dist-info directory for a target directory.
        """

        return os.path.join(target_dir,
                '{}-{}.dist-info'.format(self.name.replace('-', '_'),
                self.version_str))

    def get_dunder_init(self):
        """ Return the contents of the __init__.py to install. """

        # This default implementation will create an empty file.
        return ''

    def get_metadata_overrides(self):
        """ Return a mapping of PEP 566 metadata names and values that will
        override any corresponding values defined in the pyproject.toml file.
        A typical use is to determine a project's version dynamically.
        """

        # This default implementation does not override any metadata.
        return {}

    def get_options(self):
        """ Return the list of configurable options. """

        options = super().get_options()
        options.extend(self._options)
        options.extend(self._multibindings_options)

        return options

    def get_package_dir(self):
        """ Return the name of the package directory relative to the eventual
        target directory.  This is the directory containing the shared sip
        module (if there is one) or the target directory (if not).  It will
        normally be where the individual bindings are installed.
        """

        if self.sip_module:
            name_parts = self.sip_module.split('.')
            del name_parts[-1]

            return os.path.join(*name_parts)

        return ''

    def get_platform_tag(self):
        """ Return the platform tag to use in a wheel name.  This default
        implementation uses the platform name and applies PEP defined
        conventions depending on OS version and GLIBC version as appropriate.
        """

        platform_tag = sysconfig.get_platform()

        if self.py_platform == 'darwin' and self.minimum_macos_version:
            # We expect a three part tag so leave anything else unchanged.
            parts = platform_tag.split('-')
            if len(parts) == 3:
                min_major = int(self.minimum_macos_version[0])
                min_minor = int(self.minimum_macos_version[1])

                # For arm64 binaries enforce a valid minimum macOS version.
                if parts[2] == 'arm64' and min_major < 11:
                    min_major = 11
                    min_minor = 0

                parts[1] = '{}.{}'.format(min_major, min_minor)

                platform_tag = '-'.join(parts)

        elif self.py_platform == 'linux' and self.manylinux:
            # We expect a two part tag so leave anything else unchanged.
            parts = platform_tag.split('-')
            if len(parts) == 2:
                if self.minimum_glibc_version:
                    major, minor = self.minimum_glibc_version
                else:
                    major, minor = 2, 5

                parts[0] = 'manylinux'
                parts.insert(1, '{}.{}'.format(major, minor))

                platform_tag = '-'.join(parts)

        return platform_tag.replace('.', '_').replace('-', '_')

    def get_requires_dists(self):
        """ Return any 'Requires-Dist' to add to the project's meta-data. """

        # The only requirement is for the sip module.
        if not self.sip_module:
            return []

        requires_dist = self.metadata.get('requires-dist')
        if requires_dist is None:
            requires_dist = []
        elif isinstance(requires_dist, str):
            requires_dist = [requires_dist]

        # Ignore if the module is already defined.
        sip_project_name = self.sip_module.replace('.', '-')

        for rd in requires_dist:
            if rd.split()[0] == sip_project_name:
                return []

        next_abi_major = int(self.abi_version.split('.')[0]) + 1

        return ['{} (>={}, <{})'.format(sip_project_name, self.abi_version,
                next_abi_major)]

    def get_sip_distinfo_command_line(self, sip_distinfo, inventory,
            generator=None, wheel_tag=None, generator_version=None):
        """ Return a sequence of command line arguments to invoke sip-distinfo.
        """

        args = [
            sip_distinfo,

            '--inventory',
            inventory,

            '--project-root',
            self.root_dir,

            '--prefix',
            '\\"$(INSTALL_ROOT)\\"',
        ]

        if generator is not None:
            args.append('--generator')
            args.append(generator)

        if generator_version is not None:
            args.append('--generator-version')
            args.append(generator_version)

        if wheel_tag is not None:
            args.append('--wheel-tag')
            args.append(wheel_tag)

        for ep in self.console_scripts:
            args.append('--console-script')
            args.append(ep.replace(' ', ''))

        for ep in self.gui_scripts:
            args.append('--gui-script')
            args.append(ep.replace(' ', ''))

        for rd in self.get_requires_dists():
            args.append('--requires-dist')
            args.append('\\"{}\\"'.format(rd))

        for metadata, value in self._metadata_overrides.items():
            if value:
                metadata += '=' + value

            if ' ' in metadata:
                metadata = '\\"' + metadata + '\\"'

            args.append('--metadata')
            args.append(metadata)

        return args

    def install(self):
        """ Install the project. """

        self.builder.install()
        self._remove_build_dir()

    @property
    def minimum_glibc_version(self):
        """ The getter for the minimum GLIBC version. """

        return self._minimum_glibc_version

    @minimum_glibc_version.setter
    def minimum_glibc_version(self, value):
        """ The setter for the minimum GLIBC version. """

        # Handle the initial creation of the option value.
        if value is None and not hasattr(self, '_minimum_glibc_version'):
            self._minimum_glibc_version = None
            return

        # Make sure any minimum GLIBC version is valid and convert it to a
        # 2-tuple.
        if value:
            try:
                value = self._convert_major_minor(value)
            except ValueError:
                raise PyProjectOptionException('minimum-glibc-version',
                        "'{0}' is an invalid GLIBC version number".format(
                                value),
                        section_name='tool.sip.project')

        self._minimum_glibc_version = value

    @property
    def minimum_macos_version(self):
        """ The getter for the minimum macOS version. """

        return self._minimum_macos_version

    @minimum_macos_version.setter
    def minimum_macos_version(self, value):
        """ The setter for the minimum macOS version. """

        # Handle the initial creation of the option value.
        if value is None and not hasattr(self, '_minimum_macos_version'):
            self._minimum_macos_version = None
            return

        # Make sure any minimum macOS version is valid and convert it to a
        # 2-tuple.
        if value:
            try:
                value = self._convert_major_minor(value)
            except ValueError:
                raise PyProjectOptionException('minimum-macos-version',
                        "'{0}' is an invalid macOS version number".format(
                                value),
                        section_name='tool.sip.project')

        self._minimum_macos_version = value

    @staticmethod
    def open_for_writing(fname):
        """ Open a file for writing while handling any errors. """

        try:
            return open(fname, 'w')
        except IOError as e:
            raise UserException(
                    "There was an error creating '{0}' - make sure you have "
                    " write permission on the parent directory".format(fname),
                    detail=str(e))

    def progress(self, message):
        """ Print a progress message unless they are disabled. """

        if not self.quiet:
            if message[-1] != '.':
                message += '...'

            print(message, flush=True)

    def project_path(self, path, relative_to=None):
        """ Return a normalised version of a path.  A relative path is assumed
        to be relate to the project directory or some other provided directory.
        """

        path = os.path.normpath(path)

        if os.path.isabs(path):
            return path

        if relative_to is None:
            relative_to = self.root_dir

        return os.path.normpath(os.path.join(relative_to, path))

    def read_command_pipe(self, args, *, and_stderr=False, fatal=True):
        """ A generator for each line of a pipe from a command's stdout. """

        cmd = ' '.join(args)

        if self.verbose:
            print(cmd, flush=True)

        stderr = subprocess.STDOUT if and_stderr else subprocess.PIPE

        with subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                stdout=subprocess.PIPE, stderr=stderr) as pipe:
            for line in pipe.stdout:
                yield str(line, encoding=sys.stdout.encoding)

        if pipe.returncode != 0 and fatal:
            raise UserException(
                    "'{0}' failed returning {1}".format(cmd, pipe.returncode))

    def run_command(self, args, *, fatal=True):
        """ Run a command and display the output if requested. """

        # Read stdout and stderr until there is no more output.
        for line in self.read_command_pipe(args, and_stderr=True, fatal=fatal):
            if self.verbose:
                sys.stdout.write(line)

    def setup(self, pyproject, tool, tool_description):
        """ Complete the configuration of the project. """

        # Create any programmatically defined bindings.
        for bindings_factory in self.bindings_factories:
            bindings = bindings_factory(self)
            self.bindings[bindings.name] = bindings

        # Set the initial configuration from the pyproject.toml file.
        self._set_initial_configuration(pyproject, tool)

        # Add any tool-specific command line arguments for (so far unspecified)
        # parts of the configuration.
        self._configure_from_arguments(tool, tool_description)

        # Now that any help has been given we can report a problematic
        # pyproject.toml file.
        if pyproject.toml_error:
            raise PyProjectException(pyproject.toml_error)

        # Make sure the configuration is complete.
        self.apply_user_defaults(tool)

        # Configure the warnings module.
        if not self.verbose:
            warnings.simplefilter('ignore', UserWarning)

        # Make sure we have a clean build directory and make it current.
        if self._temp_build_dir is None:
            self.build_dir = os.path.abspath(self.build_dir)
            shutil.rmtree(self.build_dir, ignore_errors=True)
            os.mkdir(self.build_dir)

        os.chdir(self.build_dir)

        # Allow a sub-class (in a user supplied script) to make any updates to
        # the configuration.
        self.update(tool)

        os.chdir(self.root_dir)

        # Make sure the configuration is correct after any user supplied script
        # has messed with it.
        self.verify_configuration(tool)

        if tool in Option.BUILD_TOOLS and self.bindings:
            self.progress(
                    "These bindings will be built: {}.".format(
                            ', '.join(self.bindings.keys())))

    def update(self, tool):
        """ This should be re-implemented by any user supplied sub-class to
        carry out any updates to the configuration as required.  The current
        directory will be the temporary build directory.
        """

        # This default implementation calls update_buildable_bindings().
        if tool in Option.BUILD_TOOLS:
            self.update_buildable_bindings()

    def update_buildable_bindings(self):
        """ Update the list of bindings to ensure they are either buildable or
        have been explicitly enabled.
        """

        # Explicitly enabled bindings are assumed to be buildable.
        if self.enable:
            return

        for b in list(self.bindings.values()):
            if not b.is_buildable():
                del self.bindings[b.name]

        if len(self.bindings) == 0:
            raise UserException("There are no bindings that can be built")

    def verify_configuration(self, tool):
        """ Verify that the configuration is complete and consistent. """

        # Make sure any build tag is valid.
        if self.build_tag and not self.build_tag[0].isdigit():
            raise PyProjectOptionException('build-tag',
                    "'{0}' must begin with a digit".format(self.build_tag),
                    section_name='tool.sip.project')

        # Make sure relevent paths are absolute and use native separators.
        self.sip_files_dir = self.project_path(self.sip_files_dir)
        self.sip_include_dirs = [self.project_path(d)
                for d in self.sip_include_dirs]

        # Check that the targeted version of Python isn't too old.  We hope
        # that we will support newer versions automatically, but it's not
        # guaranteed.
        py_version = (self.py_major_version, self.py_minor_version)
        first_version = (3, OLDEST_SUPPORTED_MINOR)

        if py_version < first_version or self.py_major_version > 3:
            raise UserException(
                    "Python v{}.{} is not supported".format(
                            self.py_major_version, self.py_minor_version))

        # Get the supported ABI version.  (The actual version may have a later
        # minor version.)
        self.abi_version = resolve_abi_version(self.abi_version, module=False)

        # Checks for standalone projects.
        if tool in Option.BUILD_TOOLS and not self.sip_module:
            # Check there is only one set of bindings.
            if len(self.bindings) > 1:
                raise PyProjectOptionException('sip-module',
                        "must be defined when the project contains multiple "
                        "sets of bindings")

            # Make sure __init__.py is disabled.
            self.dunder_init = False

        # Check any wheel includes are valid and make sure all elements are
        # 2-tuples.
        normalised = []

        for wheel_include in self.wheel_includes:
            if isinstance(wheel_include, str):
                normalised.append((wheel_include, None))
            else:
                try:
                    wheel_include, subdir = wheel_include
                except TypeError:
                    wheel_include = subdir = None

                if isinstance(wheel_include, str) and isinstance(subdir, str):
                    normalised.append((wheel_include, subdir))
                else:
                    raise PyProjectOptionException('wheel-includes',
                            "elements must be strings or 2-tuples of strings")

        self.wheel_includes = normalised

        # Make sure that any .api directory is relative when building a wheel.
        if tool == 'wheel' and self.api_dir:
            if os.path.isabs(self.api_dir) or os.path.dirname(self.api_dir) == '..':
                raise PyProjectOptionException('api-dir',
                        "must be relative when building a wheel")

        # Verify the configuration of the builder and bindings.
        self.builder.verify_configuration(tool)

        for bindings in self.bindings.values():
            bindings.verify_configuration(tool)

    def _configure_from_arguments(self, tool, tool_description):
        """ Update the configuration from any user supplied arguments. """

        from argparse import SUPPRESS
        from .argument_parser import ArgumentParser

        parser = ArgumentParser(tool_description, argument_default=SUPPRESS)

        # Add the user configurable options to the parser.
        all_options = {}
        
        options = self.get_options()
        if len(self.bindings) < 2:
            # Remove the options that only make sense where the project has
            # multiple bindings.
            for multi in self._multibindings_options:
                options.remove(multi)

        self.add_command_line_options(parser, tool, all_options,
                options=options)

        self.builder.add_command_line_options(parser, tool, all_options)

        for bindings in self.bindings.values():
            bindings.add_command_line_options(parser, tool, all_options)

        # Parse the arguments and update the corresponding configurables.
        args = parser.parse_args(self.arguments)

        for option, configurables in all_options.items():
            for configurable in configurables:
                if hasattr(args, option.dest):
                    setattr(configurable, option.name,
                            getattr(args, option.dest))

    @staticmethod
    def _convert_major_minor(value):
        """ Convert a 'major.minor' version number to a 2-tuple of integers.
        Raise a ValueError exception if it is invalid.
        """

        parts = value.split('.')
        if len(parts) != 2:
            raise ValueError()

        return int(parts[0]), int(parts[1])

    def _enable_disable_bindings(self):
        """ Check the enabled bindings are valid and remove any disabled ones.
        """

        names = list(self.bindings.keys())

        # Check that any explicitly enabled bindings are valid.
        if self.enable:
            for enabled in self.enable:
                if enabled not in names:
                    raise UserException(
                            "Unknown enabled bindings '{0}'".format(enabled))

            # Only include explicitly enabled bindings.
            for b in list(self.bindings.values()):
                if b.name not in self.enable:
                    del self.bindings[b.name]

        # Check that any explicitly disabled bindings are valid.
        if self.disable:
            for disabled in self.disable:
                if disabled not in names:
                    raise UserException(
                            "Unknown disabled bindings '{0}'".format(disabled))

            # Remove any explicitly disabled bindings.
            for b in list(self.bindings.values()):
                if b.name in self.disable:
                    del self.bindings[b.name]

    def _remove_build_dir(self):
        """ Remove the build directory. """

        self._temp_build_dir = None

    def _set_initial_configuration(self, pyproject, tool):
        """ Set the project's initial configuration. """

        # Get the metadata and extract the version.
        self.metadata = pyproject.get_metadata()
        self._metadata_overrides = self.get_metadata_overrides()
        self.metadata.update(self._metadata_overrides)
        self.version_str = self.metadata['version']

        # Convert the version as a string to number.
        base_version = packaging.version.parse(self.version_str).base_version
        base_version = base_version.split('.')

        while len(base_version) < 3:
            base_version.append('0')

        version = 0
        for part in base_version:
            version <<= 8

            try:
                version += int(part)
            except ValueError:
                raise PyProjectOptionException('version',
                        "'{0}' is an invalid version number".format(
                                self.version_str),
                        section_name='tool.sip.metadata')

        self.version = version

        # Configure the project.
        self.configure(pyproject, 'tool.sip.project', tool)

        # Create and configure the builder.
        self.builder = self.builder_factory(self)
        self.builder.configure(pyproject, 'tool.sip.builder', tool)

        # For each set of bindings configuration make sure a bindings object
        # exists, creating it if necessary.
        bindings_sections = pyproject.get_section('tool.sip.bindings')
        if bindings_sections is not None:
            for name in bindings_sections.keys():
                if name not in self.bindings:
                    bindings = self.bindings_factory(self, name)
                    self.bindings[bindings.name] = bindings

        # Add a default set of bindings if none were defined.
        if not self.bindings:
            bindings = self.bindings_factory(self, self.metadata['name'])
            self.bindings[bindings.name] = bindings

        # Now configure each set of bindings.
        for bindings in self.bindings.values():
            bindings.configure(pyproject, 'tool.sip.bindings.' + bindings.name,
                    tool)
