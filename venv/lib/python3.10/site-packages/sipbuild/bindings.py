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
import sys

from .buildable import BuildableBindings
from .code_generator import generateCode, py2c
from .configurable import Configurable, Option
from .exceptions import UserException
from .generator import parse, resolve
from .generator.outputs import output_api, output_extract, output_pyi
from .installable import Installable
from .module import copy_nonshared_sources
from .version import SIP_VERSION


class Bindings(Configurable):
    """ The encapsulation of a module's bindings. """

    # The configurable options.
    _options = (
        # Any bindings level builder-specific settings.
        Option('builder_settings', option_type=list),

        # The list of #define names and values in the format "NAME" or
        # "NAME=VALUE".
        Option('define_macros', option_type=list),

        # Set if exception support is enabled.
        Option('exceptions', option_type=bool),

        # The list of extra compiler arguments.
        Option('extra_compile_args', option_type=list),

        # The list of extra linker arguments.
        Option('extra_link_args', option_type=list),

        # The list of extra compiled object files to link.
        Option('extra_objects', option_type=list),

        # The list of extracts to generate.
        Option('generate_extracts', option_type=list),

        # The list of additional .h files.
        Option('headers', option_type=list),

        # The list of additional C/C++ include directories to search.
        Option('include_dirs', option_type=list),

        # Set if the bindings are internal.  Internal bindings don't have their
        # .sip, .pyi or .api files installed.
        Option('internal', option_type=bool),

        # The list of library names to link against.
        Option('libraries', option_type=list),

        # The list of C/C++ library directories to search.
        Option('library_dirs', option_type=list),

        # Set to always release the Python GIL.
        Option('release_gil', option_type=bool),

        # Set to compile the module as a static library.
        Option('static', option_type=bool),

        # The name of the .sip file that specifies the bindings.  If it is
        # relative then it is relative to the project's 'sip_files_dir'.
        Option('sip_file'),

        # The filename extension to use for generated source files.
        Option('source_suffix'),

        # The list of additional C/C++ source files to compile and link.
        Option('sources', option_type=list),

        # The list of tags to enable.
        Option('tags', option_type=list),

        # The user-configurable options.  Although the use of a corresponding
        # command line option will affect all sets of bindings, putting them
        # here (as opposed to in Project) means they can have individual
        # values specified in pyproject.toml.
        Option('concatenate', option_type=int,
                help="concatenate the generated bindings into N source files",
                metavar="N"),
        Option('debug', option_type=bool, help="build with debugging symbols"),
        Option('disabled_features', option_type=list,
                help="disable the TAG feature tag", metavar="TAG"),
        Option('docstrings', option_type=bool, inverted=True,
                help="disable the generation of docstrings"),
        Option('pep484_pyi', option_type=bool,
                help="enable the generation of PEP 484 .pyi files"),
        Option('protected_is_public', option_type=bool,
                help="enable the protected/public hack (default on non-Windows)"),
        Option('protected_is_public', option_type=bool, inverted=True,
                help="disable the protected/public hack (default on Windows)"),
        Option('tracing', option_type=bool, help="build with tracing support"),
    )

    def __init__(self, project, name, **kwargs):
        """ Initialise the bindings. """

        super().__init__()

        self.project = project
        self.name = name

        self.initialise_options(kwargs)

    def apply_nonuser_defaults(self, tool):
        """ Set default values for each non-user configurable option that
        hasn't been set yet.
        """

        # Provide a default .sip file name if needed.
        if self.sip_file is None:
            self.sip_file = self.name + '.sip'

        super().apply_nonuser_defaults(tool)

    def apply_user_defaults(self, tool):
        """ Set default values for user options that haven't been set yet. """

        if self.protected_is_public is None:
            self.protected_is_public = (self.project.py_platform != 'win32')

        super().apply_user_defaults(tool)

    def generate(self):
        """ Generate the bindings source code and optional additional extracts.
        and return a BuildableBindings instance containing the details of
        everything needed to build the bindings.
        """

        project = self.project

        # The old parser had no concept of the encoding of a .sip file.  For
        # the moment we say that files should be UTF-8.  If that proves to be a
        # problem then a project-specific encoding should be able to be
        # specified in pyproject.toml which would apply to all .sip files that
        # make up the project.
        encoding = 'UTF-8'

        # Parse the input file.
        spec, modules, sip_files = parse(self.sip_file, SIP_VERSION, encoding,
                project.abi_version, self.tags, self.disabled_features,
                self.protected_is_public, self._sip_include_dirs,
                project.sip_module)

        # Resolve the types.
        resolve(spec, modules)

        pt = py2c(spec, encoding)

        module = spec.module

        uses_limited_api = module.use_limited_api or spec.is_composite

        # The details of things that will have been generated.  Note that we
        # don't include anything for .api files or generic extracts as the
        # arguments include a file name.
        buildable = BuildableBindings(self, module.fq_py_name.name,
                uses_limited_api=uses_limited_api)

        buildable.builder_settings.extend(self.builder_settings)
        buildable.debug = self.debug
        buildable.exceptions = self.exceptions
        buildable.extra_compile_args = self.extra_compile_args
        buildable.extra_link_args = self.extra_link_args
        buildable.extra_objects = self.extra_objects
        buildable.static = self.static

        # Generate any API file.
        if project.api_dir and not self.internal:
            project.progress(
                    "Generating the {0} .api file".format(buildable.target))

            output_api(spec,
                    os.path.join(project.build_dir, buildable.target + '.api'))

        # Generate any extracts.
        for extract_ref in self.generate_extracts:
            output_extract(spec, extract_ref)

        # Generate any type hints file.
        if self.pep484_pyi and not self.internal:
            project.progress(
                    "Generating the {0} .pyi file".format(buildable.target))

            pyi_path = os.path.join(buildable.build_dir,
                    buildable.target + '.pyi')

            output_pyi(spec, pyi_path)

            installable = Installable('pyi',
                    target_subdir=buildable.get_install_subdir())
            installable.files.append(pyi_path)
            buildable.installables.append(installable)

        # Generate the bindings.
        header, sources = generateCode(pt, buildable.build_dir,
                self.source_suffix, self.exceptions, self.tracing,
                self.release_gil, self.concatenate, self.tags,
                self.disabled_features, self.docstrings, project.py_debug)

        if header:
            buildable.headers.append(header)

        buildable.headers.extend(self.headers)

        buildable.sources.extend(sources)

        # Add the sip module code if it is not shared.
        buildable.include_dirs.append(buildable.build_dir)

        if project.sip_module:
            # sip.h will already be in the build directory.
            buildable.include_dirs.append(project.build_dir)

            if not self.internal:
                # Add an installable for the .sip files.
                installable = buildable.get_bindings_installable('sip')

                sip_dir = os.path.dirname(self.sip_file)

                for fn in sip_files:
                    sip_path = os.path.join(sip_dir, fn)

                    # The code generator does not report the full pathname of a
                    # .sip file (only names relative to the search directory in
                    # which it was found).  Therefore we need to check if it is
                    # actually in the directory we are installing from and
                    # ignore it if not.  This isn't really the right thing to
                    # do but is actually what we want when we have optional
                    # license .sip files.
                    if os.path.isfile(sip_path):
                        installable.files.append(sip_path)

                buildable.installables.append(installable)
        else:
            buildable.sources.extend(
                    copy_nonshared_sources(project.abi_version.split('.')[0],
                            buildable.build_dir))

        buildable.include_dirs.extend(self.include_dirs)
        buildable.sources.extend(self.sources)

        if self.protected_is_public:
            buildable.define_macros.append('SIP_PROTECTED_IS_PUBLIC')
            buildable.define_macros.append('protected=public')

        buildable.define_macros.extend(self.define_macros)

        buildable.libraries.extend(self.libraries)
        buildable.library_dirs.extend(self.library_dirs)

        return buildable

    def get_options(self):
        """ Return the list of configurable options. """

        options = super().get_options()
        options.extend(self._options)

        return options

    def is_buildable(self):
        """ Return True if the bindings are buildable.  This will not be called
        if the bindings have been explicitly enabled.
        """

        return True

    def verify_configuration(self, tool):
        """ Verify that the configuration is complete and consistent. """

        if tool not in Option.BUILD_TOOLS:
            return

        project = self.project

        super().verify_configuration(tool)

        # Make sure relevent paths are absolute and use native separators.
        self.extra_objects = [project.project_path(o)
                for o in self.extra_objects]
        self.headers = [project.project_path(h) for h in self.headers]
        self.include_dirs = [project.project_path(d)
                for d in self.include_dirs]
        self.library_dirs = [project.project_path(d)
                for d in self.library_dirs]
        self.sip_file = project.project_path(self.sip_file,
                relative_to=project.sip_files_dir)
        self.sources = [project.project_path(s) for s in self.sources]

        # Check the .sip file exists.
        if not os.path.isfile(self.sip_file):
            raise UserException(
                    "the file '{0}' for the {1} bindings does not "
                            "exist".format(self.sip_file, self.name))

        # On Windows the interpreter must be a debug build if a debug version
        # is to be built and vice versa.
        if sys.platform == 'win32':
            if self.debug:
                if not project.py_debug:
                    raise UserException(
                            "A debug version of Python must be used when "
                            "building a debug version of the {0} "
                            "bindings".format(self.name))
            elif project.py_debug:
                raise UserException(
                        "A debug version of the {0} bindings must be built "
                        "when a debug version of Python is used".format(
                                self.name))

        if not self.source_suffix:
            self.source_suffix = None
