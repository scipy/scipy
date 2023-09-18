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


from abc import abstractmethod
import glob
import os
import shutil
import stat
import sys

from .abstract_builder import AbstractBuilder
from .buildable import BuildableFromSources
from .code_generator import set_globals
from .distinfo import write_metadata
from .exceptions import UserException
from .installable import Installable
from .module import copy_sip_h, copy_sip_pyi
from .py_versions import OLDEST_SUPPORTED_MINOR
from .version import SIP_VERSION, SIP_VERSION_STR


class Builder(AbstractBuilder):
    """ The default base implementation of a project builder. """

    def build(self):
        """ Build the project in-situ. """

        self._generate_bindings()
        self._generate_scripts()

        if self.project.compile:
            self.build_project(self.project.target_dir)

    @abstractmethod
    def build_executable(self, buildable, *, fatal=True):
        """ Build an executable from a BuildableExecutable object and return
        the relative pathname of the executable.
        """

    @abstractmethod
    def build_project(self, target_dir, *, wheel_tag=None):
        """ Build the project. """

    def build_sdist(self, sdist_directory):
        """ Build an sdist for the project and return the name of the sdist
        file.
        """

        project = self.project

        # The sdist name.
        sdist_name = '{}-{}'.format(project.name.replace('-', '_'),
                project.version_str)

        # Create the sdist root directory.
        sdist_root = os.path.join(project.build_dir, sdist_name)
        os.mkdir(sdist_root)

        # Get a list of all excluded files.
        excluded = []
        for patt in project.sdist_excludes:
            excluded.extend(glob.glob(os.path.join(project.root_dir, patt)))

        for dname, dirnames, filenames in os.walk(project.root_dir):
            # Always ignore certain directories.
            if dname == project.build_dir:
                del dirnames[:]
                continue

            try:
                dirnames.remove('__pycache__')
            except ValueError:
                pass

            # Copy each file non-excluded file.
            for s_fn in filenames:
                s_fn_path = os.path.join(dname, s_fn)

                if s_fn_path in excluded:
                    continue

                d_fn_path = os.path.join(sdist_root,
                        os.path.relpath(s_fn_path, project.root_dir))
                os.makedirs(os.path.dirname(d_fn_path), exist_ok=True)

                shutil.copy2(s_fn_path, d_fn_path)

        # Create the PKG-INFO file.  This is assumed to be identical to the
        # .dist-info/METADATA file.
        write_metadata(project.metadata, project.get_requires_dists(),
                os.path.join(sdist_root, 'PKG-INFO'), project.root_dir)

        # Create the tarball.
        sdist_file = sdist_name + '.tar.gz'
        sdist_path = os.path.abspath(os.path.join(sdist_directory, sdist_file))

        saved_cwd = os.getcwd()
        os.chdir(project.build_dir)

        import tarfile

        tf = tarfile.open(sdist_path, 'w:gz', format=tarfile.PAX_FORMAT)
        tf.add(sdist_name)
        tf.close()

        os.chdir(saved_cwd)

        return sdist_file

    def build_wheel(self, wheel_directory):
        """ Build a wheel for the project and return the name of the wheel
        file.
        """

        project = self.project

        # Create a temporary directory for the wheel.
        wheel_build_dir = os.path.join(project.build_dir, 'wheel')
        os.mkdir(wheel_build_dir)

        # Build the wheel contents.
        self._generate_bindings()

        # If all buildables use the limited API then the wheel does.
        all_use_limited_api = True
        for buildable in project.buildables:
            if isinstance(buildable, BuildableFromSources):
                if not buildable.uses_limited_api:
                    all_use_limited_api = False
                    break

        # Create the wheel tag.
        wheel_tag = []

        if all_use_limited_api:
            # When the ABI tag is 'abi3' the interpreter tag is interpreted as
            # a minimum Python version.  This doesn't seem to be defined in a
            # PEP but is implemented in current pips.
            wheel_tag.append('cp3' + str(OLDEST_SUPPORTED_MINOR))
            wheel_tag.append('abi3')
        else:
            major_minor = '{}{}'.format((sys.hexversion >> 24) & 0xff,
                    (sys.hexversion >> 16) & 0xff)

            wheel_tag.append('cp{}'.format(major_minor))

            try:
                wheel_tag.append('cp' + major_minor + sys.abiflags)
            except AttributeError:
                wheel_tag.append('none')

        wheel_tag.append(project.get_platform_tag())

        wheel_tag = '-'.join(wheel_tag)

        # Add any wheel contents defined by the project.
        for nr, (patt, target_subdir) in enumerate(project.wheel_includes):
            if not os.path.isabs(patt):
                patt = os.path.join(project.root_dir, patt)

            wheel_includes = glob.glob(patt)
            if wheel_includes:
                installable = Installable(
                        'wheel_includes_{}'.format(nr) if nr else 'wheel_includes',
                        target_subdir=target_subdir)
                installable.files.extend(wheel_includes)
                project.installables.append(installable)

        # Build the project.
        self.build_project(wheel_build_dir, wheel_tag=wheel_tag)

        # Copy the wheel contents.
        self.install_project(wheel_build_dir, wheel_tag=wheel_tag)

        wheel_file = '{}-{}'.format(project.name.replace('-', '_'),
                project.version_str)

        if project.build_tag:
            wheel_file += '-{}'.format(project.build_tag)

        wheel_file += '-{}.whl'.format(wheel_tag)
        wheel_path = os.path.abspath(os.path.join(wheel_directory, wheel_file))

        # Create the .whl file.
        saved_cwd = os.getcwd()
        os.chdir(wheel_build_dir)

        from zipfile import ZipFile, ZIP_DEFLATED

        with ZipFile(wheel_path, 'w', compression=ZIP_DEFLATED) as zf:
            for dirpath, _, filenames in os.walk('.'):
                for filename in filenames:
                    # This will result in a name with no leading '.'.
                    name = os.path.relpath(os.path.join(dirpath, filename))

                    zf.write(name)

        os.chdir(saved_cwd)

        return wheel_file

    def install(self):
        """ Install the project. """

        target_dir = self.project.target_dir

        self._generate_bindings()
        self._generate_scripts()
        self.build_project(target_dir)
        self.install_project(target_dir)

    @abstractmethod
    def install_project(self, target_dir, *, wheel_tag=None):
        """ Install the project into a target directory. """

    def _generate_bindings(self):
        """ Generate the bindings for all enabled modules. """

        project = self.project

        abi_major_version, abi_minor_version = project.abi_version.split('.')

        # Get the list of directories to search for .sip files.
        sip_include_dirs = list(project.sip_include_dirs)

        if project.sip_module:
            # Add the project's sip directory.
            sip_include_dirs.append(project.sip_files_dir)

            # Add the local bindings directory to pick up the .toml files for
            # any other bindings in this package.
            local_bindings_dir = os.path.join(project.build_dir, 'bindings')
            sip_include_dirs.append(local_bindings_dir)

            # Add any bindings from previously installed packages.
            sip_include_dirs.append(
                    os.path.join(project.target_dir,
                            project.get_bindings_dir()))

            # Generate the sip.h file for the shared sip module.
            copy_sip_h(abi_major_version, project.build_dir,
                    project.sip_module, version_info=project.version_info)

        set_globals(SIP_VERSION,
                SIP_VERSION_STR if project.version_info else None,
                int(abi_major_version), int(abi_minor_version),
                project.sip_module, UserException)

        # Generate the code for each set of bindings.
        api_files = []

        for bindings in project.bindings.values():
            project.progress(
                    "Generating the {0} bindings".format(bindings.name))

            # Generate the source code.  We would prefer to pass the include
            # directories as an argument to generate but this is a fixed public
            # interface.
            bindings._sip_include_dirs = sip_include_dirs
            buildable = bindings.generate()
            del bindings._sip_include_dirs

            if not bindings.internal:
                api_files.append(
                        os.path.join(project.build_dir,
                                buildable.target + '.api'))

                # Generate the bindings configuration file.
                if project.sip_module:
                    buildable.write_configuration(local_bindings_dir)

            project.buildables.append(buildable)

        # Create __init__.py if required.
        if project.dunder_init:
            package_dir = project.get_package_dir()

            init_path = os.path.join(project.build_dir, '__init__.py')

            init_f = project.open_for_writing(init_path)
            init_f.write(project.get_dunder_init())
            init_f.close()

            installable = Installable('init', target_subdir=package_dir)
            installable.files.append(init_path)
            project.installables.append(installable)

            # Include sip.pyi if any of the bindings generate a .pyi file.
            for bindings in project.bindings.values():
                if bindings.pep484_pyi:
                    copy_sip_pyi(abi_major_version, project.build_dir)

                    installable = Installable('sip_pyi',
                            target_subdir=package_dir)
                    installable.files.append(
                            os.path.join(project.build_dir, 'sip.pyi'))
                    project.installables.append(installable)

                    # Create a PEP 561 marker file.
                    py_typed_path = os.path.join(project.build_dir, 'py.typed')
                    with open(py_typed_path, 'w') as _:
                        pass

                    installable = Installable('py_typed',
                            target_subdir=package_dir)
                    installable.files.append(py_typed_path)
                    project.installables.append(installable)

                    break

        # Create the .api file if required.
        if project.api_dir:
            api_fn = project.name + '.api'

            project.progress("Generating the {0} file".format(api_fn))

            # Concatanate the individual .api files.
            api_path = os.path.join(project.build_dir, api_fn)
            api_f = project.open_for_writing(api_path)

            for part_fn in api_files:
                with open(part_fn) as part_f:
                    api_f.write(part_f.read())

            api_f.close()

            # Add an Installable for the requested API file.
            installable = Installable('api', target_subdir=project.api_dir)
            installable.files.append(api_path)
            project.installables.append(installable)

    def _generate_scripts(self):
        """ Generate the scripts for any entry points. """

        project = self.project

        # Handle the trivial case.
        scripts = project.console_scripts + project.gui_scripts
        if not scripts:
            return

        # Create an installable for the scripts.
        installable = Installable('scripts', target_subdir=project.scripts_dir)

        for ep in scripts:
            # Parse the entry point.
            ep_parts = ep.replace(' ', '').split('=')

            if len(ep_parts) != 2:
                raise UserException(
                        "'{0}' is an invalid script specification".format(ep))

            script, module = ep_parts

            # Remove any callable name.
            module = module.split(':')[0]

            if project.py_platform == 'win32':
                script += '.bat'

            project.progress("Generating the {} script".format(script))

            script_path = os.path.join(project.build_dir, script)
            installable.files.append(script_path)

            script_f = project.open_for_writing(script_path)

            if project.py_platform == 'win32':
                script_f.write(
                        '@{} -m {} %1 %2 %3 %4 %5 %6 %7 %8 %9\n'.format(
                                sys.executable, module))
            else:
                script_f.write('#!/bin/sh\n')
                script_f.write(
                        'exec %s -m %s ${1+"$@"}\n' % (sys.executable,
                                module))

            script_f.close()

            # Make the script executable.
            os.chmod(script_path,
                    stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR|
                    stat.S_IRGRP|stat.S_IXGRP|
                    stat.S_IROTH|stat.S_IXOTH)

        project.installables.append(installable)
