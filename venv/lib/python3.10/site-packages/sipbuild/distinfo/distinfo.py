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


import base64
import hashlib
import os
import shutil
import sys

from ..exceptions import UserException
from ..pyproject import PyProject
from ..version import SIP_VERSION_STR


# The wheel format defined in PEP 427.
WHEEL_VERSION = '1.0'


def distinfo(name, console_scripts, gui_scripts, generator, generator_version,
        inventory, metadata_overrides, prefix, project_root, requires_dists,
        wheel_tag):
    """ Create and populate a .dist-info directory from an inventory file. """

    if prefix is None:
        prefix = ''

    # Read the list of installed files.
    with open(inventory) as inventory_f:
        installed_lines = inventory_f.read().strip()
        installed = installed_lines.split('\n') if installed_lines else []

    # Get the pyproject.toml file.
    saved = os.getcwd()
    os.chdir(project_root)
    pyproject = PyProject()
    os.chdir(saved)

    # Get the metadata and update it from the command line.
    metadata = pyproject.get_metadata()

    if metadata_overrides is not None:
        for oride in metadata_overrides:
            parts = oride.split('=', maxsplit=1)
            or_name = parts[0].strip()
            or_value = parts[1].strip() if len(parts) == 2 else ''
            metadata[or_name] = or_value

    # Create the directory.
    create_distinfo(name, wheel_tag, installed, metadata, requires_dists,
            project_root, console_scripts, gui_scripts, prefix_dir=prefix,
            generator=generator, generator_version=generator_version)


def create_distinfo(distinfo_dir, wheel_tag, installed, metadata,
        requires_dists, project_root, console_scripts, gui_scripts,
        prefix_dir='', generator=None, generator_version=None):
    """ Create and populate a .dist-info directory. """

    if generator is None:
        generator = 'sipbuild'

    if generator_version is None:
        generator_version = SIP_VERSION_STR

    # The prefix directory corresponds to DESTDIR or INSTALL_ROOT.
    real_distinfo_dir = prefix_dir + distinfo_dir

    # Make sure we have an empty dist-info directory.  Handle exceptions as the
    # user may be trying something silly with a system directory.
    if os.path.exists(real_distinfo_dir):
        try:
            shutil.rmtree(real_distinfo_dir)
        except Exception as e:
            raise UserException(
                    "unable remove old dist-info directory '{}'".format(
                            real_distinfo_dir),
                    str(e))

    try:
        os.mkdir(real_distinfo_dir)
    except Exception as e:
        raise UserException(
                "unable create dist-info directory '{}'".format(
                        real_distinfo_dir),
                str(e))

    # Reproducable builds.
    installed.sort()

    if wheel_tag is None:
        # Create the INSTALLER file.
        installer_fn = os.path.join(distinfo_dir, 'INSTALLER')
        installed.append(installer_fn)

        with open(prefix_dir + installer_fn, 'w') as installer_f:
            print(generator, file=installer_f)
    else:
        # Define any entry points.
        if console_scripts or gui_scripts:
            eps_fn = os.path.join(distinfo_dir, 'entry_points.txt')
            installed.append(eps_fn)

            with open(prefix_dir + eps_fn, 'w') as eps_f:
                if console_scripts:
                    eps_f.write(
                            '[console_scripts]\n' + '\n'.join(
                                    console_scripts) + '\n')

                if gui_scripts:
                    eps_f.write(
                            '[gui_scripts]\n' + '\n'.join(gui_scripts) + '\n')

        # Create the WHEEL file.
        WHEEL = '''Wheel-Version: {}
Generator: {} {}
Root-Is-Purelib: false
Tag: {}
'''

        wheel_fn = os.path.join(distinfo_dir, 'WHEEL')
        installed.append(wheel_fn)

        with open(prefix_dir + wheel_fn, 'w') as wheel_f:
            wheel_f.write(
                    WHEEL.format(WHEEL_VERSION, generator, generator_version,
                            wheel_tag))

    # Create the METADATA file.
    metadata_fn = os.path.join(distinfo_dir, 'METADATA')
    write_metadata(metadata, requires_dists, metadata_fn, project_root,
            prefix_dir=prefix_dir)
    installed.append(metadata_fn)

    # Create the RECORD file.
    record_fn = os.path.join(distinfo_dir, 'RECORD')

    distinfo_path, distinfo_base = os.path.split(distinfo_dir)
    real_distinfo_path = os.path.normcase(prefix_dir + distinfo_path)

    with open(prefix_dir + record_fn, 'w') as record_f:
        for name in installed:
            real_name = prefix_dir + name
            if os.path.isdir(real_name):
                all_fns = []

                for root, dirs, files in os.walk(real_name):
                    # Reproducable builds.
                    dirs.sort()
                    files.sort()

                    for f in files:
                        all_fns.append(os.path.join(root, f))

                    if '__pycache__' in dirs:
                        dirs.remove('__pycache__')
            else:
                all_fns = [real_name]

            for fn in all_fns:
                norm_fn = os.path.normcase(fn)

                if norm_fn.startswith(real_distinfo_path):
                    fn_name = fn[len(real_distinfo_path) + 1:].replace('\\', '/')
                elif norm_fn.startswith(prefix_dir + sys.prefix):
                    fn_name = os.path.relpath(
                            fn, real_distinfo_path).replace('\\', '/')
                else:
                    fn_name = fn[len(prefix_dir):]

                fn_f = open(fn, 'rb')
                data = fn_f.read()
                fn_f.close()

                digest = base64.urlsafe_b64encode(
                        hashlib.sha256(data).digest()).rstrip(b'=').decode('ascii')

                record_f.write(
                        '{},sha256={},{}\n'.format(fn_name, digest, len(data)))

        record_f.write('{}/RECORD,,\n'.format(distinfo_base))


def write_metadata(metadata, requires_dists, metadata_fn, project_root,
        prefix_dir=''):
    """ Write the meta-data, with additional requirements to a file. """

    if requires_dists:
        rd = metadata.get('requires-dist', [])
        if isinstance(rd, str):
            rd = [rd]

        metadata['requires-dist'] = requires_dists + rd

    with open(prefix_dir + metadata_fn, 'w') as metadata_f:
        description = None

        # Do these first for cosmetic reasons.
        for name in ('metadata-version', 'name', 'version', 'requires-python'):
            _write_metadata_item(name, metadata.pop(name), metadata_f)

        for name, value in metadata.items():
            if name == 'description-file':
                description = value
            else:
                _write_metadata_item(name, value, metadata_f)

        if description is not None:
            metadata_f.write('\n')

            # The description file uses posix separators.
            description = description.replace('/', os.sep)

            with open(os.path.join(project_root, description)) as description_f:
                metadata_f.write(description_f.read())


def _write_metadata_item(name, value, metadata_f):
    """ Write a single metadata item. """

    if isinstance(value, str):
        value = [value]

    for v in value:
        metadata_f.write('{}: {}\n'.format(name.title(), v))
