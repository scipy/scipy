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
import shutil
import subprocess
import sys

from ..version import SIP_VERSION, SIP_VERSION_STR

from .abi_version import (get_module_source_dir, get_sip_module_version,
        resolve_abi_version)


def module(sip_module, abi_version, project, sdist, setup_cfg, sip_h, sip_rst,
        target_dir):
    """ Create the various elements of a sip module. """

    # Provide some defaults.
    abi_major_version = resolve_abi_version(abi_version).split('.')[0]

    if project is None:
        project = sip_module.replace('.', '_')

    # Create the patches.
    patches = _create_patches(sip_module, abi_major_version, project)

    # The names of generated files.
    sdist_dir = project + '-' + patches['@SIP_MODULE_VERSION@']
    sip_h_fn = 'sip.h'
    sip_rst_fn = 'sip.rst'

    if target_dir:
        sdist_dir = os.path.join(target_dir, sdist_dir)
        sip_h_fn = os.path.join(target_dir, sip_h_fn)
        sip_rst_fn = os.path.join(target_dir, sip_rst_fn)

    # Generate the required files.
    if sdist:
        _create_sdist(sdist_dir, abi_major_version, patches, setup_cfg)

    if sip_h:
        _create_sip_file(sip_h_fn, abi_major_version, patches)

    if sip_rst:
        _create_sip_file(sip_rst_fn, abi_major_version, patches)


def copy_sip_h(abi_major_version, target_dir, sip_module='',
        version_info=True):
    """ Copy the sip.h file. """

    patches = _create_patches(sip_module, abi_major_version,
            version_info=version_info)
    _install_source_file('sip.h', get_module_source_dir(abi_major_version),
            target_dir, patches)


def copy_sip_pyi(abi_major_version, target_dir):
    """ Copy the sip.pyi file. """

    shutil.copy(
            os.path.join(get_module_source_dir(abi_major_version), 'sip.pyi'),
            target_dir)


def copy_nonshared_sources(abi_major_version, target_dir):
    """ Copy the module sources as a non-shared module. """

    # Copy the patched sip.h.
    copy_sip_h(abi_major_version, target_dir)

    # Copy the remaining source code.
    module_source_dir = get_module_source_dir(abi_major_version)
    sources = []

    for fn in os.listdir(module_source_dir):
        if fn.endswith('.c') or fn.endswith('.cpp') or fn.endswith('.h'):
            src_fn = os.path.join(module_source_dir, fn)
            dst_fn = os.path.join(target_dir, fn)
            shutil.copyfile(src_fn, dst_fn)

            if not fn.endswith('.h'):
                sources.append(dst_fn)

    return sources


def _create_patches(sip_module, abi_major_version, project='',
        version_info=True):
    """ Return a dict of the patches. """

    sip_module_parts = sip_module.split('.')
    sip_module_package_name = '.'.join(sip_module_parts[:-1])
    sip_module_name = sip_module_parts[-1]

    # We special case this because this should be the only package requiring
    # the support.
    legacy = (sip_module == 'PyQt5.sip')

    if version_info:
        sip_version = SIP_VERSION
        sip_version_str = SIP_VERSION_STR
    else:
        sip_version = 0
        sip_version_str = ''

    return {
        # The public patches are those that might be needed in setup.cfg or any
        # automatically generated user documentation.
        '@SIP_MODULE_FQ_NAME@':         sip_module,
        '@SIP_MODULE_PROJECT_NAME@':    project,
        '@SIP_MODULE_PACKAGE_NAME@':    sip_module_package_name,
        '@SIP_MODULE_VERSION@':         get_sip_module_version(
                                                abi_major_version),

        # These are internal to sip.h.
        '@_SIP_MODULE_FQ_NAME@':        sip_module,
        '@_SIP_MODULE_NAME@':           sip_module_name,
        '@_SIP_MODULE_SHARED@':         '1' if sip_module else '0',
        '@_SIP_MODULE_ENTRY@':          'PyInit_' + sip_module_name,
        '@_SIP_MODULE_LEGACY@':         "1" if legacy else "0",
        '@_SIP_VERSION@':               hex(sip_version),
        '@_SIP_VERSION_STR@':           sip_version_str
    }


def _create_sdist(sdist_dir, abi_version, patches, setup_cfg):
    """ Create the sdist. """

    # Remove any existing source directory.
    shutil.rmtree(sdist_dir, ignore_errors=True)

    os.mkdir(sdist_dir)

    # The source directory doesn't have sub-directories.
    module_source_dir = get_module_source_dir(abi_version)

    for name in os.listdir(module_source_dir):
        if name in ('setup.cfg.in', 'sip.pyi', 'sip.rst.in'):
            continue

        if name != 'MANIFEST.in' and name.endswith('.in'):
            name = name[:-3]

            # Don't install the default README if we are not using the default
            # setup.cfg.
            if name != 'README' or setup_cfg is None:
                _install_source_file(name, module_source_dir, sdist_dir,
                        patches)
        else:
            shutil.copy(os.path.join(module_source_dir, name), sdist_dir)

    # Write setup.cfg as required.
    if setup_cfg is None:
        setup_cfg = os.path.join(module_source_dir, 'setup.cfg.in')

    _install_file(setup_cfg, os.path.join(sdist_dir, 'setup.cfg'), patches)

    # Create the sdist file using setuptools.  This means any user supplied
    # setup.cfg should be handled correctly.
    saved_cwd = os.getcwd()
    os.chdir(sdist_dir)

    subprocess.run(
            [sys.executable, 'setup.py', '--quiet', 'sdist', '--dist-dir',
                    '..'])

    os.chdir(saved_cwd)

    # Tidy up.
    shutil.rmtree(sdist_dir)


def _create_sip_file(sip_file_fn, abi_version, patches):
    """ Create a patched file from the module source directory. """

    dname, fname = os.path.split(os.path.abspath(sip_file_fn))

    _install_source_file(fname, get_module_source_dir(abi_version), dname,
            patches)


def _install_source_file(name, module_source_dir, target_dir, patches):
    """ Install a source file in a target directory. """

    _install_file(os.path.join(module_source_dir, name) + '.in',
            os.path.join(target_dir, name), patches)


def _install_file(name_in, name_out, patches):
    """ Install a file. """

    # Read the file.
    with open(name_in) as f:
        data = f.read()

    # Patch the file.
    for patch_name, patch in patches.items():
        data = data.replace(patch_name, patch)

    # Write the file.
    with open(name_out, 'w') as f:
        f.write(data)
