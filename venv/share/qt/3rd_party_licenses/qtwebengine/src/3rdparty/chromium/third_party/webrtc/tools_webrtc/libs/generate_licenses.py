#!/usr/bin/env python

#  Copyright 2016 The WebRTC project authors. All Rights Reserved.
#
#  Use of this source code is governed by a BSD-style license
#  that can be found in the LICENSE file in the root of the source
#  tree. An additional intellectual property rights grant can be found
#  in the file PATENTS.  All contributing project authors may
#  be found in the AUTHORS file in the root of the source tree.
"""Generates license markdown for a prebuilt version of WebRTC.

Licenses are taken from dependent libraries which are determined by
GN desc command `gn desc` on all targets specified via `--target` argument.

One can see all dependencies by invoking this command:
$ gn.py desc --all --format=json <out_directory> <target> | python -m json.tool
(see "deps" subarray)

Libraries are mapped to licenses via LIB_TO_LICENSES_DICT dictionary.

"""

import sys

import argparse
import cgi
import json
import logging
import os
import re
import subprocess

# Third_party library to licences mapping. Keys are names of the libraries
# (right after the `third_party/` prefix)
LIB_TO_LICENSES_DICT = {
    'abseil-cpp': ['third_party/abseil-cpp/LICENSE'],
    'android_ndk': ['third_party/android_ndk/NOTICE'],
    'android_sdk': ['third_party/android_sdk/LICENSE'],
    'auto': ['third_party/android_deps/libs/'
             'com_google_auto_service_auto_service/LICENSE'],
    'bazel': ['third_party/bazel/LICENSE'],
    'boringssl': ['third_party/boringssl/src/LICENSE'],
    'errorprone': ['third_party/android_deps/libs/'
                   'com_google_errorprone_error_prone_core/LICENSE'],
    'fiat': ['third_party/boringssl/src/third_party/fiat/LICENSE'],
    'guava': ['third_party/guava/LICENSE'],
    'ijar': ['third_party/ijar/LICENSE'],
    'jsoncpp': ['third_party/jsoncpp/LICENSE'],
    'libaom': ['third_party/libaom/source/libaom/LICENSE'],
    'libc++': ['buildtools/third_party/libc++/trunk/LICENSE.TXT'],
    'libc++abi': ['buildtools/third_party/libc++abi/trunk/LICENSE.TXT'],
    'libevent': ['base/third_party/libevent/LICENSE'],
    'libjpeg_turbo': ['third_party/libjpeg_turbo/LICENSE.md'],
    'libsrtp': ['third_party/libsrtp/LICENSE'],
    'libvpx': ['third_party/libvpx/source/libvpx/LICENSE'],
    'libyuv': ['third_party/libyuv/LICENSE'],
    'nasm': ['third_party/nasm/LICENSE'],
    'opus': ['third_party/opus/src/COPYING'],
    'pffft': ['third_party/pffft/LICENSE'],
    'protobuf': ['third_party/protobuf/LICENSE'],
    'rnnoise': ['third_party/rnnoise/COPYING'],
    'usrsctp': ['third_party/usrsctp/LICENSE'],
    'webrtc': ['LICENSE'],
    'zlib': ['third_party/zlib/LICENSE'],
    'base64': ['rtc_base/third_party/base64/LICENSE'],
    'sigslot': ['rtc_base/third_party/sigslot/LICENSE'],
    'portaudio': ['modules/third_party/portaudio/LICENSE'],
    'fft': ['modules/third_party/fft/LICENSE'],
    'g711': ['modules/third_party/g711/LICENSE'],
    'g722': ['modules/third_party/g722/LICENSE'],
    'ooura': ['common_audio/third_party/ooura/LICENSE'],
    'spl_sqrt_floor': ['common_audio/third_party/spl_sqrt_floor/LICENSE'],

    # TODO(bugs.webrtc.org/1110): Remove this hack. This is not a lib.
    # For some reason it is listed as so in _GetThirdPartyLibraries.
    'android_deps': [],

    # Compile time dependencies, no license needed:
    'yasm': [],
    'ow2_asm': [],
    'jdk': [],
}

# Third_party library _regex_ to licences mapping. Keys are regular expression
# with names of the libraries (right after the `third_party/` prefix)
LIB_REGEX_TO_LICENSES_DICT = {
    'android_deps:android_support_annotations.*': [
        'third_party/android_deps/libs/' +
        'com_android_support_support_annotations/LICENSE'
    ],

    # Internal dependencies, licenses are already included by other dependencies
    'android_deps:com_android_support_support_annotations.*': [],
}


def FindSrcDirPath():
  """Returns the abs path to the src/ dir of the project."""
  src_dir = os.path.dirname(os.path.abspath(__file__))
  while os.path.basename(src_dir) != 'src':
    src_dir = os.path.normpath(os.path.join(src_dir, os.pardir))
  return src_dir


SCRIPT_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))
WEBRTC_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir, os.pardir))
SRC_DIR = FindSrcDirPath()
sys.path.append(os.path.join(SRC_DIR, 'build'))
import find_depot_tools

THIRD_PARTY_LIB_SIMPLE_NAME_REGEX = r'^.*/third_party/([\w\-+]+).*$'
THIRD_PARTY_LIB_REGEX_TEMPLATE = r'^.*/third_party/%s$'


class LicenseBuilder(object):

  def __init__(self,
               buildfile_dirs,
               targets,
               lib_to_licenses_dict=None,
               lib_regex_to_licenses_dict=None):
    if lib_to_licenses_dict is None:
      lib_to_licenses_dict = LIB_TO_LICENSES_DICT

    if lib_regex_to_licenses_dict is None:
      lib_regex_to_licenses_dict = LIB_REGEX_TO_LICENSES_DICT

    self.buildfile_dirs = buildfile_dirs
    self.targets = targets
    self.lib_to_licenses_dict = lib_to_licenses_dict
    self.lib_regex_to_licenses_dict = lib_regex_to_licenses_dict

    self.common_licenses_dict = self.lib_to_licenses_dict.copy()
    self.common_licenses_dict.update(self.lib_regex_to_licenses_dict)

  @staticmethod
  def _ParseLibraryName(dep):
    """Returns library name after third_party

    Input one of:
    //a/b/third_party/libname:c
    //a/b/third_party/libname:c(//d/e/f:g)
    //a/b/third_party/libname/c:d(//e/f/g:h)

    Outputs libname or None if this is not a third_party dependency.
    """
    groups = re.match(THIRD_PARTY_LIB_SIMPLE_NAME_REGEX, dep)
    return groups.group(1) if groups else None

  def _ParseLibrary(self, dep):
    """Returns library simple or regex name that matches `dep` after third_party

    This method matches `dep` dependency against simple names in
    LIB_TO_LICENSES_DICT and regular expression names in
    LIB_REGEX_TO_LICENSES_DICT keys

    Outputs matched dict key or None if this is not a third_party dependency.
    """
    libname = LicenseBuilder._ParseLibraryName(dep)

    for lib_regex in self.lib_regex_to_licenses_dict:
      if re.match(THIRD_PARTY_LIB_REGEX_TEMPLATE % lib_regex, dep):
        return lib_regex

    return libname

  @staticmethod
  def _RunGN(buildfile_dir, target):
    cmd = [
        sys.executable,
        os.path.join(find_depot_tools.DEPOT_TOOLS_PATH, 'gn.py'),
        'desc',
        '--all',
        '--format=json',
        os.path.abspath(buildfile_dir),
        target,
    ]
    logging.debug('Running: %r', cmd)
    output_json = subprocess.check_output(cmd, cwd=WEBRTC_ROOT)
    logging.debug('Output: %s', output_json)
    return output_json

  def _GetThirdPartyLibraries(self, buildfile_dir, target):
    output = json.loads(LicenseBuilder._RunGN(buildfile_dir, target))
    libraries = set()
    for described_target in output.values():
      third_party_libs = (
          self._ParseLibrary(dep) for dep in described_target['deps'])
      libraries |= set(lib for lib in third_party_libs if lib)
    return libraries

  def GenerateLicenseText(self, output_dir):
    # Get a list of third_party libs from gn. For fat libraries we must consider
    # all architectures, hence the multiple buildfile directories.
    third_party_libs = set()
    for buildfile in self.buildfile_dirs:
      for target in self.targets:
        third_party_libs |= self._GetThirdPartyLibraries(buildfile, target)
    assert len(third_party_libs) > 0

    missing_licenses = third_party_libs - set(self.common_licenses_dict.keys())
    if missing_licenses:
      error_msg = 'Missing licenses for following third_party targets: %s' % \
                  ', '.join(missing_licenses)
      logging.error(error_msg)
      raise Exception(error_msg)

    # Put webrtc at the front of the list.
    license_libs = sorted(third_party_libs)
    license_libs.insert(0, 'webrtc')

    logging.info('List of licenses: %s', ', '.join(license_libs))

    # Generate markdown.
    output_license_file = open(os.path.join(output_dir, 'LICENSE.md'), 'w+')
    for license_lib in license_libs:
      if len(self.common_licenses_dict[license_lib]) == 0:
        logging.info('Skipping compile time or internal dependency: %s',
                     license_lib)
        continue  # Compile time dependency

      output_license_file.write('# %s\n' % license_lib)
      output_license_file.write('```\n')
      for path in self.common_licenses_dict[license_lib]:
        license_path = os.path.join(WEBRTC_ROOT, path)
        with open(license_path, 'r') as license_file:
          license_text = cgi.escape(license_file.read(), quote=True)
          output_license_file.write(license_text)
          output_license_file.write('\n')
      output_license_file.write('```\n\n')

    output_license_file.close()


def main():
  parser = argparse.ArgumentParser(description='Generate WebRTC LICENSE.md')
  parser.add_argument(
      '--verbose', action='store_true', default=False, help='Debug logging.')
  parser.add_argument(
      '--target',
      required=True,
      action='append',
      default=[],
      help='Name of the GN target to generate a license for')
  parser.add_argument('output_dir', help='Directory to output LICENSE.md to.')
  parser.add_argument(
      'buildfile_dirs',
      nargs='+',
      help='Directories containing gn generated ninja files')
  args = parser.parse_args()

  logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)

  builder = LicenseBuilder(args.buildfile_dirs, args.target)
  builder.GenerateLicenseText(args.output_dir)


if __name__ == '__main__':
  sys.exit(main())
