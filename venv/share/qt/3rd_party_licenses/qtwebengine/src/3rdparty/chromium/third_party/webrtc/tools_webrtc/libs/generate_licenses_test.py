#!/usr/bin/env vpython
# pylint: disable=relative-import,protected-access,unused-argument

#  Copyright 2017 The WebRTC project authors. All Rights Reserved.
#
#  Use of this source code is governed by a BSD-style license
#  that can be found in the LICENSE file in the root of the source
#  tree. An additional intellectual property rights grant can be found
#  in the file PATENTS.  All contributing project authors may
#  be found in the AUTHORS file in the root of the source tree.

import unittest
import mock

from generate_licenses import LicenseBuilder


class TestLicenseBuilder(unittest.TestCase):

  @staticmethod
  def _FakeRunGN(buildfile_dir, target):
    return """
    {
      "target1": {
        "deps": [
          "//a/b/third_party/libname1:c",
          "//a/b/third_party/libname2:c(//d/e/f:g)",
          "//a/b/third_party/libname3/c:d(//e/f/g:h)",
          "//a/b/not_third_party/c"
        ]
      }
    }
    """

  def testParseLibraryName(self):
    self.assertEquals(
        LicenseBuilder._ParseLibraryName('//a/b/third_party/libname1:c'),
        'libname1')
    self.assertEquals(
        LicenseBuilder._ParseLibraryName('//a/b/third_party/libname2:c(d)'),
        'libname2')
    self.assertEquals(
        LicenseBuilder._ParseLibraryName('//a/b/third_party/libname3/c:d(e)'),
        'libname3')
    self.assertEquals(
        LicenseBuilder._ParseLibraryName('//a/b/not_third_party/c'), None)

  def testParseLibrarySimpleMatch(self):
    builder = LicenseBuilder([], [], {}, {})
    self.assertEquals(
        builder._ParseLibrary('//a/b/third_party/libname:c'), 'libname')

  def testParseLibraryRegExNoMatchFallbacksToDefaultLibname(self):
    lib_dict = {
        'libname:foo.*': ['path/to/LICENSE'],
    }
    builder = LicenseBuilder([], [], lib_dict, {})
    self.assertEquals(
        builder._ParseLibrary('//a/b/third_party/libname:bar_java'), 'libname')

  def testParseLibraryRegExMatch(self):
    lib_regex_dict = {
        'libname:foo.*': ['path/to/LICENSE'],
    }
    builder = LicenseBuilder([], [], {}, lib_regex_dict)
    self.assertEquals(
        builder._ParseLibrary('//a/b/third_party/libname:foo_bar_java'),
        'libname:foo.*')

  def testParseLibraryRegExMatchWithSubDirectory(self):
    lib_regex_dict = {
        'libname/foo:bar.*': ['path/to/LICENSE'],
    }
    builder = LicenseBuilder([], [], {}, lib_regex_dict)
    self.assertEquals(
        builder._ParseLibrary('//a/b/third_party/libname/foo:bar_java'),
        'libname/foo:bar.*')

  def testParseLibraryRegExMatchWithStarInside(self):
    lib_regex_dict = {
        'libname/foo.*bar.*': ['path/to/LICENSE'],
    }
    builder = LicenseBuilder([], [], {}, lib_regex_dict)
    self.assertEquals(
        builder._ParseLibrary('//a/b/third_party/libname/fooHAHA:bar_java'),
        'libname/foo.*bar.*')

  @mock.patch('generate_licenses.LicenseBuilder._RunGN', _FakeRunGN)
  def testGetThirdPartyLibrariesWithoutRegex(self):
    builder = LicenseBuilder([], [], {}, {})
    self.assertEquals(
        builder._GetThirdPartyLibraries('out/arm', 'target1'),
        set(['libname1', 'libname2', 'libname3']))

  @mock.patch('generate_licenses.LicenseBuilder._RunGN', _FakeRunGN)
  def testGetThirdPartyLibrariesWithRegex(self):
    lib_regex_dict = {
        'libname2:c.*': ['path/to/LICENSE'],
    }
    builder = LicenseBuilder([], [], {}, lib_regex_dict)
    self.assertEquals(
        builder._GetThirdPartyLibraries('out/arm', 'target1'),
        set(['libname1', 'libname2:c.*', 'libname3']))

  @mock.patch('generate_licenses.LicenseBuilder._RunGN', _FakeRunGN)
  def testGenerateLicenseTextFailIfUnknownLibrary(self):
    lib_dict = {
        'simple_library': ['path/to/LICENSE'],
    }
    builder = LicenseBuilder(['dummy_dir'], ['dummy_target'], lib_dict, {})

    with self.assertRaises(Exception) as context:
      builder.GenerateLicenseText('dummy/dir')

    self.assertEquals(
        context.exception.message,
        'Missing licenses for following third_party targets: '
        'libname1, libname2, libname3')


if __name__ == '__main__':
  unittest.main()
