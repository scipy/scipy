import re

from ..check import Check

_MIN_YEAR = 2009
_MAX_YEAR = 2030

_LICENSE_PATTERN = re.compile(
    r'\/\/ Copyright( \([cC]\))? (\w+) .+\. {1,2}All rights reserved\.[\r\n]{1,2}' +
    r'(' +
        r'\/\/ (' +
            r'This code is governed by the( BSD)? license found in the LICENSE file\.' +
            r'|' +
            r'See LICENSE for details' +
        r')' +
        r'|' +
        r'\/\/ Use of this source code is governed by a BSD-style license that can be[\r\n]{1,2}' +
        r'\/\/ found in the LICENSE file\.' +
        r'|' +
        r'\/\/ See LICENSE or https://github\.com/tc39/test262/blob/master/LICENSE' +
    r')', re.IGNORECASE)

class CheckLicense(Check):
    '''Ensure tests declare valid license information.'''
    ID = 'LICENSE'

    def run(self, name, meta, source):
        if meta and 'flags' in meta and 'generated' in meta['flags']:
            return

        match = _LICENSE_PATTERN.search(source)

        if not match:
            return 'No license information found.'

        year_str = match.group(2)
        try:
            year = int(year_str)

            if year < _MIN_YEAR or year > _MAX_YEAR:
                raise ValueError()
        except ValueError:
            return 'Invalid year: %s' % year_str
