# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, unicode_literals, print_function

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


from .. import profiling
from .. import util


class SnakevizGui(profiling.ProfilerGui):
    name = 'snakeviz'
    description = "snakeviz http://jiffyclub.github.com/snakeviz/"

    @classmethod
    def is_available(cls):
        return util.has_command('snakeviz')

    @classmethod
    def open_profiler_gui(cls, profiler_file):
        command = util.which('snakeviz')

        return util.check_call(
            [command, profiler_file],
            valid_return_codes=(0, -15), timeout=None)
