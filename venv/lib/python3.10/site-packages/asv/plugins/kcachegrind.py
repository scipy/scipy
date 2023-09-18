# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


from .. import profiling
from .. import util


class KCachegrindGui(profiling.ProfilerGui):
    name = 'kcachegrind'
    description = "kcachegrind through pyprof2calltree"

    @classmethod
    def is_available(cls):
        return (
            util.has_command("kcachegrind") and
            util.has_command("pyprof2calltree"))

    @classmethod
    def open_profiler_gui(cls, profiler_file):
        command = util.which("pyprof2calltree")

        return util.check_call(
            [command, '-i', profiler_file, '-k'], timeout=None)
