# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


class OutputPublisher(object):
    """
    A base class for pages displaying output in the JS application
    """
    name = None
    button_label = None
    description = None
    order = float('inf')

    @classmethod
    def publish(cls, conf, repo, benchmarks, graphs, revisions):
        pass
