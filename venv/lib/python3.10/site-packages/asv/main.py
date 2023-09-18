# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import sys

import six

from . import commands
from .config import Config
from .console import log
from .plugin_manager import plugin_manager
from . import util


def main():
    parser, subparsers = commands.make_argparser()

    args = parser.parse_args()

    if not hasattr(args, 'func'):
        parser.print_help()
        sys.exit(1)

    log.enable(args.verbose)

    args.config = os.path.abspath(args.config)

    # Use the path to the config file as the cwd for the remainder of
    # the run
    dirname = os.path.dirname(args.config)
    os.chdir(dirname)

    try:
        result = args.func(args)
    except util.UserError as e:
        log.error(six.text_type(e))
        sys.exit(1)
    finally:
        log.flush()

    if result is None:
        result = 0

    sys.exit(result)
