# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import re

from . import Command
from ..config import Config
from ..machine import Machine, MachineCollection
from ..results import Results, get_filename
from ..benchmarks import Benchmarks
from ..console import log
from .. import util
from .run import Run


class Update(Command):
    @classmethod
    def setup_arguments(cls, subparsers):
        parser = subparsers.add_parser(
            "update", help="Update the results and config files "
            "to the current version",
            description="Update the results and config files "
            "to the current version")

        parser.set_defaults(func=cls.run_from_args)

        return parser

    @classmethod
    def run_from_args(cls, args, _machine_file=None):
        return cls.run(args.config, _machine_file=_machine_file)

    @classmethod
    def run(cls, config_path, _machine_file=None):
        MachineCollection.update(_path=_machine_file)

        conf = Config.load(config_path)

        log.info("Updating results data...")

        for root, dirs, files in os.walk(conf.results_dir):
            for filename in files:
                path = os.path.join(root, filename)
                if filename == 'machine.json':
                    Machine.update(path)
                elif filename == "benchmarks.json":
                    pass
                elif filename.endswith('.json'):
                    Results.update(path)

                    # Rename files if necessary
                    m = re.match(r'^([0-9a-f]+)-(.*)\.json$', os.path.basename(path), re.I)
                    if m:
                        new_path = get_filename(root, m.group(1), m.group(2))
                        if new_path != path:
                            try:
                                if os.path.exists(new_path):
                                    raise OSError()
                                os.rename(path, new_path)
                            except OSError:
                                log.warning("{}: should be renamed to {}".format(path, new_path))
                    else:
                        log.warning("{}: unrecognized file name".format(path))

        # Check benchmarks.json
        log.info("Updating benchmarks.json...")
        ok = False
        try:
            Benchmarks.load(conf)
            ok = True
        except util.UserError:
            pass

        if not ok:
            # Regenerating the file is needed
            with log.indent():
                Run.run(conf, bench=['just-discover'])
