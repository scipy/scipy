# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from fnmatch import fnmatchcase
import sys

from . import Command
from .. import console
from ..console import log
from ..results import iter_results
from . import util


class Rm(Command):
    @classmethod
    def setup_arguments(cls, subparsers):
        parser = subparsers.add_parser(
            "rm", help="Remove results from the database",
            description="""
            Removes entries from the results database.
            """)

        parser.add_argument(
            'patterns', nargs='+',
            help="""Pattern(s) to match, each of the form X=Y.  X may
            be one of "benchmark", "commit_hash", "python" or any of
            the machine or environment params.  Y is a case-sensitive
            glob pattern.""")
        parser.add_argument(
            "-y", action="store_true",
            help="""Don't prompt for confirmation.""")

        parser.set_defaults(func=cls.run_from_args)

        return parser

    @classmethod
    def run_from_conf_args(cls, conf, args):
        return cls.run(conf, args.patterns, args.y)

    @classmethod
    def run(cls, conf, patterns, y=True):
        global_patterns = {}
        single_benchmark = None
        files_to_remove = set()
        count = 0

        for pattern in patterns:
            parts = pattern.split('=', 1)
            if len(parts) != 2:
                raise util.UserError("Invalid pattern '{0}'".format(pattern))

            if parts[0] == 'benchmark':
                if single_benchmark is not None:
                    raise util.UserError("'benchmark' appears more than once")
                single_benchmark = parts[1]
            else:
                if parts[0] in global_patterns:
                    raise util.UserError(
                        "'{0}' appears more than once".format(parts[0]))
                global_patterns[parts[0]] = parts[1]

        for result in iter_results(conf.results_dir):
            found = True
            for key, val in six.iteritems(global_patterns):
                if key == 'commit_hash':
                    if not util.hash_equal(result.commit_hash, val):
                        found = False
                        break
                elif key == 'python':
                    if not fnmatchcase(result.env.python, val):
                        found = False
                        break
                else:
                    if not fnmatchcase(result.params.get(key), val):
                        found = False
                        break

            if not found:
                continue

            if single_benchmark is not None:
                found = False
                for benchmark in list(result.get_all_result_keys()):
                    if fnmatchcase(benchmark, single_benchmark):
                        count += 1
                        files_to_remove.add(result)
                        result.remove_result(benchmark)
            else:
                files_to_remove.add(result)

        if single_benchmark is not None:
            log.info("Removing {0} benchmarks in {1} files".format(
                count, len(files_to_remove)))
        else:
            log.info("Removing {0} files".format(len(files_to_remove)))

        if not y:
            do = console.get_answer_default("Perform operations", "n")
            if len(do) and do.lower()[0] != 'y':
                sys.exit(0)

        if single_benchmark is not None:
            for result in files_to_remove:
                result.save(conf.results_dir)
        else:
            for result in files_to_remove:
                result.rm(conf.results_dir)
