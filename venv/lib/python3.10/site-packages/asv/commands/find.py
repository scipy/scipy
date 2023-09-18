# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math

from . import Command
from ..benchmarks import Benchmarks
from ..console import log
from ..machine import Machine
from ..repo import get_repo
from ..runner import run_benchmarks
from .. import util

from .setup import Setup

from . import common_args


def draw_graph(lo, mid, hi, total):
    nchars = 60
    scale = float(nchars) / total
    graph = ['-'] * nchars
    graph[int(lo * scale)] = '<'
    graph[int(hi * scale)] = '>'
    graph[int(mid * scale)] = 'O'
    return ''.join(graph)


class Find(Command):
    @classmethod
    def setup_arguments(cls, subparsers):
        parser = subparsers.add_parser(
            "find", help="Find commits that introduced large regressions",
            description="""Adaptively searches a range of commits for
            one that produces a large regression.  This only works well
            when the regression in the range is mostly monotonic.""")

        parser.add_argument(
            'range', type=str, metavar=('from..to',),
            help="""Range of commits to search.  For a git
            repository, this is passed as the first argument to ``git
            log``.  See 'specifying ranges' section of the
            `gitrevisions` manpage for more info.""")
        parser.add_argument(
            "bench", type=str, metavar=('benchmark_name',),
            help="""Name of benchmark to use in search.""")
        parser.add_argument(
            "--invert", "-i", action="store_true",
            help="""Search for a decrease in the benchmark value,
            rather than an increase.""")
        common_args.add_show_stderr(parser)
        common_args.add_machine(parser)
        common_args.add_environment(parser)
        common_args.add_launch_method(parser)

        parser.set_defaults(func=cls.run_from_args)

        return parser

    @classmethod
    def run_from_conf_args(cls, conf, args, **kwargs):
        return cls.run(
            conf, args.range, args.bench,
            invert=args.invert, show_stderr=args.show_stderr,
            machine=args.machine, env_spec=args.env_spec,
            launch_method=args.launch_method, **kwargs
        )

    @classmethod
    def run(cls, conf, range_spec, bench, invert=False, show_stderr=False,
            machine=None, env_spec=None, _machine_file=None, launch_method=None):
        params = {}
        machine_params = Machine.load(
            machine_name=machine,
            _path=_machine_file,
            interactive=True)
        params.update(machine_params.__dict__)
        machine_params.save(conf.results_dir)

        repo = get_repo(conf)
        repo.pull()

        commit_hashes = repo.get_hashes_from_range(range_spec)[::-1]

        if len(commit_hashes) == 0:
            log.error("No commit hashes selected")
            return 1

        environments = Setup.run(conf=conf, env_spec=env_spec)
        if len(environments) == 0:
            log.error("No environments selected")
            return 1

        benchmarks = Benchmarks.discover(conf, repo, environments,
                                         commit_hashes, regex=bench)
        if len(benchmarks) == 0:
            log.error("'{0}' benchmark not found".format(bench))
            return 1
        elif len(benchmarks) > 1:
            log.error("'{0}' matches more than one benchmark".format(bench))
            return 1

        benchmark_name, = benchmarks.keys()

        steps = int(math.log(len(commit_hashes)) / math.log(2))

        log.info(
            "Running approximately {0} benchmarks within {1} commits".format(
                steps, len(commit_hashes)))

        env = environments[0]

        results = [None] * len(commit_hashes)

        def do_benchmark(i):
            if results[i] is not None:
                return results[i]

            commit_hash = commit_hashes[i]

            commit_name = repo.get_decorated_hash(commit_hash, 8)
            log.info(
                "For {0} commit {1}:".format(
                    conf.project, commit_name))

            env.install_project(conf, repo, commit_hash)

            res = run_benchmarks(benchmarks, env,
                                 show_stderr=show_stderr,
                                 launch_method=launch_method)

            result = res.get_result_value(benchmark_name,
                                          benchmarks[benchmark_name]['params'])

            results[i] = result

            return results[i]

        def non_null_results(*results):
            """
            Whether some value is non-null in all result sets
            """
            for values in zip(*results):
                if all(x is not None for x in values):
                    return True
            return False

        def difference_3way(a, b, c):
            """
            Return largest regression (a-b, b-c).
            """
            results_ab = [0]
            results_bc = [0]
            for va, vb, vc in zip(a, b, c):
                if va is not None and vb is not None and vc is not None:
                    denom = abs(va) + abs(vb) + abs(vc)
                    if denom == 0:
                        denom = 1.0
                    results_ab.append((va - vb) / denom)
                    results_bc.append((vb - vc) / denom)
            return max(results_ab), max(results_bc)

        def do_search(lo, hi):
            if hi - lo <= 1:
                return hi

            mid = int(math.floor((hi - lo) / 2) + lo)

            log.info(
                "Testing {0}".format(
                    draw_graph(lo, mid, hi, len(commit_hashes))))

            with log.indent():
                lo_result = None
                while lo_result is None:
                    lo_result = do_benchmark(lo)
                    if not non_null_results(lo_result):
                        lo_result = None
                        lo += 1
                        if lo >= mid:
                            raise util.UserError("Too many commits failed")

                mid_result = None
                while mid_result is None:
                    mid_result = do_benchmark(mid)
                    if not non_null_results(mid_result, lo_result):
                        mid_result = None
                        mid += 1
                        if mid >= hi:
                            raise util.UserError("Too many commits failed")

                hi_result = None
                while hi_result is None:
                    hi_result = do_benchmark(hi)
                    if not non_null_results(lo_result, mid_result, hi_result):
                        hi_result = None
                        hi -= 1
                        if hi <= mid:
                            raise util.UserError("Too many commits failed")

            diff_b, diff_a = difference_3way(hi_result, mid_result, lo_result)

            if invert:
                diff_a *= -1.0
                diff_b *= -1.0

            if diff_a >= diff_b:
                return do_search(lo, mid)
            else:
                return do_search(mid, hi)

        result = do_search(0, len(commit_hashes) - 1)

        commit_name = repo.get_decorated_hash(commit_hashes[result], 8)
        log.info("Greatest regression found: {0}".format(commit_name))

        return 0
