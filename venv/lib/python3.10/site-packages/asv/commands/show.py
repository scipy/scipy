# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from collections import defaultdict

from . import Command
from ..benchmarks import Benchmarks
from ..machine import iter_machine_files
from ..results import iter_results_for_machine, iter_results_for_machine_and_hash
from ..runner import format_benchmark_result
from ..repo import get_repo, NoSuchNameError
from ..util import load_json
from ..console import log, color_print
from ..environment import get_environments
from .. import util

from . import common_args


class Show(Command):
    @classmethod
    def setup_arguments(cls, subparsers):
        parser = subparsers.add_parser(
            "show", help="Print recorded data.",
            description="""Print saved benchmark results.""")

        parser.add_argument(
            'commit', nargs='?', default=None,
            help="""The commit to show data for.""")
        parser.add_argument(
            '--details', action='store_true', default=False,
            help="""Show all result details.""")
        common_args.add_bench(parser)
        common_args.add_machine(parser)
        common_args.add_environment(parser)
        parser.set_defaults(func=cls.run_from_args)

        return parser

    @classmethod
    def run_from_conf_args(cls, conf, args, **kwargs):
        return cls.run(
            conf=conf, commit=args.commit, bench=args.bench,
            machine=args.machine, env_spec=args.env_spec,
            details=args.details, **kwargs
        )

    @classmethod
    def run(cls, conf, commit=None, bench=None, machine=None, env_spec=None,
            details=False):
        if env_spec:
            env_names = ([env.name for env in get_environments(conf, env_spec, verbose=False)]
                         + list(env_spec))
        else:
            env_names = None

        machines = []
        for path in iter_machine_files(conf.results_dir):
            d = load_json(path)
            machines.append(d['machine'])

        if len(machines) == 0:
            raise util.UserError("No results found")
        elif machine is None:
            pass
        elif machine in machines:
            machines = [machine]
        else:
            raise util.UserError(
                "Results for machine '{0} not found".format(machine))

        benchmarks = Benchmarks.load(conf, regex=bench)

        if commit is None:
            cls._print_commits(conf, machines, env_names, benchmarks)
        else:
            cls._print_results(conf, commit, machines, env_names, benchmarks,
                               show_details=details)

    @classmethod
    def _print_commits(cls, conf, machines, env_names, benchmarks):
        commits = defaultdict(lambda: {})

        for machine in machines:
            for result in iter_results_for_machine(
                    conf.results_dir, machine):
                if env_names is not None and result.env_name not in env_names:
                    continue

                if result.get_result_keys(benchmarks):
                    commits[(machine, result.env_name)][result.commit_hash] = result.date

        log.flush()

        color_print("Commits with results:")
        color_print("")

        for machine, env_name in sorted(commits.keys()):
            color_print("Machine    : {}".format(machine))
            color_print("Environment: {}".format(env_name))
            color_print("")

            cur_commits = commits[(machine, env_name)]
            commit_order = list(cur_commits.keys())
            commit_order.sort(key=lambda x: cur_commits[x])

            for commit in commit_order:
                color_print("    {}".format(commit[:conf.hash_length]))

            color_print("")

    @classmethod
    def _print_results(cls, conf, commit_hash, machines, env_names, benchmarks,
                       show_details=False):
        repo = get_repo(conf)
        try:
            commit_hash = repo.get_hash_from_name(commit_hash)
        except NoSuchNameError:
            pass

        def results_iter():
            for machine in sorted(machines):
                for result in iter_results_for_machine_and_hash(
                        conf.results_dir, machine, commit_hash):
                    if env_names is not None and result.env_name not in env_names:
                        continue
                    yield machine, result

        color_print("Commit: {}".format(repo.get_decorated_hash(commit_hash,
                                                                conf.hash_length)),
                    "blue")
        color_print("")

        for machine, result in results_iter():
            for name in sorted(result.get_result_keys(benchmarks)):
                cls._print_benchmark(machine, result, benchmarks[name],
                                     show_details=show_details)

    @classmethod
    def _print_benchmark(cls, machine, result, benchmark, show_details=False):
        color_print("{} [{}/{}]".format(benchmark['name'],
                                        machine,
                                        result.env_name),
                    'green')

        info, details = format_benchmark_result(result, benchmark)
        color_print("  {}".format(info), 'red')
        if details:
            color_print("  " + details.replace("\n", "\n  "))

        started_at = result.started_at.get(benchmark['name'])
        ended_at = result.ended_at.get(benchmark['name'])
        if started_at and ended_at:
            started_at = util.js_timestamp_to_datetime(started_at)
            ended_at = util.js_timestamp_to_datetime(ended_at)

            color_print('  started: {}, duration: {}'.format(
                started_at.strftime('%Y-%m-%d %H:%M:%S'),
                util.human_time((ended_at - started_at).total_seconds())))

        if not show_details:
            color_print("")
            return

        stats = result.get_result_stats(benchmark['name'], benchmark['params'])

        def get_stat_info(key):
            return [x.get(key) if x is not None else None for x in stats]

        for key in ['repeat', 'number', 'ci_99', 'mean', 'std', 'min', 'max']:
            values = get_stat_info(key)

            if key == 'ci_99':
                values = ["({}, {})".format(util.human_value(x[0], benchmark['unit']),
                                            util.human_value(x[1], benchmark['unit']))
                          if x is not None else None
                          for x in values]
            elif any(isinstance(x, float) for x in values):
                values = [util.human_value(x, benchmark['unit']) if x is not None else None
                          for x in values]

            if not all(x is None for x in values):
                color_print("  {}: {}".format(key, ", ".join(map(str, values))))

        samples = result.get_result_samples(benchmark['name'], benchmark['params'])
        if not all(x is None for x in samples):
            color_print("  samples: {}".format(samples))

        color_print("")
