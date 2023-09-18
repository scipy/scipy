# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import shutil
import multiprocessing
import datetime

import six

from . import Command
from ..benchmarks import Benchmarks
from ..console import log
from ..graph import GraphSet
from ..machine import iter_machine_files
from ..repo import get_repo
from ..results import iter_results
from ..publishing import OutputPublisher
from .. import statistics
from .. import util
from .. import __version__


def check_benchmark_params(name, benchmark):
    """
    Check benchmark params and param_keys items, so that the javascript can
    assume this data is valid. It is checked in benchmark.py already when it
    is generated, but best to double check in any case.
    """
    if 'params' not in benchmark:
        # Old-format benchmarks.json
        benchmark['params'] = []
        benchmark['param_names'] = []

    msg = "Information in benchmarks.json for benchmark %s is malformed" % (
        name)
    if (not isinstance(benchmark['params'], list) or
        not isinstance(benchmark['param_names'], list)):
        raise ValueError(msg)
    if len(benchmark['params']) != len(benchmark['param_names']):
        raise ValueError(msg)
    for item in benchmark['params']:
        if not isinstance(item, list):
            raise ValueError(msg)


class Publish(Command):
    @classmethod
    def setup_arguments(cls, subparsers):
        parser = subparsers.add_parser(
            "publish", help="Collate results into a website",
            description=
            """
            Collate all results into a website.  This website will be
            written to the ``html_dir`` given in the ``asv.conf.json``
            file, and may be served using any static web server.""")
        parser.add_argument(
            '--no-pull', action='store_true', dest='no_pull',
            help="Do not pull the repository")
        parser.add_argument(
            'range', nargs='?', default=None,
            help="""Optional commit range to consider""")
        parser.add_argument(
            '--html-dir', '-o', default=None, help=(
                "Optional output directory. Default is 'html_dir' "
                "from asv config"))

        parser.set_defaults(func=cls.run_from_args)

        return parser

    @classmethod
    def run_from_conf_args(cls, conf, args):
        if args.html_dir is not None:
            conf.html_dir = args.html_dir
        return cls.run(conf=conf, range_spec=args.range, pull=not args.no_pull)

    @staticmethod
    def iter_results(conf, repo, range_spec=None):
        if range_spec is not None:
            if isinstance(range_spec, list):
                hashes = range_spec
            else:
                hashes = repo.get_hashes_from_range(range_spec)
        else:
            hashes = None
        for result in iter_results(conf.results_dir):
            if hashes is None or result.commit_hash in hashes:
                yield result

    @classmethod
    def run(cls, conf, range_spec=None, pull=True):
        params = {}
        graphs = GraphSet()
        machines = {}
        benchmark_names = set()

        log.set_nitems(6 + len(list(util.iter_subclasses(OutputPublisher))))

        if os.path.exists(conf.html_dir):
            util.long_path_rmtree(conf.html_dir)

        repo = get_repo(conf)
        benchmarks = Benchmarks.load(conf)

        def copy_ignore(src, names):
            # Copy only *.js and *.css in vendor dir
            ignore = [fn for fn in names
                      if (os.path.basename(src).lower() == 'vendor' and
                          not fn.lower().endswith('.js') and
                          not fn.lower().endswith('.css'))]
            return ignore

        template_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), '..', 'www')
        shutil.copytree(template_dir, conf.html_dir, ignore=copy_ignore)

        # Ensure html_dir is writable even if template_dir is on a read-only FS
        os.chmod(conf.html_dir, 0o755)
        for (pre, ds, fs) in os.walk(conf.html_dir):
            for x in fs:
                os.chmod(os.path.join(pre, x), 0o644)
            for x in ds:
                os.chmod(os.path.join(pre, x), 0o755)

        log.step()
        log.info("Loading machine info")
        with log.indent():
            for path in iter_machine_files(conf.results_dir):
                d = util.load_json(path)
                machines[d['machine']] = d

        log.step()
        log.info("Getting params, commits, tags and branches")
        with log.indent():
            # Determine first the set of all parameters and all commits
            hash_to_date = {}
            for results in cls.iter_results(conf, repo, range_spec):
                hash_to_date[results.commit_hash] = results.date
                for key, val in six.iteritems(results.params):
                    if val is None:
                        # Backward compatibility -- null means ''
                        val = ''

                    params.setdefault(key, set())
                    params[key].add(val)

            if pull:
                repo.pull()
            tags = repo.get_tags()
            revisions = repo.get_revisions(set(hash_to_date.keys()) | set(tags.values()))

            for tag, commit_hash in list(tags.items()):
                # Map to revision number instead of commit hash and add tags to hash_to_date
                tags[tag] = revisions[tags[tag]]
                hash_to_date[commit_hash] = repo.get_date_from_name(commit_hash)

            revision_to_date = dict((r, hash_to_date[h]) for h, r in six.iteritems(revisions))

            branches = dict(
                (branch, repo.get_branch_commits(branch))
                for branch in conf.branches)

        log.step()
        log.info("Loading results")
        with log.indent():
            # Generate all graphs
            for results in cls.iter_results(conf, repo, range_spec):
                log.dot()

                branches_for_commit = [branch for branch, commits in branches.items() if
                                       results.commit_hash in commits]

                # Print a warning message if we couldn't find the branch of a commit
                if not len(branches_for_commit):
                    msg = "Couldn't find %s in %s branches"
                    log.warning(msg % (results.commit_hash, branches.keys()))

                for key in results.get_result_keys(benchmarks):
                    b = benchmarks[key]
                    b_params = b['params']

                    result = results.get_result_value(key, b_params)
                    weight = [statistics.get_weight(s)
                              for s in results.get_result_stats(key, b_params)]
                    if not b_params:
                        result = result[0]
                        weight = weight[0]

                    benchmark_names.add(key)

                    for branch in branches_for_commit:
                        cur_params = dict(results.params)
                        cur_params['branch'] = repo.get_branch_name(branch)

                        # Backward compatibility, see above
                        for param_key, param_value in list(cur_params.items()):
                            if param_value is None:
                                cur_params[param_key] = ''

                        # Fill in missing params
                        for param_key in params.keys():
                            if param_key not in cur_params:
                                cur_params[param_key] = None
                                params[param_key].add(None)

                        # Create graph
                        graph = graphs.get_graph(key, cur_params)
                        graph.add_data_point(revisions[results.commit_hash], result, weight)

            # Get the parameter sets for all graphs
            graph_param_list = []
            for path, graph in graphs:
                if 'summary' not in graph.params:
                    if graph.params not in graph_param_list:
                        graph_param_list.append(graph.params)

        log.step()
        log.info("Detecting steps")
        with log.indent():
            n_processes = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(n_processes)
            try:
                graphs.detect_steps(pool, dots=log.dot)
                pool.close()
                pool.join()
            finally:
                pool.terminate()

        log.step()
        log.info("Generating graphs")
        with log.indent():
            # Save files
            graphs.save(conf.html_dir, dots=log.dot)

        pages = []
        classes = sorted(util.iter_subclasses(OutputPublisher),
                         key=lambda cls: cls.order)
        for cls in classes:
            log.step()
            log.info("Generating output for {0}".format(cls.__name__))
            with log.indent():
                cls.publish(conf, repo, benchmarks, graphs, revisions)
                pages.append([cls.name, cls.button_label, cls.description])

        log.step()
        log.info("Writing index")
        benchmark_map = dict(benchmarks)
        for key in six.iterkeys(benchmark_map):
            check_benchmark_params(key, benchmark_map[key])
        for key, val in six.iteritems(params):
            val = list(val)
            val.sort(key=lambda x: '[none]' if x is None else str(x))
            params[key] = val
        params['branch'] = [repo.get_branch_name(branch) for branch in conf.branches]
        revision_to_hash = dict((r, h) for h, r in six.iteritems(revisions))
        util.write_json(os.path.join(conf.html_dir, "index.json"), {
            'project': conf.project,
            'project_url': conf.project_url,
            'show_commit_url': conf.show_commit_url,
            'hash_length': conf.hash_length,
            'revision_to_hash': revision_to_hash,
            'revision_to_date': revision_to_date,
            'params': params,
            'graph_param_list': graph_param_list,
            'benchmarks': benchmark_map,
            'machines': machines,
            'tags': tags,
            'pages': pages,
        }, compact=True)

        util.write_json(os.path.join(conf.html_dir, "info.json"), {
            'asv-version': __version__,
            'timestamp': util.datetime_to_js_timestamp(datetime.datetime.utcnow())
        })
