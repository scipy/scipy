# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import base64
import os
import re
import zlib
import itertools
import hashlib
import datetime

import six
from six.moves import zip as izip

from . import environment
from .console import log
from .machine import Machine
from . import statistics
from . import util


def iter_results_paths(results):
    """
    Iterate over all of the result file paths.
    """
    skip_files = set([
        'machine.json', 'benchmarks.json'
    ])
    for root, dirs, files in os.walk(results):
        # Iterate over files only if machine.json is valid json
        machine_json = os.path.join(root, "machine.json")
        try:
            data = util.load_json(machine_json, api_version=Machine.api_version)
            machine_name = data.get('machine')
            if not isinstance(machine_name, six.text_type):
                raise util.UserError("malformed {0}".format(machine_json))
        except util.UserError as err:
            machine_json_err = "Skipping results: {0}".format(six.text_type(err))
        except IOError as err:
            machine_json_err = "Skipping results: could not load {0}".format(
                machine_json)
        else:
            machine_json_err = None

        # Iterate over files
        for filename in files:
            if filename not in skip_files and filename.endswith('.json'):
                if machine_json_err is not None:
                    # Show the warning only if there are some files to load
                    log.warning(machine_json_err)
                    break

                yield (root, filename, machine_name)


def iter_results(results):
    """
    Iterate over all of the result files.
    """
    for (root, filename, machine_name) in iter_results_paths(results):
        try:
            yield Results.load(os.path.join(root, filename), machine_name=machine_name)
        except util.UserError as exc:
            log.warning(six.text_type(exc))


def iter_results_for_machine(results, machine_name):
    """
    Iterate over all of the result files for a particular machine.
    """
    return iter_results(os.path.join(results, machine_name))


def iter_results_for_machine_and_hash(results, machine_name, commit):
    """
    Iterate over all of the result files with a given hash for a
    particular machine.
    """
    full_commit = get_result_hash_from_prefix(results, machine_name, commit)

    for (root, filename, machine_name) in iter_results_paths(
            os.path.join(results, machine_name)):
        results_commit = filename.split('-')[0]
        if results_commit == full_commit:
            try:
                yield Results.load(os.path.join(root, filename), machine_name=machine_name)
            except util.UserError as exc:
                log.warning(six.text_type(exc))


def iter_existing_hashes(results):
    """
    Iterate over all of the result commit hashes and dates and yields
    commit_hash.

    May return duplicates.  Use `get_existing_hashes` if that matters.
    """
    for result in iter_results(results):
        yield result.commit_hash


def get_existing_hashes(results):
    """
    Get a list of the commit hashes that have already been tested.
    """
    log.info("Getting existing hashes")
    hashes = list(set(iter_existing_hashes(results)))
    return hashes


def get_result_hash_from_prefix(results, machine_name, commit_prefix):
    """
    Get the 8-char result commit identifier from a potentially shorter
    prefix. Only considers the set of commits that have had
    results computed.

    Returns None if there are no matches. Raises a UserError
    if the prefix is non-unique.
    """
    commits = set([])

    path = os.path.join(results, machine_name)

    for (root, filename, r_machine_name) in iter_results_paths(path):
        if r_machine_name != machine_name:
            log.warning("Skipping results '{0}': machine name is not '{1}'".format(
                os.path.join(root, filename), machine_name))
            continue

        results_commit = filename.split('-')[0]
        cmp_len = min(len(commit_prefix), len(results_commit))
        if results_commit[:cmp_len] == commit_prefix[:cmp_len]:
            commits.add(results_commit)

    if len(commits) > 1:
        commit_list_str = ', '.join(sorted(commits))
        raise util.UserError('Git hash prefix could represent one of ' +
                             'multiple commits: {0}'.format(commit_list_str))
    elif len(commits) == 1:
        return list(commits)[0]
    else:
        return None


def get_filename(machine, commit_hash, env_name):
    """
    Get the result filename for a given machine, commit_hash and
    environment.

    If the environment name is too long, use its hash instead.
    """
    if env_name and len(env_name) >= 128:
        env_name = "env-" + hashlib.md5(env_name.encode('utf-8')).hexdigest()

    return os.path.join(
        machine,
        "{0}-{1}.json".format(
            commit_hash[:8],
            env_name))


def _compatible_results(result, result_params, params):
    """
    For parameterized benchmarks, obtain values from *result* that
    are compatible with parameters of *benchmark*
    """
    if result is None:
        # All results missing, eg. build failure
        return [None for param in itertools.product(*params)]

    # Pick results for those parameters that also appear in the
    # current benchmark
    old_results = {}
    for param, value in izip(itertools.product(*result_params), result):
        old_results[param] = value

    new_results = []
    for param in itertools.product(*params):
        new_results.append(old_results.get(param))

    return new_results


class Results(object):
    """
    Manage a set of benchmark results for a single machine and commit
    hash.
    """
    api_version = 1

    def __init__(self, params, requirements, commit_hash, date, python, env_name):
        """
        Parameters
        ----------
        params : dict
            Parameters describing the environment in which the
            benchmarks were run.

        requirements : list
            Requirements of the benchmarks being run.

        commit_hash : str
            The commit hash for the benchmark run.

        date : int
            JavaScript timestamp for when the commit was merged into
            the repository.

        python : str
            A Python version specifier.

        env_name : str
            Environment name
        """
        self._params = params
        self._requirements = requirements
        self._commit_hash = commit_hash
        self._date = date
        self._results = {}
        self._samples = {}
        self._stats = {}
        self._benchmark_params = {}
        self._profiles = {}
        self._python = python
        self._env_name = env_name
        self._started_at = {}
        self._ended_at = {}
        self._benchmark_version = {}

        # Note: stderr and errcode are not saved to files
        self._stderr = {}
        self._errcode = {}

        if commit_hash is not None:
            self._filename = get_filename(
                params['machine'], self._commit_hash, env_name)
        else:
            self._filename = None

    @classmethod
    def unnamed(cls):
        return cls({}, {}, None, None, None, None)

    @property
    def commit_hash(self):
        return self._commit_hash

    @property
    def date(self):
        return self._date

    @property
    def params(self):
        return self._params

    @property
    def started_at(self):
        return self._started_at

    @property
    def ended_at(self):
        return self._ended_at

    @property
    def benchmark_version(self):
        return self._benchmark_version

    @property
    def stderr(self):
        return self._stderr

    @property
    def errcode(self):
        return self._errcode

    def get_all_result_keys(self):
        """
        Return all available result keys.
        """
        return six.iterkeys(self._results)

    def get_result_keys(self, benchmarks):
        """
        Return result keys corresponding to benchmarks.

        Parameters
        ----------
        benchmarks : Benchmarks
            Benchmarks to return results for.
            Used for checking benchmark versions.

        Returns
        -------
        keys : set
            Set of benchmark result keys

        """
        keys = set()
        for key in six.iterkeys(self._results):
            if key not in benchmarks:
                continue

            version = self._benchmark_version.get(key)
            bench_version = benchmarks[key].get('version')

            if version is not None and version != bench_version:
                continue

            keys.add(key)

        return keys

    def get_result_value(self, key, params):
        """
        Return the value of benchmark result.

        Parameters
        ----------
        key : str
            Benchmark name to return results for
        params : {list of list, None}
            Set of benchmark parameters to return values for

        Returns
        -------
        value : {float, list of float}
            Benchmark result value. If the benchmark is parameterized, return
            a list of values.
        """
        return _compatible_results(self._results[key],
                                   self._benchmark_params[key],
                                   params)

    def get_result_stats(self, key, params):
        """
        Return the statistical information of a benchmark result.

        Parameters
        ----------
        key : str
            Benchmark name to return results for
        params : {list of list, None}
            Set of benchmark parameters to return values for

        Returns
        -------
        stats : {None, dict, list of dict}
            Result statistics. If the benchmark is parameterized,
            return a list of values.
        """
        return _compatible_results(self._stats[key],
                                   self._benchmark_params[key],
                                   params)

    def get_result_samples(self, key, params):
        """
        Return the raw data points of a benchmark result.

        Parameters
        ----------
        key : str
            Benchmark name to return results for
        params : {list of list, None}
            Set of benchmark parameters to return values for

        Returns
        -------
        samples : {None, list}
            Raw result samples. If the benchmark is parameterized,
            return a list of values.

        """
        return _compatible_results(self._samples[key],
                                   self._benchmark_params[key],
                                   params)

    def get_result_params(self, key):
        """
        Return the benchmark parameters of the given result
        """
        return self._benchmark_params[key]

    def remove_result(self, key):
        """
        Remove results corresponding to a given benchmark.
        """
        del self._results[key]
        del self._benchmark_params[key]
        del self._samples[key]
        del self._stats[key]

        # Remove profiles (may be missing)
        self._profiles.pop(key, None)

        # Remove run times (may be missing in old files)
        self._started_at.pop(key, None)
        self._ended_at.pop(key, None)

        # Remove version (may be missing)
        self._benchmark_version.pop(key, None)

    def remove_samples(self, key, selected_idx=None):
        """
        Remove measurement samples from the selected benchmark.
        """
        if key not in self._results:
            raise ValueError(key)

        if selected_idx is None:
            self._samples[key] = None
        elif self._samples[key] is not None:
            for j in selected_idx:
                self._samples[key][j] = None

    def add_result(self, benchmark, result,
                   started_at=None, ended_at=None,
                   record_samples=False,
                   append_samples=False,
                   selected_idx=None):
        """
        Add benchmark result.

        Parameters
        ----------
        benchmark : dict
            Benchmark object

        result : runner.BenchmarkResult
            Result of the benchmark.

        started_at : datetime.datetime, optional
            Benchmark start time.

        ended_at : datetime.datetime, optional
            Benchmark end time.

        record_samples : bool, optional
            Whether to save samples.

        append_samples : bool, optional
            Whether to combine new samples with old.

        selected_idx : set, optional
            Which indices in a parametrized benchmark to update

        """
        new_result = list(result.result)
        new_samples = list(result.samples)
        new_number = result.number

        benchmark_name = benchmark['name']
        benchmark_version = benchmark['version']

        if started_at is None:
            started_at = datetime.datetime.utcnow()

        if ended_at is None:
            ended_at = started_at

        new_stats = [None] * len(new_result)

        if (benchmark_name in self._results and
                benchmark_version == self._benchmark_version.get(benchmark_name)):

            # Append to old samples, if requested
            if append_samples:
                old_samples = self.get_result_samples(benchmark_name, benchmark['params'])
                for j in range(len(new_samples)):
                    if old_samples[j] is not None and new_samples[j] is not None:
                        new_samples[j] = old_samples[j] + new_samples[j]

            # Retain old result where requested
            merge_idx = [j for j in range(len(new_result))
                         if selected_idx is not None and j not in selected_idx]
            if merge_idx:
                old_result = self.get_result_value(benchmark_name, benchmark['params'])
                old_samples = self.get_result_samples(benchmark_name, benchmark['params'])
                old_stats = self.get_result_stats(benchmark_name, benchmark['params'])
                for j in merge_idx:
                    new_result[j] = old_result[j]
                    new_samples[j] = old_samples[j]
                    new_stats[j] = old_stats[j]

        # Recompute stats for updated entries (and drop unnecessary data)
        for j, (r, s, n) in enumerate(zip(new_result, new_samples, new_number)):
            if util.is_na(r):
                new_samples[j] = None
                new_stats[j] = None
                continue

            if n is not None:
                new_result[j], new_stats[j] = statistics.compute_stats(s, n)

        # Compress None lists to just None
        if all(x is None for x in new_result):
            new_result = None
        if all(x is None for x in new_samples):
            new_samples = None
        if all(x is None for x in new_stats):
            new_stats = None

        # Drop samples if requested
        if not record_samples:
            new_samples = None

        # Store result
        self._results[benchmark_name] = new_result
        self._stats[benchmark_name] = new_stats
        self._samples[benchmark_name] = new_samples

        self._benchmark_params[benchmark_name] = benchmark['params'] if benchmark['params'] else []
        self._started_at[benchmark_name] = util.datetime_to_js_timestamp(started_at)
        self._ended_at[benchmark_name] = util.datetime_to_js_timestamp(ended_at)
        self._benchmark_version[benchmark_name] = benchmark_version

        self._stderr[benchmark_name] = result.stderr
        self._errcode[benchmark_name] = result.errcode

        if result.profile:
            profile_data = base64.b64encode(zlib.compress(result.profile))
            if sys.version_info[0] >= 3:
                profile_data = profile_data.decode('ascii')
            self._profiles[benchmark_name] = profile_data

    def get_profile(self, benchmark_name):
        """
        Get the profile data for the given benchmark name.

        Parameters
        ----------
        benchmark_name : str
            Name of benchmark

        Returns
        -------
        profile_data : bytes
            Raw profile data

        """
        profile_data = self._profiles[benchmark_name]
        if sys.version_info[0] >= 3:
            profile_data = profile_data.encode('ascii')
        return zlib.decompress(base64.b64decode(profile_data))

    def has_profile(self, benchmark_name):
        """
        Does the given benchmark data have profiling information?
        """
        return benchmark_name in self._profiles

    def save(self, result_dir):
        """
        Save the results to disk, replacing existing results.

        Parameters
        ----------
        result_dir : str
            Path to root of results tree.
        """
        if self._filename is None:
            raise ValueError("Cannot save unnamed Results")

        path = os.path.join(result_dir, self._filename)

        results = {}
        for key in six.iterkeys(self._samples):
            # Save omitting default values
            value = {'result': self._results[key]}
            if self._samples[key] and any(x is not None for x in self._samples[key]):
                value['samples'] = self._samples[key]
            if self._stats[key] and any(x is not None for x in self._stats[key]):
                value['stats'] = self._stats[key]
            if self._benchmark_params[key]:
                value['params'] = self._benchmark_params[key]
            if list(value.keys()) == ['result']:
                value = value['result']
                if isinstance(value, list) and len(value) == 1:
                    value = value[0]
            results[key] = value

        data = {
            'results': results,
            'params': self._params,
            'requirements': self._requirements,
            'commit_hash': self._commit_hash,
            'date': self._date,
            'env_name': self._env_name,
            'python': self._python,
            'profiles': self._profiles,
            'started_at': self._started_at,
            'ended_at': self._ended_at,
            'benchmark_version': self._benchmark_version,
        }

        util.write_json(path, data, self.api_version, compact=True)

    def load_data(self, result_dir):
        """
        Load previous results for the current parameters (if any).
        """
        if self._filename is None:
            raise ValueError("Cannot load unnamed Results")

        path = os.path.join(result_dir, self._filename)

        if os.path.isfile(path):
            old = self.load(path)
            for dict_name in ('_results', '_samples', '_stats',
                              '_benchmark_params', '_profiles', '_started_at',
                              '_ended_at', '_benchmark_version'):
                setattr(self, dict_name, getattr(old, dict_name))

    @classmethod
    def load(cls, path, machine_name=None):
        """
        Load results from disk.

        Parameters
        ----------
        path : str
            Path to results file.
        machine_name : str, optional
            If given, check that the results file is for the given machine.

        """
        d = util.load_json(path, cls.api_version)

        try:
            obj = cls(
                d['params'],
                d['requirements'],
                d['commit_hash'],
                d['date'],
                d['python'],
                d.get('env_name',
                      environment.get_env_name('', d['python'], d['requirements']))
            )

            obj._results = {}
            obj._samples = {}
            obj._stats = {}
            obj._benchmark_params = {}

            for key, value in six.iteritems(d['results']):
                # Backward compatibility
                if not isinstance(value, dict):
                    value = {'result': [value], 'samples': None,
                             'stats': None, 'params': []}

                if not isinstance(value['result'], list):
                    value['result'] = [value['result']]

                if 'stats' in value and not isinstance(value['stats'], list):
                    value['stats'] = [value['stats']]

                value.setdefault('samples', None)
                value.setdefault('stats', None)
                value.setdefault('params', [])

                # Assign results
                obj._results[key] = value['result']
                obj._samples[key] = value['samples']
                obj._stats[key] = value['stats']
                obj._benchmark_params[key] = value['params']

            if 'profiles' in d:
                obj._profiles = d['profiles']
            obj._filename = os.path.join(*path.split(os.path.sep)[-2:])

            obj._started_at = d.get('started_at', {})
            obj._ended_at = d.get('ended_at', {})
            obj._benchmark_version = d.get('benchmark_version', {})
        except KeyError as exc:
            raise util.UserError(
                "Error loading results file '{0}': missing key {1}".format(
                    path, six.text_type(exc)))

        if machine_name is not None and obj.params.get('machine') != machine_name:
            raise util.UserError(
                "Error loading results file '{0}': machine name is not '{1}'".format(
                    path, machine_name))

        return obj

    def rm(self, result_dir):
        if self._filename is None:
            raise ValueError("Cannot remove unnamed Results")

        path = os.path.join(result_dir, self._filename)
        os.remove(path)

    @classmethod
    def update(cls, path):
        util.update_json(cls, path, cls.api_version)

    @property
    def env_name(self):
        return self._env_name


def format_benchmark_result(results, benchmark):
    """
    Pretty-print a benchmark result to human-readable form.

    Parameters
    ----------
    results : Results
        Result set object
    benchmark : dict
        Benchmark dictionary

    Returns
    -------
    info : {str, None}
        One-line description of results
    details : {str, None}
        Additional details

    """
    name = benchmark['name']

    result = results.get_result_value(name, benchmark['params'])
    stats = results.get_result_stats(name, benchmark['params'])

    total_count = len(result)
    failure_count = sum(r is None for r in result)

    info = None
    details = None

    # Display status
    if failure_count > 0:
        if failure_count == total_count:
            info = "failed"
        else:
            info = "{0}/{1} failed".format(failure_count, total_count)

    # Display results
    if benchmark['params']:
        # Long format display
        if failure_count == 0:
            info = "ok"

        display_result = [(v, statistics.get_err(v, s) if s is not None else None)
                          for v, s in zip(result, stats)]
        display = _format_benchmark_result(display_result, benchmark)
        display = "\n".join(display).strip()
        details = display
    else:
        if failure_count == 0:
            # Failure already shown above
            if not result:
                display = "[]"
            else:
                if stats[0]:
                    err = statistics.get_err(result[0], stats[0])
                else:
                    err = None
                display = util.human_value(result[0], benchmark['unit'], err=err)
                if len(result) > 1:
                    display += ";..."
            info = display

    return info, details


def _format_benchmark_result(result, benchmark, max_width=None):
    """
    Format the result from a parameterized benchmark as an ASCII table
    """
    if not result:
        return ['[]']

    def do_formatting(num_column_params):
        # Fold result to a table
        if num_column_params > 0:
            column_params = benchmark['params'][-num_column_params:]
        else:
            column_params = []

        rows = []
        if column_params:
            row_params = benchmark['params'][:-len(column_params)]
            header = benchmark['param_names'][:len(row_params)]
            column_param_permutations = list(itertools.product(*column_params))
            header += [" / ".join(_format_param_value(value) for value in values)
                       for values in column_param_permutations]
            rows.append(header)
            column_items = len(column_param_permutations)
            name_header = " / ".join(benchmark['param_names'][len(row_params):])
        else:
            column_items = 1
            row_params = benchmark['params']
            name_header = ""
            header = benchmark['param_names']
            rows.append(header)

        for j, values in enumerate(itertools.product(*row_params)):
            row_results = [util.human_value(x[0], benchmark['unit'], err=x[1])
                           for x in result[j*column_items:(j+1)*column_items]]
            row = [_format_param_value(value) for value in values] + row_results
            rows.append(row)

        if name_header:
            display = util.format_text_table(rows, 1,
                                             top_header_text=name_header,
                                             top_header_span_start=len(row_params))
        else:
            display = util.format_text_table(rows, 1)

        return display.splitlines()

    # Determine how many parameters can be fit to columns
    if max_width is None:
        max_width = util.get_terminal_width() * 3//4

    text = do_formatting(0)
    for j in range(1, len(benchmark['params'])):
        new_text = do_formatting(j)
        width = max(len(line) for line in new_text)
        if width < max_width:
            text = new_text
        else:
            break

    return text


def _format_param_value(value_repr):
    """
    Format a parameter value for displaying it as test output. The
    values are string obtained via Python repr.

    """
    regexs = ["^'(.+)'$",
              "^u'(.+)'$",
              "^<class '(.+)'>$"]

    for regex in regexs:
        m = re.match(regex, value_repr)
        if m and m.group(1).strip():
            return m.group(1)

    return value_repr
