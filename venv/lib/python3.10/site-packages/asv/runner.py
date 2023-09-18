# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import io
import json
import os
import sys
import re
import time
import tempfile
import itertools
import datetime
import pstats
import socket
import struct
import threading
import traceback

import six

from .console import log
from .results import Results, format_benchmark_result
from . import statistics
from . import util


WIN = (os.name == "nt")


# Can't use benchmark.__file__, because that points to the compiled
# file, so it can't be run by another version of Python.
BENCHMARK_RUN_SCRIPT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "benchmark.py")


JSON_ERROR_RETCODE = -257


BenchmarkResult = util.namedtuple_with_doc(
    'BenchmarkResult',
    ['result', 'samples', 'number', 'errcode', 'stderr', 'profile'],
    """
    Postprocessed benchmark result

    Attributes
    ----------
    result : list of object
        List of numeric values of the benchmarks (one for each parameter
        combination).
        Values are `None` if benchmark failed or NaN if it was skipped.
    samples : list of {list, None}
        List of lists of sampled raw data points (or Nones if
        no sampling done).
    number : list of {dict, None}
        List of actual repeat counts for each sample (or Nones if
        no sampling done).
    errcode : int
        Process exit code
    stderr : str
        Process stdout/stderr output
    profile : bytes
        If `profile` is `True` and run was at least partially successful,
        this key will be a byte string containing the cProfile data.
        Otherwise, None.
    """)


def skip_benchmarks(benchmarks, env, results=None):
    """
    Mark benchmarks as skipped.

    Parameters
    ----------
    benchmarks : Benchmarks
        Set of benchmarks to skip
    env : Environment
        Environment to skip them in
    results : Results, optional
        Where to store the results.
        If omitted, stored to a new unnamed Results object.

    Returns
    -------
    results : Results
        Benchmark results.

    """
    if results is None:
        results = Results.unnamed()

    log.warning("Skipping {0}".format(env.name))
    with log.indent():
        for name, benchmark in six.iteritems(benchmarks):
            log.step()
            log.warning('{0} skipped'.format(name))

            started_at = datetime.datetime.utcnow()
            r = fail_benchmark(benchmark)
            results.add_result(benchmark, r,
                               selected_idx=benchmarks.benchmark_selection.get(name),
                               started_at=started_at,
                               ended_at=datetime.datetime.utcnow())

    return results


def run_benchmarks(benchmarks, env, results=None,
                   show_stderr=False, quick=False, profile=False,
                   extra_params=None,
                   record_samples=False, append_samples=False,
                   run_rounds=None,
                   launch_method=None):
    """
    Run all of the benchmarks in the given `Environment`.

    Parameters
    ----------
    benchmarks : Benchmarks
        Benchmarks to run
    env : Environment object
        Environment in which to run the benchmarks.
    results : Results, optional
        Where to store the results.
        If omitted, stored to a new unnamed Results object.
    show_stderr : bool, optional
        When `True`, display any stderr emitted by the benchmark.
    quick : bool, optional
        When `True`, run each benchmark function exactly once.
        This is useful to quickly find errors in the benchmark
        functions, without taking the time necessary to get
        accurate timings.
    profile : bool, optional
        When `True`, run the benchmark through the `cProfile`
        profiler.
    extra_params : dict, optional
        Override values for benchmark attributes.
    record_samples : bool, optional
        Whether to retain result samples or discard them.
    append_samples : bool, optional
        Whether to retain any previously measured result samples
        and use them in statistics computations.
    run_rounds : sequence of int, optional
        Run rounds for benchmarks with multiple processes.
        If None, run all rounds.
    launch_method : {'auto', 'spawn', 'forkserver'}, optional
        Benchmark launching method to use.

    Returns
    -------
    results : Results
        Benchmark results.

    """

    if extra_params is None:
        extra_params = {}
    else:
        extra_params = dict(extra_params)

    if quick:
        extra_params['number'] = 1
        extra_params['repeat'] = 1
        extra_params['warmup_time'] = 0
        extra_params['processes'] = 1

    if results is None:
        results = Results.unnamed()

    # Find all setup_cache routines needed
    setup_cache_timeout = {}
    benchmark_order = {}
    cache_users = {}
    max_processes = 0

    def get_processes(benchmark):
        """Get number of processes to use for a job"""
        if 'processes' in extra_params:
            return int(extra_params['processes'])
        else:
            return int(benchmark.get('processes', 1))

    for name, benchmark in sorted(six.iteritems(benchmarks)):
        key = benchmark.get('setup_cache_key')
        setup_cache_timeout[key] = max(benchmark.get('setup_cache_timeout',
                                                     benchmark['timeout']),
                                       setup_cache_timeout.get(key, 0))
        benchmark_order.setdefault(key, []).append((name, benchmark))
        max_processes = max(max_processes, get_processes(benchmark))
        cache_users.setdefault(key, set()).add(name)

    if run_rounds is None:
        run_rounds = list(range(1, max_processes + 1))

    # Interleave benchmark runs, in setup_cache order
    existing_results = results.get_result_keys(benchmarks)

    def iter_run_items():
        for run_round in run_rounds[::-1]:
            for setup_cache_key, benchmark_set in six.iteritems(benchmark_order):
                for name, benchmark in benchmark_set:
                    log.step()

                    processes = get_processes(benchmark)

                    if run_round > processes:
                        if (not append_samples and
                                run_round == run_rounds[-1] and
                                name in existing_results):
                            # We need to remove samples here so that
                            # append_samples=False has an effect on all
                            # benchmarks regardless of whether they were
                            # run this round.
                            selected_idx = benchmarks.benchmark_selection.get(name)
                            results.remove_samples(name, selected_idx)
                        continue

                    is_final = (run_round == 1)
                    yield name, benchmark, setup_cache_key, is_final

    # Run benchmarks in order
    cache_dirs = {None: None}
    failed_benchmarks = set()
    failed_setup_cache = {}

    if append_samples:
        previous_result_keys = existing_results
    else:
        previous_result_keys = set()

    log.info("Benchmarking {0}".format(env.name))

    partial_info_time = None
    indent = log.indent()
    indent.__enter__()

    spawner = get_spawner(env, benchmarks.benchmark_dir,
                          launch_method=launch_method)

    try:
        # Preimport benchmark suite (if using forkserver)
        success, out = spawner.preimport()

        if success:
            if show_stderr and out:
                log.info("Importing benchmark suite produced output:")
                with log.indent():
                    log.error(out.rstrip())
        else:
            log.warning("Importing benchmark suite failed (skipping all benchmarks).")
            if show_stderr and out:
                with log.indent():
                    log.error(out)

            stderr = 'asv: benchmark suite import failed'
            for name, benchmark, setup_cache_key, is_final in iter_run_items():
                if name in failed_benchmarks:
                    continue

                selected_idx = benchmarks.benchmark_selection.get(name)
                started_at = datetime.datetime.utcnow()
                res = fail_benchmark(benchmark, stderr=stderr)
                results.add_result(benchmark, res,
                                   selected_idx=selected_idx,
                                   started_at=started_at,
                                   ended_at=datetime.datetime.utcnow(),
                                   record_samples=record_samples)
                failed_benchmarks.add(name)
            return results

        # Run benchmarks
        for name, benchmark, setup_cache_key, is_final in iter_run_items():
            selected_idx = benchmarks.benchmark_selection.get(name)

            started_at = datetime.datetime.utcnow()

            # Don't try to rerun failed benchmarks
            if name in failed_benchmarks:
                if is_final:
                    partial_info_time = None
                    log.info(name, reserve_space=True)
                    log_benchmark_result(results, benchmark,
                                         show_stderr=show_stderr)
                continue

            # Setup cache first, if needed
            if setup_cache_key is None:
                cache_dir = None
            elif setup_cache_key in cache_dirs:
                cache_dir = cache_dirs[setup_cache_key]
            elif setup_cache_key not in failed_setup_cache:
                partial_info_time = None
                short_key = os.path.relpath(setup_cache_key, benchmarks.benchmark_dir)
                log.info("Setting up {0}".format(short_key), reserve_space=True)
                cache_dir, stderr = spawner.create_setup_cache(
                    name, setup_cache_timeout[setup_cache_key])
                if cache_dir is not None:
                    log.add_padded('ok')
                    cache_dirs[setup_cache_key] = cache_dir
                else:
                    log.add_padded('failed')
                    if stderr and show_stderr:
                        with log.indent():
                            log.error(stderr)
                    failed_setup_cache[setup_cache_key] = stderr

            if setup_cache_key in failed_setup_cache:
                # Mark benchmark as failed
                partial_info_time = None
                log.warning('{0} skipped (setup_cache failed)'.format(name))
                stderr = 'asv: setup_cache failed\n\n{}'.format(failed_setup_cache[setup_cache_key])
                res = fail_benchmark(benchmark, stderr=stderr)
                results.add_result(benchmark, res,
                                   selected_idx=selected_idx,
                                   started_at=started_at,
                                   ended_at=datetime.datetime.utcnow(),
                                   record_samples=record_samples)
                failed_benchmarks.add(name)
                continue

            # If appending to previous results, make sure to use the
            # same value for 'number' attribute.
            cur_extra_params = extra_params
            if name in previous_result_keys:
                cur_extra_params = []
                prev_stats = results.get_result_stats(name, benchmark['params'])
                for s in prev_stats:
                    if s is None or 'number' not in s:
                        p = extra_params
                    else:
                        p = dict(extra_params)
                        p['number'] = s['number']
                    cur_extra_params.append(p)

            # Run benchmark
            if is_final:
                partial_info_time = None
                log.info(name, reserve_space=True)
            elif partial_info_time is None or time.time() > partial_info_time + 30:
                partial_info_time = time.time()
                log.info('Running ({0}--)'.format(name))

            res = run_benchmark(benchmark, spawner,
                                profile=profile,
                                selected_idx=selected_idx,
                                extra_params=cur_extra_params,
                                cwd=cache_dir)

            # Save result
            results.add_result(benchmark, res,
                               selected_idx=selected_idx,
                               started_at=started_at,
                               ended_at=datetime.datetime.utcnow(),
                               record_samples=(not is_final or record_samples),
                               append_samples=(name in previous_result_keys))

            previous_result_keys.add(name)

            if all(r is None for r in res.result):
                failed_benchmarks.add(name)

            # Log result
            if is_final:
                partial_info_time = None
                log_benchmark_result(results, benchmark,
                                     show_stderr=show_stderr)
            else:
                log.add('.')

            # Cleanup setup cache, if no users left
            if cache_dir is not None and is_final:
                cache_users[setup_cache_key].remove(name)
                if not cache_users[setup_cache_key]:
                    # No users of this cache left, perform cleanup
                    util.long_path_rmtree(cache_dir, True)
                    del cache_dirs[setup_cache_key]
    finally:
        # Cleanup any dangling caches
        for cache_dir in cache_dirs.values():
            if cache_dir is not None:
                util.long_path_rmtree(cache_dir, True)
        indent.__exit__(None, None, None)
        spawner.close()

    return results


def get_spawner(env, benchmark_dir, launch_method):
    has_fork = hasattr(os, 'fork') and hasattr(socket, 'AF_UNIX')

    if launch_method in (None, 'auto'):
        # Don't use ForkServer as default on OSX, because many Apple
        # things are not fork-safe
        if has_fork and sys.platform not in ('darwin',):
            launch_method = "forkserver"
        else:
            launch_method = "spawn"

    if launch_method == "spawn":
        spawner_cls = Spawner
    elif launch_method == "forkserver":
        if not has_fork:
            raise util.UserError("'forkserver' launch method not available "
                                 "on this platform")
        spawner_cls = ForkServer
    else:
        raise ValueError("Invalid launch_method: {}".format(launch_method))

    return spawner_cls(env, benchmark_dir)


def log_benchmark_result(results, benchmark, show_stderr=False):
    info, details = format_benchmark_result(results, benchmark)

    log.add_padded(info)
    if details:
        log.info(details, color='default')

    # Dump program output
    stderr = results.stderr.get(benchmark['name'])
    if stderr and show_stderr:
        with log.indent():
            log.error(stderr)


def fail_benchmark(benchmark, stderr='', errcode=1):
    """
    Return a BenchmarkResult describing a failed benchmark.
    """
    if benchmark['params']:
        # Mark only selected parameter combinations skipped
        params = itertools.product(*benchmark['params'])
        result = [None for idx in params]
        samples = [None] * len(result)
        number = [None] * len(result)
    else:
        result = [None]
        samples = [None]
        number = [None]

    return BenchmarkResult(result=result,
                           samples=samples,
                           number=number,
                           errcode=errcode,
                           stderr=stderr,
                           profile=None)


def run_benchmark(benchmark, spawner, profile,
                  selected_idx=None,
                  extra_params=None,
                  cwd=None,
                  prev_result=None):
    """
    Run a benchmark.

    Parameters
    ----------
    benchmark : dict
        Benchmark object dict
    spawner : Spawner
        Benchmark process spawner
    profile : bool
        Whether to run with profile
    selected_idx : set, optional
        Set of parameter indices to run for.
    extra_params : {dict, list}, optional
        Additional parameters to pass to the benchmark.
        If a list, each entry should correspond to a benchmark
        parameter combination.
    cwd : str, optional
        Working directory to run the benchmark in.
        If None, run in a temporary directory.

    Returns
    -------
    result : BenchmarkResult
        Result data.

    """

    if extra_params is None:
        extra_params = {}

    result = []
    samples = []
    number = []
    profiles = []
    stderr = ''
    errcode = 0

    if benchmark['params']:
        param_iter = enumerate(itertools.product(*benchmark['params']))
    else:
        param_iter = [(0, None)]

    for param_idx, params in param_iter:
        if selected_idx is not None and param_idx not in selected_idx:
            result.append(util.nan)
            samples.append(None)
            number.append(None)
            profiles.append(None)
            continue

        if isinstance(extra_params, list):
            cur_extra_params = extra_params[param_idx]
        else:
            cur_extra_params = extra_params

        res = _run_benchmark_single_param(
            benchmark, spawner, param_idx,
            extra_params=cur_extra_params, profile=profile,
            cwd=cwd)

        result += res.result
        samples += res.samples
        number += res.number

        profiles.append(res.profile)

        if res.stderr:
            stderr += "\n\n"
            stderr += res.stderr

        if res.errcode != 0:
            errcode = res.errcode

    return BenchmarkResult(
        result=result,
        samples=samples,
        number=number,
        errcode=errcode,
        stderr=stderr.strip(),
        profile=_combine_profile_data(profiles)
    )


def _run_benchmark_single_param(benchmark, spawner, param_idx,
                                profile, extra_params, cwd):
    """
    Run a benchmark, for single parameter combination index in case it
    is parameterized

    Parameters
    ----------
    benchmark : dict
        Benchmark object dict
    spawner : Spawner
        Benchmark process spawner
    param_idx : {int, None}
        Parameter index to run benchmark for
    profile : bool
        Whether to run with profile
    extra_params : dict
        Additional parameters to pass to the benchmark
    cwd : {str, None}
        Working directory to run the benchmark in.
        If None, run in a temporary directory.

    Returns
    -------
    result : BenchmarkResult
        Result data.

    """
    name = benchmark['name']
    if benchmark['params']:
        name += '-%d' % (param_idx,)

    if profile:
        profile_fd, profile_path = tempfile.mkstemp()
        os.close(profile_fd)
    else:
        profile_path = 'None'

    params_str = json.dumps(extra_params)

    if cwd is None:
        real_cwd = tempfile.mkdtemp()
    else:
        real_cwd = cwd

    result_file = tempfile.NamedTemporaryFile(delete=False)
    try:
        result_file.close()

        out, errcode = spawner.run(
            name=name, params_str=params_str, profile_path=profile_path,
            result_file_name=result_file.name,
            timeout=benchmark['timeout'],
            cwd=real_cwd)

        if errcode != 0:
            if errcode == util.TIMEOUT_RETCODE:
                out += "\n\nasv: benchmark timed out (timeout {0}s)\n".format(benchmark['timeout'])

            result = None
            samples = None
            number = None
        else:
            with open(result_file.name, 'r') as stream:
                data = stream.read()

            try:
                data = json.loads(data)
            except ValueError as exc:
                data = None
                errcode = JSON_ERROR_RETCODE
                out += "\n\nasv: failed to parse benchmark result: {0}\n".format(exc)

            # Special parsing for timing benchmark results
            if isinstance(data, dict) and 'samples' in data and 'number' in data:
                result = True
                samples = data['samples']
                number = data['number']
            else:
                result = data
                samples = None
                number = None

        if benchmark['params'] and out:
            params, = itertools.islice(itertools.product(*benchmark['params']),
                                       param_idx, param_idx + 1)
            out = "For parameters: {0}\n{1}".format(", ".join(params), out)

        if profile:
            with io.open(profile_path, 'rb') as profile_fd:
                profile_data = profile_fd.read()
            profile_data = profile_data if profile_data else None
        else:
            profile_data = None

        return BenchmarkResult(
            result=[result],
            samples=[samples],
            number=[number],
            errcode=errcode,
            stderr=out.strip(),
            profile=profile_data)

    except KeyboardInterrupt:
        spawner.interrupt()
        raise util.UserError("Interrupted.")
    finally:
        os.remove(result_file.name)
        if profile:
            os.remove(profile_path)
        if cwd is None:
            util.long_path_rmtree(real_cwd, True)


class Spawner(object):
    """
    Manage launching individual benchmark.py commands
    """

    def __init__(self, env, benchmark_dir):
        self.env = env
        self.benchmark_dir = os.path.abspath(benchmark_dir)
        self.interrupted = False

    def interrupt(self):
        self.interrupted = True

    def create_setup_cache(self, benchmark_id, timeout):
        cache_dir = tempfile.mkdtemp()

        out, _, errcode = self.env.run(
            [BENCHMARK_RUN_SCRIPT, 'setup_cache',
             os.path.abspath(self.benchmark_dir),
             benchmark_id],
            dots=False, display_error=False,
            return_stderr=True, valid_return_codes=None,
            redirect_stderr=True,
            cwd=cache_dir, timeout=timeout)

        if errcode == 0:
            return cache_dir, None
        else:
            util.long_path_rmtree(cache_dir, True)
            return None, out.strip()

    def run(self, name, params_str, profile_path, result_file_name, timeout, cwd):
        out, _, errcode = self.env.run(
            [BENCHMARK_RUN_SCRIPT, 'run', os.path.abspath(self.benchmark_dir),
             name, params_str, profile_path, result_file_name],
            dots=False, timeout=timeout,
            display_error=False, return_stderr=True, redirect_stderr=True,
            valid_return_codes=None, cwd=cwd)
        return out, errcode

    def preimport(self):
        return True, ""

    def close(self):
        pass


class ForkServer(Spawner):
    def __init__(self, env, root):
        super(ForkServer, self).__init__(env, root)

        if not (hasattr(os, 'fork') and hasattr(os, 'setpgid')):
            raise RuntimeError("ForkServer only available on POSIX")

        self.tmp_dir = tempfile.mkdtemp(prefix='asv-forkserver-')
        self.socket_name = os.path.join(self.tmp_dir, 'socket')

        self.server_proc = env.run(
            [BENCHMARK_RUN_SCRIPT, 'run_server', self.benchmark_dir, self.socket_name],
            return_popen=True, redirect_stderr=True)

        self._server_output = None
        self.stdout_reader_thread = threading.Thread(target=self._stdout_reader)
        self.stdout_reader_thread.start()

        # Wait for the socket to appear
        while self.stdout_reader_thread.is_alive():
            if os.path.exists(self.socket_name):
                break
            time.sleep(0.05)

        if not os.path.exists(self.socket_name):
            os.rmdir(self.tmp_dir)
            raise RuntimeError("Failed to start server thread")

    def _stdout_reader(self):
        try:
            out = self.server_proc.stdout.read()
            self.server_proc.stdout.close()
            out = out.decode('utf-8', 'replace')
        except Exception as exc:
            import traceback
            out = traceback.format_exc()

        self._server_output = out

    def run(self, name, params_str, profile_path, result_file_name, timeout, cwd):
        msg = {'action': 'run',
               'benchmark_id': name,
               'params_str': params_str,
               'profile_path': profile_path,
               'result_file': result_file_name,
               'timeout': timeout,
               'cwd': cwd}
        result = self._send_command(msg)
        return result['out'], result['errcode']

    def preimport(self):
        success = True
        out = ""
        try:
            out = self._send_command({'action': 'preimport'})
        except Exception as exc:
            success = False
            out = "asv: benchmark runner crashed\n"
            if isinstance(exc, util.UserError):
                out += str(exc)
            else:
                out += traceback.format_exc()
            out = out.rstrip()

        return success, out

    def _send_command(self, msg):
        msg = json.dumps(msg)
        if sys.version_info[0] >= 3:
            msg = msg.encode('utf-8')

        # Connect (with wait+retry)
        s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        for retry in range(5, 0, -1):
            try:
                s.connect(self.socket_name)
                break
            except socket.error:
                if retry > 1:
                    time.sleep(0.2)
                else:
                    raise

        # Send command
        try:
            s.sendall(struct.pack('<Q', len(msg)))
            s.sendall(msg)

            # Read result
            read_size, = struct.unpack('<Q', util.recvall(s, 8))
            result_text = util.recvall(s, read_size)
            if sys.version_info[0] >= 3:
                result_text = result_text.decode('utf-8')
            result = json.loads(result_text)
        except Exception:
            exitcode = self.server_proc.poll()
            if exitcode is not None:
                raise util.UserError("Process exited with code {0}".format(exitcode))
            raise
        finally:
            s.close()

        return result

    def close(self):
        import signal

        # Check for termination
        if self.server_proc.poll() is None:
            util._killpg_safe(self.server_proc.pid, signal.SIGINT)

        if self.server_proc.poll() is None:
            time.sleep(0.1)

        if self.server_proc.poll() is None:
            # Kill process group
            util._killpg_safe(self.server_proc.pid, signal.SIGKILL)

        self.server_proc.wait()
        self.stdout_reader_thread.join()

        if self._server_output and not self.interrupted:
            with log.indent():
                log.error("asv: forkserver:")
                log.error(self._server_output)

        util.long_path_rmtree(self.tmp_dir)


def _combine_profile_data(datasets):
    """
    Combine a list of profile data to a single profile
    """
    datasets = [data for data in datasets if data is not None]
    if not datasets:
        return None
    elif len(datasets) == 1:
        return datasets[0]

    # Load and combine stats
    stats = None

    while datasets:
        data = datasets.pop(0)

        f = tempfile.NamedTemporaryFile(delete=False)
        try:
            f.write(data)
            f.close()
            if stats is None:
                stats = pstats.Stats(f.name)
            else:
                stats.add(f.name)
        finally:
            os.remove(f.name)

    # Write combined stats out
    f = tempfile.NamedTemporaryFile(delete=False)
    try:
        f.close()
        stats.dump_stats(f.name)
        with open(f.name, 'rb') as fp:
            return fp.read()
    finally:
        os.remove(f.name)
