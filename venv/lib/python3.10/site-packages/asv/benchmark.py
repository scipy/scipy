# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""\
Usage: python -masv.benchmark COMMAND [...]

Manage a single benchmark and, when run from the commandline, report
its runtime to a file.

commands:

  timing [...]
      Run timing benchmark for given Python statement.

internal commands:

  discover BENCHMARK_DIR RESULT_FILE
      Discover benchmarks in a given directory and store result to a file.
  setup_cache BENCHMARK_DIR BENCHMARK_ID
      Run setup_cache for given benchmark.
  run BENCHMARK_DIR BENCHMARK_ID QUICK PROFILE_PATH RESULT_FILE
      Run a given benchmark, and store result in a file.
  run_server BENCHMARK_DIR SOCKET_FILENAME
      Run a Unix socket forkserver.
"""

# !!!!!!!!!!!!!!!!!!!! NOTE !!!!!!!!!!!!!!!!!!!!
# This file, unlike most others, must be compatible with as many
# versions of Python as possible and have no dependencies outside of
# the Python standard library.  This is the only bit of code from asv
# that is imported into the benchmarking process.

# Remove asv package directory from sys.path. This script file resides
# there although it's not part of the package, and Python puts it to
# sys.path[0] on start which can shadow other modules
import sys
if __name__ == "__main__":
    _old_sys_path_head = sys.path.pop(0)
else:
    _old_sys_path_head = None

import copy
try:
    import cProfile as profile
except:
    profile = None
import ctypes
from ctypes.util import find_library
from hashlib import sha256
import errno
if sys.version_info[0] >= 3:
    import importlib.machinery
else:
    import imp
import inspect
import itertools
import json
import os
import pickle
import re
import textwrap
import timeit
import time
import tempfile
import struct
import pkgutil
import traceback
import contextlib
from importlib import import_module


# The best timer we can use is time.process_time, but it is not
# available in the Python stdlib until Python 3.3.  This is a ctypes
# backport for Pythons that don't have it.

try:
    from time import process_time
except ImportError:  # Python <3.3
    if sys.platform.startswith("linux"):
        CLOCK_PROCESS_CPUTIME_ID = 2  # time.h

        clockid_t = ctypes.c_int
        time_t = ctypes.c_long

        class timespec(ctypes.Structure):
            _fields_ = [
                ('tv_sec', time_t),         # seconds
                ('tv_nsec', ctypes.c_long)  # nanoseconds
            ]
        _clock_gettime = ctypes.CDLL(
            find_library('rt'), use_errno=True).clock_gettime
        _clock_gettime.argtypes = [clockid_t, ctypes.POINTER(timespec)]

        def process_time():
            tp = timespec()
            if _clock_gettime(CLOCK_PROCESS_CPUTIME_ID, ctypes.byref(tp)) < 0:
                err = ctypes.get_errno()
                msg = errno.errorcode[err]
                if err == errno.EINVAL:
                    msg += (
                        "The clk_id (4) specified is not supported on this system")
                raise OSError(err, msg)
            return tp.tv_sec + tp.tv_nsec * 1e-9

    elif sys.platform == 'darwin':
        RUSAGE_SELF = 0  # sys/resources.h

        time_t = ctypes.c_long
        suseconds_t = ctypes.c_int32

        class timeval(ctypes.Structure):
            _fields_ = [
                ('tv_sec', time_t),
                ('tv_usec', suseconds_t)
            ]

        class rusage(ctypes.Structure):
            _fields_ = [
                ('ru_utime', timeval),
                ('ru_stime', timeval),
                ('ru_maxrss', ctypes.c_long),
                ('ru_ixrss', ctypes.c_long),
                ('ru_idrss', ctypes.c_long),
                ('ru_isrss', ctypes.c_long),
                ('ru_minflt', ctypes.c_long),
                ('ru_majflt', ctypes.c_long),
                ('ru_nswap', ctypes.c_long),
                ('ru_inblock', ctypes.c_long),
                ('ru_oublock', ctypes.c_long),
                ('ru_msgsnd', ctypes.c_long),
                ('ru_msgrcv', ctypes.c_long),
                ('ru_nsignals', ctypes.c_long),
                ('ru_nvcsw', ctypes.c_long),
                ('ru_nivcsw', ctypes.c_long)
            ]

        _getrusage = ctypes.CDLL(find_library('c'), use_errno=True).getrusage
        _getrusage.argtypes = [ctypes.c_int, ctypes.POINTER(rusage)]

        def process_time():
            ru = rusage()
            if _getrusage(RUSAGE_SELF, ctypes.byref(ru)) < 0:
                err = ctypes.get_errno()
                msg = errno.errorcode[err]
                if err == errno.EINVAL:
                    msg += (
                        "The clk_id (0) specified is not supported on this system")
                raise OSError(err, msg)
            return float(ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6 +
                         ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6)

    else:
        # Fallback to default timer
        process_time = timeit.default_timer

wall_timer = timeit.default_timer


def get_maxrss():
    # Fallback function, in case we don't have one that works on the
    # current platform
    return None

if sys.platform.startswith('win'):
    import ctypes
    import ctypes.wintypes

    SIZE_T = ctypes.c_size_t
    class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
        _fields_ = [
            ('cb', ctypes.wintypes.DWORD),
            ('PageFaultCount', ctypes.wintypes.DWORD),
            ('PeakWorkingSetSize', SIZE_T),
            ('WorkingSetSize', SIZE_T),
            ('QuotaPeakPagedPoolUsage', SIZE_T),
            ('QuotaPagedPoolUsage', SIZE_T),
            ('QuotaPeakNonPagedPoolUsage', SIZE_T),
            ('QuotaNonPagedPoolUsage', SIZE_T),
            ('PagefileUsage', SIZE_T),
            ('PeakPagefileUsage', SIZE_T),
        ]

    GetCurrentProcess = ctypes.windll.kernel32.GetCurrentProcess
    GetCurrentProcess.argtypes = []
    GetCurrentProcess.restype = ctypes.wintypes.HANDLE

    GetProcessMemoryInfo = ctypes.windll.psapi.GetProcessMemoryInfo
    GetProcessMemoryInfo.argtypes = (ctypes.wintypes.HANDLE,
                                     ctypes.POINTER(PROCESS_MEMORY_COUNTERS),
                                     ctypes.wintypes.DWORD)
    GetProcessMemoryInfo.restype = ctypes.wintypes.BOOL

    def get_maxrss():
        proc_hnd = GetCurrentProcess()
        counters = PROCESS_MEMORY_COUNTERS()
        info = GetProcessMemoryInfo(proc_hnd, ctypes.byref(counters), ctypes.sizeof(counters))
        if info == 0:
            raise ctypes.WinError()
        return counters.PeakWorkingSetSize

    if ctypes.sizeof(ctypes.c_void_p) == ctypes.sizeof(ctypes.c_uint64):
        DWORD_PTR = ctypes.c_uint64
    elif ctypes.sizeof(ctypes.c_void_p) == ctypes.sizeof(ctypes.c_uint32):
        DWORD_PTR = ctypes.c_uint32

    SetProcessAffinityMask = ctypes.windll.kernel32.SetProcessAffinityMask
    SetProcessAffinityMask.argtypes = [ctypes.wintypes.HANDLE, DWORD_PTR]
    SetProcessAffinityMask.restype  = bool

    GetCurrentProcess = ctypes.windll.kernel32.GetCurrentProcess
    GetCurrentProcess.argtypes = []
    GetCurrentProcess.restype  = ctypes.wintypes.HANDLE

    def set_cpu_affinity(affinity_list):
        """Set CPU affinity to CPUs listed (numbered 0...n-1)"""
        mask = 0
        for num in affinity_list:
            mask |= 2**num

        # Pseudohandle, doesn't need to be closed
        handle = GetCurrentProcess()
        ok = SetProcessAffinityMask(handle, mask)
        if not ok:
            raise RuntimeError("SetProcessAffinityMask failed")
else:
    try:
        import resource

        # POSIX
        if sys.platform == 'darwin':
            def get_maxrss():
                # OSX getrusage returns maxrss in bytes
                # https://developer.apple.com/library/mac/documentation/Darwin/Reference/ManPages/man2/getrusage.2.html
                return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        else:
            def get_maxrss():
                # Linux, *BSD return maxrss in kilobytes
                return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    except ImportError:
        pass

    def set_cpu_affinity(affinity_list):
        """Set CPU affinity to CPUs listed (numbered 0...n-1)"""
        if hasattr(os, 'sched_setaffinity'):
            os.sched_setaffinity(0, affinity_list)
        else:
            import psutil
            p = psutil.Process()
            if hasattr(p, 'cpu_affinity'):
                p.cpu_affinity(affinity_list)


def recvall(sock, size):
    """
    Receive data of given size from a socket connection
    """
    data = b""
    while len(data) < size:
        s = sock.recv(size - len(data))
        data += s
        if not s:
            raise RuntimeError("did not receive data from socket "
                               "(size {}, got only {!r})".format(size, data))
    return data


def _get_attr(source, name, ignore_case=False):
    if ignore_case:
        attrs = [getattr(source, key) for key in dir(source)
                 if key.lower() == name.lower()]

        if len(attrs) > 1:
            raise ValueError(
                "{0} contains multiple {1} functions.".format(
                    source.__name__, name))
        elif len(attrs) == 1:
            return attrs[0]
        else:
            return None
    else:
        return getattr(source, name, None)


def _get_all_attrs(sources, name, ignore_case=False):
    for source in sources:
        val = _get_attr(source, name, ignore_case=ignore_case)
        if val is not None:
            yield val


def _get_first_attr(sources, name, default, ignore_case=False):
    for val in _get_all_attrs(sources, name, ignore_case=ignore_case):
        return val
    return default


def get_setup_cache_key(func):
    if func is None:
        return None
    return '{0}:{1}'.format(inspect.getsourcefile(func),
                            inspect.getsourcelines(func)[1])


def get_source_code(items):
    """
    Extract source code of given items, and concatenate and dedent it.
    """
    sources = []
    prev_class_name = None

    for func in items:

        # custom source
        if hasattr(func, 'pretty_source'):
            src = textwrap.dedent(func.pretty_source).lstrip()
        # original source
        else:
            try:
                lines, _ = inspect.getsourcelines(func)
            except TypeError:
                continue

            if not lines:
                continue

            src = "\n".join(line.rstrip() for line in lines)
            src = textwrap.dedent(src)

        class_name = None
        if inspect.ismethod(func):
            # Add class name
            if hasattr(func, 'im_class'):
                class_name = func.im_class.__name__
            elif hasattr(func, '__qualname__'):
                names = func.__qualname__.split('.')
                if len(names) > 1:
                    class_name = names[-2]

        if class_name and prev_class_name != class_name:
            src = "class {0}:\n    {1}".format(
                class_name, src.replace("\n", "\n    "))
        elif class_name:
            src = "    {1}".format(
                class_name, src.replace("\n", "\n    "))

        sources.append(src)
        prev_class_name = class_name

    return "\n\n".join(sources).rstrip()


def _get_sourceline_info(obj, basedir):
    try:
        fn = inspect.getsourcefile(obj)
        fn = os.path.relpath(fn, basedir)
        _, lineno = inspect.getsourcelines(obj)
        return " in {!s}:{!s}".format(fn, lineno)
    except Exception:
        return ""


def check_num_args(root, benchmark_name, func, min_num_args, max_num_args=None):
    if max_num_args is None:
        max_num_args = min_num_args

    try:
        if sys.version_info[0] >= 3:
            info = inspect.getfullargspec(func)
        else:
            info = inspect.getargspec(func)
    except Exception as exc:
        print("{!s}: failed to check ({!r}{!s}): {!s}".format(
            benchmark_name, func, _get_sourceline_info(func, root), exc))
        return True

    max_args = len(info.args)

    if inspect.ismethod(func):
        max_args -= 1

    if info.defaults is not None:
        min_args = max_args - len(info.defaults)
    else:
        min_args = max_args

    if info.varargs is not None:
        max_args = float('inf')

    ok = (min_args <= max_num_args) and (min_num_args <= max_args)
    if not ok:
        if min_args == max_args:
            args_str = min_args
        else:
            args_str = "{}-{}".format(min_args, max_args)
        if min_num_args == max_num_args:
            num_args_str = min_num_args
        else:
            num_args_str = "{}-{}".format(min_num_args, max_num_args)
        print("{!s}: wrong number of arguments (for {!r}{!s}): expected {}, has {}".format(
            benchmark_name, func, _get_sourceline_info(func, root),
            num_args_str, args_str))

    return ok


class Benchmark(object):
    """
    Represents a single benchmark.
    """
    # The regex of the name of function or method to be considered as
    # this type of benchmark.  The default in the base class, will
    # match nothing.
    name_regex = re.compile('^$')

    def __init__(self, name, func, attr_sources):
        self.name = name
        self.func = func
        self.pretty_name = getattr(func, "pretty_name", None)
        self._attr_sources = attr_sources
        self._setups = list(_get_all_attrs(attr_sources, 'setup', True))[::-1]
        self._teardowns = list(_get_all_attrs(attr_sources, 'teardown', True))
        self._setup_cache = _get_first_attr(attr_sources, 'setup_cache', None)
        self.setup_cache_key = get_setup_cache_key(self._setup_cache)
        self.setup_cache_timeout = _get_first_attr([self._setup_cache], "timeout", None)
        self.timeout = _get_first_attr(attr_sources, "timeout", 60.0)
        self.code = get_source_code([self.func] + self._setups + [self._setup_cache])
        if sys.version_info[0] >= 3:
            code_text = self.code.encode('utf-8')
        else:
            code_text = self.code
        code_hash = sha256(code_text).hexdigest()
        self.version = str(_get_first_attr(attr_sources, "version", code_hash))
        self.type = "base"
        self.unit = "unit"

        self._redo_setup_next = False

        self._params = _get_first_attr(attr_sources, "params", [])
        self.param_names = _get_first_attr(attr_sources, "param_names", [])
        self._current_params = ()

        # Enforce params format
        try:
            self.param_names = [str(x) for x in list(self.param_names)]
        except ValueError:
            raise ValueError("%s.param_names is not a list of strings" % (name,))

        try:
            self._params = list(self._params)
        except ValueError:
            raise ValueError("%s.params is not a list" % (name,))

        if self._params and not isinstance(self._params[0], (tuple, list)):
            # Accept a single list for one parameter only
            self._params = [self._params]
        else:
            self._params = [[item for item in entry] for entry in self._params]

        if len(self.param_names) != len(self._params):
            self.param_names = self.param_names[:len(self._params)]
            self.param_names += ['param%d' % (k+1,) for k in range(len(self.param_names),
                                                                   len(self._params))]

        # Exported parameter representations
        self.params = [[repr(item) for item in entry] for entry in self._params]

    def set_param_idx(self, param_idx):
        try:
            self._current_params, = itertools.islice(
                itertools.product(*self._params),
                param_idx, param_idx + 1)
        except ValueError:
            raise ValueError(
                "Invalid benchmark parameter permutation index: %r" % (param_idx,))

    def insert_param(self, param):
        """
        Insert a parameter at the front of the parameter list.
        """
        self._current_params = tuple([param] + list(self._current_params))

    def __repr__(self):
        return '<{0} {1}>'.format(self.__class__.__name__, self.name)

    def check(self, root):
        # Check call syntax (number of arguments only...)
        ok = True

        if self._params:
            self.set_param_idx(0)

        min_num_args = len(self._current_params)
        max_num_args = min_num_args

        if self.setup_cache_key is not None:
            ok = ok and check_num_args(root, self.name + ": setup_cache",
                                       self._setup_cache, 0)
            max_num_args += 1

        for setup in self._setups:
            ok = ok and check_num_args(root, self.name + ": setup",
                                       setup, min_num_args, max_num_args)

        ok = ok and check_num_args(root, self.name + ": call",
                                   self.func, min_num_args, max_num_args)

        for teardown in self._teardowns:
            ok = ok and check_num_args(root, self.name + ": teardown",
                                       teardown, min_num_args, max_num_args)

        return ok

    def do_setup(self):
        try:
            for setup in self._setups:
                setup(*self._current_params)
        except NotImplementedError:
            # allow skipping test
            return True
        return False

    def redo_setup(self):
        if not self._redo_setup_next:
            self._redo_setup_next = True
            return
        self.do_teardown()
        self.do_setup()

    def do_teardown(self):
        for teardown in self._teardowns:
            teardown(*self._current_params)

    def do_setup_cache(self):
        if self._setup_cache is not None:
            return self._setup_cache()

    def do_run(self):
        return self.run(*self._current_params)

    def do_profile(self, filename=None):
        def method_caller():
            run(*params)

        if profile is None:
            raise RuntimeError("cProfile could not be imported")

        if filename is not None:
            if hasattr(method_caller, 'func_code'):
                code = method_caller.func_code
            else:
                code = method_caller.__code__

            self.redo_setup()

            profile.runctx(
                code, {'run': self.func, 'params': self._current_params},
                {}, filename)


class TimeBenchmark(Benchmark):
    """
    Represents a single benchmark for timing.
    """
    name_regex = re.compile(
        '^(Time[A-Z_].+)|(time_.+)$')

    def __init__(self, name, func, attr_sources):
        Benchmark.__init__(self, name, func, attr_sources)
        self.type = "time"
        self.unit = "seconds"
        self._attr_sources = attr_sources
        self.processes = int(_get_first_attr(self._attr_sources, 'processes', 2))
        self._load_vars()

    def _load_vars(self):
        self.repeat = _get_first_attr(self._attr_sources, 'repeat', 0)
        self.min_run_count = _get_first_attr(self._attr_sources, 'min_run_count', 2)
        self.number = int(_get_first_attr(self._attr_sources, 'number', 0))
        self.sample_time = _get_first_attr(self._attr_sources, 'sample_time', 0.01)
        self.warmup_time = _get_first_attr(self._attr_sources, 'warmup_time', -1)
        self.timer = _get_first_attr(self._attr_sources, 'timer', wall_timer)

    def do_setup(self):
        result = Benchmark.do_setup(self)
        # For parameterized tests, setup() is allowed to change these
        self._load_vars()
        return result

    def run(self, *param):
        warmup_time = self.warmup_time
        if warmup_time < 0:
            if '__pypy__' in sys.modules:
                warmup_time = 1.0
            else:
                # Transient effects exist also on CPython, e.g. from
                # OS scheduling
                warmup_time = 0.1

        if param:
            func = lambda: self.func(*param)
        else:
            func = self.func

        timer = timeit.Timer(
            stmt=func,
            setup=self.redo_setup,
            timer=self.timer)

        try:
            min_repeat, max_repeat, max_time = self.repeat
        except (ValueError, TypeError):
            if self.repeat == 0:
                min_repeat = 1
                max_repeat = 10
                max_time = 20.0
                if self.processes > 1:
                    max_repeat //= 2
                    max_time /= 2.0
            else:
                min_repeat = self.repeat
                max_repeat = self.repeat
                max_time = self.timeout

        min_repeat = int(min_repeat)
        max_repeat = int(max_repeat)
        max_time = float(max_time)

        samples, number = self.benchmark_timing(timer, min_repeat, max_repeat,
                                                max_time=max_time,
                                                warmup_time=warmup_time,
                                                number=self.number,
                                                min_run_count=self.min_run_count)

        samples = [s/number for s in samples]
        return {'samples': samples, 'number': number}

    def benchmark_timing(self, timer, min_repeat, max_repeat, max_time, warmup_time,
                         number, min_run_count):

        sample_time = self.sample_time
        start_time = wall_timer()
        run_count = 0

        samples = []

        def too_slow(num_samples):
            # stop taking samples if limits exceeded
            if run_count < min_run_count:
                return False
            if num_samples < min_repeat:
                return False
            return wall_timer() > start_time + warmup_time + max_time

        if number == 0:
            # Select number & warmup.
            #
            # This needs to be done at the same time, because the
            # benchmark timings at the beginning can be larger, and
            # lead to too small number being selected.
            number = 1
            while True:
                self._redo_setup_next = False
                start = wall_timer()
                timing = timer.timeit(number)
                wall_time = wall_timer() - start
                actual_timing = max(wall_time, timing)
                run_count += number

                if actual_timing >= sample_time:
                    if wall_timer() > start_time + warmup_time:
                        break
                else:
                    try:
                        p = min(10.0, max(1.1, sample_time/actual_timing))
                    except ZeroDivisionError:
                        p = 10.0
                    number = max(number + 1, int(p * number))

            if too_slow(1):
                return [timing], number
        elif warmup_time > 0:
            # Warmup
            while True:
                self._redo_setup_next = False
                timing = timer.timeit(number)
                run_count += number
                if wall_timer() >= start_time + warmup_time:
                    break

            if too_slow(1):
                return [timing], number

        # Collect samples
        while len(samples) < max_repeat:
            timing = timer.timeit(number)
            run_count += number
            samples.append(timing)

            if too_slow(len(samples)):
                break

        return samples, number


class MemBenchmark(Benchmark):
    """
    Represents a single benchmark for tracking the memory consumption
    of an object.
    """
    name_regex = re.compile(
        '^(Mem[A-Z_].+)|(mem_.+)$')

    def __init__(self, name, func, attr_sources):
        Benchmark.__init__(self, name, func, attr_sources)
        self.type = "memory"
        self.unit = "bytes"

    def run(self, *param):
        # We can't import asizeof directly, because we haven't loaded
        # the asv package in the benchmarking process.
        path = os.path.join(
            os.path.dirname(__file__), 'extern', 'asizeof.py')
        if sys.version_info[0] >= 3:
            asizeof = importlib.machinery.SourceFileLoader('asizeof', path).load_module()
        else:
            asizeof = imp.load_source('asizeof', path)

        obj = self.func(*param)

        sizeof2 = asizeof.asizeof([obj, obj])
        sizeofcopy = asizeof.asizeof([obj, copy.copy(obj)])

        return sizeofcopy - sizeof2


class PeakMemBenchmark(Benchmark):
    """
    Represents a single benchmark for tracking the peak memory consumption
    of the whole program.
    """
    name_regex = re.compile(
        '^(PeakMem[A-Z_].+)|(peakmem_.+)$')

    def __init__(self, name, func, attr_sources):
        Benchmark.__init__(self, name, func, attr_sources)
        self.type = "peakmemory"
        self.unit = "bytes"

    def run(self, *param):
        self.func(*param)
        return get_maxrss()


class TrackBenchmark(Benchmark):
    """
    Represents a single benchmark for tracking an arbitrary value.
    """
    name_regex = re.compile(
        '^(Track[A-Z_].+)|(track_.+)$')

    def __init__(self, name, func, attr_sources):
        Benchmark.__init__(self, name, func, attr_sources)
        self.type = _get_first_attr(attr_sources, "type", "track")
        self.unit = _get_first_attr(attr_sources, "unit", "unit")

    def run(self, *param):
        return self.func(*param)


# TODO: Support the creation of custom benchmark types


benchmark_types = [
    TimeBenchmark, MemBenchmark, PeakMemBenchmark, TrackBenchmark
]


class SpecificImporter(object):
    """
    Module importer that only allows loading a given module from the
    given path.

    Using this enables importing the asv benchmark suite without
    adding its parent directory to sys.path. The parent directory can
    in principle contain anything, including some version of the
    project module (common situation if asv.conf.json is on project
    repository top level).
    """

    def __init__(self, name, root):
        self._name = name
        self._root = root

    if sys.version_info[0] >= 3:
        def find_spec(self, fullname, path, target):
            if fullname == self._name:
                if path is not None:
                    raise ValueError()
                finder = importlib.machinery.PathFinder()
                return finder.find_spec(fullname, [self._root], target)
            return None
    else:
        def find_module(self, fullname, path=None):
            if fullname == self._name:
                return self
            return None

        def load_module(self, fullname):
            file, pathname, desc = imp.find_module(fullname, [self._root])
            return imp.load_module(fullname, file, pathname, desc)


def update_sys_path(root):
    sys.meta_path.insert(0, SpecificImporter(os.path.basename(root),
                                             os.path.dirname(root)))


def _get_benchmark(attr_name, module, klass, func):
    try:
        name = func.benchmark_name
    except AttributeError:
        name = None
        search = attr_name
    else:
        search = name.split('.')[-1]

    for cls in benchmark_types:
        if cls.name_regex.match(search):
            break
    else:
        return
    # relative to benchmark_dir
    mname_parts = module.__name__.split('.', 1)[1:]
    if klass is None:
        if name is None:
            name = ".".join(mname_parts + [func.__name__])
        sources = [func, module]
    else:
        instance = klass()
        func = getattr(instance, attr_name)
        if name is None:
            name = ".".join(mname_parts + [klass.__name__, attr_name])
        sources = [func, instance, module]
    return cls(name, func, sources)


def disc_modules(module_name, ignore_import_errors=False):
    """
    Recursively import a module and all sub-modules in the package

    Yields
    ------
    module
        Imported module in the package tree

    """
    if not ignore_import_errors:
        module = import_module(module_name)
    else:
        try:
            module = import_module(module_name)
        except BaseException:
            traceback.print_exc()
            return

    yield module

    if getattr(module, '__path__', None):
        for _, name, _ in pkgutil.iter_modules(module.__path__, module_name + '.'):
            for item in disc_modules(name, ignore_import_errors=ignore_import_errors):
                yield item


def disc_benchmarks(root, ignore_import_errors=False):
    """
    Discover all benchmarks in a given directory tree, yielding Benchmark
    objects

    For each class definition, looks for any methods with a
    special name.

    For each free function, yields all functions with a special
    name.
    """

    root_name = os.path.basename(root)

    for module in disc_modules(root_name, ignore_import_errors=ignore_import_errors):
        for attr_name, module_attr in (
            (k, v) for k, v in module.__dict__.items()
            if not k.startswith('_')
        ):
            if inspect.isclass(module_attr):
                for name, class_attr in inspect.getmembers(module_attr):
                    if (inspect.isfunction(class_attr) or
                            inspect.ismethod(class_attr)):
                        benchmark = _get_benchmark(name, module, module_attr,
                                                   class_attr)
                        if benchmark is not None:
                            yield benchmark
            elif inspect.isfunction(module_attr):
                benchmark = _get_benchmark(attr_name, module, None, module_attr)
                if benchmark is not None:
                    yield benchmark


def get_benchmark_from_name(root, name, extra_params=None):
    """
    Create a benchmark from a fully-qualified benchmark name.

    Parameters
    ----------
    root : str
        Path to the root of a benchmark suite.

    name : str
        Fully-qualified name to a specific benchmark.
    """

    if '-' in name:
        try:
            name, param_idx = name.split('-', 1)
            param_idx = int(param_idx)
        except ValueError:
            raise ValueError("Benchmark id %r is invalid" % (name,))
    else:
        param_idx = None

    update_sys_path(root)
    benchmark = None

    # try to directly import benchmark function by guessing its import module
    # name
    parts = name.split('.')
    for i in [1, 2]:
        path = os.path.join(root, *parts[:-i]) + '.py'
        if not os.path.isfile(path):
            continue
        modname = '.'.join([os.path.basename(root)] + parts[:-i])
        module = import_module(modname)
        try:
            module_attr = getattr(module, parts[-i])
        except AttributeError:
            break
        if i == 1 and inspect.isfunction(module_attr):
            benchmark = _get_benchmark(parts[-i], module, None, module_attr)
            break
        elif i == 2 and inspect.isclass(module_attr):
            try:
                class_attr = getattr(module_attr, parts[-1])
            except AttributeError:
                break
            if (inspect.isfunction(class_attr) or
                    inspect.ismethod(class_attr)):
                benchmark = _get_benchmark(parts[-1], module, module_attr,
                                           class_attr)
                break

    if benchmark is None:
        for benchmark in disc_benchmarks(root):
            if benchmark.name == name:
                break
        else:
            raise ValueError(
                "Could not find benchmark '{0}'".format(name))

    if param_idx is not None:
        benchmark.set_param_idx(param_idx)

    if extra_params:
        class ExtraBenchmarkAttrs:
            pass
        for key, value in extra_params.items():
            setattr(ExtraBenchmarkAttrs, key, value)
        benchmark._attr_sources.insert(0, ExtraBenchmarkAttrs)

    return benchmark


def list_benchmarks(root, fp):
    """
    List all of the discovered benchmarks to fp as JSON.
    """
    update_sys_path(root)

    # Streaming of JSON back out to the master process

    fp.write('[')
    first = True
    for benchmark in disc_benchmarks(root):
        if not first:
            fp.write(', ')
        clean = dict(
            (k, v) for (k, v) in benchmark.__dict__.items()
            if isinstance(v, (str, int, float, list, dict, bool)) and not
               k.startswith('_'))
        json.dump(clean, fp, skipkeys=True)
        first = False
    fp.write(']')


def main_discover(args):
    benchmark_dir, result_file = args
    with open(result_file, 'w') as fp:
        list_benchmarks(benchmark_dir, fp)


def main_check(args):
    benchmark_dir, = args

    update_sys_path(benchmark_dir)

    ok = True
    for benchmark in disc_benchmarks(benchmark_dir):
        ok = ok and benchmark.check(benchmark_dir)

    sys.exit(0 if ok else 1)


def main_setup_cache(args):
    (benchmark_dir, benchmark_id) = args
    benchmark = get_benchmark_from_name(benchmark_dir, benchmark_id)
    cache = benchmark.do_setup_cache()
    with open("cache.pickle", "wb") as fd:
        pickle.dump(cache, fd)


def main_run(args):
    (benchmark_dir, benchmark_id, params_str, profile_path, result_file) = args

    extra_params = json.loads(params_str)

    if profile_path == 'None':
        profile_path = None

    affinity_list = extra_params.pop('cpu_affinity', None)
    if affinity_list is not None:
        try:
            set_cpu_affinity(affinity_list)
        except BaseException as exc:
            print("asv: setting cpu affinity {!r} failed: {!r}".format(
                affinity_list, exc))

    benchmark = get_benchmark_from_name(
        benchmark_dir, benchmark_id, extra_params=extra_params)

    if benchmark.setup_cache_key is not None:
        with open("cache.pickle", "rb") as fd:
            cache = pickle.load(fd)
        if cache is not None:
            benchmark.insert_param(cache)

    skip = benchmark.do_setup()

    try:
        if skip:
            result = float('nan')
        else:
            result = benchmark.do_run()
            if profile_path is not None:
                benchmark.do_profile(profile_path)
    finally:
        benchmark.do_teardown()

    # Write the output value
    with open(result_file, 'w') as fp:
        json.dump(result, fp)


@contextlib.contextmanager
def posix_redirect_output(filename=None, permanent=True):
    """
    Redirect stdout/stderr to a file, using posix dup2.
    """
    sys.stdout.flush()
    sys.stderr.flush()

    stdout_fd = sys.stdout.fileno()
    stderr_fd = sys.stderr.fileno()

    if not permanent:
        stdout_fd_copy = os.dup(stdout_fd)
        stderr_fd_copy = os.dup(stderr_fd)

    if filename is None:
        out_fd, filename = tempfile.mkstemp()
    else:
        out_fd = os.open(filename, os.O_WRONLY | os.O_CREAT | os.O_TRUNC)

    try:
        # Redirect stdout and stderr to file
        os.dup2(out_fd, stdout_fd)
        os.dup2(out_fd, stderr_fd)

        yield filename
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.close(out_fd)

        if not permanent:
            os.dup2(stdout_fd_copy, stdout_fd)
            os.dup2(stderr_fd_copy, stderr_fd)
            os.close(stdout_fd_copy)
            os.close(stderr_fd_copy)


def main_run_server(args):
    import io
    import signal
    import socket

    benchmark_dir, socket_name, = args

    update_sys_path(benchmark_dir)

    # Socket I/O
    s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    s.bind(socket_name)
    s.listen(1)

    # Read and act on commands from socket
    while True:
        stdout_file = None

        try:
            conn, addr = s.accept()
        except KeyboardInterrupt:
            break

        try:
            fd, stdout_file = tempfile.mkstemp()
            os.close(fd)

            # Read command
            read_size, = struct.unpack('<Q', recvall(conn, 8))
            command_text = recvall(conn, read_size)
            if sys.version_info[0] >= 3:
                command_text = command_text.decode('utf-8')

            # Parse command
            command = json.loads(command_text)
            action = command.pop('action')

            if action == 'quit':
                break
            elif action == 'preimport':
                # Import benchmark suite before forking.
                # Capture I/O to a file during import.
                with posix_redirect_output(stdout_file, permanent=False):
                    for benchmark in disc_benchmarks(benchmark_dir, ignore_import_errors=True):
                        pass

                # Report result
                with io.open(stdout_file, 'r', errors='replace') as f:
                    out = f.read()
                out = json.dumps(out)
                if sys.version_info[0] >= 3:
                    out = out.encode('utf-8')
                conn.sendall(struct.pack('<Q', len(out)))
                conn.sendall(out)
                continue

            benchmark_id = command.pop('benchmark_id')
            params_str = command.pop('params_str')
            profile_path = command.pop('profile_path')
            result_file = command.pop('result_file')
            timeout = command.pop('timeout')
            cwd = command.pop('cwd')

            if command:
                raise RuntimeError('Command contained unknown data: {!r}'.format(command_text))

            # Spawn benchmark
            run_args = (benchmark_dir, benchmark_id, params_str, profile_path, result_file)
            pid = os.fork()
            if pid == 0:
                conn.close()
                sys.stdin.close()
                exitcode = 1
                try:
                    with posix_redirect_output(stdout_file, permanent=True):
                        try:
                            os.chdir(cwd)
                            main_run(run_args)
                            exitcode = 0
                        except BaseException as ec:
                            import traceback
                            traceback.print_exc()
                finally:
                    os._exit(exitcode)

            # Wait for results
            # (Poll in a loop is simplest --- also used by subprocess.py)
            start_time = wall_timer()
            is_timeout = False
            while True:
                res, status = os.waitpid(pid, os.WNOHANG)
                if res != 0:
                    break

                if timeout is not None and wall_timer() > start_time + timeout:
                    # Timeout
                    if is_timeout:
                        os.kill(pid, signal.SIGKILL)
                    else:
                        os.kill(pid, signal.SIGTERM)
                    is_timeout = True

                time.sleep(0.05)

            # Report result
            with io.open(stdout_file, 'r', errors='replace') as f:
                out = f.read()

            info = {'out': out,
                    'errcode': -256 if is_timeout else status}

            result_text = json.dumps(info)
            if sys.version_info[0] >= 3:
                result_text = result_text.encode('utf-8')

            conn.sendall(struct.pack('<Q', len(result_text)))
            conn.sendall(result_text)
        except KeyboardInterrupt:
            break
        finally:
            conn.close()
            if stdout_file is not None:
                os.unlink(stdout_file)


def main_timing(argv):
    import argparse
    import asv.statistics
    import asv.util
    import asv.console

    if (_old_sys_path_head is not None and
        os.path.abspath(_old_sys_path_head) != os.path.abspath(os.path.dirname(__file__))):
        sys.path.insert(0, _old_sys_path_head)

    parser = argparse.ArgumentParser(usage="python -masv.benchmark timing [options] STATEMENT")
    parser.add_argument("--setup", action="store", default=(lambda: None))
    parser.add_argument("--number", action="store", type=int, default=0)
    parser.add_argument("--repeat", action="store", type=int, default=0)
    parser.add_argument("--timer", action="store", choices=("process_time", "perf_counter"),
                        default="perf_counter")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("statement")
    args = parser.parse_args(argv)

    timer_func = {
        "process_time": process_time,
        "perf_counter": timeit.default_timer,
    }[args.timer]

    class AttrSource:
        pass

    attrs = AttrSource()
    attrs.repeat = args.repeat
    attrs.number = args.number
    attrs.timer = timer_func

    bench = TimeBenchmark("tmp", args.statement, [attrs])
    bench.redo_setup = args.setup
    result = bench.run()

    value, stats = asv.statistics.compute_stats(result['samples'], result['number'])
    formatted = asv.util.human_time(value, asv.statistics.get_err(value, stats))

    if not args.json:
        asv.console.color_print(formatted, 'red')
        asv.console.color_print("", 'default')
        asv.console.color_print("\n".join("{}: {}".format(k, v) for k, v in sorted(stats.items())), 'default')
        asv.console.color_print("samples: {}".format(result['samples']), 'default')
    else:
        json.dump({'result': value,
                   'samples': result['samples'],
                   'stats': stats}, sys.stdout)


def main_help(args):
    print(__doc__)


commands = {
    'discover': main_discover,
    'setup_cache': main_setup_cache,
    'run': main_run,
    'run_server': main_run_server,
    'check': main_check,
    'timing': main_timing,
    '-h': main_help,
    '--help': main_help,
}


def main():
    if len(sys.argv) < 2:
        main_help([])
        sys.exit(1)

    mode = sys.argv[1]
    args = sys.argv[2:]

    if mode in commands:
        commands[mode](args)
        sys.exit(0)
    else:
        sys.stderr.write("Unknown mode {0}\n".format(mode))
        sys.exit(1)

if __name__ == '__main__':
    main()
