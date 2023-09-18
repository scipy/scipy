# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Various low-level utilities.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import datetime
import json
import math
import os
import re
import select
import signal
import subprocess
import struct
import sys
import time
import errno
import threading
import shutil
import stat
import shlex
import operator
import collections

import six
from six.moves import xrange

from .extern import minify_json


nan = float('nan')
inf = float('inf')

WIN = (os.name == 'nt')

if not WIN:
    try:
        from select import PIPE_BUF
    except ImportError:
        # PIPE_BUF is not available on Python 2.6
        PIPE_BUF = os.pathconf('.', os.pathconf_names['PC_PIPE_BUF'])


TIMEOUT_RETCODE = -256


class UserError(Exception):
    pass


class ParallelFailure(Exception):
    """
    Custom exception to work around a multiprocessing bug
    https://bugs.python.org/issue9400
    """
    def __new__(cls, message, exc_cls, traceback_str):
        self = Exception.__new__(cls)
        self.message = message
        self.exc_cls = exc_cls
        self.traceback_str = traceback_str
        return self

    def __reduce__(self):
        return (ParallelFailure, (self.message, self.exc_cls, self.traceback_str))

    def __str__(self):
        return "{0}: {1}\n    {2}".format(self.exc_cls.__name__,
                                          self.message,
                                          self.traceback_str.replace("\n", "\n    "))

    def reraise(self):
        if self.exc_cls is UserError:
            raise UserError(self.message)
        else:
            raise self


def human_list(l):
    """
    Formats a list of strings in a human-friendly way.
    """
    l = ["'{0}'".format(x) for x in l]

    if len(l) == 0:
        return 'nothing'
    elif len(l) == 1:
        return l[0]
    elif len(l) == 2:
        return ' and '.join(l)
    else:
        return ', '.join(l[:-1]) + ' and ' + l[-1]


def human_float(value, significant=3, truncate_small=None, significant_zeros=False):
    """
    Return a string representing a float with human friendly significant digits.
    Switches to scientific notation for too large/small numbers.
    If `truncate_small`, then leading zeros of numbers < 1 are counted as 
    significant. If not `significant_zeros`, trailing unnecessary zeros are 
    stripped.
    """
    if value == 0:
        return "0"
    elif math.isinf(value) or math.isnan(value):
        return "{}".format(value)
    elif value < 0:
        sign = "-"
        value = -value
    else:
        sign = ""

    logv = math.log10(value)
    magnitude = int(math.floor(logv)) + 1

    if truncate_small is not None:
        magnitude = max(magnitude, -truncate_small + 1)

    num_digits = significant - magnitude

    if magnitude <= -5 or magnitude >= 9:
        # Too many digits, use scientific notation
        fmt = "{{0:.{0}e}}".format(significant)
    elif value == int(value):
        value = int(round(value, num_digits))
        fmt = "{0:d}"
    elif num_digits <= 0:
        value = int(round(value, num_digits))
        fmt = "{0:d}"
    else:
        fmt = "{{0:.{0}f}}".format(num_digits)

    formatted = sign + fmt.format(value)

    if not significant_zeros and '.' in formatted and 'e' not in fmt:
        formatted = formatted.rstrip('0')
        if formatted[-1] == '.':
            formatted = formatted[:-1]

    if significant_zeros and '.' not in formatted:
        if len(formatted) < significant:
            formatted += "." + "0"*(significant - len(formatted))

    return formatted


def human_file_size(size, err=None):
    """
    Returns a human-friendly string representing a file size
    that is 2-4 characters long.

    For example, depending on the number of bytes given, can be one
    of::

        256b
        64k
        1.1G

    Parameters
    ----------
    size : int
        The size of the file (in bytes)

    Returns
    -------
    size : str
        A human-friendly representation of the size of the file
    """
    size = float(size)

    if size < 1:
        size = 0.0

    suffixes = ' kMGTPEH'
    if size == 0:
        num_scale = 0
    else:
        num_scale = int(math.floor(math.log(size) / math.log(1000)))
    if num_scale > 7:
        suffix = '?'
    else:
        suffix = suffixes[num_scale].strip()
    scale = int(math.pow(1000, num_scale))
    value = size / scale

    str_value = human_float(value, 3)

    if err is None:
        return "{0:s}{1}".format(str_value, suffix)
    else:
        str_err = human_float(err / scale, 1, truncate_small=2)
        return "{0:s}±{1:s}{2}".format(str_value, str_err, suffix)

def human_time(seconds, err=None):
    """
    Returns a human-friendly time string that is always exactly 6
    characters long.

    Depending on the number of seconds given, can be one of::

        1w 3d
        2d 4h
        1h 5m
        1m 4s
          15s

    Will be in color if console coloring is turned on.

    Parameters
    ----------
    seconds : int
        The number of seconds to represent

    Returns
    -------
    time : str
        A human-friendly representation of the given number of seconds
        that is always exactly 6 characters.
    """
    units = [
        ('ns', 0.000000001),
        ('μs', 0.000001),
        ('ms', 0.001),
        ('s', 1),
        ('m', 60),
        ('h', 60 * 60),
        ('d', 60 * 60 * 24),
        ('w', 60 * 60 * 24 * 7),
        ('y', 60 * 60 * 24 * 7 * 52),
        ('C', 60 * 60 * 24 * 7 * 52 * 100)
    ]

    seconds = float(seconds)

    for i in xrange(len(units) - 1):
        if seconds < units[i+1][1]:
            str_time = human_float(seconds / units[i][1], 3, significant_zeros=True)
            if err is None:
                return "{0:s}{1}".format(str_time, units[i][0])
            else:
                str_err = human_float(err / units[i][1], 1, truncate_small=2)
                return "{0:s}±{1:s}{2}".format(str_time, str_err, units[i][0])
    return '~0'


def human_value(value, unit, err=None):
    """
    Formats a value in a given unit in a human friendly way.

    Parameters
    ----------
    value : anything
        The value to format

    unit : str
        The unit the value is in.  Currently understands `seconds` and `bytes`.

    err : float, optional
        Std. error in the value
    """
    if isinstance(value, (int, float)):
        if value != value:
            # nan
            display = "n/a"
        elif unit == 'seconds':
            display = human_time(value, err=err)
        elif unit == 'bytes':
            display = human_file_size(value, err=err)
        else:
            display = json.dumps(value)
            if err is not None:
                display += "±{:.2g}".format(err)
    elif value is None:
        display = "failed"
    else:
        display = json.dumps(value)

    return display


def which(filename, paths=None):
    """
    Emulates the UNIX `which` command in Python.

    Raises an IOError if no result is found.
    """
    # Hide traceback from expected exceptions in pytest reports
    __tracebackhide__ = operator.methodcaller('errisinstance', IOError)

    if os.path.sep in filename:
        locations = ['']
    elif paths is not None:
        locations = paths
    else:
        locations = os.environ.get("PATH", "").split(os.pathsep)
        if WIN:
            # On windows, an entry in %PATH% may be quoted
            locations = [path[1:-1] if len(path) > 2 and path[0] == path[-1] == '"' else path
                         for path in locations]

    if WIN:
        filenames = [filename + ext for ext in ('.exe', '.bat', '.com', '')]
    else:
        filenames = [filename]

    candidates = []
    for location in locations:
        for filename in filenames:
            candidate = os.path.join(location, filename)
            if os.path.isfile(candidate) or os.path.islink(candidate):
                candidates.append(candidate)

    if len(candidates) == 0:
        if paths is None:
            loc_info = 'PATH'
        else:
            loc_info = os.pathsep.join(locations)
        raise IOError("Could not find '{0}' in {1}".format(filename, loc_info))

    return candidates[0]


def has_command(filename):
    """
    Returns `True` if the commandline utility exists.
    """
    try:
        which(filename)
    except IOError:
        return False
    else:
        return True


class ProcessError(subprocess.CalledProcessError):
    def __init__(self, args, retcode, stdout, stderr):
        self.args = args
        self.retcode = retcode
        self.stdout = stdout
        self.stderr = stderr

    def __str__(self):
        if self.retcode == TIMEOUT_RETCODE:
            return "Command '{0}' timed out".format(
                ' '.join(self.args))
        else:
            return "Command '{0}' returned non-zero exit status {1}".format(
                ' '.join(self.args), self.retcode)


def check_call(args, valid_return_codes=(0,), timeout=600, dots=True,
               display_error=True, shell=False, env=None, cwd=None):
    """
    Runs the given command in a subprocess, raising ProcessError if it
    fails.

    See `check_output` for parameters.
    """
    # Hide traceback from expected exceptions in pytest reports
    __tracebackhide__ = operator.methodcaller('errisinstance', ProcessError)

    check_output(
        args, valid_return_codes=valid_return_codes, timeout=timeout,
        dots=dots, display_error=display_error, shell=shell, env=env,
        cwd=cwd)


class DebugLogBuffer(object):
    def __init__(self, log):
        self.buf = []
        self.first = True
        self.linebreak_re = re.compile(b'.*\n')
        self.log = log

    def __call__(self, c):
        if c is None:
            text = b"".join(self.buf)
            del self.buf[:]
        elif b'\n' in c:
            m = self.linebreak_re.match(c)
            j = m.end()
            self.buf.append(c[:j])
            text = b"".join(self.buf)
            self.buf[:] = [c[j:]]
        else:
            self.buf.append(c)
            return

        text = text.decode('utf-8', 'replace')
        if text.endswith('\n'):
            text = text[:-1]

        if text:
            if self.first:
                self.log.debug('OUTPUT -------->', continued=True)
                self.first = False
            self.log.debug(text, continued=True)


def check_output(args, valid_return_codes=(0,), timeout=600, dots=True,
                 display_error=True, shell=False, return_stderr=False,
                 env=None, cwd=None, redirect_stderr=False, return_popen=False):
    """
    Runs the given command in a subprocess, raising ProcessError if it
    fails.  Returns stdout as a string on success.

    Parameters
    ----------
    valid_return_codes : list, optional
        A list of return codes to ignore. Defaults to only ignoring zero.
        Setting to None ignores all return codes.

    timeout : number, optional
        Kill the process if it does not produce any output in `timeout`
        seconds. If `None`, there is no timeout.
        Default: 10 min

    dots : bool, optional
        If `True` (default) write a dot to the console to show
        progress as the subprocess outputs content.  May also be
        a callback function to call (with no arguments) to indicate
        progress.

    display_error : bool, optional
        If `True` (default) display the stdout and stderr of the
        subprocess when the subprocess returns an error code.

    shell : bool, optional
        If `True`, run the command through the shell.  Default is
        `False`.

    return_stderr : bool, optional
        If `True`, return both the (stdout, stderr, errcode) as a
        tuple.

    env : dict, optional
        Specify environment variables for the subprocess.

    cwd : str, optional
        Specify the current working directory to use when running the
        process.

    redirect_stderr : bool, optional
        Whether to redirect stderr to stdout. In this case the returned
        ``stderr`` (when return_stderr == True) is an empty string.

    return_popen : bool, optional
        Whether to return immediately after subprocess.Popen.

    Returns
    -------
    stdout, stderr, retcode : when return_stderr == True
    stdout : otherwise
    """
    from .console import log

    # Hide traceback from expected exceptions in pytest reports
    __tracebackhide__ = operator.methodcaller('errisinstance', ProcessError)

    def get_content(header=None):
        content = []
        if header is not None:
            content.append(header)
        if redirect_stderr:
            content.extend([
                'OUTPUT -------->',
                stdout[:-1]
            ])
        else:
            content.extend([
                'STDOUT -------->',
                stdout[:-1],
                'STDERR -------->',
                stderr[:-1]
            ])
        return '\n'.join(content)

    if isinstance(args, six.string_types):
        args = [args]

    log.debug("Running '{0}'".format(' '.join(args)))

    if env and WIN and sys.version_info < (3,):
        # Environment keys and values cannot be unicode
        def _fix_env(s):
            return s.encode('mbcs') if isinstance(s, unicode) else s
        env = {_fix_env(k): _fix_env(v) for k, v in env.items()}

    kwargs = dict(shell=shell, env=env, cwd=cwd,
                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if redirect_stderr:
        kwargs['stderr'] = subprocess.STDOUT
    if WIN:
        kwargs['close_fds'] = False
        kwargs['creationflags'] = subprocess.CREATE_NEW_PROCESS_GROUP
    else:
        kwargs['close_fds'] = True
        posix = getattr(os, 'setpgid', None)
        if posix:
            # Run the subprocess in a separate process group, so that we
            # can kill it and all child processes it spawns e.g. on
            # timeouts. Note that subprocess.Popen will wait until exec()
            # before returning in parent process, so there is no race
            # condition in setting the process group vs. calls to os.killpg
            kwargs['preexec_fn'] = lambda: os.setpgid(0, 0)

    proc = subprocess.Popen(args, **kwargs)

    if return_popen:
        return proc

    last_dot_time = time.time()
    stdout_chunks = []
    stderr_chunks = []
    is_timeout = False

    if log.is_debug_enabled():
        debug_log = DebugLogBuffer(log)
    else:
        debug_log = lambda c: None

    if WIN:
        start_time = [time.time()]
        was_timeout = [False]

        def stdout_reader_run():
            while True:
                c = proc.stdout.read(1)
                if not c:
                    break
                start_time[0] = time.time()
                stdout_chunks.append(c)
                debug_log(c)

        def stderr_reader_run():
            while True:
                c = proc.stderr.read(1)
                if not c:
                    break
                start_time[0] = time.time()
                stderr_chunks.append(c)
                debug_log(c)

        def watcher_run():
            while proc.returncode is None:
                time.sleep(0.1)
                if timeout is not None and time.time() - start_time[0] > timeout:
                    was_timeout[0] = True
                    proc.send_signal(signal.CTRL_BREAK_EVENT)

        watcher = threading.Thread(target=watcher_run)
        watcher.start()

        stdout_reader = threading.Thread(target=stdout_reader_run)
        stdout_reader.start()

        if not redirect_stderr:
            stderr_reader = threading.Thread(target=stderr_reader_run)
            stderr_reader.start()

        try:
            proc.wait()
        finally:
            if proc.returncode is None:
                proc.terminate()
                proc.wait()
            watcher.join()
            if not redirect_stderr:
                stderr_reader.join()
            stdout_reader.join()

            proc.stdout.close()
            if not redirect_stderr:
                proc.stderr.close()

        is_timeout = was_timeout[0]
    else:
        try:
            if posix and is_main_thread():
                # Forward signals related to Ctrl-Z handling; the child
                # process is in a separate process group so it won't receive
                # these automatically from the terminal
                def sig_forward(signum, frame):
                    _killpg_safe(proc.pid, signum)
                    if signum == signal.SIGTSTP:
                        os.kill(os.getpid(), signal.SIGSTOP)
                signal.signal(signal.SIGTSTP, sig_forward)
                signal.signal(signal.SIGCONT, sig_forward)

            fds = {
                proc.stdout.fileno(): stdout_chunks
                }
            if not redirect_stderr:
                fds[proc.stderr.fileno()] = stderr_chunks

            while proc.poll() is None:
                try:
                    if timeout is None:
                        rlist, wlist, xlist = select.select(
                            list(fds.keys()), [], [])
                    else:
                        rlist, wlist, xlist = select.select(
                            list(fds.keys()), [], [], timeout)
                except select.error as err:
                    if err.args[0] == errno.EINTR:
                        # interrupted by signal handler; try again
                        continue
                    raise

                if len(rlist) == 0:
                    # We got a timeout
                    is_timeout = True
                    break
                for f in rlist:
                    output = os.read(f, PIPE_BUF)
                    fds[f].append(output)
                    debug_log(output)
                if dots and time.time() - last_dot_time > 0.5:
                    if dots is True:
                        log.dot()
                    elif dots:
                        dots()
                    last_dot_time = time.time()
        finally:
            if posix and is_main_thread():
                # Restore signal handlers
                signal.signal(signal.SIGTSTP, signal.SIG_DFL)
                signal.signal(signal.SIGCONT, signal.SIG_DFL)

            if proc.returncode is None:
                # Timeout or another exceptional condition occurred, and
                # the program is still running.
                if posix:
                    # Terminate the whole process group
                    _killpg_safe(proc.pid, signal.SIGTERM)

                    for j in range(10):
                        time.sleep(0.1)
                        if proc.poll() is not None:
                            break
                    else:
                        # Didn't terminate within 1 sec, so kill it
                        _killpg_safe(proc.pid, signal.SIGKILL)
                else:
                    proc.terminate()
                proc.wait()

        proc.stdout.flush()
        if not redirect_stderr:
            proc.stderr.flush()

        stdout_chunks.append(proc.stdout.read())
        if not redirect_stderr:
            stderr_chunks.append(proc.stderr.read())

        proc.stdout.close()
        if not redirect_stderr:
            proc.stderr.close()

    debug_log(None)

    stdout = b''.join(stdout_chunks)
    stderr = b''.join(stderr_chunks)

    stdout = stdout.decode('utf-8', 'replace')
    stderr = stderr.decode('utf-8', 'replace')

    if is_timeout:
        retcode = TIMEOUT_RETCODE
    else:
        retcode = proc.returncode

    if valid_return_codes is not None and retcode not in valid_return_codes:
        header = 'Error running {0} (exit status {1})'.format(' '.join(args), retcode)
        if display_error:
            if log.is_debug_enabled():
                # Output was already printed
                log.error(header)
            else:
                log.error(get_content(header))
        raise ProcessError(args, retcode, stdout, stderr)

    if return_stderr:
        return (stdout, stderr, retcode)
    else:
        return stdout


def _killpg_safe(pgid, signo):
    """
    Same as os.killpg, but deal with OSX/BSD
    """
    try:
        os.killpg(pgid, signo)
    except OSError as exc:
        if exc.errno == errno.EPERM:
            # OSX/BSD may raise EPERM on killpg if the process group
            # already terminated
            pass
        else:
            raise


def is_main_thread():
    """
    Return True if the current thread is the main thread.
    """
    if sys.version_info[0] >= 3:
        return threading.current_thread() == threading.main_thread()
    else:
        return isinstance(threading.current_thread(), threading._MainThread)


def write_json(path, data, api_version=None, compact=False):
    """
    Writes JSON to the given path, including indentation and sorting.

    Parameters
    ----------
    path : str
        File name to write
    data : object
        Data to serialize as JSON
    api_version : int, optional
        API version number
    compact : bool, optional
        Whether to produce compact, non-human readable JSON.
        Disables sorting and indentation.
    """
    path = os.path.abspath(path)

    dirname = long_path(os.path.dirname(path))
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    if api_version is not None:
        data = dict(data)
        data['version'] = api_version

    open_kwargs = {}
    if sys.version_info[0] >= 3:
        open_kwargs['encoding'] = 'utf-8'

    with long_path_open(path, 'w', **open_kwargs) as fd:
        if not compact:
            json.dump(data, fd, indent=4, sort_keys=True)
        else:
            json.dump(data, fd)


def load_json(path, api_version=None, js_comments=False):
    """
    Loads JSON from the given path.

    Parameters
    ----------
    path : str
        File name
    api_version : str or None
        API version indentifier
    js_comments : bool, optional
        Whether to allow nonstandard javascript-style comments
        in the file. Note that this slows down the loading
        significantly.
    """
    # Hide traceback from expected exceptions in pytest reports
    __tracebackhide__ = operator.methodcaller('errisinstance', UserError)

    path = os.path.abspath(path)

    open_kwargs = {}
    if sys.version_info[0] >= 3:
        open_kwargs['encoding'] = 'utf-8'

    with long_path_open(path, 'r', **open_kwargs) as fd:
        content = fd.read()

    if js_comments:
        content = minify_json.json_minify(content)
        content = content.replace(",]", "]")
        content = content.replace(",}", "}")

    try:
        d = json.loads(content)
    except ValueError as e:
        raise UserError(
            "Error parsing JSON in file '{0}': {1}".format(
                path, six.text_type(e)))

    if api_version is not None:
        if 'version' in d:
            if d['version'] < api_version:
                raise UserError(
                    "{0} is stored in an old file format.  Run "
                    "`asv update` to update it.".format(path))
            elif d['version'] > api_version:
                raise UserError(
                    "{0} is stored in a format that is newer than "
                    "what this version of asv understands.  Update "
                    "asv to use this file.".format(path))

            del d['version']
        else:
            raise UserError(
                "No version specified in {0}.".format(path))

    return d


def update_json(cls, path, api_version):
    """
    Perform JSON file format updates.

    Parameters
    ----------
    cls : object
        Object containing methods update_to_X which updates
        the given JSON tree from version X-1 to X.

    path : str
        Path to JSON file

    api_version : int
        The current API version
    """
    # Hide traceback from expected exceptions in pytest reports
    __tracebackhide__ = operator.methodcaller('errisinstance', UserError)

    d = load_json(path)
    if 'version' not in d:
        raise UserError(
            "No version specified in {0}.".format(path))

    if d['version'] < api_version:
        for x in six.moves.xrange(d['version'] + 1, api_version + 1):
            d = getattr(cls, 'update_to_{0}'.format(x), lambda x: x)(d)
        write_json(path, d, api_version)
    elif d['version'] > api_version:
        raise UserError(
            "{0} is stored in a format that is newer than "
            "what this version of asv understands. "
            "Upgrade asv in order to use or add to "
            "these results.".format(path))


def iter_chunks(s, n):
    """
    Iterator that returns elements from s in chunks of size n.
    """
    chunk = []
    for x in s:
        chunk.append(x)
        if len(chunk) == n:
            yield chunk
            chunk = []
    if len(chunk):
        yield chunk


def pick_n(items, n):
    """Pick n items, attempting to get equal index spacing.
    """
    if not (n > 0):
        raise ValueError("Invalid number of items to pick")
    spacing = max(float(len(items)) / n, 1)
    spaced = []
    i = 0
    while int(i) < len(items) and len(spaced) < n:
        spaced.append(items[int(i)])
        i += spacing
    return spaced


def get_multiprocessing(parallel):
    """
    If parallel indicates that we want to do multiprocessing,
    imports the multiprocessing module and sets the parallel
    value accordingly.
    """
    if parallel != 1:
        import multiprocessing
        if parallel <= 0:
            parallel = multiprocessing.cpu_count()
        return parallel, multiprocessing
    return parallel, None


def iter_subclasses(cls):
    """
    Returns all subclasses of a class.
    """
    for x in cls.__subclasses__():
        yield x
        for y in iter_subclasses(x):
            yield y


def hash_equal(a, b):
    """
    Returns `True` if a and b represent the same commit hash.
    """
    min_len = min(len(a), len(b))
    return a.lower()[:min_len] == b.lower()[:min_len]


def get_cpu_info():
    """
    Gets a human-friendly description of this machine's CPU.

    Returns '' if it can't be obtained.
    """
    if sys.platform.startswith('linux'):
        with open("/proc/cpuinfo", "rb") as fd:
            lines = fd.readlines()
        for line in lines:
            if b':' in line:
                key, val = line.split(b':', 1)
                key = key.strip()
                val = val.strip()
                if key == b'model name':
                    return val.decode('ascii')
    elif sys.platform.startswith('darwin'):
        sysctl = which('sysctl')
        return check_output([sysctl, '-n', 'machdep.cpu.brand_string']).strip()
    return ''


def get_memsize():
    """
    Returns the amount of physical memory in this machine.

    Returns '' if it can't be obtained.
    """
    if sys.platform.startswith('linux'):
        with open("/proc/meminfo", "rb") as fd:
            lines = fd.readlines()
        for line in lines:
            if b':' in line:
                key, val = line.split(b':', 1)
                key = key.strip()
                val = val.strip()
                if key == b'MemTotal':
                    return int(val.split()[0])
    elif sys.platform.startswith('darwin'):
        sysctl = which('sysctl')
        return int(check_output([sysctl, '-n', 'hw.memsize']).strip())
    return ''


def _get_terminal_size_fallback():
    """
    Returns a tuple (height, width) containing the height and width of
    the terminal.  Fallback for when sys.get_terminal_size() doesn't
    exist or fails.
    """
    try:
        # Unix-specific code
        import fcntl
        import termios
        s = struct.pack(str("HHHH"), 0, 0, 0, 0)
        x = fcntl.ioctl(sys.stdout, termios.TIOCGWINSZ, s)
        (lines, width, xpixels, ypixels) = struct.unpack(str("HHHH"), x)
        if lines > 12:
            lines -= 6
        if width > 10:
            width -= 1
        return (lines, width)
    except:
        # Fall back on environment variables, or if not set, (25, 80)
        try:
            return (int(os.environ.get('LINES')),
                    int(os.environ.get('COLUMNS')))
        except TypeError:
            return 25, 80


def get_terminal_width():
    """
    Return the terminal width, or an estimate thereof.
    """
    try:
        # Python 3.3 and higher: this works under Windows and Unix
        return os.get_terminal_size().columns
    except (AttributeError, OSError):
        return _get_terminal_size_fallback()[1]


def format_text_table(rows, num_headers=0,
                      top_header_span_start=0,
                      top_header_text=None):
    """
    Format rows in as a reStructuredText table, in the vein of::

       ========== ========== ==========
       --         top header text, span start 1
       ---------- ---------------------
        row0col0     r0c1      r0c2
       ========== ========== ==========
        row1col0     r1c1      r1c2
        row2col0     r2c1      r2c2
       ========== ========== ==========

    """

    # Format content
    text_rows = [["{0}".format(item).replace("\n", " ") for item in row]
                 for row in rows]

    # Ensure same number of items on all rows
    num_items = max(len(row) for row in text_rows)
    for row in text_rows:
        row.extend(['']*(num_items - len(row)))

    # Determine widths
    col_widths = [max(len(row[j]) for row in text_rows) + 2
                  for j in range(num_items)]

    # Pad content
    text_rows = [[item.center(w) for w, item in zip(col_widths, row)]
                 for row in text_rows]

    # Generate result
    headers = [" ".join(row) for row in text_rows[:num_headers]]
    content = [" ".join(row) for row in text_rows[num_headers:]]
    separator = " ".join("-"*w for w in col_widths)

    result = []
    if top_header_text is not None:
        left_span = "-".join("-"*w for w in col_widths[:top_header_span_start])
        right_span = "-".join("-"*w for w in col_widths[top_header_span_start:])
        if left_span and right_span:
            result += ["--" + " " * (len(left_span)-1) + top_header_text.center(len(right_span))]
            result += [" ".join([left_span, right_span])]
        else:
            result += [top_header_text.center(len(separator))]
            result += ["-".join([left_span, right_span])]
        result += headers
        result += [separator.replace("-", "=")]
    elif headers:
        result += headers
        result += [separator]
    result += content
    result = [separator.replace("-", "=")] + result
    result += [separator.replace("-", "=")]
    return "\n".join(result)


def _datetime_to_timestamp(dt, divisor):
    delta = dt - datetime.datetime(1970, 1, 1)
    microseconds = (delta.days * 86400 + delta.seconds) * 10**6 + delta.microseconds
    value, remainder = divmod(microseconds, divisor)
    if remainder >= divisor//2:
        value += 1
    return value


def datetime_to_timestamp(dt):
    """
    Convert a Python datetime object to a UNIX timestamp.
    """
    return _datetime_to_timestamp(dt, 10**6)


def datetime_to_js_timestamp(dt):
    """
    Convert a Python datetime object to a JavaScript timestamp.
    """
    return _datetime_to_timestamp(dt, 10**3)


def js_timestamp_to_datetime(ts):
    """
    Convert a JavaScript timestamp to a Python datetime object.
    """
    return datetime.datetime.fromtimestamp(ts / 1000)


def is_nan(x):
    """
    Returns `True` if x is a NaN value.
    """
    if isinstance(x, float):
        return x != x
    return False


def is_na(value):
    """
    Return True if value is None or NaN
    """
    return value is None or is_nan(value)


def mean_na(values):
    """
    Take a mean, with the understanding that None and NaN stand for
    missing data.
    """
    values = [x for x in values if not is_na(x)]
    if values:
        return sum(values) / len(values)
    else:
        return None


def geom_mean_na(values):
    """
    Compute geometric mean, with the understanding that None and NaN
    stand for missing data.
    """
    values = [x for x in values if not is_na(x)]
    if values:
        exponent = 1/len(values)
        prod = 1.0
        acc = 0
        for x in values:
            prod *= abs(x)**exponent
            acc += x
        return prod if acc >= 0 else -prod
    else:
        return None


def ceildiv(numerator, denominator):
    """Ceiling division"""
    return -((-numerator)//denominator)


if not WIN:
    long_path_open = open
    long_path_rmtree = shutil.rmtree
    def long_path(path):
        return path
else:
    def long_path(path):
        if path.startswith("\\\\"):
            return path
        return "\\\\?\\" + os.path.abspath(path)

    def _remove_readonly(func, path, exc_info):
        """Try harder to remove files on Windows"""

        if isinstance(exc_info[1], OSError) and exc_info[1].errno == errno.EACCES:
            # Clear read-only flag and try again
            try:
                os.chmod(path, stat.S_IWRITE | stat.S_IREAD)
                func(path)
                return
            except OSError:
                pass

        # Reraise original error
        six.reraise(*exc_info)

    def long_path_open(filename, *a, **kw):
        return open(long_path(filename), *a, **kw)

    def long_path_rmtree(path, ignore_errors=False):
        if ignore_errors:
            onerror = None
        else:
            onerror = _remove_readonly
        shutil.rmtree(long_path(path),
                      ignore_errors=ignore_errors,
                      onerror=onerror)


def sanitize_filename(filename):
    """
    Replace characters to make a string safe to use in file names.

    This is not a 1-to-1 mapping.

    The implementation needs to match www/asv.js:escape_graph_parameter
    """
    if not isinstance(filename, six.text_type):
        filename = filename.decode(sys.getfilesystemencoding())

    # ntfs & ext3
    filename = re.sub('[<>:"/\\^|?*\x00-\x1f]', '_', filename)

    # ntfs
    forbidden = ["CON", "PRN", "AUX", "NUL", "COM1", "COM2", "COM3",
                 "COM4", "COM5", "COM6", "COM7", "COM8", "COM9", "LPT1",
                 "LPT2", "LPT3", "LPT4", "LPT5", "LPT6", "LPT7", "LPT8",
                 "LPT9"]
    if filename.upper() in forbidden:
        filename = filename + "_"

    return filename


def namedtuple_with_doc(name, slots, doc):
    cls = collections.namedtuple(name, slots)
    if sys.version_info[0] >= 3:
        cls.__doc__ = doc
        return cls
    else:
        return type(str(name), (cls,), {'__doc__': doc})


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


def interpolate_command(command, variables):
    """
    Parse a command with interpolated variables to a sequence of commands.

    The command is parsed as in posix-style shell (by shlex) and split to
    parts. Additional constructs recognized:

    - ``ENVVAR=value <command>``: parsed as declaring an environment variable
      named 'ENVVAR'.
    - ``return-code=value <command>``: parsed as declaring valid return codes.
    - ``in-dir=value <command>``: parsed as declaring working directory for command.

    Parameters
    ----------
    command : str
        Command to execute, posix shell style.
    variables : dict
        Interpolation variables.

    Returns
    -------
    command : list of str
        Command arguments.
    env : dict
        Environment variables declared in the command.
    return_codes : {set, int, None}
        Valid return codes.
    cwd : {str, None}
        Current working directory for the command, if any.

    """

    parts = shlex.split(command)

    try:
        result = [c.format(**variables) for c in parts]
    except KeyError as exc:
        raise UserError("Configuration error: {{{0}}} not available "
                        "when substituting into command {1!r} "
                        "Available: {2!r}"
                        "".format(exc.args[0], command, variables))

    env = {}

    return_codes_set = False
    return_codes = {0}
    cwd = None

    while result:
        m = re.match('^([A-Za-z_][A-Za-z0-9_]*)=(.*)$', result[0])
        if m:
            env[m.group(1)] = m.group(2)
            del result[0]
            continue

        if result[0].startswith('return-code='):
            if return_codes_set:
                raise UserError("Configuration error: multiple return-code specifications "
                                "in command {0!r} "
                                "".format(command))
                break

            if result[0] == 'return-code=any':
                return_codes = None
                return_codes_set = True
                del result[0]
                continue

            m = re.match('^return-code=([0-9,]+)$', result[0])
            if m:
                try:
                    return_codes = set(int(x) for x in m.group(1).split(","))
                    return_codes_set = True
                    del result[0]
                    continue
                except ValueError as exc:
                    pass

            raise UserError("Configuration error: invalid return-code specification "
                            "{0!r} when substituting into command {1!r} "
                            "".format(result[0], command))

        if result[0].startswith('in-dir='):
            if cwd is not None:
                raise UserError("Configuration error: multiple in-dir specifications "
                                "in command {0!r} "
                                "".format(command))
                break

            cwd = result[0][7:]
            del result[0]
            continue

        break

    return result, env, return_codes, cwd


try:
    from shlex import quote as shlex_quote
except ImportError:
    _find_unsafe = re.compile(r'[^\w@%+=:,./-]').search

    def shlex_quote(s):
        """Return a shell-escaped version of the string *s*."""
        if not s:
            return "''"
        if _find_unsafe(s) is None:
            return s

        # use single quotes, and put single quotes into double quotes
        # the string $'b is then quoted as '$'"'"'b'
        return "'" + s.replace("'", "'\"'\"'") + "'"
