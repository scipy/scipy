#!/usr/bin/env python
"""Usage: python suppress_output.py COMMAND...

Run the given command and print "in progress..." messages to stderr,
one per minute, as long as the command is producing some output.  If
the command has not recently produced any output, no further messages
will be printed until it does.

When the command exits with return code 0, exit with code 0.

When the command exits with any other return code, print all output
produced by the command and exit with the same return code.

"""
from __future__ import division, absolute_import, print_function

import sys
import os
import re
import subprocess
import signal
import time
import tempfile
import shutil


TIMEOUT = 60


def main():
    command = sys.argv[1:]

    with tempfile.TemporaryFile("a+b", 0) as log:
        p = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)

        # Rather than handling some signals ourselves, forward them.

        def forward_signal(signum, frame):
            p.send_signal(signum)

        signal.signal(signal.SIGINT, forward_signal)
        signal.signal(signal.SIGTERM, forward_signal)
        signal.signal(signal.SIGQUIT, forward_signal)
        signal.signal(signal.SIGHUP, forward_signal)

        # Wait for it to finish, and print something to indicate the
        # process is alive, but only if the log file has grown (to
        # allow continuous integration environments kill a hanging
        # process accurately if it produces no output).

        start_time = time.time()
        last_blip = time.time()
        last_log_size = log.tell()
        counter = [0]

        # Just poll it -- other approaches are more complex

        try:
            sleep_time = 0.001
            while p.poll() is None:
                sleep_time = min(2*sleep_time, 0.5)
                time.sleep(sleep_time)
                if time.time() - last_blip > TIMEOUT:
                    log_size = log.tell()
                    if log_size > last_log_size:
                        msg = "    ... in progress ({0} elapsed)".format(elapsed(time.time() - start_time))
                        print(msg, file=sys.stderr)
                        sys.stderr.flush()
                        counter[0] += 1
                        last_blip = time.time()
                        last_log_size = log_size

            ret = p.wait()
        except:  # noqa: E722
            p.terminate()
            raise

        if counter[0] > 0:
            if ret == 0:
                msg = "    ... ok ({0} elapsed)".format(elapsed(time.time() - start_time))
                print(msg, file=sys.stderr)
                sys.stderr.flush()
            else:
                msg = "    ... failed ({0} elapsed, exit code {1})".format(
                    elapsed(time.time() - start_time), ret)
                print(msg, file=sys.stderr)
                sys.stderr.flush()

        if ret != 0:
            log.seek(0)
            if sys.version_info[0] >= 3:
                shutil.copyfileobj(log, sys.stdout.buffer)
            else:
                shutil.copyfileobj(log, sys.stdout)

    sys.exit(ret)


def elapsed(t):
    if t < 0:
        sgn = '-'
        t = abs(t)
    else:
        sgn = ''

    if t < 60:
        return "{0}{1:.0f} s".format(sgn, round(t))
    elif t < 3600:
        mins, secs = divmod(t, 60)
        return "{0}{1:.0f} min {2:.0f} s".format(sgn, mins, secs)
    else:
        hours, mins = divmod(t, 3600)
        mins, secs = divmod(mins, 60)
        return "{0}{1:.0f} h {2:.0f} min {3:.0f} s".format(sgn, hours, mins, secs)


def test_elapsed():
    assert elapsed(0.4) == '0 s'
    assert elapsed(30.3) == '30 s'
    assert elapsed(59.5) == '60 s'
    assert elapsed(60.5) == '1 min 0 s'
    assert elapsed(2*60 + 40.51) == '2 min 41 s'
    assert elapsed(60 * 59.999) == '59 min 60 s'
    assert elapsed(60 * 60.0) == '1 h 0 min 0 s'
    assert elapsed(266*3600 + 13*60 + 12.4243) == '266 h 13 min 12 s'


def test_exitcode():
    r0 = subprocess.call([sys.executable, __file__, sys.executable, '-c',
                          'import sys; sys.exit(0)'])
    assert r0 == 0
    r1 = subprocess.call([sys.executable, __file__, sys.executable, '-c',
                          'import sys; sys.exit(1)'])
    assert r1 == 1
    rs = subprocess.call([sys.executable, __file__, sys.executable, '-c',
                          'import os; os.kill(os.getpid(), 15)'])
    assert rs != 0


def test_suppress(tmpdir):
    p = subprocess.Popen([sys.executable, __file__, sys.executable, '-c',
                          'import sys; '
                          'sys.stdout.write("OUT"); '
                          'sys.stderr.write("ERR"); '
                          'sys.exit(0)'],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    assert p.returncode == 0
    assert out == b''
    assert err == b''

    p = subprocess.Popen([sys.executable, __file__, sys.executable, '-c',
                          'import sys; '
                          'sys.stdout.write("OUT"); '
                          'sys.stdout.flush(); '
                          'sys.stderr.write("ERR"); '
                          'sys.exit(1)'],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    assert p.returncode == 1
    assert out == b'OUTERR'
    assert err == b''


def run_script_fast(path, script):
    fn = os.path.join(path, 'suppress_output.py')

    with open(__file__, 'rb') as f:
        text = f.read()
        text = text.replace(b'TIMEOUT = 60', b'TIMEOUT = 1')

    with open(fn, 'wb') as f:
        f.write(text)

    p = subprocess.Popen([sys.executable, fn, sys.executable, '-c', script],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return p.returncode, out, err


def test_suppress_long_ok(tmpdir):
    returncode, out, err = run_script_fast(str(tmpdir),
                                           'import sys, time; '
                                           'sys.stdout.write("OUT"); '
                                           'time.sleep(1.5); '
                                           'sys.stderr.write("ERR"); '
                                           'sys.exit(0)')
    assert returncode == 0
    assert out == b''
    assert re.match(b'^    \.\.\. in progress \([0-9 sminh]* elapsed\)\n'
                    b'    \.\.\. ok \([0-9 sminh]* elapsed\)\n$', err,
                    re.S)


def test_suppress_long_failed(tmpdir):
    returncode, out, err = run_script_fast(str(tmpdir),
                                           'import sys, time; '
                                           'sys.stdout.write("OUT"); '
                                           'time.sleep(1.5); '
                                           'sys.stdout.flush(); '
                                           'sys.stderr.write("ERR"); '
                                           'sys.exit(1)')
    assert returncode == 1
    assert out == b'OUTERR'
    assert re.match(b'^    \.\.\. in progress \([0-9 sminh]* elapsed\)\n'
                    b'    \.\.\. failed \([0-9 sminh]* elapsed, exit code 1\)\n$', err,
                    re.S)


if __name__ == "__main__":
    main()
