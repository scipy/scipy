"""

    This is a tool to maintain and run scipy with valgrind against a set of suppression files.

    A suppression file is a list of known memory errors, many of which can be false positives.

    On a fixed platform (python and numpy binary)
    if we see new errors when running against an older version of suppression files, 
    then a leak or memory violation is likely.

    We use the grindmerge tool written here (the second perl script) for the merging:

        https://wiki.wxwidgets.org/Parse_valgrind_suppressions.sh

    example:

    # produce a long list of errors triggered by _lib tests:

    # generate a suppression database for _lib

    python run-valgrind.py generate --python=python3-debug/bin/python3-debug scipy/_lib/tests

    # rerun with diff to find no errors are reported; check the file name printed at the end.

    python run-valgrind.py diff --python=python3-debug/bin/python3-debug scipy/_lib/tests

    # later rerun to find new errors are reported, check for memory issues and fix them

    python run-valgrind.py diff --python=python3-debug/bin/python3-debug scipy/_lib/tests
    ...
    python run-valgrind.py diff --python=python3-debug/bin/python3-debug scipy/_lib/tests
    ...
    python run-valgrind.py diff --python=python3-debug/bin/python3-debug scipy/_lib/tests

    # update the suppresion database
    python run-valgrind.py generate --python=python3-debug/bin/python3-debug scipy/_lib/tests

"""

from argparse import ArgumentParser
import os
from subprocess import call
from tempfile import mkstemp

ap = ArgumentParser("run-valgrind")
ap.add_argument("command", choices=["diff", "generate",])
ap.add_argument("tests", nargs='+', help='a list of tests to run')
ap.add_argument("--prefix", default='valgrind',
        help='specify a python to use, usually better with python-debug')
ap.add_argument("--python", default='python3-debug',
        help='specify a python to use, usually better with python-debug')
ap.add_argument("--valgrind", default='valgrind')
ap.add_argument("--suppressions", default=[], action='append')

def main():
    ns = ap.parse_args()

    if ns.command == 'diff':
        rules(ns, diff=True)
    elif ns.command == 'generate':
        rules(ns, diff=False)

def run_valgrind(valgrind, opts, payload, suppressions):
    fid, logfilename = mkstemp()
    os.close(fid)

    valgrind_opts = opts + [
        '--log-file=%s' % logfilename,
        ]

    for fn in suppressions:
        if not os.path.exists(fn):
            raise ValueError("file %s not found" % fn)

    suppressions = [ '--suppressions=%s' % fn for fn in suppressions ]
    cmdline = valgrind + valgrind_opts + suppressions + payload
    call(cmdline)

    with open(logfilename, 'r') as ff:
        r = ff.read()

    os.unlink(logfilename)
    if 'FATAL' in r:
        print(r)
        raise RuntimeError("Valgrind failed with")
    return r

def run_grindmerge(vallog, oldrulefilename=None):
    fid, logfilename = mkstemp()
    os.close(fid)
    fid, rulefilename = mkstemp()
    os.close(fid)

    with open(logfilename, 'w') as ff:
        ff.write(vallog)

    with open(logfilename, 'r') as stdin:
        with open(rulefilename, 'w') as stdout:

            cmdline = ['valgrind/grindmerge']
            if oldrulefilename:
                cmdline = cmdline + ['-f', oldrulefilename]
            call(cmdline, stdin=stdin, stdout=stdout)

    os.unlink(logfilename)
    with open(rulefilename, 'r') as ff:
        r = ff.read()
    os.unlink(rulefilename)
    return r

def ensure_suppression_dir(prefix, path):
    try:
        os.makedirs(os.path.join(prefix, path))
    except:
        pass

def find_suppression_files(prefix, path):
    if len(path) == 0: return []
    path_suppr = os.path.join(prefix, path, 'valgrind-suppression')
    higher = find_suppression_files(prefix, os.path.dirname(path))
    print('checking' , path_suppr)
    if os.path.exists(path_suppr):
        return higher + [path_suppr]
    else:
        return higher

def find_local_suppression_file(prefix, path):
    path = os.path.join(prefix, path)
    if not os.path.isdir(path): return None
    path_suppr = os.path.join(path, 'valgrind-suppression')
    return path_suppr

def rules(ns, diff=True):
    runtests = [ ns.python, 'runtests.py' ]

    opts = [
        '--show-leak-kinds=all',
        '--leak-check=full',
        '--num-callers=80',
        '--error-limit=no',
        '--fullpath-after=',
        '--gen-suppressions=all',
    ]

    for test in ns.tests:

        ensure_suppression_dir(ns.prefix, test)

        suppressions = ns.suppressions
        all_test_suppr = find_suppression_files(ns.prefix, test)
        print("all suppression files", all_test_suppr)
        my_test_suppr = find_local_suppression_file(ns.prefix, test)
        if not diff: # remove my_test_suppr because we will generate it afresh.
            while my_test_suppr in all_test_suppr:
                all_test_suppr.remove(my_test_suppr)
        else:
            print("using current suppression in ", my_test_suppr, "for testing")

        suppressions = suppressions + all_test_suppr

        print("using suppression files", suppressions)


        print("rebuild the binaries to avoid valgrind erros during building.")
        call(runtests + ['--build'])

        print("running valgrind with the tests")
        log = run_valgrind([ns.valgrind], opts, runtests + ['-t', test], suppressions)

        with open(my_test_suppr + '.log', 'w') as ff:
            ff.write(log)

        leaked, summary = filter_log(log)
        print(summary)

        if leaked:
            print("Check the log file for possible leaks due to scipy;", my_test_suppr + '.log')
        else:
            print("No scipy error is found, merging the rules to the current rule set")

        if diff:
            rules = run_grindmerge(log, my_test_suppr)

            if not leaked:
                print("merging new errors to ", my_test_suppr)
                with open(my_test_suppr, 'w') as ff:
                    ff.write(rules)

        else:
            if leaked:
                print("the suppresion file may contain existing scipy triggered memory errors.")
            rules = run_grindmerge(log, None)
            print("writing newly generated suppression rules to ", my_test_suppr)
            with open(my_test_suppr, 'w') as ff:
                ff.write(rules)

def filter_log(log):
    any = False
    summary = []
    for i, line in enumerate(log.split('\n')):
        if i < 5: continue # by pass header
        if 'scipy/' in line:
            print('potential leak found at line :', i, line)
            any = True
        if 'LEAK SUMMARY' in line or len(summary) > 0:
            summary.append(line)
    return any, '\n'.join(summary)
main()
