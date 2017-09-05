from __future__ import print_function
"""

    This is a tool to maintain and run scipy with valgrind against a set of suppression files.

    A suppression file is a list of known memory errors, many of which can be false positives.

    We scan valgrind log to identify potential erros related to scipy, and maintain a suppression
    db for non-scipy errors.

    Example:

    # Find scipy related errors and update the non-scipy suppression db.

    python tools/run-valgrind.py --python=python3-debug/bin/python3-debug scipy/_lib/tests/test__gcutils.py

    # Find scipy related errors and replace the non-scipy suppression db.

    python tools/run-valgrind.py --update-supp=replace --python=python3-debug/bin/python3-debug scipy/_lib/tests

    # Find scipy related errors and do not update the non-scipy suppression db.

    python tools/run-valgrind.py --update-supp=no --python=python3-debug/bin/python3-debug scipy/_lib/tests


    The errors and suppression files for scipy and non-scipy entries are stored in the valgrind/ directory for
    the test case. 

    Selected rules can be manually merged into the default file valgrind-suppression in valgrind directory.

"""

from argparse import ArgumentParser
import os
from subprocess import call
from tempfile import mkstemp
import time

ap = ArgumentParser("run-valgrind")
ap.add_argument("tests", nargs='+', help='a list of tests to run')
ap.add_argument("--update-supp", choices=["merge", "replace", "no"], default='merge', help='strategy for merging non-scipy valgrind errors')
ap.add_argument("--prefix", default='valgrind',
        help='specify a python to use, usually better with python-debug')
ap.add_argument("--python", default='python3-debug',
        help='specify a python to use, usually better with python-debug')
ap.add_argument("--valgrind", default='valgrind')
ap.add_argument("--suppressions", default=[], action='append')

def main():
    ns = ap.parse_args()

    rules(ns, update_suppressions=ns.update_supp)

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
    #print("writing log to", logfilename)
    if 'FATAL' in r:
        print(r)
        raise RuntimeError("Valgrind failed with")
    return r

def ensure_suppression_dir(prefix, path):
    try:
        os.makedirs(os.path.join(prefix, path))
    except:
        pass

def find_suppression_files(prefix, path):
    path_suppr = os.path.join(prefix, path, 'valgrind-suppression')
    if len(path) == 0:
        higher = []
    else:
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

def rules(ns, update_suppressions):
    runtests = [ ns.python, 'runtests.py' ]

    opts = [
        '--show-leak-kinds=all',
        '--leak-check=full',
        '--num-callers=40',
        '--error-limit=no',
        '--fullpath-after=', # for scipy detection
        '--gen-suppressions=all', # for suppressions generation
    ]

    for test in ns.tests:
        t0 = time.time()

        ensure_suppression_dir(ns.prefix, test)

        suppressions = ns.suppressions
        all_test_suppr = [os.path.join('tools', 'scipy-master.supp')] + find_suppression_files(ns.prefix, test)
        print("all suppression files", all_test_suppr)

        suppressions = suppressions + all_test_suppr

        print("using suppression files", suppressions)

        if update_suppressions == "replace":
            local_supp = find_local_suppression_file(ns.prefix, test)
            while local_supp in suppressions:
                suppressions.remove(local_supp)

        print("rebuild the binaries to avoid valgrind erros during building.")
        if 0 != call(runtests + ['--build']):
            raise RuntimeError("building failed")

        print("running valgrind with the tests")
        log = run_valgrind([ns.valgrind], opts, runtests + ['-t', test], suppressions)

        #log = open('log-example').read()
        vlog = ValgrindLog.from_string(log)

        scipy_errors = ValgrindLog()
        non_scipy_errors = ValgrindLog()

        for section in vlog:
            if section.is_heap_summary():
                print('\n'.join(section))
            if section.is_error_summary():
                print('\n'.join(section))
            if section.is_leak_summary():
                print('\n'.join(section))

            sc = section.get_scipy_related()
            if sc is not None:
                scipy_errors.append(section)
            else:
                non_scipy_errors.append(section)

        print('Found %d valgrind anomalies that appeared to be related to scipy' % len(scipy_errors))
        if len(scipy_errors):
            print(str(scipy_errors))

        print('Found %d valgrind anomalies that appeared to be unrelated to scipy' % len(non_scipy_errors))

        with open(os.path.join(ns.prefix, test, 'scipy.log'), 'w') as ff:
            ff.write(str(scipy_errors))
        print("Scipy error log ", os.path.join(ns.prefix, test, 'scipy.log'))

        with open(os.path.join(ns.prefix, test, 'scipy.supp'), 'w') as ff:
            ff.write(str(scipy_errors.get_suppression_db()))
        print("Scipy suppression rules", os.path.join(ns.prefix, test, 'scipy.supp'))

        with open(os.path.join(ns.prefix, test, 'nonscipy.log'), 'w') as ff:
            ff.write(str(non_scipy_errors))
        print("Non-scipy error log ", os.path.join(ns.prefix, test, 'nonscipy.log'))

        with open(os.path.join(ns.prefix, test, 'nonscipy.supp'), 'w') as ff:
            ff.write(str(non_scipy_errors.get_suppression_db()))
        print("Non-scipy suppression rules ", os.path.join(ns.prefix, test, 'nonscipy.supp'))

        if update_suppressions != 'no':
            local_supp = find_local_suppression_file(ns.prefix, test)
            newdb = non_scipy_errors.get_suppression_db()

            print("Found %d suppression rules" % len(newdb))

            if update_suppressions == 'replace':
                pass
            elif update_suppressions == 'merge':
                try:
                    db = SuppressionDB.fromfile(local_supp)
                except IOError:
                    db = SuppressionDB()

                print("Merging existing %d suppression into %d rules" %( len(db), len(newdb)))

                newdb.update(db)

            print("Written %d suppression rules to %s" % (len(newdb), local_supp))
            with open(local_supp, 'w') as ff:
                ff.write(str(newdb))

        t1 = time.time()
        print("Testing %s used %g seconds" % (test, t1 - t0))

class ValgrindSection(list):
    def __init__(self):
        self.supp_rule = []

    @classmethod
    def from_list(cls, list, supp_rule):
        self = ValgrindSection()
        self.extend(list)
        self.supp_rule = supp_rule
        return self

    def is_warn(self):
        if len(self) == 0: return False
        return self[0].startswith('Warning:')

    def is_heap_summary(self):
        if len(self) == 0: return False
        for line in self:
            if line.startswith("HEAP SUMMARY:"):
                return True

    def is_leak_summary(self):
        if len(self) == 0: return False
        for line in self:
            if line.startswith("LEAK SUMMARY:"):
                return True

    def is_error_summary(self):
        if len(self) == 0: return False
        for line in self:
            if line.startswith("ERROR SUMMARY:"):
                return True

    def is_entry(self):
        return len(self.supp_rule) > 0

    def __str__(self):
        return '\n'.join(self)

    def get_scipy_related(self):
        if not self.is_entry(): return None
        r = []
        for line in self:
            if 'scipy/' in line:
                r.append(line)
        if len(r) is 0: return None
        return ValgrindSection.from_list(r, self.supp_rule)

    def format_suppression(self):
        return '\n'.join(self.supp_rule)

class ValgrindLog(list):
    def __init__(sections):
        pass

    @classmethod
    def from_string(kls, log):
        sections = kls()
        section_start = True
        section = ValgrindSection()
        for i, line in enumerate(log.split('\n')):
            if line.startswith('=='):
                pid, line = line.split(' ', 1)
            else:
                pid = None

            if pid is not None:
                if section_start:
                    sections.append(section)
                    section = ValgrindSection()
                    section_start = False

                if len(line.strip()) == 0:
                    section_start = True
                else:
                    section.append(line.strip())

            if pid is None:
                section.supp_rule.append(line.rstrip())
        return sections

    def get_suppression_db(self):
        db = SuppressionDB()
        for section in self:
            db.add(section.format_suppression())
        return db

    @classmethod
    def fromfile(cls, filename):
        return cls(open(filename).read())

    def __str__(self):
        return '\n'.join([str(i) for i in self])

class SuppressionDB(set):
    @classmethod
    def fromfile(cls, filename):
        self = cls()
        rule = None
        for i , line in enumerate(open(filename).readlines()):
            if len(line.strip()) == 0:
                continue
            if line.strip().startswith('{'):
                rule = []

            rule.append(line.rstrip())

            if line.strip().startswith('}'):
                self.add('\n'.join(rule))

        return self

    def __str__(self):
        return '\n'.join(sorted(self))

main()
"""
v = ValgrindLog.fromfile('log-example')

supp = SuppressionDB()
for s in v:
    sc = s.get_scipy_related()
    if sc:
        print(s)
    supp.add(s.format_suppression())

print(supp)
"""
