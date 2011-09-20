#!/usr/bin/env python
# -*- encoding:utf-8 -*-
"""
git-authors [OPTIONS] REV1..REV2

List the authors who contributed within a given revision interval.

"""
# Author: Pauli Virtanen <pav@iki.fi>. This script is in the public domain.

from subprocess import Popen, PIPE, call
import tempfile
import optparse
import re
import os
import subprocess

NAME_MAP = {
    u'warren.weckesser': u'Warren Weckesser',
    u'josef': u'Josef Perktold',
    u'josef-pktd': u'Josef Perktold',
    u'Gael varoquaux': u'Gaël Varoquaux',
    u'rgommers': u'Ralf Gommers',
    u'ArmstrongJ': u'Jeff Armstrong',
    u'Collin Stocks': u'Collin RM Stocks',
    u'gotgenes': u'Chris Lasher',
    u'dsimcha': u'David Simcha',
    u'cgholke': u'Christoph Gohlke',
    u'Christolph Gohlke': u'Christoph Gohlke',
    u'dhuard': u'David Huard',
    u'pierregm': u'Pierre GM',
    u'Derek Homier': u'Derek Homeier',
    u'Derek Homeir': u'Derek Homeier',
    u'87': u'Han Genuit',
    u'Han': u'Han Genuit',
    u'weathergod': u'Benjamin Root',
    u'Mark': u'Mark Wiebe',
    u'sebhaase': u'Sebastian Haase',
    u'mdroe': u'Michael Droettboom',
    u'chris.burns': u'Chris Burns',
    u'aarchiba': u'Anne Archibald',
    u'edschofield': u'Ed Schofield',
}

def main():
    p = optparse.OptionParser(__doc__.strip())
    p.add_option("-d", "--debug", action="store_true",
                 help="print debug output")
    options, args = p.parse_args()

    if len(args) != 1:
        p.error("invalid number of arguments")

    try:
        rev1, rev2 = args[0].split('..')
    except ValueError:
        p.error("argument is not a revision range")

    # Analyze log data
    all_authors = set()
    authors = set()

    def analyze_line(line, names, disp=False):
        line = line.strip().decode('utf-8')

        # Check the commit author name
        m = re.match('^@@@([^@]*)@@@', line)
        if m:
            name = m.group(1)
            line = line[m.end():]
            name = NAME_MAP.get(name, name)
            if disp:
                if name not in names:
                    print "    - Author:", name
            names.add(name)

        # Look for "thanks to" messages in the commit log
        m = re.search(ur'([Tt]hanks to|[Cc]ourtesy of) ([A-Z][A-Za-z]*? [A-Z][A-Za-z]*? [A-Z][A-Za-z]*|[A-Z][A-Za-z ]*? [A-Z][A-Za-z]*|[a-z0-9]+)($|\.| )', line)
        if m:
            name = m.group(2)
            if name not in ('this',):
                if disp:
                    print "    - Log   :", line.strip()
                name = NAME_MAP.get(name, name)
                names.add(name)

            line = line[m.end():].strip()
            line = re.sub(ur'^(and|, and|, ) ', u'Thanks to ', line)
            analyze_line(line, names)

    # Find all authors before the named range
    for line in git.pipe('log', '--pretty=@@@%an@@@%n@@@%cn@@@%n%b',
                         '%s' % (rev1,)):
        analyze_line(line, all_authors)

    # Find authors in the named range
    for line in git.pipe('log', '--pretty=@@@%an@@@%n@@@%cn@@@%n%b',
                         '%s..%s' % (rev1, rev2)):
        analyze_line(line, authors, disp=options.debug)

    # Sort
    def name_key(fullname):
        m = re.search(u' [a-z ]*[A-Za-z-]+$', fullname)
        if m:
            forename = fullname[:m.start()].strip()
            surname = fullname[m.start():].strip()
        else:
            forename = u""
            surname = fullname.strip()
        if surname.startswith('van der '):
            surname = surname[8:]
        if surname.startswith('de '):
            surname = surname[3:]
        if surname.startswith('von '):
            surname = surname[4:]
        return (surname.lower(), forename.lower())

    authors = list(authors)
    authors.sort(key=name_key)

    # Print
    print """
Authors
=======

This release contains work by the following people (contributed at least
one patch to this release, names in alphabetical order):
"""

    for author in authors:
        if author in all_authors:
            print (u"* %s" % author).encode('utf-8')
        else:
            print (u"* %s +" % author).encode('utf-8')

    print """
A total of %(count)d people contributed to this release.
People with a "+" by their names contributed a patch for the first time.
""" % dict(count=len(authors))

    print ("\nNOTE: Check this list manually! It is automatically generated "
           "and some names\n      may be missing.")

#------------------------------------------------------------------------------
# Communicating with Git
#------------------------------------------------------------------------------

class Cmd(object):
    executable = None

    def __init__(self, executable):
        self.executable = executable

    def _call(self, command, args, kw, repository=None, call=False):
        cmd = [self.executable, command] + list(args)
        cwd = None

        if repository is not None:
            cwd = os.getcwd()
            os.chdir(repository)

        try:
            if call:
                return subprocess.call(cmd, **kw)
            else:
                return subprocess.Popen(cmd, **kw)
        finally:
            if cwd is not None:
                os.chdir(cwd)

    def __call__(self, command, *a, **kw):
        ret = self._call(command, a, {}, call=True, **kw)
        if ret != 0:
            raise RuntimeError("%s failed" % self.executable)

    def pipe(self, command, *a, **kw):
        stdin = kw.pop('stdin', None)
        p = self._call(command, a, dict(stdin=stdin, stdout=subprocess.PIPE),
                      call=False, **kw)
        return p.stdout

    def read(self, command, *a, **kw):
        p = self._call(command, a, dict(stdout=subprocess.PIPE),
                      call=False, **kw)
        out, err = p.communicate()
        if p.returncode != 0:
            raise RuntimeError("%s failed" % self.executable)
        return out

    def readlines(self, command, *a, **kw):
        out = self.read(command, *a, **kw)
        return out.rstrip("\n").split("\n")

    def test(self, command, *a, **kw):
        ret = self._call(command, a, dict(stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE),
                        call=True, **kw)
        return (ret == 0)

git = Cmd("git")

#------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

