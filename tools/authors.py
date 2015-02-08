#!/usr/bin/env python
# -*- encoding:utf-8 -*-
"""
git-authors [OPTIONS] REV1..REV2

List the authors who contributed within a given revision interval.

"""
# Author: Pauli Virtanen <pav@iki.fi>. This script is in the public domain.

from __future__ import division, print_function, absolute_import

import optparse
import re
import sys
import os
import subprocess

try:
    from scipy._lib.six import u, PY3
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                    os.pardir, 'scipy', 'lib'))
    from six import u, PY3
if PY3:
    stdout_b = sys.stdout.buffer
else:
    stdout_b = sys.stdout


NAME_MAP = {
    u('87'): u('Han Genuit'),
    u('aarchiba'): u('Anne Archibald'),
    u('alex'): u('Alex Griffing'),
    u('andbo'): u('Anders Bech Borchersen'),
    u('argriffing'): u('Alex Griffing'),
    u('arichar6'): u('Steve Richardson'),
    u('ArmstrongJ'): u('Jeff Armstrong'),
    u('Benny'): u('Benny Malengier'),
    u('brettrmurphy'): u('Brett R. Murphy'),
    u('cgholke'): u('Christoph Gohlke'),
    u('cgohlke'): u('Christoph Gohlke'),
    u('chris.burns'): u('Chris Burns'),
    u('Christolph Gohlke'): u('Christoph Gohlke'),
    u('ckuster'): u('Christopher Kuster'),
    u('Collin Stocks'): u('Collin RM Stocks'),
    u('cnovak'): u('Clemens Novak'),
    u('ctokheim'): u('Collin Tokheim'),
    u('Daniel Smith'): u('Daniel B. Smith'),
    u('Dapid'): u('David Menendez Hurtado'),
    u('dellsystem'): u('Wendy Liu'),
    u('Derek Homeir'): u('Derek Homeier'),
    u('Derek Homier'): u('Derek Homeier'),
    u('DSG User'): u('Max Bolingbroke'),
    u('dhuard'): u('David Huard'),
    u('dsimcha'): u('David Simcha'),
    u('edschofield'): u('Ed Schofield'),
    u('Gael varoquaux'): u('Gaël Varoquaux'),
    u('gotgenes'): u('Chris Lasher'),
    u('Han'): u('Han Genuit'),
    u('Helder'): u('Helder Cesar'),
    u('Horta'): u('Danilo Horta'),
    u('Jake Vanderplas'): u('Jacob Vanderplas'),
    u('jamestwebber'): u('James T. Webber'),
    u('jaimefrio'): u('Jaime Fernandez del Rio'),
    u('janani'): u('Janani Padmanabhan'),
    u('Janani'): u('Janani Padmanabhan'),
    u('jesseengel'): u('Jesse Engel'),
    u('josef'): u('Josef Perktold'),
    u('josef-pktd'): u('Josef Perktold'),
    u('kat'): u('Kat Huang'),
    u('Mark'): u('Mark Wiebe'),
    u('mdroe'): u('Michael Droettboom'),
    u('patricksnape'): u('Patrick Snape'),
    u('pbrod'): u('Per Brodtkorb'),
    u('pierregm'): u('Pierre GM'),
    u('polyatail'): u('Andrew Sczesnak'),
    u('rgommers'): u('Ralf Gommers'),
    u('sebhaase'): u('Sebastian Haase'),
    u('SytseK'): u('Sytse Knypstra'),
    u('Takuya OSHIMA'): u('Takuya Oshima'),
    u('tiagopereira'): u('Tiago M.D. Pereira'),
    u('tonysyu'): u('Tony S. Yu'),
    u('Travis E. Oliphant'): u('Travis Oliphant'),
    u('warren.weckesser'): u('Warren Weckesser'),
    u('weathergod'): u('Benjamin Root'),
    u('Andreas H'): u('Andreas Hilboll'),
    u('honnorat'): u('Marc Honnorat'),
    u('lmwang'): u('Liming Wang'),
    u('wa03'): u('Josh Lawrence'),
    u('loluengo'): u('Lorenzo Luengo'),
    u('Zhenya'): u('Evgeni Burovski'),
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
        m = re.match(u('^@@@([^@]*)@@@'), line)
        if m:
            name = m.group(1)
            line = line[m.end():]
            name = NAME_MAP.get(name, name)
            if disp:
                if name not in names:
                    stdout_b.write(("    - Author: %s\n" % name).encode('utf-8'))
            names.add(name)

        # Look for "thanks to" messages in the commit log
        m = re.search(u(r'([Tt]hanks to|[Cc]ourtesy of) ([A-Z][A-Za-z]*? [A-Z][A-Za-z]*? [A-Z][A-Za-z]*|[A-Z][A-Za-z]*? [A-Z]\. [A-Z][A-Za-z]*|[A-Z][A-Za-z ]*? [A-Z][A-Za-z]*|[a-z0-9]+)($|\.| )'), line)
        if m:
            name = m.group(2)
            if name not in (u('this'),):
                if disp:
                    stdout_b.write("    - Log   : %s\n" % line.strip().encode('utf-8'))
                name = NAME_MAP.get(name, name)
                names.add(name)

            line = line[m.end():].strip()
            line = re.sub(u(r'^(and|, and|, ) '), u('Thanks to '), line)
            analyze_line(line.encode('utf-8'), names)

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
        m = re.search(u(' [a-z ]*[A-Za-z-]+$'), fullname)
        if m:
            forename = fullname[:m.start()].strip()
            surname = fullname[m.start():].strip()
        else:
            forename = ""
            surname = fullname.strip()
        if surname.startswith(u('van der ')):
            surname = surname[8:]
        if surname.startswith(u('de ')):
            surname = surname[3:]
        if surname.startswith(u('von ')):
            surname = surname[4:]
        return (surname.lower(), forename.lower())

    authors = list(authors)
    authors.sort(key=name_key)

    # Print
    stdout_b.write(b"""
Authors
=======

""")

    for author in authors:
        if author in all_authors:
            stdout_b.write(("* %s\n" % author).encode('utf-8'))
        else:
            stdout_b.write(("* %s +\n" % author).encode('utf-8'))

    stdout_b.write(("""
A total of %(count)d people contributed to this release.
People with a "+" by their names contributed a patch for the first time.
This list of names is automatically generated, and may not be fully complete.

""" % dict(count=len(authors))).encode('utf-8'))

    stdout_b.write(("\nNOTE: Check this list manually! It is automatically generated "
                    "and some names\n      may be missing.\n").encode('utf-8'))

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
