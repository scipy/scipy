#!/usr/bin/env python

__usage__ = '''
Filters filenames and commit messages from `cvs log` command
to stdout.

Usage:
  cvs log -d">2004/04/15" | /path/to/cvs_log_cleanup.py
'''

import sys
import re
line = 1
working_file = re.compile(r'Working file: (?P<filename>.*)').match
descr_start = re.compile(r'----------------------------').match
descr_end = re.compile(r'============================').match

while line:
    line = sys.stdin.readline()
    m = working_file(line)
    if m:
        filename = m.group('filename')
        descr = []
        while line and not (descr_start(line) or descr_end(line)):
            line = sys.stdin.readline()
        if descr_start(line):
            while line and not descr_end(line):
                line = sys.stdin.readline()
                line = sys.stdin.readline()
                line = sys.stdin.readline()
                descr.append('  - ')
                while line and not (descr_start(line) or descr_end(line)):
                    descr.append(line)
                    line = sys.stdin.readline()
            print '\n*',filename
            print ''.join(descr).rstrip()
