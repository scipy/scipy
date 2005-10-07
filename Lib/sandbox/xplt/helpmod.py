# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# A general purpose, crude but hopefully effective help function.
# This expects to find help files on the sys.path named "<module>.help",
# internally formatted as follows:
# 
#   <some_feature>:
#     Documentation for some_feature. The word <some_feature> begins in
#     column 1 and ends with ``:'' and a newline. To find <some_feature>
#     type ``help("some_feature")''. Notice the indentation of this text
#     inward from <some_feature>. This is not strictly necessary, but it
#     helps keep clear which things are features.
#
#   <some_other_feature>:
#     To find this particular feature, the help function would be
#     called with help("some_other_feature").  It would construct the
#     search pattern '^some_other_feature:' and pass that to the "more"
#     program.
#
# You can call help with three general forms:
#   from help import help
#   1) help("module.feature"): Look for a given feature in "module.help".
#   2) help("module."): Open "module.help" at line 1. (Notice the dot.)
#   3) help("feature"): Look for a feature in any help file, beginning
#        with the one last accessed.

import os, sys, string

_current_module = ""
_current_dir = ""

#*************
# This function should be the only external object in this module.
def help (subj = "help.help"):
  global _current_module, _current_dir
  sys.path.append(os.environ['GISTPATH'])

  if string.find (subj,'.') == -1: # subject pattern has no dot:

    rexp = "^" + subj + ":"
    p = _find_module_with_feature(rexp)

  else: # subject pattern does have a dot, e.g. help("foo.bar")

    [ module, feature ] = string.splitfields (subj, '.')
    p = _find_module(module)
    if feature == "":
      rexp = "+1 "
    else:
      rexp = "^" + feature + ":"

  if p != "" and _ispresent(rexp,p):
    _pageit (rexp, p)

  del sys.path[-1]

#*************
# Look in all help files for a particular feature
def _find_module_with_feature (rexp):
  global _current_module, _current_dir

  # Look first in the cached help file.
  p = os.path.join (_current_dir, _current_module) + ".help"
  if _ispresent(rexp, p):
    return p

  # Extend the search to all help files.
  for d in sys.path:
    xp = os.path.join (d, "*.help")
    s = "egrep -l '" + rexp + "' " + xp + " 2>/dev/null"
    pipe = os.popen (s)
    line = pipe.readline()
    pipe.close()
    line = string.strip(line)
    if line != '':	# Found one
      _current_dir = d
      _current_module = line[0:-5]
      return os.path.join(d, line)
  
  # Failure return
  print "could not find pattern", rexp, "in any help file on your sys.path"
  return ""

#*************
# Look for a particular help file "m.help" in user's sys.path
def _find_module (m):
  global _current_module, _current_dir

  # Check cache first.
  f = m + ".help"
  if m == _current_module:
    return os.path.join(_current_dir, f)

  # Extend the search over all components.
  for d in sys.path:
    p = os.path.join(d,f)
    if os.path.exists (p):
      _current_dir = d
      _current_module = m
      return p
  
  # Failure return
  print "could not find", f, "on your sys.path"
  _current_dir = _current_module = ""
  return ""

#*************
# Check a given file to see if it contains a given regular expression
def _ispresent (rexp, file):
  if rexp == "+1 ":
    return 1
  s = "egrep -l '" + rexp + "' " + file + ">/dev/null 2>&1"
  sts = os.system (s)
  if sts == 0:
    return 1
  # Failure return
  return 0

#*************
# Page the given file, starting at the given regular expression.
# This function should always succeed.
def _pageit( rexp, file):
    if rexp != "+1 ":
      rexp = "+/'" + rexp + "' "
    s = "${PAGER-more} " + rexp + file
    sts = os.system (s)
    if sts != 0:
      print "Warning: Pager returned non-zero status:", sts, "on file", file
