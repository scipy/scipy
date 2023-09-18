# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains hooks for zest.releaser for use in semi-automated releases
of asv.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import sys


def releaser_middle(data):
    """
    releaser.middle hook to monkey-patch zest.releaser to support signed
    tagging--currently this is the only way to do this.  Also monkey-patches to
    disable an annoyance where zest.releaser only creates .zip source
    distributions.  This is supposedly a workaround for a bug in Python 2.4,
    but we don't care about Python 2.4.
    """

    if data['name'] != 'asv':
        return

    from zest.releaser.git import Git
    from zest.releaser.release import Releaser

    # Copied verbatim from zest.releaser, but with the cmd string modified to
    # use the -s option to create a signed tag and add the 'v' in front of the
    # version number
    def _my_create_tag(self, version):
        version = 'v' + version
        msg = "Tagging %s" % (version,)
        cmd = 'git tag -s %s -m "%s"' % (version, msg)
        if os.path.isdir('.git/svn'):
            print("\nEXPERIMENTAL support for git-svn tagging!\n")
            cur_branch = open('.git/HEAD').read().strip().split('/')[-1]
            print("You are on branch %s." % (cur_branch,))
            if cur_branch != 'master':
                print("Only the master branch is supported for git-svn tagging.")
                print("Please tag yourself.")
                print("'git tag' needs to list tag named %s." % (version,))
                sys.exit()
            cmd = [cmd]
            local_head = open('.git/refs/heads/master').read()
            trunk = open('.git/refs/remotes/trunk').read()
            if local_head != trunk:
                print("Your local master diverges from trunk.\n")
                # dcommit before local tagging
                cmd.insert(0, 'git svn dcommit')
            # create tag in svn
            cmd.append('git svn tag -m "%s" %s' % (msg, version))
        return cmd

    # Similarly copied from zest.releaser to support use of 'v' in front
    # of the version number
    def _my_make_tag(self):
        from zest.releaser import utils

        if self.data['tag_already_exists']:
            return
        cmds = self.vcs.cmd_create_tag(self.data['version'])
        if not isinstance(cmds, list):
            cmds = [cmds]
        if len(cmds) == 1:
            print("Tag needed to proceed, you can use the following command:")
        for cmd in cmds:
            print(cmd)
            if utils.ask("Run this command"):
                print(os.system(cmd))
            else:
                # all commands are needed in order to proceed normally
                print("Please create a tag for %s yourself and rerun." % \
                        (self.data['version'],))
                sys.exit()
        if not self.vcs.tag_exists('v' + self.data['version']):
            print("\nFailed to create tag %s!" % (self.data['version'],))
            sys.exit()

    # Normally all this does is to return '--formats=zip', which is currently
    # hard-coded as an option to always add to the sdist command; they ought to
    # make this actually optional
    def _my_sdist_options(self):
        return ''

    Git.cmd_create_tag = _my_create_tag
    Releaser._make_tag = _my_make_tag
    Releaser._sdist_options = _my_sdist_options


_NEW_CHANGELOG_TEMPLATE = str("""\
New Features
^^^^^^^^^^^^

API Changes
^^^^^^^^^^^

Bug Fixes
^^^^^^^^^

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^
""".rstrip())


def postreleaser_before(data):
    """
    postreleaser.before hook to set a different dev_version_template from the
    default: By default zest.releaser uses <version>.dev0.  We want just
    <version>.dev without the mysterious 0.
    """

    if data['name'] != 'asv':
        return

    data['dev_version_template'] = '%(new_version)s.dev'
    data['nothing_changed_yet'] = _NEW_CHANGELOG_TEMPLATE
