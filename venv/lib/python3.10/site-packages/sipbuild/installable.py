# Copyright (c) 2019, Riverbank Computing Limited
# All rights reserved.
#
# This copy of SIP is licensed for use under the terms of the SIP License
# Agreement.  See the file LICENSE for more details.
#
# This copy of SIP may also used under the terms of the GNU General Public
# License v2 or v3 as published by the Free Software Foundation which can be
# found in the files LICENSE-GPL2 and LICENSE-GPL3 included in this package.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


import os
from shutil import copy2, copytree


class Installable:
    """ Encapsulate a list of files and directories that will be installed into
    a target directory.
    """

    def __init__(self, name, *, target_subdir=None):
        """ Initialise the installable.  The optional target_subdir is the
        path of a sub-directory of the eventual target where the files will be
        installed.  If target_subdir is an absolute pathname then it is used as
        the eventual target name.
        """

        self.name = name
        self.target_subdir = target_subdir
        self.files = []

    def get_full_target_dir(self, target_dir):
        """ Return the full target directory name. """

        if self.target_subdir:
            if os.path.isabs(self.target_subdir):
                target_dir = self.target_subdir
            else:
                target_dir = os.path.join(target_dir, self.target_subdir)

        return target_dir

    def install(self, target_dir, installed, *, do_install=True):
        """ Optionally install the files in a target directory and update the
        given list of installed files.
        """

        target_dir = self.get_full_target_dir(target_dir)

        if do_install:
            os.makedirs(target_dir, exist_ok=True)

        for fn in self.files:
            t_path = os.path.join(target_dir, os.path.basename(fn))
            installed.append(t_path)

            if do_install:
                copy_fn = copytree if os.path.isdir(fn) else copy2
                copy_fn(fn, t_path)
