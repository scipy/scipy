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
import sys


class UserException(Exception):
    """ An exception capturing user friendly information. """

    def __init__(self, text, *, detail=None):
        """ Initialise the exception with its user friendly text and the
        optional detail.
        """

        super().__init__()

        self.text = text
        self.detail = detail


class UserFileException(UserException):
    """ An exception related to a specific file. """

    def __init__(self, filename, text, *, detail=None):
        """ Initialise the exception. """

        super().__init__("{0}: {1}".format(filename, text), detail=detail)


class UserParseException(UserFileException):
    """ An exception related to the parsing of a specific file. """

    def __init__(self, filename, *, detail=None):
        """ Initialise the exception. """

        super().__init__(filename, "unable to parse file", detail=detail)


def handle_exception(e):
    """ Tell the user about an exception. """

    if isinstance(e, UserException):
        # An "expected" exception.
        if e.detail is not None:
            message = "{0}: {1}".format(e.text, e.detail)
        else:
            message = e.text

        print("{0}: {1}".format(os.path.basename(sys.argv[0]), message),
                file=sys.stderr)

        sys.exit(1)

    # An internal error.
    print(
            "{0}: An internal error occurred...".format(
                    os.path.basename(sys.argv[0])),
            file=sys.stderr)

    raise e
