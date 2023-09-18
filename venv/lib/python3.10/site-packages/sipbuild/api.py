# Copyright (c) 2022, Riverbank Computing Limited
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


from .abstract_project import AbstractProject
from .exceptions import handle_exception


def build_sdist(sdist_directory, config_settings=None):
    """ The PEP 517 hook for building an sdist from pyproject.toml. """

    project = AbstractProject.bootstrap('sdist',
            arguments=_convert_config_settings(config_settings))

    # pip executes this in a separate process and doesn't handle exceptions
    # very well.  However it does capture stdout and (eventually) show it to
    # the user so we use our standard exception handling.
    try:
        return project.build_sdist(sdist_directory)
    except Exception as e:
        handle_exception(e)


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    """ The PEP 517 hook for building a wheel from pyproject.toml. """

    project = AbstractProject.bootstrap('wheel',
            arguments=_convert_config_settings(config_settings))

    # pip executes this in a separate process and doesn't handle exceptions
    # very well.  However it does capture stdout and (eventually) show it to
    # the user so we use our standard exception handling.
    try:
        return project.build_wheel(wheel_directory)
    except Exception as e:
        handle_exception(e)


def _convert_config_settings(config_settings):
    """ Return any configuration settings from the frontend to a pseudo-command
    line.
    """

    if config_settings is None:
        config_settings = {}

    args = []

    for name, value in config_settings.items():
        if value:
            if not isinstance(value, list):
                value = [value]

            for m_value in value:
                args.append(name + '=' + m_value)
        else:
            args.append(name)

    return args
