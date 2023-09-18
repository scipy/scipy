# Copyright (c) 2023, Riverbank Computing Limited
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


from .exceptions import UserFileException
from .py_versions import OLDEST_SUPPORTED_MINOR
from .toml import toml_load


class PyProjectException(UserFileException):
    """ An exception related to a pyproject.toml file. """

    def __init__(self, text, *, detail=None):
        """ Initialise the exception. """

        super().__init__("pyproject.toml", text, detail=detail)


class PyProjectOptionException(PyProjectException):
    """ An exception related to a specific option of a pyproject.toml file. """

    def __init__(self, name, text, *, section_name=None, detail=None):
        """ Initialise the exception. """

        if section_name is None:
            section_name = 'tool.sip.project'

        super().__init__("'{0}.{1}': {2}".format(section_name, name, text),
                detail=detail)


class PyProjectUndefinedOptionException(PyProjectOptionException):
    """ An exception related to an undefined option of a pyproject.toml file.
    """

    def __init__(self, name, *, section_name=None):
        """ Initialise the exception. """

        super().__init__(name, "must be defined", section_name=section_name)


class PyProject:
    """ Encapsulate a parsed pyproject.toml file. """

    def __init__(self):
        """ Initialise the object. """

        self.toml_error = None

        try:
            self._pyproject = toml_load('pyproject.toml')
        except FileNotFoundError:
            self.toml_error = "there is no such file in the current directory"
        except Exception as e:
            self.toml_error = str(e)

    def get_metadata(self):
        """ Return a dict containing the PEP 566 meta-data. """

        if self.toml_error:
            # Provide a minimal default.
            return dict(name='unknown', version='0.1')

        metadata = dict()
        name = None
        version = None
        metadata_version = None
        requires_python = None

        for md_name, md_value in self.get_section('tool.sip.metadata', required=True).items():
            md_name = md_name.lower()

            # Extract specific string values.
            if md_name in ('name', 'version', 'metadata-version', 'requires-python'):
                if not isinstance(md_value, str):
                    raise PyProjectOptionException(md_name, "must be a string",
                            section_name='tool.sip')

                if md_name == 'name':
                    if not md_value.replace('-', '_').isidentifier():
                        raise PyProjectOptionException('name',
                                "'{0}' is an invalid project name".format(
                                        md_value),
                                section_name='tool.sip')

                    name = md_value
                elif md_name == 'version':
                    version = md_value
                elif md_name == 'metadata-version':
                    metadata_version = md_value
                elif md_name == 'requires-python':
                    requires_python = md_value
            else:
                # Any other value may be a string or a list of strings.
                value_list = md_value if isinstance(md_value, list) else [md_value]

                for value in value_list:
                    if not isinstance(value, str):
                        raise PyProjectOptionException(md_name,
                                "must be a string or a list of strings",
                                section_name='tool.sip')

            metadata[md_name] = md_value

        if name is None:
            raise PyProjectUndefinedOptionException('name',
                    section_name='tool.sip')

        if version is None:
            metadata['version'] = '0.1'

        if metadata_version is None:
            # Default to PEP 566.
            metadata['metadata-version'] = '2.1'

        if requires_python is None:
            # The minimal version of Python we support.
            metadata['requires-python'] = '>=3.{}'.format(
                    OLDEST_SUPPORTED_MINOR)

        return metadata

    def get_section(self, section_name, *, required=False):
        """ Return a sub-section with a dotted name. """

        if self.toml_error:
            return None

        section = self._pyproject

        for part in section_name.split('.'):
            try:
                section = section[part]
            except KeyError:
                if required:
                    raise PyProjectException(
                            "the '[{0}]' section is missing".format(
                                    section_name))

                return None

        if not self._is_section(section):
            raise PyProjectException(
                    "'{0}' is not a section".format(section_name))

        return section

    @staticmethod
    def _is_section(value):
        """ Returns True if a section value is itself a section. """

        return isinstance(value, (dict, list))
