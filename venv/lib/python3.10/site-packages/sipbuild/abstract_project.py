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


from abc import ABC, abstractmethod
import importlib
import importlib.util
import os

from .exceptions import UserException
from .pyproject import PyProject, PyProjectOptionException


class AbstractProject(ABC):
    """ This specifies the API of a project. """

    @classmethod
    def bootstrap(cls, tool, tool_description='', arguments=None):
        """ Return an AbstractProject instance fully configured for a
        particular command line tool.
        """

        # Get the contents of the pyproject.toml file.
        pyproject = PyProject()

        # Get the name of the project factory.
        project_factory_name = None

        sip_section_name = 'tool.sip'
        value_name = 'project-factory'

        sip_section = pyproject.get_section(sip_section_name)
        if sip_section is not None:
            f_name = sip_section.get(value_name)
            if f_name is not None:
                if not isinstance(f_name, str):
                    raise PyProjectOptionException(value_name,
                            "should be a 'str' and not '{0}'".format(
                                    type(f_name).__name__),
                            section_name=sip_section_name)

                project_factory_name = f_name

        # See if there is a corresponding callable.
        if project_factory_name is None:
            default_factory_py = 'project.py'

            if os.path.isfile(default_factory_py):
                project_factory = cls.import_callable(default_factory_py, cls)
            else:
                # The default project factory.
                from .project import Project as project_factory
        else:
            project_factory = cls.import_callable(project_factory_name, cls)

        project = project_factory()

        if not isinstance(project, cls):
            raise UserException(
                    "The project factory did not return an AbstractProject "
                    "object")

        # We set this as an attribute rather than change the API of the ctor or
        # setup().
        project.arguments = arguments

        # Complete the configuration of the project.
        project.setup(pyproject, tool, tool_description)

        return project

    @abstractmethod
    def build(self):
        """ Build the project in-situ. """

    @abstractmethod
    def build_sdist(self, sdist_directory):
        """ Build an sdist for the project and return the name of the sdist
        file.
        """

    @abstractmethod
    def build_wheel(self, wheel_directory):
        """ Build a wheel for the project and return the name of the wheel
        file.
        """

    @staticmethod
    def import_callable(name, base_type):
        """ Import a callable from either a .py file or specified as
        module[:object].
        """

        # See if the name refers to a .py file.
        if name.endswith('.py'):
            name = name.replace('/', os.sep)

            module_name = name[:-3]
            object_name = None

            # Try and import the .py file.
            spec = importlib.util.spec_from_file_location(module_name, name)
            module = importlib.util.module_from_spec(spec)

            try:
                spec.loader.exec_module(module)
            except Exception as e:
                raise UserException("Unable to import '{0}'".format(name),
                        detail=str(e))
        else:
            # Extract the module and any object name.
            parts = name.split(':')
            if len(parts) > 2:
                raise UserException(
                        "The callable '{0}' must be specified as "
                        "'module[:name]'".format(name))

            module_name = parts[0]
            object_name = parts[1] if len(parts) == 2 else None

            # Try and import the module.
            try:
                module = importlib.import_module(module_name)
            except ImportError as e:
                raise UserException(
                        "Unable to import '{0}'".format(module_name),
                        detail=str(e))

        # Get the callable object from the module.
        if object_name is None:
            # Look for a class that is a sub-class of the base type.
            for obj in module.__dict__.values():
                if isinstance(obj, type):
                    if issubclass(obj, base_type):
                        # Make sure the type is defined in the module and not
                        # imported by it.
                        if obj.__module__ == module_name:
                            break
            else:
                raise UserException(
                        "'{0}' does not define a {1} sub-class".format(name,
                                base_type.__name__))
        else:
            # We have the name of the callable so just get it.
            obj = getattr(module, object_name)
            if obj is None:
                raise UserException(
                        "'{0}' module has no callable '{1}'".format(
                                module_name, object_name))

        return obj

    @abstractmethod
    def install(self):
        """ Install the project. """

    @abstractmethod
    def setup(self, pyproject, tool, tool_description):
        """ Complete the configuration of the project. """
