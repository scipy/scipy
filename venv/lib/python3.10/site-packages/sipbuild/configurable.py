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


from packaging.markers import Marker

from .exceptions import UserException
from .pyproject import PyProjectOptionException


class Configurable:
    """ A base class for an object that can be configured by a pyproject.toml,
    a build script or (possibly) the user via command line options.
    """

    def add_command_line_options(self, parser, tool, all_options, options=None):
        """ Add the object's command line options to an argument parser and
        update a map of Object instances and the configurables that they should
        be applied to.
        """

        if options is None:
            options = self.get_options()

        for option in options:
            # If it has already been set explicitly then the user cannot change
            # it.
            if getattr(self, option.name) is not None:
                continue

            # If there is no help then the user can never specify it.
            if option.help is None:
                continue

            # See if the option is applicable to this tool.
            if option.tools and tool not in option.tools:
                continue

            # See if the option has already been added.
            configurables = all_options.setdefault(option, list())
            configurables.append(self)

            if len(configurables) != 1:
                continue

            # Add the option according to its type.
            argument_name = option.user_name

            if option.inverted:
                argument_name = 'no-' + argument_name

            argument_name = '--' + argument_name

            if option.option_type is bool:
                parser.add_argument(argument_name, dest=option.dest,
                        action=('store_false' if option.inverted
                                else 'store_true'),
                        help=option.help)
            elif option.option_type is list:
                # Remove any plural.
                if argument_name.endswith('s'):
                    argument_name = argument_name[:-1]

                parser.add_argument(argument_name, dest=option.dest,
                        choices=option.choices, action='append',
                        help=option.help, metavar=option.metavar)
            else:
                parser.add_argument(argument_name, dest=option.dest,
                        type=option.option_type, choices=option.choices,
                        help=option.help, metavar=option.metavar)

    def apply_nonuser_defaults(self, tool):
        """ Set default values for each non-user configurable option that
        hasn't been set yet.
        """

        self._apply_defaults(tool, user=False)

    def apply_user_defaults(self, tool):
        """ Set default values for each user configurable option that hasn't
        been set yet.
        """

        self._apply_defaults(tool, user=True)

    def configure(self, pyproject, section_name, tool):
        """ Perform the initial configuration of an object. """

        section = pyproject.get_section(section_name)

        if section is not None:
            for name, value in section.items():
                # Find the corresponding option.
                for option in self.get_options():
                    if option.user_name == name:
                        break
                else:
                    raise PyProjectOptionException(name,
                            "is not a supported option",
                            section_name=section_name)

                # Check the type of the option.
                if not isinstance(value, option.option_type):
                    raise PyProjectOptionException(name,
                            "should be of type '{0}' and not '{1}'".format(
                                    option.option_type.__name__,
                                    type(value).__name__),
                            section_name=section_name)

                # Check the option hasn't already been initialised.
                if getattr(self, option.name) is not None:
                    raise PyProjectOptionException(name,
                            "has already been set in code and cannot be "
                            "changed",
                            section_name=section_name)

                # Evaluate any environment markers if the option supports them.
                if isinstance(value, list):
                    new_value = []

                    for v in value:
                        v = self._handle_marker(v, name, section_name)
                        if v is not None:
                            new_value.append(v)

                    value = new_value

                setattr(self, option.name, value)

        self.apply_nonuser_defaults(tool)

    def get_options(self):
        """ Return a list of configurable options. """

        return list()

    def initialise_options(self, kwargs):
        """ Initialise the options. """

        # Set the value for each option from the keyword arguments or undefined
        # if not specified.
        names = []

        for option in self.get_options():
            name = option.name
            names.append(name)

            setattr(self, name, kwargs.get(name))

        # Check that all keyword arguments are valid options.
        for kw in kwargs.keys():
            if kw not in names:
                raise UserException("'{0}' is not a valid option".format(kw))

    def verify_configuration(self, tool):
        """ Verify that the configuration is complete and consistent. """

        # This default implementation does nothing.

    def _apply_defaults(self, tool, user):
        """ Set default values for each user/non-user option that hasn't been
        set yet.
        """

        for option in self.get_options():
            if user and option.help is None:
                continue

            if not user and option.help is not None:
                continue

            value = getattr(self, option.name)
            if value is None:
                value = option.default
                if value is None:
                    # Provide a default default based on the type.
                    if option.option_type is list:
                        value = []
                    elif option.option_type is float:
                        value = 0.0
                    elif option.option_type is int:
                        value = 0
                    elif option.option_type is bool:
                        value = True if option.inverted else False
                    elif option.option_type is str:
                        value = ''
                else:
                    # Make a copy of the default in case it is mutable.
                    value = option.option_type(value)

                setattr(self, option.name, value)

    @staticmethod
    def _handle_marker(value, name, section_name):
        """ Handle any environment marker in a value.  The value is returned if
        a marker evaluates to True.  None is returned if a marker evaluates to
        False.
        """

        # Handle the trivial case of there being no marker.
        if ';' not in value:
            return value

        value, marker = value.split(';', maxsplit=1)

        try:
            satisfied = Marker(marker).evaluate()
        except:
            raise PyProjectOptionException(name,
                    "has an invalid marker '{0}'".format(marker),
                    section_name=section_name)

        return value if satisfied else None


class Option:
    """ Encapsulate a configuration option.  This defines and implements an
    attribute of a Configurable object.  The value of the attribute can be set
    either by __init__(), the pyproject.toml file and by the user using a
    command line argument (in that order).  Once the value is set it cannot be
    changed subsequently.  For example, if an attribute is set
    in pyproject.toml then the user will not then be able to modify it from the
    command line.  The value can only be changed from the command line if the
    Option object has help text specified.
    """

    # The tools that will build a set of bindings.
    BUILD_TOOLS = ('build', 'install', 'wheel')

    # All the valid tools.
    _ALL_TOOLS = BUILD_TOOLS + ('sdist', )

    # This is used to make sure each option (even if they are handling the same
    # attribute) has a unique 'dest'.
    option_nr = 0

    def __init__(self, name, *, option_type=str, choices=None, default=None,
            help=None, metavar=None, inverted=False, tools=None):
        """ Initialise the option. """

        self.name = name
        self.user_name = name.replace('_', '-')
        self.option_type = option_type
        self.default = default
        self.choices = choices
        self.help = help
        self.metavar = metavar
        self.inverted = inverted

        if tools is None:
            self.tools = self.BUILD_TOOLS
        else:
            for tool in tools:
                if tool not in self._ALL_TOOLS:
                    raise UserException(
                            "'{0}' option has an invalid tools '{1}'".format(
                                    name, tool))

            self.tools = tools

        self.dest = 'd' + str(type(self).option_nr)
        type(self).option_nr += 1
