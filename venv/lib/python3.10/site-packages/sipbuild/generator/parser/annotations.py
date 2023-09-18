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


class DottedName(str):
    """ Encapsulate a dotted name.  A dedicated type is used (rather than a
    str) because we need to be able to distinguish it from a quoted string when
    used as the value of an annotation.
    """

    pass


class InvalidAnnotation(Exception):
    """ An invalid annotation. """

    def __init__(self, name, message, use):
        """ Initialise the exception. """

        self._text = "{0} {1}".format(name, message)

        # The value to use for the annotation.
        self.use = use

    def __str__(self):
        """ Return the exception as a user friendly string. """

        return self._text


class RequiredAnnotation(InvalidAnnotation):
    """ A required annotation. """

    def __init__(name, use):
        """ Initialise the exception. """

        super().__init__(name, "requires a value", use=use)


def validate_annotation_value(pm, p, symbol, name, value):
    """ Return a valid value for the annotation or raise an InvalidAnnotation
    exception.
    """

    try:
        validator = _ANNOTATION_TYPES[name]
    except KeyError:
        raise InvalidAnnotation(name, "is not a known annotation", use=None)

    return validator(pm, p, symbol, name, value)


def bind(validator, **proto_kw):
    """ Return a function that when called with a validator function and
    prototype keyword arguments will itself return a function that will create
    an annotation-specific validator.
    """

    # This takes the prototype validator-specific keyword arguments and returns
    # a function that will itself create a validator with annotation-specific
    # arguments based on the prototypes.
    def proto_validator(**bound_kw):
        # Create the annotation specific keyword arguments by taking the
        # prototypes and updating them with the ones bound to the specific
        # annotation.
        kw = proto_kw.copy()
        kw.update(bound_kw)

        # This takes the name and value of the annotation and calls the
        # validator along with the annotation-specific keyword arguments.
        def bound_validator(pm, p, symbol, name, value):
            return validator(pm, p, symbol, name, value, **kw)

        return bound_validator

    return proto_validator


def validate_boolean(pm, p, symbol, name, value):
    """ Return a valid boolean value. """

    if value is None:
        return True

    raise InvalidAnnotation(name, "must not have a value", use=False)

boolean = bind(validate_boolean)


def validate_integer(pm, p, symbol, name, value, *, optional):
    """ Return a valid, possibly optional, integer. """

    if value is None:
        if optional:
            return None

        raise RequiredAnnotation(name, use=0)

    if not isinstance(value, int):
        raise InvalidAnnotation(name, "must be an integer", use=0)

    return value

integer = bind(validate_integer, optional=False)


def validate_name(pm, p, symbol, name, value, *, allow_dots, optional):
    """ Return a valid, possibly optional, possibly dotted name. """

    if value is None:
        if optional:
            return ''

        raise RequiredAnnotation(name, use='')

    if not isinstance(value, DottedName):
        raise InvalidAnnotation(name, "must be an unquoted name", use='')

    if '.' in value and not allow_dots:
        raise InvalidAnnotation(name, "cannot contain '.'", use='')

    return value

name = bind(validate_name, allow_dots=False, optional=False)


def validate_string(pm, p, symbol, name, value):
    """ Return a valid string value. """

    if not isinstance(value, str):
        raise InvalidAnnotation(name, "must be a quoted string", use='')

    # Handle any embedded selectors.
    for part in value.split(';'):
        if ':' not in part:
            return part.strip()

        selector, subvalue = part.split(':', maxsplit=1)
        if selector.startswith('!'):
            selector = selector[1:]
            inverted = True
        else:
            inverted = False

        if pm.evaluate_feature_or_platform(p, symbol, selector, inverted):
            return subvalue.strip()

    # No value was selected so ignore the annotation completely.
    return None

string = bind(validate_string)


def validate_string_list(pm, p, symbol, name, value):
    """ Return a valid string list value. """

    if not isinstance(value, str):
        raise InvalidAnnotation(name, "must be a quoted string", use=[])

    return value.split(' ')

string_list = bind(validate_string_list)


# The annotations and the type of their values.
_ANNOTATION_TYPES = {
    '__imatmul__':              boolean(),
    '__len__':                  boolean(),
    '__matmul__':               boolean(),
    'AbortOnException':         boolean(),
    'Abstract':                 boolean(),
    'AllowNone':                boolean(),
    'Array':                    boolean(),
    'ArraySize':                boolean(),
    'AutoGen':                  name(optional=True),
    'BaseType':                 name(),
    'Capsule':                  boolean(),
    'Constrained':              boolean(),
    'Deprecated':               boolean(),
    'Default':                  boolean(),
    'DelayDtor':                boolean(),
    'DisallowNone':             boolean(),
    'ExportDerived':            boolean(),
    'External':                 boolean(),
    'Encoding':                 string(),
    'Factory':                  boolean(),
    'FileExtension':            string(),
    'GetWrapper':               boolean(),
    'HoldGIL':                  boolean(),
    'In':                       boolean(),
    'KeepReference':            integer(optional=True),
    'KeywordArgs':              string(),
    'Metatype':                 name(allow_dots=True),
    'Mixin':                    boolean(),
    'NewThread':                boolean(),
    'NoArgParser':              boolean(),
    'NoAssignmentOperator':     boolean(),
    'NoCopy':                   boolean(),
    'NoCopyCtor':               boolean(),
    'NoDefaultCtor':            boolean(),
    'NoDefaultCtors':           boolean(),
    'NoDerived':                boolean(),
    'NoRaisesPyException':      boolean(),
    'NoRelease':                boolean(),
    'NoScope':                  boolean(),
    'NoSetter':                 boolean(),
    'NoTypeHint':               boolean(),
    'NoTypeName':               boolean(),
    'NoVirtualErrorHandler':    boolean(),
    'Numeric':                  boolean(),
    'Out':                      boolean(),
    'PostHook':                 name(),
    'PreHook':                  name(),
    'PyInt':                    boolean(),
    'PyName':                   name(),
    'PyQtFlags':                integer(),
    'PyQtFlagsEnums':           string_list(),
    'PyQtInterface':            string(),
    'PyQtNoQMetaObject':        boolean(),
    'RaisesPyException':        boolean(),
    'ReleaseGIL':               boolean(),
    'ResultSize':               boolean(),
    'ScopesStripped':           integer(),
    'Sequence':                 boolean(),
    'Supertype':                name(allow_dots=True),
    'Transfer':                 boolean(),
    'TransferBack':             boolean(),
    'TransferThis':             boolean(),
    'TypeHint':                 string(),
    'TypeHintIn':               string(),
    'TypeHintOut':              string(),
    'TypeHintValue':            string(),
    'VirtualErrorHandler':      name(),
}
