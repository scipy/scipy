# Copyright (c) 2021, Riverbank Computing Limited
# All rights reserved.
#
# This copy of SIP is licensed for use under the terms of the SIP License
# Agreement.  See the file LICENSE for more details.
#
# This copy of SIP may also used under the terms of the GNU General Public
# License v2 or v3 as published by the Free Software Foundation which can be
# found in the files LICENSE-GPL2 and LICENSE-GPL3 included in this package.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ('AS IS"
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


from .specification import PySlot


# A 4-tuple of the name, the type, True if the slot requires %MethodCode and
# the number of arguments required (-1 if any number is valid).
_SLOT_DEFINITIONS = (
    ('__str__', PySlot.STR, True, 0),
    ('__int__', PySlot.INT, False, 0),
    ('__float__', PySlot.FLOAT, False, 0),
    ('__len__', PySlot.LEN, True, 0),
    ('__contains__', PySlot.CONTAINS, True, 1),
    ('__add__', PySlot.ADD, False, 1),
    ('__sub__', PySlot.SUB, False, 1),
    ('__mul__', PySlot.MUL, False, 1),
    ('__mod__', PySlot.MOD, False, 1),
    ('__floordiv__', PySlot.FLOORDIV, True, 1),
    ('__truediv__', PySlot.TRUEDIV, False, 1),
    ('__and__', PySlot.AND, False, 1),
    ('__or__', PySlot.OR, False, 1),
    ('__xor__', PySlot.XOR, False, 1),
    ('__lshift__', PySlot.LSHIFT, False, 1),
    ('__rshift__', PySlot.RSHIFT, False, 1),
    ('__iadd__', PySlot.IADD, False, 1),
    ('__isub__', PySlot.ISUB, False, 1),
    ('__imul__', PySlot.IMUL, False, 1),
    ('__imod__', PySlot.IMOD, False, 1),
    ('__ifloordiv__', PySlot.IFLOORDIV, True, 1),
    ('__itruediv__', PySlot.ITRUEDIV, False, 1),
    ('__iand__', PySlot.IAND, False, 1),
    ('__ior__', PySlot.IOR, False, 1),
    ('__ixor__', PySlot.IXOR, False, 1),
    ('__ilshift__', PySlot.ILSHIFT, False, 1),
    ('__irshift__', PySlot.IRSHIFT, False, 1),
    ('__invert__', PySlot.INVERT, False, 0),
    ('__call__', PySlot.CALL, False, -1),
    ('__getitem__', PySlot.GETITEM, False, 1),
    ('__setitem__', PySlot.SETITEM, True, 2),
    ('__delitem__', PySlot.DELITEM, True, 1),
    ('__lt__', PySlot.LT, False, 1),
    ('__le__', PySlot.LE, False, 1),
    ('__eq__', PySlot.EQ, False, 1),
    ('__ne__', PySlot.NE, False, 1),
    ('__gt__', PySlot.GT, False, 1),
    ('__ge__', PySlot.GE, False, 1),
    ('__bool__', PySlot.BOOL, True, 0),
    ('__neg__', PySlot.NEG, False, 0),
    ('__pos__', PySlot.POS, False, 0),
    ('__abs__', PySlot.ABS, True, 0),
    ('__repr__', PySlot.REPR, True, 0),
    ('__hash__', PySlot.HASH, True, 0),
    ('__index__', PySlot.INDEX, True, 0),
    ('__iter__', PySlot.ITER, True, 0),
    ('__next__', PySlot.NEXT, True, 0),
    ('__setattr__', PySlot.SETATTR, True, 2),
    ('__delattr__', PySlot.DELATTR, True, 1),
    ('__matmul__', PySlot.MATMUL, False, 1),
    ('__imatmul__', PySlot.IMATMUL, False, 1),
    ('__await__', PySlot.AWAIT, True, 0),
    ('__aiter__', PySlot.AITER, True, 0),
    ('__anext__', PySlot.ANEXT, True, 0),
)


# The mapping of slot type to slot name.
slot_type_name_map = {t: n for n, t, f, a in _SLOT_DEFINITIONS}


# The mapping of slot name to a 3-tuple of type, %MethodCode flag and argument
# count.
slot_name_detail_map = {n: (t, f, a) for n, t, f, a in _SLOT_DEFINITIONS}


def invalid_global_slot(slot):
    """ Return True if a slot cannot be specified as a global (ie. module
    level) slot.
    """

    if slot in (PySlot.NEG, PySlot.POS):
        return False

    if is_number_slot(slot):
        return False

    if is_inplace_number_slot(slot):
        return False

    if is_rich_compare_slot(slot):
        return False

    return True


def is_hash_return_slot(slot):
    """ Return True if a slot returns a Py_hash_t. """

    return slot is PySlot.HASH


_INPLACE_NUMBER_SLOTS = (PySlot.IADD, PySlot.ISUB, PySlot.IMUL, PySlot.IMOD,
        PySlot.IFLOORDIV, PySlot.ITRUEDIV, PySlot.IAND, PySlot.IOR,
        PySlot.IXOR, PySlot.ILSHIFT, PySlot.IRSHIFT, PySlot.IMATMUL)

def is_inplace_number_slot(slot):
    """ Return True if a slot is an inplace binary numeric slot. """

    return slot in _INPLACE_NUMBER_SLOTS


_INT_RETURN_SLOTS = (PySlot.BOOL, PySlot.CONTAINS)

def is_int_return_slot(slot):
    """ Return True if a slot returns an int. """

    return slot in _INT_RETURN_SLOTS


_NUMBER_SLOTS = (PySlot.ADD, PySlot.SUB, PySlot.MUL, PySlot.MOD,
        PySlot.FLOORDIV, PySlot.TRUEDIV, PySlot.AND, PySlot.OR, PySlot.XOR,
        PySlot.LSHIFT, PySlot.RSHIFT, PySlot.MATMUL)

def is_number_slot(slot):
    """ Return True if a slot is a binary numeric slot. """

    return slot in _NUMBER_SLOTS


_RICH_COMPARE_SLOTS = (PySlot.LT, PySlot.LE, PySlot.EQ, PySlot.NE, PySlot.GT,
        PySlot.GE)

def is_rich_compare_slot(slot):
    """ Return True if a slot is a rich comparision slot. """

    return slot in _RICH_COMPARE_SLOTS


def is_ssize_return_slot(slot):
    """ Return True if a slot returns a Py_ssize_t. """

    return slot is PySlot.LEN


_VOID_RETURN_SLOTS = (PySlot.SETITEM, PySlot.DELITEM, PySlot.SETATTR)

def is_void_return_slot(slot):
    """ Return True if a slot returns a void. """

    return slot in _VOID_RETURN_SLOTS


_ZERO_ARG_SLOTS = (PySlot.STR, PySlot.INT, PySlot.FLOAT, PySlot.INVERT,
        PySlot.NEG, PySlot.LEN, PySlot.BOOL, PySlot.POS, PySlot.ABS,
        PySlot.REPR, PySlot.HASH, PySlot.INDEX, PySlot.ITER, PySlot.NEXT,
        PySlot.AWAIT, PySlot.AITER, PySlot.ANEXT)

def is_zero_arg_slot(slot):
    """ Return True if a slot takes zero arguments. """

    return slot in _ZERO_ARG_SLOTS


# A map of slots and the names of their reflections.
_SLOT_REFLECTIONS = {
    PySlot.ADD: '__radd__',
    PySlot.SUB: '__rsub__',
    PySlot.MUL: '__rmul__',
    PySlot.MATMUL: '__rmatmul__',
    PySlot.TRUEDIV: '__rtruediv__',
    PySlot.FLOORDIV: '__rfloordiv__',
    PySlot.MOD: '__rmod__',
    PySlot.LSHIFT: '__rlshift__',
    PySlot.RSHIFT: '__rrshift__',
    PySlot.AND: '__rand__',
    PySlot.OR: '__ror__',
    PySlot.XOR: '__rxor__',
}

def reflected_slot(slot):
    """ Return the name of the reflected version of a slot or None if it
    doesn't have one.
    """

    return _SLOT_REFLECTIONS.get(slot)
