# Xmethods for libstdc++.

# Copyright (C) 2014-2022 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import gdb
import gdb.xmethod
import re

matcher_name_prefix = 'libstdc++::'

def get_bool_type():
    return gdb.lookup_type('bool')

def get_std_size_type():
    return gdb.lookup_type('std::size_t')

class LibStdCxxXMethod(gdb.xmethod.XMethod):
    def __init__(self, name, worker_class):
        gdb.xmethod.XMethod.__init__(self, name)
        self.worker_class = worker_class

# Xmethods for std::array

class ArrayWorkerBase(gdb.xmethod.XMethodWorker):
    def __init__(self, val_type, size):
        self._val_type = val_type
        self._size = size

    def null_value(self):
        nullptr = gdb.parse_and_eval('(void *) 0')
        return nullptr.cast(self._val_type.pointer()).dereference()

class ArraySizeWorker(ArrayWorkerBase):
    def __init__(self, val_type, size):
        ArrayWorkerBase.__init__(self, val_type, size)

    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return get_std_size_type()

    def __call__(self, obj):
        return self._size

class ArrayEmptyWorker(ArrayWorkerBase):
    def __init__(self, val_type, size):
        ArrayWorkerBase.__init__(self, val_type, size)

    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return get_bool_type()

    def __call__(self, obj):
        return (int(self._size) == 0)

class ArrayFrontWorker(ArrayWorkerBase):
    def __init__(self, val_type, size):
        ArrayWorkerBase.__init__(self, val_type, size)

    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return self._val_type

    def __call__(self, obj):
        if int(self._size) > 0:
            return obj['_M_elems'][0]
        else:
            return self.null_value()

class ArrayBackWorker(ArrayWorkerBase):
    def __init__(self, val_type, size):
        ArrayWorkerBase.__init__(self, val_type, size)

    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return self._val_type

    def __call__(self, obj):
        if int(self._size) > 0:
            return obj['_M_elems'][self._size - 1]
        else:
            return self.null_value()

class ArrayAtWorker(ArrayWorkerBase):
    def __init__(self, val_type, size):
        ArrayWorkerBase.__init__(self, val_type, size)

    def get_arg_types(self):
        return get_std_size_type()

    def get_result_type(self, obj, index):
        return self._val_type

    def __call__(self, obj, index):
        if int(index) >= int(self._size):
            raise IndexError('Array index "%d" should not be >= %d.' %
                             ((int(index), self._size)))
        return obj['_M_elems'][index]

class ArraySubscriptWorker(ArrayWorkerBase):
    def __init__(self, val_type, size):
        ArrayWorkerBase.__init__(self, val_type, size)

    def get_arg_types(self):
        return get_std_size_type()

    def get_result_type(self, obj, index):
        return self._val_type

    def __call__(self, obj, index):
        if int(self._size) > 0:
            return obj['_M_elems'][index]
        else:
            return self.null_value()

class ArrayMethodsMatcher(gdb.xmethod.XMethodMatcher):
    def __init__(self):
        gdb.xmethod.XMethodMatcher.__init__(self,
                                            matcher_name_prefix + 'array')
        self._method_dict = {
            'size': LibStdCxxXMethod('size', ArraySizeWorker),
            'empty': LibStdCxxXMethod('empty', ArrayEmptyWorker),
            'front': LibStdCxxXMethod('front', ArrayFrontWorker),
            'back': LibStdCxxXMethod('back', ArrayBackWorker),
            'at': LibStdCxxXMethod('at', ArrayAtWorker),
            'operator[]': LibStdCxxXMethod('operator[]', ArraySubscriptWorker),
        }
        self.methods = [self._method_dict[m] for m in self._method_dict]

    def match(self, class_type, method_name):
        if not re.match('^std::(__\d+::)?array<.*>$', class_type.tag):
            return None
        method = self._method_dict.get(method_name)
        if method is None or not method.enabled:
            return None
        try:
            value_type = class_type.template_argument(0)
            size = class_type.template_argument(1)
        except:
            return None
        return method.worker_class(value_type, size)

# Xmethods for std::deque

class DequeWorkerBase(gdb.xmethod.XMethodWorker):
    def __init__(self, val_type):
        self._val_type = val_type
        self._bufsize = 512 // val_type.sizeof or 1

    def size(self, obj):
        first_node = obj['_M_impl']['_M_start']['_M_node']
        last_node = obj['_M_impl']['_M_finish']['_M_node']
        cur = obj['_M_impl']['_M_finish']['_M_cur']
        first = obj['_M_impl']['_M_finish']['_M_first']
        return (last_node - first_node) * self._bufsize + (cur - first)

    def index(self, obj, idx):
        first_node = obj['_M_impl']['_M_start']['_M_node']
        index_node = first_node + int(idx) // self._bufsize
        return index_node[0][idx % self._bufsize]

class DequeEmptyWorker(DequeWorkerBase):
    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return get_bool_type()

    def __call__(self, obj):
        return (obj['_M_impl']['_M_start']['_M_cur'] ==
                obj['_M_impl']['_M_finish']['_M_cur'])

class DequeSizeWorker(DequeWorkerBase):
    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return get_std_size_type()

    def __call__(self, obj):
        return self.size(obj)

class DequeFrontWorker(DequeWorkerBase):
    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return self._val_type

    def __call__(self, obj):
        return obj['_M_impl']['_M_start']['_M_cur'][0]

class DequeBackWorker(DequeWorkerBase):
    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return self._val_type

    def __call__(self, obj):
        if (obj['_M_impl']['_M_finish']['_M_cur'] ==
            obj['_M_impl']['_M_finish']['_M_first']):
            prev_node = obj['_M_impl']['_M_finish']['_M_node'] - 1
            return prev_node[0][self._bufsize - 1]
        else:
            return obj['_M_impl']['_M_finish']['_M_cur'][-1]

class DequeSubscriptWorker(DequeWorkerBase):
    def get_arg_types(self):
        return get_std_size_type()

    def get_result_type(self, obj, subscript):
        return self._val_type

    def __call__(self, obj, subscript):
        return self.index(obj, subscript)

class DequeAtWorker(DequeWorkerBase):
    def get_arg_types(self):
        return get_std_size_type()

    def get_result_type(self, obj, index):
        return self._val_type

    def __call__(self, obj, index):
        deque_size = int(self.size(obj))
        if int(index) >= deque_size:
            raise IndexError('Deque index "%d" should not be >= %d.' %
                             (int(index), deque_size))
        else:
           return self.index(obj, index)

class DequeMethodsMatcher(gdb.xmethod.XMethodMatcher):
    def __init__(self):
        gdb.xmethod.XMethodMatcher.__init__(self,
                                            matcher_name_prefix + 'deque')
        self._method_dict = {
            'empty': LibStdCxxXMethod('empty', DequeEmptyWorker),
            'size': LibStdCxxXMethod('size', DequeSizeWorker),
            'front': LibStdCxxXMethod('front', DequeFrontWorker),
            'back': LibStdCxxXMethod('back', DequeBackWorker),
            'operator[]': LibStdCxxXMethod('operator[]', DequeSubscriptWorker),
            'at': LibStdCxxXMethod('at', DequeAtWorker)
        }
        self.methods = [self._method_dict[m] for m in self._method_dict]

    def match(self, class_type, method_name):
        if not re.match('^std::(__\d+::)?deque<.*>$', class_type.tag):
            return None
        method = self._method_dict.get(method_name)
        if method is None or not method.enabled:
            return None
        return method.worker_class(class_type.template_argument(0))

# Xmethods for std::forward_list

class ForwardListWorkerBase(gdb.xmethod.XMethodMatcher):
    def __init__(self, val_type, node_type):
        self._val_type = val_type
        self._node_type = node_type

    def get_arg_types(self):
        return None

class ForwardListEmptyWorker(ForwardListWorkerBase):
    def get_result_type(self, obj):
        return get_bool_type()

    def __call__(self, obj):
        return obj['_M_impl']['_M_head']['_M_next'] == 0

class ForwardListFrontWorker(ForwardListWorkerBase):
    def get_result_type(self, obj):
        return self._val_type

    def __call__(self, obj):
        node = obj['_M_impl']['_M_head']['_M_next'].cast(self._node_type)
        val_address = node['_M_storage']['_M_storage'].address
        return val_address.cast(self._val_type.pointer()).dereference()

class ForwardListMethodsMatcher(gdb.xmethod.XMethodMatcher):
    def __init__(self):
        matcher_name = matcher_name_prefix + 'forward_list'
        gdb.xmethod.XMethodMatcher.__init__(self, matcher_name)
        self._method_dict = {
            'empty': LibStdCxxXMethod('empty', ForwardListEmptyWorker),
            'front': LibStdCxxXMethod('front', ForwardListFrontWorker)
        }
        self.methods = [self._method_dict[m] for m in self._method_dict]

    def match(self, class_type, method_name):
        if not re.match('^std::(__\d+::)?forward_list<.*>$', class_type.tag):
            return None
        method = self._method_dict.get(method_name)
        if method is None or not method.enabled:
            return None
        val_type = class_type.template_argument(0)
        node_type = gdb.lookup_type(str(class_type) + '::_Node').pointer()
        return method.worker_class(val_type, node_type)

# Xmethods for std::list

class ListWorkerBase(gdb.xmethod.XMethodWorker):
    def __init__(self, val_type, node_type):
        self._val_type = val_type
        self._node_type = node_type

    def get_arg_types(self):
        return None

    def get_value_from_node(self, node):
        node = node.dereference()
        if node.type.fields()[1].name == '_M_data':
            # C++03 implementation, node contains the value as a member
            return node['_M_data']
        # C++11 implementation, node stores value in __aligned_membuf
        addr = node['_M_storage'].address
        return addr.cast(self._val_type.pointer()).dereference()

class ListEmptyWorker(ListWorkerBase):
    def get_result_type(self, obj):
        return get_bool_type()

    def __call__(self, obj):
        base_node = obj['_M_impl']['_M_node']
        if base_node['_M_next'] == base_node.address:
            return True
        else:
            return False

class ListSizeWorker(ListWorkerBase):
    def get_result_type(self, obj):
        return get_std_size_type()

    def __call__(self, obj):
        begin_node = obj['_M_impl']['_M_node']['_M_next']
        end_node = obj['_M_impl']['_M_node'].address
        size = 0
        while begin_node != end_node:
            begin_node = begin_node['_M_next']
            size += 1
        return size

class ListFrontWorker(ListWorkerBase):
    def get_result_type(self, obj):
        return self._val_type

    def __call__(self, obj):
        node = obj['_M_impl']['_M_node']['_M_next'].cast(self._node_type)
        return self.get_value_from_node(node)

class ListBackWorker(ListWorkerBase):
    def get_result_type(self, obj):
        return self._val_type

    def __call__(self, obj):
        prev_node = obj['_M_impl']['_M_node']['_M_prev'].cast(self._node_type)
        return self.get_value_from_node(prev_node)

class ListMethodsMatcher(gdb.xmethod.XMethodMatcher):
    def __init__(self):
        gdb.xmethod.XMethodMatcher.__init__(self,
                                            matcher_name_prefix + 'list')
        self._method_dict = {
            'empty': LibStdCxxXMethod('empty', ListEmptyWorker),
            'size': LibStdCxxXMethod('size', ListSizeWorker),
            'front': LibStdCxxXMethod('front', ListFrontWorker),
            'back': LibStdCxxXMethod('back', ListBackWorker)
        }
        self.methods = [self._method_dict[m] for m in self._method_dict]

    def match(self, class_type, method_name):
        if not re.match('^std::(__\d+::)?(__cxx11::)?list<.*>$', class_type.tag):
            return None
        method = self._method_dict.get(method_name)
        if method is None or not method.enabled:
            return None
        val_type = class_type.template_argument(0)
        node_type = gdb.lookup_type(str(class_type) + '::_Node').pointer()
        return method.worker_class(val_type, node_type)

# Xmethods for std::vector

class VectorWorkerBase(gdb.xmethod.XMethodWorker):
    def __init__(self, val_type):
        self._val_type = val_type

    def size(self, obj):
        if self._val_type.code == gdb.TYPE_CODE_BOOL:
            start = obj['_M_impl']['_M_start']['_M_p']
            finish = obj['_M_impl']['_M_finish']['_M_p']
            finish_offset = obj['_M_impl']['_M_finish']['_M_offset']
            bit_size = start.dereference().type.sizeof * 8
            return (finish - start) * bit_size + finish_offset
        else:
            return obj['_M_impl']['_M_finish'] - obj['_M_impl']['_M_start']

    def get(self, obj, index):
        if self._val_type.code == gdb.TYPE_CODE_BOOL:
            start = obj['_M_impl']['_M_start']['_M_p']
            bit_size = start.dereference().type.sizeof * 8
            valp = start + index // bit_size
            offset = index % bit_size
            return (valp.dereference() & (1 << offset)) > 0
        else:
            return obj['_M_impl']['_M_start'][index]

class VectorEmptyWorker(VectorWorkerBase):
    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return get_bool_type()

    def __call__(self, obj):
        return int(self.size(obj)) == 0

class VectorSizeWorker(VectorWorkerBase):
    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return get_std_size_type()

    def __call__(self, obj):
        return self.size(obj)

class VectorFrontWorker(VectorWorkerBase):
    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return self._val_type

    def __call__(self, obj):
        return self.get(obj, 0)

class VectorBackWorker(VectorWorkerBase):
    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return self._val_type

    def __call__(self, obj):
        return self.get(obj, int(self.size(obj)) - 1)

class VectorAtWorker(VectorWorkerBase):
    def get_arg_types(self):
        return get_std_size_type()

    def get_result_type(self, obj, index):
        return self._val_type

    def __call__(self, obj, index):
        size = int(self.size(obj))
        if int(index) >= size:
            raise IndexError('Vector index "%d" should not be >= %d.' %
                             ((int(index), size)))
        return self.get(obj, int(index))

class VectorSubscriptWorker(VectorWorkerBase):
    def get_arg_types(self):
        return get_std_size_type()

    def get_result_type(self, obj, subscript):
        return self._val_type

    def __call__(self, obj, subscript):
        return self.get(obj, int(subscript))

class VectorMethodsMatcher(gdb.xmethod.XMethodMatcher):
    def __init__(self):
        gdb.xmethod.XMethodMatcher.__init__(self,
                                            matcher_name_prefix + 'vector')
        self._method_dict = {
            'size': LibStdCxxXMethod('size', VectorSizeWorker),
            'empty': LibStdCxxXMethod('empty', VectorEmptyWorker),
            'front': LibStdCxxXMethod('front', VectorFrontWorker),
            'back': LibStdCxxXMethod('back', VectorBackWorker),
            'at': LibStdCxxXMethod('at', VectorAtWorker),
            'operator[]': LibStdCxxXMethod('operator[]',
                                           VectorSubscriptWorker),
        }
        self.methods = [self._method_dict[m] for m in self._method_dict]

    def match(self, class_type, method_name):
        if not re.match('^std::(__\d+::)?vector<.*>$', class_type.tag):
            return None
        method = self._method_dict.get(method_name)
        if method is None or not method.enabled:
            return None
        return method.worker_class(class_type.template_argument(0))

# Xmethods for associative containers

class AssociativeContainerWorkerBase(gdb.xmethod.XMethodWorker):
    def __init__(self, unordered):
        self._unordered = unordered

    def node_count(self, obj):
        if self._unordered:
            return obj['_M_h']['_M_element_count']
        else:
            return obj['_M_t']['_M_impl']['_M_node_count']

    def get_arg_types(self):
        return None

class AssociativeContainerEmptyWorker(AssociativeContainerWorkerBase):
    def get_result_type(self, obj):
        return get_bool_type()

    def __call__(self, obj):
        return int(self.node_count(obj)) == 0

class AssociativeContainerSizeWorker(AssociativeContainerWorkerBase):
    def get_result_type(self, obj):
        return get_std_size_type()

    def __call__(self, obj):
        return self.node_count(obj)

class AssociativeContainerMethodsMatcher(gdb.xmethod.XMethodMatcher):
    def __init__(self, name):
        gdb.xmethod.XMethodMatcher.__init__(self,
                                            matcher_name_prefix + name)
        self._name = name
        self._method_dict = {
            'size': LibStdCxxXMethod('size', AssociativeContainerSizeWorker),
            'empty': LibStdCxxXMethod('empty',
                                      AssociativeContainerEmptyWorker),
        }
        self.methods = [self._method_dict[m] for m in self._method_dict]

    def match(self, class_type, method_name):
        if not re.match('^std::(__\d+::)?%s<.*>$' % self._name, class_type.tag):
            return None
        method = self._method_dict.get(method_name)
        if method is None or not method.enabled:
            return None
        unordered = 'unordered' in self._name
        return method.worker_class(unordered)

# Xmethods for std::unique_ptr

class UniquePtrGetWorker(gdb.xmethod.XMethodWorker):
    "Implements std::unique_ptr<T>::get() and std::unique_ptr<T>::operator->()"

    def __init__(self, elem_type):
        self._is_array = elem_type.code == gdb.TYPE_CODE_ARRAY
        if self._is_array:
            self._elem_type = elem_type.target()
        else:
            self._elem_type = elem_type

    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return self._elem_type.pointer()

    def _supports(self, method_name):
        "operator-> is not supported for unique_ptr<T[]>"
        return method_name == 'get' or not self._is_array

    def __call__(self, obj):
        impl_type = obj.dereference().type.fields()[0].type.tag
        # Check for new implementations first:
        if re.match('^std::(__\d+::)?__uniq_ptr_(data|impl)<.*>$', impl_type):
            tuple_member = obj['_M_t']['_M_t']
        elif re.match('^std::(__\d+::)?tuple<.*>$', impl_type):
            tuple_member = obj['_M_t']
        else:
            return None
        tuple_impl_type = tuple_member.type.fields()[0].type # _Tuple_impl
        tuple_head_type = tuple_impl_type.fields()[1].type   # _Head_base
        head_field = tuple_head_type.fields()[0]
        if head_field.name == '_M_head_impl':
            return tuple_member.cast(tuple_head_type)['_M_head_impl']
        elif head_field.is_base_class:
            return tuple_member.cast(head_field.type)
        else:
            return None

class UniquePtrDerefWorker(UniquePtrGetWorker):
    "Implements std::unique_ptr<T>::operator*()"

    def __init__(self, elem_type):
        UniquePtrGetWorker.__init__(self, elem_type)

    def get_result_type(self, obj):
        return self._elem_type

    def _supports(self, method_name):
        "operator* is not supported for unique_ptr<T[]>"
        return not self._is_array

    def __call__(self, obj):
        return UniquePtrGetWorker.__call__(self, obj).dereference()

class UniquePtrSubscriptWorker(UniquePtrGetWorker):
    "Implements std::unique_ptr<T>::operator[](size_t)"

    def __init__(self, elem_type):
        UniquePtrGetWorker.__init__(self, elem_type)

    def get_arg_types(self):
        return get_std_size_type()

    def get_result_type(self, obj, index):
        return self._elem_type

    def _supports(self, method_name):
        "operator[] is only supported for unique_ptr<T[]>"
        return self._is_array

    def __call__(self, obj, index):
        return UniquePtrGetWorker.__call__(self, obj)[index]

class UniquePtrMethodsMatcher(gdb.xmethod.XMethodMatcher):
    def __init__(self):
        gdb.xmethod.XMethodMatcher.__init__(self,
                                            matcher_name_prefix + 'unique_ptr')
        self._method_dict = {
            'get': LibStdCxxXMethod('get', UniquePtrGetWorker),
            'operator->': LibStdCxxXMethod('operator->', UniquePtrGetWorker),
            'operator*': LibStdCxxXMethod('operator*', UniquePtrDerefWorker),
            'operator[]': LibStdCxxXMethod('operator[]', UniquePtrSubscriptWorker),
        }
        self.methods = [self._method_dict[m] for m in self._method_dict]

    def match(self, class_type, method_name):
        if not re.match('^std::(__\d+::)?unique_ptr<.*>$', class_type.tag):
            return None
        method = self._method_dict.get(method_name)
        if method is None or not method.enabled:
            return None
        worker = method.worker_class(class_type.template_argument(0))
        if worker._supports(method_name):
            return worker
        return None

# Xmethods for std::shared_ptr

class SharedPtrGetWorker(gdb.xmethod.XMethodWorker):
    "Implements std::shared_ptr<T>::get() and std::shared_ptr<T>::operator->()"

    def __init__(self, elem_type):
        self._is_array = elem_type.code == gdb.TYPE_CODE_ARRAY
        if self._is_array:
            self._elem_type = elem_type.target()
        else:
            self._elem_type = elem_type

    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return self._elem_type.pointer()

    def _supports(self, method_name):
        "operator-> is not supported for shared_ptr<T[]>"
        return method_name == 'get' or not self._is_array

    def __call__(self, obj):
        return obj['_M_ptr']

class SharedPtrDerefWorker(SharedPtrGetWorker):
    "Implements std::shared_ptr<T>::operator*()"

    def __init__(self, elem_type):
        SharedPtrGetWorker.__init__(self, elem_type)

    def get_result_type(self, obj):
        return self._elem_type

    def _supports(self, method_name):
        "operator* is not supported for shared_ptr<T[]>"
        return not self._is_array

    def __call__(self, obj):
        return SharedPtrGetWorker.__call__(self, obj).dereference()

class SharedPtrSubscriptWorker(SharedPtrGetWorker):
    "Implements std::shared_ptr<T>::operator[](size_t)"

    def __init__(self, elem_type):
        SharedPtrGetWorker.__init__(self, elem_type)

    def get_arg_types(self):
        return get_std_size_type()

    def get_result_type(self, obj, index):
        return self._elem_type

    def _supports(self, method_name):
        "operator[] is only supported for shared_ptr<T[]>"
        return self._is_array

    def __call__(self, obj, index):
        # Check bounds if _elem_type is an array of known bound
        m = re.match('.*\[(\d+)]$', str(self._elem_type))
        if m and index >= int(m.group(1)):
            raise IndexError('shared_ptr<%s> index "%d" should not be >= %d.' %
                             (self._elem_type, int(index), int(m.group(1))))
        return SharedPtrGetWorker.__call__(self, obj)[index]

class SharedPtrUseCountWorker(gdb.xmethod.XMethodWorker):
    "Implements std::shared_ptr<T>::use_count()"

    def __init__(self, elem_type):
        pass

    def get_arg_types(self):
        return None

    def get_result_type(self, obj):
        return gdb.lookup_type('long')

    def _supports(self, method_name):
        return True

    def __call__(self, obj):
        refcounts = obj['_M_refcount']['_M_pi']
        return refcounts['_M_use_count'] if refcounts else 0

class SharedPtrUniqueWorker(SharedPtrUseCountWorker):
    "Implements std::shared_ptr<T>::unique()"

    def __init__(self, elem_type):
        SharedPtrUseCountWorker.__init__(self, elem_type)

    def get_result_type(self, obj):
        return gdb.lookup_type('bool')

    def __call__(self, obj):
        return SharedPtrUseCountWorker.__call__(self, obj) == 1

class SharedPtrMethodsMatcher(gdb.xmethod.XMethodMatcher):
    def __init__(self):
        gdb.xmethod.XMethodMatcher.__init__(self,
                                            matcher_name_prefix + 'shared_ptr')
        self._method_dict = {
            'get': LibStdCxxXMethod('get', SharedPtrGetWorker),
            'operator->': LibStdCxxXMethod('operator->', SharedPtrGetWorker),
            'operator*': LibStdCxxXMethod('operator*', SharedPtrDerefWorker),
            'operator[]': LibStdCxxXMethod('operator[]', SharedPtrSubscriptWorker),
            'use_count': LibStdCxxXMethod('use_count', SharedPtrUseCountWorker),
            'unique': LibStdCxxXMethod('unique', SharedPtrUniqueWorker),
        }
        self.methods = [self._method_dict[m] for m in self._method_dict]

    def match(self, class_type, method_name):
        if not re.match('^std::(__\d+::)?shared_ptr<.*>$', class_type.tag):
            return None
        method = self._method_dict.get(method_name)
        if method is None or not method.enabled:
            return None
        worker = method.worker_class(class_type.template_argument(0))
        if worker._supports(method_name):
            return worker
        return None

def register_libstdcxx_xmethods(locus):
    gdb.xmethod.register_xmethod_matcher(locus, ArrayMethodsMatcher())
    gdb.xmethod.register_xmethod_matcher(locus, ForwardListMethodsMatcher())
    gdb.xmethod.register_xmethod_matcher(locus, DequeMethodsMatcher())
    gdb.xmethod.register_xmethod_matcher(locus, ListMethodsMatcher())
    gdb.xmethod.register_xmethod_matcher(locus, VectorMethodsMatcher())
    gdb.xmethod.register_xmethod_matcher(
        locus, AssociativeContainerMethodsMatcher('set'))
    gdb.xmethod.register_xmethod_matcher(
        locus, AssociativeContainerMethodsMatcher('map'))
    gdb.xmethod.register_xmethod_matcher(
        locus, AssociativeContainerMethodsMatcher('multiset'))
    gdb.xmethod.register_xmethod_matcher(
        locus, AssociativeContainerMethodsMatcher('multimap'))
    gdb.xmethod.register_xmethod_matcher(
        locus, AssociativeContainerMethodsMatcher('unordered_set'))
    gdb.xmethod.register_xmethod_matcher(
        locus, AssociativeContainerMethodsMatcher('unordered_map'))
    gdb.xmethod.register_xmethod_matcher(
        locus, AssociativeContainerMethodsMatcher('unordered_multiset'))
    gdb.xmethod.register_xmethod_matcher(
        locus, AssociativeContainerMethodsMatcher('unordered_multimap'))
    gdb.xmethod.register_xmethod_matcher(locus, UniquePtrMethodsMatcher())
    gdb.xmethod.register_xmethod_matcher(locus, SharedPtrMethodsMatcher())
