""" Attribute and method access on Python objects from C++.

    Note: std::cout type operations currently crash python...
          Not sure what is up with this...
"""
from __future__ import absolute_import, print_function

import scipy.weave as weave

#----------------------------------------------------------------------------
# get/set attribute and call methods example
#----------------------------------------------------------------------------


class Foo(object):
    def __init__(self):
        self.val = 1

    def inc(self,amount):
        self.val += amount
        return self.val
obj = Foo()
code = """
       py::tuple result(3);

       int i = obj.attr("val");
       result[0] = i;

       py::tuple args(1);
       args[0] = 2;
       i = obj.mcall("inc",args);
       result[1] = i;

       obj.set_attr("val",5);
       i = obj.attr("val");
       result[2] = i;

       return_val = result;
       """

print('initial, inc(2), set(5)/get:', weave.inline(code,['obj']))

#----------------------------------------------------------------------------
# indexing of values.
#----------------------------------------------------------------------------
from UserList import UserList
obj = UserList([1,[1,2],"hello"])
code = """
       int i;
       // find obj length and access each of its items
       //std::cout << "UserList items: ";
       //for(i = 0; i < obj.length(); i++)
       //    std::cout << obj[i].str() << " ";
       //std::cout << std::endl;
       // assign new values to each of its items
       for(i = 0; i < obj.length(); i++)
           obj[i] = "goodbye";
       """
weave.inline(code,['obj'])
print("obj with new values:", obj)
