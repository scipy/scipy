""" Test refcounting and behavior of SCXX.
"""
import time
import os,sys

from numpy.testing import *
set_package_path()
from weave import inline_tools
restore_path()

class TestObjectConstruct(NumpyTestCase):
    #------------------------------------------------------------------------
    # Check that construction from basic types is allowed and have correct
    # reference counts
    #------------------------------------------------------------------------
    def check_int(self,level=5):
        # strange int value used to try and make sure refcount is 2.
        code = """
               py::object val = 1001;
               return_val = val;
               """
        res = inline_tools.inline(code)
        assert_equal(sys.getrefcount(res),2)
        assert_equal(res,1001)
    def check_float(self,level=5):
        code = """
               py::object val = (float)1.0;
               return_val = val;
               """
        res = inline_tools.inline(code)
        assert_equal(sys.getrefcount(res),2)
        assert_equal(res,1.0)
    def check_double(self,level=5):
        code = """
               py::object val = 1.0;
               return_val = val;
               """
        res = inline_tools.inline(code)
        assert_equal(sys.getrefcount(res),2)
        assert_equal(res,1.0)
    def check_complex(self,level=5):
        code = """
               std::complex<double> num = std::complex<double>(1.0,1.0);
               py::object val = num;
               return_val = val;
               """
        res = inline_tools.inline(code)
        assert_equal(sys.getrefcount(res),2)
        assert_equal(res,1.0+1.0j)
    def check_string(self,level=5):
        code = """
               py::object val = "hello";
               return_val = val;
               """
        res = inline_tools.inline(code)
        assert_equal(sys.getrefcount(res),2)
        assert_equal(res,"hello")

    def check_std_string(self,level=5):
        code = """
               std::string s = std::string("hello");
               py::object val = s;
               return_val = val;
               """
        res = inline_tools.inline(code)
        assert_equal(sys.getrefcount(res),2)
        assert_equal(res,"hello")

class TestObjectPrint(NumpyTestCase):
    #------------------------------------------------------------------------
    # Check the object print protocol.
    #------------------------------------------------------------------------
    def check_stdout(self,level=5):
        code = """
               py::object val = "how now brown cow";
               val.print(stdout);
               """
        res = inline_tools.inline(code)
        # visual check on this one.
    def check_stringio(self,level=5):
        import cStringIO
        file_imposter = cStringIO.StringIO()
        code = """
               py::object val = "how now brown cow";
               val.print(file_imposter);
               """
        res = inline_tools.inline(code,['file_imposter'])
        print file_imposter.getvalue()
        assert_equal(file_imposter.getvalue(),"'how now brown cow'")

##    def check_failure(self,level=5):
##        code = """
##               FILE* file = 0;
##               py::object val = "how now brown cow";
##               val.print(file);
##               """
##        try:
##            res = inline_tools.inline(code)
##        except:
##            # error was supposed to occur.
##            pass


class TestObjectCast(NumpyTestCase):
    def check_int_cast(self,level=5):
        code = """
               py::object val = 1;
               int raw_val = val;
               """
        inline_tools.inline(code)
    def check_double_cast(self,level=5):
        code = """
               py::object val = 1.0;
               double raw_val = val;
               """
        inline_tools.inline(code)
    def check_float_cast(self,level=5):
        code = """
               py::object val = 1.0;
               float raw_val = val;
               """
        inline_tools.inline(code)
    def check_complex_cast(self,level=5):
        code = """
               std::complex<double> num = std::complex<double>(1.0,1.0);
               py::object val = num;
               std::complex<double> raw_val = val;
               """
        inline_tools.inline(code)
    def check_string_cast(self,level=5):
        code = """
               py::object val = "hello";
               std::string raw_val = val;
               """
        inline_tools.inline(code)

# test class used for testing python class access from C++.
class foo:
    def bar(self):
        return "bar results"
    def bar2(self,val1,val2):
        return val1, val2
    def bar3(self,val1,val2,val3=1):
        return val1, val2, val3

class str_obj:
    def __str__(self):
        return "b"

class TestObjectHasattr(NumpyTestCase):
    def check_string(self,level=5):
        a = foo()
        a.b = 12345
        code = """
               return_val = a.hasattr("b");
               """
        res = inline_tools.inline(code,['a'])
        assert res
    def check_std_string(self,level=5):
        a = foo()
        a.b = 12345
        attr_name = "b"
        code = """
               return_val = a.hasattr(attr_name);
               """
        res = inline_tools.inline(code,['a','attr_name'])
        assert res
    def check_string_fail(self,level=5):
        a = foo()
        a.b = 12345
        code = """
               return_val = a.hasattr("c");
               """
        res = inline_tools.inline(code,['a'])
        assert not res
    def check_inline(self,level=5):
        """ THIS NEEDS TO MOVE TO THE INLINE TEST SUITE
        """
        a = foo()
        a.b = 12345
        code = """
               throw_error(PyExc_AttributeError,"bummer");
               """
        try:
            before = sys.getrefcount(a)
            res = inline_tools.inline(code,['a'])
        except AttributeError:
            after = sys.getrefcount(a)
            try:
                res = inline_tools.inline(code,['a'])
            except:
                after2 = sys.getrefcount(a)
            print "after and after2 should be equal in the following"
            print 'before, after, after2:', before, after, after2
            pass

    def check_func(self,level=5):
        a = foo()
        a.b = 12345
        code = """
               return_val = a.hasattr("bar");
               """
        res = inline_tools.inline(code,['a'])
        assert res

class TestObjectAttr(NumpyTestCase):

    def generic_attr(self,code,args=['a']):
        a = foo()
        a.b = 12345

        before = sys.getrefcount(a.b)
        res = inline_tools.inline(code,args)
        assert_equal(res,a.b)
        del res
        after = sys.getrefcount(a.b)
        assert_equal(after,before)

    def check_char(self,level=5):
        self.generic_attr('return_val = a.attr("b");')

    def check_char_fail(self,level=5):
        try:
            self.generic_attr('return_val = a.attr("c");')
        except AttributeError:
            pass

    def check_string(self,level=5):
        self.generic_attr('return_val = a.attr(std::string("b"));')

    def check_string_fail(self,level=5):
        try:
            self.generic_attr('return_val = a.attr(std::string("c"));')
        except AttributeError:
            pass

    def check_obj(self,level=5):
        code = """
               py::object name = "b";
               return_val = a.attr(name);
               """
        self.generic_attr(code,['a'])

    def check_obj_fail(self,level=5):
        try:
            code = """
                   py::object name = "c";
                   return_val = a.attr(name);
                   """
            self.generic_attr(code,['a'])
        except AttributeError:
            pass

    def check_attr_call(self,level=5):
        a = foo()
        res = inline_tools.inline('return_val = a.attr("bar").call();',['a'])
        first = sys.getrefcount(res)
        del res
        res = inline_tools.inline('return_val = a.attr("bar").call();',['a'])
        second = sys.getrefcount(res)
        assert_equal(res,"bar results")
        assert_equal(first,second)

class TestObjectSetAttr(NumpyTestCase):

    def generic_existing(self, code, desired):
        args = ['a']
        a = foo()
        a.b = 12345
        res = inline_tools.inline(code,args)
        assert_equal(a.b,desired)

    def generic_new(self, code, desired):
        args = ['a']
        a = foo()
        res = inline_tools.inline(code,args)
        assert_equal(a.b,desired)

    def check_existing_char(self,level=5):
        self.generic_existing('a.set_attr("b","hello");',"hello")
    def check_new_char(self,level=5):
        self.generic_new('a.set_attr("b","hello");',"hello")
    def check_existing_string(self,level=5):
        self.generic_existing('a.set_attr("b",std::string("hello"));',"hello")
    def check_new_string(self,level=5):
        self.generic_new('a.set_attr("b",std::string("hello"));',"hello")
    def check_existing_object(self,level=5):
        code = """
               py::object obj = "hello";
               a.set_attr("b",obj);
               """
        self.generic_existing(code,"hello")
    def check_new_object(self,level=5):
        code = """
               py::object obj = "hello";
               a.set_attr("b",obj);
               """
        self.generic_new(code,"hello")
    def check_new_fail(self,level=5):
        try:
            code = """
                   py::object obj = 1;
                   a.set_attr(obj,"hello");
                   """
            self.generic_new(code,"hello")
        except:
            pass

    def check_existing_int(self,level=5):
        self.generic_existing('a.set_attr("b",1);',1)
    def check_existing_double(self,level=5):
        self.generic_existing('a.set_attr("b",1.0);',1.0)
    def check_existing_complex(self,level=5):
        code = """
               std::complex<double> obj = std::complex<double>(1,1);
               a.set_attr("b",obj);
               """
        self.generic_existing(code,1+1j)
    def check_existing_char1(self,level=5):
        self.generic_existing('a.set_attr("b","hello");',"hello")
    def check_existing_string1(self,level=5):
        code = """
               std::string obj = std::string("hello");
               a.set_attr("b",obj);
               """
        self.generic_existing(code,"hello")

class TestObjectDel(NumpyTestCase):
    def generic(self, code):
        args = ['a']
        a = foo()
        a.b = 12345
        res = inline_tools.inline(code,args)
        assert not hasattr(a,"b")

    def check_char(self,level=5):
        self.generic('a.del("b");')
    def check_string(self,level=5):
        code = """
               std::string name = std::string("b");
               a.del(name);
               """
        self.generic(code)
    def check_object(self,level=5):
        code = """
               py::object name = py::object("b");
               a.del(name);
               """
        self.generic(code)

class TestObjectCmp(NumpyTestCase):
    def check_equal(self,level=5):
        a,b = 1,1
        res = inline_tools.inline('return_val = (a == b);',['a','b'])
        assert_equal(res,(a == b))
    def check_equal_objects(self,level=5):
        class foo:
            def __init__(self,x):
                self.x = x
            def __cmp__(self,other):
                return cmp(self.x,other.x)
        a,b = foo(1),foo(2)
        res = inline_tools.inline('return_val = (a == b);',['a','b'])
        assert_equal(res,(a == b))
    def check_lt(self,level=5):
        a,b = 1,2
        res = inline_tools.inline('return_val = (a < b);',['a','b'])
        assert_equal(res,(a < b))
    def check_gt(self,level=5):
        a,b = 1,2
        res = inline_tools.inline('return_val = (a > b);',['a','b'])
        assert_equal(res,(a > b))
    def check_gte(self,level=5):
        a,b = 1,2
        res = inline_tools.inline('return_val = (a >= b);',['a','b'])
        assert_equal(res,(a >= b))
    def check_lte(self,level=5):
        a,b = 1,2
        res = inline_tools.inline('return_val = (a <= b);',['a','b'])
        assert_equal(res,(a <= b))
    def check_not_equal(self,level=5):
        a,b = 1,2
        res = inline_tools.inline('return_val = (a != b);',['a','b'])
        assert_equal(res,(a != b))
    def check_int(self,level=5):
        a = 1
        res = inline_tools.inline('return_val = (a == 1);',['a'])
        assert_equal(res,(a == 1))
    def check_int2(self,level=5):
        a = 1
        res = inline_tools.inline('return_val = (1 == a);',['a'])
        assert_equal(res,(a == 1))
    def check_unsigned_long(self,level=5):
        a = 1
        res = inline_tools.inline('return_val = (a == (unsigned long)1);',['a'])
        assert_equal(res,(a == 1))
    def check_double(self,level=5):
        a = 1
        res = inline_tools.inline('return_val = (a == 1.0);',['a'])
        assert_equal(res,(a == 1.0))
    def check_char(self,level=5):
        a = "hello"
        res = inline_tools.inline('return_val = (a == "hello");',['a'])
        assert_equal(res,(a == "hello"))
    def check_std_string(self,level=5):
        a = "hello"
        code = """
               std::string hello = std::string("hello");
               return_val = (a == hello);
               """
        res = inline_tools.inline(code,['a'])
        assert_equal(res,(a == "hello"))

class TestObjectRepr(NumpyTestCase):
    def check_repr(self,level=5):
        class foo:
            def __str__(self):
                return "str return"
            def __repr__(self):
                return "repr return"
        a = foo()
        res = inline_tools.inline('return_val = a.repr();',['a'])
        first = sys.getrefcount(res)
        del res
        res = inline_tools.inline('return_val = a.repr();',['a'])
        second = sys.getrefcount(res)
        assert_equal(first,second)
        assert_equal(res,"repr return")

class TestObjectStr(NumpyTestCase):
    def check_str(self,level=5):
        class foo:
            def __str__(self):
                return "str return"
            def __repr__(self):
                return "repr return"
        a = foo()
        res = inline_tools.inline('return_val = a.str();',['a'])
        first = sys.getrefcount(res)
        del res
        res = inline_tools.inline('return_val = a.str();',['a'])
        second = sys.getrefcount(res)
        assert_equal(first,second)
        print res
        assert_equal(res,"str return")

class TestObjectUnicode(NumpyTestCase):
    # This ain't going to win awards for test of the year...
    def check_unicode(self,level=5):
        class foo:
            def __repr__(self):
                return "repr return"
            def __str__(self):
                return "unicode"
        a= foo()
        res = inline_tools.inline('return_val = a.unicode();',['a'])
        first = sys.getrefcount(res)
        del res
        res = inline_tools.inline('return_val = a.unicode();',['a'])
        second = sys.getrefcount(res)
        assert_equal(first,second)
        assert_equal(res,"unicode")

class TestObjectIsCallable(NumpyTestCase):
    def check_true(self,level=5):
        class foo:
            def __call__(self):
                return 0
        a= foo()
        res = inline_tools.inline('return_val = a.is_callable();',['a'])
        assert res
    def check_false(self,level=5):
        class foo:
            pass
        a= foo()
        res = inline_tools.inline('return_val = a.is_callable();',['a'])
        assert not res

class TestObjectCall(NumpyTestCase):
    def check_noargs(self,level=5):
        def foo():
            return (1,2,3)
        res = inline_tools.inline('return_val = foo.call();',['foo'])
        assert_equal(res,(1,2,3))
        assert_equal(sys.getrefcount(res),2)
    def check_args(self,level=5):
        def foo(val1,val2):
            return (val1,val2)
        code = """
               py::tuple args(2);
               args[0] = 1;
               args[1] = "hello";
               return_val = foo.call(args);
               """
        res = inline_tools.inline(code,['foo'])
        assert_equal(res,(1,"hello"))
        assert_equal(sys.getrefcount(res),2)
    def check_args_kw(self,level=5):
        def foo(val1,val2,val3=1):
            return (val1,val2,val3)
        code = """
               py::tuple args(2);
               args[0] = 1;
               args[1] = "hello";
               py::dict kw;
               kw["val3"] = 3;
               return_val = foo.call(args,kw);
               """
        res = inline_tools.inline(code,['foo'])
        assert_equal(res,(1,"hello",3))
        assert_equal(sys.getrefcount(res),2)
    def check_noargs_with_args(self,level=5):
        # calling a function that does take args with args
        # should fail.
        def foo():
            return "blah"
        code = """
               py::tuple args(2);
               args[0] = 1;
               args[1] = "hello";
               return_val = foo.call(args);
               """
        try:
            first = sys.getrefcount(foo)
            res = inline_tools.inline(code,['foo'])
        except TypeError:
            second = sys.getrefcount(foo)
            try:
                res = inline_tools.inline(code,['foo'])
            except TypeError:
                third = sys.getrefcount(foo)
        # first should == second, but the weird refcount error
        assert_equal(second,third)

class TestObjectMcall(NumpyTestCase):
    def check_noargs(self,level=5):
        a = foo()
        res = inline_tools.inline('return_val = a.mcall("bar");',['a'])
        assert_equal(res,"bar results")
        first = sys.getrefcount(res)
        del res
        res = inline_tools.inline('return_val = a.mcall("bar");',['a'])
        assert_equal(res,"bar results")
        second = sys.getrefcount(res)
        assert_equal(first,second)
    def check_args(self,level=5):
        a = foo()
        code = """
               py::tuple args(2);
               args[0] = 1;
               args[1] = "hello";
               return_val = a.mcall("bar2",args);
               """
        res = inline_tools.inline(code,['a'])
        assert_equal(res,(1,"hello"))
        assert_equal(sys.getrefcount(res),2)
    def check_args_kw(self,level=5):
        a = foo()
        code = """
               py::tuple args(2);
               args[0] = 1;
               args[1] = "hello";
               py::dict kw;
               kw["val3"] = 3;
               return_val = a.mcall("bar3",args,kw);
               """
        res = inline_tools.inline(code,['a'])
        assert_equal(res,(1,"hello",3))
        assert_equal(sys.getrefcount(res),2)
    def check_std_noargs(self,level=5):
        a = foo()
        method = "bar"
        res = inline_tools.inline('return_val = a.mcall(method);',['a','method'])
        assert_equal(res,"bar results")
        first = sys.getrefcount(res)
        del res
        res = inline_tools.inline('return_val = a.mcall(method);',['a','method'])
        assert_equal(res,"bar results")
        second = sys.getrefcount(res)
        assert_equal(first,second)
    def check_std_args(self,level=5):
        a = foo()
        method = "bar2"
        code = """
               py::tuple args(2);
               args[0] = 1;
               args[1] = "hello";
               return_val = a.mcall(method,args);
               """
        res = inline_tools.inline(code,['a','method'])
        assert_equal(res,(1,"hello"))
        assert_equal(sys.getrefcount(res),2)
    def check_std_args_kw(self,level=5):
        a = foo()
        method = "bar3"
        code = """
               py::tuple args(2);
               args[0] = 1;
               args[1] = "hello";
               py::dict kw;
               kw["val3"] = 3;
               return_val = a.mcall(method,args,kw);
               """
        res = inline_tools.inline(code,['a','method'])
        assert_equal(res,(1,"hello",3))
        assert_equal(sys.getrefcount(res),2)
    def check_noargs_with_args(self,level=5):
        # calling a function that does take args with args
        # should fail.
        a = foo()
        code = """
               py::tuple args(2);
               args[0] = 1;
               args[1] = "hello";
               return_val = a.mcall("bar",args);
               """
        try:
            first = sys.getrefcount(a)
            res = inline_tools.inline(code,['a'])
        except TypeError:
            second = sys.getrefcount(a)
            try:
                res = inline_tools.inline(code,['a'])
            except TypeError:
                third = sys.getrefcount(a)
        # first should == second, but the weird refcount error
        assert_equal(second,third)

class TestObjectHash(NumpyTestCase):
    def check_hash(self,level=5):
        class foo:
            def __hash__(self):
                return 123
        a= foo()
        res = inline_tools.inline('return_val = a.hash(); ',['a'])
        print 'hash:', res
        assert_equal(res,123)

class TestObjectIsTrue(NumpyTestCase):
    def check_true(self,level=5):
        class foo:
            pass
        a= foo()
        res = inline_tools.inline('return_val = a.is_true();',['a'])
        assert_equal(res,1)
    def check_false(self,level=5):
        a= None
        res = inline_tools.inline('return_val = a.is_true();',['a'])
        assert_equal(res,0)

class TestObjectIsTrue(NumpyTestCase):
    def check_false(self,level=5):
        class foo:
            pass
        a= foo()
        res = inline_tools.inline('return_val = a.mcall("not");',['a'])
        assert_equal(res,0)
    def check_true(self,level=5):
        a= None
        res = inline_tools.inline('return_val = a.mcall("not");',['a'])
        assert_equal(res,1)

class TestObjectType(NumpyTestCase):
    def check_type(self,level=5):
        class foo:
            pass
        a= foo()
        res = inline_tools.inline('return_val = a.type();',['a'])
        assert_equal(res,type(a))

class TestObjectSize(NumpyTestCase):
    def check_size(self,level=5):
        class foo:
            def __len__(self):
                return 10
        a= foo()
        res = inline_tools.inline('return_val = a.size();',['a'])
        assert_equal(res,len(a))
    def check_len(self,level=5):
        class foo:
            def __len__(self):
                return 10
        a= foo()
        res = inline_tools.inline('return_val = a.len();',['a'])
        assert_equal(res,len(a))
    def check_length(self,level=5):
        class foo:
            def __len__(self):
                return 10
        a= foo()
        res = inline_tools.inline('return_val = a.length();',['a'])
        assert_equal(res,len(a))

from UserList import UserList
class TestObjectSetItemOpIndex(NumpyTestCase):
    def check_list_refcount(self,level=5):
        a = UserList([1,2,3])
        # temporary refcount fix until I understand why it incs by one.
        inline_tools.inline("a[1] = 1234;",['a'])
        before1 = sys.getrefcount(a)
        after1 = sys.getrefcount(a)
        assert_equal(after1,before1)
    def check_set_int(self,level=5):
        a = UserList([1,2,3])
        inline_tools.inline("a[1] = 1234;",['a'])
        assert_equal(sys.getrefcount(a[1]),2)
        assert_equal(a[1],1234)
    def check_set_double(self,level=5):
        a = UserList([1,2,3])
        inline_tools.inline("a[1] = 123.0;",['a'])
        assert_equal(sys.getrefcount(a[1]),2)
        assert_equal(a[1],123.0)
    def check_set_char(self,level=5):
        a = UserList([1,2,3])
        inline_tools.inline('a[1] = "bubba";',['a'])
        assert_equal(sys.getrefcount(a[1]),2)
        assert_equal(a[1],'bubba')
    def check_set_string(self,level=5):
        a = UserList([1,2,3])
        inline_tools.inline('a[1] = std::string("sissy");',['a'])
        assert_equal(sys.getrefcount(a[1]),2)
        assert_equal(a[1],'sissy')
    def check_set_string(self,level=5):
        a = UserList([1,2,3])
        inline_tools.inline('a[1] = std::complex<double>(1,1);',['a'])
        assert_equal(sys.getrefcount(a[1]),2)
        assert_equal(a[1],1+1j)

from UserDict import UserDict
class TestObjectSetItemOpKey(NumpyTestCase):
    def check_key_refcount(self,level=5):
        a = UserDict()
        code =  """
                py::object one = 1;
                py::object two = 2;
                py::tuple ref_counts(3);
                py::tuple obj_counts(3);
                py::tuple val_counts(3);
                py::tuple key_counts(3);
                obj_counts[0] = a.refcount();
                key_counts[0] = one.refcount();
                val_counts[0] = two.refcount();
                a[1] = 2;
                obj_counts[1] = a.refcount();
                key_counts[1] = one.refcount();
                val_counts[1] = two.refcount();
                a[1] = 2;
                obj_counts[2] = a.refcount();
                key_counts[2] = one.refcount();
                val_counts[2] = two.refcount();

                ref_counts[0] = obj_counts;
                ref_counts[1] = key_counts;
                ref_counts[2] = val_counts;
                return_val = ref_counts;
                """
        obj,key,val = inline_tools.inline(code,['a'])
        assert_equal(obj[0],obj[1])
        assert_equal(obj[1],obj[2])
        assert_equal(key[0] + 1, key[1])
        assert_equal(key[1], key[2])
        assert_equal(val[0] + 1, val[1])
        assert_equal(val[1], val[2])

    def check_set_double_exists(self,level=5):
        a = UserDict()
        key = 10.0
        a[key] = 100.0
        inline_tools.inline('a[key] = 123.0;',['a','key'])
        first = sys.getrefcount(key)
        inline_tools.inline('a[key] = 123.0;',['a','key'])
        second = sys.getrefcount(key)
        assert_equal(first,second)
        # !! I think the following should be 3
        assert_equal(sys.getrefcount(key),5)
        assert_equal(sys.getrefcount(a[key]),2)
        assert_equal(a[key],123.0)
    def check_set_double_new(self,level=5):
        a = UserDict()
        key = 1.0
        inline_tools.inline('a[key] = 123.0;',['a','key'])
        assert_equal(sys.getrefcount(key),4) # should be 3
        assert_equal(sys.getrefcount(a[key]),2)
        assert_equal(a[key],123.0)
    def check_set_complex(self,level=5):
        a = UserDict()
        key = 1+1j
        inline_tools.inline("a[key] = 1234;",['a','key'])
        assert_equal(sys.getrefcount(key),3)
        assert_equal(sys.getrefcount(a[key]),2)
        assert_equal(a[key],1234)
    def check_set_char(self,level=5):
        a = UserDict()
        inline_tools.inline('a["hello"] = 123.0;',['a'])
        assert_equal(sys.getrefcount(a["hello"]),2)
        assert_equal(a["hello"],123.0)

    def check_set_class(self,level=5):
        a = UserDict()
        class foo:
            def __init__(self,val):
                self.val = val
            def __hash__(self):
                return self.val
        key = foo(4)
        inline_tools.inline('a[key] = "bubba";',['a','key'])
        first = sys.getrefcount(key)
        inline_tools.inline('a[key] = "bubba";',['a','key'])
        second = sys.getrefcount(key)
        # I don't think we're leaking if this is true
        assert_equal(first,second)
        # !! BUT -- I think this should be 3
        assert_equal(sys.getrefcount(key),4)
        assert_equal(sys.getrefcount(a[key]),2)
        assert_equal(a[key],'bubba')
    def check_set_from_member(self,level=5):
        a = UserDict()
        a['first'] = 1
        a['second'] = 2
        inline_tools.inline('a["first"] = a["second"];',['a'])
        assert_equal(a['first'],a['second'])

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        sys.argv.extend(["--level=5"])
    NumpyTest().run()
