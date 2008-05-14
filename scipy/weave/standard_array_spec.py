from c_spec import common_base_converter
from c_spec import num_to_c_types
import numpy

num_typecode = {}
num_typecode['?'] = 'PyArray_BOOL'
num_typecode['b'] = 'PyArray_BYTE'
num_typecode['B'] = 'PyArray_UBYTE'
num_typecode['h'] = 'PyArray_SHORT'
num_typecode['H'] = 'PyArray_USHORT'
num_typecode['i'] = 'PyArray_INT'
num_typecode['I'] = 'PyArray_UINT'
num_typecode['l'] = 'PyArray_LONG'
num_typecode['L'] = 'PyArray_ULONG'
num_typecode['q'] = 'PyArray_LONGLONG'
num_typecode['Q'] = 'PyArray_ULONGLONG'
num_typecode['f'] = 'PyArray_FLOAT'
num_typecode['d'] = 'PyArray_DOUBLE'
num_typecode['g'] = 'PyArray_LONGDOUBLE'
num_typecode['F'] = 'PyArray_CFLOAT'
num_typecode['D'] = 'PyArray_CDOUBLE'
num_typecode['G'] = 'PyArray_CLONGDOUBLE'

type_check_code = \
"""
class numpy_type_handler
{
public:
    void conversion_numpy_check_type(PyArrayObject* arr_obj, int numeric_type,
                                     const char* name)
    {
        // Make sure input has correct numeric type.
        int arr_type = arr_obj->descr->type_num;
        if (PyTypeNum_ISEXTENDED(numeric_type))
        {
        char msg[80];
        sprintf(msg, "Conversion Error: extended types not supported for variable '%s'",
                name);
        throw_error(PyExc_TypeError, msg);
        }
        if (!PyArray_EquivTypenums(arr_type, numeric_type))
        {

        const char* type_names[23] = {"bool", "byte", "ubyte","short", "ushort",
                                "int", "uint", "long", "ulong", "longlong", "ulonglong",
                                "float", "double", "longdouble", "cfloat", "cdouble",
                                "clongdouble", "object", "string", "unicode", "void", "ntype",
                                "unknown"};
        char msg[500];
        sprintf(msg,"Conversion Error: received '%s' typed array instead of '%s' typed array for variable '%s'",
                type_names[arr_type],type_names[numeric_type],name);
        throw_error(PyExc_TypeError,msg);
        }
    }

    void numpy_check_type(PyArrayObject* arr_obj, int numeric_type, const char* name)
    {
        // Make sure input has correct numeric type.
        int arr_type = arr_obj->descr->type_num;
        if (PyTypeNum_ISEXTENDED(numeric_type))
        {
        char msg[80];
        sprintf(msg, "Conversion Error: extended types not supported for variable '%s'",
                name);
        throw_error(PyExc_TypeError, msg);
        }
        if (!PyArray_EquivTypenums(arr_type, numeric_type))
        {
            const char* type_names[23] = {"bool", "byte", "ubyte","short", "ushort",
                                    "int", "uint", "long", "ulong", "longlong", "ulonglong",
                                    "float", "double", "longdouble", "cfloat", "cdouble",
                                    "clongdouble", "object", "string", "unicode", "void", "ntype",
                                    "unknown"};
            char msg[500];
            sprintf(msg,"received '%s' typed array instead of '%s' typed array for variable '%s'",
                    type_names[arr_type],type_names[numeric_type],name);
            throw_error(PyExc_TypeError,msg);
        }
    }
};

numpy_type_handler x__numpy_type_handler = numpy_type_handler();
#define conversion_numpy_check_type x__numpy_type_handler.conversion_numpy_check_type
#define numpy_check_type x__numpy_type_handler.numpy_check_type

"""

size_check_code = \
"""
class numpy_size_handler
{
public:
    void conversion_numpy_check_size(PyArrayObject* arr_obj, int Ndims,
                                     const char* name)
    {
        if (arr_obj->nd != Ndims)
        {
            char msg[500];
            sprintf(msg,"Conversion Error: received '%d' dimensional array instead of '%d' dimensional array for variable '%s'",
                    arr_obj->nd,Ndims,name);
            throw_error(PyExc_TypeError,msg);
        }
    }

    void numpy_check_size(PyArrayObject* arr_obj, int Ndims, const char* name)
    {
        if (arr_obj->nd != Ndims)
        {
            char msg[500];
            sprintf(msg,"received '%d' dimensional array instead of '%d' dimensional array for variable '%s'",
                    arr_obj->nd,Ndims,name);
            throw_error(PyExc_TypeError,msg);
        }
    }
};

numpy_size_handler x__numpy_size_handler = numpy_size_handler();
#define conversion_numpy_check_size x__numpy_size_handler.conversion_numpy_check_size
#define numpy_check_size x__numpy_size_handler.numpy_check_size

"""

numeric_init_code = \
"""
Py_Initialize();
import_array();
PyImport_ImportModule("numpy");
"""

class array_converter(common_base_converter):

    def init_info(self):
        common_base_converter.init_info(self)
        self.type_name = 'numpy'
        self.check_func = 'PyArray_Check'
        self.c_type = 'PyArrayObject*'
        self.return_type = 'PyArrayObject*'
        self.to_c_return = '(PyArrayObject*) py_obj'
        self.matching_types = [numpy.ndarray]
        self.headers = ['"numpy/arrayobject.h"',
                        '<complex>','<math.h>']
        self.support_code = [size_check_code, type_check_code]
        self.module_init_code = [numeric_init_code]

    def get_var_type(self,value):
        return value.dtype.char

    def template_vars(self,inline=0):
        res = common_base_converter.template_vars(self,inline)
        if hasattr(self,'var_type'):
            res['num_type'] = num_to_c_types[self.var_type]
            res['num_typecode'] = num_typecode[self.var_type]
        res['array_name'] = self.name + "_array"
        res['cap_name'] = self.name.upper()
        return res

    def declaration_code(self,templatize = 0,inline=0):
        res = self.template_vars(inline=inline)
        cap_name = self.name.upper()
        res['cap_name'] = cap_name
        code2 = '#define %(cap_name)s1(i) (*((%(num_type)s*)(%(array_name)s->data + (i)*S%(name)s[0])))\n' \
                    '#define %(cap_name)s2(i,j) (*((%(num_type)s*)(%(array_name)s->data + (i)*S%(name)s[0] + (j)*S%(name)s[1])))\n' \
                    '#define %(cap_name)s3(i,j,k) (*((%(num_type)s*)(%(array_name)s->data + (i)*S%(name)s[0] + (j)*S%(name)s[1] + (k)*S%(name)s[2])))\n' \
                    '#define %(cap_name)s4(i,j,k,l) (*((%(num_type)s*)(%(array_name)s->data + (i)*S%(name)s[0] + (j)*S%(name)s[1] + (k)*S%(name)s[2] + (l)*S%(name)s[3])))\n'

        res['_code2_'] = code2 % res

        code = '%(py_var)s = %(var_lookup)s;\n'   \
               '%(c_type)s %(array_name)s = %(var_convert)s;\n'  \
               'conversion_numpy_check_type(%(array_name)s,%(num_typecode)s,"%(name)s");\n' \
               '%(_code2_)s' \
               'npy_intp* N%(name)s = %(array_name)s->dimensions;\n' \
               'npy_intp* S%(name)s = %(array_name)s->strides;\n' \
               'int D%(name)s = %(array_name)s->nd;\n' \
               '%(num_type)s* %(name)s = (%(num_type)s*) %(array_name)s->data;\n'
        code = code % res
        self.__doundef=1
        return code

    def cleanup_code(self):
        code = common_base_converter.cleanup_code(self)
        try:
            if self.__doundef != 1:
                return code
            cap_name = self.name.upper()
            newcode = "#undef %(cap_name)s1\n#undef %(cap_name)s2\n"\
                      "#undef %(cap_name)s3\n#undef %(cap_name)s4\n"\
                      % {'cap_name':cap_name}
            code = "%s%s" % (code, newcode)
        except AttributeError:
            pass
        return code
