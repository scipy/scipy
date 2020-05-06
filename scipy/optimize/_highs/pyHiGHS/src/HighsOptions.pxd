# distutils: language=c++
# cython: language_level=3

from libc.stdio cimport FILE

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from HConst cimport HighsOptionType

cdef extern from "HighsOptions.h" nogil:

    cdef cppclass OptionRecord:
        HighsOptionType type
        string name
        string description
        bool advanced

    cdef cppclass OptionRecordBool(OptionRecord):
        bool default_value

    cdef cppclass OptionRecordInt(OptionRecord):
        int lower_bound
        int default_value
        int upper_bound

    cdef cppclass OptionRecordDouble(OptionRecord):
        double lower_bound
        double default_value
        double upper_bound

    cdef cppclass OptionRecordString(OptionRecord):
        string default_value

    cdef cppclass HighsOptions:
        FILE * output
        int message_level
        bool mip
        string solution_file
        bool write_solution_to_file
        bool write_solution_pretty

        vector[OptionRecord*] records
