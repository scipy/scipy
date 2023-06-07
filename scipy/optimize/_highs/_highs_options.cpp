extern "C" {
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdbool.h>
} // extern "C"

#include <mutex>
#include <map>
#include <cstdio>
#include "HighsOptions.h"

std::mutex highs_options_mutex;
char buf[256];
static bool log_to_console = false;
static bool output_flag = true;
static HighsLogOptions highs_log_options = {
    // MSVC can't handle designated initializers
    /* .log_stream = */ nullptr,
    /* .output_flag = */ &output_flag,
    /* .log_to_console = */ &log_to_console,
    /* .log_dev_level = */ nullptr
};

struct cmp_str {
   bool operator()(char const *a, char const *b) const {
      return std::strcmp(a, b) < 0;
   }
};
static std::map<const char*, HighsOptionType, cmp_str> record_type_lookup;

static HighsOptions& get_highs_options() {
    static HighsOptions highs_options;
    static bool made_lookup = false;
    if (!made_lookup) {
        for (const auto record : highs_options.records) {
            record_type_lookup.emplace(record->name.c_str(), record->type);
        }
        made_lookup = true;
    }
    return highs_options;
}

template<typename OptionRecordType, typename T>
struct CheckOptionAction {
    OptionStatus operator()(const HighsLogOptions& highs_log_options,
                              OptionRecordType& record,
                              const T value) const {
        return checkOptionValue(highs_log_options, record, value);
    }
};


template<typename OptionRecordType, typename T>
static PyObject *check_option(PyObject *self, PyObject *args, const char this_fmt) {
    int res = 0;
    const char* name = nullptr;
    T value = 0;
    char fmt[] = {
        's',       // name
        this_fmt,  // value
        '\0'
    };
    if (!PyArg_ParseTuple(args, fmt, &name, &value)) {
        return nullptr;
    }

    std::FILE* tmpf = std::tmpfile();
    std::lock_guard<std::mutex> guard(highs_options_mutex);
    HighsOptions& highs_options = get_highs_options();
    buf[0] = '\0';
    highs_log_options.log_stream = tmpf;
    HighsInt idx = 0;
    const OptionStatus idx_status = getOptionIndex(highs_log_options, name, highs_options.records, idx);
    if (OptionStatus::kOk != idx_status) {
        res = 1;  // not a valid name
    }
    if (0 == res) {
        OptionRecordType& record = static_cast<OptionRecordType&>(*highs_options.records.at(idx));
        const OptionStatus check_status = CheckOptionAction<OptionRecordType, T>()(highs_log_options, record, value);
        if (OptionStatus::kIllegalValue == check_status) {
            res = 2;  // not a valid value
        }
    }

    if (0 != res) {
        std::rewind(tmpf);
        if (nullptr == std::fgets(buf, sizeof(buf), tmpf)) {
            std::sprintf(buf, "Failed to read HiGHS log output");
        }
    } else {
        std::sprintf(buf, "Check option succeeded");
    }
    std::fclose(tmpf);
    return Py_BuildValue("is", res, buf);
}

extern "C" {

static PyObject *get_option_type(PyObject *self, PyObject *args) {
    const char* name = nullptr;
    char fmt[] = {'s', '\0'};  // name
    if (!PyArg_ParseTuple(args, fmt, &name)) {
        return nullptr;
    }

    // find type
    HighsOptions& highs_options = get_highs_options();
    const auto lookup = record_type_lookup.find(name);
    if (record_type_lookup.cend() == lookup) {
        return Py_BuildValue("i", -1);
    }
    return Py_BuildValue("i", static_cast<int>(lookup->second));
}

static PyObject *check_int_option(PyObject *self, PyObject *args) {
    return check_option<OptionRecordInt, int>(self, args, 'i');
}

static PyObject *check_double_option(PyObject *self, PyObject *args) {
    return check_option<OptionRecordDouble, double>(self, args, 'd');
}

static PyObject *check_string_option(PyObject *self, PyObject *args) {
    return check_option<OptionRecordString, const char*>(self, args, 's');
}


static PyMethodDef highsOptionsMethods[] = {
    {"get_option_type", get_option_type, METH_VARARGS, ""},
    {"check_int_option", check_int_option, METH_VARARGS, ""},
    {"check_double_option", check_double_option, METH_VARARGS, ""},
    {"check_string_option", check_string_option, METH_VARARGS, ""},
    {nullptr, nullptr, 0, nullptr} /* Sentinel */
};


static struct PyModuleDef _highs_options_ext = {
    PyModuleDef_HEAD_INIT, "_highs_options", /* name of module */
    nullptr,                              /* module documentation, may be NULL */
    -1, /* size of per-interpreter state of the module, or -1 if the module
   keeps state in global variables. */
    highsOptionsMethods};


PyMODINIT_FUNC PyInit__highs_options(void) {
    return PyModule_Create(&_highs_options_ext);
}

} // extern "C"
