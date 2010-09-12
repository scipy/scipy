#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#error YOYOYO
#endif
#include <Python.h>

#if PY_VERSION_HEX >= 0x03000000
	#define PyString_Check PyUnicode_Check
	static int __pyfile_check_guard(PyObject *x)
	{
		fprintf(stderr, "PY3K error: PyFile_Check called !\n");
		return 0;
	}
	#define PyFile_Check(x) __pyfile_check_guard((x))
	static int __pyinstance_check_guard(PyObject *x)
	{
		fprintf(stderr, "PY3K error: PyInstance_Check calleed !\n");
		return 0;
	}
	#define PyInstance_Check(x) __pyinstance_check_guard((x))
#endif
