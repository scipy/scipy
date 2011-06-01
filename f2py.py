import os.path as op
import re

import numpy.f2py
import numpy.distutils.misc_util

from waflib import Task
from waflib.TaskGen import extension

CGEN_TEMPLATE   = '%smodule'
FOBJECT_FILE    = 'fortranobject.c'
FHEADER_FILE    = 'fortranobject.h'
FWRAP_TEMPLATE  = '%s-f2pywrappers.f'

# This path is relative to build directory
F2PY_TEMP_DIR = '.f2py'

# Those regex are copied from build_src in numpy.distutils.command
F2PY_MODNAME_MATCH = re.compile(r'\s*python\s*module\s*(?P<name>[\w_]+)',
                                re.I).match
F2PY_UMODNAME_MATCH = re.compile(r'\s*python\s*module\s*(?P<name>[\w_]*?'\
                                     '__user__[\w_]*)',re.I).match
# End of copy
def is_pyf(node):
    return node.name.endswith(".pyf")

def f2py_modulename(node):
    """This returns the name of the module from the pyf source file.

    source is expected to be one string, containing the whole source file
    code."""
    name = None
    for line in node.read().splitlines():
        m = F2PY_MODNAME_MATCH(line)
        if m:
            if F2PY_UMODNAME_MATCH(line): # skip *__user__* names
                continue
            name = m.group('name')
            break
    return name

@extension('.pyf')
def add_f2py_files(task_gen, node):
    ext = '.c'

    includes = task_gen.env["INCLUDES"]
    d = op.dirname(numpy.f2py.__file__)
    includes.append(op.join(d, 'src'))
    includes.extend(numpy.distutils.misc_util.get_numpy_include_dirs())
    task_gen.env["INCLUDES"] = includes

    if is_pyf(node):
        module_name = f2py_modulename(node)
        module_node = node.parent.find_or_declare("%smodule.c" % module_name)
        parent_node = module_node.parent

        # Make sure module_node is a generated file to avoid overwriting user
        # content (I really don't like find_or_declare).
        assert module_node.is_bld()
        build_dir = parent_node.bldpath()

        fortranobject_node = task_gen.bld.bldnode.find_node(op.join(F2PY_TEMP_DIR, FOBJECT_FILE))
        assert fortranobject_node is not None
        fwrapper_node = parent_node.find_or_declare(FWRAP_TEMPLATE % module_name)
        fwrapper_node.write("")

        # XXX: evil hack to make get_bld_sig work here. Find out how to do this properly
        fortranobject_node.is_bld = lambda : False
        fwrapper_node.is_bld = lambda : False
        task_gen.source.append(fortranobject_node)
        task_gen.source.append(fwrapper_node)

        task_gen.env.F2PYFLAGS.extend(["--build-dir", build_dir])

        tsk = task_gen.create_task('f2py', node, module_node)
        task_gen.source += tsk.outputs
    else:
        raise NotImplementedError("non .pyf input not supported yet.")

class f2py(Task.Task):
    run_str = '${F2PY} ${F2PYFLAGS} ${SRC}'
    color   = 'CYAN'

def configure(conf):
    if not conf.env.CC and not conf.env.CXX:
        conf.fatal('Load a C/C++ compiler first')
    if not conf.env.PYTHON:
        conf.fatal('Load the python tool first!')
    conf.find_program('f2py', var='F2PY')
    # FIXME: this has nothing to do here
    conf.env.F2PYFLAGS = ["--quiet"]

    f2py_tempdir = conf.bldnode.make_node(F2PY_TEMP_DIR)
    f2py_tempdir.mkdir()

    fobject = f2py_tempdir.make_node(FOBJECT_FILE)

    d = op.dirname(numpy.f2py.__file__)
    source_c = op.join(d, 'src', FOBJECT_FILE)
    fobject.write(open(source_c).read())

    #fheader = f2py_tempdir.make_node(FHEADER_FILE)
    #source_c = op.join(d, 'src', FHEADER_FILE)
    #fheader.write(open(source_c).read())
