import os.path as op
import re

import numpy.f2py
import numpy.distutils.misc_util

from waflib \
    import \
        Task
from waflib.TaskGen \
    import \
        extension, feature, before_method, after_method

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
F2PY_INCLUDE_MATCH = re.compile(r"^\s+\<include_file=(\S+)\>")

def is_pyf(node):
    return node.name.endswith(".pyf")

def include_pyf(node):
    includes = []
    for line in node.read().splitlines():
        m = F2PY_INCLUDE_MATCH.match(line)
        if m:
            includes.append(m.group(1))
    return includes

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

def generate_fake_interface(name, node):
    """Generate a (fake) .pyf file from another pyf file (!)."""
    content = """\
python module %(name)s
    usercode void empty_module(void) {}
    interface
    subroutine empty_module()
    intent(c) empty_module
    end subroutine empty_module
    end interface
end python module%(name)s
"""
    node.write(content % {"name": name})

@feature('f2py')
@before_method('apply_incpaths')
def apply_f2py_includes(task_gen):
    includes = task_gen.env["INCLUDES"]
    d = op.dirname(numpy.f2py.__file__)
    includes.append(op.join(d, 'src'))
    includes.extend(numpy.distutils.misc_util.get_numpy_include_dirs())
    task_gen.env["INCLUDES"] = includes

@feature('f2py')
@before_method('apply_link')
def apply_f2py_fortran_sources(task_gen):
    for s in task_gen.source:
        if is_pyf(s):
            return

    if not hasattr(task_gen, "name"):
        module_name = f2py_modulename(task_gen.source[0])
    else:
        module_name = task_gen.name

    module_node = task_gen.source[0].parent.find_or_declare("%smodule.c" % module_name)
    task_gen.create_compiled_task("c", module_node)

    build_dir = module_node.parent.bldpath()

    source = task_gen.source[:]
    tsk = task_gen.create_task('f2py_fortran', source, module_node)
    # FIXME: ask waf ML how flags sharing and co is supposed to work
    tsk.env.F2PYFLAGS = task_gen.env.F2PYFLAGS[:]
    tsk.env.F2PYFLAGS.extend(["--build-dir", build_dir])
    tsk.env.F2PYFLAGS.extend(["--lower", "-m", module_name])
    #add_f2py_extra(task_gen, module_node)
    task_gen.source += tsk.outputs

    parent_node = module_node.parent

    # Make sure module_node is a generated file to avoid overwriting user
    # content (I really don't like find_or_declare).
    assert module_node.is_bld()

    fortranobject_source_node = task_gen.bld.bldnode.find_node(op.join(F2PY_TEMP_DIR, FOBJECT_FILE))
    assert fortranobject_source_node is not None
    fortranobject_node = parent_node.make_node("fortranobject%s.c" % module_name)
    fortranobject_node.write(fortranobject_source_node.read())

    # FIXME: race condition 
    fwrapper_node = parent_node.make_node(FWRAP_TEMPLATE % module_name)
    fwrapper_node.write("")

    # XXX: evil hack to make get_bld_sig work here. Find out how to do this properly
    fwrapper_node.is_bld = lambda : False
    fortranobject_node.is_bld = lambda : False

    task_gen.create_compiled_task("c", fortranobject_node)
    task_gen.create_compiled_task("fc", fwrapper_node)

def add_f2py_extra(task_gen, module_node, module_name):
    parent_node = module_node.parent

    # Make sure module_node is a generated file to avoid overwriting user
    # content (I really don't like find_or_declare).
    assert module_node.is_bld()

    fortranobject_source_node = task_gen.bld.bldnode.find_node(op.join(F2PY_TEMP_DIR, FOBJECT_FILE))
    assert fortranobject_source_node is not None
    fortranobject_node = parent_node.make_node("%s-fortranobject.c" % module_name)
    fortranobject_node.write(fortranobject_source_node.read())

    # FIXME: race condition 
    fwrapper_node = parent_node.make_node(FWRAP_TEMPLATE % module_name)
    fwrapper_node.write("")
    # XXX: evil hack to make get_bld_sig work here. Find out how to do this properly
    fwrapper_node.is_bld = lambda : False
    fortranobject_node.is_bld = lambda : False

    task_gen.source.append(fortranobject_node)
    task_gen.source.append(fwrapper_node)

@extension('.pyf')
def add_f2py_files(task_gen, node):
    ext = '.c'

    if not is_pyf(node):
        raise ValueError("Gne ?")

    if not hasattr(task_gen, "name"):
        module_name = f2py_modulename(node)
    else:
        module_name = task_gen.name
    module_node = node.parent.find_or_declare("%smodule.c" % module_name)

    includes = include_pyf(node)
    use_interface = ("f2py_interface_gen" in task_gen.features) \
        or ("f2py_fake_interface_gen" in task_gen.features) \
        or len(includes) > 0
    if use_interface:
        from interface_gen import generate_interface

        real_pyf = node.parent.find_or_declare("%s.pyf" % module_name)
        # Guard against overwriting existing source code by accident. Did I
        # say I hate find_or_declare ?
        assert real_pyf.is_bld()
        if "f2py_fake_interface_gen" in task_gen.features:
            generate_fake_interface(module_name, real_pyf)
        else:
            generate_interface(module_name, node.abspath(), real_pyf.abspath())
        node = real_pyf

    # XXX: evil hack to make get_bld_sig work here. Find out how to do this properly
    node.is_bld = lambda : False

    add_f2py_extra(task_gen, module_node, module_name)

    tsk = task_gen.create_task('f2py', node, module_node)
    build_dir = module_node.parent.bldpath()
    # FIXME: ask waf ML how flags sharing and co is supposed to work
    tsk.env.F2PYFLAGS = task_gen.env.F2PYFLAGS[:]
    tsk.env.F2PYFLAGS.extend(["--build-dir", build_dir])
    task_gen.source += tsk.outputs

class f2py(Task.Task):
    run_str = '${F2PY} ${F2PYFLAGS} ${SRC}'
    color   = 'CYAN'

class f2py_fortran(Task.Task):
    run_str = '${F2PY} ${F2PYFLAGS} ${SRC}'
    color   = 'CYAN'
    ext_out = [".h"]

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
