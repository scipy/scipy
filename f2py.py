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

from interface_gen \
    import \
        generate_interface


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
#F2PY_INCLUDE_MATCH = re.compile(r"^\s+\<include_file=(\S+)\>")
# Copied from crackfortran.py in f2py
F2PY_INCLUDE_MATCH = re.compile(r'\s*include\s*(\'|")(?P<name>[^\'"]*)(\'|")',re.I)

def is_pyf(node):
    return node.name.endswith(".pyf")

def include_pyf(node):
    includes = []
    for line in node.read().splitlines():
        m = F2PY_INCLUDE_MATCH.match(line)
        if m:
            includes.append(m.groupdict()["name"])
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

@feature('f2py_fortran')
@before_method('process_source')
def add_fortran_tasks(task_gen):
    assert "f2py" in task_gen.features
    
    if hasattr(task_gen, "name"):
        module_name = task_gen.name
    else:
        raise ValueError("Please define a name for task_gen %r" % task_gen)
    if len(task_gen.source) < 1:
        raise ValueError("Gne ?")
    sources = task_gen.to_nodes(task_gen.source)
    module_node = sources[0].parent.find_or_declare("%smodule.c" % module_name)

    f2py_task = task_gen.create_task('f2py_fortran', sources, module_node)
    add_f2py_tasks(task_gen, f2py_task, module_name, module_node)
    f2py_task.env.F2PYFLAGS.extend(["--lower", "-m", module_name])

def fake_interface_gen_callback(task_gen, node):
    return _interface_gen_callback(task_gen, node, "f2py_fake_interface")

def interface_gen_callback(task_gen, node):
    return _interface_gen_callback(task_gen, node, "f2py_interface")

def _interface_gen_callback(task_gen, node, interface_task_name):
    if not hasattr(task_gen, "name"):
        module_name = f2py_modulename(node)
    else:
        module_name = task_gen.name

    intermediate_output = node.parent.find_or_declare("%s.pyf" % module_name)
    module_node = node.parent.find_or_declare("%smodule.c" % module_name)
    # Guard against overwriting existing source code by accident. Did I say I
    # hate find_or_declare ?
    assert intermediate_output.is_bld()
    assert module_node.is_bld()

    interface_task = task_gen.create_task(interface_task_name, node, intermediate_output)

    f2py_task = task_gen.create_task('f2py', intermediate_output, module_node)
    add_f2py_tasks(task_gen, f2py_task, module_name, module_node)

def add_f2py_tasks(task_gen, f2py_task, module_name, module_node):
    build_dir = module_node.parent.bldpath()
    # FIXME: ask waf ML how flags sharing and co is supposed to work
    f2py_task.env.F2PYFLAGS = task_gen.env.F2PYFLAGS[:]
    f2py_task.env.F2PYFLAGS.extend(["--build-dir", build_dir])
    task_gen.source += f2py_task.outputs

    fobject_node = module_node.parent.find_or_declare("%s-fortranobject.c" % module_name)
    fobject_task = task_gen.create_task("f2py_fortran_object", [], fobject_node)
    task_gen.source += fobject_task.outputs

    fwrapper_node = module_node.parent.find_or_declare(FWRAP_TEMPLATE % module_name)
    fwrapper_task = task_gen.create_task("f2py_fwrapper", [], fwrapper_node)
    task_gen.source += fwrapper_task.outputs

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

    f2py_task = task_gen.create_task('f2py', node, module_node)
    add_f2py_tasks(task_gen, f2py_task, module_name, module_node)

@extension('.ipyf')
def add_f2py_files(task_gen, node):
    ext = '.pyf'

    if not hasattr(task_gen, "name"):
        module_name = f2py_modulename(node)
    else:
        module_name = task_gen.name
    intermediate_output = node.parent.find_or_declare("%s.pyf" % module_name)
    assert intermediate_output.is_bld()

    if "f2py_interface" in task_gen.features:
        interface_task = task_gen.create_task("f2py_interface", node, intermediate_output)
    elif "f2py_fake_interface" in task_gen.features:
        interface_task = task_gen.create_task("f2py_fake_interface", node, intermediate_output)
    else:
        raise ValueError("You need to use f2py_interface or f2py_fake_interface for .ipyf !")
    task_gen.source += interface_task.outputs
    #module_node = node.parent.find_or_declare("%smodule.c" % module_name)

    #f2py_task = task_gen.create_task('f2py', node, module_node)
    #add_f2py_tasks(task_gen, f2py_task, module_name, module_node)

class _f2py_interface(Task.Task):
    pass

class f2py_fake_interface(_f2py_interface):
    def run(self):
        node = self.inputs[0]
        output = self.outputs[0]
        module_name = f2py_modulename(node)
        generate_fake_interface(module_name, output)
        return 0

class f2py_interface(_f2py_interface):
    def run(self):
        node = self.inputs[0]
        output = self.outputs[0]
        module_name = f2py_modulename(node)
        generate_interface(module_name, node.abspath(), output.abspath())
        return 0

    def scan(self):
        found = []
        missing = []

        node = self.inputs[0]
        if not hasattr(self, "name"):
            module_name = f2py_modulename(node)
        else:
            module_name = self.name

        includes = include_pyf(node)
        if includes:
            real_pyf = node.parent.find_or_declare("%s.pyf" % module_name)
            # Guard against overwriting existing source code by accident. Did I
            # say I hate find_or_declare ?
            assert real_pyf.is_bld()

            for inc in includes:
                x = node.parent.find_resource(inc)
                if x:
                    found.append(x)
                else:
                    missing.append(x)
        user_routines = node.parent.find_resource("%s_user_routines.pyf" % module_name)
        if user_routines:
            found.append(user_routines)
        return (found, missing)

class f2py_fortran_object(Task.Task):
    def run(self):
        node = self.generator.bld.bldnode.find_node(op.join(F2PY_TEMP_DIR, FOBJECT_FILE))
        output = self.outputs[0]
        output.write(node.read())
        return 0

class f2py_fwrapper(Task.Task):
    after = ["f2py"]
    def run(self):
        output = self.outputs[0]
        if not op.exists(output.abspath()):
            output.write("")
        return 0

# TODO: f2py from .pyf or from .f should be different tasks.
class _f2py_task(Task.Task):
    run_str = '${F2PY} ${F2PYFLAGS} ${SRC}'
    color   = 'CYAN'

class f2py(_f2py_task):
    ext_in = [".pyf"]
    ext_out = [".h", ".f", ".c"]

class f2py_fortran(_f2py_task):
    ext_in = [".pyf"]
    ext_out = [".h", ".f", ".c"]

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
