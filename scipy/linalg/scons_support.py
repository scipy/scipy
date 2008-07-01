from os.path import join as pjoin, splitext, basename as pbasename

from interface_gen import generate_interface

def do_generate_interface(target, source, env):
    """Generate a .pyf file from another pyf file (!)."""
    # XXX: do this correctly
    target_name = str(target[0])
    source_name = str(source[0])

    # XXX handle skip names
    name = splitext(pbasename(target_name))[0]
    generate_interface(name, source_name, target_name)
    return 0

def generate_interface_emitter(target, source, env):
    base = str(target[0])
    return (['%s.pyf' % base], source)

def do_generate_fake_interface(target, source, env):
    """Generate a (fake) .pyf file from another pyf file (!)."""
    # XXX: do this correctly
    target_name = str(target[0])
    source_name = str(source[0])

    # XXX handle skip names
    name = splitext(pbasename(target_name))[0]
    generate_interface(name, source_name, target_name)

    f = open(target_name, 'w')
    f.write('python module '+name+'\n')
    f.write('usercode void empty_module(void) {}\n')
    f.write('interface\n')
    f.write('subroutine empty_module()\n')
    f.write('intent(c) empty_module\n')
    f.write('end subroutine empty_module\n')
    f.write('end interface\nend python module'+name+'\n')
    f.close()

    return 0
