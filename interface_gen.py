#!/usr/bin/env python

import os
import re
from distutils.dir_util import mkpath

def all_subroutines(interface_in):
    # remove comments
    comment_block_exp = re.compile(r'/\*(?:\s|.)*?\*/')
    subroutine_exp = re.compile(r'subroutine (?:\s|.)*?end subroutine.*')
    function_exp = re.compile(r'function (?:\s|.)*?end function.*')

    interface = comment_block_exp.sub('',interface_in)
    subroutine_list = subroutine_exp.findall(interface)
    function_list = function_exp.findall(interface)
    subroutine_list = subroutine_list + function_list
    subroutine_list = map(lambda x: x.strip(),subroutine_list)
    return subroutine_list

def real_convert(val_string):
    return val_string

def complex_convert(val_string):
    return '(' + val_string + ',0.)'

def convert_types(interface_in,converter):
    regexp = re.compile(r'<type_convert=(.*?)>')
    interface = interface_in[:]
    while 1:
        sub = regexp.search(interface)
        if sub is None: break
        converted = converter(sub.group(1))
        interface = interface.replace(sub.group(),converted)
    return interface

def generic_expand(generic_interface,skip_names=[]):
    generic_types ={'s' :('real',            'real', real_convert,
                          'real'),
                    'd' :('double precision','double precision',real_convert,
                          'double precision'),
                    'c' :('complex',         'complex',complex_convert,
                          'real'),
                    'z' :('double complex',  'double complex',complex_convert,
                          'double precision'),
                    'cs':('complex',         'real',complex_convert,
                          'real'),
                    'zd':('double complex',  'double precision',complex_convert,
                          'double precision'),
                    'sc':('real',            'complex',real_convert,
                          'real'),
                    'dz':('double precision','double complex', real_convert,
                          'double precision')}
    generic_c_types = {'real':'float',
                       'double precision':'double',
                       'complex':'complex_float',
                       'double complex':'complex_double'}
    # cc_types is specific in ATLAS C BLAS, in particular, for complex arguments
    generic_cc_types = {'real':'float',
                       'double precision':'double',
                       'complex':'void',
                       'double complex':'void'}
    #2. get all subroutines
    subs = all_subroutines(generic_interface)
    #print len(subs)
    #loop through the subs
    type_exp = re.compile(r'<tchar=(.*?)>')
    TYPE_EXP = re.compile(r'<TCHAR=(.*?)>')
    routine_name = re.compile(r'(subroutine|function)\s*(?P<name>\w+)\s*\(')
    interface = ''
    for sub in subs:
        #3. Find the typecodes to use:
        m = type_exp.search(sub)
        if m is None:
            interface = interface + '\n\n' + sub
            continue
        type_chars = m.group(1)
        # get rid of spaces
        type_chars = type_chars.replace(' ','')
        # get a list of the characters (or character pairs)
        type_chars = type_chars.split(',')
        # Now get rid of the special tag that contained the types
        sub = re.sub(type_exp,'<tchar>',sub)
        m = TYPE_EXP.search(sub)
        if m is not None:
            sub = re.sub(TYPE_EXP,'<TCHAR>',sub)
        sub_generic = sub.strip()
        for char in type_chars:
            type_in,type_out,converter, rtype_in = generic_types[char]
            sub = convert_types(sub_generic,converter)
            function_def = sub.replace('<tchar>',char)
            function_def = function_def.replace('<TCHAR>',char.upper())
            function_def = function_def.replace('<type_in>',type_in)
            function_def = function_def.replace('<type_in_c>',
                                          generic_c_types[type_in])
            function_def = function_def.replace('<type_in_cc>',
                                          generic_cc_types[type_in])
            function_def = function_def.replace('<rtype_in>',rtype_in)
            function_def = function_def.replace('<rtype_in_c>',
                                          generic_c_types[rtype_in])
            function_def = function_def.replace('<type_out>',type_out)
            function_def = function_def.replace('<type_out_c>',
                                          generic_c_types[type_out])
            m = routine_name.match(function_def)
            if m:
                if m.group('name') in skip_names:
                    print 'Skipping',m.group('name')
                    continue
            else:
                print 'Possible bug: Failed to determine routines name'
            interface = interface + '\n\n' + function_def

    return interface

#def interface_to_module(interface_in,module_name,include_list,sdir='.'):
def interface_to_module(interface_in,module_name):
    pre_prefix = "!%f90 -*- f90 -*-\n"
    # heading and tail of the module definition.
    file_prefix = "\npython module " + module_name +" ! in\n" \
                  "!usercode '''#include \"cblas.h\"\n"\
                  "!'''\n"\
                  "    interface  \n"
    file_suffix = "\n    end interface\n" \
             "end module %s" % module_name
    return  pre_prefix + file_prefix + interface_in + file_suffix

def process_includes(interface_in,sdir='.'):
    include_exp = re.compile(r'\n\s*[^!]\s*<include_file=(.*?)>')
    include_files = include_exp.findall(interface_in)
    for filename in include_files:
        f = open(os.path.join(sdir,filename))
        interface_in = interface_in.replace('<include_file=%s>'%filename,
                                      f.read())
        f.close()
    return interface_in

def generate_interface(module_name,src_file,target_file,skip_names=[]):
    #print "generating",module_name,"interface"
    f = open(src_file)
    generic_interface = f.read()
    f.close()
    sdir = os.path.dirname(src_file)
    generic_interface = process_includes(generic_interface,sdir)
    generic_interface = generic_expand(generic_interface,skip_names)
    module_def = interface_to_module(generic_interface,module_name)
    mkpath(os.path.dirname(target_file))
    f = open(target_file,'w')
    user_routines = os.path.join(sdir,module_name+"_user_routines.pyf")
    if os.path.exists(user_routines):
        f2 = open(user_routines)
        f.write(f2.read())
        f2.close()
    f.write(module_def)
    f.close()

def process_all():
    # process the standard files.
    for name in ['fblas','cblas','clapack','flapack']:
        generate_interface(name,'generic_%s.pyf'%(name),name+'.pyf')


if __name__ == "__main__":
    process_all()
