import string, re, os

import interface_gen

this_dir,junk = os.path.split(interface_gen.__file__)

def process_imports(interface_in,sdir='.'):
    # !! lines can't contain comments, and a single import on each line
    #
    # !! expects interface to be named and live in same directory

    # insert <include_file=file> directly into interface

    include_exp = re.compile(r'<include_file=(.*?)>')
    # get a list of the import statements
    comment_block_exp = re.compile(r'/\*(?:\s|.)*?\*/')
    import_exp = re.compile(r'from\s*(.*?)\s*import\s*(.+)')
    # remove comments
    commentless_interface = comment_block_exp.sub('',interface_in) 

    include_files = include_exp.findall(commentless_interface)
    import_list = import_exp.findall(commentless_interface)
    # remove import statements
    interface = import_exp.sub('',commentless_interface)
    interface = include_exp.sub('',commentless_interface)

    external_files = {}
    for file,routine in import_list:
        if external_files.has_key(file):
            sub_list = external_files[file]
        else:
            f = open(os.path.join(sdir,file))
            sub_list = all_subroutines(f.read(-1))            
            f.close()
            external_files[file] = sub_list

        #search for this routine in the list
        search_string = routine + '('
        found = 0
        for sub in sub_list:
            if string.find(sub,search_string) != -1:
                # add subroutine to the bottom of this interface    
                interface = interface + '\n' + sub +'\n'
                found = 1
                break
        if not found:
            print '\twarning: %s not found in %s' % (routine,file)                
    return interface, include_files

def process_type(interface_in, format, reg_exp, fortran_decl, c_decl,
                 arg_count):
    if format == 'Fortran':    decl = fortran_decl
    elif format == 'C':        decl = c_decl
    else:
        raise ValueError, 'format must be "C" or "Fortran, not ' + format
    expr = re.compile(reg_exp)
    interface = interface_in[:]
    while 1:
        sub = expr.search(interface)
        if sub is None: break
        var_name = sub.group(1)
        replacements = (var_name,)*arg_count
        declaration = decl % replacements
        interface = string.replace(interface,sub.group(),declaration)
    return interface            

def process_transpose_types(interface_in,format='Fortran'):
    f_decl = "character optional, check(*%s=='N' || *%s == 'T' ||"\
             " *%s == 'C'):: %s = 'N'"                
    c_decl = "integer optional, check(%s==111 || %s == 112 ||"\
             " %s == 113):: %s = 111"
    reg_exp = r'<transpose_type=(.*?)>'
    arg_count = 4
    return process_type(interface_in,format,reg_exp,f_decl,c_decl,arg_count)
    
def process_uplo_types(interface_in,format='Fortran'):
    f_decl = "character optional, check(*%s=='L' || *%s == 'U'):: %s = 'U'"                
    c_decl = "integer optional, check(%s==121 || %s == 122):: %s = 121"
    reg_exp = r'<uplo_type=(.*?)>'
    arg_count = 3
    return process_type(interface_in,format,reg_exp,f_decl,c_decl,arg_count)
    
def process_diag_types(interface_in,format='Fortran'):
    f_decl = "character optional, check(*%s=='N' || *%s == 'U'):: %s = 'N'"
    c_decl = "integer optional, check(%s==131 || %s == 132):: %s = 131"
    reg_exp = r'<diag_type=(.*?)>'
    arg_count = 3
    return process_type(interface_in,format,reg_exp,f_decl,c_decl,arg_count)

def process_side_types(interface_in,format='Fortran'):
    f_decl = "character optional, check(*%s=='L' || *%s == 'R'):: %s = 'L'"
    c_decl = "integer optional, check(%s==141 || %s == 142):: %s = 141"
    reg_exp = r'<diag_type=(.*?)>'
    arg_count = 3
    return process_type(interface_in,format,reg_exp,f_decl,c_decl,arg_count)

def process_special_types(interface_in, format='Fortran'):
    interface = process_transpose_types(interface_in,format)
    interface = process_uplo_types(interface,format)
    interface = process_diag_types(interface,format)
    interface = process_side_types(interface,format)
    return interface
    
def insert_blas_order(interface_in):
    sub_list = all_subroutines(interface_in)
    #c_decl = "integer optional, check(order==101 || order == 102)::"\
    #         " order = 101"
    # for now, we force C arrays to be ordered as C arrays 
    #   -- novel idea...
    c_decl = "integer intent(hide) :: order = 101"
    interface = ''
    for sub in sub_list:
        blas_sub = string.replace(sub, "(" , "(order," ,1)
        #if string.find(blas_sub, "\n") == '-1': print "BROKEN!"
        #else: print string.find(blas_sub, '\n')
        blas_sub = string.replace(blas_sub, '\n' , '\n    ' + c_decl +'\n' ,1)
        interface = interface + blas_sub + '\n\n'
    return interface    

def insert_intent_c(interface_in):
    sub_list = all_subroutines(interface_in)
    c_decl = "intent(c)"
    interface = ''
    for sub in sub_list:
        new_sub = string.replace(sub, '\n' , '\n    ' + c_decl +'\n' ,1)
        interface = interface + new_sub + '\n\n'
    return interface
         
#def process_creturn(interface_in,format):
#    sub_list = all_subroutines(interface_in)
#    interface = ''
#    reg_exp = r'<c_return=(.*?)>'
#    expr = re.compile(reg_exp)    
#    for sub in sub_list:
#        creturn = expr.search(sub)
#        if creturn:
#            processed_sub = string.replace(sub,creturn.group(),'')
#            if string.lower(format) == 'c':
#                var_name = creturn.group(1)
#                #1. replace "subroutine" with "function"
#                processed_sub = string.replace(processed_sub,'subroutine',
#                                               'function')
#                reg_exp = r'<c_return=(.*?)>'
#                expr = re.compile(reg_exp)                    
#                #2. Find line with ':: var_name'
#                #   a. find func_name.
#                #   b. find type.
#                #   c. replace line with "type func_name"    
#                #3. Remove ",<ws>var_name" from function definition.
#        else:
#            processed_sub = sub
#        interface = interface + processed_sub + '\n\n'
#    return interface

def process_return_info(interface_in,format):
    sub_list = all_subroutines(interface_in)
    interface = ''
    
    for sub in sub_list:        
        if string.find(sub,'<return_info>') != -1:
            processed_sub = string.replace(sub,'<return_info>','')
            if string.lower(format) == 'c':
                # one of these two should get rif of info in the arg list
                processed_sub = string.replace(processed_sub,',info','')
                processed_sub = string.replace(processed_sub,', info','')
                #1. replace "subroutine" with "function"
                name_exp = re.compile(r'subroutine( .+?)\(')
                name = name_exp.search(processed_sub).group(1)
                processed_sub = string.replace(processed_sub,'subroutine',
                                               'function')
                sub_lines = string.split(processed_sub,'\n')
                processed_sub = ''
                for line in sub_lines: 
                    if (string.find(line,'info') != -1):
                        processed_sub =   processed_sub +  '    integer ' \
                                        + name + '\n'
                    else:
                        processed_sub = processed_sub + line + '\n'                        
                    
                #2. Find line with ':: var_name'
                #   a. find func_name.
                #   b. find type.
                #   c. replace line with "type func_name"    
                #3. Remove ",<ws>var_name" from function definition.
        else:
            processed_sub = sub

        interface = interface + processed_sub + '\n\n'
    return interface

def process_ignore_info(interface_in,format):
    sub_list = all_subroutines(interface_in)
    interface = ''
    for sub in sub_list:        
        if string.find(sub,'<ignore_info>') != -1:
            processed_sub = string.replace(sub,'<ignore_info>','')
            if string.lower(format) == 'c':
                # one of these two should get rif of info in the arg list
                processed_sub = string.replace(processed_sub,',info','')
                processed_sub = string.replace(processed_sub,', info','')

                sub_lines = string.split(processed_sub,'\n')
                processed_sub = ''
                for line in sub_lines: 
                    if (string.find(line,'info') == -1):
                        processed_sub = processed_sub + line + '\n'                                            
        else:
            processed_sub = sub
        interface = interface + processed_sub + '\n\n'
    return interface
        
def all_subroutines(interface_in):
    # remove comments
    comment_block_exp = re.compile(r'/\*(?:\s|.)*?\*/')
    subroutine_exp = re.compile(r'subroutine (?:\s|.)*?end subroutine.*')
    function_exp = re.compile(r'function (?:\s|.)*?end function.*')
    
    interface = comment_block_exp.sub('',interface_in)    
    subroutine_list = subroutine_exp.findall(interface)
    function_list = function_exp.findall(interface)
    function_list = []
    subroutine_list = subroutine_list + function_list 
    subroutine_list = map(lambda x: string.strip(x),subroutine_list)
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
        interface = string.replace(interface,sub.group(),converted)
    return interface

def convert_matrix_order(interface_in,row_major=0):
    # simple search and replace to specify which
    # axis is getting used.
    if row_major:
        first,second = ('0','1')
    else:
        first,second = ('1','0')
    order_exp = re.compile(r'<1st_axis>')
    interface = order_exp.sub(first,interface_in)
    order_exp = re.compile(r'<2nd_axis>')
    interface = order_exp.sub(second,interface)
    return interface
                
def lapack_expand(generic_interface,row_major = 0,cwrap=0):
    lapack_types ={'s' :('real',            'real', real_convert,
                         'real'), 
                   'd' :('double precision','double precision',real_convert,
                         'double precision'),
                   'c' :('complex',         'complex',complex_convert,
                         'real'), 
                   'z' :('double complex',  'double complex',complex_convert,
                         'double precision'),                   
                   'sc':('complex',         'real',real_convert,
                         'real'),
                   'dz':('double complex',  'double precision',real_convert,
                         'double precision'),
                   'cs':('real',            'complex',complex_convert,
                         'real'),
                   'zd':('double precision','double complex', complex_convert,
                         'double precision')}
    #1. Search and replace for matrix ordering (row or column).
    generic_interface = convert_matrix_order(generic_interface,row_major)
    #2. get all subroutines    
    subs = all_subroutines(generic_interface)
    #print len(subs)
    #loop through the subs
    type_exp = re.compile(r'<tchar=(.*?)>')    
    interface = ''
    for sub in subs:
        #3. Find the typecodes to use:
        type_chars = type_exp.search(sub).group(1)
        # get rid of spaces
        type_chars = string.replace(type_chars,' ','')
        # get a list of the characters (or character pairs)
        type_chars = string.split(type_chars,',')
        #print type_chars
        # Now get rid of the special tag that contained the types
        sub = re.sub(type_exp,'<tchar>',sub)
        sub_generic = string.strip(sub)
        for char in type_chars:            
            type_in,type_out,converter, rtype_in = lapack_types[char]
            sub = convert_types(sub_generic,converter)
            function_def = string.replace(sub,'<tchar>',char)    
            function_def = string.replace(function_def,'<type_in>',type_in)
            function_def = string.replace(function_def,'<rtype_in>',rtype_in) 
            function_def = string.replace(function_def,'<type_out>',type_out)
            interface = interface + '\n\n' + function_def
    
    return interface

def process_flip_dims(generic_interface, format):
    subs = all_subroutines(generic_interface)
    #print len(subs)
    #loop through the subs
    flip_exp = re.compile(r'<Cflip (.*?)>')
    interface = ''
    for sub in subs:   # replace <Cflip xxxx> with appropriate value.
        pos = 0
        newsub = ''
        while 1:
            match = flip_exp.search(sub,pos)
            if match is None:
                newsub += sub[pos:]
                break
            start,end = match.span(1)
            if format[0] in ['F','f']:
                str = sub[start:end]
                newsub += sub[pos:start-7] + str
                pos = end + 1
            elif format[0] in ['C', 'c']:
                str = sub[start:end]
                varbeg = str.find('(') + 1
                varend = str.find(')')
                if (varbeg == 0) or (varend == -1):
                    raise RuntimeError, "<Cflip > used incorrectly."
                vars = str[varbeg:varend].split(',')
                newvars = map(string.strip,vars)
                newvars.reverse()
                newstr = 'dimension(%s)' % (",".join(newvars))
                newsub += sub[pos:start-7] + newstr
                pos = end + 1
            else:
                raise ValueError, "Format invalid."    
        interface += newsub + "\n"
    return interface

def interface_to_module(interface_in,module_name,include_list,sdir='.'):
    pre_prefix = "!%f90 -*- f90 -*-\n"
    includes = ''
    for file in include_list:
        f = open(os.path.join(sdir,file))
        includes += f.read(-1)
        f.close() 
    # heading and tail of the module definition.
    file_prefix = "\npython module " + module_name +" ! in\n" \
             "    interface  \n"
    file_suffix = "\n    end interface\n" \
             "end module %s" % module_name
    return  pre_prefix + includes + file_prefix + interface_in + file_suffix

def function_name(function_interface):
    name_exp = re.compile(r'subroutine( .+?)\(')
    names = name_exp.findall(function_interface)
    name_exp = re.compile(r'function( .+?)\(')
    names = names + name_exp.findall(function_interface)
    assert(len(names)==1)
    return string.strip(names[0])

def rename_functions(interface_in,prefix,suffix):
    sub_list = all_subroutines(interface_in)    
    interface = ''
    for sub in sub_list:
        name = function_name(sub)
        new_name = prefix+name+suffix        
        c_decl = "fortranname(%s)" % new_name
        #renamed_sub = string.replace(sub, name ,new_name ,1)
        renamed_sub = sub
        renamed_sub = string.replace(renamed_sub, '\n' , 
                                                  '\n    ' + c_decl +'\n' ,1)
        interface = interface + renamed_sub + '\n\n'
    return interface    
    
def generate_flapack(sdir,output_path):
    print "generating flapack interface"
    f = open(os.path.join(sdir,'generic_lapack.pyf'))
    module_name = 'flapack'
    generic_interface = f.read(-1)
    generic_interface, include_files = process_imports(generic_interface,sdir)
    generic_interface = process_special_types(generic_interface,
                                              format='Fortran')
    interface = lapack_expand(generic_interface,row_major = 0)
    interface = process_ignore_info(interface,format='Fortran')
    interface = process_flip_dims(interface, format='Fortran')
    # must be last
    interface = process_return_info(interface,format='Fortran')
    module_definition = interface_to_module(interface,module_name, include_files,sdir)
    
    f = open(os.path.join(output_path,module_name+'.pyf'),'w')
    f.write(module_definition)
    f.close()
   
def generate_clapack(sdir,output_path):
    print "generating clapack interface"
    f = open(os.path.join(sdir,'atlas_lapack.pyf'))
    module_name = 'clapack'
    generic_interface = f.read(-1)
    generic_interface, include_files = process_imports(generic_interface,sdir)
    generic_interface = insert_blas_order(generic_interface)
    generic_interface = insert_intent_c(generic_interface)
    generic_interface = rename_functions(generic_interface,'clapack_','')
    generic_interface = process_special_types(generic_interface,format='C')

    generic_interface = lapack_expand(generic_interface,row_major=1,cwrap=1)
    
    # for now must be after lapack_expand (don't know why)
    generic_interface = process_ignore_info(generic_interface,format='C')

    generic_interface = process_flip_dims(generic_interface, format='C')

    # must be last
    generic_interface = process_return_info(generic_interface,format='C')    
    module_def = interface_to_module(generic_interface,module_name, 
                                     include_files,sdir)    
    
    # module_def = f2py_hack(module_def)
    # a bit of a cluge here on the naming - should get new name
    # form rename_functions
    f = open(os.path.join(output_path,module_name+'.pyf'),'w')
    f.write(module_def)
    f.close()

def generate_cblas_level(sdir,level):
    f = open(os.path.join(sdir,'generic_blas%d.pyf' % level))
    generic_interface = f.read(-1)
    generic_interface, include_files = process_imports(generic_interface,sdir)
    if level > 1:
        generic_interface = insert_blas_order(generic_interface)
    generic_interface = insert_intent_c(generic_interface)
    generic_interface = rename_functions(generic_interface,'cblas_','')
    generic_interface = process_special_types(generic_interface,format='C')    
    interface = lapack_expand(generic_interface,row_major=1)
    return interface, include_files
    
def generate_cblas(sdir,output_path):
    print "generating cblas interface"
    module_name = 'cblas'
    interface, include_files = generate_cblas_level(sdir,1)
    a, b = generate_cblas_level(sdir,2)
    interface += '\n' + a
    include_files += b
    a, b = generate_cblas_level(sdir,3)
    interface += '\n' + a
    include_files += b
    module_def = interface_to_module(interface,module_name,include_files,sdir)

    f = open(os.path.join(output_path,module_name+'.pyf'),'w')
    f.write(module_def)
    f.close()
    f.close()
    
def generate_fblas(sdir,output_path):
    print "generating fblas interface"
    module_name = 'fblas'
    interface = ''
    for level in (1,2,3):
        print '\tgenerating blas%d' % level
        f = open(os.path.join(sdir,'generic_blas%d.pyf' % level))
        generic_interface = f.read(-1)
        generic_interface, include_files = process_imports(generic_interface,
                                                           sdir)
        generic_interface = process_special_types(generic_interface,
                                                  format='Fortran')
        generic_interface = lapack_expand(generic_interface,row_major=0)
        interface = interface + '\n' + generic_interface

    module_definition = interface_to_module(interface,module_name, 
                                            include_files,sdir)
    f = open(os.path.join(output_path,module_name+'.pyf'),'w')
    f.write(module_definition)
    f.close()
        
def process_all():
    # process the standard files.
    generate_flapack()
    generate_clapack()
    generate_fblas()
    generate_cblas()
    #generate_fscalapack()
    
def usage():
    s = "usage: python interface_gen.py file_name module_name\n" \
        "\n" \
        "    file_name   -- file containing generic description of\n" \
        "                   lapack interface\n" \
        "    module_name -- name of module to generate\n"
    print s

if __name__ == '__main__':
    import sys
    generate_clapack('.')
    """    
    if len(sys.argv) == 2 and sys.argv[1] == 'all':
        process_all()
    elif len(sys.argv) < 3:
        usage()
    else:          
        f = open(sys.argv[1])
        module_name = sys.argv[2]
        generic_interface = f.read(-1)
        generic_interface, include_files = process_imports(generic_interface)
        all_subroutines(generic_interface)
        interface = lapack_expand(generic_interface)
        interface = rename_functions(interface,'AA','BB')
        module_defintion = interface_to_module(interface,module_name, include_files)
        f = open(module_name+'.pyf','w')
        f.write(module_definition)
        f.close()
    """