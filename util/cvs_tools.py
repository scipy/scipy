#! /usr/bin/python
import os, sys

def rmdir(dir,depth=0):
    import os
    path = os.path.abspath(dir)
    all_files = os.listdir(path)
    indent = '   ' * depth
    for i in all_files:        
        if not i == 'CVS':
            print indent, i
            if os.path.isdir(os.path.join(path,i)):
                rmdir(os.path.join(path,i),depth+1)                            
            else:    
                cmd = 'cd ' + path + ';rm -r ' + i + ';cvs rm ' + i + ';cd ..'
                print cmd
                os.system(cmd)
    
ignore = ['*.o', '*.pyc', '*.a', '*.so', 'core', '*.dll', '*.pyd']

from fnmatch import fnmatch

def allowed_file_type(file):
    for i in ignore:
        if fnmatch(file, i):
            return 0
    return 1

def adddir(dir,depth=1):
    import os
    abs_path = os.path.abspath(dir)
    path = dir
    all_files = os.listdir(path)
    all_files = filter(allowed_file_type, all_files)
    indent = '   ' * depth
    
    if os.path.isdir(dir):
        print indent[:-3], dir
        cmd = 'cvs add ' + dir
        os.system(cmd)

    for i in all_files:        
        if not i == 'CVS':
            print indent, i
            if os.path.isdir(os.path.join(path,i)):
                adddir(os.path.join(path,i),depth+1)                            
            else:    
                cmd = 'cd ' + abs_path + ';cvs add ' + i + ';cd ..'
                print cmd
                os.system(cmd)
    
