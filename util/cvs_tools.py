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
    

def adddir(dir,depth=0):
    import os
    path = os.path.abspath(dir)
    all_files = os.listdir(path)
    indent = '   ' * depth
    for i in all_files:        
        if not i == 'CVS':
            print indent, i
            if os.path.isdir(os.path.join(path,i)):
                adddir(os.path.join(path,i),depth+1)                            
            else:    
                cmd = 'cd ' + path + ';cvs add ' + i + ';cd ..'
                print cmd
                os.system(cmd)
    
