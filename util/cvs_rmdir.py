#! /usr/bin/python
import os, sys

def cvs_rmdir(dir):
    import os
    path = os.path.abspath(dir)
    all_files = os.listdir(path)
    
    for i in all_files:
        if not i == 'CVS':
            if os.path.isdir(i):
                cvs_rmdir(i)                            
            else:    
                cmd = 'cd ' + path + ';rm -r ' + i + ';cvs rm ' + i + ';cd ..'
                print cmd
                os.system(cmd)
    

if __name__ == '__main__':
    cvs_rmdir(sys.argv[1])