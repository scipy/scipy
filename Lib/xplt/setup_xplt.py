import os
from scipy_distutils.misc_util import get_path, default_config_dict
    
def configuration(parent_package=''):
    """
       gist only works with an X-windows server
       This will install *.gs and *.gp files to
       '%spython%s/site-packages/scipy/xplt' % (sys.prefix,sys.version[:3])
    """
    # Check for X11 libraries
    save_file = './saved_values.py'
    if not os.path.exists(save_file):
	save_file = '../saved_values.py'

    try:
        exec(open(save_file).read())
        try:
            X11 = X11
        except NameError:
            X11 = check_and_save()
    except IOError:
        X11 = check_and_save()
    if X11:
        config = default_config_dict()

        if parent_package:
            parent_package = parent_package + '.'    
        local_path = get_path(__name__)

        config['packages'].append(parent_package+'xplt')
        
        from scipy_distutils.core import Extension                
        sources = ['gistCmodule.c']
        sources = [os.path.join(local_path,x) for x in sources]
                                               
        ext = Extension(parent_package+'xplt.gistC',
                        sources,
                        include_dirs = ['/usr/include/X11'],
                        library_dirs = ['/usr/X11R6/lib'],
                        libraries = ['X11','m'])
        config['ext_modules'].append(ext)
        
        from glob import glob
        gist = glob(os.path.join(local_path,'gist','*.c'))
        # libraries are C static libraries
        config['libraries'].append(('gist',{'sources':gist,
                                   'macros':[('STDC_HEADERS',1)]}))
                                   
        file_ext = ['*.gs','*.gp', '*.ps', '*.help']
        xplt_files = [glob(os.path.join(local_path,x)) for x in file_ext]
	xplt_files = reduce(lambda x,y:x+y,xplt_files,[])
        xplt_path = os.path.join(local_path,'xplt')
        config['data_files'].extend( [(xplt_path,xplt_files)])
        
        return config

# placing the line X11=0 in saved_values will turn off xplt installation.
def check_and_save(file='saved_values.py'):
    import commands
    output = commands.getoutput('find /usr -name "libX11*" -print')
    X11 = (output != '')
    fid = open(file,'a')
    fid.write('X11 = %d\n' % X11)
    fid.close()
    return X11
