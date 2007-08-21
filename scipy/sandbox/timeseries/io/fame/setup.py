__version__ = '1.0'
__revision__ = "$Revision: 37 $"
__date__     = '$Date: 2006-12-08 14:30:29 -0500 (Fri, 08 Dec 2006) $'

import os, sys
from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    nxheader = join(get_numpy_include_dirs()[0],'numpy',)

    famedir = os.getenv('FAME')
    if famedir is None:
        raise EnvironmentError("FAME environment variable not found")

    if sys.platform == 'win32': msvc_flags()

    fameheader = famedir
    confgr = Configuration(parent_package=parent_package,top_path=top_path)

    sources = join('src', 'cfame.c')
    libraries = "chli"
    library_dirs = [famedir, join(famedir, "demo/hli")]
    confgr.add_extension('cfame',
                         sources=[sources],
                         include_dirs=[nxheader, fameheader, library_dirs],
                         libraries = [libraries],
                         library_dirs = [library_dirs]
                         )
    return confgr
    
def msvc_flags():
    """/DWIN32 flag is required on windows for compiling FAME
C-hli code"""

    from distutils.msvccompiler import MSVCCompiler

    # remember old initialize
    old_MSVCCompiler_initialize = MSVCCompiler.initialize

    def fame_msvccompiler_initialize(self, *args, **kws):
         apply(old_MSVCCompiler_initialize, (self,) + args, kws)
         self.compile_options.extend(['/DWIN32'])

    # "Install" new initialize
    MSVCCompiler.initialize = fame_msvccompiler_initialize

if __name__ == "__main__":

    from numpy.distutils.core import setup
    config = configuration().todict() 
    setup(**config)
    