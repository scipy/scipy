#!/usr/bin/env python

import os
import sys
import time
from scipy_distutils.misc_util   import get_path, default_config_dict, dot_join
from scipy_distutils.misc_util   import dict_append
from scipy_distutils.core        import Extension
from scipy_distutils.system_info import x11_info
from distutils.command.config    import config
from distutils.sysconfig         import get_python_lib
#from distutils.core              import Extension

# borrowed from setup.py for pygist

local_path = get_path(__name__)

cygwin = (sys.platform == 'cygwin')
for keyword in sys.argv:
    if keyword == '--x11':
        sys.argv.remove(keyword)
        cygwin = 0

windows = (sys.platform == 'win32')

x11 = not (windows or cygwin)
if 'NO_XLIB' in os.environ:
    x11 = 0

run_config = ('config' in sys.argv)
print '%%%%%%%%%%%%%%%%%%%%%%%'
print run_config
print '%%%%%%%%%%%%%%%%%%%%%%%'

#------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------

class config_pygist (config):
    def _link (self, body,
               headers, include_dirs,
               libraries, library_dirs, lang):
        (src, obj) = self._compile(body, headers, include_dirs, lang)
        prog = os.path.splitext(os.path.basename(src))[0]
        self.compiler.link_executable([obj], prog,
                                      libraries=libraries,
                                      library_dirs=library_dirs)

        # This code is added to correct error
        if self.compiler.exe_extension is None:
            self.compiler.exe_extension = ''
        # End of code different from standard 
        prog = prog + self.compiler.exe_extension
        self.temp_files.append(prog)

        return (src, obj, prog)
    
    def run (self):
        self.configfile = open(os.path.join(local_path, "pygist","Make.cfg"),'w')
        self.configfile.write('# Make.cfg from setup.py script ' + time.ctime() + '\n')
        if not windows:
            self.configfile.write('#')
            for item in os.uname():
                self.configfile.write(' '+item)
            self.configfile.write('\n')

        self.config_toplevel()
        if not windows: self.config_unix()
        if x11: self.config_x11()
        self.configfile.close()
        print 'wrote pygist/Make.cfg'
#----------------------------------------------------------------------
    def config_toplevel(self):
        print "  ============= begin top level configuration ============="

        # check alternate libm for Alpha Linux (see play/unix/README.fpu)
        if not 'MATHLIB' in os.environ:
            self.mathlib = 'm'
            testcode = """\
/* check whether libm is broken */
#include <math.h>
int main(int argc, char *argv[])
{
  return exp(-720.) > 1.0;  /* typically an IEEE denormal */
}
"""
            if self.try_link(testcode,libraries=[self.mathlib]):
                if not self.try_run(testcode,libraries=[self.mathlib]):
                    if self.try_link(testcode, libraries=['cpml']):
                        self.mathlib="cpml"
                        print "WARNING - using -lcpml instead of -lm"
                    else:
                        print "WARNING - libm broken? see play/unix/README.fpu"
                        print "  if on Alpha Linux, rerun ./configure with CC='gcc -mieee'"
            else:
                raise "math library missing; rerun setup.py after setting the MATHLIB env variable"
        else:
            self.mathlib=os.environ['MATHLIB']
        self.configfile.write('MATHLIB=-l'+self.mathlib+'\n')
        # check exp10 presence, emulate otherwise
        testcode = """\
int main(int argc, char *argv[])
{
  double x=exp10(3.);
  return (x<999.999)||(x>1000.001);
}
"""
        if self.try_link(testcode,libraries=[self.mathlib]):
            print "using exp10 found in libm"
            self.configfile.write("NO_EXP10=\n")
        else:
            print "libm does not contain exp10, will emulate"
            self.configfile.write("NO_EXP10=-DNO_EXP10\n")
#----------------------------------------------------------------------
    def config_unix(self):
        # begin play/unix configuration
        print
        print "  ============= begin play/unix configuration ============="
        print
        os.chdir(os.path.join(local_path, 'src','play','unix'))
        configfile = open('config.h','w')
        configfile.write('/* config.h used during config.sh script */\n')
        configfile.write('#ifndef CONFIG_SCRIPT\n')
        configfile.write('# error destroy this config.h and rerun configure script\n')
        configfile.write('#endif\n')
        configfile.close()
        configfile = open('config.0h','w')
        configfile.write('/* config.h from config.sh script ' + time.ctime() + '\n')
        configfile.write(' *')
        for item in os.uname():
            configfile.write(' '+item)
        configfile.write('\n')
        configfile.write(' */\n')
        self.fatality = 0
        if not cygwin:
            self.find_time(configfile)
            self.find_wall_time(configfile)
            self.find_sigfpe(configfile)
        self.find_user_name(configfile)
        self.find_tiocgpgrp(configfile)
        self.find_cwd(configfile)
        self.find_dirent(configfile)
        self.find_poll(configfile)

        if self.fatality:
            print "*** at least one play/unix component could not be configured"
            print "*** see configuration notes in play/unix/README.cfg"
        else:
            configfile.close()
            os.rename("config.0h","config.h")
            print "wrote src/play/unix/config.h"

        os.chdir(os.path.join('..','..','..'))
        print
        print "  ============== end play/unix configuration =============="
        
    def find_time(self,configfile):
        # find CPU time function (getrusage is best if present)
        testcode = """\
/* check settings of: USE_GETRUSAGE USE_TIMES */
#define CONFIG_SCRIPT
#include "timeu.c"
int
main(int argc, char *argv[])
{
  double s;
  double t = p_cpu_secs(&s);
  return 0;
}
"""
        if self.try_link("#define USE_GETRUSAGE\n"+testcode,include_dirs=[".."]):
            print "using getrusage() (CPU timer)"
            configfile.write('#define USE_GETRUSAGE\n')
        elif self.try_link("#define USE_TIMES\n"+testcode,include_dirs=[".."]):
            print "using times() (CPU timer)"
            configfile.write('#define USE_TIMES\n')
        elif self.try_link(testcode,include_dirs=[".."]):
            print "fallback to clock(), getrusage() and times() missing (CPU timer)"
        else:
            print "FATAL getrusage(), times(), and clock() all missing (timeu.c)"
            self.fatality = 1

    def find_wall_time(self,configfile):
        # find wall time function (gettimeofday is best if present)
        testcode = """\
/* check settings of: USE_GETTIMEOFDAY */
#define CONFIG_SCRIPT
#include "timew.c"
int
main(int argc, char *argv[])
{
  double t = p_wall_secs();
  return 0;
}
"""
        if self.try_link("#define USE_GETTIMEOFDAY\n"+testcode,include_dirs=[".."]):
            print "using gettimeofday() (wall timer)"
            configfile.write('#define USE_GETTIMEOFDAY\n')
        elif self.try_link(testcode,include_dirs=[".."]):
            print "fallback to time()+difftime(), gettimeofday() missing (wall timer)"
        else:
            print "FATAL gettimeofday(), and time() or difftime() missing (timew.c)"
            self.fatality = 1

    def find_sigfpe(self,configfile):
#----------------------------------------------------------------------
# try to figure out how to get SIGFPE delivered
#----------------------------------------------------------------------
        if sys.platform=='cygwin':
            # No SIGFPE delivery on Windows.
            print 'hardwiring #define FPU_IGNORE'
            configfile.write('#define FPU_IGNORE\n')
            return

        testcode="""\
#define CONFIG_SCRIPT
#include "fputest.c"
"""

        testgcc="""\
#define CONFIG_SCRIPT
#include "config.c"
"""

        fpedef=""
        fpelib=""
        fpelibm=""
        if self.try_link("#define FPU_DIGITAL\n"+testcode):
            print "using FPU_DIGITAL (SIGFPE delivery)"
            configfile.write('#define FPU_DIGITAL\n')
            fpedef = "-DFPU_DIGITAL"
        elif self.try_link("#define FPU_AIX\n"+testcode):
            print "using FPU_AIX (SIGFPE delivery)"
            configfile.write('#define FPU_AIX\n')
            fpedef = "-DFPU_AIX"
        elif self.try_link("#define FPU_HPUX\n"+testcode):
            print "using FPU_HPUX (SIGFPE delivery)"
            configfile.write('#define FPU_HPUX\n')
            fpedef = "-DFPU_HPUX"
            fpelibm = self.mathlib
        elif self.try_link("#define FPU_SOLARIS\n"+testcode):
            print "using FPU_SOLARIS (SIGFPE delivery)"
            configfile.write('#define FPU_SOLARIS\n')
            fpedef="-DFPU_SOLARIS"
            # note this works under IRIX 6.3, while FPU_IRIX does not??
        elif self.try_link("#define FPU_SUN4\n"+testcode,libraries=[self.mathlib]):
            print "using FPU_SUN4 (-lm) (SIGFPE delivery)"
            configfile.write('#define FPU_SUN4\n')
            fpedef="-DFPU_SUN4"
            fpelibm=self.mathlib
        elif self.try_link("#define FPU_SUN4\n"+testcode,libraries=["sunmath"]):
            print "using FPU_SUN4 (-lsunmath) (SIGFPE delivery)"
            configfile.write('#define FPU_SUN4\n')
            fpedef="-DFPU_SUN4"
            fpelib="sunmath"
        elif self.try_link("#define FPU_IRIX\n"+testcode,libraries=["fpe"]):
            # FPU_SOLARIS seems to work better??
            print "using FPU_IRIX (SIGFPE delivery)"
            configfile.write('#define FPU_IRIX\n')
            fpedef="-DFPU_IRIX"
            fpelib="fpe"
        elif self.try_link("#define FPU_IRIX\n"+testcode):
            print "using FPU_IRIX (SIGFPE delivery), but no libfpe??"
            configfile.write('#define FPU_IRIX\n')
            fpedef="-DFPU_IRIX"
        elif self.try_link("#define FPU_MACOSX\n"+testcode):
            print "using FPU_MACOSX (SIGFPE delivery)"
            configfile.write('#define FPU_MACOSX\n')
            fpedef="-DFPU_MACOSX"
        elif self.try_compile("#define TEST_GCC\n"+testgcc):
            if self.try_link("#define FPU_ALPHA_LINUX\n" + testcode):
                print "using FPU_ALPHA_LINUX (SIGFPE delivery)"
                configfile.write('#define FPU_ALPHA_LINUX\n')
                fpedef="-DFPU_ALPHA_LINUX"
                print "...libm may be broken -- read play/unix/README.fpu for more"
                print "...fputest failure may not mean that pygist itself is broken"
                # CC="$CC -mfp-trap-mode=su -mtrap-precision=i"
            elif self.try_link("#define FPU_GCC_I86\n" + testcode):
                print "using FPU_GCC_I86 (SIGFPE delivery)"
                configfile.write('#define FPU_GCC_I86\n')
                fpedef="-DFPU_GCC_I86"
            elif self.try_link("#define FPU_GCC_SPARC\n" + testcode):
                print "using FPU_GCC_SPARC (SIGFPE delivery)"
                configfile.write('#define FPU_GCC_SPARC\n')
                fpedef="-DFPU_GCC_SPARC"
            elif self.try_link("#define FPU_GCC_M86K\n" + testcode):
                print "using FPU_GCC_M68K (SIGFPE delivery)"
                configfile.write('#define FPU_GCC_M68K\n')
                fpedef="-DFPU_GCC_M68K"
            elif self.try_link("#define FPU_GCC_POWERPC\n" + testcode):
                print "using FPU_GCC_POWERPC (SIGFPE delivery)"
                configfile.write('#define FPU_GCC_POWERPC\n')
                fpedef="-DFPU_GCC_POWERPC"
            elif self.try_link("#define FPU_GCC_ARM\n" + testcode):
                print "using FPU_GCC_ARM (SIGFPE delivery)"
                configfile.write('#define FPU_GCC_ARM\n')
                fpedef="-DFPU_GCC_ARM"
        elif self.try_link("#define FPU_GNU_FENV \n" + testcode, libraries=[self.mathlib]):
	    raise 'kutje'
            print "using FPU_GNU_FENV (SIGFPE delivery)"
            configfile.write('#define FPU_GNU_FENV\n')
            fpedef="-DFPU_GNU_FENV"
            fpelibm=self.mathlib
        elif self.try_link("#define FPU_UNICOS\n" + testcode, libraries=[self.mathlib]):
            print "using FPU_UNICOS (SIGFPE delivery)"
            self.configfile.write('#define FPU_UNICOS\n')
            fpedef="-DFPU_UNICOS"
            fpelibm=self.mathlib

        if "fpedef" in os.environ:
            if self.try_link("#define FPU_IGNORE\n" + testcode, libraries=[self.mathlib]):
                print "using FPU_IGNORE (SIGFPE delivery)"
                configfile.write('#define FPU_IGNORE\n')
                fpedef="-DFPU_IGNORE"
            else:
                print "FATAL unable to build SIGFPE fputest? (fputest.c, fpuset.c)"
                self.fatality=1

        if fpelib:
            self.configfile.write("FPELIB="+fpelib+"\n")
        else:
            self.configfile.write("FPELIB=\n")

        if fpelibm:
            self.configfile.write("FPELIBM="+fpelibm+"\n")
        else:
            self.configfile.write("FPELIBM=\n")

        if fpedef:
            # on IRIX be sure that TRAP_FPE environment variable is turned off
            if "TRAP_FPE" in os.environ: del os.environ["TRAP_FPE"]
            testcode = "#define " + fpedef[2:] + "\n" + testcode
            libraries = []
            if fpelib: libraries.append(fpelib)
            if fpelibm: libraries.append(fpelibm)
            if not self.try_run(testcode,libraries=libraries):
                print
                print "*************************WARNING***************************"
                print "*** play/unix configuration failed to get SIGFPE delivered"
                print "*** read the notes in play/unix/README.fpu"
                print "*************************WARNING***************************"
                print

    def find_user_name(self,configfile):
        # find function to get user name
        testcode = """\
/* check settings of: NO_PASSWD */
#define CONFIG_SCRIPT
#include "usernm.c"
int
main(int argc, char *argv[])
{
  int value;
  char *u = p_getuser(); \
  value = (u!=0);
  return value;
}
"""
        if self.try_link(testcode,include_dirs=[".."]):
            print "using POSIX getlogin(), getpwuid(), getuid() functions"
        elif self.try_link("#define USE_PASSWD\n"+testcode,include_dirs=[".."]):
            print "fallback to cuserid(), POSIX getlogin() family missing"
            configfile.write('#define NO_PASSWD\n')
        else:
            print "FATAL cuserid(), POSIX getlogin() family both missing (usernm.c)"
            self.fatality = 1

    def find_tiocgpgrp(self, configfile):
        # find function to get controlling terminal process group
        testcode = """\
/* check settings of: USE_TIOCGPGRP_IOCTL */
#define CONFIG_SCRIPT
#include "uinbg.c"
int
main(int argc, char *argv[])
{
  int value;
  value = u_in_background();
  return value;
}
"""
        if self.try_compile(testcode,include_dirs=[".."]):
            pass
        elif self.try_compile("#define USE_POSIX_GETPGRP\n"+testcode,include_dirs=[".."]):
            print "using strict POSIX getpgrp prototype"
            testcode = "#define USE_POSIX_GETPGRP\n" + testcode
            configfile.write('#define USE_POSIX_GETPGRP\n')
        if self.try_link(testcode,include_dirs=[".."]):
            print "using POSIX tcgetpgrp() function"
        elif self.try_link('#define USE_TIOCGPGRP_IOCTL=<sys/termios.h>\n'+testcode,include_dirs=[".."]):
            print "fallback to TIOCGPGRP in sys/termios.h, POSIX tcgetpgrp() missing"
            configfile.write('#define USE_TIOCGPGRP_IOCTL <sys/termios.h>\n')
        elif self.try_link('#define USE_TIOCGPGRP_IOCTL=<sgtty.h>\n'+testcode,include_dirs=[".."]):
            print "fallback to TIOCGPGRP in sgtty.h, POSIX tcgetpgrp() missing"
            configfile.write('#define USE_TIOCGPGRP_IOCTL <sgtty.h>\n')
        else:
            print "FATAL unable to find TIOCGPGRP ioctl header (uinbg.c)"
            print "  (you can patch config.0h by hand if you know header)"
            configfile.write('#define USE_TIOCGPGRP_IOCTL <???>\n')
            self.fatality = 1

    def find_cwd(self, configfile):
        # find function to get current working directory
        testcode = """\
/* check settings of: USE_GETWD */
#define CONFIG_SCRIPT
#include <unistd.h>
static char dirbuf[1024];
#ifdef USE_GETWD
#define getcwd(x,y) getwd(x)
#endif
int
main(int argc, char *argv[])
{
  int value;
  char *u = getcwd(dirbuf, 1024);
  value = (u!=0);
  return value;
}
"""
        if self.try_link(testcode):
            print "using POSIX getcwd() function"
        elif self.try_link("#define USE_GETWD\n"+testcode):
            print "fallback to getwd(), POSIX getcwd() missing"
            configfile.write('#define USE_GETWD\n')
        else:
            print "FATAL getcwd(), getwd() both missing (dir.c)"
            self.fatality=1

    def find_dirent(self, configfile):
        # find headers required to read directories
        testcode = """\
/* check settings of: DIRENT_HEADER USE_GETWD */
#define TEST_DIRENT
#define CONFIG_SCRIPT
#include "dir.c"
p_twkspc p_wkspc;
int
main(int argc, char *argv[])
{
  int value;
  p_dir *d = p_dopen("no/such/thing");
  char *l = p_dnext(d, &value);
  value = p_chdir(l) || p_rmdir(l) || p_mkdir(l);
  return value;
}
"""
        if self.try_link(testcode,include_dirs=[".."]):
            print "using POSIX dirent.h header for directory ops"
        elif self.try_link('#define DIRENT_HEADER=<sys/dir.h>\n'+testcode,include_dirs=[".."]):
            print "using sys/dir.h header for directory ops"
            configfile.write('#define DIRENT_HEADER <sys/dir.h>\n')
        elif self.try_link('#define DIRENT_HEADER=<sys/ndir.h>\n'+testcode,include_dirs=[".."]):
            print "using sys/ndir.h header for directory ops"
            configfile.write('#define DIRENT_HEADER <sys/ndir.h>')
        elif self.try_link('#define DIRENT_HEADER=<ndir.h>\n'+testcode,include_dirs=[".."]):
            print "using ndir.h header for directory ops"
            configfile.write('#define DIRENT_HEADER <ndir.h>')
        else:
            print "FATAL dirent.h, sys/dir.h, sys/ndir.h, ndir.h all missing (dir.c)"
            self.fatality=1

    def find_poll(self,configfile):
        # find headers and functions required for poll/select functionality
        testcode = """\
/* check settings of: USE_SYS_POLL_H USE_SELECT HAVE_SYS_SELECT_H
                      NO_SYS_TIME_H NEED_SELECT_PROTO */
#define CONFIG_SCRIPT
#define TEST_POLL
#include "uevent.c"
int
main(int argc, char *argv[])
{
  int p = u_poll(1000); \
  return 0;
}
"""
        if self.try_link(testcode,include_dirs=[".."]):
            print "using poll(), poll.h header"
        elif self.try_link('#define USE_SYS_POLL_H\n'+testcode,include_dirs=[".."]):
            print "using poll(), sys/poll.h header"
            configfile.write('#define USE_SYS_POLL_H\n')
        elif self.try_link('#define USE_SELECT\n'+'#define HAVE_SYS_SELECT_H\n'+testcode,include_dirs=[".."]):
            print "using select(), sys/select.h header"
            configfile.write('#define USE_SELECT\n')
            configfile.write('#define HAVE_SYS_SELECT_H\n')
        elif self.try_link('#define USE_SELECT\n'+'#define NEED_SELECT_PROTO\n'+testcode,include_dirs=[".."]):
            print "using select(), sys/time.h, sys/types.h headers"
            configfile.write('#define USE_SELECT\n')
            configfile.write('#define NEED_SELECT_PROTO\n')
        elif self.try_link('#define USE_SELECT\n'+'#define NO_SYS_TIME_H\n'+'#define NEED_SELECT_PROTO\n'+testcode,include_dirs=[".."]):
            print "using select(), time.h, sys/types.h headers"
            configfile.write('#define USE_SELECT\n')
            configfile.write('#define NO_SYS_TIME_H\n')
            configfile.write('#define NEED_SELECT_PROTO\n')
        elif  self.try_link('#define USE_SELECT\n'+testcode,include_dirs=[".."]):
            print "using select(), sys/time.h header"
            configfile.write('#define USE_SELECT')
        elif  self.try_link('#define USE_SELECT\n'+'#define NO_SYS_TIME_H\n'+testcode,include_dirs=[".."]):
            print "using select(), time.h header"
            configfile.write('#define USE_SELECT\n')
            configfile.write('#define NO_SYS_TIME_H\n')
        else:
            print "FATAL neither poll() nor select() usable? (uevent.c, upoll.c)"
            self.fatality=1
#----------------------------------------------------------------------
    def config_x11(self):
        print
        print "  ============= begin play/x11 configuration =============="
        print
        from string import replace
        self.fatality=0

        # figure out directories to compile and load with X11
        if 'X11BASE' in os.environ:
            X11BASE=os.environ['X11BASE']
        else:
            X11BASE="/no/suggested/x11dir"

        # directory list is from autoconf, except openwin promoted near top
        xlist = [X11BASE+"/include",
                 "/usr/X11R6/include",
                 "/usr/X11R5/include",
                 "/usr/X11R4/include",
                 "/usr/include/X11R6",
                 "/usr/include/X11R5",
                 "/usr/include/X11R4",
                 "/usr/openwin/include",
                 "/usr/openwin/share/include",
                 "/usr/local/X11R6/include",
                 "/usr/local/X11R5/include",
                 "/usr/local/X11R4/include",
                 "/usr/local/include/X11R6",
                 "/usr/local/include/X11R5",
                 "/usr/local/include/X11R4",
                 "/usr/X11/include",
                 "/usr/include/X11",
                 "/usr/local/X11/include"
                 "/usr/local/include/X11",
                 "/usr/X386/include",
                 "/usr/x386/include",
                 "/usr/XFree86/include/X11",
                 "/usr/include",
                 "/usr/local/include",
                 "/usr/unsupported/include",
                 "/usr/athena/include",
                 "/usr/local/x11r5/include",
                 "/usr/lpp/Xamples/include"]
        testcode = """\
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>
extern XVisualInfo xvi;
XVisualInfo xvi;             /* XVisualInfo declared in Xutil.h */
int main(int argc, char *argv[])
{
  Display *dpy = XOpenDisplay("nosuchserver:0.0");
  xvi.screen = XK_BackSpace;  /* XK_BackSpace defined in keysym.h */
  xvi.red_mask = XA_PRIMARY;  /* XA_PRIMARY defined in Xatom.h */
  xvi.depth = XC_crosshair;   /* XC_crosshair defined in cursorfont.h */
  XCloseDisplay(dpy);
  return 0;
}
"""
        xinc = ""
        xlib = ""
        xfound = 0
        if self.try_compile(testcode):
            xfound=1
        else:
            for d in xlist:
                if os.path.exists(d+'/X11/Xlib.h'):
                    if self.try_compile(testcode,include_dirs=[d]):
                        xfound = 1
                        xinc = d
                        xlib = replace(d,"include","lib")
                        break
        if not xfound:
            print "FATAL unable to find X11 includes (play/x11) $xinc"
            self.fatality=1
        else:
            if self.try_link(testcode,include_dirs=[xinc],libraries=['X11']):
                xlib=""
                xfound=2
            elif xlib and self.try_link(testcode,include_dirs=[xinc],library_dirs=[xlib],libraries=['X11']):
                xfound=2
            else:
                xlist = [replace(directory,'include','lib') for directory in xlist]
		if xinc:
                    for d in xlist:
                        if self.try_link(testcode,include_dirs=[xinc],library_dirs=[d],libraries=['X11']):
                            xlib = d
                            xfound=2
                            break
		else:
                    for d in xlist:
                        if self.try_link(testcode,library_dirs=[d],libraries=['X11']):
                            xlib = d
                            xfound=2
                            break
        if xfound:
            print "found X Window System, X11 headers and libraries"
	    if xinc:
            	print "  - using X11 header switch -I"+xinc
	    else:
            	print "  - using X11 header switch [none]"
            print "  - using X11 loader switch -L"+xlib
        else:
            print "FATAL unable to find X11 libraries (play/x11) $xlib"
            self.fatality=1
        if xinc:
	    self.configfile.write("XINC=-I"+xinc+"\n")
	else:
	    self.configfile.write("XINC=\n")
        self.configfile.write("XLIB=-L"+xlib+"\n")

        print "appended to ../../Make.cfg"
        print
        print "  ============== end play/x11 configuration ==============="


def get_special_dirs(plat):
    if plat in ['aix4', 'aix5','sunos5']:
        return ['/usr/local/lib']
    elif plat in ['linux2','cygwin']:
        return ['/usr/lib']
    return []

gistsource = ["src/gist/gist.c",
              "src/gist/tick.c",
              "src/gist/tick60.c",
              "src/gist/engine.c",
              "src/gist/gtext.c",
              "src/gist/draw.c",
              "src/gist/draw0.c",
              "src/gist/clip.c",
              "src/gist/gread.c",
              "src/gist/gcntr.c",
              "src/gist/hlevel.c",
              "src/gist/ps.c",
              "src/gist/cgm.c",
              "src/gist/xfancy.c",
              "src/gist/xbasic.c"]

if cygwin:
    unixsource = ["src/play/unix/dir.c",
                  "src/play/unix/files.c",
                  "src/play/unix/pathnm.c",
                  "src/play/unix/slinks.c",
                  "src/play/unix/stdinit.c",
                  "src/play/unix/uevent.c",
                  "src/play/unix/uinbg.c",
                  "src/play/unix/usernm.c"]
elif not windows:
    unixsource = ["src/play/unix/dir.c",
                  "src/play/unix/files.c",
                  "src/play/unix/fpuset.c",
                  "src/play/unix/pathnm.c",
                  "src/play/unix/timew.c",
                  "src/play/unix/uevent.c",
                  "src/play/unix/ugetc.c",
                  "src/play/unix/umain.c",
                  "src/play/unix/usernm.c",
                  "src/play/unix/slinks.c"]

if not (windows or cygwin):
    x11source = ["src/play/x11/colors.c",
                 "src/play/x11/connect.c",
                 "src/play/x11/cursors.c",
                 "src/play/x11/errors.c",
                 "src/play/x11/events.c",
                 "src/play/x11/fills.c",
                 "src/play/x11/fonts.c",
                 "src/play/x11/images.c",
                 "src/play/x11/lines.c",
                 "src/play/x11/pals.c",
                 "src/play/x11/pwin.c",
                 "src/play/x11/resource.c",
                 "src/play/x11/rgbread.c",
                 "src/play/x11/textout.c",
                 "src/play/x11/rect.c",
                 "src/play/x11/clips.c",
                 "src/play/x11/points.c"]

if windows:
    winsource = ["src/play/win/pscr.c",
                 "src/play/win/pals.c",
                 "src/play/win/ptext.c",
                 "src/play/win/pfill.c",
                 "src/play/win/pcell.c",
                 "src/play/win/pmin.c",
                 "src/play/win/plines.c",
                 "src/play/win/prect.c",
                 "src/play/win/points.c",
                 "src/play/win/cursors.c",
                 "src/play/win/pwin.c",
                 "src/play/win/timew.c",
                 "src/play/win/clips.c",
                 "src/play/win/getdc.c",
                 "src/play/win/files.c",
                 "src/play/win/usernm.c",
                 "src/play/win/pathnm.c"]
elif cygwin:
    winsource = ["src/play/win/pscr.c",
                 "src/play/win/pals.c",
                 "src/play/win/ptext.c",
                 "src/play/win/pfill.c",
                 "src/play/win/pcell.c",
                 "src/play/win/pmin.c",
                 "src/play/win/plines.c",
                 "src/play/win/prect.c",
                 "src/play/win/points.c",
                 "src/play/win/cursors.c",
                 "src/play/win/pwin.c",
                 "src/play/win/timew.c",
                 "src/play/win/clips.c",
                 "src/play/win/getdc.c"]

allsource = ["src/play/all/hash.c",
             "src/play/all/hash0.c",
             "src/play/all/mm.c",
             "src/play/all/alarms.c",
             "src/play/all/pstrcpy.c",
             "src/play/all/pstrncat.c",
             "src/play/all/p595.c",
             "src/play/all/bitrev.c",
             "src/play/all/bitlrot.c",
             "src/play/all/bitmrot.c"]

def getallparams(gistpath,local):
    if windows:
        extra_compile_args = ['-DGISTPATH="\\"' + gistpath + '\\""' ]
    else:
        extra_compile_args = ['-DGISTPATH="' + gistpath + '"' ]

    extra_link_args = []
    if windows or cygwin:
        extra_compile_args.append("-DWINDOWS")
        extra_compile_args.append("-mwindows")
        extra_link_args.append("-mwindows")
    if cygwin:
        extra_compile_args.append("-DCYGWIN")

    include_dirs = [ 'src/gist', 'src/play', 'src/play/unix' ]

    if windows or cygwin:
        libraries = []
    else:
        libraries = ['X11']

    library_dirs = [os.path.join(local,x) for x in ['.','src']]
    library_dirs.extend(get_special_dirs(sys.platform))

    include_dirs = [os.path.join(local,x) for x in include_dirs]

    if not run_config:
        inputfile = open(os.path.join(local,"pygist","Make.cfg"))
        lines = inputfile.readlines()
        inputfile.close()
        for line in lines:
            if line[:8]=="MATHLIB=":
                mathlib = line[8:-1] #removing the \n
                # remove the -l
                mathlib = mathlib[2:]
                libraries.append(mathlib)
            if line[:9]=="NO_EXP10=":
                no_exp10 = line[9:-1] # removing \n
                if no_exp10: extra_compile_args.append(no_exp10)
            if line[:5]=="XINC=":
                xinc = line[5:-1] # removing \n
                if xinc and not (windows or cygwin):
                    # remove the -I
                    xinc = xinc[2:]
                    if xinc: include_dirs.append(xinc)
            if line[:5]=="XLIB=":
                xlib = line[5:-1] # removing \n
                if xlib and not (windows or cygwin):
                    # remove the -L
                    xlib = xlib[2:]
                    library_dirs.append(xlib)    

    return include_dirs, library_dirs, libraries, \
               extra_compile_args, extra_link_args
    

def configuration(parent_package=''):
    """
       This will install *.gs and *.gp files to
       'site-packages/scipy/xplt/gistdata' 
    """

    package = 'xplt'
    print '***********************************************'
    print local_path,parent_package
    print '***********************************************'    
        
    xplt_path = os.path.join(parent_package,'xplt')
    config = default_config_dict(package,parent_package)

    if windows:
        playsource = winsource + allsource
    elif cygwin:
        playsource = unixsource + winsource + allsource
    else:
        playsource = unixsource + x11source + allsource

    gistpath = os.path.join(get_python_lib(1),xplt_path,"gistdata")
    gistpath = gistpath.replace('\\','\\\\\\\\')

    sources = ['pygist/gistCmodule.c'] + gistsource + playsource
    sources = [os.path.join(local_path,x) for x in sources]

    include_dirs, library_dirs, libraries, \
                  extra_compile_args, extra_link_args = getallparams(gistpath,\
                                                                     local_path)

    ext_arg = {'name':dot_join(parent_package,package,'gistC'),
               'sources':sources,
               'include_dirs':include_dirs,
               'library_dirs':library_dirs,
               'libraries':libraries,
               'extra_compile_args':extra_compile_args,
               'extra_link_args':extra_link_args}
    ext = Extension (**ext_arg)
    config['ext_modules'].append(ext)


    from glob import glob
    file_ext = ['*.gs','*.gp', '*.ps', '*.help']
    xplt_files = [glob(os.path.join(local_path,'gistdata',x)) for x in file_ext]
    xplt_files += [glob(os.path.join(local_path,'src','g',x)) for x in file_ext]
    xplt_files = reduce(lambda x,y:x+y,xplt_files,[])
    data_path = os.path.join(xplt_path,'gistdata')
    config['data_files'].extend( [(data_path,xplt_files)])

    config['cmdclass'] = {'config' : config_pygist}
    
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
