
import os
import sys
import time
from distutils.command.config    import config
from scipy_distutils.system_info import get_info

#------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------

class config_pygist(config):

    def __init__(self,local_path,config_path):
        from distutils.dist import Distribution
        config.__init__(self, Distribution())
        self.local_path = local_path
        self.config_path = config_path
        self.config_h = os.path.join(self.config_path, "config.h")

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
        fn = os.path.join(self.config_path, "Make.cfg")
        if os.path.isfile(fn) and os.path.isfile(self.config_h):
            print '*'*70
            print 'Files\n%s\n%s\n exist.' % (fn,self.config_h)
            print 'Skipping pygist configuration'\
                  ' (remove %s to force reconfiguration).' % fn
            print '*'*70
            return

        from distutils.ccompiler import new_compiler
        self.compiler = new_compiler(compiler=self.compiler,
                                     verbose=1)
        from distutils.sysconfig import customize_compiler
        customize_compiler(self.compiler)

        self.configfile = open(fn,'w')
        self.configfile.write('# Make.cfg from setup.py script ' + time.ctime() + '\n')
        if sys.platform != 'win32':
            self.configfile.write('#')
            for item in os.uname():
                self.configfile.write(' '+item)
            self.configfile.write('\n')

        self.config_toplevel()
        self.config_unix()
        self.config_x11()
        self.configfile.close()
        print 'wrote',fn

    if sys.version[:3]<'2.2':
        def try_run (self, body,
                     headers=None, include_dirs=None,
                     libraries=None, library_dirs=None,
                     lang="c"):
            """Try to compile, link to an executable, and run a program
            built from 'body' and 'headers'.  Return true on success, false
            otherwise.
            """
            from distutils.ccompiler import CompileError, LinkError,\
                 DistutilsExecError
            self._check_compiler()
            try:
                src,obj,exe=self._link(body, headers, include_dirs,
                                       libraries, library_dirs, lang)
                self.spawn([os.path.join('.',exe)])
                ok = 1
            except (CompileError, LinkError, DistutilsExecError):
                ok = 0

            self.announce(ok and "success!" or "failure.")
            self._clean()
            return ok

        def try_compile (self, body, headers=None, include_dirs=None, lang="c"):
            """Try to compile a source file built from 'body' and 'headers'.
            Return true on success, false otherwise.
            """
            from distutils.ccompiler import CompileError
            self._check_compiler()
            try:
                self._compile(body, headers, include_dirs, lang)
                ok = 1
            except CompileError:
                ok = 0

            self.announce(ok and "success!" or "failure.")
            self._clean()
            return ok


#----------------------------------------------------------------------
    def config_toplevel(self):
        print "  ============= begin top level configuration ============="

        testcode = """\
/* check whether libm is broken */
#include <math.h>
int main(int argc, char *argv[])
{
  return exp(-720.) > 1.0;  /* typically an IEEE denormal */
}
"""
        # First try without libm to workaround MSVC6.0 compiler.
        if self.try_link(testcode) and self.try_run(testcode):
            print "No libm needed."
            self.mathlib = None
        # check alternate libm for Alpha Linux (see play/unix/README.fpu)
        elif not os.environ.has_key('MATHLIB'):
            self.mathlib = 'm'
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
        libraries = []
        if self.mathlib:
            self.configfile.write('MATHLIB=-l'+self.mathlib+'\n')
            libraries.append(self.mathlib)
        # check exp10 presence, emulate otherwise
        testcode = """\
int main(int argc, char *argv[])
{
  double x=exp10(3.);
  return (x<999.999)||(x>1000.001);
}
"""
        if self.try_link(testcode,libraries=libraries):
            print "using exp10 found in libm"
            self.configfile.write("NO_EXP10=\n")
        else:
            print "libm does not contain exp10, will emulate"
            self.configfile.write("NO_EXP10=-DNO_EXP10\n")
#----------------------------------------------------------------------
    def config_unix(self):
        if sys.platform == 'win32':
            configfile = open(self.config_h,'w')
            configfile.write('\n')
            configfile.close()
            return
        # begin play/unix configuration
        print
        print "  ============= begin play/unix configuration ============="
        print
        configfile = open(self.config_h,'w')
        configfile.write('/* config.h used during config.sh script */\n')
        configfile.write('#ifndef CONFIG_SCRIPT\n')
        configfile.write('# error destroy this config.h and rerun configure script\n')
        configfile.write('#endif\n')
        configfile.close()
        configfile = open(self.config_h+'0','w')
        configfile.write('/* config.h from config.sh script ' + time.ctime() + '\n')
        configfile.write(' *')
        for item in os.uname():
            configfile.write(' '+item)
        configfile.write('\n')
        configfile.write(' */\n')
        self.fatality = 0

        self.unix_include_dirs = [os.path.abspath(os.path.dirname(self.config_h)),
                                  os.path.abspath(os.path.join(self.local_path,
                                                               'src','play'))]
        _olddir = os.getcwd()
        os.chdir(os.path.join(self.local_path, 'src','play','unix'))
        if sys.platform!='cygwin':
            self.find_time(configfile)
            self.find_wall_time(configfile)
            self.find_sigfpe(configfile)
        self.find_user_name(configfile)
        self.find_tiocgpgrp(configfile)
        self.find_cwd(configfile)
        self.find_dirent(configfile)
        self.find_poll(configfile)
        os.chdir(_olddir)

        if self.fatality:
            print "*** at least one play/unix component could not be configured"
            print "*** see configuration notes in play/unix/README.cfg"
        else:
            configfile.close()
            os.rename(self.config_h+'0',self.config_h)
            print "wrote config.h"


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
        if self.try_link("#define USE_GETRUSAGE\n"+testcode,
                         include_dirs=self.unix_include_dirs):
            print "using getrusage() (CPU timer)"
            configfile.write('#define USE_GETRUSAGE\n')
        elif self.try_link("#define USE_TIMES\n"+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using times() (CPU timer)"
            configfile.write('#define USE_TIMES\n')
        elif self.try_link(testcode,include_dirs=self.unix_include_dirs):
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
        if self.try_link("#define USE_GETTIMEOFDAY\n"+testcode,
                         include_dirs=self.unix_include_dirs):
            print "using gettimeofday() (wall timer)"
            configfile.write('#define USE_GETTIMEOFDAY\n')
        elif self.try_link(testcode,include_dirs=self.unix_include_dirs):
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
        if self.try_link("#define FPU_DIGITAL\n"+testcode,
                         include_dirs=self.unix_include_dirs):
            print "using FPU_DIGITAL (SIGFPE delivery)"
            configfile.write('#define FPU_DIGITAL\n')
            fpedef = "-DFPU_DIGITAL"
        elif self.try_link("#define FPU_AIX\n"+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using FPU_AIX (SIGFPE delivery)"
            configfile.write('#define FPU_AIX\n')
            fpedef = "-DFPU_AIX"
        elif self.try_link("#define FPU_HPUX\n"+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using FPU_HPUX (SIGFPE delivery)"
            configfile.write('#define FPU_HPUX\n')
            fpedef = "-DFPU_HPUX"
            fpelibm = self.mathlib
        elif self.try_link("#define FPU_SOLARIS\n"+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using FPU_SOLARIS (SIGFPE delivery)"
            configfile.write('#define FPU_SOLARIS\n')
            fpedef="-DFPU_SOLARIS"
            # note this works under IRIX 6.3, while FPU_IRIX does not??
        elif self.try_link("#define FPU_SUN4\n"+testcode,
                           libraries=[self.mathlib],
                           include_dirs=self.unix_include_dirs):
            print "using FPU_SUN4 (-lm) (SIGFPE delivery)"
            configfile.write('#define FPU_SUN4\n')
            fpedef="-DFPU_SUN4"
            fpelibm=self.mathlib
        elif self.try_link("#define FPU_SUN4\n"+testcode,
                           libraries=["sunmath"],
                           include_dirs=self.unix_include_dirs):
            print "using FPU_SUN4 (-lsunmath) (SIGFPE delivery)"
            configfile.write('#define FPU_SUN4\n')
            fpedef="-DFPU_SUN4"
            fpelib="sunmath"
        elif self.try_link("#define FPU_IRIX\n"+testcode,
                           libraries=["fpe"],
                           include_dirs=self.unix_include_dirs):
            # FPU_SOLARIS seems to work better??
            print "using FPU_IRIX (SIGFPE delivery)"
            configfile.write('#define FPU_IRIX\n')
            fpedef="-DFPU_IRIX"
            fpelib="fpe"
        elif self.try_link("#define FPU_IRIX\n"+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using FPU_IRIX (SIGFPE delivery), but no libfpe??"
            configfile.write('#define FPU_IRIX\n')
            fpedef="-DFPU_IRIX"
        elif self.try_link("#define FPU_MACOSX\n"+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using FPU_MACOSX (SIGFPE delivery)"
            configfile.write('#define FPU_MACOSX\n')
            fpedef="-DFPU_MACOSX"
        elif self.try_compile("#define TEST_GCC\n"+testgcc,
                              include_dirs=self.unix_include_dirs):
            if self.try_link("#define FPU_ALPHA_LINUX\n" + testcode,
                             include_dirs=self.unix_include_dirs):
                print "using FPU_ALPHA_LINUX (SIGFPE delivery)"
                configfile.write('#define FPU_ALPHA_LINUX\n')
                fpedef="-DFPU_ALPHA_LINUX"
                print "...libm may be broken -- read play/unix/README.fpu for more"
                print "...fputest failure may not mean that pygist itself is broken"
                # CC="$CC -mfp-trap-mode=su -mtrap-precision=i"
            elif self.try_link("#define FPU_GCC_I86\n" + testcode,
                               include_dirs=self.unix_include_dirs):
                print "using FPU_GCC_I86 (SIGFPE delivery)"
                configfile.write('#define FPU_GCC_I86\n')
                fpedef="-DFPU_GCC_I86"
            elif self.try_link("#define FPU_GCC_SPARC\n" + testcode,
                               include_dirs=self.unix_include_dirs):
                print "using FPU_GCC_SPARC (SIGFPE delivery)"
                configfile.write('#define FPU_GCC_SPARC\n')
                fpedef="-DFPU_GCC_SPARC"
            elif self.try_link("#define FPU_GCC_M86K\n" + testcode,
                               include_dirs=self.unix_include_dirs):
                print "using FPU_GCC_M68K (SIGFPE delivery)"
                configfile.write('#define FPU_GCC_M68K\n')
                fpedef="-DFPU_GCC_M68K"
            elif self.try_link("#define FPU_GCC_POWERPC\n" + testcode,
                               include_dirs=self.unix_include_dirs):
                print "using FPU_GCC_POWERPC (SIGFPE delivery)"
                configfile.write('#define FPU_GCC_POWERPC\n')
                fpedef="-DFPU_GCC_POWERPC"
            elif self.try_link("#define FPU_GCC_ARM\n" + testcode,
                               include_dirs=self.unix_include_dirs):
                print "using FPU_GCC_ARM (SIGFPE delivery)"
                configfile.write('#define FPU_GCC_ARM\n')
                fpedef="-DFPU_GCC_ARM"
        elif self.try_link("#define FPU_GNU_FENV \n" + testcode,
                           libraries=[self.mathlib],
                           include_dirs=self.unix_include_dirs):
	    raise 'kutje'
            print "using FPU_GNU_FENV (SIGFPE delivery)"
            configfile.write('#define FPU_GNU_FENV\n')
            fpedef="-DFPU_GNU_FENV"
            fpelibm=self.mathlib
        elif self.try_link("#define FPU_UNICOS\n" + testcode,
                           libraries=[self.mathlib],
                           include_dirs=self.unix_include_dirs):
            print "using FPU_UNICOS (SIGFPE delivery)"
            self.configfile.write('#define FPU_UNICOS\n')
            fpedef="-DFPU_UNICOS"
            fpelibm=self.mathlib

        if os.environ.has_key('fpedef'):
            if self.try_link("#define FPU_IGNORE\n" + testcode,
                             libraries=[self.mathlib],
                             include_dirs=self.unix_include_dirs):
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
            if os.environ.has_key('TRAP_FPE'): del os.environ["TRAP_FPE"]
            testcode = "#define " + fpedef[2:] + "\n" + testcode
            libraries = []
            if fpelib: libraries.append(fpelib)
            if fpelibm: libraries.append(fpelibm)
            if not self.try_run(testcode,
                                libraries=libraries,
                                include_dirs=self.unix_include_dirs):
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
        if self.try_link(testcode,include_dirs=self.unix_include_dirs):
            print "using POSIX getlogin(), getpwuid(), getuid() functions"
        elif self.try_link("#define USE_PASSWD\n"+testcode,
                           include_dirs=self.unix_include_dirs):
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
        if self.try_compile(testcode,include_dirs=self.unix_include_dirs):
            pass
        elif self.try_compile("#define USE_POSIX_GETPGRP\n"+testcode,
                              include_dirs=self.unix_include_dirs):
            print "using strict POSIX getpgrp prototype"
            testcode = "#define USE_POSIX_GETPGRP\n" + testcode
            configfile.write('#define USE_POSIX_GETPGRP\n')
        if self.try_link(testcode,include_dirs=self.unix_include_dirs):
            print "using POSIX tcgetpgrp() function"
        elif self.try_link('#define USE_TIOCGPGRP_IOCTL=<sys/termios.h>\n'+testcode,
                           include_dirs=self.unix_include_dirs):
            print "fallback to TIOCGPGRP in sys/termios.h, POSIX tcgetpgrp() missing"
            configfile.write('#define USE_TIOCGPGRP_IOCTL <sys/termios.h>\n')
        elif self.try_link('#define USE_TIOCGPGRP_IOCTL=<sgtty.h>\n'+testcode,
                           include_dirs=self.unix_include_dirs):
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
        if self.try_link(testcode,include_dirs=self.unix_include_dirs):
            print "using POSIX getcwd() function"
        elif self.try_link("#define USE_GETWD\n"+testcode,
                           include_dirs=self.unix_include_dirs):
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
        if self.try_link(testcode,include_dirs=self.unix_include_dirs):
            print "using POSIX dirent.h header for directory ops"
        elif self.try_link('#define DIRENT_HEADER=<sys/dir.h>\n'+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using sys/dir.h header for directory ops"
            configfile.write('#define DIRENT_HEADER <sys/dir.h>\n')
        elif self.try_link('#define DIRENT_HEADER=<sys/ndir.h>\n'+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using sys/ndir.h header for directory ops"
            configfile.write('#define DIRENT_HEADER <sys/ndir.h>')
        elif self.try_link('#define DIRENT_HEADER=<ndir.h>\n'+testcode,
                           include_dirs=self.unix_include_dirs):
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
        if self.try_link(testcode,include_dirs=self.unix_include_dirs):
            print "using poll(), poll.h header"
        elif self.try_link('#define USE_SYS_POLL_H\n'+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using poll(), sys/poll.h header"
            configfile.write('#define USE_SYS_POLL_H\n')
        elif self.try_link('#define USE_SELECT\n'+'#define HAVE_SYS_SELECT_H\n'+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using select(), sys/select.h header"
            configfile.write('#define USE_SELECT\n')
            configfile.write('#define HAVE_SYS_SELECT_H\n')
        elif self.try_link('#define USE_SELECT\n'+'#define NEED_SELECT_PROTO\n'+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using select(), sys/time.h, sys/types.h headers"
            configfile.write('#define USE_SELECT\n')
            configfile.write('#define NEED_SELECT_PROTO\n')
        elif self.try_link('#define USE_SELECT\n'+'#define NO_SYS_TIME_H\n'+'#define NEED_SELECT_PROTO\n'+testcode,
                           include_dirs=self.unix_include_dirs):
            print "using select(), time.h, sys/types.h headers"
            configfile.write('#define USE_SELECT\n')
            configfile.write('#define NO_SYS_TIME_H\n')
            configfile.write('#define NEED_SELECT_PROTO\n')
        elif  self.try_link('#define USE_SELECT\n'+testcode,
                            include_dirs=self.unix_include_dirs):
            print "using select(), sys/time.h header"
            configfile.write('#define USE_SELECT')
        elif  self.try_link('#define USE_SELECT\n'+'#define NO_SYS_TIME_H\n'+testcode,
                            include_dirs=self.unix_include_dirs):
            print "using select(), time.h header"
            configfile.write('#define USE_SELECT\n')
            configfile.write('#define NO_SYS_TIME_H\n')
        else:
            print "FATAL neither poll() nor select() usable? (uevent.c, upoll.c)"
            self.fatality=1
#----------------------------------------------------------------------
    def config_x11(self):
        if sys.platform in ['win32','cygwin']:
            return
        if os.environ.has_key('NO_XLIB'):
            return
        print
        print "  ============= begin play/x11 configuration =============="
        print
        from string import replace
        self.fatality=0

        x11_info = get_info('x11')
        if x11_info:
            print "found X Window System, X11 headers and libraries"
            xinc = x11_info.get('include_dirs','')
            xinc = ' '.join(['-I'+d for d in xinc])
            xlib = x11_info.get('library_dirs','')
            xlib = ' '.join(['-L'+d for d in xlib])
	    if xinc:
            	print "  - using X11 header switch "+xinc
            print "  - using X11 loader switch "+xlib
            self.configfile.write("XINC="+xinc+"\n")
            self.configfile.write("XLIB="+xlib+"\n")
            print "appended to Make.cfg"
        else:
            print "FATAL unable to find X11 includes (play/x11) $xinc"
            self.fatality=1

        print
        print "  ============== end play/x11 configuration ==============="
        return
        # figure out directories to compile and load with X11
        X11BASE=os.environ.get('X11BASE','/no/suggested/x11dir')

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

        print "appended to Make.cfg"
        print
        print "  ============== end play/x11 configuration ==============="
