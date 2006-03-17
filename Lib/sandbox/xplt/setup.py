## Automatically adapted for scipy Oct 31, 2005 by

#!/usr/bin/env python

import os
import sys
from distutils import dir_util
from distutils.sysconfig import get_python_lib


def get_special_dirs(plat):
    if plat in ['aix4', 'aix5','sunos5']:
        return ['/usr/local/lib']
    elif plat in ['linux2','cygwin']:
        return ['/usr/lib']
    return []


cygwin = 0
if sys.platform=='cygwin':
    cygwin = 1

macosx = 0
if sys.platform=='darwin':
    macosx = 1

for keyword in sys.argv:
    if keyword=='--x11':
        sys.argv.remove(keyword)
        cygwin = 0
        macosx = 0

windows = 0
if sys.platform=='win32':
    windows = 1

x11 = 0
if not (windows or cygwin or macosx):
    x11 = 1
if 'NO_XLIB' in os.environ:
    x11 = 0



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
              "src/gist/eps.c",
              "src/gist/style.c",
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
elif macosx:
    unixsource = ["src/play/unix/dir.c",
                  "src/play/unix/files.c",
                  "src/play/unix/pathnm.c",
                  "src/play/unix/timew.c",
                  "src/play/unix/slinks.c",
                  "src/play/unix/stdinit.c",
                  "src/play/unix/uevent.c",
                  "src/play/unix/uinbg.c",
                  "src/play/unix/usernm.c"]
elif not (windows):
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

if not (windows or cygwin or macosx):
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
elif macosx:
    macsource = ["src/play/mac/pscr.m",
                 "src/play/mac/pals.m",
                 "src/play/mac/text.m",
                 "src/play/mac/cell.m",
                 "src/play/mac/bitblt.m",
                 "src/play/mac/points.m",
                 "src/play/mac/cursors.m",
                 "src/play/mac/pwin.m",
                 "src/play/mac/clips.m",
                 "src/play/mac/pen.m",
                 "src/play/mac/color.m",
                 "src/play/mac/font.m"]

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



def getallparams(gistpath,local_path,config_path):
    from numpy.distutils.system_info import get_info
    x11_info = get_info('x11')
    extra_compile_args = ['-DGISTPATH="\\"' + gistpath + '\\""' ]

    extra_link_args = []
    if windows or cygwin:
        extra_compile_args.append("-DWINDOWS")
        extra_compile_args.append("-mwindows")
        extra_link_args.append("-mwindows")
        libraries = []
    else:
        libraries = x11_info.get('libraries',['X11'])
    if cygwin:
        extra_compile_args.append("-DCYGWIN")
    if macosx:
        extra_compile_args.append("-DMACOSX")
        extra_link_args.append('-framework')
        extra_link_args.append('Cocoa')


    include_dirs = ['src/gist', 'src/play', 'src/play/unix' ]

    library_dirs = [os.path.join(local_path,x) for x in ['.','src']]
    library_dirs.extend(x11_info.get('library_dirs',[]))
    library_dirs.extend(get_special_dirs(sys.platform))

    include_dirs = [os.path.join(local_path,x) for x in include_dirs]
    include_dirs.extend(x11_info.get('include_dirs',[]))

    if 1:
        inputfile = open(os.path.join(config_path,"Make.cfg"))
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
                if xinc and sys.platform not in ['cygwin','win32']:
                    # remove the -I
                    xinc = xinc[2:]
                    if xinc: include_dirs.append(xinc)
            if line[:5]=="XLIB=":
                xlib = line[5:-1] # removing \n
                if xlib and sys.platform not in ['cygwin','win32']:
                    # remove the -L
                    xlib = xlib[2:]
                    library_dirs.append(xlib)

    return include_dirs, library_dirs, libraries, \
               extra_compile_args, extra_link_args


def configuration(parent_package='',top_path=None):
    """
       This will install *.gs and *.gp files to
       'site-packages/scipy/xplt/gistdata'
    """
    from numpy.distutils.misc_util   import Configuration
    from config_pygist import config_pygist
    config = Configuration('xplt',parent_package, top_path)
    local_path = config.local_path

    all_playsource = [os.path.join('src','play','*','*.c'),
                      os.path.join('src','play','*.h')
                      ]

    gistpath = os.path.join(get_python_lib(1),config.path_in_package,"gistdata")
    gistpath = gistpath.replace("\\",r"\\\\")

    def get_playsource(extension,build_dir):
        if windows:
            playsource = winsource + allsource
        elif cygwin:
            playsource = unixsource + winsource + allsource
        elif macosx:
            playsource = unixsource + macsource + allsource
        else:
            playsource = unixsource + x11source + allsource
        sources = [os.path.join(local_path,n) for n in playsource]

        config_path = os.path.join(build_dir,'config_pygist')
        dir_util.mkpath(config_path)
        conf = config_pygist(local_path,config_path)
        # Look to see if compiler is set on command line and add it
        #    This is repeating code, but I'm not sure how to avoid it
        #    As this gets run before overall setup does.
        #    This is needed so that compiler can be over-ridden from the
        #    platform default in the configuration section of xplt.
        for arg in sys.argv[1:]:
            if arg[:11] == '--compiler=':
                conf.compiler = arg[11:]
                break
            if arg[:2] == '-c':
                conf.compiler = arg[2:]
                break
        # Generate Make.cfg and config.h:
        conf.run()

        include_dirs, library_dirs, libraries, \
                      extra_compile_args, extra_link_args \
                      = getallparams(gistpath,local_path,config_path)
        include_dirs.insert(0,os.path.dirname(conf.config_h))

        extension.include_dirs.extend(include_dirs)
        extension.library_dirs.extend(library_dirs)
        extension.libraries.extend(libraries)
        extension.extra_compile_args.extend(extra_compile_args)
        extension.extra_link_args.extend(extra_link_args)
        return sources




    gistC = os.path.join('pygist','gistCmodule.c')
    sources = gistsource
    sources = [gistC] + sources + [get_playsource]

    config.add_extension('gistC',
                         sources,
                         depends = ['src']
                         )
    config.add_extension('gistfuncs',
                         [os.path.join('pygist','gistfuncsmodule.c')])


    file_ext = ['*.gs','*.gp', '*.ps', '*.help']
    xplt_files = [os.path.join('gistdata',x) for x in file_ext]
    xplt_files += [os.path.join('src','g',x) for x in file_ext]

    config.add_data_dir('gistdata')
    config.add_data_dir((os.path.join(config.path_in_package,'gistdata'),
                         os.path.abspath(config.paths('src/g')[0])))

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
