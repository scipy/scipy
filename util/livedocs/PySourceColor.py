"""
#############################################################################
# PySourceColor.py
#############################################################################
# A python source to colorized html converter.
# Hacked by M.E.Farmer Jr. 2004
# Python license
#############################################################################
# Now supports two types of output markup:
# - HTML markup does not create w3c valid html, but it works on every
#   browser i've tried so far.(I.E.,Mozilla/Firefox,Opera,wxHTML).
# - CSS markup is w3c valid html 4.01 strict,
#   but will not render correctly everywhere!
#   Too bad CSS is not supported on all browsers.
#############################################################################
# Features:
# -Can seperate and colorize:
#    12 types of strings
#    2 comment types
#    numbers
#    operators
#    class / name
#    def / name
#    decorator / name
#    keywords
#    arguments class/def/deco
#    text
# -Eight colorschemes built-in:
#    null
#    mono
#    dark (default)
#    dark2
#    lite
#    idle
#    viewcvs
#    pythonwin
# -Two types of markup:
#    html (default)
#    css/html
# -Any combination of four text styles:
#    none (default)
#    bold
#    italic
#    underline
#############################################################################
# Example usage:
#############################################################################
# # import
# import PySourceColor as psc
# psc.convert('c:/Python22/PySourceColor.py', colors=psc.idle, show=1)
#----------------------------------------------------------------------------
# # from module import *
# from PySourceColor import *
# convert('c:/Python22/Lib', colors=lite, markup="css")
#----------------------------------------------------------------------------
# # How to use a custom colorscheme.
# from PySourceColor import *
# new = {
#    ERRORTOKEN:             ('bui','#FF8080',''),
#    DECORATOR_NAME:         ('','#AACBBC',''),
#    DECORATOR:              ('','#333333',''),
#    NAME:                   ('','#1133AA','#DDFF22'),
#    NUMBER:                 ('','#236676','#FF5555'),
#    OPERATOR:               ('b','#435345','#BBBB11'),
#    COMMENT:                ('','#545366','#AABBFF'),
#    DOUBLECOMMENT:          ('','#553455','#FF00FF'),
#    CLASS_NAME:             ('','#000000','#FFFFFF'),
#    DEF_NAME:               ('u','#897845','#000022'),
#    KEYWORD:                ('','#345345','#FFFF22'),
#    SINGLEQUOTE:            ('','#223344','#AADDCC'),
#    SINGLEQUOTE_R:          ('','#344522',''),
#    SINGLEQUOTE_U:          ('','#234234',''),
#    DOUBLEQUOTE:            ('','#334421',''),
#    DOUBLEQUOTE_R:          ('','#345345',''),
#    DOUBLEQUOTE_U:          ('','#678673',''),
#    TRIPLESINGLEQUOTE:      ('','#FFFFFF','#000000'),
#    TRIPLESINGLEQUOTE_R:    ('bu','#443256','#DDFFDA'),
#    TRIPLESINGLEQUOTE_U:    ('','#423454','#DDFFDA'),
#    TRIPLEDOUBLEQUOTE:      ('ib','#000000','#FFFFFF'),
#    TRIPLEDOUBLEQUOTE_R:    ('ub','#000000','#FFFFFF'),
#    TRIPLEDOUBLEQUOTE_U:    ('', '#CCAABB','#FFFAFF'),
#    PAGEBACKGROUND:         '#FFFAAA',
#    }
# if __name__ == '__main__':
#     import sys
#     convert(sys.argv[1], './css.html',colors=new, markup='css', show=1)
#     convert(sys.argv[1], './html.html',colors=new, markup='html', show=1)
#############################################################################
"""
__all__ = ['ERRORTOKEN','DECORATOR_NAME', 'DECORATOR', 'ARGS',
           'NAME', 'NUMBER', 'OPERATOR', 'COMMENT',
           'DOUBLECOMMENT', 'CLASS_NAME', 'DEF_NAME', 'KEYWORD',
           'SINGLEQUOTE','SINGLEQUOTE_R','SINGLEQUOTE_U','DOUBLEQUOTE',
           'DOUBLEQUOTE_R', 'DOUBLEQUOTE_U', 'TRIPLESINGLEQUOTE',
           'TRIPLESINGLEQUOTE_R', 'TRIPLESINGLEQUOTE_U', 'TRIPLEDOUBLEQUOTE',
           'TRIPLEDOUBLEQUOTE_R', 'TRIPLEDOUBLEQUOTE_U', 'PAGEBACKGROUND',
           'null', 'mono', 'lite', 'dark','dark2', 'pythonwin','idle',
           'viewcvs', 'Usage', 'cli', 'str2stdout', 'path2stdout', 'Parser',
           'str2file', 'str2html', 'path2file', 'path2html', 'convert',
           'walkdir', 'defaultColors', 'showpage']
__title__ = 'PySourceColor'
__version__ = "1.9.10"
__date__ = '27 November 2004'
__author__ = "M.E.Farmer Jr."
__credits__ = '''This was originally based on a python recipe
submitted by Jürgen Hermann to ASPN.
M.E.Farmer 2004
Python license
'''
import os
import cgi
import sys
import time
import glob
import getopt
import keyword
import token
import tokenize
import cStringIO
import traceback
import webbrowser

# Do not edit
NAME = token.NAME
NUMBER = token.NUMBER
COMMENT = tokenize.COMMENT
OPERATOR = token.OP
ERRORTOKEN = token.ERRORTOKEN
ARGS = token.NT_OFFSET + 1
DOUBLECOMMENT = token.NT_OFFSET + 2
CLASS_NAME = token.NT_OFFSET + 3
DEF_NAME = token.NT_OFFSET + 4
KEYWORD = token.NT_OFFSET + 5
SINGLEQUOTE = token.NT_OFFSET + 6
SINGLEQUOTE_R = token.NT_OFFSET + 7
SINGLEQUOTE_U = token.NT_OFFSET + 8
DOUBLEQUOTE = token.NT_OFFSET + 9
DOUBLEQUOTE_R = token.NT_OFFSET + 10
DOUBLEQUOTE_U = token.NT_OFFSET + 11
TRIPLESINGLEQUOTE = token.NT_OFFSET + 12
TRIPLESINGLEQUOTE_R = token.NT_OFFSET + 13
TRIPLESINGLEQUOTE_U = token.NT_OFFSET + 14
TRIPLEDOUBLEQUOTE = token.NT_OFFSET + 15
TRIPLEDOUBLEQUOTE_R = token.NT_OFFSET + 16
TRIPLEDOUBLEQUOTE_U = token.NT_OFFSET + 17
PAGEBACKGROUND = token.NT_OFFSET + 18
DECORATOR = token.NT_OFFSET + 19
DECORATOR_NAME = token.NT_OFFSET + 20

# Do not edit (markup classname lookup)
MARKUPDICT = {
        ERRORTOKEN:             'err',
        DECORATOR_NAME:         'decn',
        DECORATOR:              'dec',
        ARGS:                   'args',
        NAME:                   'name',
        NUMBER:                 'num',
        OPERATOR:               'op',
        COMMENT:                'comm',
        DOUBLECOMMENT:          'dcom',
        CLASS_NAME:             'clsn',
        DEF_NAME:               'defn',
        KEYWORD:                'key',
        SINGLEQUOTE:            'sq',
        SINGLEQUOTE_R:          'sqr',
        SINGLEQUOTE_U:          'squ',
        DOUBLEQUOTE:            'dq',
        DOUBLEQUOTE_R:          'dqr',
        DOUBLEQUOTE_U:          'dqu',
        TRIPLESINGLEQUOTE:      'tsq',
        TRIPLESINGLEQUOTE_R:    'tsqr',
        TRIPLESINGLEQUOTE_U:    'tsqu',
        TRIPLEDOUBLEQUOTE:      'tdq',
        TRIPLEDOUBLEQUOTE_R:    'tdqr',
        TRIPLEDOUBLEQUOTE_U:    'tdqu',
        }

######################################################################
# Edit colors and styles to taste
# Create your own scheme, just copy one below , rename and edit.
# fore and back color are rgb hex and must be specified. #RRGGBB
# Styles are optional: b=bold, i=italic, u=underline
# Cutom styles must at least define NAME, ERRORTOKEN, PAGEBACKGROUND,
# all missing elements will default to NAME.
######################################################################
# Copy null and use it as a starter colorscheme.
null = {# tokentype:            ('tags', 'textforecolor', 'textbackcolor')
        ERRORTOKEN:             ('','#FF8080',''),# Error token
        DECORATOR_NAME:         ('','#000000',''),# Decorator name
        DECORATOR:              ('','#000000',''),# @
        ARGS:                   ('','#000000',''),
        NAME:                   ('','#000000',''),# All other text
        NUMBER:                 ('','#000000',''),# 0->10
        OPERATOR:               ('','#000000',''),# ()<>=!.:;^>%, etc.
        COMMENT:                ('','#000000',''),# Single comment
        DOUBLECOMMENT:          ('','#000000',''),## Double comment
        CLASS_NAME:             ('','#000000',''),# Class name
        DEF_NAME:               ('','#000000',''),# Def name
        KEYWORD:                ('','#000000',''),# Python keywords
        SINGLEQUOTE:            ('','#000000',''),# 'SINGLEQUOTE'
        SINGLEQUOTE_R:          ('','#000000',''),# r'SINGLEQUOTE'
        SINGLEQUOTE_U:          ('','#000000',''),# u'SINGLEQUOTE'
        DOUBLEQUOTE:            ('','#000000',''),# "DOUBLEQUOTE"
        DOUBLEQUOTE_R:          ('','#000000',''),# r"DOUBLEQUOTE"
        DOUBLEQUOTE_U:          ('','#000000',''),# u"DOUBLEQUOTE"
        TRIPLESINGLEQUOTE:      ('','#000000',''),# '''TRIPLESINGLEQUOTE'''
        TRIPLESINGLEQUOTE_R:    ('','#000000',''),# r'''TRIPLESINGLEQUOTE'''
        TRIPLESINGLEQUOTE_U:    ('','#000000',''),# u'''TRIPLESINGLEQUOTE'''
        TRIPLEDOUBLEQUOTE:      ('','#000000',''),# """TRIPLEDOUBLEQUOTE"""
        TRIPLEDOUBLEQUOTE_R:    ('','#000000',''),# r"""TRIPLEDOUBLEQUOTE"""
        TRIPLEDOUBLEQUOTE_U:    ('','#000000',''),# u"""TRIPLEDOUBLEQUOTE"""
        PAGEBACKGROUND:         '#FFFFFF'# set the page background
        }

mono = {
        ERRORTOKEN:             ('','#FF8080',''),
        DECORATOR_NAME:         ('bu','#000000',''),
        DECORATOR:              ('b','#000000',''),
        ARGS:                   ('b','#000000',''),
        NAME:                   ('','#000000',''),
        NUMBER:                 ('b','#000000',''),
        OPERATOR:               ('b','#000000',''),
        COMMENT:                ('i','#000000',''),
        DOUBLECOMMENT:          ('b','#000000',''),
        CLASS_NAME:             ('bu','#000000',''),
        DEF_NAME:               ('b','#000000',''),
        KEYWORD:                ('b','#000000',''),
        SINGLEQUOTE:            ('','#000000',''),
        SINGLEQUOTE_R:          ('','#000000',''),
        SINGLEQUOTE_U:          ('','#000000',''),
        DOUBLEQUOTE:            ('','#000000',''),
        DOUBLEQUOTE_R:          ('','#000000',''),
        DOUBLEQUOTE_U:          ('','#000000',''),
        TRIPLESINGLEQUOTE:      ('','#000000',''),
        TRIPLESINGLEQUOTE_R:    ('','#000000',''),
        TRIPLESINGLEQUOTE_U:    ('','#000000',''),
        TRIPLEDOUBLEQUOTE:      ('i','#000000',''),
        TRIPLEDOUBLEQUOTE_R:    ('i','#000000',''),
        TRIPLEDOUBLEQUOTE_U:    ('i','#000000',''),
        PAGEBACKGROUND:         '#FFFFFF'
        }

dark = {
        ERRORTOKEN:             ('','#FF8080',''),
        DECORATOR_NAME:         ('b','#FFBBAA',''),
        DECORATOR:              ('b','#CC5511',''),
        ARGS:                   ('b','#CCCCDD',''),
        NAME:                   ('','#FFFFFF',''),
        NUMBER:                 ('','#FF0000',''),
        OPERATOR:               ('b','#FAF785',''),
        COMMENT:                ('','#45FCA0',''),
        DOUBLECOMMENT:          ('i','#A7C7A9',''),
        CLASS_NAME:             ('b','#B599FD',''),
        DEF_NAME:               ('b','#EBAE5C',''),
        KEYWORD:                ('b','#8680FF',''),
        SINGLEQUOTE:            ('','#F8BAFE',''),
        SINGLEQUOTE_R:          ('','#F8BAFE',''),
        SINGLEQUOTE_U:          ('','#F8BAFE',''),
        DOUBLEQUOTE:            ('','#FF80C0',''),
        DOUBLEQUOTE_R:          ('','#FF80C0',''),
        DOUBLEQUOTE_U:          ('','#FF80C0',''),
        TRIPLESINGLEQUOTE:      ('','#FF9595',''),
        TRIPLESINGLEQUOTE_R:    ('','#FF9595',''),
        TRIPLESINGLEQUOTE_U:    ('','#FF9595',''),
        TRIPLEDOUBLEQUOTE:      ('','#B3FFFF',''),
        TRIPLEDOUBLEQUOTE_R:    ('','#B3FFFF',''),
        TRIPLEDOUBLEQUOTE_U:    ('','#B3FFFF',''),
        PAGEBACKGROUND:         '#000000'
        }

dark2 = {
        ERRORTOKEN:             ('','#FF8080',''),
        DECORATOR_NAME:         ('b','#FFBBAA',''),
        DECORATOR:              ('b','#CC5511',''),
        ARGS:                   ('b','#FFFFFF',''),
        NAME:                   ('','#c0c0c0',''),
        NUMBER:                 ('b','#00FF00',''),
        OPERATOR:               ('b','#FF090F',''),
        COMMENT:                ('i','#F0F709','#844200'),
        DOUBLECOMMENT:          ('i','#F0F709','#844200'),
        CLASS_NAME:             ('b','#7E58C7',''),
        DEF_NAME:               ('b','#FF8040',''),
        KEYWORD:                ('b','#4726E1',''),
        SINGLEQUOTE:            ('','#8080C0',''),
        SINGLEQUOTE_R:          ('','#8080C0',''),
        SINGLEQUOTE_U:          ('','#8080C0',''),
        DOUBLEQUOTE:            ('','#ADB9F1',''),
        DOUBLEQUOTE_R:          ('','#ADB9F1',''),
        DOUBLEQUOTE_U:          ('','#ADB9F1',''),
        TRIPLESINGLEQUOTE:      ('','#00C1C1',''),
        TRIPLESINGLEQUOTE_R:    ('','#00C1C1',''),
        TRIPLESINGLEQUOTE_U:    ('','#00C1C1',''),
        TRIPLEDOUBLEQUOTE:      ('','#33e3e3',''),
        TRIPLEDOUBLEQUOTE_R:    ('','#33e3e3',''),
        TRIPLEDOUBLEQUOTE_U:    ('','#33e3e3',''),
        PAGEBACKGROUND:         '#000000'
        }

lite = {
        ERRORTOKEN:             ('','#FF8080',''),
        DECORATOR_NAME:         ('b','#BB4422',''),
        DECORATOR:              ('b','#3333af',''),
        ARGS:                   ('','#000000',''),
        NAME:                   ('','#000000',''),
        NUMBER:                 ('','#FF2200',''),
        OPERATOR:               ('b','#303000',''),
        COMMENT:                ('','#007F00',''),
        DOUBLECOMMENT:          ('','#606060',''),
        CLASS_NAME:             ('','#0000FF',''),
        DEF_NAME:               ('b','#9C7A00',''),
        KEYWORD:                ('b','#0000AF',''),
        SINGLEQUOTE:            ('','#600080',''),
        SINGLEQUOTE_R:          ('','#600080',''),
        SINGLEQUOTE_U:          ('','#600080',''),
        DOUBLEQUOTE:            ('','#A0008A',''),
        DOUBLEQUOTE_R:          ('','#A0008A',''),
        DOUBLEQUOTE_U:          ('','#A0008A',''),
        TRIPLESINGLEQUOTE:      ('','#337799',''),
        TRIPLESINGLEQUOTE_R:    ('','#337799',''),
        TRIPLESINGLEQUOTE_U:    ('','#337799',''),
        TRIPLEDOUBLEQUOTE:      ('','#1188AA',''),
        TRIPLEDOUBLEQUOTE_R:    ('','#1188AA',''),
        TRIPLEDOUBLEQUOTE_U:    ('','#1188AA',''),
        PAGEBACKGROUND:         '#FFFFFF'
        }

idle = {
        ERRORTOKEN:             ('','#FF8080',''),
        DECORATOR_NAME:         ('','#900090',''),
        DECORATOR:              ('','#000000',''),
        NAME:                   ('','#000000',''),
        NUMBER:                 ('','#000000',''),
        OPERATOR:               ('','#000000',''),
        COMMENT:                ('','#DD0000',''),
        DOUBLECOMMENT:          ('','#DD0000',''),
        CLASS_NAME:             ('','#0000FF',''),
        DEF_NAME:               ('','#0000FF',''),
        KEYWORD:                ('','#FF7700',''),
        SINGLEQUOTE:            ('','#00AA00',''),
        SINGLEQUOTE_R:          ('','#00AA00',''),
        SINGLEQUOTE_U:          ('','#00AA00',''),
        DOUBLEQUOTE:            ('','#00AA00',''),
        DOUBLEQUOTE_R:          ('','#00AA00',''),
        DOUBLEQUOTE_U:          ('','#00AA00',''),
        TRIPLESINGLEQUOTE:      ('','#00AA00',''),
        TRIPLESINGLEQUOTE_R:    ('','#00AA00',''),
        TRIPLESINGLEQUOTE_U:    ('','#00AA00',''),
        TRIPLEDOUBLEQUOTE:      ('','#00AA00',''),
        TRIPLEDOUBLEQUOTE_R:    ('','#00AA00',''),
        TRIPLEDOUBLEQUOTE_U:    ('','#00AA00',''),
        PAGEBACKGROUND:         '#FFFFFF'
        }

pythonwin = {
        ERRORTOKEN:             ('','#FF8080',''),
        DECORATOR_NAME:         ('b','#303030',''),
        DECORATOR:              ('b','#DD0080',''),
        ARGS:                   ('','#000000',''),
        NAME:                   ('','#303030',''),
        NUMBER:                 ('','#008080',''),
        OPERATOR:               ('','#000000',''),
        COMMENT:                ('','#007F00',''),
        DOUBLECOMMENT:          ('','#7F7F7F',''),
        CLASS_NAME:             ('b','#0000FF',''),
        DEF_NAME:               ('b','#007F7F',''),
        KEYWORD:                ('b','#000080',''),
        SINGLEQUOTE:            ('','#808000',''),
        SINGLEQUOTE_R:          ('','#808000',''),
        SINGLEQUOTE_U:          ('','#808000',''),
        DOUBLEQUOTE:            ('','#808000',''),
        DOUBLEQUOTE_R:          ('','#808000',''),
        DOUBLEQUOTE_U:          ('','#808000',''),
        TRIPLESINGLEQUOTE:      ('','#808000',''),
        TRIPLESINGLEQUOTE_R:    ('','#808000',''),
        TRIPLESINGLEQUOTE_U:    ('','#808000',''),
        TRIPLEDOUBLEQUOTE:      ('','#808000',''),
        TRIPLEDOUBLEQUOTE_R:    ('','#808000',''),
        TRIPLEDOUBLEQUOTE_U:    ('','#808000',''),
        PAGEBACKGROUND:         '#FFFFFF'
        }

viewcvs = {
        ERRORTOKEN:             ('b','#FF8080',''),
        DECORATOR_NAME:         ('','#000000',''),
        DECORATOR:              ('','#000000',''),
        ARGS:                   ('','#000000',''),
        NAME:                   ('','#000000',''),
        NUMBER:                 ('','#000000',''),
        OPERATOR:               ('','#000000',''),
        COMMENT:                ('i','#b22222',''),
        DOUBLECOMMENT:          ('i','#b22222',''),
        CLASS_NAME:             ('','#000000',''),
        DEF_NAME:               ('b','#0000ff',''),
        KEYWORD:                ('b','#a020f0',''),
        SINGLEQUOTE:            ('b','#bc8f8f',''),
        SINGLEQUOTE_R:          ('b','#bc8f8f',''),
        SINGLEQUOTE_U:          ('b','#bc8f8f',''),
        DOUBLEQUOTE:            ('b','#bc8f8f',''),
        DOUBLEQUOTE_R:          ('b','#bc8f8f',''),
        DOUBLEQUOTE_U:          ('b','#bc8f8f',''),
        TRIPLESINGLEQUOTE:      ('b','#bc8f8f',''),
        TRIPLESINGLEQUOTE_R:    ('b','#bc8f8f',''),
        TRIPLESINGLEQUOTE_U:    ('b','#bc8f8f',''),
        TRIPLEDOUBLEQUOTE:      ('b','#bc8f8f',''),
        TRIPLEDOUBLEQUOTE_R:    ('b','#bc8f8f',''),
        TRIPLEDOUBLEQUOTE_U:    ('b','#bc8f8f',''),
        PAGEBACKGROUND:         '#FFFFFF'
        }

defaultColors = dark

def Usage():
    """
 -----------------------------------------------------------------------------
  PySourceColor.py ver: %s
 -----------------------------------------------------------------------------
  Module summary:
     This module is designed to colorize python source code.
         Standalone:
             This module will work from the command line with options.
             This module will work with redirected stdio.
         Imported:
             This module can be imported and used directly in your code.
 -----------------------------------------------------------------------------
  Command line options:
     -h, --help
         Optional-> Display this help message.
     -t, --test
         Optional-> Will ignore all others flags but  --profile
             test all schemes and markup combinations
     -p, --profile
         Optional-> Works only with --test or -t
             runs profile.py and makes the test work in quiet mode.
     -i, --in, --input
         Use any of these for the current dir (.,cwd)
         Input can be file or dir.
         Input from stdin use one of the following (-,stdin)
         If stdin is used as input stdout is output unless specified.
     -o, --out, --output
         Optional-> output dir for the colorized source.
             default: output dir is the input dir.
         To output html to stdout use one of the following (-,stdout)
         Stdout can be used without stdin if you give a file as input.
     -c, --color
         Optional-> null, mono, dark, dark2, lite, idle, pythonwin, viewcvs
             default: dark
     -s, --show
         Optional-> Show webpage after creation.
             default: no show
     -m, --markup
         Optional-> html, css
             default: HTML
 -----------------------------------------------------------------------------
  Option usage:
   # Test and show pages
      python PySourceColor.py -t
   # Test and only show profile results
      python PySAourceColor.py -t -p
   # Colorize all .py,.pyw files in cwdir you can also use: (.,cwd)
      python PySourceColor.py -i .
   # Using long options w/ =
      python PySourceColor.py --in=c:/myDir/my.py --color=lite --show
   # Using short options w/out =
      python PySourceColor.py -i c:/myDir/  -c idle -m css
   # Using any mix
      python PySourceColor.py --in . -o=c:/myDir --show
 -----------------------------------------------------------------------------
  Stdio usage:
   # Stdio using no options
      python PySourceColor.py < c:/MyFile.py >> c:/tmp/MyFile.html
   # Using stdin alone automatically uses stdout for output: (stdin,-)
      python PySourceColor.py -i- < c:/MyFile.py >> c:/tmp/myfile.html
   # Stdout can also be written to directly from a file instead of stdin
      python PySourceColor.py -i c:/MyFile.py -m css -o- >> c:/tmp/myfile.html
   # Stdin can be used as input , but output can still be specified
      python PySourceColor.py -i- -o c:/pydoc.py.html -s < c:/Python22/my.py
 _____________________________________________________________________________
 """
    print Usage.__doc__% (__version__)
    sys.exit(1)

###################################################### Command line interface

def cli():
    """Handle command line args and redirections"""
    try:
        # try to get command line args
        opts, args = getopt.getopt(sys.argv[1:],
              "hsqtpi:o:c:m:",["help", "show", "quiet", "profile", "test",
              "input=", "output=", "color=", "markup="])
    except getopt.GetoptError:
        # on error print help information and exit:
        Usage()

    # init some names
    input = None
    output = None
    colorscheme = None
    markup = 'html'
    show = 0
    quiet = 0
    test = 0
    profile = 0
    # if we have args then process them
    for o, a in opts:
        if o in ["-h", "--help"]:
            Usage()
            sys.exit()
        if o in ["-o", "--output", "--out"]:
            output = a
        if o in ["-i", "--input", "--in"]:
            input = a
            if input in [".", "cwd"]:
                input = os.getcwd()
        if o in ["-s", "--show"]:
            show = 1
        if o in ["-q", "--quiet"]:
            quiet = 1
        if o in ["-t", "--test"]:
            test = 1
        if o in ["-m", "--markup"]:
            markup = str(a)
        if o in ["-p", "--profile"]:
            profile = 1
        if o in ["-c", "--color"]:
            try:
                colorscheme = globals().get(a.lower())
            except:
                traceback.print_exc()
                Usage()
    if test:
        if profile:
            import profile
            profile.run('_test(show=%s, quiet=%s)'%(show,quiet))
        else:
            # Parse this script in every possible colorscheme and markup
            _test(show,quiet)
    elif input in [None, "-", "stdin"] or output in ["-", "stdout"]:
        # determine if we are going to use stdio
        if input not in [None, "-", "stdin"]:
            if os.path.isfile(input) :
                path2stdout(input)
            else:
                raise PathError, 'File does not exists!'
        else:
            try:
                if sys.stdin.isatty():
                    raise InputError, 'Please check input!'
                else:
                    if output in [None,"-","stdout"]:
                        str2stdout(sys.stdin.read())
                    else:
                        str2file(sys.stdin.read(), output, show)
            except:
                traceback.print_exc()
                Usage()
    else:
        if os.path.exists(input):
            # if there was at least an input given we can proceed
            convert(input, output, colorscheme, show, markup, quiet)
        else:
            raise PathError, 'File does not exists!'
            Usage()

######################################################### Simple markup tests

def _test(show=0, quiet=0):
    """Test the parser and most of the functions.

       There are 15 _test total(seven colorschemes in two diffrent markups,
       and a str2file test. Most functions are tested by this.
    """
    fi = sys.argv[0]
    if not fi.endswith('.exe'):# Do not test if frozen as an archive
        path2file(fi, '/tmp/null.html', null, show=show, quiet=quiet)
        path2file(fi, '/tmp/null_css.html', null, show=show,
                  markup='css', quiet=quiet)
        path2file(fi, '/tmp/mono.html', mono, show=show, quiet=quiet)
        path2file(fi, '/tmp/mono_css.html', mono, show=show,
                  markup='css', quiet=quiet)
        path2file(fi, '/tmp/lite.html', lite, show=show, quiet=quiet)
        path2file(fi, '/tmp/lite_css.html', lite, show=show,
                  markup='css', quiet=quiet)
        path2file(fi, '/tmp/dark.html', dark, show=show, quiet=quiet)
        path2file(fi, '/tmp/dark_css.html', dark, show=show,
                  markup='css', quiet=quiet)
        path2file(fi, '/tmp/dark2.html', dark2, show=show, quiet=quiet)
        path2file(fi, '/tmp/dark2_css.html', dark2, show=show,
                  markup='css', quiet=quiet)
        path2file(fi, '/tmp/idle.html', idle, show=show, quiet=quiet)
        path2file(fi, '/tmp/idle_css.html', idle, show=show,
                  markup='css', quiet=quiet)
        path2file(fi, '/tmp/viewcvs.html', viewcvs, show=show, quiet=quiet)
        path2file(fi, '/tmp/viewcvs_css.html', viewcvs, show=show,
                  markup='css', quiet=quiet)
        path2file(fi, '/tmp/pythonwin.html', pythonwin, show=show,
                  quiet=quiet)
        path2file(fi, '/tmp/pythonwin_css.html', pythonwin, show=show,
                  markup='css', quiet=quiet)
        teststr=r'''"""This is a test of decorators and other things"""
@whatever(arg,arg2)
def my(arg,arg2):
   """This is a docstring"""
   print 'fgfdgfgfdgfdgfdgfdg'
# Comment before decorators
@staticmethod## Double comments are kewl
def hgfghf():
   """this is a docstring"""
   """Second doc -->

===================
Just for giggles
===================
## comments in a second docstring
>>> This = 123
"""
   ##double comments------------>
   w = {'sdfdfdf':'1','dfdfdf':'2'}
   p = r'dfgfg\nfgfdgfgfdfgfg\dfgfgfdgfdgfd\dfgfd\g'
   j = p.split(" ")
   Ur"""asdasdasdsadsadsadsadsadsadsadsa"""
   U"""asdasdsadsadsadsadsad"""
   r'sdfsdfdsfdsfdsfdfdf'
   rU"sdfdsfdsfdsfdsfdfd"
   u'gdfgfdgfgfdgfgfgfdg'
   import os
   os.path.split(r'/tmp/null')
   e = u'ügdg'+r"""ytryt="""
   g = """dfgfdgfdgdfgdfgdfgfdg"""+"ghgfhgfhghghh"
   # who said that
   d = 'fgfdgfdgfdgfdgfdgfgfgfgfdgfdgfgfdgfdgfdgfddffdgfdgfdg\
dfgdfgfdgfdgdfgdfgfgfgfdgfgdgdfgfg'
   print d'''
        htmlPath = os.path.abspath('/tmp/strtest.html')
        str2file(teststr, htmlPath, colors=dark, show=show)
        _printinfo("  wrote %s" % htmlPath, quiet)
    else:
        Usage()

####################################################### User level funtctions

def str2stdout(sourcestring, colors=None, form=None):
    """Converts a code(string) to colorized HTML. Writes to stdout.

       form='code',or'snip' (for "<pre>yourcode</pre>" only)
       colors=null,mono,lite,dark,idle,or pythonwin
    """
    Parser(sourcestring, colors).format(form)

def path2stdout(sourcepath, colors=None, form=None):
    """Converts code(file) to colorized HTML. Writes to stdout.

       form='code',or'snip' (for "<pre>yourcode</pre>" only)
       colors=null,mono,lite,dark,idle,or pythonwin
    """
    sourcestring = open(sourcepath).read()
    Parser(sourcestring, colors, title=sourcepath).format(form)

def str2html(sourcestring, colors=None, form=None):
    """Converts a code(string) to colorized HTML. Returns an HTML string.

       form='code',or'snip' (for "<pre>yourcode</pre>" only)
       colors=null,mono,lite,dark,idle,or pythonwin
    """
    stringIO = cStringIO.StringIO()
    Parser(sourcestring, colors, out=stringIO, markup='html').format(form)
    stringIO.seek(0)
    return stringIO.read()

def str2file(sourcestring, outfile, colors=None, show=0):
    """Converts a string to a file.

       makes no attempt at correcting bad pathnames
    """
    html = str2html(sourcestring, colors)
    f = open(outfile,'wt')
    f.writelines(html)
    f.close()
    if show:
        showpage(outfile)

def path2html(sourcepath, colors=None, form=None):
    """Converts code(file) to colorized HTML. Returns an HTML string.

       form='code',or'snip' (for "<pre>yourcode</pre>" only)
       colors=null,mono,lite,dark,idle,or pythonwin
    """
    stringIO = cStringIO.StringIO()
    sourcestring = open(sourcepath).read()
    Parser(sourcestring, colors, title=sourcepath, out=stringIO,
           markup='html').format(form)
    stringIO.seek(0)
    return stringIO.read()

def convert(source, outdir=None, colors=None,
              show=0, markup='html', quiet=0):
    """Takes a file or dir as input and places the html in the outdir.

       If outdir is none it defaults to the input dir
    """
    c=0
    # If it is a filename then path2file
    if not os.path.isdir(source):
        if os.path.isfile(source):
            c+=1
            path2file(source, outdir, colors, show, markup, quiet)
        else:
            raise PathError, 'File does not exists!'
    # If we pass in a dir we need to walkdir for files.
    # Then we need to colorize them with path2file
    else:
        fileList = walkdir(source)
        if fileList != None:
            # make sure outdir is a dir
            if outdir != None:
                if os.path.splitext(outdir)[1] != '':
                    outdir = os.path.split(outdir)[0]
            for item in fileList:
                c+=1
                path2file(item, outdir, colors, show, markup, quiet)
            _printinfo('Completed colorizing %s files.'%str(c), quiet)
        else:
            _printinfo("No files to convert in dir.", quiet)

def path2file(sourcePath, out=None, colors=None, show=0,
                markup='html', quiet=0):
    """ Converts python source to html file"""
    # If no outdir is given we use the sourcePath
    if out == None:#this is a guess
        htmlPath = sourcePath + '.html'
    else:
        # If we do give an out_dir, and it does
        # not exist , it will be created.
        if os.path.splitext(out)[1] == '':
            if not os.path.isdir(out):
                os.makedirs(out)
            sourceName = os.path.basename(sourcePath)
            htmlPath = os.path.join(out,sourceName)+'.html'
        # If we do give an out_name, and its dir does
        # not exist , it will be created.
        else:
            outdir = os.path.split(out)[0]
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
            htmlPath = out
    htmlPath = os.path.abspath(htmlPath)
    # Open the text and do the parsing.
    source = open(sourcePath).read()
    Parser(source, colors, sourcePath, open(htmlPath, 'wt'),
           markup).format()
    _printinfo("  wrote %s" % htmlPath, quiet)
    if show:
        # load HTML page into the default web browser.
        showpage(htmlPath)
    return htmlPath

def walkdir(dir):
    """Return a list of .py and .pyw files from a given directory.

       This function can be written as a generator Python 2.3, or a genexp
       in Python 2.4. But 2.2 and 2.1 would be left out....
    """
    # Get a list of files that match *.py*
    GLOB_PATTERN = os.path.join(dir, "*.[p][y]*")
    pathlist = glob.glob(GLOB_PATTERN)
    # Now filter out all but py and pyw
    filterlist = [x for x in pathlist
                        if x.endswith('.py')
                        or x.endswith('.pyw')]
    if filterlist != []:
        # if we have a list send it
        return filterlist
    else:
        return None

def showpage(path):
    """Helper function to open webpages"""
    try:
        webbrowser.open_new(os.path.abspath(path))
    except:
        traceback.print_exc()

def _printinfo(message, quiet):
    """Helper to print messages"""
    if not quiet:
        print message

########################################################### Custom Exceptions

class PySourceColorError(Exception):
    pass# base for custom error

class PathError(PySourceColorError):
    pass# custom error

class InputError(PySourceColorError):
    pass# custom error

########################################################## Python code parser

class Parser:

    """MoinMoin python parser heavily chopped :)"""

    def __init__(self, raw, colors=None, title='', out=sys.stdout,
                   markup='html', header=0, footer=0):
        """Store the source text & set some flags"""
        if colors == None:
            colors = defaultColors
        self.raw = raw.expandtabs().strip()
        self.title = os.path.basename(title)
        self.out = out
        self.argFlag = 0
        self.classFlag = 0
        self.defFlag = 0
        self.decoratorFlag = 0
        self.markup = markup.upper()
        self.colors = colors
        self.header = header
        self.footer = footer

    def format(self, form=None):
        """Parse and send the colorized source"""
        if form in ('snip','code'):
            self.addEnds = 0
        else:
            self.addEnds = 1
        # Store line offsets in self.lines
        self.lines = [0, 0]
        pos = 0

        # Gather lines
        while 1:
            pos = self.raw.find('\n', pos) + 1
            if not pos: break
            self.lines.append(pos)
        self.lines.append(len(self.raw))

        # Wrap text in a filelike object
        self.pos = 0
        text = cStringIO.StringIO(self.raw)

        # Markup start
        if self.addEnds:
            self._doPageStart()
        else:
            self._doSnippetStart()

        ## Tokenize calls the __call__
        ## function for each token till done.
        # Parse the source and write out the results.
        try:
            tokenize.tokenize(text.readline, self)
        except tokenize.TokenError, ex:
            msg = ex[0]
            line = ex[1][0]
            self.out.write("<h3>ERROR: %s</h3>%s\n"%
                            (msg, self.raw[self.lines[line]:]))
            traceback.print_exc()

        # Markup end
        if self.addEnds:
            self._doPageEnd()
        else:
            self._doSnippetEnd()

    def __call__(self, toktype, toktext, (srow,scol), (erow,ecol), line):
        """Token handler. Order is important do not rearrange."""
        style = ''
        # Calculate new positions
        oldpos = self.pos
        newpos = self.lines[srow] + scol
        self.pos = newpos + len(toktext)

        # Handle newlines
        if toktype in (token.NEWLINE, tokenize.NL):
            self.out.write('\n')
            return

        # Send the original whitespace, if needed
        if newpos > oldpos:
            self.out.write(self.raw[oldpos:newpos])

        # Skip indenting tokens
        if toktype in (token.INDENT, token.DEDENT):
            self.pos = newpos
            return

        # Look for operators
        if token.LPAR <= toktype and toktype <= token.OP:
            # Trap decorators py2.4>
            if toktext == '@':
                toktype = DECORATOR
                # Set a flag if this was the decorator start so
                # the decorator name and arguments can be identified
                self.decoratorFlag = self.argFlag = 1
            else:
                # Find the start for arguments
                if toktext == '(' and self.argFlag:
                    self.argFlag = 2
                # Find the end for arguments
                elif toktext == ':':
                    self.argFlag = 0
                toktype = token.OP

        # Look for keywords
        elif toktype == token.NAME and keyword.iskeyword(toktext):
            toktype = KEYWORD
            # Set a flag if this was the class / def start so
            # the class / def name and arguments can be identified
            if toktext =='class':
                self.classFlag = self.argFlag = 1
            elif toktext == 'def':
                self.defFlag = self.argFlag = 1

        # Look for class, def, decorator name
        elif self.classFlag or self.defFlag or self.decoratorFlag:
            if self.classFlag:
                self.classFlag = 0
                toktype = CLASS_NAME
            elif self.defFlag:
                self.defFlag = 0
                toktype = DEF_NAME
            elif self.decoratorFlag:
                self.decoratorFlag = 0
                toktype = DECORATOR_NAME

        # Seperate the diffrent comment types
        # Order of evaluation is important do not change.
        elif toktype == token.STRING:
            text = toktext.lower()
            # TRIPLE DOUBLE QUOTE's
            if (text[:3].lower() == '"""'):
                toktype = TRIPLEDOUBLEQUOTE
            elif (text[:4] == 'r"""'):
                toktype = TRIPLEDOUBLEQUOTE_R
            elif (text[:4] == 'u"""' or
                   text[:5] == 'ur"""' or
                   text[:5] == 'ru"""'):
                toktype = TRIPLEDOUBLEQUOTE_U
            # DOUBLE QUOTE's
            elif (text[:1] == '"'):
                toktype = DOUBLEQUOTE
            elif (text[:2] == 'r"'):
                toktype = DOUBLEQUOTE_R
            elif (text[:2] == 'u"' or
                   text[:3] == 'ur"' or
                   text[:3] == 'ru"'):
                toktype = DOUBLEQUOTE_U
            # TRIPLE SINGLE QUOTE's
            elif (text[:3] == "'''"):
                toktype = TRIPLESINGLEQUOTE
            elif (text[:4] == "r'''"):
                toktype = TRIPLESINGLEQUOTE_R
            elif (text[:4] == "u'''" or
                   text[:5] == "ur'''" or
                   text[:5] == "ru'''"):
                toktype = TRIPLESINGLEQUOTE_U
            # SINGLE QUOTE's
            elif (text[:1] == "'"):
                toktype = SINGLEQUOTE
            elif (text[:2] == "r'"):
                toktype = SINGLEQUOTE_R
            elif (text[:2] == "u'" or
                   text[:3] == "ur'" or
                   text[:3] == "ru'" ):
                toktype = SINGLEQUOTE_U

        # Seperate the diffrent comment types
        elif toktype == tokenize.COMMENT:
            if toktext[:2] == "##":
                toktype = DOUBLECOMMENT

        # Seperate errors from decorators
        elif toktype == token.ERRORTOKEN:
            # trap decorators...<py2.4
            if toktext == '@':
                toktype = DECORATOR
                # Set a flag if this was the decorator start so
                # the decorator name and arguments can be identified
                self.decoratorFlag = self.argFlag = 1

        # Seperate args from names
        elif (self.argFlag == 2 and
              toktype == token.NAME and
              toktext != 'None'):
            toktype = ARGS

        # Send text for any markup
        getattr(self, '_send%sText'%(self.markup))(toktype, toktext)
        return

    ################################################################# Helpers

    def _doSnippetStart(self):
        # Start of html snippet
        self.out.write('<pre>\n')

    def _doSnippetEnd(self):
        # Start of html snippet
        self.out.write('</pre>\n')

    ######################################################## markup selectors

    def _doPageStart(self):
        getattr(self, '_do%sStart'%(self.markup))()

    def _doPageHeader(self):
        getattr(self, '_do%sHeader'%(self.markup))()

    def _doPageFooter(self):
        getattr(self, '_do%sFooter'%(self.markup))()

    def _doPageEnd(self):
        getattr(self, '_do%sEnd'%(self.markup))()

    ################################################### color/style retrieval

    def _getTags(self, key):
        # style tags
        return self.colors.get(key, self.colors[NAME])[0]

    def _getForeColor(self, key):
        # get text foreground color, if not set to black
        color = self.colors.get(key, self.colors[NAME])[1]
        if color[:1] != '#':
            color = '#000000'
        return color

    def _getBackColor(self, key):
        # get text background color
        return self.colors.get(key, self.colors[NAME])[2]

    def _getPageColor(self):
        # get page background color
        return self.colors.get(PAGEBACKGROUND, '#FFFFFF')

    def _getStyle(self, key):
        # get the token style from the color dictionary
        return self.colors.get(key, self.colors[NAME])

    def _getMarkupClass(self, key):
        # get the markup class name from the markup dictionary
        return MARKUPDICT.get(key, MARKUPDICT[NAME])

    def _setDocumentCreatedBy(self):
        return '<!--This document created by %s ver.%s on: %s-->\n'%(
                  __title__,__version__,time.ctime())

    ################################################### HTML markup functions

    def _doHTMLStart(self):
        # Start of html page
        self.out.write('<!DOCTYPE html PUBLIC \
"-//W3C//DTD HTML 4.01//EN">\n')
        self.out.write('<html><head><title>%s</title>\n'%(self.title))
        self.out.write(self._setDocumentCreatedBy())
        self.out.write('<meta http-equiv="Content-Type" \
content="text/html;charset=iso-8859-1">\n')
        # Get background
        self.out.write('</head><body bgcolor="%s">\n'%self._getPageColor())
        # Write a little info at the top.
        if self.header:
            self._doPageHeader()
        self.out.write('<pre>')

    def _sendHTMLText(self, toktype, toktext):
        # If it is an error set a red box around the bad tokens
        # older browsers will ignore it
        if toktype == ERRORTOKEN:
            style = ' style="border: solid 1.5pt #FF0000;"'
        else:
            style = ''
        # Get styles
        tags, color = self.colors.get(toktype, self.colors[NAME])[:2]#
        tagstart=[]
        tagend=[]
        # check for styles and set them if needed.
        if 'b' in tags:#Bold
            tagstart.append('<b>')
            tagend.append('</b>')
        if 'i' in tags:#Italics
            tagstart.append('<i>')
            tagend.append('</i>')
        if 'u' in tags:#Underline
            tagstart.append('<u>')
            tagend.append('</u>')
        # HTML tags should be paired like so : <b><i><u>Doh!</u></i></b>
        tagend.reverse()
        starttag="".join(tagstart)
        endtag="".join(tagend)
        sendtext = cgi.escape(toktext)
        # send text
        ## Output optimization
        # skip font tag if black text, but styles will still be sent. (b,u,i)
        if color !='#000000':
            startfont = '<font color="%s"%s>'%(color, style)
            endfont = '</font>'
        else:
            startfont, endfont = ('','')
        self.out.write(''.join([startfont,starttag,sendtext,endtag,endfont]))
        return

    def _doHTMLHeader(self):
        # Optional
        color = self._getForeColor(token.NAME)
        self.out.write(' <h3><b><font color="%s">#%s %s</font></b></h3>\n'%
                        (color, self.title, time.ctime()))

    def _doHTMLFooter(self):
        # Optional
        color = self._getForeColor(token.NAME)
        self.out.write(' <h3><b><font color="%s">#%s %s</font></b></h3>\n'%
                        (color, self.title,time.ctime()))

    def _doHTMLEnd(self):
        # End of html page
        self.out.write('</pre>\n')
        # Write a little info at the bottom
        if self.footer:
            self._doPageFooter()
        self.out.write('</body></html>\n')

    #################################################### CSS markup functions

    def _getCSSStyle(self, key):
        # Get the tags and colors from the dictionary
        tags, forecolor, backcolor = self._getStyle(key)
        style=[]
        if tags:
            if 'b' in tags:# Bold
                style.append('font-weight:bold;')
            if 'i' in tags:# Italic
                style.append('font-style:italic;')
            if 'u' in tags:# Underline
                style.append('text-decoration:underline;')
        style.append('color:%s;'% forecolor)
        if backcolor:
            style.append('background-color:%s;'%backcolor)
        return (self._getMarkupClass(key),''.join(style))

    def _doCSSStart(self):
        # Start of css/html page
        self.out.write('<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN">')
        self.out.write('<html><head><title>%s</title>\n'%(self.title))
        self.out.write(self._setDocumentCreatedBy())
        self.out.write('<meta http-equiv="Content-Type" \
content="text/html;charset=iso-8859-1">\n')
        self.out.write('<style type="text/css">\n')
        # Get page background color and write styles ignore all but b,i,u
        self.out.write('<!--\n')
        self.out.write('body {background:%s;}\n'%self._getPageColor())
        # write out the various css styles
        for key in MARKUPDICT:
            if key == ERRORTOKEN:
                # If it is an errortoken set a red box around it
                self.out.write('.%s {border: solid 1.5pt #FF0000;%s}\n'
                                %self._getCSSStyle(ERRORTOKEN))
            else:
                # set the styles for all but errortokens
                self.out.write('.%s {%s}\n'%self._getCSSStyle(key))
        self.out.write('--></style>\n</head>\n<body>\n')
        # Write a little info at the top.
        if self.header:
            self._doPageHeader()
        self.out.write('<pre>')

    def _sendCSSText(self, toktype, toktext):
        # Get styles
        ##Output optimization
        # send text plain if black and no tags.
        tags, color = self.colors.get(toktype, self.colors[NAME])[:2]#
        if color == '#000000' and not tags:
            startspan, endspan = ('','')
        else:
            markupclass = MARKUPDICT.get(toktype, MARKUPDICT[NAME])#
            startspan = '<span class="%s">'%(markupclass)
            endspan = '</span>'
        sendtext = cgi.escape(toktext)
        self.out.write(''.join([startspan,sendtext,endspan]))
        return

    def _doCSSHeader(self):
        # Optional
        self.out.write('<h3><span class="name">#%s %s</span></h3>\n'%
                        (self.title, time.ctime()))

    def _doCSSFooter(self):
        # Optional
        self.out.write('<h3><span class="name">#%s %s</span></h3>\n'%
                        (self.title, time.ctime()))

    def _doCSSEnd(self):
        # End of css/html page
        self.out.write('</pre>\n')
        # Write a little info at the bottom
        if self.footer:
            self._doPageFooter()
        self.out.write('</body></html>\n')

#############################################################################

if __name__ == '__main__':
    cli()

#############################################################################
# 2004 M.E.Farmer Jr.
# Python license
