/*
 * BROWSER.C
 *
 * $Id$
 *
 * Main for GIST CGM viewer
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include <string.h>
#include "pstdlib.h"
#include "pstdio.h"

#include "xbasic.h"
#include "cgm.h"

/* OpenCGM, ReadCGM, CGMRelative, cgmLandscape, bg0fg1 in cgmin.h */
#include "cgmin.h"
/* GpEPSEngine declared in eps.h */
#include "eps.h"

extern long strtol(const char *, char **, int);

extern int cgmScaleToFit;
extern int gtDoEscapes;

/* List of commands in alphabetical order, 0 terminated */
static char *commandList[]= {
  "cgm", "display", "draw", "eps", "free", "help", "info", "open",
  "ps", "quit", "send", 0 };

/* Corresponding functions */
static int MakeCGM(int help);
static int MakeX(int help);
static int Draw(int help);
static int EPS(int help);
static int FreeAny(int help);
static int Help(int help);
static int Info(int help);
static int OpenIn(int help);
static int MakePS(int help);
static int Quit(int help);
static int Send(int help);
static int Special(int help);  /* must correspond to 0 in commandList */
static int nPrefix, cSuffix;   /* arguments for Special */
typedef int (*Command)(int help);
static Command CommandList[]= {
  &MakeCGM, &MakeX, &Draw, &EPS, &FreeAny, &Help, &Info, &OpenIn,
  &MakePS, &Quit, &Send, &Special };
static int GetCommand(char *token);

int amBatch= 0;

int warningCount= 0;
int no_warnings= 0;

static int CreateCGM(int device, char *name);
static int CreatePS(int device, char *name);
static int CreateX(int device, char *name);

static int GetPageGroup(char *token);
static int CheckEOL(char *command);
static void DrawSend(int ds, char *token);
static int GetDeviceList(int ds, char *command);
static int FindDevice(void);
static int SaveName(int device, char *name);
static void DoSpecial(int nPrefix, int cSuffix);

static int xPrefix= 0;
static void HandleExpose(Engine *engine, Drauing *drauing, int xy[]);
static void HandleOther(Engine *engine, int k, int md);
static int HelpAndExit(void);
static int MessageAndExit(char *msg);

extern int on_idle(void);
extern int on_quit(void);
extern void on_stdin(char *line);
static int did_startup = 0;
extern int do_startup(void);

int defaultDPI;

char *outNames[8];
int outLengths[8];
int outDraw[8], outSend[8], outMark[8];

int nPageGroups= 0;
int mPage[32], nPage[32], sPage[32];
int n_active_groups = 0;

static int nOut, noDisplay, x_only;
static char *inName;

int
on_launch(int argc, char *argv[])
{
  int i, j;
  char *arg;

  nOut= 0;
  for (i=0 ; i<8 ; i++) {
    outNames[i]= 0;
    outEngines[i]= 0;
    outDraw[i]= outSend[i]= outLengths[i]= 0;
  }
  noDisplay= amBatch= x_only= 0;
  defaultDPI= 100;
  inName= 0;

  p_quitter(&on_quit);
  p_idler(&on_idle);
  p_stdinit(&on_stdin);
  g_stdout = p_stdout;

  for (i=1 ; i<argc ; i++) {
    arg= argv[i];

    if (arg[0]=='-') {
      int fileType= -1;
      arg++;
      if (strcmp(arg, "cgm")==0) {
        fileType= 0;
        i++;
        if (i>=argc) return MessageAndExit("Missing file or display name");
        arg= argv[i];
      } else if (strcmp(arg, "ps")==0) {
        fileType= 1;
        i++;
        if (i>=argc) return MessageAndExit("Missing file or display name");
        arg= argv[i];
      } else if (strcmp(arg, "in")==0) {
        i++;
        if (i>=argc) return MessageAndExit("Missing file or display name");
        if (inName) return HelpAndExit();
        else inName= argv[i];
      } else if (strcmp(arg, "display")==0 || strcmp(arg, "d")==0) {
        fileType= 2;
        i++;
        if (i>=argc) return MessageAndExit("Missing file or display name");
        arg= argv[i];
      } else if (strcmp(arg, "f")==0) {
        amBatch= 1;
        fileType= 1;
        arg= "*stdout*";
      } else if (strcmp(arg, "nd")==0) noDisplay= 1;
      else if (strcmp(arg, "b")==0) amBatch= 1;
      else if (strcmp(arg, "nowarn")==0) no_warnings= 1;
      else if (strcmp(arg, "geometry")==0) {
        char *suffix;
        int w=0,h=0;
        i++;
        if (i>=argc) MessageAndExit("Missing geometry");
        arg = argv[i];
        w = (int)strtol(arg, &suffix, 10);
        if (suffix!=arg && *suffix=='x') {
          arg = suffix+1;
          h = (int)strtol(arg, &suffix, 10);
        }
        if (w < 10 || h < 10) MessageAndExit("Invalid geometry");
        gx75width = gx100width = w;
        gx75height = gx100height = h;
      } else if (strcmp(arg, "75")==0) defaultDPI= 75;
      else if (strcmp(arg, "100")==0) defaultDPI= 100;
      else if (strcmp(arg, "gks")==0) {
        gx75width= gx75height= 600;     /* 8x8 X window size */
        gx100width= gx100height= 800;   /* 8x8 X window size */
        cgmScaleToFit= 1;               /* 8x8 PostScript plotting area */
        gtDoEscapes= 0;
      } else if (strcmp(arg, "x")==0) x_only= 1;
      else if (strcmp(arg, "gks")==0) {
        cgmScaleToFit= 1;
        gtDoEscapes= 0;
      }
      else if (strcmp(arg, "fmbug")==0) epsFMbug= 1;
      else if (strcmp(arg, "bg0fg1")==0) bg0fg1= 1;
      else if (strcmp(arg, "esc0")==0) gtDoEscapes= 0;
      else if (strcmp(arg, "esc1")==0) gtDoEscapes= 1;
      else return HelpAndExit();

      if (fileType>=0) {
        if (nOut>=8)
          return MessageAndExit("At most 8 output files/displays allowed");
        for (j=0 ; j<nOut ; j++) if (strcmp(outNames[j], arg)==0)
          return MessageAndExit("Duplicate output filenames not allowed");
        outNames[nOut]= arg;
        outTypes[nOut]= fileType;
        nOut++;
      }

    } else if (arg[0]<='9' && arg[0]>='0') {
      if (GetPageGroup(arg)) return MessageAndExit("Try again");

    } else if (!inName) {
      inName= arg;

    } else {
      return HelpAndExit();
    }
  }

  if (inName && OpenCGM(inName)) inName= 0;

  if (amBatch) {
    if (!inName)
      return MessageAndExit("Must specify an input CGM file to use -b or -f");
    if (!nOut)
      return MessageAndExit("Must specify some output file to use -b");
    noDisplay= 1;
  }

  /* Create CGM and PostScript engines */
  for (i=0 ; i<nOut ; i++) {
    if (outTypes[i]==0) CreateCGM(i, outNames[i]);
    else if (outTypes[i]==1) CreatePS(i, outNames[i]);
  }

  /* If page list specified, do implied send command */
  if (amBatch && nPageGroups<=0) {
    mPage[0]= 1;
    nPage[0]= 32767;
    sPage[0]= 1;
    nPageGroups= 1;
  }
  if (nPageGroups>0) {
    for (i=0 ; i<8 ; i++) {
      if (!outSend[i]) GpDeactivate(outEngines[i]);
      if (outSend[i] && !GpActivate(outEngines[i])) n_active_groups++;
    }
  }

  g_initializer(&argc, argv);

  return 0;
}

int
do_startup(void)
{
  int i;
  did_startup = 1;

  if (!noDisplay) {
    int noX= 1;
    for (i=0 ; i<nOut ; i++) if (outTypes[i]==2) {
      if (!CreateX(i, outNames[i])) noX= 0;
    }
    if (noX && nOut<8) {
      if (!CreateX(nOut, 0)) {
        nOut++;
        noX= 0;
      }
    }
    noDisplay= noX;
  }
  if (noDisplay && x_only)
    return MessageAndExit("Must be an X display to use -x");
  if (x_only && !inName)
    return MessageAndExit("Must specify an input CGM file to use -x");

  if (n_active_groups)
    ReadCGM(mPage, nPage, sPage, nPageGroups);
  return 0;
}

static int need_prompt = 1;

int
on_idle(void)
{
  int flag;
  if (!did_startup && do_startup()) return 0;
  if (need_prompt && !x_only && !amBatch) p_stdout("gist> ");
  need_prompt = 0;
  flag = CatalogCGM();
  if (!flag && amBatch) p_quit();
  return flag;
}

int
on_quit(void)
{
  int i;
  for (i=0 ; i<8 ; i++) {
    if (outEngines[i]) {
      GpDeactivate(outEngines[i]);
      GpKillEngine(outEngines[i]);
    }
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

static int HelpAndExit(void)
{
  p_stderr("Usage:   gist [cgminput] [options] [page list]\n\n");
  p_stderr("   options:  -in cgminput  (open cgminput for input)\n");
  p_stderr("             -cgm cgmout   (open cgmout for output)\n");
  p_stderr("             -ps psout     (open psout for output)\n");
  p_stderr("             -display host:n.m  (connect to X server)\n");
  p_stderr("             -75  (default to 75 dpi)\n");
  p_stderr("             -100 (default to 100 dpi)\n");
  p_stderr("             -gks (default to 8x8 inches)\n");
  p_stderr("             -geometry WxH (initial window size in pixels)\n");
  p_stderr("             -nd  (not even default X server)\n");
  p_stderr("             -b   (batch mode, implies -nd)\n");
  p_stderr("             -x   (graphics window only, no keyboard)\n");
  p_stderr("             -nowarn (only one warning message printed)\n");
  p_stderr("             -f   (PostScript to stdout, implies -b)\n");
  p_stderr("             -fmbug (if EPS files are for FrameMaker)\n");
  p_stderr("             -bg0fg1 (force color 0 to bg, 1 to fg)\n");
  p_stderr("             -esc0 (skip ^_! escapes-assumed if -gks)\n");
  p_stderr("             -esc1 (handle ^_! text escapes-default)\n");
  p_stderr("   page list:  [page group] [page group] ...\n");
  p_stderr("   page group:  n     (output page n)\n");
  p_stderr("   page group:  n-m   (output pages n->m inclusive)\n");
  p_stderr("   page group:  n-m-s (output every s-th page, n-m)\n");
  p_quit();
  return 1;
}

static int MessageAndExit(char *msg)
{
  p_stderr("gist: ");
  p_stderr(msg);
  p_stderr("\ngist: try    gist -help    for usage\n");
  p_quit();
  return 1;
}

void Warning(char *general, char *particular)
{
  if (++warningCount > 5) return;
  if (no_warnings) {
    if (no_warnings==1)
      p_stderr("gist: (WARNING) rerun without -nowarn to see warnings");
    no_warnings= 2;
    return;
  }
  p_stderr("gist: (WARNING) ");
  p_stderr(general);
  p_stderr(particular);
  p_stderr("\n");
}

/* ------------------------------------------------------------------------ */

static int GetCommand(char *token)
{
  char **command= commandList;
  char *cmd;
  char *now= token;
  int quitOK= 0;

  if (!now) return -1;

  /* Quit has two synonyms, ? is a synonym for help, and f, b, and g
     are special command forms.  */
  if (strcmp(token, "exit")==0 || strcmp(token, "end")==0 ||
      strcmp(token, "quit")==0) {
    /* Quit and its synonyms cannot be abbreviated */
    now= "quit";
    quitOK= 1;
  } else if (strcmp(token, "?")==0) now= "help";
  else if (strcmp(token, "f")==0) now= "1f";
  else if (strcmp(token, "b")==0) now= "1b";
  else if (strcmp(token, "g")==0) now= "1g";
  else if (strcmp(token, "s")==0) now= "send";

  /* Check for nf, nb, ng */
  if (*now<='9' && *now>='0') {
    char *suffix;
    nPrefix= (int)strtol(now, &suffix, 10);
    cSuffix= *suffix;
    if ((cSuffix=='f' || cSuffix=='b' || cSuffix=='g') &&
        *(suffix+1)=='\0') {
      while (*command) command++;
      return command-commandList;
    }
  }

  cmd= *command;
  while (cmd && *now>*cmd) cmd= *(++command);
  while (cmd && *now==*cmd) {
    while (*cmd && *now==*cmd) { now++;  cmd++; }
    if (!*now) break;  /* token matches cmd */
    now= token;
    cmd= *(++command);
  }
  if (!cmd || *now) {
    p_stderr("gist: Unknown command: ");
    p_stderr(token);
    p_stderr("\n      Type help for help.\n");
    return -1;
  }

  if (*cmd) {
    now= token;
    cmd= *(command+1);
    if (cmd && *now++==*cmd++) {
      while (*cmd && *now==*cmd) { now++;  cmd++; }
      if (!*now) {
        p_stderr("gist: Ambiguous command: ");
        p_stderr(token);
        p_stderr("\n      Type help for help.\n");
        return -1;
      }
    }
  }

  if (strcmp(*command, "quit")==0 && !quitOK) {
    p_stderr("gist: The quit command cannot be abbreviated.\n");
    return -1;
  }

  return command-commandList;
}

static int GetPageGroup(char *token)
{
  int n, m, s;
  char *suffix, *tok;

  s= 1;
  n= m= (int)strtol(token, &suffix, 10);
  if (suffix!=token && *suffix=='-') {
    tok= suffix+1;
    n= (int)strtol(tok, &suffix, 10);
    if (suffix!=tok && *suffix=='-') {
      tok= suffix+1;
      s= (int)strtol(tok, &suffix, 10);
      if (suffix==tok) suffix= tok-1;  /* this is an error */
    }
  }

  if (*suffix) {
    p_stderr("gist: (SYNTAX) ");
    p_stderr(token);
    p_stderr(" is not a legal page number group\n");
    return 1;
  }

  if (nPageGroups>=32) {
    p_stderr("gist: (SYNTAX) too many page number groups (32 max)\n");
    return 1;
  }

  mPage[nPageGroups]= m;
  nPage[nPageGroups]= n;
  sPage[nPageGroups]= s;
  nPageGroups++;
  return 0;
}

static int CheckEOL(char *command)
{
  if (strtok(0, " \t\n")) {
    p_stderr("gist: (SYNTAX) garbage after ");
    p_stderr(command);
    p_stderr(" command\n");
    return 1;
  } else {
    return 0;
  }
}

static int Help(int help)
{
  int cmd;

  if (help) {
    p_stderr("gist: help command syntax:\n     help [command]\n");
    p_stderr("  Print help message (for given command).\n");
    return 0;
  }

  cmd= GetCommand(strtok(0, " \t\n"));

  if (cmd<0) {
    int len;
    char line[80], **command= commandList;
    p_stderr("gist: Input Syntax:\n     command [arg1] [arg2] ...\n");
    strcpy(line, "  Available commands are:  ");
    len= 27;
    for (;;) {
      for (;;) {
        if (len+strlen(*command) > 72) {
          strcpy(line+len, "\n");
          p_stderr(line);
          break;
        }
        strcpy(line+len, *command);
        len+= strlen(*command++);
        if (!*command) {
          strcpy(line+len, "\n");
          p_stderr(line);
          break;
        }
        strcpy(line+len, ", ");
        len+= 2;
      }
      if (!*command) break;
      strcpy(line, "     ");
      len= 5;
    }
    p_stderr("  The quit command has two synonyms:  exit, end\n");
    p_stderr("  Any command except quit may be abbreviated by\n");
    p_stderr("  omitting characters from the right.\n");
    p_stderr("  The help command explains specific syntax, e.g.:\n");
    p_stderr("     help cgm\n");
    p_stderr("  describes the syntax of the cgm command.\n");
    p_stderr("  Five commands can be typed to a gist X window:\n");
    p_stderr("     nf   - forward n pages and draw (default 1)\n");
    p_stderr("     nb   - backward n pages and draw (default 1)\n");
    p_stderr("     ng   - go to page n and draw (default 1)\n");
    p_stderr("     s    - send current page\n");
    p_stderr("     q    - quit\n");

  } else {
    CommandList[cmd](1);
  }

  return 0;
}

static char line[256];

void
on_stdin(char *lin)
{
  int cmd;

  line[0] = '\0';
  strncat(line, lin, 255);
  cmd= GetCommand(strtok(line, " \t\n"));

  warningCount= xPrefix= 0;
  if (cmd>=0 && CommandList[cmd](0))
    p_quit();

  need_prompt = 1;
}

static int Quit(int help)
{
  if (help) {
    p_stderr("gist: quit command syntax:\n     quit\n");
    p_stderr("  Finish and close any output files, then exit.\n");
    p_stderr("  Synonyms:  exit, end   (no abbreviations allowed)\n");
    return 0;
  }
  return 1;
}

/* ------------------------------------------------------------------------ */

static int SaveName(int device, char *name)
{
  int len;

  if (name == outNames[device]) return 0;

  len= name? (int)strlen(name) : 0;
  if (len>outLengths[device]) {
    if (outLengths[device]>0) p_free(outNames[device]);
    outNames[device]= (char *)p_malloc(len+1);
    if (!outNames[device]) {
      p_stderr("gist: (SEVERE) memory manager failed in SaveName\n");
      outLengths[device]= 0;
      return 1;
    }
    outLengths[device]= len;
  }

  if (name) strcpy(outNames[device], name);
  else if (outLengths[device]>0) outNames[device][0]= '\0';
  else outNames[device]= 0;
  return 0;
}

static int CreateCGM(int device, char *name)
{
  if (SaveName(device, name)) return 1;

  outEngines[device]=
    GpCGMEngine("CGM Viewer", cgmLandscape, 0, name);
  if (!outEngines[device]) {
    Warning(gistError, "");
    Warning("Unable to create CGM engine ", name);
    return 1;
  }

  outTypes[device]= 0;
  outDraw[device]= 0;
  outSend[device]= 1;

  return 0;
}

static int CreatePS(int device, char *name)
{
  if (SaveName(device, name)) return 1;

  outEngines[device]=
    GpPSEngine("CGM Viewer", cgmLandscape, 0, name);
  if (!outEngines[device]) {
    Warning(gistError, "");
    Warning("Unable to create PostScript engine ", name);
    return 1;
  }

  outTypes[device]= 1;
  outDraw[device]= 0;
  outSend[device]= 1;

  return 0;
}

#ifndef NO_XLIB
static int CreateX(int device, char *name)
{
  if (SaveName(device, name)) return 1;

  gist_input_hint = 1;
  outEngines[device]=
    GpBXEngine("CGM Viewer", cgmLandscape, defaultDPI, name);
  if (!outEngines[device]) {
    Warning(gistError, "");
    Warning("Unable to open X server ", name? name : "(default)");
    return 1;
  }

  outTypes[device]= 2;
  outDraw[device]= 1;
  outSend[device]= 0;

  GxInput(outEngines[device], &HandleExpose, 0, 0, &HandleOther);

  return 0;
}
#endif

/* ------------------------------------------------------------------------ */

static int OpenIn(int help)
{
  char *token;

  if (help) {
    p_stderr("gist: open command syntax:\n     open cgminput\n");
    p_stderr("  Closes the current CGM input file, then opens cgminput.\n");
    p_stderr("  Only a Gist-compliant binary CGM file is legal.\n");
    p_stderr("  The cgminput may be the first file of a family.\n");
    p_stderr("  Subsequent page numbers refer to this input file.\n");
    return 0;
  }

  token= strtok(0, " \t\n");
  if (!token) {
    p_stderr("gist: (SYNTAX) cgminput name missing in open command\n");
    return 0;
  }
  if (CheckEOL("open")) return 0;

  if (no_warnings) no_warnings= 1;  /* one warning per file family */
  OpenCGM(token);
  return 0;
}

static int FindDevice(void)
{
  int i;
  for (i=0 ; i<8 ; i++) if (!outEngines[i]) break;
  if (i>=8)
    Warning("8 devices already open for output, command ignored", "");
  return i;
}

static int MakeCGM(int help)
{
  char *token, *cgmout;
  long size= 0;
  int device;

  if (help) {
    p_stderr("gist: cgm command syntax:\n     cgm cgmout [size]\n");
    p_stderr("  Opens a CGM file cgmout for output.\n");
    p_stderr("  The size (default 1000000) is the maximum size of a\n");
    p_stderr("  single file in the output family, in bytes.\n");
    p_stderr("  Subsequent send commands will write to cgmout,\n");
    p_stderr("  unless the send to list is modified (see send).\n");
    return 0;
  }

  token= strtok(0, " \t\n");
  if (!token) {
    p_stderr("gist: (SYNTAX) cgmout name missing in cgm command\n");
    return 0;
  }
  cgmout= token;
  token= strtok(0, " \t\n");
  if (token) {
    char *suffix;
    size= strtol(token, &suffix, 0);
    if (*suffix) {
      p_stderr("gist: (SYNTAX) size unintelligble in cgm command\n");
      return 0;
    }
    if (CheckEOL("cgm")) return 0;
  }

  device= FindDevice();
  if (device>=8) return 0;

  if (!CreateCGM(device, cgmout) &&
      size>0) ((CGMEngine *)outEngines[device])->fileSize= size;

  return 0;
}

static int MakePS(int help)
{
  char *token;
  int device;

  if (help) {
    p_stderr("gist: ps command syntax:\n     ps psout\n");
    p_stderr("  Opens a PostScript file psout for output.\n");
    p_stderr("  Subsequent send commands will write to psout,\n");
    p_stderr("  unless the send to list is modified (see send).\n");
    return 0;
  }

  token= strtok(0, " \t\n");
  if (!token) {
    p_stderr("gist: (SYNTAX) psout name missing in ps command\n");
    return 0;
  }
  if (CheckEOL("ps")) return 0;

  device= FindDevice();
  if (device>=8) return 0;

  CreatePS(device, token);

  return 0;
}

static int MakeX(int help)
{
  char *token, *server;
  int dpi, device, defDPI;

  if (help) {
    p_stderr("gist: display command syntax:\n     ");
    p_stderr("display host:server.screen [dpi]\n");
    p_stderr("  Connects to the specified X server.\n");
    p_stderr("  Subsequent draw commands will write to server,\n");
    p_stderr("  unless the draw to list is modified (see draw).\n");
    p_stderr("  If specified, 40<=dpi<=200 (default 100).\n");
    return 0;
  }

  token= strtok(0, " \t\n");
  if (!token) {
    p_stderr("gist: (SYNTAX) cgmoutput name missing in cgm command\n");
    return 0;
  }
  server= token;
  token= strtok(0, " \t\n");
  if (token) {
    char *suffix;
    dpi= (int)strtol(token, &suffix, 0);
    if (*suffix) {
      p_stderr("gist: (SYNTAX) dpi unintelligble in display command\n");
      return 0;
    }
    if (dpi<40 && dpi>200) {
      p_stderr(
        "gist: (SYNTAX) dpi not between 40 and 200 in display command\n");
      return 0;
    }
    if (CheckEOL("display")) return 0;
  } else {
    dpi= 100;
  }

  device= FindDevice();
  if (device>=8) return 0;

  defDPI= defaultDPI;
  defaultDPI= dpi;
  CreateX(device, server);
  defaultDPI= defDPI;

  return 0;
}

static int EPS(int help)
{
  char *token;
  int device;

  if (help) {
    p_stderr("gist: eps command syntax:\n     eps epsout\n");
    p_stderr("  Open an Encapsulated PostScript file epsout, write\n");
    p_stderr("  the current page to it, then close epsout.\n");
    p_stderr("  (Note that an EPS file can have only a single page.)\n");
    return 0;
  }

  token= strtok(0, " \t\n");
  if (!token) {
    p_stderr("gist: (SYNTAX) epsout name missing in eps command\n");
    return 0;
  }
  if (CheckEOL("eps")) return 0;

  device= FindDevice();
  if (device>=8) return 0;

  device= FindDevice();
  if (device>=8) return 0;

  outEngines[device]=
    GpPSEngine("CGM Viewer", cgmLandscape, 0, "_tmp.eps");
  if (!outEngines[device]) {
    Warning(gistError, "");
    Warning("Unable to create PostScript engine ", token);
    return 0;
  }

  GpPreempt(outEngines[device]);

  nPage[0]= mPage[0]= CGMRelative(0);
  sPage[0]= 1;
  nPageGroups= 1;

  /* First read does PS part, second computes EPS preview */
  if (!ReadCGM(mPage, nPage, sPage, nPageGroups)) {
    GpPreempt(0);
    outEngines[device]= EPSPreview(outEngines[device], token);
    if (outEngines[device]) {
      GpPreempt(outEngines[device]);
      ReadCGM(mPage, nPage, sPage, nPageGroups);
    } else {
      Warning("memory manager failed creating EPS engine ", token);
      return 0;
    }
  }

  if (outEngines[device]) {
    GpPreempt(0);

    GpKillEngine(outEngines[device]);
    outEngines[device]= 0;
  }

  return 0;
}

static int FreeAny(int help)
{
  int i;

  if (help) {
    p_stderr("gist: free command syntax:\n     free [device# ...]\n");
    p_stderr("  Finish and close the device#(s).  If none given,\n");
    p_stderr("  frees all send devices,\n");
    p_stderr("  (Use the info command to describe current device numbers.)\n");
    return 0;
  }

  if (GetDeviceList(1, "free")) return 0;

  for (i=0 ; i<8 ; i++) {
    if (outMark[i]) {
      GpKillEngine(outEngines[i]);
      outEngines[i]= 0;
    }
  }

  return 0;
}

static char *yorn[2]= { "No  ", "Yes " };
static char *tname[3]= { "CGM", "PS ", "X  " };

extern void CGMinfo(void);  /* cgmin.c */

static int Info(int help)
{
  int i;

  if (help) {
    p_stderr("gist: info command syntax:\n     info\n");
    p_stderr("  Print descriptions of all current output files.\n");
    return 0;
  }

  p_stdout("\n");
  CGMinfo();
  p_stdout("Number Draw Send Type   Name\n");
  for (i=0 ; i<8 ; i++) {
    if (outEngines[i]) {
      char msg[80];
      sprintf(msg, "%3d    %s %s %s  ", i, yorn[outDraw[i]], yorn[outSend[i]],
             tname[outTypes[i]]);
      p_stdout(msg);
      p_stdout(outNames[i]? outNames[i] : "");
      p_stdout("\n");
    }
  }
  p_stdout("\n");

  return 0;
}

/* ------------------------------------------------------------------------ */

static int Draw(int help)
{
  int i, n;
  char *token;

  if (help) {
    p_stderr("gist: draw command syntax:\n     draw [page list]\n");
    p_stderr(
    "  Copy the page(s) (default current page) from the current CGM input\n");
    p_stderr("  to all display output devices.\n");
    p_stderr("  By default, these are all X windows.\n");
    p_stderr("  Use alternate syntax:\n     draw to [device#1 ...]\n");
    p_stderr("  to specify a particular list of devices to be used\n");
    p_stderr("  by the draw command.  Without any device numbers,\n");
    p_stderr("  draw to restores the default list of devices.\n");
    p_stderr("  (Use the info command to describe current device numbers.)\n");
    p_stderr("  Page list syntax:   group1 [group2 ...]\n");
    p_stderr("  Page group syntax:   n   - just page n\n");
    p_stderr("                     m-n   - pages n thru m\n");
    p_stderr("                   m-n-s   - pages n thru m, step s\n");
    return 0;
  }

  token= strtok(0, " \t\n");
  if (token && strcmp(token, "to")==0) {
    DrawSend(1, token);
    return 0;
  }

  n= 0;
  for (i=0 ; i<8 ; i++) {
    if (!outDraw[i]) GpDeactivate(outEngines[i]);
    if (outDraw[i] && !GpActivate(outEngines[i])) n++;
  }

  if (!n && (i= FindDevice())<8) {
    if (!CreateX(i, 0) && !GpActivate(outEngines[i])) n++;
  }

  if (n) DrawSend(0, token);
  else Warning("no devices active for draw command", "");
  return 0;
}

static int Send(int help)
{
  int i, n;
  char *token;

  if (help) {
    p_stderr("gist: send command syntax:\n     send [page list]\n");
    p_stderr(
    "  Copy the page(s) (default current page) from the current CGM input\n");
    p_stderr("  to all display output devices.\n");
    p_stderr("  By default, these are all X windows.\n");
    p_stderr("  Use alternate syntax:\n     send to [device#1] ...\n");
    p_stderr("  to specify a particular list of devices to be used\n");
    p_stderr("  by the send command.  Without any device numbers,\n");
    p_stderr("  send to restores the default list of devices.\n");
    p_stderr("  (Use the info command to describe current device numbers.)\n");
    p_stderr("  Page list syntax:   group1 [group2 ...]\n");
    p_stderr("  Page group syntax:   n   - just page n\n");
    p_stderr("                     m-n   - pages n thru m\n");
    p_stderr("                   m-n-s   - pages n thru m, step s\n");
    return 0;
  }

  token= strtok(0, " \t\n");
  if (token && strcmp(token, "to")==0) {
    DrawSend(1, token);
    return 0;
  }

  n= 0;
  for (i=0 ; i<8 ; i++) {
    if (!outSend[i]) GpDeactivate(outEngines[i]);
    if (outSend[i] && !GpActivate(outEngines[i])) n++;
  }

  if (n) DrawSend(1, token);
  else Warning("no devices active for send command", "");
  return 0;
}

static int GetDeviceList(int ds, char *command)
{
  int device;
  char *token= strtok(0, " \t\n");

  if (token) {
    char *suffix;
    for (device=0 ; device<8 ; device++) outMark[device]= 0;
    do {
      device= (int)strtol(token, &suffix, 10);
      if (*suffix || device<0 || device>7) {
        p_stderr("gist: (SYNTAX) ");
        p_stderr(token);
        p_stderr(" not a legal device# in ");
        p_stderr(command);
        p_stderr(" command\n");
        return 1;
      }
      if (!outEngines[device])
        Warning("ignoring non-existent device# ", token);
      else
        outMark[device]= 1;
    } while ((token= strtok(0, " \t\n")));
    if (ds==0) {
      for (device=0 ; device<8 ; device++) if (outMark[device]) break;
      if (device>=8) {
        Warning(command, " command with no legal devices, no action taken");
        return 1;
      }
    }

  } else if (ds==0) {
    for (device=0 ; device<8 ; device++) outMark[device]= 0;
  } else {
    for (device=0 ; device<8 ; device++) outMark[device]= outSend[device];
  }

  return 0;
}

static char *dsName[]= { "draw to", "send to" };

static void DrawSend(int ds, char *token)
{
  nPageGroups= 0;

  if (token && strcmp(token, "to")==0) {
    /* draw to  or send to  merely resets outDraw or outSend list */
    int i, n= 0;
    int *outDS= ds? outSend : outDraw;
    if (GetDeviceList(0, dsName[ds])) return;
    for (i=0 ; i<8 ; i++) if ((outDS[i]= outMark[i])) n++;
    if (!n) for (i=0 ; i<8 ; i++)
      outDS[i]= outEngines[i]? (ds? outTypes[i]<2 : outTypes[i]==2) : 0;
    return;

  } else if (token) {
    do {
      if (GetPageGroup(token)) return;
    } while ((token= strtok(0, " \t\n")));
  } else {
    nPage[0]= mPage[0]= CGMRelative(0);
    sPage[0]= 1;
    nPageGroups= 1;
  }

  ReadCGM(mPage, nPage, sPage, nPageGroups);
}

static int Special(int help)
{
  if (help) {
    char msg[80];
    sprintf(msg, "gist: n%c command syntax:\n     n%c\n",
            cSuffix, cSuffix);
    p_stderr(msg);
    if (cSuffix=='f')
      p_stderr("  Forward n (default 1) pages, then draw\n");
    else if (cSuffix=='b')
      p_stderr("  Backward n (default 1) pages, then draw\n");
    else
      p_stderr("  Go to page n (default 1), then draw\n");
    return 0;
  }

  if (CheckEOL("nf, nb, or ng")) return 0;

  DoSpecial(nPrefix, cSuffix);
  return 0;
}

static void DoSpecial(int nPrefix, int cSuffix)
{
  int i, n;

  n= 0;
  for (i=0 ; i<8 ; i++) {
    if (!outDraw[i]) GpDeactivate(outEngines[i]);
    if (outDraw[i] && !GpActivate(outEngines[i])) n++;
  }

  if (cSuffix=='f') mPage[0]= CGMRelative(nPrefix);
  else if (cSuffix=='b') mPage[0]= CGMRelative(-nPrefix);
  else mPage[0]= nPrefix;
  nPage[0]= mPage[0];
  sPage[0]= 1;
  nPageGroups= 1;

  if (n) ReadCGM(mPage, nPage, sPage, nPageGroups);
  else Warning("no devices active for nf, nb, or ng command", "");
}

/* ------------------------------------------------------------------------ */

/* ARGSUSED */
static void HandleExpose(Engine *engine, Drauing *drauing, int xy[])
{
  XEngine *xEngine= GisXEngine(engine);

  if (!xEngine) return;

  /* Redraw current picture on this engine only */

  GpPreempt(engine);

  nPage[0]= mPage[0]= CGMRelative(0);
  sPage[0]= 1;
  nPageGroups= 1;

  ReadCGM(mPage, nPage, sPage, nPageGroups);

  GpPreempt(0);
  GxExpose(engine, drauing, xy);
}

static int cSuffices[]= { 'f', 'b', 'g' };

/* ARGSUSED */
static void HandleOther(Engine *engine, int k, int md)
{
  XEngine *xEngine= GisXEngine(engine);
  int go, bad;

  if (!xEngine) return;

  go= bad= 0;

  if (k == '0') xPrefix= 10*xPrefix;
  else if (k == '1') xPrefix= 10*xPrefix+1;
  else if (k == '2') xPrefix= 10*xPrefix+2;
  else if (k == '3') xPrefix= 10*xPrefix+3;
  else if (k == '4') xPrefix= 10*xPrefix+4;
  else if (k == '5') xPrefix= 10*xPrefix+5;
  else if (k == '6') xPrefix= 10*xPrefix+6;
  else if (k == '7') xPrefix= 10*xPrefix+7;
  else if (k == '8') xPrefix= 10*xPrefix+8;
  else if (k == '9') xPrefix= 10*xPrefix+9;
  else if (k=='f' || k=='F' || (k=='+' && (md&P_KEYPAD))) go= 1;
  else if (k=='b' || k=='B' || (k=='-' && (md&P_KEYPAD))) go= 2;
  else if (k=='g' || k=='G' || k=='\r') go= 3;
  else if (k=='s' || k=='S' || (k=='=' && (md&P_KEYPAD))) go= 4;
  else if (k=='q' || k=='Q') go= 5;
  else bad= 1;

  if ((go==4||go==5) && xPrefix!=0) bad= 1;
  if (go && !bad) {
    if (go<4) {
      if (xPrefix==0) xPrefix= 1;
      DoSpecial(xPrefix, cSuffices[go-1]);
    } else if (go==4) {
      int i, n= 0;
      for (i=0 ; i<8 ; i++) {
        if (!outSend[i]) GpDeactivate(outEngines[i]);
        if (outSend[i] && !GpActivate(outEngines[i])) n++;
      }

      nPage[0]= mPage[0]= CGMRelative(0);
      sPage[0]= 1;
      nPageGroups= 1;

      if (n) ReadCGM(mPage, nPage, sPage, nPageGroups);
      else Warning("no devices active for send command", "");
    } else if (go==5) {
      p_quit();
    }
    xPrefix= 0;
    warningCount= 0;
  } else if (bad) {
    p_feep(xEngine->win);
    xPrefix= 0;
  }
}

/* ------------------------------------------------------------------------ */
