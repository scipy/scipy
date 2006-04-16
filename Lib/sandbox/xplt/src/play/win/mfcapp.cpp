/*
 * mfcapp.cpp -- $Id$
 * MFC implementation of play, boss and worker thread classes
 *
 * Copyright (c) 2000.  See accompanying LEGAL file for details.
 */

#include "mfcapp.h"
#include "mfcterm.h"
#include "mfcres.h"

#include "playw.h"
#include "pstdlib.h"

static int w_result = 0;

static CWinThread *the_worker = 0;
static CMultiDocTemplate *mfc_template = 0;

static void on_view(mfc_term_view *view);
static void on_update_view(mfc_term_view *view, CCmdUI *ui);
static int mfc_showing = 0;
static void mfc_show();
static char *mfc_cpy(HANDLE heap, const char *text, long len);
static int mfc_run_hack();
static int (*w_on_launch)(int, char **)= 0;

static void mfc_quit(void);
static int mfc_quitting = 0;
static HWND mfc_parent(int width, int height, char *title, int hints);

static void
mfc_quit(void)
{
  if (!mfc_quitting)
    the_boss.m_pMainWnd->SendMessage(WM_CLOSE, 0, 0);
  mfc_quitting = 1;
}

/*------------------------------------------------------------------------*/
/* boss thread (main application) */

class mfc_worker : public CWinThread {
  DECLARE_DYNCREATE(mfc_worker)
public:
  mfc_worker();
  virtual BOOL InitInstance();
  virtual int Run();
  virtual BOOL PreTranslateMessage(MSG* pmsg);
  virtual BOOL OnIdle(LONG count);
  virtual int ExitInstance();
  int run_hack();
};

BEGIN_MESSAGE_MAP(mfc_boss, CWinApp)
  ON_COMMAND(ID_APP_ABOUT, on_about)
  ON_COMMAND(ID_SIGINT, on_sigint)
  ON_COMMAND(ID_VIEW_TERM, on_view_term)
  ON_COMMAND(ID_VIEW_HIST, on_view_hist)
  ON_UPDATE_COMMAND_UI(ID_VIEW_TERM, on_update_view_term)
  ON_UPDATE_COMMAND_UI(ID_VIEW_HIST, on_update_view_hist)
  ON_COMMAND(ID_FILE_NEW, CWinApp::OnFileNew)
  ON_COMMAND(ID_FILE_OPEN, CWinApp::OnFileOpen)
  ON_COMMAND(ID_FILE_PRINT_SETUP, CWinApp::OnFilePrintSetup)
END_MESSAGE_MAP()

static int w_argc = 0;
static char **w_argv = 0;

mfc_boss::mfc_boss(int (*on_launch)(int, char **))
{
  w_on_launch = on_launch;
}

class mfc_frame_wnd : public CMDIFrameWnd
{
  DECLARE_DYNCREATE(mfc_frame_wnd)
public:
  mfc_frame_wnd();
  virtual LRESULT WindowProc(UINT message, WPARAM wParam, LPARAM lParam);
};

IMPLEMENT_DYNCREATE(mfc_frame_wnd, CMDIFrameWnd)

mfc_frame_wnd::mfc_frame_wnd()
{
}
static int oops = 0;
LRESULT
mfc_frame_wnd::WindowProc(UINT message, WPARAM wParam, LPARAM lParam) 
{
  LRESULT result = 0;
  if (message==ID_CALL_FUNC) {
    void (*func)(void *arg)= (void (*)(void *))lParam;
    void *arg = (void *)wParam;
    func(arg);
  } else {
    if (message==WM_QUERYNEWPALETTE && w_screen.active)
      return ::SendMessage(w_screen.active->w, message, wParam, lParam);
    result = CMDIFrameWnd::WindowProc(message, wParam, lParam);
  }
  return result;
}

BOOL
mfc_boss::InitInstance()
{
  if (!AfxOleInit()) {  /* all this for richedit... */
    AfxMessageBox(IDP_OLE_INIT_FAILED);
    return 0;
  }

#ifdef _AFXDLL
  Enable3dControls();
#else
  Enable3dControlsStatic();
#endif

  SetRegistryKey(m_pszAppName); // AFX_IDS_APP_TITLE, considered IDR_MAINFRAME
  LoadStdProfileSettings(10); // Load standard INI file options (including MRU)

  mfc_template =
    new CMultiDocTemplate(IDR_EDITFRAME,
                          RUNTIME_CLASS(mfc_edit_doc),
                          RUNTIME_CLASS(mfc_edit_child),
                          RUNTIME_CLASS(mfc_edit_view));
  AddDocTemplate(mfc_template);

  /* create the main window, but don't show it until
   * the terminal window or a graphics child window is created
   * -- MFC Run method kills process if no main window */
  CMDIFrameWnd* fw = new mfc_frame_wnd;
  if (!fw->LoadFrame(IDR_MAINFRAME)) return 0;
  m_pMainWnd = fw;
  m_pMainWnd->DragAcceptFiles();

  EnableShellOpen();
  RegisterShellFileTypes(0);  /* do not register for printing */

  /* some command line stuff needs to be merged with w_init... */
  CCommandLineInfo cmd_info;
  ParseCommandLine(cmd_info);
  /* prevent opening an empty untitled document at startup */
  if (cmd_info.m_nShellCommand == CCommandLineInfo::FileNew)
    cmd_info.m_nShellCommand = CCommandLineInfo::FileNothing;
  if (!ProcessShellCommand(cmd_info)) return 0;

  the_worker=
    AfxBeginThread(RUNTIME_CLASS(mfc_worker), THREAD_PRIORITY_BELOW_NORMAL);
  return 1;
}

int
mfc_boss::ExitInstance()
{
  the_worker->PostThreadMessage(WM_QUIT,0,0);
  Sleep(200);
  if (!mfc_quitting) {
    mfc_quitting = 1;
    w_sigint(0);
    Sleep(200);
  }
  CWinApp::ExitInstance();
  return w_result;
}

void
mfc_boss::on_about()
{
  CDialog about(IDD_ABOUTBOX);
  about.DoModal();
}

void
mfc_boss::on_sigint()
{
  mfc_reset_stdin();
  if (!w_sigint(1)) mfc_quit();
}

void
mfc_boss::on_view_term()
{
  on_view(term_view);
}

void
mfc_boss::on_view_hist()
{
  on_view(hist_view);
}

void
mfc_boss::on_update_view_term(CCmdUI *ui)
{
  on_update_view(term_view, ui);
}

void
mfc_boss::on_update_view_hist(CCmdUI *ui)
{
  on_update_view(hist_view, ui);
}

static void
on_view(mfc_term_view *view)
{
  if (view && IsWindow(view->m_hWnd)) {
    CMDIFrameWnd *fw = (CMDIFrameWnd *)(the_boss.m_pMainWnd);
    if (view->is_visible) {
      CMDIChildWnd *cw = (CMDIChildWnd *)fw->GetActiveFrame();
      if (view == cw->GetActiveView()) fw->MDINext();
      view->GetParent()->ShowWindow(SW_HIDE);
    } else {
      view->GetParent()->ShowWindow(SW_SHOW);
          fw->MDIActivate(view->GetParent());
    }
    view->is_visible = !view->is_visible;
    SendMessage(fw->m_hWndMDIClient, WM_MDIREFRESHMENU, 0, 0);
  }
}

static void
on_update_view(mfc_term_view *view, CCmdUI *ui)
{
  if (view && IsWindow(view->m_hWnd)) {
    ui->Enable(1);
    ui->SetCheck(view->is_visible);
  } else {
    ui->Enable(0);
  }
}

/*------------------------------------------------------------------------*/
/* worker thread */

IMPLEMENT_DYNCREATE(mfc_worker, CWinThread)

mfc_worker::mfc_worker()
{
}

BOOL
mfc_worker::InitInstance()
{
  w_initialize(the_boss.m_hInstance, the_boss.m_pMainWnd->m_hWnd,
               mfc_quit, mfc_stdinit, mfc_parent);

  /* crack command line into argc, argv */
  HANDLE heap = GetProcessHeap();
  LPSTR cmd_line = the_boss.m_lpCmdLine;
  char module_name[1028];
  DWORD len = GetModuleFileName((HMODULE)0, module_name, 1024);
  module_name[len] = '\0';
  w_argc = 0;
  w_argv = (char**)HeapAlloc(heap, HEAP_GENERATE_EXCEPTIONS, sizeof(char *)*9);
  w_argv[w_argc++] = mfc_cpy(heap, w_unixpath(module_name), len);
  if (cmd_line) {
    char *c = cmd_line;
    char delim;
    for (;;) {
      while (c[0]==' ' || c[0]=='\t' || c[0]=='\r' || c[0]=='\n') c++;
      delim = c[0];
      if (!delim) break;
      cmd_line = c;
      if (delim=='"' || delim=='\'') {
        cmd_line = ++c;
        while (c[0] && c[0]!=delim) c++;
      } else {
        while (c[0] && c[0]!=' ' && c[0]!='\t' &&
               c[0]!='\r' && c[0]!='\n') c++;
        delim = 'x';
      }
      if (w_argc>1 || cmd_line[0]!='-'||cmd_line[1]!='n'||cmd_line[2]!='o'||
          cmd_line[3]!='m'||cmd_line[4]!='d'||cmd_line[5]!='i') {
        if (!(w_argc&7))
          w_argv = (char **)HeapReAlloc(heap, HEAP_GENERATE_EXCEPTIONS,
                                       w_argv, sizeof(char *)*(2*w_argc+1));
        w_argv[w_argc++] = mfc_cpy(heap, cmd_line, c - cmd_line);
      } else {
        w_no_mdi = 1;
      }
      if (c[0] == delim) c++;
      cmd_line = c;
    }
  }
  w_argv[w_argc] = 0;

  return 1;
}

static char *
mfc_cpy(HANDLE heap, const char *text, long len)
{
  if (!len) while (text[len]) len++;
  char *buf = (char *)HeapAlloc(heap, HEAP_GENERATE_EXCEPTIONS, len+1);
  char *buffer = buf;
  while (len--) *buf++= *text++;
  buf[0] = '\0';
  return buffer;
}

int
mfc_worker::Run()
{
  return w_protect(mfc_run_hack);
}

int
mfc_worker::run_hack()
{
  if (w_on_launch) {
    int (*on_launch)(int, char **)= w_on_launch;
    w_on_launch = 0;
    p_mminit();
    w_pollinit();
    int result = on_launch(w_argc, w_argv);
    if (result) return result;
  }
  return CWinThread::Run();
}

static int
mfc_run_hack()
{
  return ((mfc_worker *)the_worker)->run_hack();
}

BOOL
mfc_worker::PreTranslateMessage(MSG* pmsg)
{
  return w_app_msg(pmsg) || CWinThread::PreTranslateMessage(pmsg);
}

BOOL
mfc_worker::OnIdle(LONG count)
{
  return CWinThread::OnIdle(count) || w_work_idle();
}

int
mfc_worker::ExitInstance()
{
  w_result = w_on_quit();
  CWinThread::ExitInstance();
  AfxEndThread(w_result, 1);  // simply returning doesn't delete the_worker
  return w_result;            // not reached
}

/*------------------------------------------------------------------------*/

static void
mfc_show()
{
  if (!(mfc_showing&1)) {
    CMDIFrameWnd *fw = (CMDIFrameWnd *)(the_boss.m_pMainWnd);
    fw->ShowWindow(the_boss.m_nCmdShow);
    fw->UpdateWindow();
    mfc_showing |= 1;
  }
}

/* ARGSUSED */
void
mfc_term_init(void *arg)
{
  mfc_show();

  if (!(mfc_showing&2)) {
    mfc_edit_doc *doc = new mfc_edit_doc(mfc_template, 0);
    term_view = (mfc_term_view *)doc->GetView();
    term_view->is_visible = 1;
    mfc_template->InitialUpdateFrame((CFrameWnd*)term_view->GetParent(),
                                     doc, 1);
    doc = new mfc_edit_doc(mfc_template, 1);
    hist_view = (mfc_term_view *)doc->GetView();
    hist_view->is_visible = 0;
    mfc_template->InitialUpdateFrame((CFrameWnd*)hist_view->GetParent(),
                                     doc, 0);

    mfc_showing |= 2;
  }
}

/*------------------------------------------------------------------------*/

class mfc_child : public CMDIChildWnd
{
  DECLARE_DYNCREATE(mfc_child)
public:
  mfc_child();
  int width, height, hints;
  HWND child;

protected:
  int precreated;
  virtual BOOL PreCreateWindow( CREATESTRUCT& cs );
  virtual ~mfc_child();

  afx_msg void OnParentNotify(UINT msg, LPARAM lp);
  afx_msg void OnDestroy();
  afx_msg void OnSetFocus(CWnd *w);
  afx_msg void OnSize(UINT type, int cx, int cy);
  afx_msg void OnMDIActivate(BOOL activate, CWnd* aw, CWnd* dw);

  DECLARE_MESSAGE_MAP()
};

struct w_parent_args {
  char *title;
  int width, height, hints;
  HWND handle;
};
static void mfc_window(struct w_parent_args *args);

static HWND
mfc_parent(int width, int height, char *title, int hints)
{
  if (mfc_showing) {
    struct w_parent_args args;
    args.title = title;
    args.width = width;
    args.height = height;
    args.hints = hints;
    args.handle = 0;
    the_boss.m_pMainWnd->SendMessage(ID_CALL_FUNC, (WPARAM)&args,
                                     (LPARAM)mfc_window);
    return args.handle;
  } else {
    return 0;
  }
}

static void
mfc_window(struct w_parent_args *args)
{
  mfc_child *gw = new mfc_child;
  gw->width = args->width;
  gw->height = args->height;
  gw->hints = args->hints;
  gw->child = 0;
  gw->Create(0, args->title);
  args->handle = gw->m_hWnd;
}

BEGIN_MESSAGE_MAP(mfc_child, CMDIChildWnd)
  ON_WM_PARENTNOTIFY()
  ON_WM_DESTROY()
  ON_WM_SETFOCUS()
  ON_WM_SIZE()
  ON_WM_MDIACTIVATE()
END_MESSAGE_MAP()

IMPLEMENT_DYNCREATE(mfc_child, CMDIChildWnd)

mfc_child::mfc_child()
{
  precreated = 0;
}

mfc_child::~mfc_child()
{
}

BOOL
mfc_child::PreCreateWindow(CREATESTRUCT& cs)
{
  if (!precreated) {
    RECT rect;
    rect.left = 0;
    rect.top = 0;
    rect.right = cs.cx = width;
    rect.bottom = cs.cy = height;
    if (AdjustWindowRect(&rect, cs.style, 0)) {
      cs.cx = rect.right - rect.left + 4;  // is 4 really some variable?
      cs.cy = rect.bottom - rect.top + 4;
    }
    if (hints & P_NORESIZE)
      cs.style &= ~(WS_THICKFRAME | WS_MAXIMIZEBOX);
    precreated = 1;
  }
  return CMDIChildWnd::PreCreateWindow(cs);
}

void
mfc_child::OnDestroy() 
{
  CMDIChildWnd::OnDestroy();
  CMDIFrameWnd *fw = (CMDIFrameWnd *)the_boss.m_pMainWnd;
  ::SendMessage(fw->m_hWndMDIClient, WM_MDIREFRESHMENU, 0, 0);
}

void
mfc_child::OnSetFocus(CWnd *w) 
{
  /* first call to OnSetFocus occurs before child is created */
  if (child) {
    ::SetFocus(child);
    CMDIFrameWnd *fw = (CMDIFrameWnd *)the_boss.m_pMainWnd;
    ::SendMessage(fw->m_hWndMDIClient, WM_MDIREFRESHMENU, 0, 0);
  } else {
    CWnd::OnSetFocus(w);
  }
  /* apparently AttachThreadInput in mfc_worker::InitInstance unneeded? */
}

void
mfc_child::OnParentNotify(UINT msg, LPARAM lp)
{
  if (msg==WM_CREATE) {
    child = (HWND)lp;
    ::SetFocus(child);  /* mfc_child (this) already has focus */
  }
  CWnd::OnParentNotify(msg, lp);
}

void
mfc_child::OnSize(UINT type, int cx, int cy)
{
  RECT r;
  GetClientRect(&r);
  ::MoveWindow(child, 0,0, r.right,r.bottom, 1);
  CMDIChildWnd::OnSize(type, cx, cy);
}

void
mfc_child::OnMDIActivate(BOOL activate, CWnd* aw, CWnd* dw)
{
  if (child) {
    p_win *pw = child? (p_win *)GetWindowLong(child, GWL_USERDATA) : 0;
    if (pw && activate && pw->palette) pw->s->active = pw;
  }
  CMDIChildWnd::OnMDIActivate(activate, aw, dw);
}

/*------------------------------------------------------------------------*/
