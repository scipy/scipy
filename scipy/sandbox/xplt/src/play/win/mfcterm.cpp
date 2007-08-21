/*
 * mfcterm.cpp
 * richedit class for play MDI development environment
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "mfcapp.h"
#include "mfcterm.h"
#include "mfcres.h"
#include "playw.h"

#include <afxole.h>

#include <malloc.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

static HANDLE w_sending = 0;
static void mfc_sender(void *context);
static void mfc_stdout(char *output_line, long len);
static void mfc_stderr(char *output_line, long len);
static void mfc_bstdout(char *output_line);
static void mfc_bstderr(char *output_line);
static int mfc_deliver(const char *buf, long len);

static char *w_deprompt(char *text);

IMPLEMENT_DYNCREATE(mfc_edit_view, CRichEditView)

BEGIN_MESSAGE_MAP(mfc_edit_view, CRichEditView)
  ON_COMMAND(ID_GOTO_LINE, on_goto_line)
  ON_COMMAND(ID_GOTO_OUT, on_goto_out)
  ON_COMMAND(ID_SELECT_OUT, on_select_out)
  ON_COMMAND(ID_OPEN_LINE, on_open_line)
  ON_COMMAND(ID_CHOOSE_FILE, on_choose_file)
  ON_COMMAND(ID_EDIT_REDO, on_redo)
  ON_UPDATE_COMMAND_UI(ID_GOTO_OUT, on_update_term)
  ON_UPDATE_COMMAND_UI(ID_SELECT_OUT, on_update_term)
  ON_UPDATE_COMMAND_UI(ID_OPEN_LINE, on_update_term)
  ON_UPDATE_COMMAND_UI(ID_CHOOSE_FILE, on_update_term)
  ON_UPDATE_COMMAND_UI(ID_EDIT_REDO, on_update_redo)
  ON_WM_CREATE()

  ON_COMMAND(ID_FILE_PRINT, CRichEditView::OnFilePrint)
  ON_COMMAND(ID_FILE_PRINT_DIRECT, CRichEditView::OnFilePrint)
  ON_COMMAND(ID_FILE_PRINT_PREVIEW, CRichEditView::OnFilePrintPreview)
END_MESSAGE_MAP()

mfc_edit_view::mfc_edit_view()
{
  m_nWordWrap = WrapNone;
  is_term = riched2 = 0;
}

mfc_edit_view::~mfc_edit_view()
{
}

BOOL
mfc_edit_view::PreCreateWindow(CREATESTRUCT& cs)
{
  BOOL result = CRichEditView::PreCreateWindow(cs);

#if (_RICHEDIT_VER < 0x200) && !defined(NO_RICHEDIT2)
// The following lines will turn on richedit 2.0 or 3.0 controls.
// However, the compatibility with the rest of the class is botched so that
//   many of the MFC functions only work with 1.0 semantics (_RICHEDIT_VER is
//   set in afxwin.h to force 1.0 in riched.h).
// I prefer not to risk this with the terminal and history windows, but
//   the lure of multilevel undo is irresistable for editor windows.
  static int re2_initialized = 0;
  if (result) {
    if (!re2_initialized) {
      if (!LoadLibraryA("riched20.dll"))  // riched32.dll is 1.0 version
        re2_initialized = 2;
      else
        re2_initialized = 1;
    }
    if (re2_initialized==1) {
      cs.lpszClass = RICHEDIT_CLASSA;  // RICHEDIT_CLASS set to 1.0 version
                                      // unless _RICHEDIT_VER >= 0x200
      riched2 = 1;
    }
  }
#endif

  return result;
}

void
mfc_edit_view::OnInitialUpdate()
{
  CRichEditView::OnInitialUpdate();
  // set the printing margins (720 twips = 1/2 inch).
  SetMargins(CRect(720, 720, 720, 720));
}

BOOL
mfc_edit_view::OnPreparePrinting(CPrintInfo* pInfo)
{
  // default preparation
  return DoPreparePrinting(pInfo);
}

void
mfc_edit_view::OnDestroy()
{
  CRichEditView::OnDestroy();

  // Deactivate the item on destruction; this is important
  // when a splitter view is being used.
  COleClientItem* item = GetDocument()->GetInPlaceActiveItem(this);
  if (item && item->GetActiveView()==this)
    item->Deactivate();
}

static CFont fixed_font;
static int font_init = 0;

int
mfc_edit_view::OnCreate(LPCREATESTRUCT lpcs)
{
  if (CCtrlView::OnCreate(lpcs) != 0)
    return -1;
  GetRichEditCtrl().LimitText(lMaxSize);
  GetRichEditCtrl().SetEventMask(ENM_SELCHANGE | ENM_CHANGE | ENM_SCROLL);
  VERIFY(GetRichEditCtrl().SetOLECallback(&m_xRichEditOleCallback));
  m_lpRichEditOle = GetRichEditCtrl().GetIRichEditOle();
  DragAcceptFiles(0);
  GetRichEditCtrl().SetOptions(ECOOP_OR, ECO_AUTOWORDSELECTION);

  if (!font_init) {
    HFONT f = (HFONT)GetStockObject(ANSI_FIXED_FONT);
        LOGFONT lf;
        if (::GetObject((HGDIOBJ)f, sizeof(LOGFONT), &lf)) {
      fixed_font.CreateFontIndirect(&lf);
    } else {
      return 1;
    }
    font_init = 1;
  }
  GetRichEditCtrl().SetFont(&fixed_font,0);

  WrapChanged();
  ASSERT(m_lpRichEditOle != NULL);
  return 0;
}

// without this, idiot will paste files and other objects into window
HRESULT
mfc_edit_view::QueryAcceptData(LPDATAOBJECT lpdataobj,
                               CLIPFORMAT* lpcfFormat, DWORD reco,
                               BOOL fReally, HGLOBAL hMetaPict)
{
  if (*lpcfFormat == CF_TEXT) return S_OK;
  COleDataObject dataobj;
  dataobj.Attach(lpdataobj, FALSE);
  if (*lpcfFormat==0 && dataobj.IsDataAvailable(CF_TEXT)) {
    *lpcfFormat = CF_TEXT;
    return S_OK;
  }
  return S_FALSE;
}

class mfc_goto : public CDialog
{
public:
  mfc_goto(mfc_edit_view *v);
  virtual void OnOK();
  mfc_edit_view *view;
  DECLARE_MESSAGE_MAP()
};

BEGIN_MESSAGE_MAP(mfc_goto, CDialog)
END_MESSAGE_MAP()

mfc_goto::mfc_goto(mfc_edit_view *v) : CDialog(IDD_GOTO_LINE)
{
  view = v;
}

void
mfc_goto::OnOK()
{
  CEdit *ce = (CEdit *)GetDlgItem(IDC_GOTO_LINE);
  char line[16];
  int n = ce->GetLine(0, line, 12);
  if (n>0) {
    line[n] = '\0';
    n = atoi(line)-1;
    if (n < 0) n = 0;
    n = view->GetRichEditCtrl().LineIndex(n);
    view->GetRichEditCtrl().SetSel(n,n);
  }
  CDialog::OnOK();
}

void
mfc_edit_view::on_goto_line()
{
  mfc_goto goto_dlg(this);
  goto_dlg.DoModal();
}

void
mfc_edit_view::on_redo()
{
  ASSERT(::IsWindow(m_hWnd));
  ::SendMessage(m_hWnd, EM_REDO, 0, 0);
  m_bSyncCharFormat = m_bSyncParaFormat = 1;
}

void
mfc_edit_view::on_update_term(CCmdUI *ui)
{
  ui->Enable(is_term);
}

void
mfc_edit_view::on_update_redo(CCmdUI *ui)
{
  if (riched2) {
    ASSERT(::IsWindow(m_hWnd));
    ui->Enable((BOOL)::SendMessage(m_hWnd, EM_CANREDO, 0, 0));
  } else {
    ui->Enable(0);
  }
}

/*------------------------------------------------------------------------*/

IMPLEMENT_DYNCREATE(mfc_edit_child, CMDIChildWnd)

BEGIN_MESSAGE_MAP(mfc_edit_child, CMDIChildWnd)
  ON_WM_CLOSE()
END_MESSAGE_MAP()

mfc_edit_child::mfc_edit_child()
{
}

mfc_edit_child::~mfc_edit_child()
{
}

void
mfc_edit_child::OnClose() 
{
  if (this == term_view->GetParent())
    AfxGetApp()->m_pMainWnd->SendMessage(WM_COMMAND, ID_VIEW_TERM, 0);
  else if (this == hist_view->GetParent())
    AfxGetApp()->m_pMainWnd->SendMessage(WM_COMMAND, ID_VIEW_HIST, 0);
  else
    CMDIChildWnd::OnClose();
}

/*------------------------------------------------------------------------*/

IMPLEMENT_DYNCREATE(mfc_edit_doc, CRichEditDoc)

BEGIN_MESSAGE_MAP(mfc_edit_doc, CRichEditDoc)
  //{{AFX_MSG_MAP(mfc_edit_doc)
  ON_COMMAND(ID_FILE_CLOSE, OnFileClose)
  //}}AFX_MSG_MAP
  // Enable default OLE container implementation
  ON_UPDATE_COMMAND_UI(ID_OLE_EDIT_LINKS, CRichEditDoc::OnUpdateEditLinksMenu)
  ON_COMMAND(ID_OLE_EDIT_LINKS, CRichEditDoc::OnEditLinks)
  ON_UPDATE_COMMAND_UI_RANGE(ID_OLE_VERB_FIRST, ID_OLE_VERB_LAST, CRichEditDoc::OnUpdateObjectVerbMenu)
END_MESSAGE_MAP()

mfc_edit_doc::mfc_edit_doc()
{
}

mfc_edit_doc::mfc_edit_doc(CMultiDocTemplate *mdt, int hist)
{
  mdt->AddDocument(this);
  m_pDocTemplate = mdt;
  mfc_edit_child *frame = new mfc_edit_child;
  CCreateContext context;
  context.m_pCurrentFrame = 0;
  context.m_pCurrentDoc = this;
  context.m_pNewViewClass = RUNTIME_CLASS(mfc_term_view);
  context.m_pNewDocTemplate = mdt;
  frame->LoadFrame(IDR_EDITFRAME, WS_OVERLAPPEDWINDOW | FWS_ADDTOTITLE,
                   0, &context);
  SetTitle(hist? "*history*" : "*terminal*");
  OnNewDocument();
  if (!hist)
    ((mfc_term_view *)GetView())->is_term = 1;
}

mfc_edit_doc::~mfc_edit_doc()
{
}

CRichEditCntrItem*
mfc_edit_doc::CreateClientItem(REOBJECT* preo) const
{
  return 0;
}

void
mfc_edit_doc::Serialize(CArchive& ar)
{
  if (ar.IsStoring()) {
    // TODO: add storing code here
  } else {
    // TODO: add loading code here
  }

  // Calling the base class CRichEditDoc enables serialization
  //  of the container document's COleClientItem objects.
  m_bRTF = 0;
  CRichEditDoc::Serialize(ar);
}

BOOL
mfc_edit_doc::OnNewDocument()
{
  if (!CRichEditDoc::OnNewDocument())
    return FALSE;

  // TODO: add reinitialization code here
  // (SDI documents will reuse this document)

  return TRUE;
}

BOOL
mfc_edit_doc::OnOpenDocument(LPCTSTR lpszPathName) 
{
  if (!CRichEditDoc::OnOpenDocument(lpszPathName))
    return FALSE;

  // TODO: Add your specialized creation code here

  return TRUE;
}

BOOL
mfc_edit_doc::SaveModified() 
{
  mfc_edit_view *view = (mfc_edit_view *)GetView();
  if (view==term_view || view==hist_view)
    return 1;
  else
    return CRichEditDoc::SaveModified();
}

void
mfc_edit_doc::OnFileClose() 
{
  mfc_edit_view *view = (mfc_edit_view *)GetView();
  if (view == term_view)
    AfxGetApp()->m_pMainWnd->SendMessage(WM_COMMAND, ID_VIEW_TERM, 0);
  else if (view == hist_view)
    AfxGetApp()->m_pMainWnd->SendMessage(WM_COMMAND, ID_VIEW_HIST, 0);
  else
    COleServerDoc::OnCloseDocument();
}

/*------------------------------------------------------------------------*/

IMPLEMENT_DYNCREATE(mfc_term_view, mfc_edit_view)

BEGIN_MESSAGE_MAP(mfc_term_view, mfc_edit_view)
  ON_COMMAND(ID_EDIT_UNDO, on_undo)
END_MESSAGE_MAP()

mfc_term_view::mfc_term_view()
{
  mark = smin = smax = len = 0;
  is_visible = recalling = recursing = 0;
}

mfc_term_view::~mfc_term_view()
{
}

BOOL
mfc_term_view::PreCreateWindow(CREATESTRUCT& cs)
{
  return CRichEditView::PreCreateWindow(cs);
}

void
mfc_term_view::OnInitialUpdate()
{
  mfc_edit_view::OnInitialUpdate();
}

LRESULT
mfc_term_view::WindowProc(UINT message, WPARAM wParam, LPARAM lParam) 
{
  LRESULT result = 0;
  int key_unshifted = 0;
  long len0 = 0;
  int top = !recursing;

  if (top) {
    if (message==WM_DESTROY) top = 0;
    recursing = 1;
  }

  switch (message) {
  case WM_KEYDOWN:
  case WM_CHAR:
    key_unshifted = (GetKeyState(VK_SHIFT)>=0 && GetKeyState(VK_CONTROL)>=0 &&
                    GetKeyState(VK_MENU)>=0);
    if (key_unshifted) {
      if (wParam == '\r') {
        if (message == WM_KEYDOWN) {
          if (this == term_view) send_or_copy();
          else recall_line();
        }
        if (top) recursing = 0;
        return 0;
      } else if (this==term_view && message==WM_KEYDOWN) {
        if ((wParam==VK_UP || wParam==VK_DOWN) && at_eob()) {
          if (wParam == VK_UP)  recall_prev();
          else recall_next();
          recalling = 1;
          if (top) recursing = 0;
          return 0;
        }
      }
    }
    /* drop through */
  case WM_CUT:
  case WM_PASTE:
  case EM_REPLACESEL:
    get_state();
    result = CRichEditView::WindowProc(message, wParam, lParam);
    if (smin<mark ||
        (message==WM_KEYDOWN && wParam=='\b' && smin==smax && smin==mark)) {
      if (smax > mark) mark = smax;
      mark += eob() - len;
    }
    if (this==term_view && message==WM_KEYDOWN && wParam==VK_HOME)
      home_mark();
    recalling = 0;
    break;

  default:
    if (top && recalling) get_state();
    result = CRichEditView::WindowProc(message, wParam, lParam);
    if (top && recalling) {
      long smin0=smin, smax0=smax, len0=len;
      get_state();
      recalling = (smin0==smin && smax0==smax && len0==len);
    }
  }

  if (top) recursing = 0;
  return result;
}

void
mfc_term_view::on_undo()
{
  /* note: WM_UNDO never sent or never reaches WindowProc
   * -- neither does ID_EDIT_UNDO WM_COMMAND message?? */
  long len0 = eob();
  CRichEditView::OnEditUndo();
  get_state();
  if (smin < mark)
    mark += len - len0;
  recalling = 0;
}

void
mfc_term_view::get_state()
{
  len = eob();
  GetRichEditCtrl().GetSel(smin, smax);
}

long
mfc_term_view::bol(int offset)
{
  long line = GetRichEditCtrl().LineFromChar(-1);
  if (offset) {
    line += offset;
    if (line < 0)
      line = 0;
    else if (line >= GetRichEditCtrl().GetLineCount())
      line = GetRichEditCtrl().GetLineCount() - 1;
  }
  return GetRichEditCtrl().LineIndex(line);
}

long
mfc_term_view::eol(int offset)
{
  long i = bol(offset);
  return i + GetRichEditCtrl().LineLength(i);
}

long
mfc_term_view::eob(int beg)
{
  long i = GetRichEditCtrl().LineIndex(GetRichEditCtrl().GetLineCount() - 1);
  if (!beg) i += GetRichEditCtrl().LineLength(i);
  return i;
}

int
mfc_term_view::at_eob()
{
  long mn, mx;
  GetRichEditCtrl().GetSel(mn, mx);
  return (mn==mx && mx==eob());
}

CString
mfc_term_view::grab_line(int offset) const
{
  long line = GetRichEditCtrl().LineFromChar(-1);
  if (offset) {
    line += offset;
    if (line<0 || line>=GetRichEditCtrl().GetLineCount())
      return (const char *)0;
  }
  int len = GetRichEditCtrl().LineLength(GetRichEditCtrl().LineIndex(line));
  LPSTR s = (char*)_alloca((len + 2)*2);   // copied from MFC GetSelText
  GetRichEditCtrl().GetLine(line, s, len+2);
  s[len] = '\0';
  return s;
}

void
mfc_term_view::select_lines(int off1, int off2)
{
  GetRichEditCtrl().SetSel(bol(off1), bol(off2));
}

void
mfc_term_view::save_line(const char *txt)    // to end of hist_view
{
  long mn0, mx0;
  GetRichEditCtrl().GetSel(mn0, mx0);
  long i = eob(1);
  if (mx0 > i) mn0 = mx0 = i;
  GetRichEditCtrl().SetSel(i, eob(0));
  GetRichEditCtrl().ReplaceSel(txt);
  GetRichEditCtrl().SetSel(mn0, mx0);
}

void
mfc_term_view::set_command(const char *txt)  // at end of term_view
{
  GetRichEditCtrl().SetSel(mark, eob(0));
  GetRichEditCtrl().ReplaceSel(txt);
  long mn, mx;
  GetRichEditCtrl().GetSel(mn, mx);
  GetRichEditCtrl().SetSel(mx, mx);
}

void
mfc_term_view::add_output(const char *txt)   // at end of term_view
{
  if (::IsWindow(m_hWnd)) {
    long mn0, mx0, i, m = mark;
    GetRichEditCtrl().GetSel(mn0, mx0);
    GetRichEditCtrl().SetSel(mark, mark);
    GetRichEditCtrl().ReplaceSel(txt);
    GetRichEditCtrl().GetSel(i, mark);
    if (mx0>m || (mx0==m && mn0==mx0)) {
      if (mn0 < m) mn0 = mx0;
      mn0 += mark-m, mx0 += mark-m;
    }
    GetRichEditCtrl().SetSel(mn0, mx0);
  }
}

void
mfc_term_view::send_or_copy()  // bound to Enter in term_view
{
  long i, m;
  GetRichEditCtrl().GetSel(i, m);
  m = GetRichEditCtrl().LineIndex(GetRichEditCtrl().LineFromChar(mark));
  if (i >= m) {  // send new input as command line
    m = eob();
    //if (GetRichEditCtrl().LineLength(m)) {
      GetRichEditCtrl().SetSel(m, m);
      GetRichEditCtrl().ReplaceSel("\r\n");
      GetRichEditCtrl().GetSel(i, m);
    //}
    GetRichEditCtrl().SetSel(mark, m);
    CString txt = GetRichEditCtrl().GetSelText();
    if (mfc_deliver(txt, txt.GetLength())) {
      hist_view->save_line(txt);
      mark = m;
    }
    GetRichEditCtrl().SetSel(m, m);

  } else {       // copy current line as pending new command line
    CString txt = grab_line();
    set_command(w_deprompt((char *)(const char *)txt));
  }
}

void
mfc_term_view::recall_prev()   // bound to VK_UP in term_view
{
  long i;
  if (!recalling) {
    i = hist_view->eob(1);
    hist_view->GetRichEditCtrl().SetSel(i, i);
  } else {
    i = hist_view->bol();
  }
  long j = hist_view->bol(-1);
  if (j < i) {
    hist_view->GetRichEditCtrl().SetSel(j, j);
    CString txt = hist_view->grab_line();
    set_command(txt);
  } else {
    MessageBeep(MB_OK);
  }
}

void
mfc_term_view::recall_next()   // bound to VK_DOWN in term_view
{
  long i = hist_view->bol(1);
  if (i < hist_view->eob(1)) {
    hist_view->GetRichEditCtrl().SetSel(i, i);
    CString txt = hist_view->grab_line();
    set_command(txt);
  } else {
    MessageBeep(MB_OK);
  }
}

void
mfc_term_view::home_mark()     // after VK_HOME in term_view
{
  long i, m;
  GetRichEditCtrl().GetSel(i, m);
  if (i < mark) {  // adjust to mark
    if (GetRichEditCtrl().LineFromChar(mark) ==
        GetRichEditCtrl().LineFromChar(i))
      GetRichEditCtrl().SetSel(mark, mark);
  }
}

void
mfc_term_view::recall_line()   // bound to Enter in hist_view
{
  CString txt = grab_line();
  term_view->set_command(txt);
  long i = bol(0);     // get first char of current line
  if (i == bol(1)) {  // this line typed at end: send it
    term_view->send_or_copy();
    i = eob();
  } else if (term_view->is_visible) {
    CMDIFrameWnd *fw = (CMDIFrameWnd *)(the_boss.m_pMainWnd);
        fw->MDIActivate(term_view->GetParent());
  }
  GetRichEditCtrl().SetSel(i, i);
}

mfc_term_view *term_view = 0;
mfc_term_view *hist_view = 0;

/*------------------------------------------------------------------------*/

int
mfc_stdinit(void (**wout)(char*,long), void (**werr)(char*,long))
{ /* we are in worker thread here */
  int result = 1;
  *wout = *werr = 0;
  w_sending = CreateEvent(0, 1, 0, 0);
  if (w_sending) {
    the_boss.m_pMainWnd->SendMessage(ID_CALL_FUNC, (WPARAM)0,
                                     (LPARAM)&mfc_term_init);
    w_add_input(w_sending, mfc_sender, 0);
    *wout = mfc_stdout;
    *werr = mfc_stderr;
    result = 0;
  }
  return result;
}

/* ARGSUSED */
static void
mfc_sender(void *context)
{ /* we are in worker thread here */
  ResetEvent(w_sending);  /* set in boss thread to trigger this callback */
  w_deliver(w_sendbuf(-1));
}

static void
mfc_stdout(char *output_line, long len)
{ /* we are in worker thread here */
  if (the_boss.m_pMainWnd)
    the_boss.m_pMainWnd->SendMessage(ID_CALL_FUNC, (WPARAM)output_line,
                                     (LPARAM)&mfc_bstdout);
}

static void
mfc_stderr(char *output_line, long len)
{ /* we are in worker thread here */
  if (the_boss.m_pMainWnd)
    the_boss.m_pMainWnd->SendMessage(ID_CALL_FUNC, (WPARAM)output_line,
                                     (LPARAM)&mfc_bstderr);
}

static void
mfc_bstdout(char *output_line)
{ /* we are in boss thread here */
  term_view->add_output(output_line);
}

static void
mfc_bstderr(char *output_line)
{ /* we are in boss thread here */
  term_view->add_output(output_line);
}

static int
mfc_deliver(const char *buf, long len)
{ /* we are in boss thread here */
  int ready = (WaitForSingleObject(w_sending,0) == WAIT_TIMEOUT);
  if (ready && len>=0) {
    char *line = w_sendbuf(len);
    while (len--) *line++= *buf++;
    line[0] = '\0';
    SetEvent(w_sending);
  }
  return ready;
}

void
mfc_reset_stdin()
{ /* abandon any pending input if SIGINT received
   * -- without this, deadlock is possible
   * -- actually, should be called from wpoll when clearing queue */
  ResetEvent(w_sending);
}

/*------------------------------------------------------------------------*/

void
mfc_edit_view::on_choose_file()
{
  CFileDialog dlg(1);
  CString name;
  dlg.m_ofn.lpstrTitle = "Insert Filename";
  dlg.m_ofn.lpstrFile = name.GetBuffer(1025);
  if (dlg.DoModal() == IDOK)
    GetRichEditCtrl().ReplaceSel(name);
  name.ReleaseBuffer();
}

void
mfc_edit_view::on_goto_out()
{
  int len;
  char s[80];
  long line = GetRichEditCtrl().LineFromChar(-1);
  for (line-- ; line>0 ; line--) {
    len = GetRichEditCtrl().LineLength(GetRichEditCtrl().LineIndex(line));
    if (len>78) len = 78;
    GetRichEditCtrl().GetLine(line, s, len);
    s[len] = '\0';
    if (w_deprompt(s) != s) break;
  }
  line = GetRichEditCtrl().LineIndex(line+1);
  GetRichEditCtrl().SetSel(line, line);
}

void
mfc_edit_view::on_select_out()
{
  int len;
  char s[80];
  long line = GetRichEditCtrl().LineFromChar(-1);
  long linemx = GetRichEditCtrl().GetLineCount();
  long line0 = line;
  long i, j;
  for (line-- ; line>=0 ; line--) {
    len = GetRichEditCtrl().LineLength(GetRichEditCtrl().LineIndex(line));
    if (len>78) len = 78;
    GetRichEditCtrl().GetLine(line, s, len);
    s[len] = '\0';
    if (w_deprompt(s) != s) break;
  }
  i = GetRichEditCtrl().LineIndex(line+1);
  for (line=line0 ; line<linemx ; line++) {
    len = GetRichEditCtrl().LineLength(GetRichEditCtrl().LineIndex(line));
    if (len>78) len = 78;
    GetRichEditCtrl().GetLine(line, s, len);
    s[len] = '\0';
    if (w_deprompt(s) != s) break;
  }
  j = GetRichEditCtrl().LineIndex(line);
  GetRichEditCtrl().SetSel(i, j);
}

void
mfc_edit_view::on_open_line()
{
  int len, i, n, oline=0;
  char s[1076];
  long line = GetRichEditCtrl().LineFromChar(-1);
  for (n=0 ; n<10 ; n++,line--) {
    len = GetRichEditCtrl().LineLength(GetRichEditCtrl().LineIndex(line));
    if (len>1075) len = 1075;
    GetRichEditCtrl().GetLine(line, s, len);
    s[len] = '\0';
    for (i=4 ; i<len-8 ; i++) {
      if (s[i]==':' &&
          s[i-1]=='E' && s[i-2]=='N' && s[i-3]=='I' && s[i-4]=='L') {
        for (i++ ; i<len-7 ; i++) if (s[i]!=' ' && s[i]!='\t') break;
        while (s[i]>='0' && s[i]<='9') oline = 10*oline + (s[i++]-'0');
        for (i+=4 ; i<len-1 ; i++) {
          if (s[i]==':' &&
              s[i-1]=='E' && s[i-2]=='L' && s[i-3]=='I' && s[i-4]=='F') {
            for (i++ ; i<len ; i++) if (s[i]!=' ' && s[i]!='\t') break;
            if (i < len) {
              /* document name s line number oline */
              mfc_edit_doc *doc=
                (mfc_edit_doc *)the_boss.OpenDocumentFile(&s[i]);
              if (doc) {
                mfc_edit_view *view = (mfc_edit_view *)doc->GetView();
                long ll = view->GetRichEditCtrl().LineIndex(oline-1);
                view->GetRichEditCtrl().SetSel(ll, ll);
              } else {
                MessageBeep(MB_OK);
              }
              return;
            }
          }
        }
      }
    }
  }
  MessageBeep(MB_OK);
}

static char *
w_deprompt(char *text)
{
  char *t = text;
  /* remove regexp "^[a-zA-Z0-9_-]*> *" from text */
  while ((t[0]>='a' && t[0]<='z') || (t[0]>='A' && t[0]<='Z') ||
         (t[0]>='0' && t[0]<='9') || t[0]=='_' || t[0]=='-') t++;
  if (t[0]=='>') {
    t++;
    while (t[0]==' ' || t[0]=='\t') t++;
  } else {
    t = text;
  }
  return t;
}
