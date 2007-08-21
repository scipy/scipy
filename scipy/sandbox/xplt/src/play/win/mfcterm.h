/*
 * mfcterm.h
 * richedit class for play MDI development environment
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

/*
 * mfc_edit_view -- editor windows for yorick source code, data files
 * 1. uses monospaced font (OnCreate)
 * 2. prevents MFC from allowing files and other OLE/COM objects
 *    from being pasted into files (QueryAcceptData)
 * 3. turns off wordwrap to number lines correctly (constructor)
 * 4. intercepts WM_DESTROY incomprehensibly ala MSVC++ (OnDestroy)
 * X. experiment with richedit 2.0/3.0 (PreCreateWindow)
 *
 * mfc_edit_child -- MDIChild parent of mfc_edit_view
 * 1. intercepts WM_CLOSE to prevent killing terminal, history windows
 *
 * mfc_edit_doc -- editor documents
 * 1. forces files to be text format, not RTF (Serialize)
 * 2. prevents save dialog box for terminal, history windows (SaveModified)
 * 3. intercepts ID_FILE_CLOSE to prevent killing terminal, history windows
 *
 * mfc_term_view -- terminal and history windows
 * 1. maintains "process mark" separating old input and output
 *    from new input in terminal window - allows typeahead
 *   -intercepts EN_SELCHANGE for this purpose (on_selchange, WindowProc)
 * 2. new bindings for Enter, Up and Down arrows (latter only at eob),
 *    extra processing for Delete and Home keys, cut, paste, undo, replace
 *    messages (on_selchange, WindowProc)
 */

/*------------------------------------------------------------------------*/

/* assumes mfcapp.h included first */
#include <afxrich.h>

class mfc_edit_doc;

class mfc_edit_view : public CRichEditView
{
  DECLARE_DYNCREATE(mfc_edit_view)
public:
  mfc_edit_view();
  virtual ~mfc_edit_view();
  int is_term, riched2;

  mfc_edit_doc* GetDocument() { return (mfc_edit_doc*)m_pDocument; }
  virtual BOOL PreCreateWindow(CREATESTRUCT& cs);

protected:
  virtual void OnInitialUpdate(); // called first time after construct
  virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
  virtual HRESULT QueryAcceptData(LPDATAOBJECT lpdataobj,
                                  CLIPFORMAT* lpcfFormat, DWORD reco,
                                  BOOL fReally, HGLOBAL hMetaPict);
  int OnCreate(LPCREATESTRUCT lpcs);  // on WM_CREATE message

  afx_msg void on_choose_file();
  afx_msg void on_goto_line();
  afx_msg void on_goto_out();
  afx_msg void on_select_out();
  afx_msg void on_open_line();
  afx_msg void on_redo();
  afx_msg void on_update_term(CCmdUI *ui);
  afx_msg void on_update_redo(CCmdUI *ui);
  afx_msg void OnDestroy();

  DECLARE_MESSAGE_MAP()
};

/*------------------------------------------------------------------------*/

class mfc_edit_child : public CMDIChildWnd
{
  DECLARE_DYNCREATE(mfc_edit_child)
public:
  mfc_edit_child();
  virtual ~mfc_edit_child();

protected:
  afx_msg void OnClose();

  DECLARE_MESSAGE_MAP()
};

/*------------------------------------------------------------------------*/

class mfc_edit_doc : public CRichEditDoc
{
  DECLARE_DYNCREATE(mfc_edit_doc)
public:
  mfc_edit_doc();
  virtual ~mfc_edit_doc();
  mfc_edit_doc(CMultiDocTemplate *mdt, int hist);

  virtual BOOL OnNewDocument();
  virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
  virtual void Serialize(CArchive& ar);
protected:
  virtual BOOL SaveModified();
  virtual CRichEditCntrItem* CreateClientItem(REOBJECT* preo) const;

  afx_msg void OnFileClose();

  DECLARE_MESSAGE_MAP()
};

/*------------------------------------------------------------------------*/

class mfc_term_view : public mfc_edit_view
{
  DECLARE_DYNCREATE(mfc_term_view)
public:
  mfc_term_view();
  virtual ~mfc_term_view();

  mfc_edit_doc* GetDocument() { return (mfc_edit_doc*)m_pDocument; }

  int is_visible;

  long bol(int offset = 0);
  long eol(int offset = 0);
  long eob(int beg = 0);
  int at_eob();
  CString grab_line(int offset = 0) const;
  void select_lines(int off1, int off2);
  void get_state();
  long smin, smax, len;  // set by get_state
  int recursing;         // set if WindowProc recursing
  int recalling;         // set during command recall sequences

  void save_line(const char *txt);   // to end of hist_view

  void set_command(const char *txt); // at end of term_view
  void add_output(const char *txt);  // at end of term_view
  void send_or_copy();  // bound to Enter in term_view
  void recall_prev();   // bound to VK_UP in term_view
  void recall_next();   // bound to VK_DOWN in term_view
  void home_mark();     // called after VK_HOME in term_view

  void recall_line();   // bound to Enter in hist_view

  // TODO functions to:
  // 1. select output block (back to previous prompt)
  // 2. extract LINE/FILE, open and goto
  // 3. goto line (in mfc_edit_view)

public:
  virtual BOOL PreCreateWindow(CREATESTRUCT& cs);

protected:
  virtual void OnInitialUpdate(); // called first time after construct
  virtual LRESULT WindowProc(UINT message, WPARAM wParam, LPARAM lParam);

  long mark;  // process mark -- insertion point for p_stdout

  afx_msg void on_undo();
  DECLARE_MESSAGE_MAP()
};

/*------------------------------------------------------------------------*/

extern mfc_term_view *term_view, *hist_view;
extern void mfc_term_init(void *arg);
extern int mfc_stdinit(void (**wout)(char*,long), void (**werr)(char*,long));
extern void mfc_reset_stdin();
