/*
 * mfcapp.cpp -- $Id$
 * MFC implementation of play, boss and worker thread declarations
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#define VC_EXTRALEAN

#include <afxwin.h>
#include <afxcmn.h>

class mfc_boss : public CWinApp {
public:
  mfc_boss(int (*on_launch)(int, char **));
  virtual BOOL InitInstance();
  virtual int ExitInstance();

  void on_sigint();
  void on_view_term();
  void on_view_hist();
  void on_update_view_term(CCmdUI *ui);
  void on_update_view_hist(CCmdUI *ui);

  afx_msg void on_about();
  DECLARE_MESSAGE_MAP()
};

extern mfc_boss the_boss;
