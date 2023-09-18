// Example for use of GNU gettext.
// This file is in the public domain.

// Source code of the C++ program.

#include <wx/wx.h>
#include <wx/intl.h>

/* Get getpid() declaration.  */
#if HAVE_UNISTD_H
# include <unistd.h>
#endif

class MyApp: public wxApp
{
public:
  virtual bool OnInit();
private:
  // wxWidgets has the concept of a "current locale". It is the one returned
  // by wxGetLocale() and implicitly used by wxGetTranslation.
  // But there is no way to explicitly set this current locale! Rather, it is
  // always set to the last constructed locale(!), and is modified when a
  // locale is destroyed. In such a way that the current locale points to
  // invalid memory after you do
  //    wxLocale *a = new wxLocale;
  //    wxLocale *b = new wxLocale;
  //    delete a;
  //    delete b;
  // So, to avoid problems, we use exactly one instance of wxLocale, and keep
  // it alive for the entire application lifetime.
  wxLocale appLocale;
};

class MyFrame: public wxFrame
{
public:
  MyFrame();
};

// This defines the main() function.
IMPLEMENT_APP(MyApp)

bool MyApp::OnInit()
{
  // First, register the base directory where to look up .mo files.
  wxLocale::AddCatalogLookupPathPrefix(wxT(LOCALEDIR));
  // Second, initialize the locale and set the application-wide message domain.
  appLocale.Init();
  appLocale.AddCatalog(wxT("hello-c++-wxwidgets"));
  // Now wxGetLocale() is initialized appropriately.

  // Then only start building the GUI elements of the application.

  // Create the main frame window.
  MyFrame *frame = new MyFrame();

  // Show the frame.
  frame->Show(true);
  SetTopWindow(frame);

  return true;
}

MyFrame::MyFrame()
  : wxFrame(NULL, wxID_ANY, _T("Hello example"))
{
  wxStaticText *label1 =
    new wxStaticText(this, wxID_ANY, _("Hello, world!"));

  wxString label2text =
    wxString::Format(_("This program is running as process number %d."),
                     getpid());
  wxStaticText *label2 =
    new wxStaticText(this, wxID_ANY, label2text);

  wxBoxSizer *topSizer = new wxBoxSizer(wxVERTICAL);
  topSizer->Add(label1);
  topSizer->Add(label2);
  SetSizer(topSizer);
}
