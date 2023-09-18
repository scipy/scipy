// Example for use of GNU gettext.
// This file is in the public domain.
//
// Source code of the C#/Forms program.

using System; /* String, EventHandler */
using GNU.Gettext; /* GettextResourceManager */
using System.Diagnostics; /* Process */
using System.Threading; /* Thread */
using System.Drawing; /* Point, Size */
using System.Windows.Forms; /* Application, Form, Label, Button */

public class Hello {

  private static GettextResourceManager catalog =
    new GettextResourceManager("hello-csharp-forms");

  class HelloWindow : Form {

    private int border;
    private Label label1;
    private Label label2;
    private Button ok;

    public HelloWindow () {
      border = 2;

      label1 = new Label();
      label1.Text = catalog.GetString("Hello, world!");
      label1.ClientSize = new Size(label1.PreferredWidth, label1.PreferredHeight);
      Controls.Add(label1);

      label2 = new Label();
      label2.Text =
        String.Format(
            catalog.GetString("This program is running as process number {0}."),
            Process.GetCurrentProcess().Id);
      label2.ClientSize = new Size(label2.PreferredWidth, label2.PreferredHeight);
      Controls.Add(label2);

      ok = new Button();
      Label okLabel = new Label();
      ok.Text = okLabel.Text = "OK";
      ok.ClientSize = new Size(okLabel.PreferredWidth + 12, okLabel.PreferredHeight + 4);
      ok.Click += new EventHandler(Quit);
      Controls.Add(ok);

      Size total = ComputePreferredSizeWithoutBorder();
      LayoutControls(total.Width, total.Height);
      ClientSize = new Size(border + total.Width + border, border + total.Height + border);
    }

    protected override void OnResize(EventArgs ev) {
      LayoutControls(ClientSize.Width - border - border, ClientSize.Height - border - border);
      base.OnResize(ev);
    }

    // Layout computation, part 1: The preferred size of this panel.
    private Size ComputePreferredSizeWithoutBorder () {
      int totalWidth = Math.Max(Math.Max(label1.PreferredWidth, label2.PreferredWidth),
                                ok.Width);
      int totalHeight = label1.PreferredHeight + label2.PreferredHeight + 6 + ok.Height;
      return new Size(totalWidth, totalHeight);
    }

    // Layout computation, part 2: Determine where to put the sub-controls.
    private void LayoutControls (int totalWidth, int totalHeight) {
      label1.Location = new Point(border, border);
      label2.Location = new Point(border, border + label1.PreferredHeight);
      ok.Location = new Point(border + totalWidth - ok.Width, border + totalHeight - ok.Height);
    }

    private void Quit (Object sender, EventArgs ev) {
      Application.Exit();
    }
  }

  public static void Main () {
    Application.Run(new HelloWindow());
  }
}
