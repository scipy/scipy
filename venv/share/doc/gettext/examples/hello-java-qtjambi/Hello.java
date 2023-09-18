// Example for use of GNU gettext.
// This file is in the public domain.
//
// Source code of the Java/QtJambi program.

import java.util.*;
import java.io.*;
import java.text.*;
import com.trolltech.qt.core.*;
import com.trolltech.qt.gui.*;
import gnu.gettext.*;

public class Hello {
  public static void main (String[] args) {
    ResourceBundle catalog = ResourceBundle.getBundle("hello-java-qtjambi");

    QApplication.initialize(args);

    QMainWindow window = new QMainWindow();
    window.setWindowTitle("Hello example");

    QWidget panel = new QWidget();
    QVBoxLayout panelLayout = new QVBoxLayout();
    panelLayout.setSpacing(2);

    QLabel label1 =
      new QLabel(GettextResource.gettext(catalog,"Hello, world!"));
    panelLayout.addWidget(label1);

    QLabel label2 =
      new QLabel(
          MessageFormat.format(
              GettextResource.gettext(catalog,
                  "This program is running as process number {0}."),
              new Object[] { getPid() }));
    panelLayout.addWidget(label2);

    QWidget buttonBar = new QWidget();
    QHBoxLayout buttonBarLayout = new QHBoxLayout();
    QWidget filler = new QWidget(); // makes the button right-aligned
    buttonBarLayout.addWidget(filler);
    QPushButton button = new QPushButton("OK");
    button.setMaximumWidth(button.sizeHint().width()+20);
    button.clicked.connect(window, "close()");
    buttonBarLayout.addWidget(button);
    buttonBar.setLayout(buttonBarLayout);
    panelLayout.addWidget(buttonBar);

    panel.setLayout(panelLayout);

    window.setCentralWidget(panel);

    window.show();

    QApplication.exec();
  }

  /* Return the process ID of the current process.  */
  private static String getPid () {
    try {
      String[] args = new String[] { "/bin/sh", "-c", "echo $PPID" };
      Process p = Runtime.getRuntime().exec(args);
      InputStream p_out = p.getInputStream();
      String s = (new BufferedReader(new InputStreamReader(p_out))).readLine();
      p.destroy();
      if (s != null)
        return s;
    } catch (IOException e) {
      e.printStackTrace();
    }
    return "???";
  }
}
