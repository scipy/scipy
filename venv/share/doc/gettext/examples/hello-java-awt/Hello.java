// Example for use of GNU gettext.
// This file is in the public domain.
//
// Source code of the Java/AWT program.

import java.util.*;
import java.io.*;
import java.text.*;
import java.awt.*;
import java.awt.event.*;
import gnu.gettext.*;

public class Hello {
  public static void main (String[] args) {
    ResourceBundle catalog = ResourceBundle.getBundle("hello-java-awt");
    Frame frame = new Frame("Hello example");
    frame.addWindowListener(
        new WindowAdapter() {
          public void windowClosing (WindowEvent event) {
            System.exit(0);
          }
        });
    Label label1 = new Label(GettextResource.gettext(catalog,"Hello, world!"));
    Label label2 =
      new Label(
          MessageFormat.format(
              GettextResource.gettext(catalog,
                  "This program is running as process number {0}."),
              new Object[] { getPid() }));
    Button button = new Button("OK");
    button.addActionListener(
        new ActionListener() {
          public void actionPerformed (ActionEvent event) {
            System.exit(0);
          }
        });
    Container labels = new Container();
    labels.setLayout(new GridLayout(2, 1));
    labels.add(label1);
    labels.add(label2);
    Container buttons = new Container();
    buttons.setLayout(new FlowLayout(FlowLayout.RIGHT));
    buttons.add(button);
    frame.setLayout(new BorderLayout());
    frame.add(labels, BorderLayout.CENTER);
    frame.add(buttons, BorderLayout.SOUTH);
    frame.pack();
    frame.setVisible(true);
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
