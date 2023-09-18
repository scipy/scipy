// Example for use of GNU gettext.
// This file is in the public domain.
//
// Source code of the Java/Swing program.

import java.util.*;
import java.io.*;
import java.text.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import gnu.gettext.*;

public class Hello {
  public static void main (String[] args) {
    ResourceBundle catalog = ResourceBundle.getBundle("hello-java-swing");
    JFrame frame = new JFrame("Hello example");
    frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    JLabel label1 =
      new JLabel(GettextResource.gettext(catalog,"Hello, world!"));
    JLabel label2 =
      new JLabel(
          MessageFormat.format(
              GettextResource.gettext(catalog,
                  "This program is running as process number {0}."),
              new Object[] { getPid() }));
    JButton button = new JButton("OK");
    button.addActionListener(
        new ActionListener() {
          public void actionPerformed (ActionEvent event) {
            System.exit(0);
          }
        });
    JPanel labels = new JPanel();
    labels.setLayout(new GridLayout(2, 1));
    labels.add(label1);
    labels.add(label2);
    JPanel buttons = new JPanel();
    buttons.setLayout(new FlowLayout(FlowLayout.RIGHT));
    buttons.add(button);
    frame.getContentPane().setLayout(new BorderLayout());
    frame.getContentPane().add(labels, BorderLayout.CENTER);
    frame.getContentPane().add(buttons, BorderLayout.SOUTH);
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
