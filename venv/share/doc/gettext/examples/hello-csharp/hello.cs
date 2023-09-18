// Example for use of GNU gettext.
// This file is in the public domain.
//
// Source code of the C# program.

using System; /* String, Console */
using GNU.Gettext; /* GettextResourceManager */
using System.Diagnostics; /* Process */

public class Hello {
  public static void Main (String[] args) {
    GettextResourceManager catalog =
      new GettextResourceManager("hello-csharp");
    Console.WriteLine(catalog.GetString("Hello, world!"));
    Console.WriteLine(
        String.Format(
            catalog.GetString("This program is running as process number {0}."),
            Process.GetCurrentProcess().Id));
  }
}
