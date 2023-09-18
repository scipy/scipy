// Test for Java 1.5 or newer.
// This file is in the public domain.
import java.util.*;
/**
 * @author Bruno Haible
 */
public class Test15 {
  public static void main (String[] args) {
    try {
      foo();
    } catch (Throwable e) {
      System.exit(1);
    }
    // Check the JVM version is at least 1.5.
    String version = System.getProperty("java.specification.version");
    int i = 0;
    while (i < version.length()
           && (Character.isDigit(version.charAt(i)) || version.charAt(i)=='.'))
      i++;
    float fversion = Float.valueOf(version.substring(0,i));
    if (!(fversion >= 1.5f)) System.exit(1);
    // Check the VM is not GNU libgcj.
    String vm = System.getProperty("java.vm.name");
    if (vm.startsWith("GNU")) System.exit(1);
    System.exit(0);
  }
  private static List<Integer> foo() {
    return new ArrayList<Integer>();
  }
}
