// mcs -out:cstest call_highs_from_csharp.cs -r:highscslib.dll

using System;

class Program {
   static void Main(string[] args) {
      double[] cc = {1, -2};
      double[] cl = {0, 0};
      double[] cu = {10, 10};
      double[] rl = {0, 0};
      double[] ru = {2, 1};
      int[] astart = {0, 2, 4};
      int[] aindex = {0, 1, 0, 1};
      double[] avalue = {1, 2, 1, 3};

      HighsModel model = new HighsModel(cc, cl, cu, rl, ru, astart, aindex, avalue);

      HighsLpSolver solver = new HighsLpSolver();

      HighsStatus status = solver.passLp(model);
      status = solver.run();
      HighsSolution sol = solver.getSolution();
      HighsBasis bas = solver.getBasis();
      HighsModelStatus modelStatus = solver.GetModelStatus();
      
      Console.WriteLine("Status: " + status);
      Console.WriteLine("Modelstatus: " + modelStatus);
   
      for (int i=0; i<sol.rowvalue.Length; i++) {
         Console.WriteLine("Activity for row " + i + " = " + sol.rowvalue[i]);
      }
      for (int i=0; i<sol.coldual.Length; i++) {
         Console.WriteLine("Reduced cost x[" + i + "] = " + sol.coldual[i]);
      }
      for (int i=0; i<sol.rowdual.Length; i++) {
         Console.WriteLine("Dual value for row " + i + " = " + sol.rowdual[i]);
      }
      for (int i=0; i<sol.colvalue.Length; i++) {
         Console.WriteLine("x" + i + " = " + sol.colvalue[i] + " is " + bas.colbasisstatus[i]);
      }
   }
}