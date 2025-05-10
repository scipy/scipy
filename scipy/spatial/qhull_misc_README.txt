The file scipy/spatial/qhull_misc.c contains a function
"qh_new_qhull_scipy" derived from Qhull sources, via the following
patch:

--- a/subprojects/qhull_r/libqhull_r/user_r.c
+++ b/subprojects/qhull_r/libqhull_r/user_r.c
@@ -122,7 +122,7 @@
     An example of using qh_new_qhull is user_eg_r.c
 */
 int qh_new_qhull(qhT *qh, int dim, int numpoints, coordT *points, boolT ismalloc,
-                char *qhull_cmd, FILE *outfile, FILE *errfile) {
+                char *qhull_cmd, FILE *outfile, FILE *errfile, coordT* feaspoint) {
   /* gcc may issue a "might be clobbered" warning for dim, points, and ismalloc [-Wclobbered].
      These parameters are not referenced after a longjmp() and hence not clobbered.
      See http://stackoverflow.com/questions/7721854/what-sense-do-these-clobbered-variable-warnings-make */
@@ -159,6 +159,26 @@ int qh_new_qhull(qhT *qh, int dim, int numpoints, coordT *points, boolT ismalloc
          the last coordinate of each halfspace is its offset */
       hulldim= dim-1;
-      qh_setfeasible(qh, hulldim);
+      if(feaspoint)
+      {
+        coordT* coords;
+        coordT* value;
+        int i;
+        if (!(qh->feasible_point= (pointT*)qh_malloc(hulldim * sizeof(coordT)))) {
+          qh_fprintf(qh, qh->ferr, 6079, "qhull error: insufficient memory for 'Hn,n,n'\n");
+          qh_errexit(qh, qh_ERRmem, NULL, NULL);
+        }
+        coords = qh->feasible_point;
+        value = feaspoint;
+        for(i = 0; i < hulldim; ++i)
+        {
+          *(coords++) = *(value++);
+        }
+      }
+      else
+      {
+        qh_setfeasible(qh, hulldim);
+      }
       new_points= qh_sethalfspace_all(qh, dim, numpoints, points, qh->feasible_point);
       new_ismalloc= True;
       if (ismalloc)
