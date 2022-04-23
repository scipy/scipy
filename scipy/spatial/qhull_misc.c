#include "qhull_src/src/qhull_ra.h"


/* This is a patched version of qhull_src/src/user_r.c:qh_new_qhull,
   with an additional "feaspoint" argument.

   See qhull_src/README.scipy
*/
int qh_new_qhull_scipy(qhT *qh, int dim, int numpoints, coordT *points, boolT ismalloc,
                char *qhull_cmd, FILE *outfile, FILE *errfile, coordT* feaspoint) {
  /* gcc may issue a "might be clobbered" warning for dim, points, and ismalloc [-Wclobbered].
     These parameters are not referenced after a longjmp() and hence not clobbered.
     See http://stackoverflow.com/questions/7721854/what-sense-do-these-clobbered-variable-warnings-make */
  int exitcode, hulldim;
  boolT new_ismalloc;
  coordT *new_points;

  if(!errfile){
    errfile= stderr;
  }
  if (!qh->qhmem.ferr) {
    qh_meminit(qh, errfile);
  } else {
    qh_memcheck(qh);
  }
  if (strncmp(qhull_cmd, "qhull ", (size_t)6) && strcmp(qhull_cmd, "qhull") != 0) {
    qh_fprintf(qh, errfile, 6186, "qhull error (qh_new_qhull): start qhull_cmd argument with \"qhull \" or set to \"qhull\"\n");
    return qh_ERRinput;
  }
  qh_initqhull_start(qh, NULL, outfile, errfile);
  if(numpoints==0 && points==NULL){
      trace1((qh, qh->ferr, 1047, "qh_new_qhull: initialize Qhull\n"));
      return 0;
  }
  trace1((qh, qh->ferr, 1044, "qh_new_qhull: build new Qhull for %d %d-d points with %s\n", numpoints, dim, qhull_cmd));
  exitcode= setjmp(qh->errexit);
  if (!exitcode){
    qh->NOerrexit= False;
    qh_initflags(qh, qhull_cmd);
    if (qh->DELAUNAY)
      qh->PROJECTdelaunay= True;
    if (qh->HALFspace) {
      /* points is an array of halfspaces,
         the last coordinate of each halfspace is its offset */
      hulldim= dim-1;
      if(feaspoint)
      {
        coordT* coords;
        coordT* value;
        int i;
        if (!(qh->feasible_point= (pointT*)qh_malloc(hulldim * sizeof(coordT)))) {
          qh_fprintf(qh, qh->ferr, 6079, "qhull error: insufficient memory for 'Hn,n,n'\n");
          qh_errexit(qh, qh_ERRmem, NULL, NULL);
        }
        coords = qh->feasible_point;
        value = feaspoint;
        for(i = 0; i < hulldim; ++i)
        {
          *(coords++) = *(value++);
        }
      }
      else
      {
        qh_setfeasible(qh, hulldim);
      }
      new_points= qh_sethalfspace_all(qh, dim, numpoints, points, qh->feasible_point);
      new_ismalloc= True;
      if (ismalloc)
        qh_free(points);
    }else {
      hulldim= dim;
      new_points= points;
      new_ismalloc= ismalloc;
    }
    qh_init_B(qh, new_points, numpoints, hulldim, new_ismalloc);
    qh_qhull(qh);
    qh_check_output(qh);
    if (outfile) {
      qh_produce_output(qh);
    }else {
      qh_prepare_output(qh);
    }
    if (qh->VERIFYoutput && !qh->FORCEoutput && !qh->STOPadd && !qh->STOPcone && !qh->STOPpoint)
      qh_check_points(qh);
  }
  qh->NOerrexit= True;
  return exitcode;
} /* new_qhull */
