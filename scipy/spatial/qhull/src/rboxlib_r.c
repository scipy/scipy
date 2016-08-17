/*<html><pre>  -<a                             href="index_r.htm#TOC"
  >-------------------------------</a><a name="TOP">-</a>

   rboxlib_r.c
     Generate input points

   notes:
     For documentation, see prompt[] of rbox_r.c
     50 points generated for 'rbox D4'

   WARNING:
     incorrect range if qh_RANDOMmax is defined wrong (user_r.h)
*/

#include "libqhull_r.h"  /* First for user_r.h */
#include "random_r.h"

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <setjmp.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _MSC_VER  /* Microsoft Visual C++ */
#pragma warning( disable : 4706)  /* assignment within conditional expression. */
#pragma warning( disable : 4996)  /* this function (strncat,sprintf,strcpy) or variable may be unsafe. */
#endif

#define MAXdim 200
#define PI 3.1415926535897932384

/* ------------------------------ prototypes ----------------*/
int qh_roundi(qhT *qh, double a);
void qh_out1(qhT *qh, double a);
void qh_out2n(qhT *qh, double a, double b);
void qh_out3n(qhT *qh, double a, double b, double c);
void qh_outcoord(qhT *qh, int iscdd, double *coord, int dim);
void qh_outcoincident(qhT *qh, int coincidentpoints, double radius, int iscdd, double *coord, int dim);

void    qh_fprintf_rbox(qhT *qh, FILE *fp, int msgcode, const char *fmt, ... );
void    qh_free(void *mem);
void   *qh_malloc(size_t size);
int     qh_rand(qhT *qh);
void    qh_srand(qhT *qh, int seed);

/*-<a                             href="qh-qhull_r.htm#TOC"
  >-------------------------------</a><a name="rboxpoints">-</a>

  qh_rboxpoints(qh, rbox_command )
    Generate points to qh->fout according to rbox options
    Report errors on qh->ferr

  returns:
    0 (qh_ERRnone) on success
    1 (qh_ERRinput) on input error
    4 (qh_ERRmem) on memory error
    5 (qh_ERRqhull) on internal error

  notes:
    To avoid using stdio, redefine qh_malloc, qh_free, and qh_fprintf_rbox (user_r.c)

  design:
    Straight line code (consider defining a struct and functions):

    Parse arguments into variables
    Determine the number of points
    Generate the points
*/
int qh_rboxpoints(qhT *qh, char* rbox_command) {
  int i,j,k;
  int gendim;
  int coincidentcount=0, coincidenttotal=0, coincidentpoints=0;
  int cubesize, diamondsize, seed=0, count, apex;
  int dim=3 , numpoints= 0, totpoints, addpoints=0;
  int issphere=0, isaxis=0,  iscdd= 0, islens= 0, isregular=0, iswidth=0, addcube=0;
  int isgap=0, isspiral=0, NOcommand= 0, adddiamond=0;
  int israndom=0, istime=0;
  int isbox=0, issimplex=0, issimplex2=0, ismesh=0;
  double width=0.0, gap=0.0, radius=0.0, coincidentradius=0.0;
  double coord[MAXdim], offset, meshm=3.0, meshn=4.0, meshr=5.0;
  double *coordp, *simplex= NULL, *simplexp;
  int nthroot, mult[MAXdim];
  double norm, factor, randr, rangap, lensangle= 0, lensbase= 1;
  double anglediff, angle, x, y, cube= 0.0, diamond= 0.0;
  double box= qh_DEFAULTbox; /* scale all numbers before output */
  double randmax= qh_RANDOMmax;
  char command[200], seedbuf[200];
  char *s= command, *t, *first_point= NULL;
  time_t timedata;
  int exitcode;

  exitcode= setjmp(qh->rbox_errexit);
  if (exitcode) {
    /* same code for error exit and normal return, qh->NOerrexit is set */
    if (simplex)
        qh_free(simplex);
    return exitcode;
  }

  *command= '\0';
  strncat(command, rbox_command, sizeof(command)-strlen(command)-1);

  while (*s && !isspace(*s))  /* skip program name */
    s++;
  while (*s) {
    while (*s && isspace(*s))
      s++;
    if (*s == '-')
      s++;
    if (!*s)
      break;
    if (isdigit(*s)) {
      numpoints= qh_strtol(s, &s);
      continue;
    }
    /* ============= read flags =============== */
    switch (*s++) {
    case 'c':
      addcube= 1;
      t= s;
      while (isspace(*t))
        t++;
      if (*t == 'G')
        cube= qh_strtod(++t, &s);
      break;
    case 'd':
      adddiamond= 1;
      t= s;
      while (isspace(*t))
        t++;
      if (*t == 'G')
        diamond= qh_strtod(++t, &s);
      break;
    case 'h':
      iscdd= 1;
      break;
    case 'l':
      isspiral= 1;
      break;
    case 'n':
      NOcommand= 1;
      break;
    case 'r':
      isregular= 1;
      break;
    case 's':
      issphere= 1;
      break;
    case 't':
      istime= 1;
      if (isdigit(*s)) {
        seed= qh_strtol(s, &s);
        israndom= 0;
      }else
        israndom= 1;
      break;
    case 'x':
      issimplex= 1;
      break;
    case 'y':
      issimplex2= 1;
      break;
    case 'z':
      qh->rbox_isinteger= 1;
      break;
    case 'B':
      box= qh_strtod(s, &s);
      isbox= 1;
      break;
    case 'C':
      if (*s)
        coincidentpoints=  qh_strtol(s, &s);
      if (*s == ',') {
        ++s;
        coincidentradius=  qh_strtod(s, &s);
      }
      if (*s == ',') {
        ++s;
        coincidenttotal=  qh_strtol(s, &s);
      }
      if (*s && !isspace(*s)) {
        qh_fprintf_rbox(qh, qh->ferr, 7080, "rbox error: arguments for 'Cn,r,m' are not 'int', 'float', and 'int'.  Remaining string is '%s'\n", s);
        qh_errexit_rbox(qh, qh_ERRinput);
      }
      if (coincidentpoints==0){
        qh_fprintf_rbox(qh, qh->ferr, 6268, "rbox error: missing arguments for 'Cn,r,m' where n is the number of coincident points, r is the radius (default 0.0), and m is the number of points\n");
        qh_errexit_rbox(qh, qh_ERRinput);
      }
      if (coincidentpoints<0 || coincidenttotal<0 || coincidentradius<0.0){
        qh_fprintf_rbox(qh, qh->ferr, 6269, "rbox error: negative arguments for 'Cn,m,r' where n (%d) is the number of coincident points, m (%d) is the number of points, and r (%.2g) is the radius (default 0.0)\n", coincidentpoints, coincidenttotal, coincidentradius);
        qh_errexit_rbox(qh, qh_ERRinput);
      }
      break;
    case 'D':
      dim= qh_strtol(s, &s);
      if (dim < 1
      || dim > MAXdim) {
        qh_fprintf_rbox(qh, qh->ferr, 6189, "rbox error: dimension, D%d, out of bounds (>=%d or <=0)", dim, MAXdim);
        qh_errexit_rbox(qh, qh_ERRinput);
      }
      break;
    case 'G':
      if (isdigit(*s))
        gap= qh_strtod(s, &s);
      else
        gap= 0.5;
      isgap= 1;
      break;
    case 'L':
      if (isdigit(*s))
        radius= qh_strtod(s, &s);
      else
        radius= 10;
      islens= 1;
      break;
    case 'M':
      ismesh= 1;
      if (*s)
        meshn= qh_strtod(s, &s);
      if (*s == ',') {
        ++s;
        meshm= qh_strtod(s, &s);
      }else
        meshm= 0.0;
      if (*s == ',') {
        ++s;
        meshr= qh_strtod(s, &s);
      }else
        meshr= sqrt(meshn*meshn + meshm*meshm);
      if (*s && !isspace(*s)) {
        qh_fprintf_rbox(qh, qh->ferr, 7069, "rbox warning: assuming 'M3,4,5' since mesh args are not integers or reals\n");
        meshn= 3.0, meshm=4.0, meshr=5.0;
      }
      break;
    case 'O':
      qh->rbox_out_offset= qh_strtod(s, &s);
      break;
    case 'P':
      if (!first_point)
        first_point= s-1;
      addpoints++;
      while (*s && !isspace(*s))   /* read points later */
        s++;
      break;
    case 'W':
      width= qh_strtod(s, &s);
      iswidth= 1;
      break;
    case 'Z':
      if (isdigit(*s))
        radius= qh_strtod(s, &s);
      else
        radius= 1.0;
      isaxis= 1;
      break;
    default:
      qh_fprintf_rbox(qh, qh->ferr, 7070, "rbox error: unknown flag at %s.\nExecute 'rbox' without arguments for documentation.\n", s);
      qh_errexit_rbox(qh, qh_ERRinput);
    }
    if (*s && !isspace(*s)) {
      qh_fprintf_rbox(qh, qh->ferr, 7071, "rbox error: missing space between flags at %s.\n", s);
      qh_errexit_rbox(qh, qh_ERRinput);
    }
  }

  /* ============= defaults, constants, and sizes =============== */
  if (qh->rbox_isinteger && !isbox)
    box= qh_DEFAULTzbox;
  if (addcube) {
    cubesize= (int)floor(ldexp(1.0,dim)+0.5);
    if (cube == 0.0)
      cube= box;
  }else
    cubesize= 0;
  if (adddiamond) {
    diamondsize= 2*dim;
    if (diamond == 0.0)
      diamond= box;
  }else
    diamondsize= 0;
  if (islens) {
    if (isaxis) {
        qh_fprintf_rbox(qh, qh->ferr, 6190, "rbox error: can not combine 'Ln' with 'Zn'\n");
        qh_errexit_rbox(qh, qh_ERRinput);
    }
    if (radius <= 1.0) {
        qh_fprintf_rbox(qh, qh->ferr, 6191, "rbox error: lens radius %.2g should be greater than 1.0\n",
               radius);
        qh_errexit_rbox(qh, qh_ERRinput);
    }
    lensangle= asin(1.0/radius);
    lensbase= radius * cos(lensangle);
  }

  if (!numpoints) {
    if (issimplex2)
        ; /* ok */
    else if (isregular + issimplex + islens + issphere + isaxis + isspiral + iswidth + ismesh) {
        qh_fprintf_rbox(qh, qh->ferr, 6192, "rbox error: missing count\n");
        qh_errexit_rbox(qh, qh_ERRinput);
    }else if (adddiamond + addcube + addpoints)
        ; /* ok */
    else {
        numpoints= 50;  /* ./rbox D4 is the test case */
        issphere= 1;
    }
  }
  if ((issimplex + islens + isspiral + ismesh > 1)
  || (issimplex + issphere + isspiral + ismesh > 1)) {
    qh_fprintf_rbox(qh, qh->ferr, 6193, "rbox error: can only specify one of 'l', 's', 'x', 'Ln', or 'Mn,m,r' ('Ln s' is ok).\n");
    qh_errexit_rbox(qh, qh_ERRinput);
  }
  if (coincidentpoints>0 && (numpoints == 0 || coincidenttotal > numpoints)) {
    qh_fprintf_rbox(qh, qh->ferr, 6270, "rbox error: 'Cn,r,m' requested n coincident points for each of m points.  Either there is no points or m (%d) is greater than the number of points (%d).\n", coincidenttotal, numpoints);
    qh_errexit_rbox(qh, qh_ERRinput);
  }
  if (coincidenttotal == 0)
    coincidenttotal= numpoints;

  /* ============= print header with total points =============== */
  if (issimplex || ismesh)
    totpoints= numpoints;
  else if (issimplex2)
    totpoints= numpoints+dim+1;
  else if (isregular) {
    totpoints= numpoints;
    if (dim == 2) {
        if (islens)
          totpoints += numpoints - 2;
    }else if (dim == 3) {
        if (islens)
          totpoints += 2 * numpoints;
      else if (isgap)
        totpoints += 1 + numpoints;
      else
        totpoints += 2;
    }
  }else
    totpoints= numpoints + isaxis;
  totpoints += cubesize + diamondsize + addpoints;
  totpoints += coincidentpoints*coincidenttotal;

  /* ============= seed randoms =============== */
  if (istime == 0) {
    for (s=command; *s; s++) {
      if (issimplex2 && *s == 'y') /* make 'y' same seed as 'x' */
        i= 'x';
      else
        i= *s;
      seed= 11*seed + i;
    }
  }else if (israndom) {
    seed= (int)time(&timedata);
    sprintf(seedbuf, " t%d", seed);  /* appends an extra t, not worth removing */
    strncat(command, seedbuf, sizeof(command)-strlen(command)-1);
    t= strstr(command, " t ");
    if (t)
      strcpy(t+1, t+3); /* remove " t " */
  } /* else, seed explicitly set to n */
  qh_RANDOMseed_(qh, seed);

  /* ============= print header =============== */

  if (iscdd)
      qh_fprintf_rbox(qh, qh->fout, 9391, "%s\nbegin\n        %d %d %s\n",
      NOcommand ? "" : command,
      totpoints, dim+1,
      qh->rbox_isinteger ? "integer" : "real");
  else if (NOcommand)
      qh_fprintf_rbox(qh, qh->fout, 9392, "%d\n%d\n", dim, totpoints);
  else
      /* qh_fprintf_rbox special cases 9393 to append 'command' to the RboxPoints.comment() */
      qh_fprintf_rbox(qh, qh->fout, 9393, "%d %s\n%d\n", dim, command, totpoints);

  /* ============= explicit points =============== */
  if ((s= first_point)) {
    while (s && *s) { /* 'P' */
      count= 0;
      if (iscdd)
        qh_out1(qh, 1.0);
      while (*++s) {
        qh_out1(qh, qh_strtod(s, &s));
        count++;
        if (isspace(*s) || !*s)
          break;
        if (*s != ',') {
          qh_fprintf_rbox(qh, qh->ferr, 6194, "rbox error: missing comma after coordinate in %s\n\n", s);
          qh_errexit_rbox(qh, qh_ERRinput);
        }
      }
      if (count < dim) {
        for (k=dim-count; k--; )
          qh_out1(qh, 0.0);
      }else if (count > dim) {
        qh_fprintf_rbox(qh, qh->ferr, 6195, "rbox error: %d coordinates instead of %d coordinates in %s\n\n",
                  count, dim, s);
        qh_errexit_rbox(qh, qh_ERRinput);
      }
      qh_fprintf_rbox(qh, qh->fout, 9394, "\n");
      while ((s= strchr(s, 'P'))) {
        if (isspace(s[-1]))
          break;
      }
    }
  }

  /* ============= simplex distribution =============== */
  if (issimplex+issimplex2) {
    if (!(simplex= (double*)qh_malloc( dim * (dim+1) * sizeof(double)))) {
      qh_fprintf_rbox(qh, qh->ferr, 6196, "rbox error: insufficient memory for simplex\n");
      qh_errexit_rbox(qh, qh_ERRmem); /* qh_ERRmem */
    }
    simplexp= simplex;
    if (isregular) {
      for (i=0; i<dim; i++) {
        for (k=0; k<dim; k++)
          *(simplexp++)= i==k ? 1.0 : 0.0;
      }
      for (k=0; k<dim; k++)
        *(simplexp++)= -1.0;
    }else {
      for (i=0; i<dim+1; i++) {
        for (k=0; k<dim; k++) {
          randr= qh_RANDOMint;
          *(simplexp++)= 2.0 * randr/randmax - 1.0;
        }
      }
    }
    if (issimplex2) {
        simplexp= simplex;
      for (i=0; i<dim+1; i++) {
        if (iscdd)
          qh_out1(qh, 1.0);
        for (k=0; k<dim; k++)
          qh_out1(qh, *(simplexp++) * box);
        qh_fprintf_rbox(qh, qh->fout, 9395, "\n");
      }
    }
    for (j=0; j<numpoints; j++) {
      if (iswidth)
        apex= qh_RANDOMint % (dim+1);
      else
        apex= -1;
      for (k=0; k<dim; k++)
        coord[k]= 0.0;
      norm= 0.0;
      for (i=0; i<dim+1; i++) {
        randr= qh_RANDOMint;
        factor= randr/randmax;
        if (i == apex)
          factor *= width;
        norm += factor;
        for (k=0; k<dim; k++) {
          simplexp= simplex + i*dim + k;
          coord[k] += factor * (*simplexp);
        }
      }
      for (k=0; k<dim; k++)
        coord[k] *= box/norm;
      qh_outcoord(qh, iscdd, coord, dim);
      if(coincidentcount++ < coincidenttotal)
        qh_outcoincident(qh, coincidentpoints, coincidentradius, iscdd, coord, dim);
    }
    isregular= 0; /* continue with isbox */
    numpoints= 0;
  }

  /* ============= mesh distribution =============== */
  if (ismesh) {
    nthroot= (int)(pow((double)numpoints, 1.0/dim) + 0.99999);
    for (k=dim; k--; )
      mult[k]= 0;
    for (i=0; i < numpoints; i++) {
      coordp= coord;
      for (k=0; k < dim; k++) {
        if (k == 0)
          *(coordp++)= mult[0] * meshn + mult[1] * (-meshm);
        else if (k == 1)
          *(coordp++)= mult[0] * meshm + mult[1] * meshn;
        else
          *(coordp++)= mult[k] * meshr;
      }
      qh_outcoord(qh, iscdd, coord, dim);
      if(coincidentcount++ < coincidenttotal)
        qh_outcoincident(qh, coincidentpoints, coincidentradius, iscdd, coord, dim);
      for (k=0; k < dim; k++) {
        if (++mult[k] < nthroot)
          break;
        mult[k]= 0;
      }
    }
  }
  /* ============= regular points for 's' =============== */
  else if (isregular && !islens) {
    if (dim != 2 && dim != 3) {
      qh_free(simplex);
      qh_fprintf_rbox(qh, qh->ferr, 6197, "rbox error: regular points can be used only in 2-d and 3-d\n\n");
      qh_errexit_rbox(qh, qh_ERRinput);
    }
    if (!isaxis || radius == 0.0) {
      isaxis= 1;
      radius= 1.0;
    }
    if (dim == 3) {
      if (iscdd)
        qh_out1(qh, 1.0);
      qh_out3n(qh, 0.0, 0.0, -box);
      if (!isgap) {
        if (iscdd)
          qh_out1(qh, 1.0);
        qh_out3n(qh, 0.0, 0.0, box);
      }
    }
    angle= 0.0;
    anglediff= 2.0 * PI/numpoints;
    for (i=0; i < numpoints; i++) {
      angle += anglediff;
      x= radius * cos(angle);
      y= radius * sin(angle);
      if (dim == 2) {
        if (iscdd)
          qh_out1(qh, 1.0);
        qh_out2n(qh, x*box, y*box);
      }else {
        norm= sqrt(1.0 + x*x + y*y);
        if (iscdd)
          qh_out1(qh, 1.0);
        qh_out3n(qh, box*x/norm, box*y/norm, box/norm);
        if (isgap) {
          x *= 1-gap;
          y *= 1-gap;
          norm= sqrt(1.0 + x*x + y*y);
          if (iscdd)
            qh_out1(qh, 1.0);
          qh_out3n(qh, box*x/norm, box*y/norm, box/norm);
        }
      }
    }
  }
  /* ============= regular points for 'r Ln D2' =============== */
  else if (isregular && islens && dim == 2) {
    double cos_0;

    angle= lensangle;
    anglediff= 2 * lensangle/(numpoints - 1);
    cos_0= cos(lensangle);
    for (i=0; i < numpoints; i++, angle -= anglediff) {
      x= radius * sin(angle);
      y= radius * (cos(angle) - cos_0);
      if (iscdd)
        qh_out1(qh, 1.0);
      qh_out2n(qh, x*box, y*box);
      if (i != 0 && i != numpoints - 1) {
        if (iscdd)
          qh_out1(qh, 1.0);
        qh_out2n(qh, x*box, -y*box);
      }
    }
  }
  /* ============= regular points for 'r Ln D3' =============== */
  else if (isregular && islens && dim != 2) {
    if (dim != 3) {
      qh_free(simplex);
      qh_fprintf_rbox(qh, qh->ferr, 6198, "rbox error: regular points can be used only in 2-d and 3-d\n\n");
      qh_errexit_rbox(qh, qh_ERRinput);
    }
    angle= 0.0;
    anglediff= 2* PI/numpoints;
    if (!isgap) {
      isgap= 1;
      gap= 0.5;
    }
    offset= sqrt(radius * radius - (1-gap)*(1-gap)) - lensbase;
    for (i=0; i < numpoints; i++, angle += anglediff) {
      x= cos(angle);
      y= sin(angle);
      if (iscdd)
        qh_out1(qh, 1.0);
      qh_out3n(qh, box*x, box*y, 0.0);
      x *= 1-gap;
      y *= 1-gap;
      if (iscdd)
        qh_out1(qh, 1.0);
      qh_out3n(qh, box*x, box*y, box * offset);
      if (iscdd)
        qh_out1(qh, 1.0);
      qh_out3n(qh, box*x, box*y, -box * offset);
    }
  }
  /* ============= apex of 'Zn' distribution + gendim =============== */
  else {
    if (isaxis) {
      gendim= dim-1;
      if (iscdd)
        qh_out1(qh, 1.0);
      for (j=0; j < gendim; j++)
        qh_out1(qh, 0.0);
      qh_out1(qh, -box);
      qh_fprintf_rbox(qh, qh->fout, 9398, "\n");
    }else if (islens)
      gendim= dim-1;
    else
      gendim= dim;
    /* ============= generate random point in unit cube =============== */
    for (i=0; i < numpoints; i++) {
      norm= 0.0;
      for (j=0; j < gendim; j++) {
        randr= qh_RANDOMint;
        coord[j]= 2.0 * randr/randmax - 1.0;
        norm += coord[j] * coord[j];
      }
      norm= sqrt(norm);
      /* ============= dim-1 point of 'Zn' distribution ========== */
      if (isaxis) {
        if (!isgap) {
          isgap= 1;
          gap= 1.0;
        }
        randr= qh_RANDOMint;
        rangap= 1.0 - gap * randr/randmax;
        factor= radius * rangap / norm;
        for (j=0; j<gendim; j++)
          coord[j]= factor * coord[j];
      /* ============= dim-1 point of 'Ln s' distribution =========== */
      }else if (islens && issphere) {
        if (!isgap) {
          isgap= 1;
          gap= 1.0;
        }
        randr= qh_RANDOMint;
        rangap= 1.0 - gap * randr/randmax;
        factor= rangap / norm;
        for (j=0; j<gendim; j++)
          coord[j]= factor * coord[j];
      /* ============= dim-1 point of 'Ln' distribution ========== */
      }else if (islens && !issphere) {
        if (!isgap) {
          isgap= 1;
          gap= 1.0;
        }
        j= qh_RANDOMint % gendim;
        if (coord[j] < 0)
          coord[j]= -1.0 - coord[j] * gap;
        else
          coord[j]= 1.0 - coord[j] * gap;
      /* ============= point of 'l' distribution =============== */
      }else if (isspiral) {
        if (dim != 3) {
          qh_free(simplex);
          qh_fprintf_rbox(qh, qh->ferr, 6199, "rbox error: spiral distribution is available only in 3d\n\n");
          qh_errexit_rbox(qh, qh_ERRinput);
        }
        coord[0]= cos(2*PI*i/(numpoints - 1));
        coord[1]= sin(2*PI*i/(numpoints - 1));
        coord[2]= 2.0*(double)i/(double)(numpoints-1) - 1.0;
      /* ============= point of 's' distribution =============== */
      }else if (issphere) {
        factor= 1.0/norm;
        if (iswidth) {
          randr= qh_RANDOMint;
          factor *= 1.0 - width * randr/randmax;
        }
        for (j=0; j<dim; j++)
          coord[j]= factor * coord[j];
      }
      /* ============= project 'Zn s' point in to sphere =============== */
      if (isaxis && issphere) {
        coord[dim-1]= 1.0;
        norm= 1.0;
        for (j=0; j<gendim; j++)
          norm += coord[j] * coord[j];
        norm= sqrt(norm);
        for (j=0; j<dim; j++)
          coord[j]= coord[j] / norm;
        if (iswidth) {
          randr= qh_RANDOMint;
          coord[dim-1] *= 1 - width * randr/randmax;
        }
      /* ============= project 'Zn' point onto cube =============== */
      }else if (isaxis && !issphere) {  /* not very interesting */
        randr= qh_RANDOMint;
        coord[dim-1]= 2.0 * randr/randmax - 1.0;
      /* ============= project 'Ln' point out to sphere =============== */
      }else if (islens) {
        coord[dim-1]= lensbase;
        for (j=0, norm= 0; j<dim; j++)
          norm += coord[j] * coord[j];
        norm= sqrt(norm);
        for (j=0; j<dim; j++)
          coord[j]= coord[j] * radius/ norm;
        coord[dim-1] -= lensbase;
        if (iswidth) {
          randr= qh_RANDOMint;
          coord[dim-1] *= 1 - width * randr/randmax;
        }
        if (qh_RANDOMint > randmax/2)
          coord[dim-1]= -coord[dim-1];
      /* ============= project 'Wn' point toward boundary =============== */
      }else if (iswidth && !issphere) {
        j= qh_RANDOMint % gendim;
        if (coord[j] < 0)
          coord[j]= -1.0 - coord[j] * width;
        else
          coord[j]= 1.0 - coord[j] * width;
      }
      /* ============= scale point to box =============== */
      for (k=0; k<dim; k++)
          coord[k]= coord[k] * box;

      /* ============= write output =============== */
      qh_outcoord(qh, iscdd, coord, dim);
      if(coincidentcount++ < coincidenttotal)
        qh_outcoincident(qh, coincidentpoints, coincidentradius, iscdd, coord, dim);
    }
  }

  /* ============= write cube vertices =============== */
  if (addcube) {
    for (j=0; j<cubesize; j++) {
      if (iscdd)
        qh_out1(qh, 1.0);
      for (k=dim-1; k>=0; k--) {
        if (j & ( 1 << k))
          qh_out1(qh, cube);
        else
          qh_out1(qh, -cube);
      }
      qh_fprintf_rbox(qh, qh->fout, 9400, "\n");
    }
  }

  /* ============= write diamond vertices =============== */
  if (adddiamond) {
    for (j=0; j<diamondsize; j++) {
      if (iscdd)
        qh_out1(qh, 1.0);
      for (k=dim-1; k>=0; k--) {
        if (j/2 != k)
          qh_out1(qh, 0.0);
        else if (j & 0x1)
          qh_out1(qh, diamond);
        else
          qh_out1(qh, -diamond);
      }
      qh_fprintf_rbox(qh, qh->fout, 9401, "\n");
    }
  }

  if (iscdd)
    qh_fprintf_rbox(qh, qh->fout, 9402, "end\nhull\n");

  /* same code for error exit and normal return */
  qh_free(simplex);
  return qh_ERRnone;
} /* rboxpoints */

/*------------------------------------------------
outxxx - output functions for qh_rboxpoints
*/
int qh_roundi(qhT *qh, double a) {
  if (a < 0.0) {
    if (a - 0.5 < INT_MIN) {
      qh_fprintf_rbox(qh, qh->ferr, 6200, "rbox input error: negative coordinate %2.2g is too large.  Reduce 'Bn'\n", a);
      qh_errexit_rbox(qh, qh_ERRinput);
    }
    return (int)(a - 0.5);
  }else {
    if (a + 0.5 > INT_MAX) {
      qh_fprintf_rbox(qh, qh->ferr, 6201, "rbox input error: coordinate %2.2g is too large.  Reduce 'Bn'\n", a);
      qh_errexit_rbox(qh, qh_ERRinput);
    }
    return (int)(a + 0.5);
  }
} /* qh_roundi */

void qh_out1(qhT *qh, double a) {

  if (qh->rbox_isinteger)
    qh_fprintf_rbox(qh, qh->fout, 9403, "%d ", qh_roundi(qh, a+qh->rbox_out_offset));
  else
    qh_fprintf_rbox(qh, qh->fout, 9404, qh_REAL_1, a+qh->rbox_out_offset);
} /* qh_out1 */

void qh_out2n(qhT *qh, double a, double b) {

  if (qh->rbox_isinteger)
    qh_fprintf_rbox(qh, qh->fout, 9405, "%d %d\n", qh_roundi(qh, a+qh->rbox_out_offset), qh_roundi(qh, b+qh->rbox_out_offset));
  else
    qh_fprintf_rbox(qh, qh->fout, 9406, qh_REAL_2n, a+qh->rbox_out_offset, b+qh->rbox_out_offset);
} /* qh_out2n */

void qh_out3n(qhT *qh, double a, double b, double c) {

  if (qh->rbox_isinteger)
    qh_fprintf_rbox(qh, qh->fout, 9407, "%d %d %d\n", qh_roundi(qh, a+qh->rbox_out_offset), qh_roundi(qh, b+qh->rbox_out_offset), qh_roundi(qh, c+qh->rbox_out_offset));
  else
    qh_fprintf_rbox(qh, qh->fout, 9408, qh_REAL_3n, a+qh->rbox_out_offset, b+qh->rbox_out_offset, c+qh->rbox_out_offset);
} /* qh_out3n */

void qh_outcoord(qhT *qh, int iscdd, double *coord, int dim) {
    double *p= coord;
    int k;

    if (iscdd)
      qh_out1(qh, 1.0);
    for (k=0; k < dim; k++)
      qh_out1(qh, *(p++));
    qh_fprintf_rbox(qh, qh->fout, 9396, "\n");
} /* qh_outcoord */

void qh_outcoincident(qhT *qh, int coincidentpoints, double radius, int iscdd, double *coord, int dim) {
  double *p;
  double randr, delta;
  int i,k;
  double randmax= qh_RANDOMmax;

  for (i= 0; i<coincidentpoints; i++) {
    p= coord;
    if (iscdd)
      qh_out1(qh, 1.0);
    for (k=0; k < dim; k++) {
      randr= qh_RANDOMint;
      delta= 2.0 * randr/randmax - 1.0; /* -1..+1 */
      delta *= radius;
      qh_out1(qh, *(p++) + delta);
    }
    qh_fprintf_rbox(qh, qh->fout, 9410, "\n");
  }
} /* qh_outcoincident */

/*------------------------------------------------
   Only called from qh_rboxpoints or qh_fprintf_rbox
   qh_fprintf_rbox is only called from qh_rboxpoints
*/
void qh_errexit_rbox(qhT *qh, int exitcode)
{
    longjmp(qh->rbox_errexit, exitcode);
} /* qh_errexit_rbox */

