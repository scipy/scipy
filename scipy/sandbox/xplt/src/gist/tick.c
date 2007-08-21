/*
 * TICK.C
 *
 * $Id$
 *
 * Implement ticks and labels for GIST coordinate systems
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

/* NOTE-- probably still a few bugs, especially with log scales.
   I wish the log scale coding were cleaner.  The only bug I've
   noticed is fairly minor:  With limits 45 to 280 (e.g.) on
   a log scale, the (unlabeled) tick at 45 is missing.  Since
   this is the bottommost tick, the scale still doesn't look
   bad.
 */

#include "gist.h"

/* Use sprintf function heavily */
#include <stdio.h>

/* The nice unit routine defined here is also used in draw.c */
extern GpReal GpNiceUnit(GpReal finest, int *base, int *power);

#define TICK_ANY (TICK_L | TICK_U | TICK_C)
#define LABEL_ANY (LABEL_L | LABEL_U)

/* WC->NDC mapping stored by GaTicks, used by Draw?Ticks
   and Draw?Labels routines */
static GpReal scalx, offx, scaly, offy;

/* The main output from FindTicks is a cumulative count of the
   number of ticks at each level, nChangeLevel[0:nLevel]
   (nLevel<TICK_LEVELS), and the floating point tick values
   themselves, ticks[0:nTotal-1] */
static long nTotal;
static int nLevel;
static int nChangeLevel[TICK_LEVELS];
static GpReal *ticks;

/* Scratch space maintained by GaGetScratchP */
extern int GaGetScratchP(long n);
extern GpReal *gaxScratch, *gayScratch;

/* Additional stuff required for sub-decade tick labels */
static GpReal subTick0[3], subUnit[3];
static int subPower[3], nSubtotal[3];
static GpReal minSubdecade;

/* FindTicks also assures that there is enough working space for
   the scratch arrays (px,py), (qx,qy) used by Draw?Ticks */
GpReal *px, *py, *qx, *qy;

static GaAltTicks *altticks= 0;
static GaAltLabel *altlabel= 0;

/* use fancy exponents if escapes recognized */
extern int gtDoEscapes;

/* ------------------------------------------------------------------------ */

extern char *strcpy(char *, const char *);

extern double floor(double);
extern double ceil(double);
extern double fabs(double);
extern double log10(double);

#ifndef NO_EXP10
  extern double exp10(double);
#else
# define exp10(x) pow(10.,x)
  extern double pow(double,double);
#endif

#define LOG2 0.301029996
#define LOG5 0.698970004
#define LOG10o9 0.0457574905
#define LOG6o5 0.0791812460
#define LOG5o4 0.0969100130

static void EvenlySpace(GpReal start, GpReal unit, GpReal stop);
static GpReal Subdivide(GpReal unit, int nDivs, GpReal lo, GpReal hi,
                        GpReal *start);
static void FirstSubDecade(GpReal *unit, GpReal lo, GpReal hi,
                           GpReal zSmallest);
static int LogDivide(GpReal *unit, int *base, GpReal *fine,
                     GpReal lo, GpReal hi, GpReal zSmallest);
static void ScopeOutTicks(GpReal *loin, GpReal *hiin,
                          GaAxisStyle *style, int isLog,
                          int *useLog, GpReal *nMajor, GpReal *nMinor,
                          GpReal *jUnit, int *jBase, int *jPower,
                          GpReal *itick0);
static int FindOmitted(GpReal lo, GpReal hi);
static GpReal FindOrigin(GpReal lo, GpReal hi, GaAxisStyle *style, int isLog);
static void FindTicks(GpReal lo, GpReal hi, GaAxisStyle *style, int isLog);
static void TickEndpoints(GpReal place, GpReal len, int lo, int hi,
                          GpReal *p0, GpReal *p1);
static void DrawXTicks(GpReal place, GpReal *lens, int inout, int upper,
                       GpLineAttribs *attribs);
static void DrawYTicks(GpReal place, GpReal *lens, int inout, int upper,
                       GpLineAttribs *attribs);
static void DrawOrigin(GpReal x0, GpReal x1, GpReal y0, GpReal y1,
                       GpLineAttribs *attribs);
static void InitLabels(int nDigits, char exponent[]);
static void NextLabel(char label[]);
static void NiceLogLabel(char *label, char *format, GpReal value, GpReal tick);
static void GtText(GpReal x0, GpReal y0, const char *text);
static int DrawXLabels(int isLog, GpReal place, int upper, int nDigits,
                       GpTextAttribs *attribs);
static int DrawYLabels(int isLog, GpReal place, int upper, int nDigits,
                       GpTextAttribs *attribs);
static void DrawOverflow(GpReal x, GpReal y);

/* ------------------------------------------------------------------------ */

GpReal GpNiceUnit(GpReal finest, int *base, int *power)
{
  int b, p;
  GpReal unit;
  if (finest==0.0) finest= 1.0e-6;
  p= (int)floor(log10(fabs(finest)));
  unit= exp10((double)p);
  finest= finest/unit;
  if (finest>5.0) {
    p++;
    unit*= 10.0;
    b= 1;
  } else if (finest>2.0) b= 5;
  else if (finest>1.0) b= 2;
  else b= 1;
  *base= b;
  *power= p;
  return b*unit;
}

static void EvenlySpace(GpReal start, GpReal unit, GpReal stop)
{
  if (start<=stop) {
    GpReal v, vo;  /* avoid infinite loop if stop-start and unit  << start */
    ticks[nTotal++]= vo= v= start;
    for (v+=unit ; v<=stop && v>vo ; v+=unit) ticks[nTotal++]= vo= v;
  }
}

static GpReal Subdivide(GpReal unit, int nDivs, GpReal lo, GpReal hi,
                        GpReal *start)
{
  GpReal newUnit= unit/(GpReal)nDivs;
  GpReal value= *start;
  GpReal test= lo+unit;
  int j;
  for (j=1 ; j<nDivs ; j++) {
    value+= newUnit;
    if (value>=test) {
      value-= unit;
      *start= value;
    }
    EvenlySpace(value, unit, hi);
  }
  return newUnit;
}

static GpReal zi[3]= {1.0, 2.0, 5.0};
static GpReal zf[3]= {2.0, 5.0, 10.0};
static GpReal zlogi[3]= {0.0, LOG2, LOG5};
static GpReal zlogf[3]= {LOG2, LOG5, 1.0};

#define EPS_LOG 0.00001

static GpReal zfSigh[3]= {2.0+EPS_LOG, 5.0+EPS_LOG, 10.0-EPS_LOG};

static void FirstSubDecade(GpReal *unit, GpReal lo, GpReal hi,
                           GpReal zSmallest)
{
  int ii, iiset= -1;
  GpReal z, zlog, oldUnit;
  GpReal decade= ceil(lo);
  GpReal zDecTest= lo+1.0;
  GpReal interval0or1;

  if (hi <= decade) decade-= 1.0;

  /* Outer loop is on interval (1,2], (2,5], or (5,10). */
  for (ii=0 ; ii<3 ; ii++) {
    oldUnit= unit[ii];
    if (ii>0) nSubtotal[ii-1]= nTotal;
    interval0or1= ii<2? EPS_LOG : 0.0;

  retry:
    if (zlogf[ii]+interval0or1 < lo-decade) continue;
    zlog= zlogi[ii]+decade;
    if (zlog+interval0or1 < lo) {
      /* Some ticks are present in no decades at all (this can only happen
         if hi-lo<1.0).  Try to avoid computing log10s by jumping directly
         to 1st tick.  */
      z= oldUnit*ceil(zSmallest/oldUnit);
    } else {
      if (zlog > hi) {
        /* This can only happen if hi-lo < 1.0. */
        if (decade < lo) goto done;
        decade-= 1.0;
        goto retry;
      }
      /* All ticks are present, but must deal with case that oldUnit
         is larger than interval without duplicating ticks... */
      z= oldUnit*ceil(zi[ii]/oldUnit);
    }

    /* Lower boundary of all three intervals must be skipped */
    if (ii>iiset) {
      /* never want to do this after retry if some ticks placed */
      if (z<=zi[ii]) z+= oldUnit;
      subTick0[ii]= z;
      subUnit[ii]= oldUnit;
    }

    /* Next loop is on unit divisions. */
    for ( ; z<zfSigh[ii]+interval0or1 ; z+=oldUnit) {
      zlog= log10(z)+decade;
      if (zlog >= zDecTest) { /* zDecTest= lo+1.0 */
        zlog-= 1.0;
        decade-= 1.0;
      }
      if (zlog > hi) {
        /* This can only happen if hi-lo < 1.0. */
        if (decade < lo) goto done;
        decade-= 1.0;
        goto retry;
      }
      /* Inner loop is on decades. */
      EvenlySpace(zlog, 1.0, hi);
      iiset= ii;
      if (zlog < minSubdecade) minSubdecade= zlog;
    }
  }
  ii= 2;
 done:
  for ( ; ii<3 ; ii++) nSubtotal[ii]= nTotal;
  minSubdecade= floor(minSubdecade);
}

static int LogDivide(GpReal *unit, int *base, GpReal *fine,
                     GpReal lo, GpReal hi, GpReal zSmallest)
{
  int i, ii, nDivs, newBase;
  GpReal z, zu, zlog, ztop, newUnit, oldUnit;
  GpReal du, decade= ceil(lo);
  GpReal zDecTest= lo+1.0;
  GpReal zTest;
  int anyDivisions= 0;

  if (hi <= decade) decade-= 1.0;

  /* Outer loop is on interval (1,2], (2,5], or (5,10). */
  for (ii=0 ; ii<3 ; ii++) {
    oldUnit= unit[ii];
    newBase= base[ii];
    if (newBase==1) {
      /* Use halves if it can be subdivided again, else try fifths. */
      if (oldUnit>=10.0*fine[ii] || oldUnit<5.0*fine[ii]) nDivs= 2;
      else nDivs= 5;
      newBase= 7-nDivs;
    } else {
      nDivs= newBase;
      newBase= 1;
    }
    newUnit= oldUnit/nDivs;
    if (newUnit<fine[ii]) continue;
    unit[ii]= newUnit;
    base[ii]= newBase;
    anyDivisions= 1;

  retry:
    if (zlogf[ii]+decade < lo) continue;
    zlog= zlogi[ii]+decade;
    if (zlog < lo) {
      /* Some ticks are present in no decades at all (this can only happen
         if hi-lo<1.0).  Try to avoid computing log10s by jumping directly
         to 1st tick.  */
      z= oldUnit*ceil(zSmallest/oldUnit);
      zTest= zSmallest+oldUnit;
    } else {
      if (zlog > hi) {
        if (decade < lo) goto done;
        decade-= 1.0;
        goto retry;
      }
      /* All ticks are present, but must deal with case that oldUnit
         is larger than interval without duplicating ticks... */
      /* z= oldUnit*ceil(zi[ii]/oldUnit); */
      z= zi[ii];
      zTest= zf[ii];
    }
    ztop= zf[ii]-EPS_LOG;  /* normally, don't tick top of interval */

    /* Obnoxius special cases... */
    if (newBase==1 && newUnit>0.5) { /* i.e.- newUnit==1.0 */
      if (ii==0) continue;
      if (ii==1) {
        if (oldUnit<3.0) ztop= 5.5;     /* i.e.- oldUnit==2.0, do tick at 5 */
      } else {
        if (oldUnit<3.0 && z<6.0) z= 6.0; /* i.e.- oldUnit==2.0, begin at 7 */
      }
    }

    /* Next loop is on subdivisions. */
    du= decade;
    for (i=1 ; i<nDivs ; i++) {
      z+= newUnit;
      if (z > zTest) { /* this can only happen on zSmallest branch above */
        z-= oldUnit;
        if (z <= zi[ii]+EPS_LOG) break; /* allows for no ticks */
      }
      decade= du;
      /* Next loop is on unit divisions. */
      for (zu=z ; zu<ztop ; zu+=oldUnit) {
        zlog= log10(zu)+decade;
        if (zlog >= zDecTest) { /* zDecTest= lo+1.0 */
          zlog-= 1.0;
          decade-= 1.0;
        }
        if (zlog > hi) {
          /* This can only happen if hi-lo < 1.0. */
          if (zu!=z) break;  /* log(z) may still be less than hi */
          if (decade < lo) goto done;
          decade-= 1.0;
          goto retry;
        }
        /* Inner loop is on decades. */
        EvenlySpace(zlog, 1.0, hi);
      }
    }

  }
 done:

  return anyDivisions;
}

/* ------------------------------------------------------------------------ */

/* Things set by ScopeOutTicks() utility routine */
static int useLog;        /* set if log ticks must be generated */
static GpReal nMajor;     /* extracted from GaAxis style, adjusted if... */
static GpReal nMinor;     /* ...this is log scale reverted to linear */
static GpReal jUnit;      /* nice number just above delta/nMajor */
static int jBase, jPower; /* jUnit= jBase*10^jPower, jBase= 1, 2, or 5 */
static GpReal itick0;     /* integer "number of major unit" for 1st tick */

/* After FindTicks, first tick is at tick0, and first major tick
   has index itick0 (tick is at itick0*jUnit) */
static GpReal tick0;
static int subDecadeTicks; /* set if there are sub-decade major ticks */

static void ScopeOutTicks(GpReal *loin, GpReal *hiin,
                          GaAxisStyle *style, int isLog,
                          int *useLog, GpReal *nMajor, GpReal *nMinor,
                          GpReal *jUnit, int *jBase, int *jPower,
                          GpReal *itick0)
{
  GpReal lo= *loin;
  GpReal hi= *hiin;
  GpReal finest, delta;

  *nMajor= style->nMajor;
  *nMinor= style->nMinor;
  if (*nMinor<*nMajor) *nMinor= *nMajor;

  if (hi<lo) { GpReal tmp= lo; lo= hi; hi= tmp; }
  delta= hi-lo;
  if (delta==0.0) {
    if (lo>0.0) delta= lo;
    else if (lo<0.0) delta= -lo;
    else delta= 0.01;
  }

  /* Adjust lo and hi slightly to avoid rounding errors that can cause
     ticks and or labels falling exactly at the endpoints of the
     interval to be missed.  This is a little heavy handed, but it
     does the desirable thing in most cases.  */
  lo-= 0.0001*delta;
  hi+= 0.0001*delta;
  delta= hi-lo;

  if (isLog) {
    *useLog= delta > LOG2;
    if (!*useLog) {
      /* Revert to linear scale if full scale less than factor of 2... */
      lo= exp10(lo);
      hi= exp10(hi);
      /* ...but adjust number of divisions so that finest linear division
         is no finer than finest log division would have been */
      *nMajor= (1.-(lo/hi))/(1.-exp10(-delta/(*nMajor)));
      *nMinor= (1.-(lo/hi))/(1.-exp10(-delta/(*nMinor)));
      delta= hi-lo;
    }
  } else {
    *useLog= 0;
  }

  finest= delta/(*nMajor);
  if (!(*useLog) || finest>1.0) *jUnit= GpNiceUnit(finest, jBase, jPower);
  else { *jUnit= 1.0;   *jBase= 1;   *jPower= 0; }

  *itick0= ceil(lo/(*jUnit));

  *loin= lo;
  *hiin= hi;
}

static int FindOmitted(GpReal lo, GpReal hi)
{
  GpReal z= itick0*jUnit;
  int nTicks= nChangeLevel[0];
  if (nTicks>0) {
    GpReal multiple= 100./(GpReal)jBase; /* multiple of jUnit ends in "00" */
    multiple*= ceil(itick0/multiple);    /* itick0 rounded up to multiple */
    nTicks--;
    if (multiple <= itick0+(GpReal)nTicks) return (int)(multiple-itick0);
    else return 0;
  } else {
    /* This is a log scale with no decades (sigh).
       Know that jUnit=1.0, but hi-lo>LOG2, since it hasn't
       reverted to a linear scale.  Therefore, either 2, 5, or 4
       must be in the interval-- take the first one that is.
       Initially, z=itick0*jUnit is the decade above the interval.  */
    int i;
    if (z+(LOG2-1.0) >= lo) z+= LOG2-1.0;         /* 2 */
    else if (z-LOG2 <= hi) z-= LOG2;              /* 5 */
    else z+= 2.0*LOG2-1.0;                        /* 4 */
    /* Since there are no decade ticks, the ticks are in strictly
       increasing order.  */
    for (i=0 ; i<nChangeLevel[1] && ticks[i]+EPS_LOG<z ; i++);
    return i;
  }
}

static GpReal FindOrigin(GpReal lo, GpReal hi, GaAxisStyle *style, int isLog)
{
  int useLog;
  GpReal nMajor, nMinor, jUnit;
  int jBase, jPower;
  GpReal itick0;

  GpReal z, multiple, origin;
  int nTicks;

  ScopeOutTicks(&lo, &hi, style, isLog, &useLog, &nMajor, &nMinor,
                &jUnit, &jBase, &jPower, &itick0);

  /* Count number of ticks using same algorithm as EvenlySpace */
  origin= z= itick0*jUnit;
  for (nTicks=0 ; z<=hi ; z+=jUnit) nTicks++;

  if (nTicks>0) {
    multiple= 100./(GpReal)jBase; /* multiple of jUnit that ends in "00" */
    multiple*= ceil(itick0/multiple);   /* itick0 rounded up to multiple */
    nTicks--;
    if (multiple <= itick0+(GpReal)nTicks) origin= multiple*jUnit;
    /* else origin remains itick0*jUnit */
  } else {
    /* This is a log scale with no decades (sigh).
       Know that jUnit=1.0, but hi-lo>LOG2, since it hasn't
       reverted to a linear scale.  Therefore, either 2, 5, or 4
       must be in the interval-- take the first one that is.
       Initially, z=itick0*jUnit is the decade above the interval.  */
    if (z+(LOG2-1.0) >= lo) origin= z+(LOG2-1.0);       /* 2 */
    else if (z-LOG2 <= hi) origin= z-LOG2;              /* 5 */
    else origin= z+(2.0*LOG2-1.0);                      /* 4 */
  }

  if (isLog && !useLog) origin= log10(origin);

  return origin;
}

static void FindTicks(GpReal lo, GpReal hi, GaAxisStyle *style, int isLog)
{
  GpReal finest, finest2, finest5, finest10;
  GpReal nUnit = 0;
  int /*nBase,*/ nPower;
  int reqdSpace;

  /* Swap limits if necessary to force lo<hi, and exponentiate if this
       scale will revert from log to linear (set useLog appropriately),
     extract nMajor, nMinor from style,
       possibly adjusting for reverted log scales,
     compute the nice major unit (jUnit, jBase, jPower),
       and itick0, the integer "unit number" of the first major tick  */
  ScopeOutTicks(&lo, &hi, style, isLog, &useLog, &nMajor, &nMinor,
                &jUnit, &jBase, &jPower, &itick0);

  /* Be sure there is enough scratch space to hold the ticks */
  reqdSpace= (int)(useLog? style->logAdjMinor*nMinor : nMinor) + 16;
  if (GaGetScratchP(3*reqdSpace)) {
    nTotal= 0;
    nLevel= 0;
    nChangeLevel[0]= 0;
    return; /* memory manager failure */
  } else {
    px= gaxScratch;
    qx= px+reqdSpace;
    py= gayScratch;
    qy= py+reqdSpace;
    ticks= qx+reqdSpace;
  }

  if (altticks && !isLog) {
    int i;
    for (i=0 ; i<TICK_LEVELS ; i++) nChangeLevel[i]= 0;
    if (!altticks(lo, hi, nMajor, nMinor, ticks, nChangeLevel)) {
      nTotal= nChangeLevel[0];
      for (nLevel=TICK_LEVELS-1 ; nLevel>0 ; nLevel--)
        if ((nTotal= nChangeLevel[nLevel])) break;
      /* attempt to set jUnit, jPower, jBase */
      if (nChangeLevel[0]>0) {
        if (nChangeLevel[0]>1) {
          GpReal pwr;
          jUnit= fabs(ticks[1]-ticks[0]);
          jPower= (int)floor(log10(jUnit)+0.000001);
          pwr= exp10((double)jPower);
          for (i=4 ; i ; i--) {
            jBase= (int)(jUnit/pwr+0.000001);
            if (fabs(jBase-jUnit/pwr)<0.0001) break;
            jPower--;
            pwr/= 10;
          }
        }
        /* note: itick0*jUnit must equal 1st tick value
                 *and* other ticks must be at itick0+n*jBase */
        itick0= ticks[0]/jUnit;
      }
      tick0= itick0*jUnit;
      goto done;
    }
    altticks= 0;  /* can be used as a flag later if necessary */
  }

  /* Descend hierarchy from major divisions to minor divisions.
     The minor divisions are nUnit= nBase*10^nPower.
     nLevel+1= the number of levels (major is level 0,
                                     minor is level nLevel) */
  nPower= jPower;
  nUnit= jUnit;
  finest= (hi-lo)/nMinor;
  finest2= finest*2.0;
  finest5= finest*5.0;
  finest10= finest*10.0;

  tick0= itick0*jUnit;
  nLevel= 0;
  nTotal= 0;
  EvenlySpace(tick0, jUnit, hi);
  nChangeLevel[0]= nTotal;

  if (jBase==2) {
    if (nUnit < finest2) goto done;
    nUnit= Subdivide(nUnit, 2, lo, hi, &tick0);
    nChangeLevel[++nLevel]= nTotal;
  } else if (jBase==5) goto mid;

  for (;;) {  /* Loop to finer and finer decades */
    /* nUnit is 1.0eNN */
    if (nUnit<finest2 || (useLog&&nPower==0) || nLevel>=TICK_LEVELS-1) {
      /*nBase= 1;*/
      break;
    }
    nPower--;

    if (nUnit>=finest5 && nUnit<finest10) {
      /* Subdivide into 2.0eNN if and only if no further subdivisions */
      nUnit= Subdivide(nUnit, 5, lo, hi, &tick0);
      nChangeLevel[++nLevel]= nTotal;
      /*nBase= 2;*/
      break;
    }
    nUnit= Subdivide(nUnit, 2, lo, hi, &tick0);
    nChangeLevel[++nLevel]= nTotal;

  mid:
    /* nUnit is 5.0eNN */
    if (nUnit<finest5 || nLevel>=TICK_LEVELS-1) {
      /*nBase= 5;*/
      break;
    }
    nUnit= Subdivide(nUnit, 5, lo, hi, &tick0);
    nChangeLevel[++nLevel]= nTotal;
  }
 done:

  if (isLog && !useLog) {
    /* must adjust ticks to produce log scale (lo and hi were exp10'd) */
    int i;
    for (i=0 ; i<nTotal ; i++) ticks[i]= log10(ticks[i]);
  }

  /*------------ Do sub-decade log ticks if necessary -------------*/
  subDecadeTicks= 0;
  if (useLog && nUnit<1.1 && nLevel<TICK_LEVELS-1) {
    GpReal delta= hi-lo;
    GpReal logAdjMajor= style->logAdjMajor;
    GpReal logAdjMinor= style->logAdjMinor;
    GpReal fine[3];  /* 0 is interval (1,2], 1 is (2,5], 2 is (5,10) */
    GpReal unit[3];
    int base[3];
    GpReal zSmallest;

    /* Apply subdecade adjustment factors */
    nMajor*= logAdjMajor;
    nMinor*= logAdjMinor;

    /* Compute separate finest divisions for each of the 3 parts of
       a decade-- (1,2], (2,5], and (5,10).  The formulas are the same
       as for the revert-to-linear adjustment if isLog && !useLog. */
    finest= (1.-exp10(-delta/nMinor));
    fine[0]= 2.0*finest;
    fine[1]= 5.0*finest;
    fine[2]= 10.0*finest;

    /* Compute the major units from which to begin subdivision.
       The "natural" progression of subdivision is adjusted to give
       scales which "look logarithmic", i.e.- so that the tick spacing
       at each level varies more than it would if the 3 intervals
       were divided like linear scales. */

    finest= delta/nMajor;
    if (finest>LOG2) {
      finest= delta/nMinor;
      if (finest>LOG2) return;
      subDecadeTicks= 0;  /* no major sub-decade ticks */
    } else {
      subDecadeTicks= 1;  /* set if there are sub-decade major ticks */
    }
    zSmallest= hi-lo<1.0? exp10(lo-floor(lo)) : 1.0;

    if (subDecadeTicks && (delta<LOG5 || finest<=LOG10o9)) {
      /* Choose linear units best suited to each of the 3 intervals.
         These scales tend to "look linear".  If there are to be no
         major ticks, never need to choose this. */
      finest= (1.-exp10(-finest));  /* see note above */
      unit[0]= GpNiceUnit(2.0*finest, base+0, subPower+0);
      unit[1]= GpNiceUnit(5.0*finest, base+1, subPower+1);
      unit[2]= GpNiceUnit(10.0*finest, base+2, subPower+2);
      if (finest>LOG6o5 && base[1]==1) {
        /* log10(5/4) >= finest > log10(6/5) needs correction. */
        base[1]= 2;
        unit[1]= 2.0;
      }
    } else if (finest<=LOG5o4) {
      /* Subdecade ticks at 2, 4, 6, and 8 only. */
      unit[0]= 2; /* skips next level (subdivision into 1s) */
      unit[1]= 2; /* next level is 1s */
      unit[2]= 2; /* next level is 1s */
      base[0]= base[1]= base[2]= 2;
      subPower[0]= subPower[1]= subPower[2]= 0;
    } else {
      /* Subdecade ticks at 2 and 5 only. */
      unit[0]= 2; /* skips next level (subdivision into 1s) */
      unit[1]= 5; /* next level is 1s */
      unit[2]= 5; /* next level is 1s */
      base[0]= 2;
      base[1]= base[2]= 5;
      subPower[0]= subPower[1]= subPower[2]= 0;
    }

    minSubdecade= ceil(lo);
    FirstSubDecade(unit, lo, hi, zSmallest);
    nChangeLevel[++nLevel]= nTotal;

    /* Loop on tick level.  Each of the 3 intervals is subdivided once
       (or not at all) per pass.  */
    while (nLevel<TICK_LEVELS-1 &&
           LogDivide(unit, base, fine, lo, hi, zSmallest))
      nChangeLevel[++nLevel]= nTotal;

    /* Correct for possibility that final pass produced no ticks.  */
    if (nLevel>0 && nChangeLevel[nLevel-1]==nChangeLevel[nLevel]) nLevel--;
  }
}

/* ------------------------------------------------------------------------ */

static void TickEndpoints(GpReal place, GpReal len, int lo, int hi,
                          GpReal *p0, GpReal *p1)
{
  if (lo) {
    if (hi) {   /* symmetric ticks */
      *p0= place-0.5*len;
      *p1= place+0.5*len;
    } else {    /* right or top justified ticks */
      *p0= place-len;
      *p1= place;
    }
  } else {    /* left or bottom justified ticks */
    *p0= place;
    *p1= place+len;
  }
}

static void DrawXTicks(GpReal place, GpReal *lens, int inout, int upper,
                       GpLineAttribs *attribs)
{
  int i, level, llev;
  GpReal y0, y1;
  int lo= upper? (inout & TICK_IN) : (inout & TICK_OUT);
  int hi= upper? (inout & TICK_OUT) : (inout & TICK_IN);

  level= llev= 0;
  TickEndpoints(place, lens[level], lo, hi, &y0, &y1);
  for (i=0 ; i<nTotal ; i++) {
    if (i==nChangeLevel[level]) {
      level++;
      if (i>0) TickEndpoints(place, lens[++llev], lo, hi, &y0, &y1);
    }
    py[i]= y0;
    qy[i]= y1;
    qx[i]= px[i]= scalx*ticks[i]+offx;
  }

  gistA.l= *attribs;
  GpDisjoint(nTotal, px, py, qx, qy);
}

static void DrawYTicks(GpReal place, GpReal *lens, int inout, int upper,
                       GpLineAttribs *attribs)
{
  int i, level, llev;
  GpReal x0, x1;
  int lo= upper? (inout & TICK_IN) : (inout & TICK_OUT);
  int hi= upper? (inout & TICK_OUT) : (inout & TICK_IN);

  level= llev= 0;
  TickEndpoints(place, lens[level], lo, hi, &x0, &x1);
  for (i=0 ; i<nTotal ; i++) {
    if (i==nChangeLevel[level]) {
      level++;
      if (i>0) TickEndpoints(place, lens[++llev], lo, hi, &x0, &x1);
    }
    px[i]= x0;
    qx[i]= x1;
    qy[i]= py[i]= scaly*ticks[i]+offy;
  }

  gistA.l= *attribs;
  GpDisjoint(nTotal, px, py, qx, qy);
}

static void DrawOrigin(GpReal x0, GpReal x1, GpReal y0, GpReal y1,
                       GpLineAttribs *attribs)
{
  gistA.l= *attribs;
  GpDisjoint(1L, &x0, &y0, &x1, &y1);
}

/* ------------------------------------------------------------------------ */

static char overflow[32];
static char fixedFormat[16];   /* e.g.- " %+05.0f" */
static int decimalPoint;
static int niceDecades;
static char *niceDecs[]=
       { "0.001", "0.01", "0.1", "1.0", "10.0", "100.0", "1000." };
static int overflowChar;
static char integerFormat[]= "% .0f";
static char decadeFormat[]= "E%+.0f";
static char nicedFormat[]= "10^%+.0f";

static GpReal iValue;  /* integer value of current label */
static int omitX, omitY;  /* indices of origin label for TICK_C */

static void InitLabels(int nDigits, char exponent[])
{
  int nLabels= nChangeLevel[0];
  GpReal maxAbs;
  int iPower, xPower;

  if (nLabels<=0) {
    int subDigits= (int)minSubdecade + subPower[0];
    /* log axis with no decade labels only (know jBase=1) */
    if (subDigits>=-3 && itick0<=3.0) {
      exponent[0]= '\0';
      niceDecades= 1;
    } else {
      sprintf(exponent, "E%+.0f", itick0-1.0);
      niceDecades= 0;
    }
    return;
  }

  /* iValue is the INTEGER (not normalized) value of the first label
     -- the first tick is at iValue*10^jPower */
  iValue= itick0*((GpReal)jBase);

  /* iPower is the decade of the largest integer in the list of tick values
     -- the decade of the largest label is iPower+jPower */
  maxAbs= fabs(itick0+(GpReal)(nLabels-1));
  if (nLabels>1 && maxAbs<fabs(itick0)) maxAbs= fabs(itick0);
  iPower= (int)log10(maxAbs*((GpReal)jBase)+0.5);
  xPower= iPower+jPower;

  if (nDigits<2) nDigits= 2;
  decimalPoint= 0;
  exponent[0]= '\0';
  niceDecades= 0;
  overflow[0]= '\0';
  overflowChar= 0;

  if (useLog) {
    /* Labels are decades on a log scale */
    int subDigits= (int)minSubdecade + subPower[0];
    int jp= jPower;
    GpReal fValue= iValue + (GpReal)((nLabels-1)*jBase);
    while (jp--) { iValue*= 10.0;  jBase*= 10; }
    if (subDecadeTicks) niceDecades= (subDigits>=-3 && fValue<=3.0);
    else niceDecades= (iValue>=-3.0 && fValue<=3.0);

  } else if ((xPower>3 || jPower<-3 || iPower>3) && iPower<nDigits) {
    /* Labels which are 10,000 and up, or which involve units less
       than 0.001, or which require more than 4 digits, will be printed
       in scientific notation */
    sprintf(exponent, "E%+d", xPower);
    sprintf(fixedFormat, " %%+0%d.0f", 2+iPower);
    decimalPoint= 2;  /* sign and leading digit move left for decimal pt */

  } else if (jPower>=0 && (iPower<nDigits || xPower<3)) {
    /* Labels which are multiples of 1.0 up to 9999 will be
       printed as integers */
    /* Make iValue, jBase the actual integer values */
    int value= (int)iValue;
    while (jPower--) {                /* WARNING- clobbers jBase, jPower */
      value*= 10;
      jBase*= 10;
    }
    iValue= (GpReal)value;

  } else if (iPower<nDigits || (jPower>-3 && iPower<3)) {
    /* Labels with units less than 1.0 but greater than 0.001, which
       can be written with 4 digits or fewer will be printed in
       fixed point notation */
    int iDigits= xPower<0? -jPower : iPower;
    sprintf(fixedFormat, " %%+0%d.0f", 2+iDigits);
    decimalPoint= 2+iDigits+jPower;

  } else {
    /* These labels are too wide to fit within the nDigit field */
    GpReal origin= 100.0*ceil(iValue/100.0); /* round iValue up to 100s */
    sprintf(fixedFormat, "%%+0%d.0f", 2+iPower);
    overflowChar= iPower;   /* first iPower-1 digits are overflow, last
                               two are used on individual labels */
    if (origin-iValue > (GpReal)(jBase*(nLabels-1))) origin= iValue;
    /* Assume origin has same number of digits as maxAbs (iPower) */
    sprintf(overflow, "x0= %+.0fE%+d", origin, xPower);
    /* Want, e.g.-  "x0=-1.2345??E+12"   (iPower=6) */
    if (origin<0.0) overflow[3]= '-';
    overflow[4]= overflow[5];
    overflow[5]= '.';
    overflow[4+iPower]= overflow[5+iPower]= '?';
  }
}

static GpReal positiveZero= 0.0; /* Yecch */

static void NextLabel(char label[])
{
  if (iValue==0.0) iValue= positiveZero; /* if sign bit set, prints "-0" */
  if (decimalPoint) {
    /* scientific or fixed point (starts with - or blank) */
    int i;
    sprintf(label, fixedFormat, iValue);
    for (i=0 ; i<decimalPoint ; i++) label[i]= label[i+1];
    label[decimalPoint]= '.';
    /* get rid of leading zero sometimes */
    for (i=0 ; i<decimalPoint-1 ; i++) {
      if (label[i]<'0' || label[i]>'9') continue;
      if (label[i]=='0') {
        do {
          label[i]= label[i+1];
          i++;
        } while (label[i]);
      }
      break;
    }
  } else if (overflowChar) {
    /* overflow (beyond nDigits) */
    sprintf(label, fixedFormat, iValue);
    label[0]= '?';
    label[1]= label[overflowChar];
    label[2]= label[overflowChar+1];
    label[3]= '\0';
  } else if (niceDecades) {
    /* All decades between 0.001 and 1000. inclusive */
    strcpy(label, niceDecs[3+(int)iValue]);
  } else if (useLog) {
    /* E-12 style decade */
    sprintf(label, gtDoEscapes? nicedFormat:decadeFormat, iValue);
  } else {
    /* integer (starts with - or blank) */
    sprintf(label, integerFormat, iValue);
  }

  /* increment label value for next tick */
  iValue+= (GpReal)jBase;
}

static void NiceLogLabel(char *label, char *format, GpReal value, GpReal tick)
{
       /* "0.001", "0.01", "0.1", "1.0", "10.0", "100.0", "1000." */
  int decade= (int)floor(tick);
  char *src, *dst, scratch[32];

  sprintf(scratch, format, value);
  src= scratch;
  dst= label;

  if (decade<0) {
    *dst++= '0';   *dst++= '.';
    while (++decade) *dst++= '0';
    while ((*dst++= *src++)) if (*src=='.') src++;
    if (dst[-2]=='0') dst[-2]= '\0';

  } else {  /* decade>=0 */
    int digits= 0;
    while ((*dst++= *src++)) {
      if (*src=='.') src++;
      if (digits==decade) *dst++= '.';
      digits++;
    }
    if (digits<=decade) {
      dst--;
      while (digits<=decade) { *dst++= '0';  digits++; }
      *dst++= '.';
      *dst= '\0';
    }
  }
}

static void GtText(GpReal x0, GpReal y0, const char *text)
{
  if (*text=='+') text++;
  GpText(x0, y0, text);
}

static int DrawXLabels(int isLog, GpReal place, int upper, int nDigits,
                       GpTextAttribs *attribs)
{
  int i, nLabel= nChangeLevel[0];
  char label[32], expspace[16], *exponent;
  GpReal tick0= 0.0;
  int altflag;

  gistA.t= *attribs;
  gistA.t.alignH= TH_CENTER;
  gistA.t.alignV= upper? TV_BASE : TV_CAP;

  exponent= expspace+2;
  InitLabels(nDigits, exponent); /* also sets overflow */
  if (gtDoEscapes && exponent[0]) {
    exponent= expspace;
    exponent[0]= '1';
    exponent[1]= '0';
    exponent[2]= '^';
  }

  altflag= (altlabel && !overflow[0] && !isLog);
  for (i=0 ; altflag && i<nLabel ; i++)
    altflag= !altlabel((char *)0, ticks[i]);

  for (i=0 ; i<nLabel ; i++) {
    if (!altflag) NextLabel(label);
    else altlabel(label, ticks[i]);
    if (i==0 && nLabel==1) {
      tick0= scalx*ticks[i]+offx;
      if (useLog && subDecadeTicks) omitX= -1;
    }
    if (i!=omitX) GtText(scalx*ticks[i]+offx, place, label);
  }

  if (useLog && subDecadeTicks) {
    GpReal value, unit, test= 0.0;
    int j;
    for (j=0 ; j<3 ; j++) if (nSubtotal[j]>nLabel) break;
    sprintf(fixedFormat, "%%.%df", -subPower[j]);
    nLabel= nSubtotal[2];   /* actual total count of labels */
    for (j=0 ; j<3 ; j++) { /* loop on (1,2] (2,5] (5,10) intervals */
      value= subTick0[j];
      unit= subUnit[j];
      if (i<nSubtotal[j]) {
        if (niceDecades) NiceLogLabel(label, fixedFormat, value, ticks[i]);
        else sprintf(label, fixedFormat, value);
        test= ticks[i];
      }
      for ( ; i<nSubtotal[j] ; i++) {
        if (test > ticks[i]+EPS_LOG) { /* i.e.- test!=ticks[i] */
          value+= unit;
          if (!niceDecades) sprintf(label, fixedFormat, value);
          test= ticks[i];
        }
        if (niceDecades) NiceLogLabel(label, fixedFormat, value, ticks[i]);
        if (i==0 && nLabel==1) tick0= scalx*ticks[i]+offx;
        if (i!=omitX) GtText(scalx*ticks[i]+offx, place, label);
        test+= 1.0;   /* this is what EvenlySpace originally did... */
      }
    }
  }

  if (exponent[0]) {
    /* if present, exponent goes full line below X labels */
    GpReal exploc;
    GpReal nudged= place + (upper? +gistA.t.height : -gistA.t.height);
    if (nLabel>1) {
      exploc= 0.5*(ticks[nLabel-1]+ticks[nLabel-2]);
      exploc= scalx*exploc+offx;
    } else {
      if (fabs(tick0-gistT.viewport.xmin) >
          fabs(tick0-gistT.viewport.xmax))
        exploc= 0.5*(tick0+gistT.viewport.xmin);
      else
        exploc= 0.5*(tick0+gistT.viewport.xmax);
    }
    GpText(exploc, nudged, exponent);
  }

  return overflow[0]!='\0';
}

static int DrawYLabels(int isLog, GpReal place, int upper, int nDigits,
                       GpTextAttribs *attribs)
{
  int i, nLabel= nChangeLevel[0];
  char label[32], expspace[16], *exponent;
  GpReal tick0= 0.0;
  int altflag;

  gistA.t= *attribs;
  gistA.t.alignH= upper? TH_LEFT : TH_RIGHT;
  gistA.t.alignV= TV_HALF;

  exponent= expspace+2;
  InitLabels(nDigits, exponent); /* also sets overflow */
  if (gtDoEscapes && exponent[0]) {
    exponent= expspace;
    exponent[0]= '1';
    exponent[1]= '0';
    exponent[2]= '^';
  }

  if (overflow[0]) {
    overflow[0]= 'y';
    /* put special marker where exponent would go --
       this is on y-axis only, since overflow text is likely to be
       "far away" from the y-axis labels, but "close to" x-axis labels */
    strcpy(exponent, "y0+?");
  }

  altflag= (altlabel && !overflow[0] && !isLog);
  for (i=0 ; altflag && i<nLabel ; i++)
    altflag= !altlabel((char *)0, ticks[i]);

  for (i=0 ; i<nLabel ; i++) {
    if (!altflag) NextLabel(label);
    else altlabel(label, ticks[i]);
    if (i==0 && nLabel==1) {
      tick0= scaly*ticks[i]+offy;
      if (useLog && subDecadeTicks) omitY= -1;
    }
    if (i!=omitY) GtText(place, scaly*ticks[i]+offy, label);
  }

  if (useLog && subDecadeTicks) {
    GpReal value, unit, test= 0.0;
    int j;
    for (j=0 ; j<3 ; j++) if (nSubtotal[j]>nLabel) break;
    sprintf(fixedFormat, "%%.%df", -subPower[j]);
    nLabel= nSubtotal[2];   /* actual total count of labels */
    for (j=0 ; j<3 ; j++) { /* loop on (1,2] (2,5] (5,10) intervals */
      value= subTick0[j];
      unit= subUnit[j];
      if (i<nSubtotal[j]) {
        if (niceDecades) NiceLogLabel(label, fixedFormat, value, ticks[i]);
        else sprintf(label, fixedFormat, value);
        test= ticks[i];
      }
      for ( ; i<nSubtotal[j] ; i++) {
        if (test > ticks[i]+EPS_LOG) { /* i.e.- test!=ticks[i] */
          value+= unit;
          if (!niceDecades) sprintf(label, fixedFormat, value);
          test= ticks[i];
        }
        if (niceDecades) NiceLogLabel(label, fixedFormat, value, ticks[i]);
        if (i==0 && nLabel==1) tick0= scalx*ticks[i]+offx;
        if (i!=omitY) GtText(place, scaly*ticks[i]+offy, label);
        test+= 1.0;   /* this is what EvenlySpace originally did... */
      }
    }
  }

  if (exponent[0]) {
    /* if present, exponent goes slightly left of Y labels */
    GpReal exploc;
    GpReal nudged= place + (upper? +0.4 : -0.4)*gistA.t.height;
    if (nLabel>1) {
      exploc= 0.5*(ticks[nLabel-1]+ticks[nLabel-2]);
      exploc= scaly*exploc+offy;
    } else {
      if (fabs(tick0-gistT.viewport.ymin) >
          fabs(tick0-gistT.viewport.ymax))
        exploc= 0.5*(tick0+gistT.viewport.ymin);
      else
        exploc= 0.5*(tick0+gistT.viewport.ymax);
    }
    GpText(nudged, exploc, exponent);
  }

  return overflow[0]!='\0';
}

static void DrawOverflow(GpReal x, GpReal y)
{
  /* Just finished drawing labels (DrawXLabels, DrawYLabels), so text
     attributes should all be set properly except for alignment */
  gistA.t.alignH= TH_NORMAL;
  gistA.t.alignV= TV_NORMAL;
  GpText(x, y, overflow);
}

/* ------------------------------------------------------------------------ */

static GpTransform saveMap;
static GpTransform unitMap= {{0., 2., 0., 2.}, {0., 2., 0., 2.}};

int GaTicks(GaTickStyle *ticks, int xIsLog, int yIsLog)
{
  return GaAltTick(ticks, xIsLog, yIsLog,
                   (GaAltTicks *)0, (GaAltLabel *)0,
                   (GaAltTicks *)0, (GaAltLabel *)0);
}

int GaAltTick(GaTickStyle *ticks, int xIsLog, int yIsLog,
              GaAltTicks *xtick, GaAltLabel *xlabel,
              GaAltTicks *ytick, GaAltLabel *ylabel)
{
       /* Draws a system of tick marks and labels for the current
          transformation (gistT), according to the specified
          tick style.  */

  GpReal xmin= gistT.viewport.xmin;
  GpReal xmax= gistT.viewport.xmax;
  GpReal ymin= gistT.viewport.ymin;
  GpReal ymax= gistT.viewport.ymax;

  GpReal wxmin= gistT.window.xmin;
  GpReal wxmax= gistT.window.xmax;
  GpReal wymin= gistT.window.ymin;
  GpReal wymax= gistT.window.ymax;

  scalx= (xmax-xmin)/(wxmax-wxmin);
  offx= xmin - scalx*wxmin;
  scaly= (ymax-ymin)/(wymax-wymin);
  offy= ymin - scaly*wymin;

  if (xmin>xmax) { GpReal tmp= xmin; xmin= xmax; xmax= tmp; }
  if (ymin>ymax) { GpReal tmp= ymin; ymin= ymax; ymax= tmp; }

  saveMap= gistT;
  GpSetTrans(&unitMap);
  gistClip= 0;

  if (ticks->horiz.flags & TICK_ANY) {
    altticks= (ticks->horiz.flags & ALT_TICK)? xtick : 0;
    altlabel= (ticks->horiz.flags & ALT_LABEL)? xlabel : 0;
    FindTicks(wxmin, wxmax, &ticks->horiz, xIsLog);

    if (ticks->horiz.flags & TICK_C) {
      /* Axis and ticks to be drawn in middle of viewport */
      GpReal yOrigin= FindOrigin(wymin, wymax, &ticks->vert, yIsLog);
      yOrigin= scaly*yOrigin+offy;

      DrawXTicks(yOrigin, ticks->horiz.tickLen, TICK_IN | TICK_OUT, 0,
                 &ticks->horiz.tickStyle);

      if ((ticks->horiz.flags & LABEL_ANY)) {
        omitX= FindOmitted(wxmin, wxmax);
        if (DrawXLabels(xIsLog, yOrigin-ticks->horiz.labelOff, 0,
                      ticks->horiz.nDigits, &ticks->horiz.textStyle))
          DrawOverflow(ticks->horiz.xOver, ticks->horiz.yOver);
      }

    } else {
      /* Axis and ticks to be drawn around edges of viewport */
      int ovfl, jp, jb;

      if (ticks->horiz.flags & TICK_L)
        DrawXTicks(ymin-ticks->horiz.tickOff, ticks->horiz.tickLen,
                   ticks->horiz.flags, 0, &ticks->horiz.tickStyle);
      if (ticks->horiz.flags & TICK_U)
        DrawXTicks(ymax+ticks->horiz.tickOff, ticks->horiz.tickLen,
                   ticks->horiz.flags, 1, &ticks->horiz.tickStyle);

      omitX= -1;
      jp= jPower;  jb= jBase;
      ovfl= ((ticks->horiz.flags & LABEL_L) &&
             DrawXLabels(xIsLog, ymin-ticks->horiz.labelOff, 0,
                         ticks->horiz.nDigits, &ticks->horiz.textStyle));
      jPower= jp;  jBase= jb;
      ovfl|= ((ticks->horiz.flags & LABEL_U) &&
              DrawXLabels(xIsLog, ymax+ticks->horiz.labelOff, 1,
                          ticks->horiz.nDigits, &ticks->horiz.textStyle));
      if (ovfl)
        DrawOverflow(ticks->horiz.xOver, ticks->horiz.yOver);
    }
  }

  if (ticks->horiz.flags & GRID_F) {
    GpReal gridLen[TICK_LEVELS];
    GpReal delta= ymax-ymin;
    int i;
    if (nLevel>ticks->horiz.gridLevel)
      nTotal= nChangeLevel[ticks->horiz.gridLevel];
    for (i=0 ; i<=nLevel ; i++) gridLen[i]= delta;
    DrawXTicks(ymin, gridLen, TICK_IN, 0, &ticks->horiz.gridStyle);
  } else if (ticks->horiz.flags & GRID_O) {
    GpReal xOrigin= FindOrigin(wxmin, wxmax, &ticks->horiz, xIsLog);
    xOrigin= scalx*xOrigin+offx;
    DrawOrigin(xOrigin, xOrigin, ymin, ymax, &ticks->horiz.gridStyle);
  }

  if (ticks->vert.flags & TICK_ANY) {
    altticks= (ticks->horiz.flags & ALT_TICK)? ytick : 0;
    altlabel= (ticks->horiz.flags & ALT_LABEL)? ylabel : 0;
    FindTicks(wymin, wymax, &ticks->vert, yIsLog);

    if (ticks->vert.flags & TICK_C) {
      /* Axis and ticks to be drawn in middle of viewport */
      GpReal xOrigin= FindOrigin(wxmin, wxmax, &ticks->horiz, xIsLog);
      xOrigin= scalx*xOrigin+offx;

      DrawYTicks(xOrigin, ticks->vert.tickLen, TICK_IN | TICK_OUT, 0,
                 &ticks->vert.tickStyle);

      if ((ticks->vert.flags & LABEL_ANY)) {
        omitY= FindOmitted(wymin, wymax);
        if (DrawYLabels(yIsLog, xOrigin-ticks->vert.labelOff, 0,
                        ticks->vert.nDigits, &ticks->vert.textStyle))
          DrawOverflow(ticks->vert.xOver, ticks->vert.yOver);
      }

    } else {
      /* Axis and ticks to be drawn around edges of viewport */
      int ovfl, jp, jb;

      if (ticks->vert.flags & TICK_L)
        DrawYTicks(xmin-ticks->vert.tickOff, ticks->vert.tickLen,
                   ticks->vert.flags, 0, &ticks->vert.tickStyle);
      if (ticks->vert.flags & TICK_U)
        DrawYTicks(xmax+ticks->vert.tickOff, ticks->vert.tickLen,
                   ticks->vert.flags, 1, &ticks->vert.tickStyle);

      omitY= -1;
      jp= jPower;  jb= jBase;
      ovfl= ((ticks->vert.flags & LABEL_L) &&
             DrawYLabels(xIsLog, xmin-ticks->vert.labelOff, 0,
                         ticks->vert.nDigits, &ticks->vert.textStyle));
      jPower= jp;  jBase= jb;
      ovfl|= ((ticks->vert.flags & LABEL_U) &&
              DrawYLabels(xIsLog, xmax+ticks->vert.labelOff, 1,
                          ticks->vert.nDigits, &ticks->vert.textStyle));
      if (ovfl)
        DrawOverflow(ticks->vert.xOver, ticks->vert.yOver);
    }
  }

  if (ticks->vert.flags & GRID_F) {
    GpReal gridLen[TICK_LEVELS];
    GpReal delta= xmax-xmin;
    int i;
    if (nLevel>ticks->vert.gridLevel)
      nTotal= nChangeLevel[ticks->vert.gridLevel];
    for (i=0 ; i<=nLevel ; i++) gridLen[i]= delta;
    DrawYTicks(xmin, gridLen, TICK_IN, 0, &ticks->vert.gridStyle);
  } else if (ticks->vert.flags & GRID_O) {
    GpReal yOrigin= FindOrigin(wymin, wymax, &ticks->vert, yIsLog);
    yOrigin= scaly*yOrigin+offy;
    DrawOrigin(xmin, xmax, yOrigin, yOrigin, &ticks->vert.gridStyle);
  }

  if (ticks->frame) {
    GpReal x0[4], y0[4], x1[4], y1[4];
    x0[0]= x1[2]= x0[3]= x1[3]= xmin-ticks->vert.tickOff;
    x1[0]= x0[1]= x1[1]= x0[2]= xmax+ticks->vert.tickOff;
    y0[0]= y1[0]= y0[1]= y1[3]= ymin-ticks->horiz.tickOff;
    y1[1]= y0[2]= y1[2]= y0[3]= ymax+ticks->horiz.tickOff;
    gistA.l= ticks->frameStyle;
    GpDisjoint(4L, x0, y0, x1, y1);
  }

  GpSetTrans(&saveMap);
  return 0;
}

/* ------------------------------------------------------------------------ */
