
 /************************************************************************/
 /**  Functions to compute cumulative distributions and their inverses  **/
 /**  for the NIfTI-1 statistical types.  Much of this code is taken    **/
 /**  from other sources.  In particular, the cdflib functions by       **/
 /**  Brown and Lovato make up the bulk of this file.  That code        **/
 /**  was placed in the public domain.  The code by K. Krishnamoorthy   **/
 /**  is also released for unrestricted use.  Finally, the other parts  **/
 /**  of this file (by RW Cox) are released to the public domain.       **/
 /**                                                                    **/
 /**  Most of this file comprises a set of "static" functions, to be    **/
 /**  called by the user-level functions at the very end of the file.   **/
 /**  At the end of the file is a simple main program to drive these    **/
 /**  functions.                                                        **/
 /**                                                                    **/
 /**  To find the user-level functions, search forward for the string   **/
 /**  "nifti_", which will be at about line 11000.                      **/
 /************************************************************************/
 /*****==============================================================*****/
 /***** Neither the National Institutes of Health (NIH), the DFWG,   *****/
 /***** nor any of the members or employees of these institutions    *****/
 /***** imply any warranty of usefulness of this material for any    *****/
 /***** purpose, and do not assume any liability for damages,        *****/
 /***** incidental or otherwise, caused by any use of this document. *****/
 /***** If these conditions are not acceptable, do not use this!     *****/
 /*****==============================================================*****/
 /************************************************************************/

#include <stdlib.h>

/****************************************************************************
 Statistical codes implemented :

     NIFTI_INTENT_CORREL     = correlation statistic
     NIFTI_INTENT_TTEST      = t statistic (central)
     NIFTI_INTENT_FTEST      = F statistic (central)
     NIFTI_INTENT_ZSCORE     = N(0,1) statistic
     NIFTI_INTENT_CHISQ      = Chi-squared (central)
     NIFTI_INTENT_BETA       = Beta variable (central)
     NIFTI_INTENT_BINOM      = Binomial variable
     NIFTI_INTENT_GAMMA      = Gamma distribution
     NIFTI_INTENT_POISSON    = Poisson distribution
     NIFTI_INTENT_FTEST_NONC = noncentral F statistic
     NIFTI_INTENT_CHISQ_NONC = noncentral chi-squared
     NIFTI_INTENT_TTEST_NONC = noncentral t statistic
     NIFTI_INTENT_CHI        = Chi statistic (central)
     NIFTI_INTENT_INVGAUSS   = inverse Gaussian variable
     NIFTI_INTENT_WEIBULL    = Weibull distribution
     NIFTI_INTENT_EXTVAL     = Extreme value type I
     NIFTI_INTENT_NORMAL     = N(mu,variance) normal
     NIFTI_INTENT_LOGISTIC   = Logistic distribution
     NIFTI_INTENT_LAPLACE    = Laplace distribution
     NIFTI_INTENT_UNIFORM    = Uniform distribution
     NIFTI_INTENT_PVAL       = "p-value"
     NIFTI_INTENT_LOGPVAL    = -ln(p)
     NIFTI_INTENT_LOG10PVAL  = -log10(p)
*****************************************************************************/

extern char *inam[];

int nifti_intent_code( char *name );
double nifti_stat2cdf( double val, int code, double p1,double p2,double p3 );
double nifti_stat2rcdf( double val, int code, double p1,double p2,double p3 );
double nifti_cdf2stat( double p , int code, double p1,double p2,double p3 );
#if defined(__COMPILE_UNUSED_FUNCTIONS__)
double nifti_rcdf2stat( double q , int code, double p1,double p2,double p3 );
#endif/*(__COMPILE_UNUSED_FUNCTIONS__)*/
double nifti_stat2zscore( double val , int code, double p1,double p2,double p3);
double nifti_stat2hzscore( double val, int code, double p1,double p2,double p3);

/** Prototypes for cdflib functions **/

double algdiv(double*,double*);
double alngam(double*);
double alnrel(double*);
double apser(double*,double*,double*,double*);
double basym(double*,double*,double*,double*);
double bcorr(double*,double*);
double betaln(double*,double*);
double bfrac(double*,double*,double*,double*,double*,double*);
void bgrat(double*,double*,double*,double*,double*,double*,int*i);
double bpser(double*,double*,double*,double*);
void bratio(double*,double*,double*,double*,double*,double*,int*);
double brcmp1(int*,double*,double*,double*,double*);
double brcomp(double*,double*,double*,double*);
double bup(double*,double*,double*,double*,int*,double*);
void cdfbet(int*,double*,double*,double*,double*,double*,double*,
                   int*,double*);
void cdfbin(int*,double*,double*,double*,double*,double*,double*,
                   int*,double*);
void cdfchi(int*,double*,double*,double*,double*,int*,double*);
void cdfchn(int*,double*,double*,double*,double*,double*,int*,double*);
void cdff(int*,double*,double*,double*,double*,double*,int*,double*);
void cdffnc(int*,double*,double*,double*,double*,double*,double*,
                   int*s,double*);
void cdfgam(int*,double*,double*,double*,double*,double*,int*,double*);
#if defined(__COMPILE_UNUSED_FUNCTIONS__)
void cdfnbn(int*,double*,double*,double*,double*,double*,double*,
                   int*,double*);
void cdfnor(int*,double*,double*,double*,double*,double*,int*,double*);
#endif /*defined(__COMPILE_UNUSED_FUNCTIONS__)*/
void cdfpoi(int*,double*,double*,double*,double*,int*,double*);
void cdft(int*,double*,double*,double*,double*,int*,double*);
void cumbet(double*,double*,double*,double*,double*,double*);
void cumbin(double*,double*,double*,double*,double*,double*);
void cumchi(double*,double*,double*,double*);
void cumchn(double*,double*,double*,double*,double*);
void cumf(double*,double*,double*,double*,double*);
void cumfnc(double*,double*,double*,double*,double*,double*);
void cumgam(double*,double*,double*,double*);
#if defined(__COMPILE_UNUSED_FUNCTIONS__)
void cumnbn(double*,double*,double*,double*,double*,double*);
#endif /*defined(__COMPILE_UNUSED_FUNCTIONS__)*/
void cumnor(double*,double*,double*);
void cumpoi(double*,double*,double*,double*);
void cumt(double*,double*,double*,double*);
#if defined(__COMPILE_UNUSED_FUNCTIONS__)
double dbetrm(double*,double*);
#endif /*defined(__COMPILE_UNUSED_FUNCTIONS__)*/
double devlpl(double [],int*,double*);
#if defined(__COMPILE_UNUSED_FUNCTIONS__)
double dexpm1(double*);
double dinvnr(double *p,double *q);
#endif /*defined(__COMPILE_UNUSED_FUNCTIONS__)*/
void E0000(int,int*,double*,double*,unsigned long*,
                  unsigned long*,double*,double*,double*,
                  double*,double*,double*,double*);
void dinvr(int*,double*,double*,unsigned long*,unsigned long*);
void dstinv(double*,double*,double*,double*,double*,double*,
                   double*);
#if defined(__COMPILE_UNUSED_FUNCTIONS__)
double dlanor(double*);
double dln1mx(double*);
double dln1px(double*);
double dlnbet(double*,double*);
double dlngam(double*);
double dstrem(double*);
#endif /*defined(__COMPILE_UNUSED_FUNCTIONS__)*/
double dt1(double*,double*,double*);
void E0001(int,int*,double*,double*,double*,double*,
                  unsigned long*,unsigned long*,double*,double*,
                  double*,double*);
void dzror(int*,double*,double*,double*,double *,
                  unsigned long*,unsigned long*);
void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl);
double erf1(double*);
double erfc1(int*,double*);
double esum(int*,double*);
double exparg(int*);
double fpser(double*,double*,double*,double*);
double gam1(double*);
void gaminv(double*,double*,double*,double*,double*,int*);
double gamln(double*);
double gamln1(double*);
double Xgamm(double*);
void grat1(double*,double*,double*,double*,double*,double*);
void gratio(double*,double*,double*,double*,int*);
double gsumln(double*,double*);
double psi(double*);
double rcomp(double*,double*);
double rexp(double*);
double rlog(double*);
double rlog1(double*);
double spmpar(int*);
double stvaln(double*);
double fifdint(double);
double fifdmax1(double,double);
double fifdmin1(double,double);
double fifdsign(double,double);
long fifidint(double);
long fifmod(long,long);
void ftnstop(char*);
int ipmpar(int*);

/** end: prototypes for cdflib functions **/

