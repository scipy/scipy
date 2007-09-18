#include "nifti1.h"   /* for the NIFTI_INTENT_* constants */
#include "nifticdf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main( int argc , char *argv[] )
{
   double val , p , p1=0.0,p2=0.0,p3=0.0 ;
   double vbot,vtop,vdel ;
   int code , iarg=1 , doq=0 , dod=0 , doi=0 , doz=0 , doh=0 ;

   /*-- print some help for the pitiful user --*/

   if( argc < 3 || strstr(argv[1],"help") != NULL ){
    int ii ;
    printf("\n") ;
    printf("Demo program for computing NIfTI statistical functions.\n") ;
    printf("Usage: nifti_stats [-q|-d|-1|-z] val CODE [p1 p2 p3]\n") ;
    printf(" val can be a single number or in the form bot:top:step.\n") ;
    printf(" default ==> output p = Prob(statistic < val).\n") ;
    printf("  -q     ==> output is 1-p.\n") ;
    printf("  -d     ==> output is density.\n") ;
    printf("  -1     ==> output is x such that Prob(statistic < x) = val.\n") ;
    printf("  -z     ==> output is z such that Normal cdf(z) = p(val).\n") ;
    printf("  -h     ==> output is z such that 1/2-Normal cdf(z) = p(val).\n");
    printf(" Allowable CODEs:\n") ;
    for( ii=NIFTI_FIRST_STATCODE ; ii <= NIFTI_LAST_STATCODE ; ii++ ){
     printf("  %-10s",inam[ii]); if((ii-NIFTI_FIRST_STATCODE)%6==5)printf("\n");
    }
    printf("\n") ;
    printf(" Following CODE are distributional parameters, as needed.\n");
    printf("\n") ;
    printf("Results are written to stdout, 1 number per output line.\n") ;
    printf("Example (piping output into AFNI program 1dplot):\n") ;
    printf(" nifti_stats -d 0:4:.001 INVGAUSS 1 3 | 1dplot -dx 0.001 -stdin\n");
    printf("\n") ;
    printf("Author - RW Cox - SSCC/NIMH/NIH/DHHS/USA/EARTH - March 2004\n") ;
    printf("\n") ;
    exit(0) ;
   }

   /*-- check first arg to see if it is an output option;
        if so, set the appropriate output flag to determine what to compute --*/

        if( strcmp(argv[iarg],"-q") == 0 ){ doq = 1 ; iarg++ ; }
   else if( strcmp(argv[iarg],"-d") == 0 ){ dod = 1 ; iarg++ ; }
   else if( strcmp(argv[iarg],"-1") == 0 ){ doi = 1 ; iarg++ ; }
   else if( strcmp(argv[iarg],"-z") == 0 ){ doz = 1 ; iarg++ ; }
   else if( strcmp(argv[iarg],"-h") == 0 ){ doh = 1 ; iarg++ ; }

   /*-- get the value(s) to process --*/

   vbot=vtop=vdel = 0.0 ;
   sscanf( argv[iarg++] , "%lf:%lf:%lf" , &vbot,&vtop,&vdel ) ;
   if( vbot >= vtop ) vdel = 0.0 ;
   if( vdel <= 0.0  ) vtop = vbot ;

   /*-- decode the CODE into the integer signifying the distribution --*/

   code = nifti_intent_code(argv[iarg++]) ;
     if( code < 0 ){ fprintf(stderr,"illegal code=%s\n",argv[iarg-1]); exit(1); }

   /*-- get the parameters, if present (defaults are 0) --*/

   if( argc > iarg ) p1 = strtod(argv[iarg++],NULL) ;
   if( argc > iarg ) p2 = strtod(argv[iarg++],NULL) ;
   if( argc > iarg ) p3 = strtod(argv[iarg++],NULL) ;

   /*-- loop over input value(s), compute output, write to stdout --*/

   for( val=vbot ; val <= vtop ; val += vdel ){
     if( doq )                                        /* output = 1-cdf */
       p = nifti_stat2rcdf( val , code,p1,p2,p3 ) ;
     else if( dod )                                   /* output = density */
       p = 1000.0*( nifti_stat2cdf(val+.001,code,p1,p2,p3)
                   -nifti_stat2cdf(val     ,code,p1,p2,p3)) ;
     else if( doi )                                   /* output = inverse */
       p = nifti_cdf2stat( val , code,p1,p2,p3 ) ;
     else if( doz )                                   /* output = z score */
       p = nifti_stat2zscore( val , code,p1,p2,p3 ) ;
     else if( doh )                                   /* output = halfz score */
       p = nifti_stat2hzscore( val , code,p1,p2,p3 ) ;
     else                                              /* output = cdf */
       p = nifti_stat2cdf( val , code,p1,p2,p3 ) ;

     printf("%.9g\n",p) ;
     if( vdel <= 0.0 ) break ;  /* the case of just 1 value */
   }

   /*-- terminus est --*/

   exit(0) ;
}
