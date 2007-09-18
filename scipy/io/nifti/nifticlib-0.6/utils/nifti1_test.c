#include <nifti1_io.h>   /* directly include I/O library functions */

/*-----------------------------------------------*/
/*    cc -o nifti1_test -O2 nifti1_test.c -lm    */
/*-----------------------------------------------*/

/****************************************************************************/

int main( int argc , char *argv[] )
{
   nifti_image *nim ;
   int iarg=1 , outmode=1 , ll , argn=1, usegzip=0;
   char *tmpstr;
   
   if( argc < 2 || strcmp(argv[1],"-help") == 0 ){
     printf("Usage: nifti1_test [-n2|-n1|-na|-a2] infile [prefix]\n"
            "\n"
            " If prefix is given, then the options mean:\n"
            "  -a2 ==> write an ANALYZE 7.5 file pair: prefix.hdr/prefix.img\n"
            "  -n2 ==> write a NIFTI-1 file pair: prefix.hdr/prefix.img\n"
            "  -n1 ==> write a NIFTI-1 single file: prefix.nii\n"
            "  -na ==> write a NIFTI-1 ASCII+binary file: prefix.nia\n");
     printf("  -za2 => write an ANALYZE 7.5 file pair:\n"
            "          prefix.hdr.gz/prefix.img.gz\n"
            "  -zn2 => write a NIFTI-1 file pair: prefix.hdr.gz/prefix.img.gz\n"
            "  -zn1 => write a NIFTI-1 single file: prefix.nii.gz\n"
            "  -zna => write a NIFTI-1 ASCII+binary file: prefix.nia.gz\n");
     printf(" The default is '-n1'.\n"
            "\n"
            " If prefix is not given, then the header info from infile\n"
            " file is printed to stdout.\n"
            "\n"
            " Please note that the '.nia' format is NOT part of the\n"
            " NIFTI-1 specification, but is provided mostly for ease\n"
            " of visualization (e.g., you can edit a .nia file and\n"
            " change some header fields, then rewrite it as .nii)\n"
           ) ;
     printf("\nsizeof(nifti_1_header)=%u\n",
            (unsigned int)sizeof(nifti_1_header)) ;
     exit(0) ;
   }

   if( argv[1][0] == '-' ){
     argn=1;
     if (argv[1][1] == 'z' ) {
       usegzip = 1;
       argn=2;
     }
     if( argv[1][argn] == 'a' ){
       outmode = 0 ;
     } else if( argv[1][argn] == 'n' ){
       switch( argv[1][argn+1] ){
         case '1': outmode = 1 ; break ;
         default:  outmode = 2 ; break ;
         case 'a': outmode = 3 ; break ;
       }
     }
     iarg++ ;
   }

   if( iarg >= argc ){
     fprintf(stderr,"** ERROR: no input file on command line!?\n"); exit(1);
   }

   nim = nifti_image_read( argv[iarg++] , 1 ) ;
   if( nim == NULL ) exit(1) ;

   if( iarg >= argc ){ nifti_image_infodump(nim); exit(0); }

   nim->nifti_type = outmode ;
   if( nim->fname != NULL ) free(nim->fname) ;
   if( nim->iname != NULL ) free(nim->iname) ;

   ll = strlen(argv[iarg]) ;
   tmpstr = nifti_makebasename(argv[iarg]);
   nim->fname = (char *)calloc(1,ll+8) ; strcpy(nim->fname,tmpstr) ;
   nim->iname = (char *)calloc(1,ll+8) ; strcpy(nim->iname,tmpstr) ;
   free(tmpstr);
   if( nim->nifti_type == 1 ){
     strcat(nim->fname,".nii") ;
     strcat(nim->iname,".nii") ;
   } else if ( nim->nifti_type == 3 ){
     strcat(nim->fname,".nia") ;
     strcat(nim->iname,".nia") ;
   } else {
     strcat(nim->fname,".hdr") ;
     strcat(nim->iname,".img") ;
   }
   if (usegzip) {
     strcat(nim->fname,".gz");
     strcat(nim->iname,".gz");
   }
   nifti_image_write( nim ) ;
   nifti_image_free( nim ) ;
   exit(0) ;
}
