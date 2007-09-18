/*
 * test program for NIFTI lib.
 */
#include <nifti1_io.h>
enum NIFTITEST_BOOL {
  NIFTITEST_TRUE=1,
  NIFTITEST_FALSE=0
};

void _PrintTest(const int line,const char * message,const int FailureOccured, const enum NIFTITEST_BOOL isFatal,int *ErrorAccum)
{
  if(FailureOccured==NIFTITEST_TRUE)  /* This line can be commented out for a more verbose output */
    {
    char const * const PREFIX= (FailureOccured)?"==========ERROR":"..........SUCCESS";
    char const * const ISFATALPREFIX= (isFatal && FailureOccured)?" FATAL":"";
    printf("%s%s (LINE %d): %s\n",PREFIX,ISFATALPREFIX,line,message);
    fflush(stdout);
    *ErrorAccum+=FailureOccured;
    if(isFatal==NIFTITEST_TRUE && FailureOccured==NIFTITEST_TRUE)
      {
      printf("\n\nTOTAL ERRORS=%d\n",*ErrorAccum);
      exit( *ErrorAccum);
      }
    }
  return;
}
#define PrintTest(message,failure,isfailure,errorcount) \
  _PrintTest(__LINE__,message,failure,isfailure,errorcount)

nifti_image * generate_reference_image( const char * write_image_filename , int * const Errors)
{
  nifti_1_header reference_header;
  memset(&reference_header,0,sizeof(reference_header));
  reference_header.sizeof_hdr=sizeof(reference_header);
  reference_header.regular='r';
  reference_header.extents=16384;

  reference_header.dim[0]=5;
  reference_header.dim[1]=23;
  reference_header.dim[2]=17;
  reference_header.dim[3]=11;
  reference_header.dim[4]=7;
  reference_header.dim[5]=3;
  reference_header.dim[6]=1; //This MUST be 1 anything else is invalid due to code that usees huristics to fix other possible problems;
  reference_header.dim[7]=1; //This MUST be 1 anything else is invalid due to code that usees huristics to fix other possible problems;
  reference_header.intent_p1=10101010.101F;
  reference_header.intent_p2=987654321.0F;
  reference_header.intent_p3=-1234.0F;
  reference_header.intent_code=NIFTI_INTENT_ESTIMATE;
  reference_header.datatype=DT_INT32;
  reference_header.pixdim[0]=-1.0F; /* this is really qfac */
  reference_header.pixdim[1]=0.25F;
  reference_header.pixdim[2]=0.5F;
  reference_header.pixdim[3]=1.0F;
  reference_header.pixdim[4]=2.0F;
  reference_header.pixdim[5]=4.0F;
  reference_header.pixdim[6]=-2.0e10F;
  reference_header.pixdim[7]=-2.0e10F;
  reference_header.vox_offset=0;
  reference_header.scl_slope=0.25;
  reference_header.scl_inter=128;
  reference_header.qform_code=NIFTI_XFORM_SCANNER_ANAT;
  reference_header.quatern_b=-0.5F;
  reference_header.quatern_c= 0.5F;
  reference_header.quatern_d=-0.5F;
  reference_header.qoffset_x=reference_header.dim[1]/2.0F;
  reference_header.qoffset_y=reference_header.dim[2]/2.0F;
  reference_header.qoffset_z=reference_header.dim[3]/2.0F;
  reference_header.sform_code=NIFTI_XFORM_SCANNER_ANAT;
  reference_header.srow_x[0]=0.5;
  reference_header.srow_x[1]=0.0;
  reference_header.srow_x[2]=0.0;
  reference_header.srow_x[3]=0.0;
  reference_header.srow_y[0]=0.0;
  reference_header.srow_y[1]=1.0;
  reference_header.srow_y[2]=0.0;
  reference_header.srow_y[3]=0.0;
  reference_header.srow_z[0]=0.0;
  reference_header.srow_z[1]=0.0;
  reference_header.srow_z[2]=2.0;
  reference_header.srow_z[3]=0.0;
  reference_header.magic[0]='n';
  reference_header.magic[1]='+';
  reference_header.magic[2]='1';
  reference_header.magic[3]='\0';
  /* String is purposfully too long */
  strncpy(reference_header.intent_name,"PHANTOM_DATA to be used for regression testing the nifti reader/writer",16);
  strncpy(reference_header.descrip,"This is a very long dialog here to use up more than 80 characters of space to test to see if the code is robust enough to deal appropriatly with very long and obnoxious lines.",80);

  {
  int nbyper;
  int swapsize;
  nifti_datatype_sizes(reference_header.datatype ,&nbyper,&swapsize);
  reference_header.bitpix=nbyper*8;
  }

  nifti_image * reference_image=nifti_convert_nhdr2nim(reference_header,write_image_filename);
  {
  const unsigned int NumVoxels=reference_image->nx*reference_image->ny*reference_image->nz*reference_image->nt*reference_image->nu;
  reference_image->data=(signed int *)calloc(NumVoxels,sizeof(signed int)) ; /*!< pointer to data: nbyper*nvox bytes     */
  PrintTest("Checking memory allocation",reference_image->data ==0 ,NIFTITEST_TRUE,Errors);
  {
  signed int i=0;
  for(; i < (signed int)NumVoxels ; i++)
    {
    ((signed int *)(reference_image->data))[i]=i;
    }
  }
  }
  PrintTest("Setting filenames",nifti_set_filenames( reference_image,write_image_filename, 0, 0 ) != 0, NIFTITEST_TRUE,Errors);
  PrintTest("Setting type from names",nifti_set_type_from_names( reference_image ) != 0, NIFTITEST_TRUE,Errors);
  /*   PrintTest("Checking type and names",nifti_type_and_names_match( reference_image , 1 ) != 1, NIFTITEST_TRUE,Errors); */
  PrintTest("Check reference_image data is non null",(reference_image->data==0),NIFTITEST_TRUE,Errors);
  return reference_image;
}


void compare_reference_image_values(nifti_image const * const reference_image, nifti_image const * const reloaded_image, int * const Errors)
{
  PrintTest("Checking nifti_type",(reference_image->nifti_type!=reloaded_image->nifti_type),NIFTITEST_FALSE,Errors);
  PrintTest("Checking fname",(strcmp(reference_image->fname,reloaded_image->fname)),NIFTITEST_FALSE,Errors);
  PrintTest("Checking iname",(strcmp(reference_image->iname,reloaded_image->iname)),NIFTITEST_FALSE,Errors);
  PrintTest("Checking ndim",(reference_image->ndim!=reloaded_image->ndim),NIFTITEST_FALSE,Errors);
  PrintTest("Checking nx",(reference_image->nx!=reloaded_image->nx),NIFTITEST_FALSE,Errors);
  PrintTest("Checking ny",(reference_image->ny!=reloaded_image->ny),NIFTITEST_FALSE,Errors);
  PrintTest("Checking nz",(reference_image->nz!=reloaded_image->nz),NIFTITEST_FALSE,Errors);
  PrintTest("Checking nt",(reference_image->nt!=reloaded_image->nt),NIFTITEST_FALSE,Errors);
  PrintTest("Checking nu",(reference_image->nu!=reloaded_image->nu),NIFTITEST_FALSE,Errors);
  PrintTest("Checking dx",(reference_image->dx!=reloaded_image->dx),NIFTITEST_FALSE,Errors);
  PrintTest("Checking dy",(reference_image->dy!=reloaded_image->dy),NIFTITEST_FALSE,Errors);
  PrintTest("Checking dz",(reference_image->dz!=reloaded_image->dz),NIFTITEST_FALSE,Errors);
  PrintTest("Checking dt",(reference_image->dt!=reloaded_image->dt),NIFTITEST_FALSE,Errors);
  PrintTest("Checking du",(reference_image->du!=reloaded_image->du),NIFTITEST_FALSE,Errors);
  PrintTest("Checking datatype",(reference_image->datatype!=reloaded_image->datatype),NIFTITEST_FALSE,Errors);
  {
  const unsigned int NumVoxels=reference_image->nx*reference_image->ny*reference_image->nz*reference_image->nt*reference_image->nu;
  PrintTest("Check loaded data is non null",(reloaded_image->data==0),NIFTITEST_TRUE,Errors);
  PrintTest("Check reference_image data is non null",(reference_image->data==0),NIFTITEST_TRUE,Errors);
  {
  unsigned int CurrVoxel=0;
  for(; CurrVoxel < NumVoxels ; CurrVoxel++)
    {
    /*printf("%d ",CurrVoxel); fflush(stdout);*/
    if( ((int *)(reference_image->data))[CurrVoxel] !=  ((int *)(reloaded_image->data))[CurrVoxel]) 
      {
      PrintTest("Incorrect Pixel Value Found",0,NIFTITEST_FALSE,Errors);
      }
    }
  }
  }
  PrintTest("Checking xyz_units",(reference_image->xyz_units!=reloaded_image->xyz_units),NIFTITEST_FALSE,Errors);
  PrintTest("Checking time_units",(reference_image->time_units!=reloaded_image->time_units),NIFTITEST_FALSE,Errors);
  PrintTest("Checking intent_code",(reference_image->intent_code!=reloaded_image->intent_code),NIFTITEST_FALSE,Errors);
  PrintTest("Checking intent_name",(strncmp(reference_image->intent_name,reloaded_image->intent_name,16) )!=0,NIFTITEST_FALSE,Errors);
  PrintTest("Checking description",(strncmp(reference_image->descrip,reloaded_image->descrip,80))!=0,NIFTITEST_FALSE,Errors);
  return ;
}

int main (int argc, char *argv[])
{
  char TEMP_STR[256];
  nifti_set_debug_level(3);
  int Errors=0;
  {
  PrintTest("NOT REALLY AN ERROR, JUST TESTING THE ERROR TEST REPORTING MECHANISM",1,NIFTITEST_FALSE,&Errors);
  PrintTest("NOT REALLY AN ERROR, JUST TESTING THE ERROR COUNTING MECHANISM",Errors==1,NIFTITEST_FALSE,&Errors);
  Errors=0;
  }
  {
  const char write_image_filename[6][64]={
    "ATestReferenceImageForReadingAndWriting.nii",
    "ATestReferenceImageForReadingAndWriting.hdr",
    "ATestReferenceImageForReadingAndWriting.img",
    "ATestReferenceImageForReadingAndWriting.nii.gz",
    "ATestReferenceImageForReadingAndWriting.hdr.gz",
    "ATestReferenceImageForReadingAndWriting.img.gz"
  };
  printf("======= Testing All Nifti Valid Names ======\n");
  fflush(stdout);
  unsigned int filenameindex;
  for(filenameindex=0;filenameindex<6; filenameindex++)
    {
    char buf[512];
    int CompressedTwoFile = strstr(write_image_filename[filenameindex],".img.gz") != 0 ||
      strstr(write_image_filename[filenameindex],".hdr.gz") != 0;

    printf("======= Testing with filename: %s ======\n",write_image_filename[filenameindex]);
    fflush(stdout);
    nifti_image * reference_image = generate_reference_image(write_image_filename[filenameindex],&Errors);
    /*
     * Add an extension to test extension reading
     */
    {
    static char ext[] = "THIS IS A TEST";
    sprintf(buf,"nifti_add_extension %s",write_image_filename[filenameindex]);
    PrintTest(buf,
              nifti_add_extension(reference_image,
                                  ext,sizeof(ext),
                                  NIFTI_ECODE_COMMENT) == -1,
              NIFTITEST_FALSE,&Errors);
    sprintf(buf,"valid_nifti_extension %s",write_image_filename[filenameindex]);
    PrintTest("valid_nifti_extensions",
              valid_nifti_extensions(reference_image) == 0,
              NIFTITEST_FALSE,&Errors);
    }
    PrintTest("Create reference image",reference_image==0,NIFTITEST_TRUE,&Errors);
    nifti_image_write   ( reference_image ) ;
    /*
     * test nifti_copy_extension
     */
    {
      nifti_image *nim = nifti_simple_init_nim();
      PrintTest("nifti_copy_extension",
                nifti_copy_extensions(nim,reference_image),
                NIFTITEST_FALSE,&Errors);
      
      nifti_image_free(nim);
      nim = nifti_copy_nim_info(reference_image);
      PrintTest("nifti_copy_nim_info",
                nim == 0,
                NIFTITEST_FALSE,&Errors);
      PrintTest("nifti_nim_is_valid",
                nifti_nim_is_valid(nim,0) == 0,
                NIFTITEST_FALSE,&Errors);
      
      
      nifti_image_free(nim);
      
    }
    {
    nifti_image * reloaded_image = nifti_image_read(reference_image->fname,1);
    PrintTest("Reload of image ",reloaded_image==0,NIFTITEST_TRUE,&Errors);
    
    {
    /*
     * if the file is named '.img', '.hdr', '.img.gz', or '.hdr.gz', then
     * the header extensions won't be saved with the file.
     * The test will fail if it finds an extension in a 2-file NIfTI, or
     * fails to find one in a '.nii' or '.nii.gz' file.
     */
    int result = valid_nifti_extensions(reloaded_image);
    sprintf(buf,"reload valid_nifti_extensions %s",write_image_filename[filenameindex]);
      PrintTest(buf,
                CompressedTwoFile ? result != 0 : result == 0,
                NIFTITEST_FALSE,&Errors);
    }
    nifti_image_infodump(reloaded_image);
    compare_reference_image_values(reference_image,reloaded_image,&Errors);
    nifti_image_free(reloaded_image);
    }

    {
    nifti_brick_list   NB_orig, NB_select;

    nifti_image      * nim_orig, * nim_select;

    int                blist[5] = { 7, 0, 5, 5, 9 };

    /*
     * test some error paths in the nifti_image_read_bricks
     */
    nim_orig = nifti_image_read_bricks(reference_image->fname,0,blist, &NB_orig);
    PrintTest("invalid arg bricked image read 1",nim_orig != 0,NIFTITEST_FALSE,&Errors);

    nim_orig   = nifti_image_read_bricks(reference_image->fname, 0, NULL,  &NB_orig);
    PrintTest("Reload of bricked image",nim_orig == 0,NIFTITEST_FALSE,&Errors);
    nifti_free_NBL(&NB_orig);
    nifti_image_free(nim_orig);
    
    nim_select = nifti_image_read_bricks(reference_image->fname, 5, blist, &NB_select);
    PrintTest("Reload of bricked image with blist",nim_orig == 0,NIFTITEST_FALSE,&Errors);
    nifti_free_NBL(&NB_select);
    nifti_image_free(nim_select);

    }
    /*
     * test nifti_update_dims_from_array
     */
    PrintTest("nifti_update_dims_from_array -- valid dims",
              nifti_update_dims_from_array(reference_image) != 0,
              NIFTITEST_FALSE,&Errors);
    reference_image->dim[0] = 8;
    PrintTest("nifti_update_dims_from_array -- invalid dims",
              nifti_update_dims_from_array(reference_image) == 0,
              NIFTITEST_FALSE,&Errors);
    {
    nifti_1_header x = nifti_convert_nim2nhdr(reference_image);
    char buf[512];
    sprintf(buf,"nifti_hdr_looks_good %s",reference_image->fname);
    PrintTest(buf,
              !nifti_hdr_looks_good(&x),
              NIFTITEST_FALSE,&Errors);
    }

    nifti_image_free(reference_image);
    }
  /*
   * check nifti_findimgname
   */
  {
  char *imgname = nifti_findimgname("ATestReferenceImageForReadingAndWriting.hdr",2);
  PrintTest("nifti_findimgname",
            imgname == 0 || 
            strcmp(imgname,"ATestReferenceImageForReadingAndWriting.img") != 0,
            NIFTITEST_FALSE,&Errors);
  free(imgname);
  }
  {
  int IsNiftiFile;
  IsNiftiFile = is_nifti_file(write_image_filename[0]);
  PrintTest("is_nifti_file0",
            IsNiftiFile != 1,NIFTITEST_FALSE,&Errors);
  IsNiftiFile = is_nifti_file(write_image_filename[1]);
  PrintTest("is_nifti_file1",
            IsNiftiFile != 2,NIFTITEST_FALSE,&Errors);
  IsNiftiFile = is_nifti_file(write_image_filename[3]);
  PrintTest("is_nifti_file2",
            IsNiftiFile != 1,NIFTITEST_FALSE,&Errors);
  IsNiftiFile = is_nifti_file(write_image_filename[4]);
  PrintTest("is_nifti_file2",
            IsNiftiFile != 2,NIFTITEST_FALSE,&Errors);
  }
  
  }
  {
  /*
   * test writing and reading an ascii file
   */
  nifti_image * reference_image = 
    generate_reference_image("TestAsciiImage.nia",&Errors);
  reference_image->nifti_type = 3;
  nifti_image_write(reference_image);
  nifti_image * reloaded_image = nifti_image_read("TestAsciiImage.nia",1);
  PrintTest("Read/Write Ascii image",
            reloaded_image == 0,NIFTITEST_FALSE,&Errors);
  nifti_image_free(reference_image);
  nifti_image_free(reloaded_image);
  }

  {
  enum { NUM_FILE_NAMES=8 };
  const char * FILE_NAMES[NUM_FILE_NAMES]={
    "myimage",
    "myimage.tif",
    "myimage.tif.gz",
    "myimage.nii",
    "myimage.img.gz",
    ".nii",
    ".myhiddenimage",
    ".myhiddenimage.nii"
  };
  const char * KNOWN_FILE_BASENAMES[NUM_FILE_NAMES]={
    "myimage",
    "myimage.tif",
    "myimage.tif.gz",
    "myimage",
    "myimage",
    "",
    ".myhiddenimage",
    ".myhiddenimage"
  };
  const int KNOWN_nifti_validfilename[NUM_FILE_NAMES]={
    1,
    1,
    1,
    1,
    1,
    0,
    1,
    1
  };
  const int KNOWN_nifti_is_complete_filename[NUM_FILE_NAMES]={
    0,
    0,
    0,
    1,
    1,
    0,
    0,
    1
  };
  unsigned int fni;
  for(fni=0;fni<NUM_FILE_NAMES;fni++)
    {
    printf("\nTesting \"%s\" filename\n",FILE_NAMES[fni]);
    {
    int KnownValid=nifti_validfilename(FILE_NAMES[fni]);
    snprintf(TEMP_STR,256,"nifti_validfilename(\"%s\")=%d",FILE_NAMES[fni],KnownValid);
    PrintTest(TEMP_STR,KnownValid != KNOWN_nifti_validfilename[fni],NIFTITEST_FALSE,&Errors);
    }
    {
    int KnownValid=nifti_is_complete_filename(FILE_NAMES[fni]);
    snprintf(TEMP_STR,256,"nifti_is_complete_filename(\"%s\")=%d",FILE_NAMES[fni],KnownValid);
    PrintTest(TEMP_STR,KnownValid != KNOWN_nifti_is_complete_filename[fni],NIFTITEST_FALSE,&Errors);
    }

    {
    char * basename=nifti_makebasename(FILE_NAMES[fni]);
    snprintf(TEMP_STR,256,"nifti_makebasename(\"%s\")=\"%s\"",FILE_NAMES[fni],basename);
    PrintTest(TEMP_STR,strcmp(basename,KNOWN_FILE_BASENAMES[fni]) != 0,NIFTITEST_FALSE,&Errors);
    free(basename);
   
    }
    }
  /*
   * the following 2 calls aren't tested, because all they do is display
   * compile-time information -- no way to fail unless writing to stdout fails.
   */
  nifti_disp_lib_hist();
  nifti_disp_lib_version();
  /*
   * the following exercises error path code in nifti_image_read_bricks
   */
  PrintTest(
            "nifti_image_read_bricks 1",
            nifti_image_read_bricks((char *)0,-1,(const int *)0,(nifti_brick_list *)0) != 0,
            NIFTITEST_FALSE,
            &Errors);
  PrintTest(
            "nifti_image_read_bricks 1",
            nifti_image_read_bricks("NOFILE.NOFILE",-1,(const int *)0,(nifti_brick_list *)0) != 0,
            NIFTITEST_FALSE,
            &Errors);
  }
  /*
   * call nifti_datatype_string with all possible values
   */
#define nifti_datatype_test(constant,string)                            \
  {                                                                     \
  char buf[64];                                                         \
  sprintf(buf,"nifti_datatype_string %s",string);                       \
  PrintTest(                                                   \
            buf,                                                        \
            strcmp(nifti_datatype_string(constant),string) != 0,        \
            NIFTITEST_FALSE,                                            \
            &Errors);                                                   \
  }
  nifti_datatype_test(DT_UNKNOWN,"UNKNOWN"); 
  nifti_datatype_test(DT_BINARY, "BINARY");
  nifti_datatype_test(DT_INT8, "INT8");
  nifti_datatype_test(DT_UINT8, "UINT8");
  nifti_datatype_test(DT_INT16, "INT16");
  nifti_datatype_test(DT_UINT16, "UINT16");
  nifti_datatype_test(DT_INT32, "INT32");
  nifti_datatype_test(DT_UINT32, "UINT32");
  nifti_datatype_test(DT_INT64, "INT64");
  nifti_datatype_test(DT_UINT64, "UINT64");
  nifti_datatype_test(DT_FLOAT32, "FLOAT32");
  nifti_datatype_test(DT_FLOAT64, "FLOAT64");
  nifti_datatype_test(DT_FLOAT128, "FLOAT128");
  nifti_datatype_test(DT_COMPLEX64, "COMPLEX64");
  nifti_datatype_test(DT_COMPLEX128, "COMPLEX128");
  nifti_datatype_test(DT_COMPLEX256, "COMPLEX256");
  nifti_datatype_test(DT_RGB24, "RGB24");
#define nifti_is_inttype_test(constant,rval)            \
  {                                                     \
  char buf[64];                                         \
  sprintf(buf,"nifti_datatype_string %d",constant);     \
  PrintTest(                                   \
            buf,                                        \
            nifti_is_inttype(constant) != rval,         \
            NIFTITEST_FALSE,                            \
            &Errors);                                   \
  }
  nifti_is_inttype_test(DT_UNKNOWN,0);
  nifti_is_inttype_test(DT_BINARY,0);
  nifti_is_inttype_test(DT_INT8,1);
  nifti_is_inttype_test(DT_UINT8,1);
  nifti_is_inttype_test(DT_INT16,1);
  nifti_is_inttype_test(DT_UINT16,1);
  nifti_is_inttype_test(DT_INT32,1);
  nifti_is_inttype_test(DT_UINT32,1);
  nifti_is_inttype_test(DT_INT64,1);
  nifti_is_inttype_test(DT_UINT64,1);
  nifti_is_inttype_test(DT_FLOAT32,0);
  nifti_is_inttype_test(DT_FLOAT64,0);
  nifti_is_inttype_test(DT_FLOAT128,0);
  nifti_is_inttype_test(DT_COMPLEX64,0);
  nifti_is_inttype_test(DT_COMPLEX128,0);
  nifti_is_inttype_test(DT_COMPLEX256,0);
  nifti_is_inttype_test(DT_RGB24,1);
#define nifti_units_string_test(constant,string)                \
  {                                                             \
  char buf[64];                                                 \
  sprintf(buf,"nifti_units_string_test %s",string);             \
  PrintTest(                                           \
            buf,                                                \
            strcmp(nifti_units_string(constant),string) != 0,   \
            NIFTITEST_FALSE,                                    \
            &Errors);                                           \
  }
  nifti_units_string_test(NIFTI_UNITS_METER,"m");
  nifti_units_string_test(NIFTI_UNITS_MM,"mm");
  nifti_units_string_test(NIFTI_UNITS_MICRON,"um");
  nifti_units_string_test(NIFTI_UNITS_SEC,"s");
  nifti_units_string_test(NIFTI_UNITS_MSEC,"ms");
  nifti_units_string_test(NIFTI_UNITS_USEC,"us");
  nifti_units_string_test(NIFTI_UNITS_HZ,"Hz");
  nifti_units_string_test(NIFTI_UNITS_PPM,"ppm");
  nifti_units_string_test(NIFTI_UNITS_RADS,"rad/s");
#define nifti_intent_string_test(constant,string)               \
  {                                                             \
  char buf[64];                                                 \
  sprintf(buf,"nifti_intent_string %s",string);                 \
  PrintTest(                                           \
            buf,                                                \
            strcmp(nifti_intent_string(constant),string) != 0,  \
            NIFTITEST_FALSE,                                    \
            &Errors);                                           \
  }
  nifti_intent_string_test(NIFTI_INTENT_CORREL,"Correlation statistic");
  nifti_intent_string_test(NIFTI_INTENT_TTEST,"T-statistic");
  nifti_intent_string_test(NIFTI_INTENT_FTEST,"F-statistic");
  nifti_intent_string_test(NIFTI_INTENT_ZSCORE,"Z-score");
  nifti_intent_string_test(NIFTI_INTENT_CHISQ,"Chi-squared distribution");
  nifti_intent_string_test(NIFTI_INTENT_BETA,"Beta distribution");
  nifti_intent_string_test(NIFTI_INTENT_BINOM,"Binomial distribution");
  nifti_intent_string_test(NIFTI_INTENT_GAMMA,"Gamma distribution");
  nifti_intent_string_test(NIFTI_INTENT_POISSON,"Poisson distribution");
  nifti_intent_string_test(NIFTI_INTENT_NORMAL,"Normal distribution");
  nifti_intent_string_test(NIFTI_INTENT_FTEST_NONC,"F-statistic noncentral");
  nifti_intent_string_test(NIFTI_INTENT_CHISQ_NONC,"Chi-squared noncentral");
  nifti_intent_string_test(NIFTI_INTENT_LOGISTIC,"Logistic distribution");
  nifti_intent_string_test(NIFTI_INTENT_LAPLACE,"Laplace distribution");
  nifti_intent_string_test(NIFTI_INTENT_UNIFORM,"Uniform distribition");
  nifti_intent_string_test(NIFTI_INTENT_TTEST_NONC,"T-statistic noncentral");
  nifti_intent_string_test(NIFTI_INTENT_WEIBULL,"Weibull distribution");
  nifti_intent_string_test(NIFTI_INTENT_CHI,"Chi distribution");
  nifti_intent_string_test(NIFTI_INTENT_INVGAUSS,"Inverse Gaussian distribution");
  nifti_intent_string_test(NIFTI_INTENT_EXTVAL,"Extreme Value distribution");
  nifti_intent_string_test(NIFTI_INTENT_PVAL,"P-value");
  nifti_intent_string_test(NIFTI_INTENT_LOGPVAL,"Log P-value");
  nifti_intent_string_test(NIFTI_INTENT_LOG10PVAL,"Log10 P-value");
  nifti_intent_string_test(NIFTI_INTENT_ESTIMATE,"Estimate");
  nifti_intent_string_test(NIFTI_INTENT_LABEL,"Label index");
  nifti_intent_string_test(NIFTI_INTENT_NEURONAME,"NeuroNames index");
  nifti_intent_string_test(NIFTI_INTENT_GENMATRIX,"General matrix");
  nifti_intent_string_test(NIFTI_INTENT_SYMMATRIX,"Symmetric matrix");
  nifti_intent_string_test(NIFTI_INTENT_DISPVECT,"Displacement vector");
  nifti_intent_string_test(NIFTI_INTENT_VECTOR,"Vector");
  nifti_intent_string_test(NIFTI_INTENT_POINTSET,"Pointset");
  nifti_intent_string_test(NIFTI_INTENT_TRIANGLE,"Triangle");
  nifti_intent_string_test(NIFTI_INTENT_QUATERNION,"Quaternion");
  nifti_intent_string_test(NIFTI_INTENT_DIMLESS,"Dimensionless number");
  nifti_intent_string_test(-200,"Unknown");

#define nifti_slice_string_test(constant,string)                \
  {                                                             \
  char buf[64];                                                 \
  sprintf(buf,"nifti_slice_string_test %s",string);             \
  PrintTest(                                           \
            buf,                                                \
            strcmp(nifti_slice_string(constant),string) != 0,   \
            NIFTITEST_FALSE,                                    \
            &Errors);                                           \
  }
  nifti_slice_string_test(NIFTI_SLICE_SEQ_INC,"sequential_increasing");
  nifti_slice_string_test(NIFTI_SLICE_SEQ_DEC,"sequential_decreasing");
  nifti_slice_string_test(NIFTI_SLICE_ALT_INC,"alternating_increasing");
  nifti_slice_string_test(NIFTI_SLICE_ALT_DEC,"alternating_decreasing");
  nifti_slice_string_test(NIFTI_SLICE_ALT_INC2,"alternating_increasing_2");
  nifti_slice_string_test(NIFTI_SLICE_ALT_DEC2,"alternating_decreasing_2");
#define nifti_orientation_string_test(constant,string)                  \
  {                                                                     \
  char buf[64];                                                         \
  sprintf(buf,"nifti_orientation_string_test %s",string);               \
  PrintTest(                                                   \
            buf,                                                        \
            strcmp(nifti_orientation_string(constant),string) != 0,     \
            NIFTITEST_FALSE,                                            \
            &Errors);                                                   \
  }
  nifti_orientation_string_test(NIFTI_L2R,"Left-to-Right");
  nifti_orientation_string_test(NIFTI_R2L,"Right-to-Left");
  nifti_orientation_string_test(NIFTI_P2A,"Posterior-to-Anterior");
  nifti_orientation_string_test(NIFTI_A2P,"Anterior-to-Posterior");
  nifti_orientation_string_test(NIFTI_I2S,"Inferior-to-Superior");
  nifti_orientation_string_test(NIFTI_S2I,"Superior-to-Inferior");

#define nifti_datatype_sizes_test(constant,Nbyper,Swapsize)     \
  {                                                             \
  int nbyper;                                                   \
  int swapsize;                                                 \
  char buf[64];                                                 \
  sprintf(buf,"nifti_datatype_sizes_test %d",constant);         \
  nifti_datatype_sizes(constant,&nbyper,&swapsize);             \
  PrintTest(                                           \
            buf,                                                \
            nbyper != Nbyper || swapsize != Swapsize,           \
            NIFTITEST_FALSE,                                    \
            &Errors);                                           \
  }

  nifti_datatype_sizes_test(DT_UINT8,1,0);
  nifti_datatype_sizes_test(DT_UINT16,2,2);
  nifti_datatype_sizes_test(DT_RGB24,3,0);
  nifti_datatype_sizes_test(DT_FLOAT32,4,4);
  nifti_datatype_sizes_test(DT_COMPLEX64,8,4);
  nifti_datatype_sizes_test(DT_UINT64,8,8);
  nifti_datatype_sizes_test(DT_FLOAT128,16,16);
  nifti_datatype_sizes_test(DT_COMPLEX128,16,8);
  nifti_datatype_sizes_test(DT_COMPLEX256,32,16);

  {
  mat44 R;
  unsigned i,j;
  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
      R.m[i][j] = (i == j ? 1 : 0);
  float qb;
  float qc;
  float qd;
  float qx;
  float qy;
  float qz;
  float dx;
  float dy;
  float dz;
  float qfac;
  nifti_mat44_to_quatern(R,&qb,&qc,&qd,&qx,&qy,&qz,&dx,&dy,&dz,&qfac);
  PrintTest("nifti_mat44_to_quatern",
            qb != 0.000000 || qc != 0.000000 || qd != 0.000000 ||  
            qx != 0.000000 || qy != 0.000000 || qd != 0.000000 ||  
            dx != 1.000000 || dy != 1.000000 || dz != 1.000000 ||  
            qfac != 1.000000,
            NIFTITEST_FALSE,&Errors);
  }
  {
  mat44 x = nifti_make_orthog_mat44(0.14,0.0,0.0,
                                    0.0,0.9,0.0,
                                    0.0,0.0,1.1);
  
  PrintTest("nifti_make_orthog_mat44",
            x.m[0][0] != 1.0 || x.m[1][1] != 1.0 ||
            x.m[2][2] != 1.0 || x.m[3][3] != 1.0,
            NIFTITEST_FALSE,&Errors);
  }
  {
  static char x[16] = { 'a','b','c','d','e','f','g','h',
                        'H','G','F','E','D','C','B','A' };
  nifti_swap_Nbytes(1,16,(void *)x);
  PrintTest("nifti_swap_16bytes",
            x[0] != 'A' || x[1] != 'B' || x[2] != 'C' || x[3] != 'D' ||
            x[4] != 'E' || x[5] != 'F' || x[6] != 'G' || x[7] != 'H' ||
            x[8] != 'h' || x[9] != 'g' || x[10] != 'f' || x[11] != 'e' ||
            x[12] != 'd' || x[13] != 'c' || x[14] != 'b' || x[15] != 'a',
            NIFTITEST_FALSE,&Errors);
    
  }
  {
  static char x[8] = { 'a','b','c','d','D','C','B','A' };
  nifti_swap_Nbytes(1,8,(void *)x);
  PrintTest("nifti_swap_8bytes",
            x[0] != 'A' || x[1] != 'B' || x[2] != 'C' || x[3] != 'D' ||
            x[4] != 'd' || x[5] != 'c' || x[6] != 'b' || x[7] != 'a',
            NIFTITEST_FALSE,&Errors);
    
  }
  {
  /*
   * test nifti_simple_init_nim
   */
  nifti_image *nim = nifti_simple_init_nim();
  PrintTest("nifti_simple_init_nim",
            nim == 0,NIFTITEST_FALSE,&Errors);
  nifti_image_free(nim); nim = 0;
  /*
   * test nifti_image_open
   */
  znzFile f = nifti_image_open("ATestReferenceImageForReadingAndWriting.hdr","r",&nim);
  PrintTest("nifti_image_open",
            nim == 0 || f == 0,
            NIFTITEST_FALSE,&Errors);
  PrintTest("nifti_image_load",
            nifti_image_load(nim) == -1,
            NIFTITEST_FALSE,&Errors);
  nifti_image_unload(nim);
  PrintTest("nifti_image_unload",
            nim->data != 0,
            NIFTITEST_FALSE,&Errors);

  znzclose(f);
  nifti_image_free(nim);
  }
  /*
   * call various functions from nifti_stats
   */
  printf("\n\nTOTAL ERRORS=%d\n",Errors);
  return Errors;
}
