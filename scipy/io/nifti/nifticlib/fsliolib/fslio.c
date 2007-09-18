/*
    fslio.c  (Input and output routines for images in FSL)

    Mark Jenkinson
    FMRIB Image Analysis Group
*/


    
/*

    The fslio.c file was originally part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fslio.c has now been placed in the public domain.
   
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
*/

/** \file fslio.c 
    \brief Main collection of FSL i/o routines, written by Mark Jenkinson, FMRIB

    - updates by Rick Reynolds, SSCC, NIMH
 */

#include "fslio.h"
#include "assert.h"

static int FslIgnoreMFQ=0;
static int FslOverrideOutputType=-1;

#define FSLIOERR(x) { fprintf(stderr,"Error:: %s\n",(x)); fflush(stderr); exit(EXIT_FAILURE); }


/************************************************************
 * FslFileTypeString
 ************************************************************/
/*! \fn char* FslFileTypeString(int filetype)
    \brief  Return a string describing the format of the dataset
    \param filetype  FSL data format code.  Legal values are as defined
        by FSL_TYPE.
    \return  A string with the data format name, e.g. "ANALYZE-7.5"
    \sa FSL_TYPE
*/
char* FslFileTypeString(int filetype)
{
  if (filetype==FSL_TYPE_ANALYZE)          return "ANALYZE-7.5";
  if (filetype==FSL_TYPE_NIFTI)            return "NIFTI-1+";
  if (filetype==FSL_TYPE_NIFTI_PAIR)       return "NIFTI-1";
  if (filetype==FSL_TYPE_ANALYZE_GZ)       return "ANALYZE-7.5";
  if (filetype==FSL_TYPE_NIFTI_GZ)         return "NIFTI-1+";
  if (filetype==FSL_TYPE_NIFTI_PAIR_GZ)    return "NIFTI-1";
  return "UNKNOWN";
}


int FslIsValidFileType(int filetype)
{
  if ( (filetype!=FSL_TYPE_ANALYZE)    && (filetype!=FSL_TYPE_ANALYZE_GZ) && 
       (filetype!=FSL_TYPE_NIFTI)      && (filetype!=FSL_TYPE_NIFTI_GZ) && 
       (filetype!=FSL_TYPE_NIFTI_PAIR) && (filetype!=FSL_TYPE_NIFTI_PAIR_GZ) &&
       (filetype!=FSL_TYPE_MINC)       && (filetype!=FSL_TYPE_MINC_GZ) ) {
    fprintf(stderr,"Error: unrecognised file type: %d\n",filetype);
    return 0;
  }
  return 1;
}


int FslBaseFileType(int filetype)
{
  /* returns -1 to indicate error - unrecognised filetype */
  if ( (filetype==FSL_TYPE_ANALYZE_GZ)    || (filetype==FSL_TYPE_ANALYZE) )
    return FSL_TYPE_ANALYZE;
  if ( (filetype==FSL_TYPE_NIFTI_GZ)      || (filetype==FSL_TYPE_NIFTI) )
    return FSL_TYPE_NIFTI;
  if ( (filetype==FSL_TYPE_NIFTI_PAIR_GZ) || (filetype==FSL_TYPE_NIFTI_PAIR) )
    return FSL_TYPE_NIFTI_PAIR;
  if ( (filetype==FSL_TYPE_MINC_GZ)       || (filetype==FSL_TYPE_MINC) )  
    return FSL_TYPE_MINC;
  fprintf(stderr,"Error: unrecognised file type (%d)\n",filetype);
  return -1;
}


int FslGetFileType2(const FSLIO *fslio, int quiet)
{
  FSLIO *mutablefslio;
  if (fslio==NULL)  FSLIOERR("FslGetFileType: Null pointer passed for FSLIO");
  if ( (fslio->file_mode==FSL_TYPE_MINC) || (fslio->file_mode==FSL_TYPE_MINC_GZ) ) {
    return fslio->file_mode;
  }
  if ( !FslIsValidFileType(fslio->file_mode) )    return -1;

  if (fslio->niftiptr!=NULL) {   /* check that it is nifti_type and filetype are consistent */
    if (fslio->niftiptr->nifti_type != FslBaseFileType(fslio->file_mode)) {
      if (!quiet) {
        fprintf(stderr,"Warning: nifti structure and fsl structure disagree on file type\n");
        fprintf(stderr,"nifti = %d and fslio = %d\n",fslio->niftiptr->nifti_type,fslio->file_mode);
      }
      mutablefslio = (FSLIO *) fslio;  /* dodgy and will generate warnings */
      mutablefslio->niftiptr->nifti_type = FslBaseFileType(fslio->file_mode);
      return fslio->file_mode;
    }
 }
  return fslio->file_mode;
}

int FslGetFileType(const FSLIO *fslio)
{ 
  return FslGetFileType2(fslio,0);
}



void FslSetFileType(FSLIO *fslio, int filetype)
{
  if (fslio==NULL)  FSLIOERR("FslSetFileType: Null pointer passed for FSLIO");
  if ( (filetype==FSL_TYPE_MINC) || (filetype==FSL_TYPE_MINC_GZ) ) {
    fslio->file_mode = filetype;
    return;
  }
  if (! FslIsValidFileType(filetype)) { return; } 
  fslio->file_mode = filetype;  /* indicates general nifti - details in niftiptr */
  if (fslio->niftiptr!=NULL) { 
    fslio->niftiptr->nifti_type = FslBaseFileType(filetype); 
    nifti_set_iname_offset(fslio->niftiptr);
  }
}



int FslIsSingleFileType(int filetype)
{
  if ( (filetype==FSL_TYPE_NIFTI) || (filetype==FSL_TYPE_NIFTI_GZ) || 
       (filetype==FSL_TYPE_MINC)  || (filetype==FSL_TYPE_MINC_GZ) )
    return 1;
  return 0;
}


int FslIsCompressedFileType(int filetype)
{
  if ( filetype >=100 ) return 1;
  return 0;
}


int FslGetWriteMode(const FSLIO *fslio)
{
  if (fslio==NULL)  FSLIOERR("FslGetWriteMode: Null pointer passed for FSLIO");
  return fslio->write_mode;
} 


void FslSetWriteMode(FSLIO *fslio, int mode)
{
  if (fslio==NULL)  FSLIOERR("FslSetWriteMode: Null pointer passed for FSLIO");
  fslio->write_mode = mode;
}


int FslGetEnvOutputType(void)
{
  /* return type is one of FSL_TYPE_* or -1 to indicate error */
  char *otype;
  if (FslOverrideOutputType>=0)  return FslOverrideOutputType;
  otype = getenv("FSLOUTPUTTYPE");
  if (otype == NULL) {
    fprintf(stderr,"ERROR:: Environment variable FSLOUTPUTTYPE is not set!\n");
    fprintf(stderr,"Please make sure that the appropriate configuration file is sourced by your shell (e.g. by putting it in .profile).\n");
    fprintf(stderr,"e.g. bash or sh users add the line \". ${FSLDIR}/etc/fslconf/fsl.sh\"\n");
    fprintf(stderr,"e.g. tcsh or csh users add the line \"source ${FSLDIR}/etc/fslconf/fsl.csh\"\n");
    exit(EXIT_FAILURE);
  }
  if (strcmp(otype,"ANALYZE")==0) { return FSL_TYPE_ANALYZE; }
  if (strcmp(otype,"ANALYZE_GZ")==0) { return FSL_TYPE_ANALYZE_GZ; }
  if (strcmp(otype,"NIFTI")==0) { return FSL_TYPE_NIFTI; }
  if (strcmp(otype,"NIFTI_GZ")==0) { return FSL_TYPE_NIFTI_GZ; }
  if (strcmp(otype,"NIFTI_PAIR")==0) { return FSL_TYPE_NIFTI_PAIR; }
  if (strcmp(otype,"NIFTI_PAIR_GZ")==0) { return FSL_TYPE_NIFTI_PAIR_GZ; }
  if (strcmp(otype,"MINC")==0) { return FSL_TYPE_MINC; }
  if (strcmp(otype,"MINC_GZ")==0) { return FSL_TYPE_MINC_GZ; }
  fprintf(stderr,"ERROR:: Unrecognised value (%s) of environment variable FSLOUTPUT\n",otype);
  fprintf(stderr,"Legal values are: ANALYZE, NIFTI, NIFTI_PAIR, MINC, ANALYZE_GZ, NIFTI_GZ, NIFTI_PAIR_GZ, MINC_GZ\n");
  exit(EXIT_FAILURE);
  return -1;
}
    

int FslFileType(const char* fname) 
{
  /* return type is FSL_TYPE_* or -1 to indicate undetermined */
  /* use name as first priority but if that is ambiguous then resolve using environment */
  int flen;
  int retval=-1;
  if (fname==NULL) return retval;
  flen = strlen(fname);
  if (flen<5) return retval;  /* smallest name + extension is a.nii */
  if (strcmp(fname + flen - 4,".nii")==0)  retval=FSL_TYPE_NIFTI;
  if (strcmp(fname + flen - 7,".nii.gz")==0)  retval=FSL_TYPE_NIFTI_GZ;
  if (strcmp(fname + flen - 4,".mnc")==0)  retval=FSL_TYPE_MINC;
  if (strcmp(fname + flen - 7,".mnc.gz")==0)  retval=FSL_TYPE_MINC;
  if (strcmp(fname + flen - 4,".hdr")==0)  retval=FSL_TYPE_NIFTI_PAIR;
  if (strcmp(fname + flen - 4,".img")==0)  retval=FSL_TYPE_NIFTI_PAIR;
  if (strcmp(fname + flen - 7,".hdr.gz")==0)  retval=FSL_TYPE_NIFTI_PAIR_GZ;
  if (strcmp(fname + flen - 7,".img.gz")==0)  retval=FSL_TYPE_NIFTI_PAIR_GZ;
  if ( (retval==FSL_TYPE_NIFTI_PAIR) || (retval==FSL_TYPE_NIFTI_PAIR_GZ) ) {
    /* If it was hdr or img, check if Analyze was requested by environment */
    if ( (FslGetEnvOutputType() == FSL_TYPE_ANALYZE) && (retval == FSL_TYPE_NIFTI_PAIR) ) 
      retval=FSL_TYPE_ANALYZE;
    if ( (FslGetEnvOutputType() == FSL_TYPE_ANALYZE_GZ) && (retval == FSL_TYPE_NIFTI_PAIR_GZ) )
      retval=FSL_TYPE_ANALYZE_GZ;
  }
  return retval;
}


/************************************************************
 * FslGetReadFileType
 ************************************************************/
/*! \fn int FslGetReadFileType(const FSLIO *fslio)
    \brief  return the best estimate of the true file type 

  This function is used to return the best estimate of the true file type once
   a simple open has occurred - for now it is used after a nifti open call is made 

    \param  fslio data structure
    \return FSL_TYPE filetype code
    \sa FSL_TYPE
*/
int FslGetReadFileType(const FSLIO *fslio)
{
  int filetype=FSL_TYPE_ANALYZE;  /* unused default */
  if (fslio==NULL)  FSLIOERR("FslReadGetFileType: Null pointer passed for FSLIO");
  /* Don't use fslio->file_mode as it hasn't been set yet */
  if (fslio->niftiptr!=NULL) {   
    /* use the nifti_type and hdr or img name to determine the actual type */
    if (fslio->niftiptr->nifti_type == FSL_TYPE_ANALYZE) {
      if (FslIsCompressedFileType(FslFileType(fslio->niftiptr->iname))) {
        filetype = FSL_TYPE_ANALYZE_GZ;
      } else {
        filetype = FSL_TYPE_ANALYZE;
      }
    }
    if (fslio->niftiptr->nifti_type == FSL_TYPE_NIFTI_PAIR) {
      if (FslIsCompressedFileType(FslFileType(fslio->niftiptr->iname))) {
        filetype = FSL_TYPE_NIFTI_PAIR_GZ;
      } else {
        filetype = FSL_TYPE_NIFTI_PAIR;
      }
    }
    if (fslio->niftiptr->nifti_type == FSL_TYPE_NIFTI) {
      if (FslIsCompressedFileType(FslFileType(fslio->niftiptr->fname))) {
        filetype = FSL_TYPE_NIFTI_GZ;
      } else {
        filetype = FSL_TYPE_NIFTI;
      }
    }
    
  }
  if (fslio->mincptr!=NULL) { 
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
    filetype = FSL_TYPE_MINC;
  }
  return filetype;
}


int FslFileExists(const char *filename)
{ 
  /* return 1 if file(s) exists, otherwise return 0 */
  char *hdrname = nifti_findhdrname(filename);
  char *imgname = NULL;
  if (hdrname!=NULL){
      imgname = nifti_findimgname(filename,
                                  FslBaseFileType(FslFileType(hdrname)));
      free(hdrname);
      if (imgname != NULL) { free(imgname);  return 1; }
  }
  return 0;
}

char *FslMakeBaseName(const char *fname)
{
  char *basename;
  int blen;
  basename = nifti_makebasename(fname);
  blen = strlen(basename);
#ifdef HAVE_ZLIB
  if (strcmp(basename + blen-7,".mnc.gz") == 0) 
    { basename[blen-7]='\0'; return basename; }
#endif
  if (strcmp(basename + blen-4,".mnc") == 0) 
    { basename[blen-4]='\0'; return basename; }
  return basename;
}


void FslGetHdrImgNames(const char* filename, const FSLIO* fslio, 
                       char** hdrname, char** imgname)
{
  char *basename;
  int filetype;
  basename = FslMakeBaseName(filename);
  *hdrname = (char *)calloc(sizeof(char),strlen(basename)+8);
  *imgname = (char *)calloc(sizeof(char),strlen(basename)+8);
  strcpy(*hdrname,basename);
  strcpy(*imgname,basename);
  filetype = FslGetFileType(fslio);
  if (filetype==FSL_TYPE_NIFTI_GZ) {
    strcat(*hdrname,".nii.gz");
    strcat(*imgname,".nii.gz");
    free(basename);
    return;
  }
  if (filetype==FSL_TYPE_NIFTI) {
    strcat(*hdrname,".nii");
    strcat(*imgname,".nii");
    free(basename);
    return;
  }
  if (filetype==FSL_TYPE_MINC_GZ) {
    strcat(*hdrname,".mnc.gz");
    strcat(*imgname,".mnc.gz");
    free(basename);
    return;
  }
  if (filetype==FSL_TYPE_MINC) {
    strcat(*hdrname,".mnc");
    strcat(*imgname,".mnc");
    free(basename);
    return;
  }
  if ( (filetype==FSL_TYPE_NIFTI_PAIR_GZ) || (filetype==FSL_TYPE_ANALYZE_GZ) ) {
    strcat(*hdrname,".hdr.gz");
    strcat(*imgname,".img.gz");
    free(basename);
    return;
  }
  if ( (filetype==FSL_TYPE_NIFTI_PAIR) || (filetype==FSL_TYPE_ANALYZE) ) {
    strcat(*hdrname,".hdr");
    strcat(*imgname,".img");
    free(basename);
    return;
  }

  fprintf(stderr,"Error: Unrecognised filetype (%d)\n",FslGetFileType(fslio));
  free(basename);
  
  /* Failure */
  *hdrname = NULL;
  *imgname = NULL;
}



/***************************************************************
 * FslInit()
 ***************************************************************/
/*! \fn FSLIO *FslInit()
    \brief allocate space for the FSLIO struct and set some sensible defaults 
    \return  A pointer to an initialized FSLIO data structure
 */
FSLIO *FslInit(void)
{
  FSLIO *fslio;
  fslio = (FSLIO *) calloc(1,sizeof(FSLIO));
  FslSetInit(fslio);
  return fslio;
}

void FslSetInit(FSLIO* fslio)
{
  /* set some sensible defaults */
  fslio->niftiptr = NULL;
  fslio->mincptr  = NULL;
  FslSetFileType(fslio,FslGetEnvOutputType());
  FslSetWriteMode(fslio,0);
  fslio->written_hdr = 0;
}



void FslInit4Write(FSLIO* fslio, const char* filename, int ft)
{
  /* ft determines filetype if ft>=0*/ 
  int imgtype;

  FslSetWriteMode(fslio,1);

  /* Determine file type from image name (first priority) or environment (default) */
  imgtype = FslFileType(filename);
  if (imgtype<0)  imgtype = FslGetEnvOutputType();

  if (ft >= 0) imgtype = ft;

  if (!FslIsValidFileType(imgtype)) {
    fprintf(stderr,"Error: Failed to determine file type for writing in FslOpen()\n");
    exit(EXIT_FAILURE);
  }
  
  if ( (FslBaseFileType(imgtype)!=FSL_TYPE_MINC) ) {
    FslInitHeader(fslio, NIFTI_TYPE_FLOAT32,
                  1, 1, 1, 3,  0.0, 0.0, 0.0, 0.0,  4, "mm");
    
    FslSetFileType(fslio,imgtype);  /* this is after InitHeader as niftiptr set there */
    
    /* determine the header and image filename */
    FslGetHdrImgNames(filename,fslio,&(fslio->niftiptr->fname),&(fslio->niftiptr->iname));
    if ( (fslio->niftiptr->fname == NULL) || (fslio->niftiptr->iname == NULL) ) { 
      fprintf(stderr,"Error: cannot find filenames for %s\n",filename); 
    }

  } else if (FslBaseFileType(imgtype)==FSL_TYPE_MINC) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
    return;
  } else {
    fprintf(stderr,"Error:: unrecognised image type requested\n");
    return;
  }

  return;
}



void FslInitHeader(FSLIO *fslio, short t, 
                   size_t x, size_t y, size_t z, size_t v,
                   float vx, float vy, float vz, float tr,
                   size_t dim,
                   const char* units)
{
  /* NB: This function does not set the file type or write mode*/

  if (fslio==NULL)  FSLIOERR("FslInitHeader: Null pointer passed for FSLIO");
  
  fslio->niftiptr = nifti_simple_init_nim();
  /* make nifti type consistent with fslio */
  fslio->niftiptr->nifti_type = FslBaseFileType(fslio->file_mode);

  fslio->mincptr = NULL;

  FslSetDataType(fslio,t);
  FslSetDim(fslio,x,y,z,v);
  FslSetVoxDim(fslio,vx,vy,vz,tr);
  FslSetTimeUnits(fslio,"s");
  FslSetDimensionality(fslio,dim);
}


void FslCloneHeader(FSLIO *dest, const FSLIO *src)
{
  /* only clone the information that is stored in the disk version of the header */
  /*  - therefore _not_ the filenames, output type, write mode, etc */ 

  char *fname=NULL, *iname=NULL;
  void *data=NULL;
  int filetype, writemode;
  int preserve_nifti_values = 0;
  if (dest==NULL)  FSLIOERR("FslCloneHeader: Null pointer passed for FSLIO");
  if (src==NULL)   FSLIOERR("FslCloneHeader: Null pointer passed for FSLIO");

  if (src->niftiptr!=NULL) {
    /* preserve the filenames, output type and write mode */
    if (dest->niftiptr != NULL) {
      fname = dest->niftiptr->fname;
      iname = dest->niftiptr->iname;
      data = dest->niftiptr->data;
      preserve_nifti_values = 1;
    }
    filetype = FslGetFileType2(dest,1);
    writemode = FslGetWriteMode(dest);

    /* copy _all_ info across */
    dest->niftiptr = nifti_copy_nim_info(src->niftiptr);

    /* restore old values */
    if (preserve_nifti_values) {
      dest->niftiptr->data = data; 
    } else { 
        /* destroy the values that the nifti copy creates */
      free(dest->niftiptr->fname);
      free(dest->niftiptr->iname);
      nifti_free_extensions(dest->niftiptr);

      dest->niftiptr->fname = NULL;
      dest->niftiptr->iname = NULL; 
      dest->niftiptr->data = NULL;   /* should already be NULL */
    }
    FslSetFileType(dest,filetype);
    FslSetWriteMode(dest,writemode);
  }

  if (src->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }

}


int  fsl_fileexists(const char* fname)
{
   znzFile fp;
   fp = znzopen( fname , "rb" , 1 ) ;
   if( !znz_isnull(fp) )  { znzclose(fp);  return 1; }
   return 0;
}


int FslCheckForMultipleFileNames(const char* filename)
{
  char *basename, *tmpname;
  int singlecount=0, hdrcount=0, imgcount=0, ambiguous=0;
  basename =  nifti_makebasename(filename);
  tmpname = (char *)calloc(strlen(basename) + 10,sizeof(char));

  strcpy(tmpname,basename);
  strcat(tmpname,".nii"); 
  if (fsl_fileexists(tmpname)) { singlecount++; }
  strcpy(tmpname,basename);
  strcat(tmpname,".nii.gz"); 
  if (fsl_fileexists(tmpname)) { singlecount++; }
  strcpy(tmpname,basename);
  strcat(tmpname,".mnc"); 
  if (fsl_fileexists(tmpname)) { singlecount++; }
  strcpy(tmpname,basename);
  strcat(tmpname,".mnc.gz"); 
  if (fsl_fileexists(tmpname)) { singlecount++; }

  strcpy(tmpname,basename);
  strcat(tmpname,".img"); 
  if (fsl_fileexists(tmpname)) { imgcount++; }
  strcpy(tmpname,basename);
  strcat(tmpname,".img.gz"); 
  if (fsl_fileexists(tmpname)) { imgcount++; }

  strcpy(tmpname,basename);
  strcat(tmpname,".hdr"); 
  if (fsl_fileexists(tmpname)) { hdrcount++; }
  strcpy(tmpname,basename);
  strcat(tmpname,".hdr.gz"); 
  if (fsl_fileexists(tmpname)) { hdrcount++; }
 
  ambiguous = 1;
  if ( (hdrcount==1) && (imgcount==1) && (singlecount==0) )  { ambiguous=0; }
  if ( (hdrcount==0) && (imgcount==0) && (singlecount==1) )  { ambiguous=0; }

  /* treat no image found as not ambiguous - want opening errors instead */
  if ( (hdrcount==0) && (imgcount==0) && (singlecount==0) )  { ambiguous=0; }

  free(tmpname);
  free(basename);
  return ambiguous;
}



int check_for_multiple_filenames(const char* filename)
{
  char *basename, *tmpname;
  char *otype;
  if (FslCheckForMultipleFileNames(filename))
    {  /* take action */
      basename =  nifti_makebasename(filename);
      tmpname = (char *)calloc(strlen(basename) + 10,sizeof(char));
      fprintf(stderr,"\n\n\nWARNING!!!! Multiple image files detected:\n");
      /* list the offending files */
      strcpy(tmpname,basename);
      strcat(tmpname,".nii"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".nii.gz"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".mnc"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".mnc.gz"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".img"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".img.gz"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".hdr"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".hdr.gz"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      fprintf(stderr,"\n\n");

      if (!FslIgnoreMFQ) {
        otype = getenv("FSLMULTIFILEQUIT");
        if (otype!=NULL) {
          fprintf(stderr,"STOPPING PROGRAM\n");
          exit(EXIT_FAILURE);
        }
      }
      return 1;
    }
  return 0;
}



/***************************************************************
 * FslOpen
 ***************************************************************/
/*! \fn FSLIO *FslOpen(const char *filename, const char *opts)
    \brief Opens a file for either reading or writing. 

        The format of the output dataset is determined automatically by 
        passing filetype -1 to FslXOpen.
    \sa FslXOpen
 */
FSLIO *FslOpen(const char *filename, const char *opts)
{
  /* Note: -1 for filetype indicates that FslXOpen should determine filetype for itself */
  return FslXOpen(filename,opts,-1);
}


/***************************************************************
 * FslXOpen
 ***************************************************************/
/*! \fn FSLIO *FslXOpen(const char *filename, const char *opts, int filetype)
    \brief Opens a file for either reading or writing

        Files to be read are automatically read whether 
        compressed or not.  Also, reading uses the file extension 
        and will fail if that file does not exist.
        For a more robust read, pass the basename in as then all types 
        will be tried.
    \param filename Name (or basename) of the file to open
    \param opts Flags for fopen() of dataset, eg "r", "wb", etc.
    \param filetype specifies the type of file to be written. Legal
        values are as defined by FSL_TYPE.  If filetype is less than 
        zero, then it is ignored and the type is determined by the 
        filename extension or, failing that, the environment default.
    \return pointer to FSLIO dataset datastructure
    \sa FSLIO
    \sa FSL_TYPE
 */
FSLIO *FslXOpen(const char *filename, const char *opts, int filetype)
{

  FSLIO *fslio;
  char bopts[1024];
  size_t i, bi;
  int imgtype;

  fslio = FslInit();

  bi=0;
  for(i=0;i<strlen(opts);i++) {
    if (opts[i]=='w') { FslSetWriteMode(fslio,1); }
    if (opts[i]!='b' && opts[i]!='t') { bopts[bi++]=opts[i]; }
  }
  /* add in 'b' (at the end) for windows compatibility */
  bopts[bi++]='b';
  bopts[bi]='\0';
  

  if (FslGetWriteMode(fslio)==1) {
    
    /** ====================== Open file for writing ====================== **/
   
    FslInit4Write(fslio,filename,filetype);
    imgtype = FslGetFileType(fslio);
    fslio->written_hdr = 0;

    /* open the image file - not the header */
    fslio->fileptr = znzopen(fslio->niftiptr->iname,bopts,FslIsCompressedFileType(imgtype));
    if (znz_isnull(fslio->fileptr)) { 
      fprintf(stderr,"Error: failed to open file %s\n",fslio->niftiptr->iname); 
      return NULL;
    }

    if (!FslIsSingleFileType(imgtype)) {
      /* set up pointer at end of iname_offset for dual file formats (not singles) */
      FslSeekVolume(fslio,0);
    }
    return fslio;

  }



  /** ======================== Open file for reading ====================== **/

  check_for_multiple_filenames(filename);

  /* see if the extension indicates a minc file */
  imgtype = FslFileType(filename);
  if ((imgtype>=0) && (FslBaseFileType(imgtype)==FSL_TYPE_MINC)) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
    return NULL;
  }

  /* otherwise open nifti file: read header and open img file (may be same file) */
  fslio->fileptr = nifti_image_open(filename,bopts,&(fslio->niftiptr));
  if (znz_isnull(fslio->fileptr)) { 
    fprintf(stderr,"Error: failed to open file %s\n",filename); 
    return NULL;
  }

  /* set the file type given what has been read - it uses nifti_type and filenames */
  imgtype = FslGetReadFileType(fslio);
  FslSetFileType(fslio,imgtype);
  FslSetWriteMode(fslio,0);

  if (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_ANALYZE) {
    /* For the ANALYZE case in FSL, must cheat and grab the originator field! */
    /* Note that the header file is always separate here and closed by now */
    struct dsr ahdr;
    short orig[5];
    FslReadRawHeader(&ahdr,fslio->niftiptr->fname);
    if (fslio->niftiptr->byteorder != nifti_short_order()) {
      AvwSwapHeader(&ahdr);
    }
    /* Read the origin and set the sform up (if origin is non-zero) */
    /* Note that signed pixdims are passed in to set the LR orientation */
    memcpy(orig,&(ahdr.hist.originator),10);
    FslSetAnalyzeSform(fslio, orig, fslio->niftiptr->pixdim[1],
                       fslio->niftiptr->pixdim[2], fslio->niftiptr->pixdim[3]);
  }

  /* from now on force all vox dims to be positive - LR info is in sform */
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->dx = fabs(fslio->niftiptr->dx);
    fslio->niftiptr->dy = fabs(fslio->niftiptr->dy);
    fslio->niftiptr->dz = fabs(fslio->niftiptr->dz);
    fslio->niftiptr->pixdim[1] = fabs(fslio->niftiptr->pixdim[1]);
    fslio->niftiptr->pixdim[2] = fabs(fslio->niftiptr->pixdim[2]);
    fslio->niftiptr->pixdim[3] = fabs(fslio->niftiptr->pixdim[3]);
  }
  /* set up pointer at end of iname_offset , ready for reading */
  FslSeekVolume(fslio,0);  

  return fslio;

}



/***************************************************************
 * FslReadAllVolumes
 ***************************************************************/
/*! \fn void* FslReadAllVolumes(FSLIO* fslio, char* filename)
    \brief Read the header and all data into the FSLIO structure

        There is no need for FslOpen or FslClose calls when FslReadAllVolumes()
        is called.  
        <br>This routine allocates the buffer to hold the entire dataset. 
        <br>The data block returned will contain the data in whatever
        datatype it is stored as on disk (therefore it is a void *).
        <br>The data buffer will be byteswapped to native-endian.
        <br>The data buffer will not be scaled. 
        <br>The best call to make before this is FslInit() or a calloc() for 
        fslio.  (??? why calloc if this allocates the buffer ???)

    \param fslio pointer to an open dataset
    \param filename Name of the dataset to read.
    \return A pointer to the data block buffer (allocated by this function).
        <br> Return Null on error ??? is this true ???
        <ul>
        <li>Note this pointer is also in the FSLIO structure as 
        fslio->niftiptr->data.</li>
        <li>Note a void pointer is returned, as the datablock is of
        variable datatype.</li>
        </ul>

 */
void* FslReadAllVolumes(FSLIO* fslio, char* filename)
{

  int imgtype;
  if (fslio==NULL)  FSLIOERR("FslReadAllVolumes: Null pointer passed for FSLIO");

  /* see if the extension indicates a minc file */
  imgtype = FslFileType(filename);
  if ((imgtype>=0) && (FslBaseFileType(imgtype)==FSL_TYPE_MINC)) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
    return NULL;
  }

  /** otherwise it is a nifti file - so read it! **/
  fslio->mincptr = NULL;
  /* make sure an FslOpen hasn't locked the file */
  if (!znz_isnull(fslio->fileptr)) FslClose(fslio);  
  
  fslio->niftiptr = nifti_image_read(filename,1);

  /* check for failure, from David Akers */
  if (fslio->niftiptr == NULL) {
        FSLIOERR("FslReadAllVolumes: error reading NIfTI image");
        return(NULL);
  }

  FslSetFileType(fslio,fslio->niftiptr->nifti_type);
  FslSetWriteMode(fslio,0);
  return fslio->niftiptr->data;
}



/***************************************************************
 * FslReadVolumes
 ***************************************************************/
/*! \fn size_t FslReadVolumes(FSLIO *fslio, void *buffer, size_t nvols)
    \brief Read the first nvols Volumes from a 4D dataset

    \param fslio pointer to open dataset
    \param buffer buffer to read data into, allocated by ???
    \param nvols  number of volumes to read
    \return Number of volumes read.
 */
size_t FslReadVolumes(FSLIO *fslio, void *buffer, size_t nvols)
{
  int volbytes;
  size_t retval=0;
  if (fslio==NULL)  FSLIOERR("FslReadVolumes: Null pointer passed for FSLIO");
  if (znz_isnull(fslio->fileptr))  FSLIOERR("FslReadVolumes: Null file pointer");
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->data = buffer;
    volbytes = FslGetVolSize(fslio)  * fslio->niftiptr->nbyper;
    retval = nifti_read_buffer(fslio->fileptr,fslio->niftiptr->data,nvols*volbytes,fslio->niftiptr);
    retval /= volbytes;
  }

  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return retval;
}



/***************************************************************
 * FslWriteAllVolumes
 ***************************************************************/
/*! \fn void FslWriteAllVolumes(FSLIO *fslio, const void *buffer)
    \brief  Writes all data from buffer (using size info from fslio) to file.
        
        Dimension and datatype of buffer are as is specified in nifti_image structure
        fslio->niftiptr.
        Note: If file format is Analyze (not nifti) and in Neurological order then 
        SWAP DATA into Radiological order.

    \param fslio pointer to open dataset
    \param buffer pointer to data array. Size and datatype of this buffer  
 */
void FslWriteAllVolumes(FSLIO *fslio, const void *buffer)
{
  short x,y,z,t;

  if (fslio==NULL)  FSLIOERR("FslWriteAllVolumes: Null pointer passed for FSLIO");

  FslGetDim(fslio,&x,&y,&z,&t);
  FslWriteHeader(fslio);
  FslWriteVolumes(fslio,buffer,t);
  return;
}



/***************************************************************
 * FslWriteVolumes
 ***************************************************************/
/*! \fn size_t FslWriteVolumes(FSLIO *fslio, const void *buffer, size_t nvols)
    \brief Write the first nvols volumes in buffer to disk.  

        Dimension and datatype of buffer are as is specified in nifti_image structure
        fslio->niftiptr.
        Note: If file format is Analyze (not nifti) and in Neurological order then 
        SWAP DATA into Radiological order.

        
    \param fslio        pointer to open dataset
    \param buffer       pointer to data array. Size and datatype of this buffer  
    \param nvols        number of volumes to write
    \return ??? looks like return of retval is missing ???  0 on error.
 */
size_t FslWriteVolumes(FSLIO *fslio, const void *buffer, size_t nvols)
{
  /* The dimensions and datatype must be set before calling this function */
  int retval;
  if (fslio==NULL)  FSLIOERR("FslWriteVolumes: Null pointer passed for FSLIO");
  if ( (!fslio->written_hdr) && (FslIsSingleFileType(FslGetFileType(fslio))) &&
       (FslIsCompressedFileType(FslGetFileType(fslio))) )
    { FSLIOERR("FslWriteVolumes: header must be written before data for single compressed file types"); }
  
  if (fslio->niftiptr!=NULL) {
    long int nbytes, bpv;
    bpv = fslio->niftiptr->nbyper;  /* bytes per voxel */
    nbytes = nvols * FslGetVolSize(fslio) * bpv;

    if ( (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_ANALYZE)
         && (FslGetLeftRightOrder(fslio)==FSL_NEUROLOGICAL) ) {
      /* If it is Analyze and Neurological order then SWAP DATA into Radiological order */
      /* This is nasty - but what else can be done?!? */
      char *tmpbuf, *inbuf;
      long int x, b, n, nrows;
      short nx, ny, nz, nv;
      inbuf = (char *) buffer;
      tmpbuf = (char *)calloc(nbytes,1);
      FslGetDim(fslio,&nx,&ny,&nz,&nv);
      nrows = nbytes / (nx * bpv);
      for (n=0; n<nrows; n++) {
        for (x=0; x<nx; x++) {
          for (b=0; b<bpv; b++) {
            tmpbuf[b +  bpv * (n*nx + nx - 1 - x)] = inbuf[b + bpv * (n*nx + x)];
          }
        }
      }
      retval = nifti_write_buffer(fslio->fileptr, tmpbuf, nbytes);
      free(tmpbuf);
    } else {
      retval = nifti_write_buffer(fslio->fileptr, buffer, nbytes);
    }
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;  /* failure */
}


/***************************************************************
 * FslWriteHeader
 ***************************************************************/
/*! \fn void FslWriteHeader(FSLIO *fslio)
    \brief Writes nifti/anz header and opens img file ready for writing
        
    \param fslio        pointer to open dataset
 */
void FslWriteHeader(FSLIO *fslio)
{
  /* writes header and opens img file ready for writing */
  if (fslio==NULL)  FSLIOERR("FslWriteHeader: Null pointer passed for FSLIO");

  if (fslio->niftiptr!=NULL) {
    fslio->written_hdr = 1;
    if (znz_isnull(fslio->fileptr)) FSLIOERR("FslWriteHeader: no file opened!");

    strcpy(fslio->niftiptr->descrip,"FSL3.2beta");
    if (FslIsSingleFileType(FslGetFileType(fslio))) {
      /* write header info but don't close the file */
      nifti_image_write_hdr_img2(fslio->niftiptr,2,"wb",fslio->fileptr,NULL);
      /* set up pointer at end of iname_offset for single files only */
      FslSeekVolume(fslio,0);
    } else {
      /* open a new hdr file, write it and close it */
      nifti_image_write_hdr_img(fslio->niftiptr,0,"wb");
    }
  }

  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return;
}


/***************************************************************
 * FslReadSliceSeries
 ***************************************************************/
/*! \fn size_t FslReadSliceSeries(FSLIO *fslio, void *buffer, short slice, size_t nvols)
    \brief Read one slice from each of the first nvols volumes in the dataset, ie get an xyt buffer.

        Dimension and datatype of buffer are as is specified in nifti_image structure
        fslio->niftiptr.
        Note: filepointer in file data array is restored to its initial position.
        
    \param fslio        pointer to open dataset
    \param buffer       buffer large enough to hold 1 slice from each volume
    \param slice        slice number (0 based) to read  [0 z-1]
    \param nvols        number of volumes to read a slice from 
    \return             Number of volumes from which a slice was successfully read. 0 on error.
 */
size_t FslReadSliceSeries(FSLIO *fslio, void *buffer, short slice, size_t nvols)
{
  size_t slbytes,volbytes;
  size_t n, orig_offset;
  short x,y,z,v,type;

  if (fslio==NULL)  FSLIOERR("FslReadSliceSeries: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    
    FslGetDim(fslio,&x,&y,&z,&v);
    
    if ((slice<0) || (slice>=z)) FSLIOERR("FslReadSliceSeries: slice outside valid range");
    
    slbytes = x * y * (FslGetDataType(fslio, &type) / 8);
    volbytes = slbytes * z;
    
    orig_offset = znztell(fslio->fileptr);
    znzseek(fslio->fileptr, slbytes*slice, SEEK_CUR);
    
    for (n=0; n<nvols; n++) {
      if (n>0) znzseek(fslio->fileptr, volbytes - slbytes, SEEK_CUR);
      if (znzread((char *)buffer+n*slbytes, 1, slbytes, fslio->fileptr) != slbytes)
        FSLIOERR("FslReadSliceSeries: failed to read values");
     if (fslio->niftiptr->byteorder != nifti_short_order())
        nifti_swap_Nbytes(slbytes / fslio->niftiptr->swapsize,
                          fslio->niftiptr->swapsize, (char *)buffer+n*slbytes);
     }
    
    
    /* restore file pointer to original position */
    znzseek(fslio->fileptr,orig_offset,SEEK_SET);
    return n;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
}


/***************************************************************
 * FslReadRowSeries
 ***************************************************************/
/*! \fn size_t FslReadRowSeries(FSLIO *fslio, void *buffer, short row, short slice, size_t nvols)
    \brief Read one row from one slice for first nvols volumes in dataset; ie get an xt buffer.

        Dimension and datatype of buffer are as is specified in nifti_image structure
        fslio->niftiptr.
        Note: filepointer in file data array is restored to its initial position.
        
    \param fslio        pointer to open dataset
    \param buffer       buffer to hold one row from each volume.
    \param row          row number (0 based) to read [0 y-1]
    \param slice        slice number (0 based) to read
    \param nvols        number of volumes to read a row from 
    \return Number of volumes from which a row was successfully read. 0 on error.
 */
size_t FslReadRowSeries(FSLIO *fslio, void *buffer, short row, short slice, size_t nvols)
{
  size_t rowbytes,slbytes,volbytes;
  size_t n, orig_offset;
  short x,y,z,v,type;
  
  if (fslio==NULL)  FSLIOERR("FslReadRowSeries: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    
    FslGetDim(fslio,&x,&y,&z,&v);
    
    if ((slice<0) || (slice>=z)) FSLIOERR("FslReadRowSeries: slice outside valid range");
    if ((row<0) || (row>=y)) FSLIOERR("FslReadRowSeries: row outside valid range");
    
    rowbytes = x * (FslGetDataType(fslio, &type)) / 8;
    slbytes = rowbytes * y;
    volbytes = slbytes * z;
    
    orig_offset = znztell(fslio->fileptr);
    znzseek(fslio->fileptr, rowbytes*row + slbytes*slice, SEEK_CUR);
    
    for (n=0; n<nvols; n++){
      if (n>0) znzseek(fslio->fileptr, volbytes - rowbytes, SEEK_CUR);
      if (znzread((char *)buffer+n*rowbytes, 1, rowbytes, fslio->fileptr) != rowbytes)
        FSLIOERR("FslReadRowSeries: failed to read values");
      if (fslio->niftiptr->byteorder != nifti_short_order())
        nifti_swap_Nbytes(rowbytes / fslio->niftiptr->swapsize,
                          fslio->niftiptr->swapsize, (char *)buffer+n*rowbytes);
    }
    
    /* restore file pointer to original position */
    znzseek(fslio->fileptr,orig_offset,SEEK_SET);
    return n;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
}


/***************************************************************
 * FslReadTimeSeries
 ***************************************************************/
/*! \fn size_t FslReadTimeSeries(FSLIO *fslio, void *buffer, short xVox, short yVox, short zVox, size_t nvols)
    \brief Read one voxel (xyz location) from first nvols volumes in dataset; ie get a t  dim buffer.

        Dimension and datatype of buffer are as is specified in nifti_image structure
        fslio->niftiptr.
        Note: filepointer in file data array is restored to its initial position.
        
    \param fslio        pointer to open dataset
    \param buffer       buffer to hold one timeseries vector
    \param xVox         x voxel [0 x-1]
    \param yVox         y voxel [0 y-1]
    \param zVox         z voxel [0 z-1]
    \param nvols        number of volumes to read a voxel from
    \return Number of volumes from which a voxel was successfully read. 0 on error.
 */
size_t FslReadTimeSeries(FSLIO *fslio, void *buffer, short xVox, short yVox, short zVox, 
                         size_t nvols)
{
  size_t volbytes, offset, orig_offset;
  size_t n;
  short xdim,ydim,zdim,v,wordsize;

  if (fslio==NULL)  FSLIOERR("FslReadTimeSeries: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {

    FslGetDim(fslio,&xdim,&ydim,&zdim,&v);
    
    if ((xVox<0) || (xVox >=xdim)) FSLIOERR("FslReadTimeSeries: voxel outside valid range");
    if ((yVox<0) || (yVox >=ydim)) FSLIOERR("FslReadTimeSeries: voxel outside valid range");
    if ((zVox<0) || (zVox >=zdim)) FSLIOERR("FslReadTimeSeries: voxel outside valid range");
    
    wordsize = fslio->niftiptr->nbyper;
    volbytes = xdim * ydim * zdim * wordsize;
    
    orig_offset = znztell(fslio->fileptr);
    offset = ((ydim * zVox + yVox) * xdim + xVox) * wordsize;
    znzseek(fslio->fileptr,offset,SEEK_CUR);
    
    for (n=0; n<nvols; n++) {
      if (n>0) znzseek(fslio->fileptr, volbytes - wordsize, SEEK_CUR);
      if (znzread((char *)buffer+(n*wordsize), 1, wordsize,fslio->fileptr) != wordsize)
        FSLIOERR("FslReadTimeSeries: failed to read values"); 
      if (fslio->niftiptr->byteorder != nifti_short_order())
        nifti_swap_Nbytes(1,fslio->niftiptr->swapsize,
                          (char *)buffer+(n*wordsize));
    }
    
    /* restore file pointer to original position */
    znzseek(fslio->fileptr,orig_offset,SEEK_SET);
    return n;

  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
}


size_t FslReadCplxVolumes(FSLIO *fslio, void *buffer, size_t nvols, char mode)
{
  if (fslio==NULL)  FSLIOERR("FslReadCplxVolumes: Null pointer passed for FSLIO");
  fprintf(stderr,"Warning:: FslReadCplxVolumes is not yet supported\n");
  return 0;
}

size_t FslWriteCplxVolumes(FSLIO *fslio, void *buffer, size_t nvols, char mode)
{
  if (fslio==NULL)  FSLIOERR("FslWriteCplxVolumes: Null pointer passed for FSLIO");
  fprintf(stderr,"Warning:: FslWriteCplxVolumes is not yet supported\n");
  return 0;
}

int FslSeekVolume(FSLIO *fslio, size_t vols)
{
  int offset;
  if (fslio==NULL)  FSLIOERR("FslSeekVolume: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    offset = fslio->niftiptr->iname_offset + 
      vols * FslGetVolSize(fslio) * fslio->niftiptr->nbyper;
    if (znz_isnull(fslio->fileptr)) FSLIOERR("FslSeekVolume: Null file pointer");
    return znzseek(fslio->fileptr,offset,SEEK_SET);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
}


size_t FslGetVolSize(FSLIO *fslio)
{
  /* returns number of voxels per 3D volume */
  if (fslio==NULL)  FSLIOERR("FslGetVolSize: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    return (fslio->niftiptr->nx * fslio->niftiptr->ny * fslio->niftiptr->nz);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
}


void FslSetDim(FSLIO *fslio, short x, short y, short z, short v)
{
  int ndim;
  if (fslio==NULL)  FSLIOERR("FslSetDim: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {

    ndim=4;
    if (v<=1) {ndim--; if (z<=1) {ndim--; if (y<=1) {ndim--; if (x<=1) {ndim--;}}}}

    fslio->niftiptr->ndim = ndim;

    if (x>=1) fslio->niftiptr->nx = x; else fslio->niftiptr->nx=1;
    if (y>=1) fslio->niftiptr->ny = y; else fslio->niftiptr->ny=1;
    if (z>=1) fslio->niftiptr->nz = z; else fslio->niftiptr->nz=1;
    if (v>=1) fslio->niftiptr->nt = v; else fslio->niftiptr->nt=1;
    fslio->niftiptr->nu = 1;
    fslio->niftiptr->nv = 1;
    fslio->niftiptr->nw = 1;

    /* deal with stupid redundancies */
    fslio->niftiptr->dim[0] = fslio->niftiptr->ndim ;
    fslio->niftiptr->dim[1] = fslio->niftiptr->nx;
    fslio->niftiptr->dim[2] = fslio->niftiptr->ny;
    fslio->niftiptr->dim[3] = fslio->niftiptr->nz;
    fslio->niftiptr->dim[4] = fslio->niftiptr->nt;
    fslio->niftiptr->dim[5] = fslio->niftiptr->nu;
    fslio->niftiptr->dim[6] = fslio->niftiptr->nv;
    fslio->niftiptr->dim[7] = fslio->niftiptr->nw;

    fslio->niftiptr->nvox =  fslio->niftiptr->nx * fslio->niftiptr->ny * fslio->niftiptr->nz
      * fslio->niftiptr->nt * fslio->niftiptr->nu * fslio->niftiptr->nv * fslio->niftiptr->nw ;

  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetDim(FSLIO *fslio, short *x, short *y, short *z, short *v)
{
  if (fslio==NULL)  FSLIOERR("FslGetDim: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *x = fslio->niftiptr->nx;
    *y = fslio->niftiptr->ny;
    *z = fslio->niftiptr->nz;
    *v = fslio->niftiptr->nt;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetDimensionality(FSLIO *fslio, size_t dim)
{
  if (fslio==NULL)  FSLIOERR("FslSetDimensionality: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->ndim = dim;
    fslio->niftiptr->dim[0] = dim;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetDimensionality(FSLIO *fslio, size_t *dim)
{
  if (fslio==NULL)  FSLIOERR("FslGetDimensionality: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *dim = fslio->niftiptr->ndim;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetVoxDim(FSLIO *fslio, float x, float y, float z, float tr)
{
  if (fslio==NULL)  FSLIOERR("FslSetVoxDim: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->dx = fabs(x);
    fslio->niftiptr->dy = fabs(y);
    fslio->niftiptr->dz = fabs(z);
    fslio->niftiptr->dt = fabs(tr);
    fslio->niftiptr->pixdim[1] = fabs(x);
    fslio->niftiptr->pixdim[2] = fabs(y);
    fslio->niftiptr->pixdim[3] = fabs(z);
    fslio->niftiptr->pixdim[4] = fabs(tr);
    /* set the units to mm and seconds */
    fslio->niftiptr->xyz_units  = NIFTI_UNITS_MM;
    fslio->niftiptr->time_units = NIFTI_UNITS_SEC;
 }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetVoxDim(FSLIO *fslio, float *x, float *y, float *z, float *tr)
{
  if (fslio==NULL)  FSLIOERR("FslGetVoxDim: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *x = fabs(fslio->niftiptr->dx);
    *y = fabs(fslio->niftiptr->dy);
    *z = fabs(fslio->niftiptr->dz);
    *tr = fabs(fslio->niftiptr->dt);
    /* now check the units and convert to mm and sec */
    if (fslio->niftiptr->xyz_units == NIFTI_UNITS_METER) 
    { *x *= 1000.0;   *y *= 1000.0;   *z *= 1000.0; }
    if (fslio->niftiptr->xyz_units == NIFTI_UNITS_MICRON) 
    { *x /= 1000.0;   *y /= 1000.0;   *z /= 1000.0; }
    if (fslio->niftiptr->xyz_units == NIFTI_UNITS_MSEC) 
    { *tr /= 1000.0; }
    if (fslio->niftiptr->xyz_units == NIFTI_UNITS_USEC) 
    { *tr /= 1000000.0; }
    /* if it is Hz or other frequency then leave it */
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetCalMinMax(FSLIO *fslio, float *min, float *max)
{
  if (fslio==NULL)  FSLIOERR("FslGetCalMinMax: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *min = fslio->niftiptr->cal_min;
    *max = fslio->niftiptr->cal_max;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetCalMinMax(FSLIO *fslio, float  min, float  max)
{
  if (fslio==NULL)  FSLIOERR("FslSetCalMinMax: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->cal_min = min;
    fslio->niftiptr->cal_max = max;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetAuxFile(FSLIO *fslio,char *aux_file)
{
  if (fslio==NULL)  FSLIOERR("FslGetAuxFile: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    strncpy(aux_file,fslio->niftiptr->aux_file, 24);
    aux_file[23] = '\0';
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetAuxFile(FSLIO *fslio,const char *aux_file)
{
  if (fslio==NULL)  FSLIOERR("FslSetAuxFile: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    strncpy(fslio->niftiptr->aux_file, aux_file, 24);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetVoxUnits(FSLIO *fslio, const char *units)
{
  int unitcode=0;
  if (fslio==NULL)  FSLIOERR("FslSetVoxUnits: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    if (strcmp(units,nifti_units_string(NIFTI_UNITS_METER))==0) {
      unitcode = NIFTI_UNITS_METER;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_MM))==0) {
      unitcode = NIFTI_UNITS_MM;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_MICRON))==0) {
      unitcode = NIFTI_UNITS_MICRON;
    }
    fslio->niftiptr->xyz_units  = unitcode;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetVoxUnits(FSLIO *fslio, char *units)
{
  if (fslio==NULL)  FSLIOERR("FslGetVoxUnits: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    strcpy(units,nifti_units_string(fslio->niftiptr->xyz_units));
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}

void FslSetTimeUnits(FSLIO *fslio, const char *units)
{
  int unitcode=0;
  if (fslio==NULL)  FSLIOERR("FslSetTimeUnits: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    if (strcmp(units,nifti_units_string(NIFTI_UNITS_HZ))==0) {
      unitcode = NIFTI_UNITS_HZ;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_PPM))==0) {
      unitcode = NIFTI_UNITS_PPM;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_RADS))==0) {
      unitcode = NIFTI_UNITS_RADS;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_SEC))==0) {
      unitcode = NIFTI_UNITS_SEC;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_MSEC))==0) {
        fprintf(stderr,"Warning::Setting time units to msec is not fully recommended in fslio\n");
        unitcode = NIFTI_UNITS_MSEC;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_USEC))==0) {
        fprintf(stderr,"Warning::Setting time units to msec is not fully recommended in fslio\n");
        unitcode = NIFTI_UNITS_USEC;
    }
    fslio->niftiptr->time_units = unitcode;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetTimeUnits(FSLIO *fslio, char *units)
{
  if (fslio==NULL)  FSLIOERR("FslGetTimeUnits: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    strcpy(units,nifti_units_string(fslio->niftiptr->time_units));
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetDataType(FSLIO *fslio, short t)
{
  int nbytepix=0, ss=0;
  if (fslio==NULL)  FSLIOERR("FslSetDataType: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->datatype = t;
    nifti_datatype_sizes(t,&nbytepix,&ss);
    fslio->niftiptr->nbyper = nbytepix;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}

size_t FslGetDataType(FSLIO *fslio, short *t)
{
    /* returns bits per pixel */
  int nbytepix=32, ss=0;
  if (fslio==NULL)  FSLIOERR("FslGetDataType: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *t = fslio->niftiptr->datatype;
    nifti_datatype_sizes(*t,&nbytepix,&ss);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return (size_t) 8 * nbytepix;
}


void FslGetMMCoord(mat44 stdmat, float voxx, float voxy, float voxz, 
                   float *mmx, float *mmy, float *mmz) 
{
    *mmx = stdmat.m[0][0] * voxx + stdmat.m[0][1] * voxy + stdmat.m[0][2] * voxz 
        + stdmat.m[0][3];
    *mmy = stdmat.m[1][0] * voxx + stdmat.m[1][1] * voxy + stdmat.m[1][2] * voxz 
        + stdmat.m[1][3];
    *mmz = stdmat.m[2][0] * voxx + stdmat.m[2][1] * voxy + stdmat.m[2][2] * voxz 
        + stdmat.m[2][3];
}


void FslGetVoxCoord(mat44 stdmat, float mmx, float mmy, float mmz, 
                   float *voxx, float *voxy, float *voxz) 
{
  mat44 mm2vox;

  mm2vox = nifti_mat44_inverse(stdmat);
    *voxx = mm2vox.m[0][0] * mmx + mm2vox.m[0][1] * mmy + mm2vox.m[0][2] * mmz 
        + mm2vox.m[0][3];
    *voxy = mm2vox.m[1][0] * mmx + mm2vox.m[1][1] * mmy + mm2vox.m[1][2] * mmz 
        + mm2vox.m[1][3];
    *voxz = mm2vox.m[2][0] * mmx + mm2vox.m[2][1] * mmy + mm2vox.m[2][2] * mmz 
        + mm2vox.m[2][3];
}


void FslSetStdXform(FSLIO *fslio, short sform_code, mat44 stdmat)
{
    /* NB: stdmat must point to a 4x4 array */
  if (fslio==NULL)  FSLIOERR("FslSetStdXform: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
      fslio->niftiptr->sform_code = sform_code;
      fslio->niftiptr->sto_xyz.m[0][0] = stdmat.m[0][0];
      fslio->niftiptr->sto_xyz.m[0][1] = stdmat.m[0][1];
      fslio->niftiptr->sto_xyz.m[0][2] = stdmat.m[0][2];
      fslio->niftiptr->sto_xyz.m[0][3] = stdmat.m[0][3];
      fslio->niftiptr->sto_xyz.m[1][0] = stdmat.m[1][0];
      fslio->niftiptr->sto_xyz.m[1][1] = stdmat.m[1][1];
      fslio->niftiptr->sto_xyz.m[1][2] = stdmat.m[1][2];
      fslio->niftiptr->sto_xyz.m[1][3] = stdmat.m[1][3];
      fslio->niftiptr->sto_xyz.m[2][0] = stdmat.m[2][0];
      fslio->niftiptr->sto_xyz.m[2][1] = stdmat.m[2][1];
      fslio->niftiptr->sto_xyz.m[2][2] = stdmat.m[2][2];
      fslio->niftiptr->sto_xyz.m[2][3] = stdmat.m[2][3];
      fslio->niftiptr->sto_xyz.m[3][0] = 0;
      fslio->niftiptr->sto_xyz.m[3][1] = 0;
      fslio->niftiptr->sto_xyz.m[3][2] = 0;
      fslio->niftiptr->sto_xyz.m[3][3] = 1;
      fslio->niftiptr->sto_ijk = nifti_mat44_inverse(fslio->niftiptr->sto_xyz);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


short FslGetStdXform(FSLIO *fslio, mat44 *stdmat)
{
    /* returns sform code  (NB: stdmat must point to a 4x4 array) */
    float dx,dy,dz,tr;
    if (fslio==NULL)  FSLIOERR("FslGetStdXform: Null pointer passed for FSLIO");
    if (fslio->niftiptr!=NULL) {
        stdmat->m[0][0] = fslio->niftiptr->sto_xyz.m[0][0];
        stdmat->m[0][1] = fslio->niftiptr->sto_xyz.m[0][1];
        stdmat->m[0][2] = fslio->niftiptr->sto_xyz.m[0][2];
        stdmat->m[0][3] = fslio->niftiptr->sto_xyz.m[0][3];
        stdmat->m[1][0] = fslio->niftiptr->sto_xyz.m[1][0];
        stdmat->m[1][1] = fslio->niftiptr->sto_xyz.m[1][1];
        stdmat->m[1][2] = fslio->niftiptr->sto_xyz.m[1][2];
        stdmat->m[1][3] = fslio->niftiptr->sto_xyz.m[1][3];
        stdmat->m[2][0] = fslio->niftiptr->sto_xyz.m[2][0];
        stdmat->m[2][1] = fslio->niftiptr->sto_xyz.m[2][1];
        stdmat->m[2][2] = fslio->niftiptr->sto_xyz.m[2][2];
        stdmat->m[2][3] = fslio->niftiptr->sto_xyz.m[2][3];
        stdmat->m[3][0] = 0.0;
        stdmat->m[3][1] = 0.0;
        stdmat->m[3][2] = 0.0;
        stdmat->m[3][3] = 1.0;
        
        /* the code below gives a default but it really should never be used */
        if (fslio->niftiptr->sform_code == NIFTI_XFORM_UNKNOWN) {
            FslGetVoxDim(fslio,&dx,&dy,&dz,&tr);
            stdmat->m[0][0] = -dx;  /* default Radiological convention */
            stdmat->m[0][1] = 0;
            stdmat->m[0][2] = 0;
            stdmat->m[0][3] = 0;
            stdmat->m[1][0] = 0;
            stdmat->m[1][1] = dy;
            stdmat->m[1][2] = 0;
            stdmat->m[1][3] = 0;
            stdmat->m[2][0] = 0;
            stdmat->m[2][1] = 0;
            stdmat->m[2][2] = dz;
            stdmat->m[2][3] = 0;
            stdmat->m[3][0] = 0.0;
            stdmat->m[3][1] = 0.0;
            stdmat->m[3][2] = 0.0;
            stdmat->m[3][3] = 1.0;
        }
        return fslio->niftiptr->sform_code;
    }
    if (fslio->mincptr!=NULL) {
        fprintf(stderr,"Warning:: Minc is not yet supported\n");
    }
    return NIFTI_XFORM_UNKNOWN;
}


void FslSetRigidXform(FSLIO *fslio, short qform_code, mat44 rigidmat)
{
    /* NB: rigidmat must point to an allocated mat44 */
  float dx, dy, dz;
  if (fslio==NULL)  FSLIOERR("FslSetRigidXform: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
      fslio->niftiptr->qform_code = qform_code;
      fslio->niftiptr->qto_xyz.m[0][0] = rigidmat.m[0][0];
      fslio->niftiptr->qto_xyz.m[0][1] = rigidmat.m[0][1];
      fslio->niftiptr->qto_xyz.m[0][2] = rigidmat.m[0][2];
      fslio->niftiptr->qto_xyz.m[0][3] = rigidmat.m[0][3];
      fslio->niftiptr->qto_xyz.m[1][0] = rigidmat.m[1][0];
      fslio->niftiptr->qto_xyz.m[1][1] = rigidmat.m[1][1];
      fslio->niftiptr->qto_xyz.m[1][2] = rigidmat.m[1][2];
      fslio->niftiptr->qto_xyz.m[1][3] = rigidmat.m[1][3];
      fslio->niftiptr->qto_xyz.m[2][0] = rigidmat.m[2][0];
      fslio->niftiptr->qto_xyz.m[2][1] = rigidmat.m[2][1];
      fslio->niftiptr->qto_xyz.m[2][2] = rigidmat.m[2][2];
      fslio->niftiptr->qto_xyz.m[2][3] = rigidmat.m[2][3];
      fslio->niftiptr->qto_xyz.m[3][0] = 0;
      fslio->niftiptr->qto_xyz.m[3][1] = 0;
      fslio->niftiptr->qto_xyz.m[3][2] = 0;
      fslio->niftiptr->qto_xyz.m[3][3] = 1;
      nifti_mat44_to_quatern(
            fslio->niftiptr->qto_xyz,&(fslio->niftiptr->quatern_b),
            &(fslio->niftiptr->quatern_c),&(fslio->niftiptr->quatern_d),
            &(fslio->niftiptr->qoffset_x),&(fslio->niftiptr->qoffset_y),
            &(fslio->niftiptr->qoffset_z),&dx,&dy,&dz,&(fslio->niftiptr->qfac));
      fslio->niftiptr->qto_ijk = nifti_mat44_inverse(fslio->niftiptr->qto_xyz);

  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


short FslGetRigidXform(FSLIO *fslio, mat44 *rigidmat)
{
    /* returns qform code  (NB: rigidmat must point to an allocated mat44) */
    float dx,dy,dz,tr;
    if (fslio==NULL)  FSLIOERR("FslGetRigidXform: Null pointer passed for FSLIO");
    if (fslio->niftiptr!=NULL) {
        rigidmat->m[0][0] = fslio->niftiptr->qto_xyz.m[0][0];
        rigidmat->m[0][1] = fslio->niftiptr->qto_xyz.m[0][1];
        rigidmat->m[0][2] = fslio->niftiptr->qto_xyz.m[0][2];
        rigidmat->m[0][3] = fslio->niftiptr->qto_xyz.m[0][3];
        rigidmat->m[1][0] = fslio->niftiptr->qto_xyz.m[1][0];
        rigidmat->m[1][1] = fslio->niftiptr->qto_xyz.m[1][1];
        rigidmat->m[1][2] = fslio->niftiptr->qto_xyz.m[1][2];
        rigidmat->m[1][3] = fslio->niftiptr->qto_xyz.m[1][3];
        rigidmat->m[2][0] = fslio->niftiptr->qto_xyz.m[2][0];
        rigidmat->m[2][1] = fslio->niftiptr->qto_xyz.m[2][1];
        rigidmat->m[2][2] = fslio->niftiptr->qto_xyz.m[2][2];
        rigidmat->m[2][3] = fslio->niftiptr->qto_xyz.m[2][3];
        rigidmat->m[3][0] = 0.0;
        rigidmat->m[3][1] = 0.0;
        rigidmat->m[3][2] = 0.0;
        rigidmat->m[3][3] = 1.0;
        
        /* the code gives a default but it should never really be used */
        if (fslio->niftiptr->sform_code == NIFTI_XFORM_UNKNOWN) {
          FslGetVoxDim(fslio,&dx,&dy,&dz,&tr);
          rigidmat->m[0][0] = dx;
          rigidmat->m[0][1] = 0;
          rigidmat->m[0][2] = 0;
          rigidmat->m[0][3] = 0;
          rigidmat->m[1][0] = 0;
          rigidmat->m[1][1] = dy;
          rigidmat->m[1][2] = 0;
          rigidmat->m[1][3] = 0;
          rigidmat->m[2][0] = 0;
          rigidmat->m[2][1] = 0;
          rigidmat->m[2][2] = dz;
          rigidmat->m[3][0] = 0.0;
          rigidmat->m[3][1] = 0.0;
          rigidmat->m[3][2] = 0.0;
          rigidmat->m[3][3] = 1.0;
        }
        return fslio->niftiptr->qform_code;
    }
    if (fslio->mincptr!=NULL) {
        fprintf(stderr,"Warning:: Minc is not yet supported\n");
    }
    return NIFTI_XFORM_UNKNOWN;
}


void FslSetIntent(FSLIO *fslio, short intent_code, float p1, float p2, float p3)
{
  if (fslio==NULL)  FSLIOERR("FslSetIntent: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
      fslio->niftiptr->intent_code = intent_code;
      fslio->niftiptr->intent_p1 = p1;
      fslio->niftiptr->intent_p2 = p2;
      fslio->niftiptr->intent_p3 = p3;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


short FslGetIntent(FSLIO *fslio, short *intent_code, float *p1, float *p2,
                   float *p3)
{
  /* also returns intent code */
  if (fslio==NULL)  FSLIOERR("FslGetIntent: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
      *intent_code = fslio->niftiptr->intent_code;
      *p1 = fslio->niftiptr->intent_p1;
      *p2 = fslio->niftiptr->intent_p2;
      *p3 = fslio->niftiptr->intent_p3;
      return fslio->niftiptr->intent_code;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return NIFTI_INTENT_NONE;
}




void FslSetIntensityScaling(FSLIO *fslio, float slope, float intercept)
{
  if (fslio==NULL)  FSLIOERR("FslSetIntensityScaling: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
      fslio->niftiptr->scl_slope = slope;
      fslio->niftiptr->scl_inter = intercept;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


int FslGetIntensityScaling(FSLIO *fslio, float *slope, float *intercept)
{
  /* returns 1 if scaling required or 0 otherwise */
  if (fslio==NULL)  FSLIOERR("FslGetIntensityScaling: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *slope = fslio->niftiptr->scl_slope;
    *intercept = fslio->niftiptr->scl_inter;
    if (fabs(*slope)<1e-30) {
      *slope = 1.0;
      *intercept = 0.0;
      return 0;
    }
    if ( (fabs(*slope - 1.0)>1e-30) || (fabs(*intercept)>1e-30) ) {
      return 1;
    } else {
      return 0;
    }
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
 
}


mat33 mat44_to_mat33(mat44 x)
{
  mat33 y;
  int i,j;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      y.m[i][j] = x.m[i][j];
    }
  }
  return y;
}


int FslGetLeftRightOrder(FSLIO *fslio)
{
  /* Determines if the image is stored in neurological or radiological convention */
  int order=FSL_RADIOLOGICAL, sform_code, qform_code;
  float det=-1.0;
  mat44 sform44, qform44;
  mat33 sform33, qform33;
  if (fslio==NULL)  FSLIOERR("FslGetLeftRightOrder: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    sform_code = FslGetStdXform(fslio,&sform44);
    qform_code = FslGetRigidXform(fslio,&qform44);
    if (sform_code!=NIFTI_XFORM_UNKNOWN) { 
      sform33 = mat44_to_mat33(sform44);
      det = nifti_mat33_determ(sform33);
    } else if (qform_code!=NIFTI_XFORM_UNKNOWN) { 
      qform33 = mat44_to_mat33(qform44);
      det = nifti_mat33_determ(qform33); 
    }

    if (det<0.0) order=FSL_RADIOLOGICAL;
    else order=FSL_NEUROLOGICAL;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return order;
}



void FslSetAnalyzeSform(FSLIO *fslio, const short *orig,
                        float dx, float dy, float dz)
{
  /* Creates an sform matrix for an Analyze file */
  /* THIS ALWAYS CREATES A RADIOLOGICAL ORDERED SFORM */
  /* NB: the origin passed in here is in Analyze convention - starting at 1, not 0 */
  float x, y, z;
  if (fslio==NULL)  FSLIOERR("FslSetAnalyzeSform: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    if (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_ANALYZE) {
      /* default case */
      fslio->niftiptr->sform_code = NIFTI_XFORM_UNKNOWN;
    }
    /* ignore all zero origins - really all serious coord stuff should
       be done via the FslSetStdCoord call */
    if ((orig[0]!=0) || (orig[1]!=0) || (orig[2]!=0))
      {
        short origx=0, origy=0, origz=0;
        if ((orig[0]!=0) || (orig[1]!=0) || (orig[2]!=0)) {
          /* convert to nifti conventions (start at 0 not 1) */
          origx = orig[0] - 1;
          origy = orig[1] - 1;
          origz = orig[2] - 1;
        }
        if ( dx * dy * dz > 0 ) {
          /* change neurological convention to radiological if necessary */
          dx = -dx;
        }
        if ( (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_ANALYZE) 
             || (fslio->niftiptr->sform_code == NIFTI_XFORM_UNKNOWN) ) {
          /* make a default transform with the requested origin at xyz=000 */ 
          fslio->niftiptr->sform_code = NIFTI_XFORM_ALIGNED_ANAT;
          fslio->niftiptr->sto_xyz.m[0][0] = dx;
          fslio->niftiptr->sto_xyz.m[0][1] = 0;
          fslio->niftiptr->sto_xyz.m[0][2] = 0;
          fslio->niftiptr->sto_xyz.m[0][3] = -(origx)*(dx);
          fslio->niftiptr->sto_xyz.m[1][0] = 0;
          fslio->niftiptr->sto_xyz.m[1][1] = dy;
          fslio->niftiptr->sto_xyz.m[1][2] = 0;
          fslio->niftiptr->sto_xyz.m[1][3] = -(origy)*(dy);
          fslio->niftiptr->sto_xyz.m[2][0] = 0;
          fslio->niftiptr->sto_xyz.m[2][1] = 0;
          fslio->niftiptr->sto_xyz.m[2][2] = dz;
          fslio->niftiptr->sto_xyz.m[2][3] = -(origz)*(dz);
          fslio->niftiptr->sto_xyz.m[3][0] = 0;
          fslio->niftiptr->sto_xyz.m[3][1] = 0;
          fslio->niftiptr->sto_xyz.m[3][2] = 0;
          fslio->niftiptr->sto_xyz.m[3][3] = 1;
          fslio->niftiptr->sto_ijk =
                 nifti_mat44_inverse(fslio->niftiptr->sto_xyz);
        } else {
          /* update the existing origin */
          /* find out what the existing xyz of the requested origin is */
          x = fslio->niftiptr->sto_xyz.m[0][0] * origx
            + fslio->niftiptr->sto_xyz.m[0][1] * origy
            + fslio->niftiptr->sto_xyz.m[0][2] * origz
            + fslio->niftiptr->sto_xyz.m[0][3];
          y = fslio->niftiptr->sto_xyz.m[1][0] * origx
            + fslio->niftiptr->sto_xyz.m[1][1] * origy
            + fslio->niftiptr->sto_xyz.m[1][2] * origz
            + fslio->niftiptr->sto_xyz.m[1][3];
          z = fslio->niftiptr->sto_xyz.m[2][0] * origx
            + fslio->niftiptr->sto_xyz.m[2][1] * origy
            + fslio->niftiptr->sto_xyz.m[2][2] * origz
            + fslio->niftiptr->sto_xyz.m[2][3];
          /* subtract off whatever is currently the xyz of the origin */
          fslio->niftiptr->sto_xyz.m[0][3] -= x;
          fslio->niftiptr->sto_xyz.m[1][3] -= y;
          fslio->niftiptr->sto_xyz.m[2][3] -= z;
          fslio->niftiptr->sto_ijk =
                 nifti_mat44_inverse(fslio->niftiptr->sto_xyz);
        }
        
      }
    
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetAnalyzeOrigin(FSLIO *fslio, short orig[5])
{
  /* NB: orig returned here is in Analyze convention - starting at 1, not 0 */
  if (fslio==NULL)  FSLIOERR("FslGetAnalyzeOrigin: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    /* Use sform or qform to determine the origin - default is zero */
      orig[0]=0; 
      orig[1]=0; 
      orig[2]=0; 
      orig[3]=0; 
      orig[4]=0;

      if (fslio->niftiptr->qform_code != NIFTI_XFORM_UNKNOWN) {
        orig[0]=(short) fslio->niftiptr->qto_ijk.m[0][3] + 1;
        orig[1]=(short) fslio->niftiptr->qto_ijk.m[1][3] + 1;
        orig[2]=(short) fslio->niftiptr->qto_ijk.m[2][3] + 1;
      } 

      if (fslio->niftiptr->sform_code != NIFTI_XFORM_UNKNOWN) {
        orig[0]=(short) fslio->niftiptr->sto_ijk.m[0][3] + 1;
        orig[1]=(short) fslio->niftiptr->sto_ijk.m[1][3] + 1;
        orig[2]=(short) fslio->niftiptr->sto_ijk.m[2][3] + 1;
      } 
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}



/***************************************************************
 * FslClose
 ***************************************************************/
/*! \fn int FslClose(FSLIO *fslio)
    \brief  Write header and image data if this dataset was open for
        writing.  Close the dataset header and data files.

        
    \param fslio  pointer to FSLIO data structure
    \return  -1 on error, 0 OK ???.
 */
int FslClose(FSLIO *fslio)
{
  int retval=0, filetype;
  struct dsr *hdr;
  znzFile hptr=NULL;

  if (fslio==NULL)   return 0;

  /* close the (data) file */
  if (!znz_isnull(fslio->fileptr)) retval=znzclose(fslio->fileptr);

  /** ----- if writing the image, need to worry about the header bit ----- **/

  if ( (fslio->niftiptr!=NULL) && (FslGetWriteMode(fslio)==1) 
       && (fslio->written_hdr==0) ) {

    /* ensure that the type is set correctly */
    fslio->niftiptr->nifti_type = FslBaseFileType(FslGetFileType(fslio));

    /* must write the header now */
    filetype = FslGetFileType(fslio);
    strcpy(fslio->niftiptr->descrip,"FSL3.2beta");
    if (!FslIsSingleFileType(filetype)) {
      /* for file pairs - open new header file and write it */
      nifti_image_write_hdr_img(fslio->niftiptr,0,"wb");
    } else {
      /* for single files it is more complicated */
      if (!FslIsCompressedFileType(filetype)) {
        /* noncompressed -> reopen this file in r+ mode and write the header part again */
        nifti_image_write_hdr_img(fslio->niftiptr,0,"r+b");
      } else {
        /* compressed mode -> not possible! */
        fprintf(stderr,"Error:: header must be written before writing any other data.\n");
        return -1;
      }
    }
  }
    
  /* --- nasty hack to write the origin in Analyze files --- */

  if ( (FslGetWriteMode(fslio)==1) && (fslio->niftiptr!=NULL) && 
       (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_ANALYZE) ) {
 
    /* read in the old header, change the origin and write it out again */
    hdr = (struct dsr *) calloc(1,sizeof(struct dsr));
    FslReadRawHeader(hdr,fslio->niftiptr->fname);
    if (fslio->niftiptr->byteorder != nifti_short_order()) {AvwSwapHeader(hdr);}
    
    /* calculate origin from sform (if set) */
    {
      short blah[5];
      FslGetAnalyzeOrigin(fslio,blah);
      memcpy(hdr->hist.originator,blah,5*sizeof(short));
    
      /* Write out in radiological order if origin is non-zero */
      /* set negative pixdim if needed to keep LR orientation consistent */
      if ( (blah[0]!=0) || (blah[1]!=0) || (blah[2]!=0) ) {
        if (hdr->dime.pixdim[1] * hdr->dime.pixdim[2] * hdr->dime.pixdim[3] > 0) {
          hdr->dime.pixdim[1] = - hdr->dime.pixdim[1]; 
        }
      }
    }

    /* swap back byte order and write out */
    if (fslio->niftiptr->byteorder != nifti_short_order()) {AvwSwapHeader(hdr);}
    hptr = znzopen(fslio->niftiptr->fname,"wb",FslIsCompressedFileType(FslGetFileType(fslio)));
    if (znz_isnull(hptr)) {     
      fprintf(stderr,"Error:: Could not write origin data to header file %s.\n",
              fslio->niftiptr->fname);
      return -1;
    };
    
    znzwrite(hdr,1,sizeof(struct dsr),hptr);
    znzclose(hptr);
    free(hdr);
  }

  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
    return -1;
  }

  return retval;
}
  

void AvwSwapHeader(struct dsr *avw)
{
  char *ptr;

  ptr = (char *) &(avw->hk);
  nifti_swap_4bytes(1,ptr);             /* sizeof_hdr */
  ptr += 32;
  nifti_swap_4bytes(1,ptr);             /* extents */
  ptr += 4;
  nifti_swap_2bytes(1,ptr);             /* session_error */
  
  ptr = (char *) &(avw->dime);
  nifti_swap_2bytes(8,ptr);             /* dims */
  ptr += 28;
  nifti_swap_2bytes(4,ptr);             /* unused1, datatype, bitpix, dim_un0 */
  ptr += 8;
  nifti_swap_4bytes(18,ptr);            /* pixdim, vox_offset, ... */
                                        /* cal_min, compressed, ... glmin */

  ptr = (char *) &(avw->hist);
  ptr += 105;
  nifti_swap_2bytes(5,ptr);             /* originator (used to store origin) */
  ptr += 63;
  nifti_swap_4bytes(8,ptr);             /* views, ... smin */
}


int FslReadRawHeader(void *buffer, const char* filename)
{
  znzFile fp;
  int retval;
  fp = znzopen(filename,"rb",1);
  if (znz_isnull(fp)) {
    fprintf(stderr,"Could not open header %s\n",filename);
    return 0;
  }
  retval = znzread(buffer,1,348,fp);
  znzclose(fp);
  if (retval != 348) {
    fprintf(stderr,"Could not read header %s\n",filename);
    return retval;
  }
  return retval;
}

void FslSetOverrideOutputType(int type)
{
  if ( (type==-1) || (FslIsValidFileType(type)) ) {
    FslOverrideOutputType=type;
  } else {
    fprintf(stderr,"Invalid file type (%d) requested - ignoring this\n",type);
  }
}

int FslGetOverrideOutputType(void)
{
  return FslOverrideOutputType;
}

void FslSetIgnoreMFQ(int flag)
{
  assert((flag==0) || (flag==1));
  FslIgnoreMFQ=flag;
}


int FslGetIgnoreMFQ(void)
{
  return FslIgnoreMFQ;
}

/***************************************************************
 * FslReadHeader
 ***************************************************************/
/*! \fn FSLIO * FslReadHeader(char *fname)
    \brief Reads nifti/anz header, no data is read
        
    \param fname        filename specification (could be .img,.hdr,.nii, or no ext
    \return FSLIO data structure with the nifti_image structure fields filled 
            as per fname header.
            NULL on error 
 */
FSLIO * FslReadHeader(char *fname)
{
   char *hdrname, *imgname;
   FSLIO *fslio;


   fslio = FslInit();
  
  /** get header file name */
  FslGetHdrImgNames(fname, fslio, &hdrname, &imgname);

  /** read header information */
  fslio->niftiptr = nifti_image_read(hdrname, 0);

  if (fslio->niftiptr == NULL) {
        FSLIOERR("FslReadHeader: error reading header information");
        return(NULL);
  }

  fslio->file_mode = FslGetReadFileType(fslio);

  return(fslio);
}


/***************************************************************
 * FslGetVolumeAsScaledDouble
 ***************************************************************/
/*! \fn double *** FslGetVolumeAsScaledDouble(FSLIO *fslio, int vol)
    \brief Return volume #vol (0-based) as a 3D array of scaled doubles. 

        Volume Array is indexed as [0..zdim-1][0..ydim-1][0..xdim-1].  
        <br>The array will be byteswapped to native-endian.
        <br>Array values are scaled as per fslio header slope and intercept fields.

    \param fslio pointer to open dataset
    \param vol volume number to read (legal range [0..tdim-1])
    \return Pointer to 3D double array, NULL on error
 */
double ***FslGetVolumeAsScaledDouble(FSLIO *fslio, int vol)
{
  double ***newbuf;
  void *diskbuf;
  int xx,yy,zz;
  int ret;
  float inter, slope;
  int dims_to_get[8];
  int i;

  if (fslio==NULL)  FSLIOERR("FslGetVolumeAsScaledDouble: Null pointer passed for FSLIO");

  if ((fslio->niftiptr->dim[0] < 3) || (fslio->niftiptr->dim[0] > 4))
        FSLIOERR("FslGetVolumeAsScaledDouble: Incorrect dataset dimension, 3D-4D needed");

  /***** nifti dataset */
  if (fslio->niftiptr!=NULL) {
        xx = (fslio->niftiptr->nx == 0 ? 1 : (long)fslio->niftiptr->nx);
        yy = (fslio->niftiptr->ny == 0 ? 1 : (long)fslio->niftiptr->ny);
        zz = (fslio->niftiptr->nz == 0 ? 1 : (long)fslio->niftiptr->nz);

        if (fslio->niftiptr->scl_slope == 0) {
                slope = 1.0;
                inter = 0.0;
        }
        else {
                slope = fslio->niftiptr->scl_slope;
                inter = fslio->niftiptr->scl_inter;
        }
        

    /** allocate new 3D buffer */
    newbuf = d3matrix(zz-1,yy-1,xx-1);


    /** read in the data in disk format */
    dims_to_get[0] = 0;
    for (i=1; i<8; i++)
        dims_to_get[i] = -1;
    dims_to_get[4] = vol;


    diskbuf = NULL;
    ret = nifti_read_collapsed_image(fslio->niftiptr, dims_to_get, &diskbuf );
    if (ret <= 0) {
        fprintf(stderr,"ERROR:: read of disk buffer for volume %d from %s failed.\n",vol,fslio->niftiptr->iname);
        return(NULL);
    }


    /** cvt disk buffer to scaled double buffer  */
    ret = convertBufferToScaledDouble(newbuf[0][0], diskbuf, (long)(xx*yy*zz), slope, inter, fslio->niftiptr->datatype);

    free(diskbuf);

    if (ret == 0)
        return(newbuf);
    else
        return(NULL);

  } /* nifti data */


  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }

  return(NULL);
}




/***************************************************************
 * FslGetBufferAsScaledDouble
 ***************************************************************/
/*! \fn double **** FslGetBufferAsScaledDouble(FSLIO *fslio)
    \brief Return the fslio data buffer of a 1-4D dataset as a 4D array of 
        scaled doubles. 

        Array is indexed as buf[0..tdim-1][0..zdim-1][0..ydim-1][0..xdim-1].  
        <br>The array will be byteswapped to native-endian.
        <br>Array values are scaled as per fslio header slope and intercept fields.

    \param fslio pointer to open dataset
    \return Pointer to 4D double array, NULL on error
 */
double ****FslGetBufferAsScaledDouble(FSLIO *fslio)
{
  double ****newbuf;
  int xx,yy,zz,tt;
  int ret;
  float inter, slope;

  if (fslio==NULL)  FSLIOERR("FslGetBufferAsScaledDouble: Null pointer passed for FSLIO");

  if ((fslio->niftiptr->dim[0] <= 0) || (fslio->niftiptr->dim[0] > 4))
        FSLIOERR("FslGetBufferAsScaledDouble: Incorrect dataset dimension, 1-4D needed");

  /***** nifti dataset */
  if (fslio->niftiptr!=NULL) {
        xx = (fslio->niftiptr->nx == 0 ? 1 : (long)fslio->niftiptr->nx);
        yy = (fslio->niftiptr->ny == 0 ? 1 : (long)fslio->niftiptr->ny);
        zz = (fslio->niftiptr->nz == 0 ? 1 : (long)fslio->niftiptr->nz);
        tt = (fslio->niftiptr->nt == 0 ? 1 : (long)fslio->niftiptr->nt);

        if (fslio->niftiptr->scl_slope == 0) {
                slope = 1.0;
                inter = 0.0;
        }
        else {
                slope = fslio->niftiptr->scl_slope;
                inter = fslio->niftiptr->scl_inter;
        }
        

    /** allocate new 4D buffer */
    newbuf = d4matrix(tt-1,zz-1,yy-1,xx-1);

    /** cvt it */
    ret = convertBufferToScaledDouble(newbuf[0][0][0], fslio->niftiptr->data, (long)(xx*yy*zz*tt), slope, inter, fslio->niftiptr->datatype);

    if (ret == 0)
        return(newbuf);
    else
        return(NULL);

  } /* nifti data */


  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }

  return(NULL);
}

/***************************************************************
 * convertBufferToScaledDouble
 ***************************************************************/
/*! \fn int  convertBufferToScaledDouble(double *outbuf, void *inbuf, long len, float slope, float inter, int nifti_datatype )
    \brief allocate a 4D buffer, use 1 contiguous buffer for the data 

        Array is indexed as buf[0..th-1][0..zh-1][0..yh-1][0..xh-1].  
        <br>To access all elements as a vector, use buf[0][0][0][i] where
        i can range from 0 to th*zh*yh*xh - 1.

    \param outbuf pointer to array of doubles of size len
    \param inbuf void pointer to an array of len items of datatype nifti_datatype
    \param len number of elements in outbuf and inbuf
    \param slope slope term of scaling to be applied
    \param inter intercept term of scaling to be applied:  out = (in*slope)+inter
    \param nifti_datatype NIFTI datatype code for the datatype of the elements in inbuf
    \return error code: 0=OK -1=error
 */
int  convertBufferToScaledDouble(double *outbuf, void *inbuf, long len, float slope, float inter, int nifti_datatype ) 
{

        long i;


    /** fill the buffer */
    for (i=0; i<len; i++)
        switch(nifti_datatype) {
            case NIFTI_TYPE_UINT8:
                outbuf[i] = (double) ( *((THIS_UINT8 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_INT8:
                outbuf[i] = (double) ( *((THIS_INT8 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_UINT16:
                outbuf[i] = (double) ( *((THIS_UINT16 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_INT16:
                outbuf[i] = (double) ( *((THIS_INT16 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_UINT64:
                outbuf[i] = (double) ( *((THIS_UINT64 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_INT64:
                outbuf[i] = (double) ( *((THIS_INT64 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_UINT32:
                outbuf[i] = (double) ( *((THIS_UINT32 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_INT32:
                outbuf[i] = (double) ( *((THIS_INT32 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_FLOAT32:
                outbuf[i] = (double) ( *((THIS_FLOAT32 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_FLOAT64:
                outbuf[i] = (double) ( *((THIS_FLOAT64 *)(inbuf)+i) * slope + inter);
                break;

            case NIFTI_TYPE_FLOAT128:
            case NIFTI_TYPE_COMPLEX128:
            case NIFTI_TYPE_COMPLEX256:
            case NIFTI_TYPE_COMPLEX64:
            default:
                fprintf(stderr, "\nWarning, cannot support %s yet.\n",nifti_datatype_string(nifti_datatype));
                return(-1);
        }

return(0);
}

/***************************************************************
 * d3matrix
 ***************************************************************/
/*! \fn double ****d3matrix(int zh,  int yh, int xh)
    \brief allocate a 3D buffer, use 1 contiguous buffer for the data 

        Array is indexed as buf[0..zh][0..yh][0..xh].  
        <br>To access all elements as a vector, use buf[0][0][i] where
        i can range from 0 to zh*yh*xh - 1.
        Adaptation of Numerical Recipes in C nrutil.c allocation routines. 

    \param zh slowest changing dimension
    \param yh 2nd fastest changing dimension
    \param xh fastest changing dimension
    \return Pointer to 3D double array
 */
double ***d3matrix(int zh,  int yh, int xh)
{

        int j;
        int nslice = zh+1;
        int nrow = yh+1;
        int ncol = xh+1;
        double ***t;


        /** allocate pointers to slices */
        t=(double ***) malloc((size_t)((nslice)*sizeof(double**)));
        if (!t) FSLIOERR("d3matrix: allocation failure");

        /** allocate pointers for ydim */
        t[0]=(double **) malloc((size_t)((nslice*nrow)*sizeof(double*)));
        if (!t[0]) FSLIOERR("d3matrix: allocation failure");


        /** allocate the data blob */
        t[0][0]=(double *) malloc((size_t)((nslice*nrow*ncol)*sizeof(double)));
        if (!t[0][0]) FSLIOERR("d3matrix: allocation failure");


        /** point everything to the data blob */
        for(j=1;j<nrow*nslice;j++) t[0][j]=t[0][j-1]+ncol;
        for(j=1;j<nslice;j++) t[j]=t[j-1]+nrow;

        return t;
}


/***************************************************************
 * d4matrix
 ***************************************************************/
/*! \fn double ****d4matrix(int th, int zh,  int yh, int xh)
    \brief allocate a 4D buffer, use 1 contiguous buffer for the data 

        Array is indexed as buf[0..th][0..zh][0..yh][0..xh].  
        <br>To access all elements as a vector, use buf[0][0][0][i] where
        i can range from 0 to th*zh*yh*xh - 1.
        Adaptation of Numerical Recipes in C nrutil.c allocation routines. 

    \param th slowest changing dimension
    \param zh 2nd slowest changing dimension
    \param yh 2nd fastest changing dimension
    \param xh fastest changing dimension
    \return Pointer to 4D double array
 */
double ****d4matrix(int th, int zh,  int yh, int xh)
{

        int j;
        int nvol = th+1;
        int nslice = zh+1;
        int nrow = yh+1;
        int ncol = xh+1;
        double ****t;


        /** allocate pointers to vols */
        t=(double ****) malloc((size_t)((nvol)*sizeof(double***)));
        if (!t) FSLIOERR("d4matrix: allocation failure");

        /** allocate pointers to slices */
        t[0]=(double ***) malloc((size_t)((nvol*nslice)*sizeof(double**)));
        if (!t[0]) FSLIOERR("d4matrix: allocation failure");

        /** allocate pointers for ydim */
        t[0][0]=(double **) malloc((size_t)((nvol*nslice*nrow)*sizeof(double*)));
        if (!t[0][0]) FSLIOERR("d4matrix: allocation failure");


        /** allocate the data blob */
        t[0][0][0]=(double *) malloc((size_t)((nvol*nslice*nrow*ncol)*sizeof(double)));
        if (!t[0][0][0]) FSLIOERR("d4matrix: allocation failure");


        /** point everything to the data blob */
        for(j=1;j<nrow*nslice*nvol;j++) t[0][0][j]=t[0][0][j-1]+ncol;
        for(j=1;j<nslice*nvol;j++) t[0][j]=t[0][j-1]+nrow;
        for(j=1;j<nvol;j++) t[j]=t[j-1]+nslice;

        return t;
}

