/** \file znzlib.c
    \brief Low level i/o interface to compressed and noncompressed files.
        Written by Mark Jenkinson, FMRIB

This library provides an interface to both compressed (gzip/zlib) and
uncompressed (normal) file IO.  The functions are written to have the
same interface as the standard file IO functions.

To use this library instead of normal file IO, the following changes
are required:
 - replace all instances of FILE* with znzFile
 - change the name of all function calls, replacing the initial character
   f with the znz  (e.g. fseek becomes znzseek)
   one exception is rewind() -> znzrewind()
 - add a third parameter to all calls to znzopen (previously fopen)
   that specifies whether to use compression (1) or not (0)
 - use znz_isnull rather than any (pointer == NULL) comparisons in the code
   for znzfile types (normally done after a return from znzopen)
 
NB: seeks for writable files with compression are quite restricted

 */

#include "znzlib.h"

/*
znzlib.c  (zipped or non-zipped library)

*****            This code is released to the public domain.            *****

*****  Author: Mark Jenkinson, FMRIB Centre, University of Oxford       *****
*****  Date:   September 2004                                           *****

*****  Neither the FMRIB Centre, the University of Oxford, nor any of   *****
*****  its employees imply any warranty of usefulness of this software  *****
*****  for any purpose, and do not assume any liability for damages,    *****
*****  incidental or otherwise, caused by any use of this document.     *****

*/


/* Note extra argument (use_compression) where 
   use_compression==0 is no compression
   use_compression!=0 uses zlib (gzip) compression
*/

znzFile znzopen(const char *path, const char *mode, int use_compression)
{
  znzFile file;
  file = (znzFile) calloc(1,sizeof(struct znzptr));
  if( file == NULL ){
     fprintf(stderr,"** ERROR: znzopen failed to alloc znzptr\n");
     return NULL;
  }

  file->nzfptr = NULL;

#ifdef HAVE_ZLIB
  file->zfptr = NULL;

  if (use_compression) {
    file->withz = 1;
    if((file->zfptr = gzopen(path,mode)) == NULL) {
        free(file);
        file = NULL;
    }
  } else {
#endif

    file->withz = 0;
    if((file->nzfptr = fopen(path,mode)) == NULL) {
      free(file);
      file = NULL;
    }

#ifdef HAVE_ZLIB
  }
#endif

  return file;
}


znzFile znzdopen(int fd, const char *mode, int use_compression)
{
  znzFile file;
  file = (znzFile) calloc(1,sizeof(struct znzptr));
  if( file == NULL ){
     fprintf(stderr,"** ERROR: znzdopen failed to alloc znzptr\n");
     return NULL;
  }
#ifdef HAVE_ZLIB
  if (use_compression) {
    file->withz = 1;
    file->zfptr = gzdopen(fd,mode);
    file->nzfptr = NULL;
  } else {
#endif
    file->withz = 0;
#ifdef HAVE_FDOPEN
    file->nzfptr = fdopen(fd,mode);
#endif
#ifdef HAVE_ZLIB
    file->zfptr = NULL;
  };
#endif
  return file;
}


int Xznzclose(znzFile * file)
{
  int retval = 0;
  if (*file!=NULL) {
#ifdef HAVE_ZLIB
    if ((*file)->zfptr!=NULL)  { retval = gzclose((*file)->zfptr); }
#endif
    if ((*file)->nzfptr!=NULL) { retval = fclose((*file)->nzfptr); }
                                                                                
    free(*file);
    *file = NULL;
  }
  return retval;
}


size_t znzread(void* buf, size_t size, size_t nmemb, znzFile file)
{
  if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
  if (file->zfptr!=NULL) 
    return (size_t) (gzread(file->zfptr,buf,((int) size)*((int) nmemb)) / size);
#endif
  return fread(buf,size,nmemb,file->nzfptr);
}

size_t znzwrite(const void* buf, size_t size, size_t nmemb, znzFile file)
{
  if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
  if (file->zfptr!=NULL)
      {
      /*  NOTE:  We must typecast const away from the buffer because
          gzwrite does not have complete const specification */
    return (size_t) ( gzwrite(file->zfptr,(void *)buf,size*nmemb) / size );
      }
#endif
  return fwrite(buf,size,nmemb,file->nzfptr);
}

long znzseek(znzFile file, long offset, int whence)
{
  if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
  if (file->zfptr!=NULL) return (long) gzseek(file->zfptr,offset,whence);
#endif
  return fseek(file->nzfptr,offset,whence);
}

int znzrewind(znzFile stream)
{
  if (stream==NULL) { return 0; }
#ifdef HAVE_ZLIB
  /* On some systems, gzrewind() fails for uncompressed files.
     Use gzseek(), instead.               10, May 2005 [rickr]

     if (stream->zfptr!=NULL) return gzrewind(stream->zfptr);
  */

  if (stream->zfptr!=NULL) return (int)gzseek(stream->zfptr, 0L, SEEK_SET);
#endif
  rewind(stream->nzfptr);
  return 0;
}

long znztell(znzFile file)
{
  if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
  if (file->zfptr!=NULL) return (long) gztell(file->zfptr);
#endif
  return ftell(file->nzfptr);
}

int znzputs(const char * str, znzFile file)
{
  if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
  if (file->zfptr!=NULL) return gzputs(file->zfptr,str);
#endif
  return fputs(str,file->nzfptr);
}


char * znzgets(char* str, int size, znzFile file)
{
  if (file==NULL) { return NULL; }
#ifdef HAVE_ZLIB
  if (file->zfptr!=NULL) return gzgets(file->zfptr,str,size);
#endif
  return fgets(str,size,file->nzfptr);
}


int znzflush(znzFile file)
{
  if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
  if (file->zfptr!=NULL) return gzflush(file->zfptr,Z_SYNC_FLUSH);
#endif
  return fflush(file->nzfptr);
}


int znzeof(znzFile file)
{
  if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
  if (file->zfptr!=NULL) return gzeof(file->zfptr);
#endif
  return feof(file->nzfptr);
}


int znzputc(int c, znzFile file)
{
  if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
  if (file->zfptr!=NULL) return gzputc(file->zfptr,c);
#endif
  return fputc(c,file->nzfptr);
}


int znzgetc(znzFile file)
{
  if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
  if (file->zfptr!=NULL) return gzgetc(file->zfptr);
#endif
  return fgetc(file->nzfptr);
}

#if !defined (WIN32)
int znzprintf(znzFile stream, const char *format, ...)
{
  int retval=0;
  char *tmpstr;
  va_list va;
  if (stream==NULL) { return 0; }
  va_start(va, format);
#ifdef HAVE_ZLIB
  if (stream->zfptr!=NULL) {
    int size;  /* local to HAVE_ZLIB block */
    size = strlen(format) + 1000000;  /* overkill I hope */
    tmpstr = (char *)calloc(1, size);
    if( tmpstr == NULL ){
       fprintf(stderr,"** ERROR: znzprintf failed to alloc %d bytes\n", size);
       return retval;
    }
    vsprintf(tmpstr,format,va);
    retval=gzprintf(stream->zfptr,"%s",tmpstr);
    free(tmpstr);
  } else 
#endif
  {
   retval=vfprintf(stream->nzfptr,format,va);
  }
  va_end(va);
  return retval;
}

#endif

