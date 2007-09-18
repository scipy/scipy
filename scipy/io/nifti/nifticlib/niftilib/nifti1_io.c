#define _NIFTI1_IO_C_

#include "nifti1_io.h"   /* typedefs, prototypes, macros, etc. */

/*****===================================================================*****/
/*****     Sample functions to deal with NIFTI-1 and ANALYZE files       *****/
/*****...................................................................*****/
/*****            This code is released to the public domain.            *****/
/*****...................................................................*****/
/*****  Author: Robert W Cox, SSCC/DIRP/NIMH/NIH/DHHS/USA/EARTH          *****/
/*****  Date:   August 2003                                              *****/
/*****...................................................................*****/
/*****  Neither the National Institutes of Health (NIH), nor any of its  *****/
/*****  employees imply any warranty of usefulness of this software for  *****/
/*****  any purpose, and do not assume any liability for damages,        *****/
/*****  incidental or otherwise, caused by any use of this document.     *****/
/*****===================================================================*****/

/** \file nifti1_io.c
    \brief main collection of nifti1 i/o routines
           - written by Bob Cox, SSCC NIMH
           - revised by Mark Jenkinson, FMRIB
           - revised by Rick Reynolds, SSCC, NIMH
           - revised by Kate Fissell, University of Pittsburgh

        The library history can be viewed via "nifti_tool -nifti_hist".
    <br>The library version can be viewed via "nifti_tool -nifti_ver".
 */

/*! global history and version strings, for printing */
static char * gni_history[] = 
{
  "----------------------------------------------------------------------\n"
  "history (of nifti library changes):\n"
  "\n",
  "0.0  August, 2003 [rwcox]\n"
  "     (Robert W Cox of the National Institutes of Health, SSCC/DIRP/NIMH)\n"
  "   - initial version\n"
  "\n",
  "0.1  July/August, 2004 [Mark Jenkinson]\n"
  "     (FMRIB Centre, University of Oxford, UK)\n"
  "   - Mainly adding low-level IO and changing things to allow gzipped\n"
  "     files to be read and written\n"
  "   - Full backwards compatability should have been maintained\n"
  "\n",
  "0.2  16 Nov 2004 [rickr]\n"
  "     (Rick Reynolds of the National Institutes of Health, SSCC/DIRP/NIMH)\n"
  "   - included Mark's changes in the AFNI distribution (including znzlib/)\n"
  "     (HAVE_ZLIB is commented out for the standard distribution)\n"
  "   - modified nifti_validfilename() and nifti_makebasename()\n"
  "   - added nifti_find_file_extension()\n"
  "\n",
  "0.3  3 Dec 2004 [rickr]\n"
  "   - note: header extensions are not yet checked for\n"
  "   - added formatted history as global string, for printing\n"
  "   - added nifti_disp_lib_hist(), to display the nifti library history\n"
  "   - added nifti_disp_lib_version(), to display the nifti library history\n",
  "   - re-wrote nifti_findhdrname()\n"
  "       o used nifti_find_file_extension()\n"
  "       o changed order of file tests (default is .nii, depends on input)\n"
  "       o free hdrname on failure\n"
  "   - made similar changes to nifti_findimgname()\n"
  "   - check for NULL return from nifti_findhdrname() calls\n",
  "   - removed most of ERREX() macros\n"
  "   - modified nifti_image_read()\n"
  "       o added debug info and error checking (on gni_debug > 0, only)\n"
  "       o fail if workingname is NULL\n"
  "       o check for failure to open header file\n"
  "       o free workingname on failure\n"
  "       o check for failure of nifti_image_load()\n"
  "       o check for failure of nifti_convert_nhdr2nim()\n",
  "   - changed nifti_image_load() to int, and check nifti_read_buffer return\n"
  "   - changed nifti_read_buffer() to fail on short read, and to count float\n"
  "     fixes (to print on debug)\n"
  "   - changed nifti_image_infodump to print to stderr\n"
  "   - updated function header comments, or moved comments above header\n"
  "   - removed const keyword\n"
  "   - added LNI_FERR() macro for error reporting on input files\n"
  "\n",
  "0.4  10 Dec 2004 [rickr]  - added header extensions\n"
  "   - in nifti1_io.h:\n"
  "       o added num_ext and ext_list to the definition of nifti_image\n"
  "       o made many functions static (more to follow)\n"
  "       o added LNI_MAX_NIA_EXT_LEN, for max nifti_type 3 extension length\n",
  "   - added __DATE__ to version output in nifti_disp_lib_version()\n"
  "   - added nifti_disp_matrix_orient() to print orientation information\n"
  "   - added '.nia' as a valid file extension in nifti_find_file_extension()\n"
  "   - added much more debug output\n"
  "   - in nifti_image_read(), in the case of an ASCII header, check for\n"
  "     extensions after the end of the header\n",
  "   - added nifti_read_extensions() function\n"
  "   - added nifti_read_next_extension() function\n"
  "   - added nifti_add_exten_to_list() function\n"
  "   - added nifti_check_extension() function\n"
  "   - added nifti_write_extensions() function\n"
  "   - added nifti_extension_size() function\n"
  "   - in nifti_set_iname_offest():\n"
  "       o adjust offset by the extension size and the extender size\n",
  "       o fixed the 'ceiling modulo 16' computation\n"
  "   - in nifti_image_write_hdr_img2(): \n"
  "       o added extension writing\n"
  "       o check for NULL return from nifti_findimgname()\n"
  "   - include number of extensions in nifti_image_to_ascii() output\n"
  "   - in nifti_image_from_ascii():\n"
  "       o return bytes_read as a parameter, computed from the final spos\n"
  "       o extract num_ext from ASCII header\n"
  "\n",
  "0.5  14 Dec 2004 [rickr]  - added sub-brick reading functions\n"
  "   - added nifti_brick_list type to nifti1_io.h, along with new prototypes\n"
  "   - added main nifti_image_read_bricks() function, with description\n"
  "   - added nifti_image_load_bricks() - library function (requires nim)\n"
  "   - added valid_nifti_brick_list() - library function\n"
  "   - added free_NBL() - library function\n",
  "   - added update_nifti_image_for_brick_list() for dimension update\n"
  "   - added nifti_load_NBL_bricks(), nifti_alloc_NBL_mem(),\n"
  "           nifti_copynsort() and force_positive() (static functions)\n"
  "   - in nifti_image_read(), check for failed load only if read_data is set\n"
  "   - broke most of nifti_image_load() into nifti_image_load_prep()\n"
  "\n",
  "0.6  15 Dec 2004 [rickr]  - added sub-brick writing functionality\n"
  "   - in nifti1_io.h, removed znzlib directory from include - all nifti\n"
  "       library files are now under the nifti directory\n"
  "   - nifti_read_extensions(): print no offset warning for nifti_type 3\n"
  "   - nifti_write_all_data():\n"
  "       o pass nifti_brick_list * NBL, for optional writing\n"
  "       o if NBL, write each sub-brick, sequentially\n",
  "   - nifti_set_iname_offset(): case 1 must have sizeof() cast to int\n"
  "   - pass NBL to nifti_image_write_hdr_img2(), and allow NBL or data\n"
  "   - added nifti_image_write_bricks() wrapper for ...write_hdr_img2()\n"
  "   - included compression abilities\n"
  "\n",
  "0.7  16 Dec 2004 [rickr] - minor changes to extension reading\n"
  "\n",
  "0.8  21 Dec 2004 [rickr] - restrict extension reading, and minor changes\n"
  "   - in nifti_image_read(), compute bytes for extensions (see remaining)\n"
  "   - in nifti_read_extensions(), pass 'remain' as space for extensions,\n"
  "        pass it to nifti_read_next_ext(), and update for each one read \n"
  "   - in nifti_check_extension(), require (size <= remain)\n",
  "   - in update_nifti_image_brick_list(), update nvox\n"
  "   - in nifti_image_load_bricks(), make explicit check for nbricks <= 0\n"
  "   - in int_force_positive(), check for (!list)\n"
  "   - in swap_nifti_header(), swap sizeof_hdr, and reorder to struct order\n"
  "   - change get_filesize functions to signed ( < 0 is no file or error )\n",
  "   - in nifti_validfilename(), lose redundant (len < 0) check\n"
  "   - make print_hex_vals() static\n"
  "   - in disp_nifti_1_header, restrict string field widths\n"
  "\n",
  "0.9  23 Dec 2004 [rickr] - minor changes\n"
  "   - broke ASCII header reading out of nifti_image_read(), into new\n"
  "        functions has_ascii_header() and read_ascii_image()\n",
  "   - check image_read failure and znzseek failure\n"
  "   - altered some debug output\n"
  "   - nifti_write_all_data() now returns an int\n"
  "\n",
  "0.10 29 Dec 2004 [rickr]\n"
  "   - renamed nifti_valid_extension() to nifti_check_extension()\n"
  "   - added functions nifti_makehdrname() and nifti_makeimgname()\n"
  "   - added function valid_nifti_extensions()\n"
  "   - in nifti_write_extensions(), check for validity before writing\n",
  "   - rewrote nifti_image_write_hdr_img2():\n"
  "       o set write_data and leave_open flags from write_opts\n"
  "       o add debug print statements\n"
  "       o use nifti_write_ascii_image() for the ascii case\n"
  "       o rewrote the logic of all cases to be easier to follow\n",
  "   - broke out code as nifti_write_ascii_image() function\n"
  "   - added debug to top-level write functions, and free the znzFile\n"
  "   - removed unused internal function nifti_image_open()\n"
  "\n",
  "0.11 30 Dec 2004 [rickr] - small mods\n"
  "   - moved static function prototypes from header to C file\n"
  "   - free extensions in nifti_image_free()\n"
  "\n",
  "1.0  07 Jan 2005 [rickr] - INITIAL RELEASE VERSION\n"
  "   - added function nifti_set_filenames()\n"
  "   - added function nifti_read_header()\n"
  "   - added static function nhdr_looks_good()\n"
  "   - added static function need_nhdr_swap()\n"
  "   - exported nifti_add_exten_to_list symbol\n",
  "   - fixed #bytes written in nifti_write_extensions()\n"
  "   - only modify offset if it is too small (nifti_set_iname_offset)\n"
  "   - added nifti_type 3 to nifti_makehdrname and nifti_makeimgname\n"
  "   - added function nifti_set_filenames()\n"
  "\n",
  "1.1  07 Jan 2005 [rickr]\n"
  "   - in nifti_read_header(), swap if needed\n"
  "\n",
  "1.2  07 Feb 2005 [kate fissell c/o rickr] \n"
  "   - nifti1.h: added doxygen comments for main struct and #define groups\n"
  "   - nifti1_io.h: added doxygen comments for file and nifti_image struct\n"
  "   - nifti1_io.h: added doxygen comments for file and some functions\n"
  "   - nifti1_io.c: changed nifti_copy_nim_info to use memcpy\n"
  "\n",
  "1.3  09 Feb 2005 [rickr]\n"
  "   - nifti1.h: added doxygen comments for extension structs\n"
  "   - nifti1_io.h: put most #defines in #ifdef _NIFTI1_IO_C_ block\n"
  "   - added a doxygen-style description to every exported function\n"
  "   - added doxygen-style comments within some functions\n"
  "   - re-exported many znzFile functions that I had made static\n"
  "   - re-added nifti_image_open (sorry, Mark)\n"
  "   - every exported function now has 'nifti' in the name (19 functions)\n",
  "   - made sure every alloc() has a failure test\n"
  "   - added nifti_copy_extensions function, for use in nifti_copy_nim_info\n"
  "   - nifti_is_gzfile: added initial strlen test\n"
  "   - nifti_set_filenames: added set_byte_order parameter option\n"
  "     (it seems appropriate to set the BO when new files are associated)\n"
  "   - disp_nifti_1_header: prints to stdout (a.o.t. stderr), with fflush\n"
  "\n",
  "1.4  23 Feb 2005 [rickr] - sourceforge merge\n"
  "   - merged into the nifti_io CVS directory structure at sourceforge.net\n"
  "   - merged in 4 changes by Mark, and re-added his const keywords\n"
  "   - cast some pointers to (void *) for -pedantic compile option\n"
  "   - added nifti_free_extensions()\n"
  "\n",
  "1.5  02 Mar 2005 [rickr] - started nifti global options\n"
  "   - gni_debug is now g_opts.debug\n"
  "   - added validity check parameter to nifti_read_header\n"
  "   - need_nhdr_swap no longer does test swaps on the stack\n"
  "\n",
  "1.6  05 April 2005 [rickr] - validation and collapsed_image_read\n"
  "   - added nifti_read_collapsed_image(), an interface for reading partial\n"
  "     datasets, specifying a subset of array indices\n"
  "   - for read_collapsed_image, added static functions: rci_read_data(),\n"
  "     rci_alloc_mem(), and make_pivot_list()\n",
  "   - added nifti_nim_is_valid() to check for consistency (more to do)\n"
  "   - added nifti_nim_has_valid_dims() to do many dimensions tests\n"
  "\n",
  "1.7  08 April 2005 [rickr]\n"
  "   - added nifti_update_dims_from_array() - to update dimensions\n"
  "   - modified nifti_makehdrname() and nifti_makeimgname():\n"
  "       if prefix has a valid extension, use it (else make one up)\n"
  "   - added nifti_get_intlist - for making an array of ints\n"
  "   - fixed init of NBL->bsize in nifti_alloc_NBL_mem()  {thanks, Bob}\n"
  "\n",
  "1.8  14 April 2005 [rickr]\n"
  "   - added nifti_set_type_from_names(), for nifti_set_filenames()\n"
  "     (only updates type if number of files does not match it)\n"
  "   - added is_valid_nifti_type(), just to be sure\n"
  "   - updated description of nifti_read_collapsed_image() for *data change\n"
  "     (if *data is already set, assume memory exists for results)\n"
  "   - modified rci_alloc_mem() to allocate only if *data is NULL\n"
  "\n",
  "1.9  19 April 2005 [rickr]\n"
  "   - added extension codes NIFTI_ECODE_COMMENT and NIFTI_ECODE_XCEDE\n"
  "   - added nifti_type codes NIFTI_MAX_ECODE and NIFTI_MAX_FTYPE\n"
  "   - added nifti_add_extension() {exported}\n"
  "   - added nifti_fill_extension() as a static function\n"
  "   - added nifti_is_valid_ecode() {exported}\n",
  "   - nifti_type values are now NIFTI_FTYPE_* file codes\n"
  "   - in nifti_read_extensions(), decrement 'remain' by extender size, 4\n"
  "   - in nifti_set_iname_offset(), case 1, update if offset differs\n"
  "   - only output '-d writing nifti file' if debug > 1\n"
  "\n",
  "1.10 10 May 2005 [rickr]\n"
  "   - files are read using ZLIB only if they end in '.gz'\n"
  "\n",
  "1.11 12 August 2005 [kate fissell]\n"
  "   - Kate's 0.2 release packaging, for sourceforge\n"
  "\n",
  "1.12 17 August 2005 [rickr] - comment (doxygen) updates\n"
  "   - updated comments for most functions (2 updates from Cinly Ooi)\n"
  "   - added nifti_type_and_names_match()\n"
  "\n",
  "1.12a 24 August 2005 [rickr] - remove all tabs from Clibs/*/*.[ch]\n",
  "1.12b 25 August 2005 [rickr] - changes by Hans Johnson\n",
  "1.13  25 August 2005 [rickr]\n",
  "   - finished changes by Hans for Insight\n"
  "   - added const in all appropraite parameter locations (30-40)\n"
  "     (any pointer referencing data that will not change)\n"
  "   - shortened all string constants below 509 character limit\n"
  "1.14  28 October 2005 [HJohnson]\n",
  "   - use nifti_set_filenames() in nifti_convert_nhdr2nim()\n"
  "1.15  02 November 2005 [rickr]\n",
  "   - added skip_blank_ext to nifti_global_options\n"
  "   - added nifti_set_skip_blank_ext(), to set option\n"
  "   - if skip_blank_ext and no extensions, do not read/write extender\n"
  "1.16 18 November 2005 [rickr]\n",
  "   - removed any test or access of dim[i], i>dim[0]\n"
  "   - do not set pixdim for collapsed dims to 1.0, leave them as they are\n"
  "   - added magic and dim[i] tests in nifti_hdr_looks_good()\n"
  "   - added 2 size_t casts\n"
  "1.17 22 November 2005 [rickr]\n",
  "   - in hdr->nim, for i > dim[0], pass 0 or 1, else set to 1\n"
  "1.18 02 March 2006 [rickr]\n",
  "   - in nifti_alloc_NBL_mem(), fixed nt=0 case from 1.17 change\n"
  "1.19 23 May 2006 [HJohnson,rickr]\n",
  "   - nifti_write_ascii_image(): free(hstr)\n"
  "   - nifti_copy_extensions(): clear num_ext and ext_list\n"
  "1.20 27 Jun 2006 [rickr]\n",
  "   - nifti_findhdrname(): fixed assign of efirst to match stated logic\n"
  "     (problem found by Atle Bjørnerud)\n"
  "1.21 05 Sep 2006 [rickr] update for nifticlib-0.4 release\n",
  "   - was reminded to actually add nifti_set_skip_blank_ext()\n"
  "   - init g_opts.skip_blank_ext to 0\n"
  "1.22 01 Jun 2007 nifticlib-0.5 release\n",
  "1.23 05 Jun 2007 nifti_add_exten_to_list: revert on failure, free old list\n"
  "1.24 07 Jun 2007 nifti_copy_extensions: use esize-8 for data size\n"
  "1.25 12 Jun 2007 [rickr] EMPTY_IMAGE creation\n",
  "   - added nifti_make_new_header() - to create from dims/dtype\n"
  "   - added nifti_make_new_nim() - to create from dims/dtype/fill\n"
  "   - added nifti_is_valid_datatype(), and more debug info\n",
  "1.26 27 Jul 2007 [rickr] handle single volumes > 2^31 bytes (but < 2^32)\n",
  "1.27 28 Jul 2007 [rickr] nim->nvox, NBL-bsize are now type size_t\n"
  "1.28 30 Jul 2007 [rickr] size_t updates\n"
  "1.29 08 Aug 2007 [rickr] for list, valid_nifti_brick_list requires 3 dims\n"
  "----------------------------------------------------------------------\n"
};
static char gni_version[] = "nifti library version 1.29 (08 Aug, 2007)";

/*! global nifti options structure */
static nifti_global_options g_opts = { 1, 0 };

/*---------------------------------------------------------------------------*/
/* prototypes for internal functions - not part of exported library          */

/* extension routines */
static int  nifti_read_extensions( nifti_image *nim, znzFile fp, int remain );
static int  nifti_read_next_extension( nifti1_extension * nex, nifti_image *nim,                                       int remain, znzFile fp );
static int  nifti_check_extension(nifti_image *nim, int size,int code, int rem);
static void update_nifti_image_for_brick_list(nifti_image * nim , int nbricks);
static int  nifti_add_exten_to_list(nifti1_extension *  new_ext,
                                    nifti1_extension ** list, int new_length);
static int  nifti_fill_extension(nifti1_extension * ext, const char * data,
                                 int len, int ecode);

/* NBL routines */
static int  nifti_load_NBL_bricks(nifti_image * nim , int * slist, int * sindex,                                  nifti_brick_list * NBL, znzFile fp );
static int  nifti_alloc_NBL_mem(  nifti_image * nim, int nbricks,
                                  nifti_brick_list * nbl);
static int  nifti_copynsort(int nbricks, const int *blist, int **slist,
                            int **sindex);

/* for nifti_read_collapsed_image: */
static int  rci_read_data(nifti_image *nim, int *pivots, int *prods, int nprods,
                  const int dims[], char *data, znzFile fp, size_t base_offset);
static int  rci_alloc_mem(void ** data, int prods[8], int nprods, int nbyper );
static int  make_pivot_list(nifti_image * nim, const int dims[], int pivots[],
                            int prods[], int * nprods );

/* misc */
static int   need_nhdr_swap    (short dim0, int hdrsize);
static int   print_hex_vals    (const char * data, int nbytes, FILE * fp);
static int   unescape_string   (char *str);  /* string utility functions */
static char *escapize_string   (const char *str);

/* internal I/O routines */
static znzFile nifti_image_load_prep( nifti_image *nim );
static int     has_ascii_header(znzFile fp);
/*---------------------------------------------------------------------------*/


/* for calling from some main program */

/*----------------------------------------------------------------------*/
/*! display the nifti library module history (via stdout)
*//*--------------------------------------------------------------------*/
void nifti_disp_lib_hist( void )
{
   int c, len = sizeof(gni_history)/sizeof(char *);
   for( c = 0; c < len; c++ )
       fputs(gni_history[c], stdout);
}

/*----------------------------------------------------------------------*/
/*! display the nifti library version (via stdout)
*//*--------------------------------------------------------------------*/
void nifti_disp_lib_version( void )
{
   printf("%s, compiled %s\n", gni_version, __DATE__);
}


/*----------------------------------------------------------------------*/
/*! nifti_image_read_bricks        - read nifti data as array of bricks
 *
 *                                   13 Dec 2004 [rickr]
 * 
 *  \param  hname    - filename of dataset to read (must be valid)
 *  \param  nbricks  - number of sub-bricks to read
 *                     (if blist is valid, nbricks must be > 0)
 *  \param  blist    - list of sub-bricks to read
 *                     (can be NULL; if NULL, read complete dataset)
 *  \param  NBL      - pointer to empty nifti_brick_list struct
 *                     (must be a valid pointer)
 *
 *  \return
 *     <br> nim      - same as nifti_image_read, but nim->data will be NULL
 *     <br> NBL      - filled with data
 *
 * By default, this function will read the nifti dataset and break the data
 * into a list of nt*nu*nv*nw sub-bricks, each having size nx*ny*nz elements.
 * That is to say, instead of reading the entire dataset as a single array,
 * break it up into sub-bricks, each of size nx*ny*nz elements.
 *
 * If 'blist' is valid, it is taken to be a list of sub-bricks, of length
 * 'nbricks'.  The data will still be separated into sub-bricks of size
 * nx*ny*nz elements, but now 'nbricks' sub-bricks will be returned, of the
 * caller's choosing via 'blist'.
 *
 * E.g. consider a dataset with 12 sub-bricks (numbered 0..11), and the
 * following code:
 *
 * <pre>
 * { nifti_brick_list   NB_orig, NB_select;
 *   nifti_image      * nim_orig, * nim_select;
 *   int                blist[5] = { 7, 0, 5, 5, 9 };
 *
 *   nim_orig   = nifti_image_read_bricks("myfile.nii", 0, NULL,  &NB_orig);
 *   nim_select = nifti_image_read_bricks("myfile.nii", 5, blist, &NB_select);
 * }
 * </pre>
 *
 * Here, nim_orig gets the entire dataset, where NB_orig.nbricks = 11.  But
 * nim_select has NB_select.nbricks = 5.
 *
 * Note that the first case is not quite the same as just calling the
 * nifti_image_read function, as here the data is separated into sub-bricks.
 *
 * Note that valid blist elements are in [0..nt*nu*nv*nw-1],
 * or written [ 0 .. (dim[4]*dim[5]*dim[6]*dim[7] - 1) ].
 *
 * Note that, as is the case with all of the reading functions, the
 * data will be allocated, read in, and properly byte-swapped, if
 * necessary.
 *
 * \sa nifti_image_load_bricks, nifti_free_NBL, valid_nifti_brick_list,
       nifti_image_read
*//*----------------------------------------------------------------------*/
nifti_image *nifti_image_read_bricks(const char * hname, int nbricks,
                                     const int * blist, nifti_brick_list * NBL)
{
   nifti_image * nim;

   if( !hname || !NBL ){
      fprintf(stderr,"** nifti_image_read_bricks: bad params (%p,%p)\n",
              hname, (void *)NBL);
      return NULL;
   }

   if( blist && nbricks <= 0 ){
      fprintf(stderr,"** nifti_image_read_bricks: bad nbricks, %d\n", nbricks);
      return NULL;
   }

   nim = nifti_image_read(hname, 0);  /* read header, but not data */

   if( !nim ) return NULL;   /* errors were already printed */

   /* if we fail, free image and return */
   if( nifti_image_load_bricks(nim, nbricks, blist, NBL) <= 0 ){
      nifti_image_free(nim);
      return NULL;
   }

   if( blist ) update_nifti_image_for_brick_list(nim, nbricks);

   return nim;
}


/*----------------------------------------------------------------------
 * update_nifti_image_for_brick_list  - update nifti_image
 *
 * When loading a specific brick list, the distinction between
 * nt, nu, nv and nw is lost.  So put everything in t, and set
 * dim[0] = 4.
 *----------------------------------------------------------------------*/
static void update_nifti_image_for_brick_list( nifti_image * nim , int nbricks )
{
   int ndim;

   if( g_opts.debug > 2 ){
      fprintf(stderr,"+d updating image dimensions for %d bricks in list\n",
              nbricks);
      fprintf(stderr,"   ndim = %d\n",nim->ndim);
      fprintf(stderr,"   nx,ny,nz,nt,nu,nv,nw: (%d,%d,%d,%d,%d,%d,%d)\n",
              nim->nx, nim->ny, nim->nz, nim->nt, nim->nu, nim->nv, nim->nw);
   }

   nim->nt = nbricks;
   nim->nu = nim->nv = nim->nw = 1;
   nim->dim[4] = nbricks;
   nim->dim[5] = nim->dim[6] = nim->dim[7] = 1;

   /* compute nvox                                                       */
   /* do not rely on dimensions above dim[0]         16 Nov 2005 [rickr] */
   for( nim->nvox = 1, ndim = 1; ndim <= nim->dim[0]; ndim++ )
      nim->nvox *= nim->dim[ndim];

   /* update the dimensions to 4 or lower */
   for( ndim = 4; (ndim > 1) && (nim->dim[ndim] <= 1); ndim-- )
       ;

   if( g_opts.debug > 2 ){
      fprintf(stderr,"+d ndim = %d -> %d\n",nim->ndim, ndim);
      fprintf(stderr," --> (%d,%d,%d,%d,%d,%d,%d)\n",
              nim->nx, nim->ny, nim->nz, nim->nt, nim->nu, nim->nv, nim->nw);
   }

   nim->dim[0] = nim->ndim = ndim;
}


/*----------------------------------------------------------------------*/
/*! nifti_update_dims_from_array  - update nx, ny, ... from nim->dim[]

    Fix all the dimension information, based on a new nim->dim[].

    Note: we assume that dim[0] will not increase.

    Check for updates to pixdim[], dx,...,  nx,..., nvox, ndim, dim[0].
*//*--------------------------------------------------------------------*/
int nifti_update_dims_from_array( nifti_image * nim )
{
   int c, ndim;

   if( !nim ){
      fprintf(stderr,"** update_dims: missing nim\n");
      return 1;
   }

   if( g_opts.debug > 2 ){
      fprintf(stderr,"+d updating image dimensions given nim->dim:");
      for( c = 0; c < 8; c++ ) fprintf(stderr," %d", nim->dim[c]);
      fputc('\n',stderr);
   }

   /* verify dim[0] first */
   if(nim->dim[0] < 1 || nim->dim[0] > 7){
      fprintf(stderr,"** invalid dim[0], dim[] = ");
      for( c = 0; c < 8; c++ ) fprintf(stderr," %d", nim->dim[c]);
      fputc('\n',stderr);
      return 1;
   }

   /* set nx, ny ..., dx, dy, ..., one by one */

   /* less than 1, set to 1, else copy */
   if(nim->dim[1] < 1) nim->nx = nim->dim[1] = 1;
   else                nim->nx = nim->dim[1];
   nim->dx = nim->pixdim[1];

   /* if undefined, or less than 1, set to 1 */
   if(nim->dim[0] < 2 || (nim->dim[0] >= 2 && nim->dim[2] < 1))
      nim->ny = nim->dim[2] = 1;
   else
      nim->ny = nim->dim[2];
   /* copy delta values, in any case */
   nim->dy = nim->pixdim[2];

   if(nim->dim[0] < 3 || (nim->dim[0] >= 3 && nim->dim[3] < 1))
      nim->nz = nim->dim[3] = 1;
   else /* just copy vals from arrays */
      nim->nz = nim->dim[3];
   nim->dz = nim->pixdim[3];

   if(nim->dim[0] < 4 || (nim->dim[0] >= 4 && nim->dim[4] < 1))
      nim->nt = nim->dim[4] = 1;
   else /* just copy vals from arrays */
      nim->nt = nim->dim[4];
   nim->dt = nim->pixdim[4];

   if(nim->dim[0] < 5 || (nim->dim[0] >= 5 && nim->dim[5] < 1))
      nim->nu = nim->dim[5] = 1;
   else /* just copy vals from arrays */
      nim->nu = nim->dim[5];
   nim->du = nim->pixdim[5];

   if(nim->dim[0] < 6 || (nim->dim[0] >= 6 && nim->dim[6] < 1))
      nim->nv = nim->dim[6] = 1;
   else /* just copy vals from arrays */
      nim->nv = nim->dim[6];
   nim->dv = nim->pixdim[6];

   if(nim->dim[0] < 7 || (nim->dim[0] >= 7 && nim->dim[7] < 1))
      nim->nw = nim->dim[7] = 1;
   else /* just copy vals from arrays */
      nim->nw = nim->dim[7];
   nim->dw = nim->pixdim[7];

   for( c = 1, nim->nvox = 1; c <= nim->dim[0]; c++ )
      nim->nvox *= nim->dim[c];

   /* compute ndim, assuming it can be no larger than the old one */
   for( ndim = nim->dim[0]; (ndim > 1) && (nim->dim[ndim] <= 1); ndim-- )
       ;

   if( g_opts.debug > 2 ){
      fprintf(stderr,"+d ndim = %d -> %d\n",nim->ndim, ndim);
      fprintf(stderr," --> (%d,%d,%d,%d,%d,%d,%d)\n",
              nim->nx, nim->ny, nim->nz, nim->nt, nim->nu, nim->nv, nim->nw);
   }

   nim->dim[0] = nim->ndim = ndim;

   return 0;
}


/*----------------------------------------------------------------------*/
/*! Load the image data from disk into an already-prepared image struct.
 *
 * \param    nim      - initialized nifti_image, without data
 * \param    nbricks  - the length of blist (must be 0 if blist is NULL)
 * \param    blist    - an array of xyz volume indices to read (can be NULL)
 * \param    NBL      - pointer to struct where resulting data will be stored
 *
 * If blist is NULL, read all sub-bricks.
 * 
 * \return the number of loaded bricks (NBL->nbricks),
 *    0 on failure, < 0 on error
 *
 * NOTE: it is likely that another function will copy the data pointers
 *       out of NBL, in which case the only pointer the calling function
 *       will want to free is NBL->bricks (not each NBL->bricks[i]).
*//*--------------------------------------------------------------------*/
int nifti_image_load_bricks( nifti_image * nim , int nbricks,
                             const int * blist, nifti_brick_list * NBL )
{
   int     * slist = NULL, * sindex = NULL, rv;
   znzFile   fp;

   /* we can have blist == NULL */
   if( !nim || !NBL ){
      fprintf(stderr,"** nifti_image_load_bricks, bad params (%p,%p)\n",
              (void *)nim, (void *)NBL);
      return -1;
   }

   if( blist && nbricks <= 0 ){
      if( g_opts.debug > 1 )
         fprintf(stderr,"-d load_bricks: received blist with nbricks = %d,"
                        "ignoring blist\n", nbricks);
      blist = NULL; /* pretend nothing was passed */
   }

   if( blist && ! valid_nifti_brick_list(nim, nbricks, blist, g_opts.debug>0) )
      return -1;

   /* for efficiency, let's read the file in order */
   if( blist && nifti_copynsort( nbricks, blist, &slist, &sindex ) != 0 )
      return -1;

   /* open the file and position the FILE pointer */
   fp = nifti_image_load_prep( nim );
   if( !fp ){
      if( g_opts.debug > 0 )
         fprintf(stderr,"** nifti_image_load_bricks, failed load_prep\n");
      if( blist ){ free(slist); free(sindex); }
      return -1;
   }

   /* this will flag to allocate defaults */
   if( !blist ) nbricks = 0;
   if( nifti_alloc_NBL_mem( nim, nbricks, NBL ) != 0 ){
      if( blist ){ free(slist); free(sindex); }
      znzclose(fp);
      return -1;
   }

   rv = nifti_load_NBL_bricks(nim, slist, sindex, NBL, fp);

   if( rv != 0 ){
      nifti_free_NBL( NBL );  /* failure! */
      NBL->nbricks = 0; /* repetative, but clear */
   }

   if( slist ){ free(slist); free(sindex); }

   znzclose(fp);

   return NBL->nbricks;
}


/*----------------------------------------------------------------------*/
/*! nifti_free_NBL      - free all pointers and clear structure
 *
 * note: this does not presume to free the structure pointer
*//*--------------------------------------------------------------------*/
void nifti_free_NBL( nifti_brick_list * NBL )
{
   int c;

   if( NBL->bricks ){
      for( c = 0; c < NBL->nbricks; c++ )
         if( NBL->bricks[c] ) free(NBL->bricks[c]);
      free(NBL->bricks);
      NBL->bricks = NULL;
   }

   NBL->bsize = NBL->nbricks = 0;
}


/*----------------------------------------------------------------------
 * nifti_load_NBL_bricks      - read the file data into the NBL struct
 *
 * return 0 on success, -1 on failure
 *----------------------------------------------------------------------*/
static int nifti_load_NBL_bricks( nifti_image * nim , int * slist, int * sindex,
                                  nifti_brick_list * NBL, znzFile fp )
{
   size_t oposn, fposn;      /* orig and current file positions */
   size_t rv;
   long   test;
   int    c;
   int    prev, isrc, idest; /* previous and current sub-brick, and new index */

   test = znztell(fp);  /* store current file position */
   if( test < 0 ){
      fprintf(stderr,"** load bricks: ztell failed??\n");
      return -1;
   }
   fposn = oposn = test;

   /* first, handle the default case, no passed blist */
   if( !slist ){
      for( c = 0; c < NBL->nbricks; c++ ) {
         rv = nifti_read_buffer(fp, NBL->bricks[c], NBL->bsize, nim);
         if( rv != NBL->bsize ){
            fprintf(stderr,"** load bricks: cannot read brick %d from '%s'\n",
                    c, nim->iname ? nim->iname : nim->fname);
            return -1;
         }
      }
      if( g_opts.debug > 1 )
         fprintf(stderr,"+d read %d default %u-byte bricks from file %s\n",
                 NBL->nbricks, (unsigned int)NBL->bsize,
                 nim->iname ? nim->iname:nim->fname );
      return 0;
   }

   if( !sindex ){
      fprintf(stderr,"** load_NBL_bricks: missing index list\n");
      return -1;
   }

   prev = -1;   /* use prev for previous sub-brick */
   for( c = 0; c < NBL->nbricks; c++ ){
       isrc = slist[c];   /* this is original brick index (c is new one) */
       idest = sindex[c]; /* this is the destination index for this data */

       /* if this sub-brick is not the previous, we must read from disk */
       if( isrc != prev ){

          /* if we are not looking at the correct sub-brick, scan forward */
          if( fposn != (oposn + isrc*NBL->bsize) ){
             fposn = oposn + isrc*NBL->bsize;
             if( znzseek(fp, fposn, SEEK_SET) < 0 ){
                fprintf(stderr,"** failed to locate brick %d in file '%s'\n",
                        isrc, nim->iname ? nim->iname : nim->fname);
                return -1;
             }
          }

          /* only 10,000 lines later and we're actually reading something! */
          rv = nifti_read_buffer(fp, NBL->bricks[idest], NBL->bsize, nim);
          if( rv != NBL->bsize ){
             fprintf(stderr,"** failed to read brick %d from file '%s'\n",
                     isrc, nim->iname ? nim->iname : nim->fname);
             return -1;
          }
          fposn += NBL->bsize;
       } else {
          /* we have already read this sub-brick, just copy the previous one */
          /* note that this works because they are sorted */
          memcpy(NBL->bricks[idest], NBL->bricks[sindex[c-1]], NBL->bsize);
       }

       prev = isrc;  /* in any case, note the now previous sub-brick */
   }

   return 0;
}


/*----------------------------------------------------------------------
 * nifti_alloc_NBL_mem      - allocate memory for bricks
 *
 * return 0 on success, -1 on failure
 *----------------------------------------------------------------------*/
static int nifti_alloc_NBL_mem(nifti_image * nim, int nbricks,
                               nifti_brick_list * nbl)
{
   int c;

   /* if nbricks is not specified, use the default */
   if( nbricks > 0 ) nbl->nbricks = nbricks;
   else {  /* I missed this one with the 1.17 change    02 Mar 2006 [rickr] */
      nbl->nbricks = 1;
      for( c = 4; c <= nim->ndim; c++ )
          nbl->nbricks *= nim->dim[c];
   }

   nbl->bsize   = (size_t)nim->nx * nim->ny * nim->nz * nim->nbyper;/* bytes */
   nbl->bricks  = (void **)malloc(nbl->nbricks * sizeof(void *));

   if( ! nbl->bricks ){
      fprintf(stderr,"** NANM: failed to alloc %d void ptrs\n",nbricks);
      return -1;
   }

   for( c = 0; c < nbl->nbricks; c++ ){
      nbl->bricks[c] = (void *)malloc(nbl->bsize);
      if( ! nbl->bricks[c] ){
         fprintf(stderr,"** NANM: failed to alloc %u bytes for brick %d\n",
                 (unsigned int)nbl->bsize, c);
         /* so free and clear everything before returning */
         while( c > 0 ){
            c--;
            free(nbl->bricks[c]);
         }
         free(nbl->bricks);
         nbl->bricks = NULL;
         nbl->bsize = nbl->nbricks = 0;
         return -1;
      }
   }

   if( g_opts.debug > 2 )
      fprintf(stderr,"+d NANM: alloc'd %d bricks of %u bytes for NBL\n",
              nbl->nbricks, (unsigned int)nbl->bsize);

   return 0;
}


/*----------------------------------------------------------------------
 * nifti_copynsort      - copy int list, and sort with indices
 *
 * 1. duplicate the incoming list
 * 2. create an sindex list, and init with 0..nbricks-1
 * 3. do a slow insertion sort on the small slist, along with sindex list
 * 4. check results, just to be positive
 * 
 * So slist is sorted, and sindex hold original positions.
 *
 * return 0 on success, -1 on failure
 *----------------------------------------------------------------------*/
static int nifti_copynsort(int nbricks, const int * blist, int ** slist,
                           int ** sindex)
{
   int * stmp, * itmp;   /* for ease of typing/reading */
   int   c1, c2, spos, tmp;

   *slist  = (int *)malloc(nbricks * sizeof(int));
   *sindex = (int *)malloc(nbricks * sizeof(int));

   if( !*slist || !*sindex ){
      fprintf(stderr,"** NCS: failed to alloc %d ints for sorting\n",nbricks);
      if(*slist)  free(*slist);   /* maybe one succeeded */
      if(*sindex) free(*sindex);
      return -1;
   }

   /* init the lists */
   memcpy(*slist, blist, nbricks*sizeof(int));
   for( c1 = 0; c1 < nbricks; c1++ ) (*sindex)[c1] = c1;

   /* now actually sort slist */
   stmp = *slist;
   itmp = *sindex;
   for( c1 = 0; c1 < nbricks-1; c1++ ) {
      /* find smallest value, init to current */
      spos = c1;
      for( c2 = c1+1; c2 < nbricks; c2++ )
         if( stmp[c2] < stmp[spos] ) spos = c2;
      if( spos != c1 ) /* swap: fine, don't maintain sub-order, see if I care */
      {
         tmp        = stmp[c1];      /* first swap the sorting values */
         stmp[c1]   = stmp[spos];
         stmp[spos] = tmp;

         tmp        = itmp[c1];      /* then swap the index values */
         itmp[c1]   = itmp[spos];
         itmp[spos] = tmp;
      }
   }

   if( g_opts.debug > 2 ){
      fprintf(stderr,  "+d sorted indexing list:\n");
      fprintf(stderr,  "  orig   : ");
      for( c1 = 0; c1 < nbricks; c1++ ) fprintf(stderr,"  %d",blist[c1]);
      fprintf(stderr,"\n  new    : ");
      for( c1 = 0; c1 < nbricks; c1++ ) fprintf(stderr,"  %d",stmp[c1]);
      fprintf(stderr,"\n  indices: ");
      for( c1 = 0; c1 < nbricks; c1++ ) fprintf(stderr,"  %d",itmp[c1]);
      fputc('\n', stderr);
   }

   /* check the sort (why not?  I've got time...) */
   for( c1 = 0; c1 < nbricks-1; c1++ ){
       if( (stmp[c1] > stmp[c1+1]) || (blist[itmp[c1]] != stmp[c1]) ){
          fprintf(stderr,"** sorting screw-up, way to go, rick!\n");
          free(stmp); free(itmp); *slist = NULL; *sindex = NULL;
          return -1;
       }
   }

   if( g_opts.debug > 2 ) fprintf(stderr,"-d sorting is okay\n");

   return 0;
}


/*----------------------------------------------------------------------*/
/*! valid_nifti_brick_list      - check sub-brick list for image
 *
 * This function verifies that nbricks and blist are appropriate
 * for use with this nim, based on the dimensions.
 *
 * \param nim        nifti_image to check against
 * \param nbricks    number of brick indices in blist
 * \param blist      list of brick indices to check in nim
 * \param disp_error if this flag is set, report errors to user
 *
 * \return 1 if valid, 0 if not
*//*--------------------------------------------------------------------*/
int valid_nifti_brick_list(nifti_image * nim , int nbricks,
                           const int * blist, int disp_error)
{
   int c, nsubs;

   if( !nim ){
      if( disp_error || g_opts.debug > 0 )
         fprintf(stderr,"** valid_nifti_brick_list: missing nifti image\n");
      return 0;
   }

   if( nbricks <= 0 || !blist ){
      if( disp_error || g_opts.debug > 1 )
         fprintf(stderr,"** valid_nifti_brick_list: no brick list to check\n");
      return 0;
   }

   if( nim->dim[0] < 3 ){
      if( disp_error || g_opts.debug > 1 )
         fprintf(stderr,"** cannot read explict brick list from %d-D dataset\n",
                 nim->dim[0]);
      return 0;
   }

   /* nsubs sub-brick is nt*nu*nv*nw */
   for( c = 4, nsubs = 1; c <= nim->dim[0]; c++ )
      nsubs *= nim->dim[c];

   if( nsubs <= 0 ){
      fprintf(stderr,"** VNBL warning: bad dim list (%d,%d,%d,%d)\n",
                     nim->dim[4], nim->dim[5], nim->dim[6], nim->dim[7]);
      return 0;
   }

   for( c = 0; c < nbricks; c++ )
      if( (blist[c] < 0) || (blist[c] >= nsubs) ){
         if( disp_error || g_opts.debug > 1 )
            fprintf(stderr,
               "-d ** bad sub-brick chooser %d (#%d), valid range is [0,%d]\n",
               blist[c], c, nsubs-1);
         return 0;
      }

   return 1;  /* all is well */
}

#if 0
/* set any non-positive values to 1 */
static int int_force_positive( int * list, int nel )
{
   int c;
   if( !list || nel < 0 ){
      if( g_opts.debug > 0 )
         fprintf(stderr,"** int_force_positive: bad params (%p,%d)\n",
                 (void *)list,nel);
      return -1;
   }
   for( c = 0; c < nel; c++ )
      if( list[c] <= 0 ) list[c] = 1;
   return 0;
}
#endif
/* end of new nifti_image_read_bricks() functionality */

/*----------------------------------------------------------------------*/
/*! display the orientation from the quaternian fields
 *
 * \param mesg if non-NULL, display this message first
 * \param mat  the matrix to convert to "nearest" orientation
 *
 * \return -1 if results cannot be determined, 0 if okay
*//*--------------------------------------------------------------------*/
int nifti_disp_matrix_orient( const char * mesg, mat44 mat )
{
   int i, j, k;

   if ( mesg ) fputs( mesg, stderr );  /* use stdout? */

   nifti_mat44_to_orientation( mat, &i,&j,&k );
   if ( i <= 0 || j <= 0 || k <= 0 ) return -1;

   /* so we have good codes */
   fprintf(stderr, "  i orientation = '%s'\n"
                   "  j orientation = '%s'\n"
                   "  k orientation = '%s'\n",
                   nifti_orientation_string(i),
                   nifti_orientation_string(j),
                   nifti_orientation_string(k) );
   return 0;
}


/*----------------------------------------------------------------------*/
/*! duplicate the given string (alloc length+1)
 *
 * \return allocated pointer (or NULL on failure)
*//*--------------------------------------------------------------------*/
char *nifti_strdup(const char *str)
{
  char *dup= (char *)malloc( strlen(str)+1 );
  if (dup) strcpy(dup,str);
  return dup;
}


/*---------------------------------------------------------------------------*/
/*! Return a pointer to a string holding the name of a NIFTI datatype.

    \param dt NIfTI-1 datatype

    \return pointer to static string holding the datatype name

    \warning Do not free() or modify this string!
             It points to static storage.

    \sa NIFTI1_DATATYPES group in nifti1.h
*//*-------------------------------------------------------------------------*/
char *nifti_datatype_string( int dt )
{
   switch( dt ){
     case DT_UNKNOWN:    return "UNKNOWN"    ;
     case DT_BINARY:     return "BINARY"     ;
     case DT_INT8:       return "INT8"       ;
     case DT_UINT8:      return "UINT8"      ;
     case DT_INT16:      return "INT16"      ;
     case DT_UINT16:     return "UINT16"     ;
     case DT_INT32:      return "INT32"      ;
     case DT_UINT32:     return "UINT32"     ;
     case DT_INT64:      return "INT64"      ;
     case DT_UINT64:     return "UINT64"     ;
     case DT_FLOAT32:    return "FLOAT32"    ;
     case DT_FLOAT64:    return "FLOAT64"    ;
     case DT_FLOAT128:   return "FLOAT128"   ;
     case DT_COMPLEX64:  return "COMPLEX64"  ;
     case DT_COMPLEX128: return "COMPLEX128" ;
     case DT_COMPLEX256: return "COMPLEX256" ;
     case DT_RGB24:      return "RGB24"      ;
   }
   return "**ILLEGAL**" ;
}

/*----------------------------------------------------------------------*/
/*! Determine if the datatype code dt is an integer type (1=YES, 0=NO).

    \return whether the given NIfTI-1 datatype code is valid

    \sa     NIFTI1_DATATYPES group in nifti1.h
*//*--------------------------------------------------------------------*/
int nifti_is_inttype( int dt )
{
   switch( dt ){
     case DT_UNKNOWN:    return 0 ;
     case DT_BINARY:     return 0 ;
     case DT_INT8:       return 1 ;
     case DT_UINT8:      return 1 ;
     case DT_INT16:      return 1 ;
     case DT_UINT16:     return 1 ;
     case DT_INT32:      return 1 ;
     case DT_UINT32:     return 1 ;
     case DT_INT64:      return 1 ;
     case DT_UINT64:     return 1 ;
     case DT_FLOAT32:    return 0 ;
     case DT_FLOAT64:    return 0 ;
     case DT_FLOAT128:   return 0 ;
     case DT_COMPLEX64:  return 0 ;
     case DT_COMPLEX128: return 0 ;
     case DT_COMPLEX256: return 0 ;
     case DT_RGB24:      return 1 ;
   }
   return 0 ;
}

/*---------------------------------------------------------------------------*/
/*! Return a pointer to a string holding the name of a NIFTI units type.

    \param  uu NIfTI-1 unit code

    \return pointer to static string for the given unit type

    \warning Do not free() or modify this string!
             It points to static storage.

    \sa     NIFTI1_UNITS group in nifti1.h
*//*-------------------------------------------------------------------------*/
char *nifti_units_string( int uu )
{
   switch( uu ){
     case NIFTI_UNITS_METER:  return "m" ;
     case NIFTI_UNITS_MM:     return "mm" ;
     case NIFTI_UNITS_MICRON: return "um" ;
     case NIFTI_UNITS_SEC:    return "s" ;
     case NIFTI_UNITS_MSEC:   return "ms" ;
     case NIFTI_UNITS_USEC:   return "us" ;
     case NIFTI_UNITS_HZ:     return "Hz" ;
     case NIFTI_UNITS_PPM:    return "ppm" ;
     case NIFTI_UNITS_RADS:   return "rad/s" ;
   }
   return "Unknown" ;
}

/*---------------------------------------------------------------------------*/
/*! Return a pointer to a string holding the name of a NIFTI transform type.

    \param  xx NIfTI-1 xform code

    \return pointer to static string describing xform code

    \warning Do not free() or modify this string!
             It points to static storage.

    \sa     NIFTI1_XFORM_CODES group in nifti1.h
*//*-------------------------------------------------------------------------*/
char *nifti_xform_string( int xx )
{
   switch( xx ){
     case NIFTI_XFORM_SCANNER_ANAT:  return "Scanner Anat" ;
     case NIFTI_XFORM_ALIGNED_ANAT:  return "Aligned Anat" ;
     case NIFTI_XFORM_TALAIRACH:     return "Talairach" ;
     case NIFTI_XFORM_MNI_152:       return "MNI_152" ;
   }
   return "Unknown" ;
}

/*---------------------------------------------------------------------------*/
/*! Return a pointer to a string holding the name of a NIFTI intent type.

    \param  ii NIfTI-1 intent code

    \return pointer to static string describing code

    \warning Do not free() or modify this string!
             It points to static storage.

    \sa     NIFTI1_INTENT_CODES group in nifti1.h
*//*-------------------------------------------------------------------------*/
char *nifti_intent_string( int ii )
{
   switch( ii ){
     case NIFTI_INTENT_CORREL:     return "Correlation statistic" ;
     case NIFTI_INTENT_TTEST:      return "T-statistic" ;
     case NIFTI_INTENT_FTEST:      return "F-statistic" ;
     case NIFTI_INTENT_ZSCORE:     return "Z-score"     ;
     case NIFTI_INTENT_CHISQ:      return "Chi-squared distribution" ;
     case NIFTI_INTENT_BETA:       return "Beta distribution" ;
     case NIFTI_INTENT_BINOM:      return "Binomial distribution" ;
     case NIFTI_INTENT_GAMMA:      return "Gamma distribution" ;
     case NIFTI_INTENT_POISSON:    return "Poisson distribution" ;
     case NIFTI_INTENT_NORMAL:     return "Normal distribution" ;
     case NIFTI_INTENT_FTEST_NONC: return "F-statistic noncentral" ;
     case NIFTI_INTENT_CHISQ_NONC: return "Chi-squared noncentral" ;
     case NIFTI_INTENT_LOGISTIC:   return "Logistic distribution" ;
     case NIFTI_INTENT_LAPLACE:    return "Laplace distribution" ;
     case NIFTI_INTENT_UNIFORM:    return "Uniform distribition" ;
     case NIFTI_INTENT_TTEST_NONC: return "T-statistic noncentral" ;
     case NIFTI_INTENT_WEIBULL:    return "Weibull distribution" ;
     case NIFTI_INTENT_CHI:        return "Chi distribution" ;
     case NIFTI_INTENT_INVGAUSS:   return "Inverse Gaussian distribution" ;
     case NIFTI_INTENT_EXTVAL:     return "Extreme Value distribution" ;
     case NIFTI_INTENT_PVAL:       return "P-value" ;

     case NIFTI_INTENT_LOGPVAL:    return "Log P-value" ;
     case NIFTI_INTENT_LOG10PVAL:  return "Log10 P-value" ;

     case NIFTI_INTENT_ESTIMATE:   return "Estimate" ;
     case NIFTI_INTENT_LABEL:      return "Label index" ;
     case NIFTI_INTENT_NEURONAME:  return "NeuroNames index" ;
     case NIFTI_INTENT_GENMATRIX:  return "General matrix" ;
     case NIFTI_INTENT_SYMMATRIX:  return "Symmetric matrix" ;
     case NIFTI_INTENT_DISPVECT:   return "Displacement vector" ;
     case NIFTI_INTENT_VECTOR:     return "Vector" ;
     case NIFTI_INTENT_POINTSET:   return "Pointset" ;
     case NIFTI_INTENT_TRIANGLE:   return "Triangle" ;
     case NIFTI_INTENT_QUATERNION: return "Quaternion" ;

     case NIFTI_INTENT_DIMLESS:    return "Dimensionless number" ;
   }
   return "Unknown" ;
}

/*---------------------------------------------------------------------------*/
/*! Return a pointer to a string holding the name of a NIFTI slice_code.

    \param  ss NIfTI-1 slice order code

    \return pointer to static string describing code

    \warning Do not free() or modify this string!
             It points to static storage.

    \sa     NIFTI1_SLICE_ORDER group in nifti1.h
*//*-------------------------------------------------------------------------*/
char *nifti_slice_string( int ss )
{
   switch( ss ){
     case NIFTI_SLICE_SEQ_INC:  return "sequential_increasing"    ;
     case NIFTI_SLICE_SEQ_DEC:  return "sequential_decreasing"    ;
     case NIFTI_SLICE_ALT_INC:  return "alternating_increasing"   ;
     case NIFTI_SLICE_ALT_DEC:  return "alternating_decreasing"   ;
     case NIFTI_SLICE_ALT_INC2: return "alternating_increasing_2" ;
     case NIFTI_SLICE_ALT_DEC2: return "alternating_decreasing_2" ;
   }
   return "Unknown" ;
}

/*---------------------------------------------------------------------------*/
/*! Return a pointer to a string holding the name of a NIFTI orientation.

    \param ii orientation code

    \return pointer to static string holding the orientation information

    \warning Do not free() or modify the return string!
             It points to static storage.

    \sa  NIFTI_L2R in nifti1_io.h
*//*-------------------------------------------------------------------------*/
char *nifti_orientation_string( int ii )
{
   switch( ii ){
     case NIFTI_L2R: return "Left-to-Right" ;
     case NIFTI_R2L: return "Right-to-Left" ;
     case NIFTI_P2A: return "Posterior-to-Anterior" ;
     case NIFTI_A2P: return "Anterior-to-Posterior" ;
     case NIFTI_I2S: return "Inferior-to-Superior" ;
     case NIFTI_S2I: return "Superior-to-Inferior" ;
   }
   return "Unknown" ;
}

/*--------------------------------------------------------------------------*/
/*! Given a datatype code, set number of bytes per voxel and the swapsize.

    \param datatype nifti1 datatype code
    \param nbyper   pointer to return value: number of bytes per voxel
    \param swapsize pointer to return value: size of swap blocks

    \return appropriate values at nbyper and swapsize

    The swapsize is set to 0 if this datatype doesn't ever need swapping.

    \sa NIFTI1_DATATYPES in nifti1.h
*//*------------------------------------------------------------------------*/
void nifti_datatype_sizes( int datatype , int *nbyper, int *swapsize )
{
   int nb=0, ss=0 ;
   switch( datatype ){
     case DT_INT8:
     case DT_UINT8:       nb =  1 ; ss =  0 ; break ;

     case DT_INT16:
     case DT_UINT16:      nb =  2 ; ss =  2 ; break ;

     case DT_RGB24:       nb =  3 ; ss =  0 ; break ;

     case DT_INT32:
     case DT_UINT32:
     case DT_FLOAT32:     nb =  4 ; ss =  4 ; break ;

     case DT_COMPLEX64:   nb =  8 ; ss =  4 ; break ;

     case DT_FLOAT64:
     case DT_INT64:
     case DT_UINT64:      nb =  8 ; ss =  8 ; break ;

     case DT_FLOAT128:    nb = 16 ; ss = 16 ; break ;

     case DT_COMPLEX128:  nb = 16 ; ss =  8 ; break ;

     case DT_COMPLEX256:  nb = 32 ; ss = 16 ; break ;
   }

   ASSIF(nbyper,nb) ; ASSIF(swapsize,ss) ; return ;
}

/*---------------------------------------------------------------------------*/
/*! Given the quaternion parameters (etc.), compute a transformation matrix.

   See comments in nifti1.h for details.
     - qb,qc,qd = quaternion parameters
     - qx,qy,qz = offset parameters
     - dx,dy,dz = grid stepsizes (non-negative inputs are set to 1.0)
     - qfac     = sign of dz step (< 0 is negative; >= 0 is positive)

   <pre>
   If qx=qy=qz=0, dx=dy=dz=1, then the output is a rotation matrix.
   For qfac >= 0, the rotation is proper.
   For qfac <  0, the rotation is improper.
   </pre>

   \see "QUATERNION REPRESENTATION OF ROTATION MATRIX" in nifti1.h
   \see nifti_mat44_to_quatern, nifti_make_orthog_mat44,
       nifti_mat44_to_orientation

*//*-------------------------------------------------------------------------*/
mat44 nifti_quatern_to_mat44( float qb, float qc, float qd,
                              float qx, float qy, float qz,
                              float dx, float dy, float dz, float qfac )
{
   mat44 R ;
   double a,b=qb,c=qc,d=qd , xd,yd,zd ;

   /* last row is always [ 0 0 0 1 ] */

   R.m[3][0]=R.m[3][1]=R.m[3][2] = 0.0 ; R.m[3][3]= 1.0 ;

   /* compute a parameter from b,c,d */

   a = 1.0l - (b*b + c*c + d*d) ;
   if( a < 1.e-7l ){                   /* special case */
     a = 1.0l / sqrt(b*b+c*c+d*d) ;
     b *= a ; c *= a ; d *= a ;        /* normalize (b,c,d) vector */
     a = 0.0l ;                        /* a = 0 ==> 180 degree rotation */
   } else{
     a = sqrt(a) ;                     /* angle = 2*arccos(a) */
   }

   /* load rotation matrix, including scaling factors for voxel sizes */

   xd = (dx > 0.0) ? dx : 1.0l ;       /* make sure are positive */
   yd = (dy > 0.0) ? dy : 1.0l ;
   zd = (dz > 0.0) ? dz : 1.0l ;

   if( qfac < 0.0 ) zd = -zd ;         /* left handedness? */

   R.m[0][0] =        (a*a+b*b-c*c-d*d) * xd ;
   R.m[0][1] = 2.0l * (b*c-a*d        ) * yd ;
   R.m[0][2] = 2.0l * (b*d+a*c        ) * zd ;
   R.m[1][0] = 2.0l * (b*c+a*d        ) * xd ;
   R.m[1][1] =        (a*a+c*c-b*b-d*d) * yd ;
   R.m[1][2] = 2.0l * (c*d-a*b        ) * zd ;
   R.m[2][0] = 2.0l * (b*d-a*c        ) * xd ;
   R.m[2][1] = 2.0l * (c*d+a*b        ) * yd ;
   R.m[2][2] =        (a*a+d*d-c*c-b*b) * zd ;

   /* load offsets */

   R.m[0][3] = qx ; R.m[1][3] = qy ; R.m[2][3] = qz ;

   return R ;
}

/*---------------------------------------------------------------------------*/
/*! Given the 3x4 upper corner of the matrix R, compute the quaternion
   parameters that fit it.

     - Any NULL pointer on input won't get assigned (e.g., if you don't want
       dx,dy,dz, just pass NULL in for those pointers).
     - If the 3 input matrix columns are NOT orthogonal, they will be
       orthogonalized prior to calculating the parameters, using
       the polar decomposition to find the orthogonal matrix closest
       to the column-normalized input matrix.
     - However, if the 3 input matrix columns are NOT orthogonal, then
       the matrix produced by nifti_quatern_to_mat44 WILL have orthogonal
       columns, so it won't be the same as the matrix input here.
       This "feature" is because the NIFTI 'qform' transform is
       deliberately not fully general -- it is intended to model a volume
       with perpendicular axes.
     - If the 3 input matrix columns are not even linearly independent,
       you'll just have to take your luck, won't you?

   \see "QUATERNION REPRESENTATION OF ROTATION MATRIX" in nifti1.h

   \see nifti_quatern_to_mat44, nifti_make_orthog_mat44,
       nifti_mat44_to_orientation
*//*-------------------------------------------------------------------------*/
void nifti_mat44_to_quatern( mat44 R ,
                             float *qb, float *qc, float *qd,
                             float *qx, float *qy, float *qz,
                             float *dx, float *dy, float *dz, float *qfac )
{
   double r11,r12,r13 , r21,r22,r23 , r31,r32,r33 ;
   double xd,yd,zd , a,b,c,d ;
   mat33 P,Q ;

   /* offset outputs are read write out of input matrix  */

   ASSIF(qx,R.m[0][3]) ; ASSIF(qy,R.m[1][3]) ; ASSIF(qz,R.m[2][3]) ;

   /* load 3x3 matrix into local variables */

   r11 = R.m[0][0] ; r12 = R.m[0][1] ; r13 = R.m[0][2] ;
   r21 = R.m[1][0] ; r22 = R.m[1][1] ; r23 = R.m[1][2] ;
   r31 = R.m[2][0] ; r32 = R.m[2][1] ; r33 = R.m[2][2] ;

   /* compute lengths of each column; these determine grid spacings  */

   xd = sqrt( r11*r11 + r21*r21 + r31*r31 ) ;
   yd = sqrt( r12*r12 + r22*r22 + r32*r32 ) ;
   zd = sqrt( r13*r13 + r23*r23 + r33*r33 ) ;

   /* if a column length is zero, patch the trouble */

   if( xd == 0.0l ){ r11 = 1.0l ; r21 = r31 = 0.0l ; xd = 1.0l ; }
   if( yd == 0.0l ){ r22 = 1.0l ; r12 = r32 = 0.0l ; yd = 1.0l ; }
   if( zd == 0.0l ){ r33 = 1.0l ; r13 = r23 = 0.0l ; zd = 1.0l ; }

   /* assign the output lengths */

   ASSIF(dx,xd) ; ASSIF(dy,yd) ; ASSIF(dz,zd) ;

   /* normalize the columns */

   r11 /= xd ; r21 /= xd ; r31 /= xd ;
   r12 /= yd ; r22 /= yd ; r32 /= yd ;
   r13 /= zd ; r23 /= zd ; r33 /= zd ;

   /* At this point, the matrix has normal columns, but we have to allow
      for the fact that the hideous user may not have given us a matrix
      with orthogonal columns.

      So, now find the orthogonal matrix closest to the current matrix.

      One reason for using the polar decomposition to get this
      orthogonal matrix, rather than just directly orthogonalizing
      the columns, is so that inputting the inverse matrix to R
      will result in the inverse orthogonal matrix at this point.
      If we just orthogonalized the columns, this wouldn't necessarily hold. */

   Q.m[0][0] = r11 ; Q.m[0][1] = r12 ; Q.m[0][2] = r13 ; /* load Q */
   Q.m[1][0] = r21 ; Q.m[1][1] = r22 ; Q.m[1][2] = r23 ;
   Q.m[2][0] = r31 ; Q.m[2][1] = r32 ; Q.m[2][2] = r33 ;

   P = nifti_mat33_polar(Q) ;  /* P is orthog matrix closest to Q */

   r11 = P.m[0][0] ; r12 = P.m[0][1] ; r13 = P.m[0][2] ; /* unload */
   r21 = P.m[1][0] ; r22 = P.m[1][1] ; r23 = P.m[1][2] ;
   r31 = P.m[2][0] ; r32 = P.m[2][1] ; r33 = P.m[2][2] ;

   /*                            [ r11 r12 r13 ]               */
   /* at this point, the matrix  [ r21 r22 r23 ] is orthogonal */
   /*                            [ r31 r32 r33 ]               */

   /* compute the determinant to determine if it is proper */

   zd = r11*r22*r33-r11*r32*r23-r21*r12*r33
       +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;  /* should be -1 or 1 */

   if( zd > 0 ){             /* proper */
     ASSIF(qfac,1.0) ;
   } else {                  /* improper ==> flip 3rd column */
     ASSIF(qfac,-1.0) ;
     r13 = -r13 ; r23 = -r23 ; r33 = -r33 ;
   }

   /* now, compute quaternion parameters */

   a = r11 + r22 + r33 + 1.0l ;

   if( a > 0.5l ){                /* simplest case */
     a = 0.5l * sqrt(a) ;
     b = 0.25l * (r32-r23) / a ;
     c = 0.25l * (r13-r31) / a ;
     d = 0.25l * (r21-r12) / a ;
   } else {                       /* trickier case */
     xd = 1.0 + r11 - (r22+r33) ;  /* 4*b*b */
     yd = 1.0 + r22 - (r11+r33) ;  /* 4*c*c */
     zd = 1.0 + r33 - (r11+r22) ;  /* 4*d*d */
     if( xd > 1.0 ){
       b = 0.5l * sqrt(xd) ;
       c = 0.25l* (r12+r21) / b ;
       d = 0.25l* (r13+r31) / b ;
       a = 0.25l* (r32-r23) / b ;
     } else if( yd > 1.0 ){
       c = 0.5l * sqrt(yd) ;
       b = 0.25l* (r12+r21) / c ;
       d = 0.25l* (r23+r32) / c ;
       a = 0.25l* (r13-r31) / c ;
     } else {
       d = 0.5l * sqrt(zd) ;
       b = 0.25l* (r13+r31) / d ;
       c = 0.25l* (r23+r32) / d ;
       a = 0.25l* (r21-r12) / d ;
     }
     if( a < 0.0l ){ b=-b ; c=-c ; d=-d; a=-a; }
   }

   ASSIF(qb,b) ; ASSIF(qc,c) ; ASSIF(qd,d) ;
   return ;
}

/*---------------------------------------------------------------------------*/
/*! Compute the inverse of a bordered 4x4 matrix.

     <pre>
   - Some numerical code fragments were generated by Maple 8.
   - If a singular matrix is input, the output matrix will be all zero.
   - You can check for this by examining the [3][3] element, which will
     be 1.0 for the normal case and 0.0 for the bad case.

     The input matrix should have the form:
        [ r11 r12 r13 v1 ]
        [ r21 r22 r23 v2 ]
        [ r31 r32 r33 v3 ]
        [  0   0   0   1 ]
     </pre>
*//*-------------------------------------------------------------------------*/
mat44 nifti_mat44_inverse( mat44 R )
{
   double r11,r12,r13,r21,r22,r23,r31,r32,r33,v1,v2,v3 , deti ;
   mat44 Q ;
                                                       /*  INPUT MATRIX IS:  */
   r11 = R.m[0][0]; r12 = R.m[0][1]; r13 = R.m[0][2];  /* [ r11 r12 r13 v1 ] */
   r21 = R.m[1][0]; r22 = R.m[1][1]; r23 = R.m[1][2];  /* [ r21 r22 r23 v2 ] */
   r31 = R.m[2][0]; r32 = R.m[2][1]; r33 = R.m[2][2];  /* [ r31 r32 r33 v3 ] */
   v1  = R.m[0][3]; v2  = R.m[1][3]; v3  = R.m[2][3];  /* [  0   0   0   1 ] */

   deti = r11*r22*r33-r11*r32*r23-r21*r12*r33
         +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;

   if( deti != 0.0l ) deti = 1.0l / deti ;

   Q.m[0][0] = deti*( r22*r33-r32*r23) ;
   Q.m[0][1] = deti*(-r12*r33+r32*r13) ;
   Q.m[0][2] = deti*( r12*r23-r22*r13) ;
   Q.m[0][3] = deti*(-r12*r23*v3+r12*v2*r33+r22*r13*v3
                     -r22*v1*r33-r32*r13*v2+r32*v1*r23) ;

   Q.m[1][0] = deti*(-r21*r33+r31*r23) ;
   Q.m[1][1] = deti*( r11*r33-r31*r13) ;
   Q.m[1][2] = deti*(-r11*r23+r21*r13) ;
   Q.m[1][3] = deti*( r11*r23*v3-r11*v2*r33-r21*r13*v3
                     +r21*v1*r33+r31*r13*v2-r31*v1*r23) ;

   Q.m[2][0] = deti*( r21*r32-r31*r22) ;
   Q.m[2][1] = deti*(-r11*r32+r31*r12) ;
   Q.m[2][2] = deti*( r11*r22-r21*r12) ;
   Q.m[2][3] = deti*(-r11*r22*v3+r11*r32*v2+r21*r12*v3
                     -r21*r32*v1-r31*r12*v2+r31*r22*v1) ;

   Q.m[3][0] = Q.m[3][1] = Q.m[3][2] = 0.0l ;
   Q.m[3][3] = (deti == 0.0l) ? 0.0l : 1.0l ; /* failure flag if deti == 0 */

   return Q ;
}

/*---------------------------------------------------------------------------*/
/*! Input 9 floats and make an orthgonal mat44 out of them.

   Each row is normalized, then nifti_mat33_polar() is used to orthogonalize
   them.  If row #3 (r31,r32,r33) is input as zero, then it will be taken to
   be the cross product of rows #1 and #2.

   This function can be used to create a rotation matrix for transforming
   an oblique volume to anatomical coordinates.  For this application:
    - row #1 (r11,r12,r13) is the direction vector along the image i-axis
    - row #2 (r21,r22,r23) is the direction vector along the image j-axis
    - row #3 (r31,r32,r33) is the direction vector along the slice direction
      (if available; otherwise enter it as 0's)

   The first 2 rows can be taken from the DICOM attribute (0020,0037)
   "Image Orientation (Patient)".

   After forming the rotation matrix, the complete affine transformation from
   (i,j,k) grid indexes to (x,y,z) spatial coordinates can be computed by
   multiplying each column by the appropriate grid spacing:
    - column #1 (R.m[0][0],R.m[1][0],R.m[2][0]) by delta-x
    - column #2 (R.m[0][1],R.m[1][1],R.m[2][1]) by delta-y
    - column #3 (R.m[0][2],R.m[1][2],R.m[2][2]) by delta-z

   and by then placing the center (x,y,z) coordinates of voxel (0,0,0) into
   the column #4 (R.m[0][3],R.m[1][3],R.m[2][3]).

   \sa nifti_quatern_to_mat44, nifti_mat44_to_quatern,
       nifti_mat44_to_orientation
*//*-------------------------------------------------------------------------*/
mat44 nifti_make_orthog_mat44( float r11, float r12, float r13 ,
                               float r21, float r22, float r23 ,
                               float r31, float r32, float r33  )
{
   mat44 R ;
   mat33 Q , P ;
   double val ;

   R.m[3][0] = R.m[3][1] = R.m[3][2] = 0.0l ; R.m[3][3] = 1.0l ;

   Q.m[0][0] = r11 ; Q.m[0][1] = r12 ; Q.m[0][2] = r13 ; /* load Q */
   Q.m[1][0] = r21 ; Q.m[1][1] = r22 ; Q.m[1][2] = r23 ;
   Q.m[2][0] = r31 ; Q.m[2][1] = r32 ; Q.m[2][2] = r33 ;

   /* normalize row 1 */

   val = Q.m[0][0]*Q.m[0][0] + Q.m[0][1]*Q.m[0][1] + Q.m[0][2]*Q.m[0][2] ;
   if( val > 0.0l ){
     val = 1.0l / sqrt(val) ;
     Q.m[0][0] *= val ; Q.m[0][1] *= val ; Q.m[0][2] *= val ;
   } else {
     Q.m[0][0] = 1.0l ; Q.m[0][1] = 0.0l ; Q.m[0][2] = 0.0l ;
   }

   /* normalize row 2 */

   val = Q.m[1][0]*Q.m[1][0] + Q.m[1][1]*Q.m[1][1] + Q.m[1][2]*Q.m[1][2] ;
   if( val > 0.0l ){
     val = 1.0l / sqrt(val) ;
     Q.m[1][0] *= val ; Q.m[1][1] *= val ; Q.m[1][2] *= val ;
   } else {
     Q.m[1][0] = 0.0l ; Q.m[1][1] = 1.0l ; Q.m[1][2] = 0.0l ;
   }

   /* normalize row 3 */

   val = Q.m[2][0]*Q.m[2][0] + Q.m[2][1]*Q.m[2][1] + Q.m[2][2]*Q.m[2][2] ;
   if( val > 0.0l ){
     val = 1.0l / sqrt(val) ;
     Q.m[2][0] *= val ; Q.m[2][1] *= val ; Q.m[2][2] *= val ;
   } else {
     Q.m[2][0] = Q.m[0][1]*Q.m[1][2] - Q.m[0][2]*Q.m[1][1] ;  /* cross */
     Q.m[2][1] = Q.m[0][2]*Q.m[1][0] - Q.m[0][0]*Q.m[1][2] ;  /* product */
     Q.m[2][2] = Q.m[0][0]*Q.m[1][1] - Q.m[0][1]*Q.m[1][0] ;
   }

   P = nifti_mat33_polar(Q) ;  /* P is orthog matrix closest to Q */

   R.m[0][0] = P.m[0][0] ; R.m[0][1] = P.m[0][1] ; R.m[0][2] = P.m[0][2] ;
   R.m[1][0] = P.m[1][0] ; R.m[1][1] = P.m[1][1] ; R.m[1][2] = P.m[1][2] ;
   R.m[2][0] = P.m[2][0] ; R.m[2][1] = P.m[2][1] ; R.m[2][2] = P.m[2][2] ;

   R.m[0][3] = R.m[1][3] = R.m[2][3] = 0.0 ; return R ;
}

/*----------------------------------------------------------------------*/
/*! compute the inverse of a 3x3 matrix
*//*--------------------------------------------------------------------*/
mat33 nifti_mat33_inverse( mat33 R )   /* inverse of 3x3 matrix */
{
   double r11,r12,r13,r21,r22,r23,r31,r32,r33 , deti ;
   mat33 Q ;
                                                       /*  INPUT MATRIX:  */
   r11 = R.m[0][0]; r12 = R.m[0][1]; r13 = R.m[0][2];  /* [ r11 r12 r13 ] */
   r21 = R.m[1][0]; r22 = R.m[1][1]; r23 = R.m[1][2];  /* [ r21 r22 r23 ] */
   r31 = R.m[2][0]; r32 = R.m[2][1]; r33 = R.m[2][2];  /* [ r31 r32 r33 ] */

   deti = r11*r22*r33-r11*r32*r23-r21*r12*r33
         +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;

   if( deti != 0.0l ) deti = 1.0l / deti ;

   Q.m[0][0] = deti*( r22*r33-r32*r23) ;
   Q.m[0][1] = deti*(-r12*r33+r32*r13) ;
   Q.m[0][2] = deti*( r12*r23-r22*r13) ;

   Q.m[1][0] = deti*(-r21*r33+r31*r23) ;
   Q.m[1][1] = deti*( r11*r33-r31*r13) ;
   Q.m[1][2] = deti*(-r11*r23+r21*r13) ;

   Q.m[2][0] = deti*( r21*r32-r31*r22) ;
   Q.m[2][1] = deti*(-r11*r32+r31*r12) ;
   Q.m[2][2] = deti*( r11*r22-r21*r12) ;

   return Q ;
}

/*----------------------------------------------------------------------*/
/*! compute the determinant of a 3x3 matrix
*//*--------------------------------------------------------------------*/
float nifti_mat33_determ( mat33 R )   /* determinant of 3x3 matrix */
{
   double r11,r12,r13,r21,r22,r23,r31,r32,r33 ;
                                                       /*  INPUT MATRIX:  */
   r11 = R.m[0][0]; r12 = R.m[0][1]; r13 = R.m[0][2];  /* [ r11 r12 r13 ] */
   r21 = R.m[1][0]; r22 = R.m[1][1]; r23 = R.m[1][2];  /* [ r21 r22 r23 ] */
   r31 = R.m[2][0]; r32 = R.m[2][1]; r33 = R.m[2][2];  /* [ r31 r32 r33 ] */

   return r11*r22*r33-r11*r32*r23-r21*r12*r33
         +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;
}

/*----------------------------------------------------------------------*/
/*! compute the max row norm of a 3x3 matrix
*//*--------------------------------------------------------------------*/
float nifti_mat33_rownorm( mat33 A )  /* max row norm of 3x3 matrix */
{
   float r1,r2,r3 ;

   r1 = fabs(A.m[0][0])+fabs(A.m[0][1])+fabs(A.m[0][2]) ;
   r2 = fabs(A.m[1][0])+fabs(A.m[1][1])+fabs(A.m[1][2]) ;
   r3 = fabs(A.m[2][0])+fabs(A.m[2][1])+fabs(A.m[2][2]) ;
   if( r1 < r2 ) r1 = r2 ;
   if( r1 < r3 ) r1 = r3 ;
   return r1 ;
}

/*----------------------------------------------------------------------*/
/*! compute the max column norm of a 3x3 matrix
*//*--------------------------------------------------------------------*/
float nifti_mat33_colnorm( mat33 A )  /* max column norm of 3x3 matrix */
{
   float r1,r2,r3 ;

   r1 = fabs(A.m[0][0])+fabs(A.m[1][0])+fabs(A.m[2][0]) ;
   r2 = fabs(A.m[0][1])+fabs(A.m[1][1])+fabs(A.m[2][1]) ;
   r3 = fabs(A.m[0][2])+fabs(A.m[1][2])+fabs(A.m[2][2]) ;
   if( r1 < r2 ) r1 = r2 ;
   if( r1 < r3 ) r1 = r3 ;
   return r1 ;
}

/*----------------------------------------------------------------------*/
/*! multiply 2 3x3 matrices
*//*--------------------------------------------------------------------*/
mat33 nifti_mat33_mul( mat33 A , mat33 B )  /* multiply 2 3x3 matrices */
{
   mat33 C ; int i,j ;
   for( i=0 ; i < 3 ; i++ )
    for( j=0 ; j < 3 ; j++ )
      C.m[i][j] =  A.m[i][0] * B.m[0][j]
                 + A.m[i][1] * B.m[1][j]
                 + A.m[i][2] * B.m[2][j] ;
   return C ;
}

/*---------------------------------------------------------------------------*/
/*! polar decomposition of a 3x3 matrix

   This finds the closest orthogonal matrix to input A
   (in both the Frobenius and L2 norms).

   Algorithm is that from NJ Higham, SIAM J Sci Stat Comput, 7:1160-1174.
*//*-------------------------------------------------------------------------*/
mat33 nifti_mat33_polar( mat33 A )
{
   mat33 X , Y , Z ;
   float alp,bet,gam,gmi , dif=1.0 ;
   int k=0 ;

   X = A ;

   /* force matrix to be nonsingular */

   gam = nifti_mat33_determ(X) ;
   while( gam == 0.0 ){        /* perturb matrix */
     gam = 0.00001 * ( 0.001 + nifti_mat33_rownorm(X) ) ;
     X.m[0][0] += gam ; X.m[1][1] += gam ; X.m[2][2] += gam ;
     gam = nifti_mat33_determ(X) ;
   }

   while(1){
     Y = nifti_mat33_inverse(X) ;
     if( dif > 0.3 ){     /* far from convergence */
       alp = sqrt( nifti_mat33_rownorm(X) * nifti_mat33_colnorm(X) ) ;
       bet = sqrt( nifti_mat33_rownorm(Y) * nifti_mat33_colnorm(Y) ) ;
       gam = sqrt( bet / alp ) ;
       gmi = 1.0 / gam ;
     } else {
       gam = gmi = 1.0 ;  /* close to convergence */
     }
     Z.m[0][0] = 0.5 * ( gam*X.m[0][0] + gmi*Y.m[0][0] ) ;
     Z.m[0][1] = 0.5 * ( gam*X.m[0][1] + gmi*Y.m[1][0] ) ;
     Z.m[0][2] = 0.5 * ( gam*X.m[0][2] + gmi*Y.m[2][0] ) ;
     Z.m[1][0] = 0.5 * ( gam*X.m[1][0] + gmi*Y.m[0][1] ) ;
     Z.m[1][1] = 0.5 * ( gam*X.m[1][1] + gmi*Y.m[1][1] ) ;
     Z.m[1][2] = 0.5 * ( gam*X.m[1][2] + gmi*Y.m[2][1] ) ;
     Z.m[2][0] = 0.5 * ( gam*X.m[2][0] + gmi*Y.m[0][2] ) ;
     Z.m[2][1] = 0.5 * ( gam*X.m[2][1] + gmi*Y.m[1][2] ) ;
     Z.m[2][2] = 0.5 * ( gam*X.m[2][2] + gmi*Y.m[2][2] ) ;

     dif = fabs(Z.m[0][0]-X.m[0][0])+fabs(Z.m[0][1]-X.m[0][1])
          +fabs(Z.m[0][2]-X.m[0][2])+fabs(Z.m[1][0]-X.m[1][0])
          +fabs(Z.m[1][1]-X.m[1][1])+fabs(Z.m[1][2]-X.m[1][2])
          +fabs(Z.m[2][0]-X.m[2][0])+fabs(Z.m[2][1]-X.m[2][1])
          +fabs(Z.m[2][2]-X.m[2][2])                          ;

     k = k+1 ;
     if( k > 100 || dif < 3.e-6 ) break ;  /* convergence or exhaustion */
     X = Z ;
   }

   return Z ;
}

/*---------------------------------------------------------------------------*/
/*! compute the (closest) orientation from a 4x4 ijk->xyz tranformation matrix

   <pre>
   Input:  4x4 matrix that transforms (i,j,k) indexes to (x,y,z) coordinates,
           where +x=Right, +y=Anterior, +z=Superior.
           (Only the upper-left 3x3 corner of R is used herein.)
   Output: 3 orientation codes that correspond to the closest "standard"
           anatomical orientation of the (i,j,k) axes.
   Method: Find which permutation of (x,y,z) has the smallest angle to the
           (i,j,k) axes directions, which are the columns of the R matrix.
   Errors: The codes returned will be zero.

   For example, an axial volume might get return values of
     *icod = NIFTI_R2L   (i axis is mostly Right to Left)
     *jcod = NIFTI_P2A   (j axis is mostly Posterior to Anterior)
     *kcod = NIFTI_I2S   (k axis is mostly Inferior to Superior)
   </pre>

   \see "QUATERNION REPRESENTATION OF ROTATION MATRIX" in nifti1.h

   \see nifti_quatern_to_mat44, nifti_mat44_to_quatern,
        nifti_make_orthog_mat44
*//*-------------------------------------------------------------------------*/
void nifti_mat44_to_orientation( mat44 R , int *icod, int *jcod, int *kcod )
{
   float xi,xj,xk , yi,yj,yk , zi,zj,zk , val,detQ,detP ;
   mat33 P , Q , M ;
   int i,j,k=0,p,q,r , ibest,jbest,kbest,pbest,qbest,rbest ;
   float vbest ;

   if( icod == NULL || jcod == NULL || kcod == NULL ) return ; /* bad */

   *icod = *jcod = *kcod = 0 ; /* error returns, if sh*t happens */

   /* load column vectors for each (i,j,k) direction from matrix */

   /*-- i axis --*/ /*-- j axis --*/ /*-- k axis --*/

   xi = R.m[0][0] ; xj = R.m[0][1] ; xk = R.m[0][2] ;
   yi = R.m[1][0] ; yj = R.m[1][1] ; yk = R.m[1][2] ;
   zi = R.m[2][0] ; zj = R.m[2][1] ; zk = R.m[2][2] ;

   /* normalize column vectors to get unit vectors along each ijk-axis */

   /* normalize i axis */

   val = sqrt( xi*xi + yi*yi + zi*zi ) ;
   if( val == 0.0 ) return ;                 /* stupid input */
   xi /= val ; yi /= val ; zi /= val ;

   /* normalize j axis */

   val = sqrt( xj*xj + yj*yj + zj*zj ) ;
   if( val == 0.0 ) return ;                 /* stupid input */
   xj /= val ; yj /= val ; zj /= val ;

   /* orthogonalize j axis to i axis, if needed */

   val = xi*xj + yi*yj + zi*zj ;    /* dot product between i and j */
   if( fabs(val) > 1.e-4 ){
     xj -= val*xi ; yj -= val*yi ; zj -= val*zi ;
     val = sqrt( xj*xj + yj*yj + zj*zj ) ;  /* must renormalize */
     if( val == 0.0 ) return ;              /* j was parallel to i? */
     xj /= val ; yj /= val ; zj /= val ;
   }

   /* normalize k axis; if it is zero, make it the cross product i x j */

   val = sqrt( xk*xk + yk*yk + zk*zk ) ;
   if( val == 0.0 ){ xk = yi*zj-zi*yj; yk = zi*xj-zj*xi ; zk=xi*yj-yi*xj ; }
   else            { xk /= val ; yk /= val ; zk /= val ; }

   /* orthogonalize k to i */

   val = xi*xk + yi*yk + zi*zk ;    /* dot product between i and k */
   if( fabs(val) > 1.e-4 ){
     xk -= val*xi ; yk -= val*yi ; zk -= val*zi ;
     val = sqrt( xk*xk + yk*yk + zk*zk ) ;
     if( val == 0.0 ) return ;      /* bad */
     xk /= val ; yk /= val ; zk /= val ;
   }

   /* orthogonalize k to j */

   val = xj*xk + yj*yk + zj*zk ;    /* dot product between j and k */
   if( fabs(val) > 1.e-4 ){
     xk -= val*xj ; yk -= val*yj ; zk -= val*zj ;
     val = sqrt( xk*xk + yk*yk + zk*zk ) ;
     if( val == 0.0 ) return ;      /* bad */
     xk /= val ; yk /= val ; zk /= val ;
   }

   Q.m[0][0] = xi ; Q.m[0][1] = xj ; Q.m[0][2] = xk ;
   Q.m[1][0] = yi ; Q.m[1][1] = yj ; Q.m[1][2] = yk ;
   Q.m[2][0] = zi ; Q.m[2][1] = zj ; Q.m[2][2] = zk ;

   /* at this point, Q is the rotation matrix from the (i,j,k) to (x,y,z) axes */

   detQ = nifti_mat33_determ( Q ) ;
   if( detQ == 0.0 ) return ; /* shouldn't happen unless user is a DUFIS */

   /* Build and test all possible +1/-1 coordinate permutation matrices P;
      then find the P such that the rotation matrix M=PQ is closest to the
      identity, in the sense of M having the smallest total rotation angle. */

   /* Despite the formidable looking 6 nested loops, there are
      only 3*3*3*2*2*2 = 216 passes, which will run very quickly. */

   vbest = -666.0 ; ibest=pbest=qbest=rbest=1 ; jbest=2 ; kbest=3 ;
   for( i=1 ; i <= 3 ; i++ ){     /* i = column number to use for row #1 */
    for( j=1 ; j <= 3 ; j++ ){    /* j = column number to use for row #2 */
     if( i == j ) continue ;
      for( k=1 ; k <= 3 ; k++ ){  /* k = column number to use for row #3 */
       if( i == k || j == k ) continue ;
       P.m[0][0] = P.m[0][1] = P.m[0][2] =
        P.m[1][0] = P.m[1][1] = P.m[1][2] =
         P.m[2][0] = P.m[2][1] = P.m[2][2] = 0.0 ;
       for( p=-1 ; p <= 1 ; p+=2 ){    /* p,q,r are -1 or +1      */
        for( q=-1 ; q <= 1 ; q+=2 ){   /* and go into rows #1,2,3 */
         for( r=-1 ; r <= 1 ; r+=2 ){
           P.m[0][i-1] = p ; P.m[1][j-1] = q ; P.m[2][k-1] = r ;
           detP = nifti_mat33_determ(P) ;           /* sign of permutation */
           if( detP * detQ <= 0.0 ) continue ;  /* doesn't match sign of Q */
           M = nifti_mat33_mul(P,Q) ;

           /* angle of M rotation = 2.0*acos(0.5*sqrt(1.0+trace(M)))       */
           /* we want largest trace(M) == smallest angle == M nearest to I */

           val = M.m[0][0] + M.m[1][1] + M.m[2][2] ; /* trace */
           if( val > vbest ){
             vbest = val ;
             ibest = i ; jbest = j ; kbest = k ;
             pbest = p ; qbest = q ; rbest = r ;
           }
   }}}}}}

   /* At this point ibest is 1 or 2 or 3; pbest is -1 or +1; etc.

      The matrix P that corresponds is the best permutation approximation
      to Q-inverse; that is, P (approximately) takes (x,y,z) coordinates
      to the (i,j,k) axes.

      For example, the first row of P (which contains pbest in column ibest)
      determines the way the i axis points relative to the anatomical
      (x,y,z) axes.  If ibest is 2, then the i axis is along the y axis,
      which is direction P2A (if pbest > 0) or A2P (if pbest < 0).

      So, using ibest and pbest, we can assign the output code for
      the i axis.  Mutatis mutandis for the j and k axes, of course. */

   switch( ibest*pbest ){
     case  1: i = NIFTI_L2R ; break ;
     case -1: i = NIFTI_R2L ; break ;
     case  2: i = NIFTI_P2A ; break ;
     case -2: i = NIFTI_A2P ; break ;
     case  3: i = NIFTI_I2S ; break ;
     case -3: i = NIFTI_S2I ; break ;
   }

   switch( jbest*qbest ){
     case  1: j = NIFTI_L2R ; break ;
     case -1: j = NIFTI_R2L ; break ;
     case  2: j = NIFTI_P2A ; break ;
     case -2: j = NIFTI_A2P ; break ;
     case  3: j = NIFTI_I2S ; break ;
     case -3: j = NIFTI_S2I ; break ;
   }

   switch( kbest*rbest ){
     case  1: k = NIFTI_L2R ; break ;
     case -1: k = NIFTI_R2L ; break ;
     case  2: k = NIFTI_P2A ; break ;
     case -2: k = NIFTI_A2P ; break ;
     case  3: k = NIFTI_I2S ; break ;
     case -3: k = NIFTI_S2I ; break ;
   }

   *icod = i ; *jcod = j ; *kcod = k ; return ;
}

/*---------------------------------------------------------------------------*/
/* Routines to swap byte arrays in various ways:
    -  2 at a time:  ab               -> ba               [short]
    -  4 at a time:  abcd             -> dcba             [int, float]
    -  8 at a time:  abcdDCBA         -> ABCDdcba         [long long, double]
    - 16 at a time:  abcdefghHGFEDCBA -> ABCDEFGHhgfedcba [long double]
-----------------------------------------------------------------------------*/

typedef struct { unsigned char a,b ; } twobytes ;

/*----------------------------------------------------------------------*/
/*! swap each byte pair from the given list of n pairs
*//*--------------------------------------------------------------------*/
void nifti_swap_2bytes( int n , void *ar )    /* 2 bytes at a time */
{
   register int ii ;
   register twobytes *tb = (twobytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].b ; tb[ii].b = tt ;
   }
   return ;
}

/*---------------------------------------------------------------------------*/

typedef struct { unsigned char a,b,c,d ; } fourbytes ;

/*----------------------------------------------------------------------*/
/*! swap 4 bytes at a time from the given list of n sets of 4 bytes
*//*--------------------------------------------------------------------*/
void nifti_swap_4bytes( int n , void *ar )    /* 4 bytes at a time */
{
   register int ii ;
   register fourbytes *tb = (fourbytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].d ; tb[ii].d = tt ;
     tt = tb[ii].b ; tb[ii].b = tb[ii].c ; tb[ii].c = tt ;
   }
   return ;
}

/*---------------------------------------------------------------------------*/

typedef struct { unsigned char a,b,c,d , D,C,B,A ; } eightbytes ;

/*----------------------------------------------------------------------*/
/*! swap 8 bytes at a time from the given list of n sets of 8 bytes
*//*--------------------------------------------------------------------*/
void nifti_swap_8bytes( int n , void *ar )    /* 8 bytes at a time */
{
   register int ii ;
   register eightbytes *tb = (eightbytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].A ; tb[ii].A = tt ;
     tt = tb[ii].b ; tb[ii].b = tb[ii].B ; tb[ii].B = tt ;
     tt = tb[ii].c ; tb[ii].c = tb[ii].C ; tb[ii].C = tt ;
     tt = tb[ii].d ; tb[ii].d = tb[ii].D ; tb[ii].D = tt ;
   }
   return ;
}

/*---------------------------------------------------------------------------*/

typedef struct { unsigned char a,b,c,d,e,f,g,h ,
                               H,G,F,E,D,C,B,A  ; } sixteenbytes ;

/*----------------------------------------------------------------------*/
/*! swap 16 bytes at a time from the given list of n sets of 16 bytes
*//*--------------------------------------------------------------------*/
void nifti_swap_16bytes( int n , void *ar )    /* 16 bytes at a time */
{
   register int ii ;
   register sixteenbytes *tb = (sixteenbytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].A ; tb[ii].A = tt ;
     tt = tb[ii].b ; tb[ii].b = tb[ii].B ; tb[ii].B = tt ;
     tt = tb[ii].c ; tb[ii].c = tb[ii].C ; tb[ii].C = tt ;
     tt = tb[ii].d ; tb[ii].d = tb[ii].D ; tb[ii].D = tt ;

     tt = tb[ii].e ; tb[ii].e = tb[ii].E ; tb[ii].E = tt ;
     tt = tb[ii].f ; tb[ii].f = tb[ii].F ; tb[ii].F = tt ;
     tt = tb[ii].g ; tb[ii].g = tb[ii].G ; tb[ii].G = tt ;
     tt = tb[ii].h ; tb[ii].h = tb[ii].H ; tb[ii].H = tt ;
   }
   return ;
}

/*---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*! based on siz, call the appropriate nifti_swap_Nbytes() function
*//*--------------------------------------------------------------------*/
void nifti_swap_Nbytes( int n , int siz , void *ar )  /* subsuming case */
{
   switch( siz ){
     case 2:  nifti_swap_2bytes ( n , ar ) ; break ;
     case 4:  nifti_swap_4bytes ( n , ar ) ; break ;
     case 8:  nifti_swap_8bytes ( n , ar ) ; break ;
     case 16: nifti_swap_16bytes( n , ar ) ; break ;
   }
   return ;
}


/*-------------------------------------------------------------------------*/
/*! Byte swap NIFTI-1 file header in various places and ways.

    If is_nifti is nonzero, will also swap the NIFTI-specific
    components of the header; otherwise, only the components
    common to NIFTI and ANALYZE will be swapped.
*//*---------------------------------------------------------------------- */
void swap_nifti_header( struct nifti_1_header *h , int is_nifti )
{

#if 0                /* ANALYZE fields not used by this software */
   swap_4(h->sizeof_hdr) ;
   swap_4(h->extents) ;
   swap_2(h->session_error) ;
   swap_4(h->compressed) ;
   swap_4(h->glmax) ; swap_4(h->glmin) ;
#endif

   /* this stuff is always present, for ANALYZE and NIFTI */

   swap_4(h->sizeof_hdr) ;
   nifti_swap_2bytes( 8 , h->dim ) ;
   nifti_swap_4bytes( 8 , h->pixdim ) ;

   swap_2(h->datatype) ;
   swap_2(h->bitpix) ;

   swap_4(h->vox_offset); swap_4(h->cal_max); swap_4(h->cal_min);

   /* this stuff is NIFTI specific */

   if( is_nifti ){
     swap_4(h->intent_p1); swap_4(h->intent_p2); swap_4(h->intent_p3);
     swap_2(h->intent_code);

     swap_2(h->slice_start);    swap_2(h->slice_end);
     swap_4(h->scl_slope);      swap_4(h->scl_inter);
     swap_4(h->slice_duration); swap_4(h->toffset);

     swap_2(h->qform_code); swap_2(h->sform_code);
     swap_4(h->quatern_b); swap_4(h->quatern_c); swap_4(h->quatern_d);
     swap_4(h->qoffset_x); swap_4(h->qoffset_y); swap_4(h->qoffset_z);
     nifti_swap_4bytes(4,h->srow_x);
     nifti_swap_4bytes(4,h->srow_y);
     nifti_swap_4bytes(4,h->srow_z);
   }
   return ;
}


#define USE_STAT
#ifdef  USE_STAT
/*---------------------------------------------------------------------------*/
/* Return the file length (0 if file not found or has no contents).
   This is a Unix-specific function, since it uses stat().
-----------------------------------------------------------------------------*/
#include <sys/types.h>
#include <sys/stat.h>

/*---------------------------------------------------------------------------*/
/*! return the size of a file, in bytes

    \return size of file on success, -1 on error or no file

    changed to return int, -1 means no file or error      20 Dec 2004 [rickr]
*//*-------------------------------------------------------------------------*/
int nifti_get_filesize( const char *pathname )
{
   struct stat buf ; int ii ;

   if( pathname == NULL || *pathname == '\0' ) return -1 ;
   ii = stat( pathname , &buf ); if( ii != 0 ) return -1 ;
   return (unsigned int)buf.st_size ;
}

#else  /*---------- non-Unix version of the above, less efficient -----------*/

int nifti_get_filesize( const char *pathname )
{
   znzFile fp ; int len ;

   if( pathname == NULL || *pathname == '\0' ) return -1 ;
   fp = znzopen(pathname,"rb",0); if( znz_isnull(fp) ) return -1 ;
   znzseek(fp,0L,SEEK_END) ; len = znztell(fp) ;
   znzclose(fp) ; return len ;
}

#endif /* USE_STAT */


/*----------------------------------------------------------------------*/
/*! return the total volume size, in bytes

    This is computed as nvox * nbyper.
*//*--------------------------------------------------------------------*/
size_t nifti_get_volsize(const nifti_image *nim)
{
   return (size_t)(nim->nbyper) * (size_t)(nim->nvox) ; /* total bytes */
}


/*--------------------------------------------------------------------------*/
/* Support functions for filenames in read and write
   - allows for gzipped files
*/


/*----------------------------------------------------------------------*/
/*! simple check for file existence

    \return 1 on existence, 0 otherwise
*//*--------------------------------------------------------------------*/
int nifti_fileexists(const char* fname)
{
   znzFile fp;
   fp = znzopen( fname , "rb" , 1 ) ;
   if( !znz_isnull(fp) )  { znzclose(fp);  return 1; }
   return 0; /* fp is NULL */
}

/*----------------------------------------------------------------------*/
/*! return whether the filename is valid

    The name is considered valid if the file basename has length greater than
    zero, AND one of the valid nifti extensions is provided.
    fname input          | return |
    ===============================
    "myimage"            |  0     |
    "myimage.tif"        |  0     |
    "myimage.tif.gz"     |  0     |
    "myimage.nii"        |  1     |
    ".nii"               |  0     |
    ".myhiddenimage"     |  0     |
    ".myhiddenimage.nii" |  1     |
*//*--------------------------------------------------------------------*/
int nifti_is_complete_filename(const char* fname)
{
   char * ext;

   /* check input file(s) for sanity */
   if( fname == NULL || *fname == '\0' ){
      if ( g_opts.debug > 1 )
         fprintf(stderr,"-- empty filename in nifti_validfilename()\n");
      return 0;
   }

   ext = nifti_find_file_extension(fname);
   if ( ext == NULL ) { /*Invalid extension given */
      if ( g_opts.debug > 0 )
         fprintf(stderr,"-- no nifti valid extension for filename '%s'\n", fname);
       return 0;
   }

   if ( ext && ext == fname ) {   /* then no filename prefix */
      if ( g_opts.debug > 0 )
         fprintf(stderr,"-- no prefix for filename '%s'\n", fname);
      return 0;
   }
   return 1;
}

/*----------------------------------------------------------------------*/
/*! return whether the filename is valid

    The name is considered valid if its length is positive, excluding
    any nifti filename extension.
    fname input         |  return | result of nifti_makebasename
    ====================================================================
    "myimage"           |  1      | "myimage"
    "myimage.tif"       |  1      | "myimage.tif"
    "myimage.tif.gz"    |  1      | "myimage.tif"
    "myimage.nii"       |  1      | "myimage"
    ".nii"              |  0      | <ERROR - basename has zero length>
    ".myhiddenimage"    |  1      | ".myhiddenimage"
    ".myhiddenimage.nii |  1      | ".myhiddenimage"
*//*--------------------------------------------------------------------*/
int nifti_validfilename(const char* fname)
{
   char * ext;

   /* check input file(s) for sanity */
   if( fname == NULL || *fname == '\0' ){
      if ( g_opts.debug > 1 )
         fprintf(stderr,"-- empty filename in nifti_validfilename()\n");
      return 0;
   }

   ext = nifti_find_file_extension(fname);

   if ( ext && ext == fname ) {   /* then no filename prefix */
      if ( g_opts.debug > 0 )
         fprintf(stderr,"-- no prefix for filename '%s'\n", fname);
      return 0;
   }

   return 1;
}

/*----------------------------------------------------------------------*/
/*! check the end of the filename for a valid nifti extension

    Valid extensions are currently .nii, .hdr, .img, .nia,
    or any of them followed by .gz.  Note that '.' is part of
    the extension.

    \return a pointer to the extension (within the filename), or NULL
*//*--------------------------------------------------------------------*/
char * nifti_find_file_extension( const char * name )
{
   char * ext;
   int    len;

   if ( ! name ) return NULL;

   len = strlen(name);
   if ( len < 4 ) return NULL;

   ext = (char *)name + len - 4;

   if ( (strcmp(ext, ".hdr") == 0) || (strcmp(ext, ".img") == 0) ||
        (strcmp(ext, ".nia") == 0) || (strcmp(ext, ".nii") == 0) )
      return ext;

#ifdef HAVE_ZLIB
   if ( len < 7 ) return NULL;

   ext = (char *)name + len - 7;

   if ( (strcmp(ext, ".hdr.gz") == 0) || (strcmp(ext, ".img.gz") == 0) ||
        (strcmp(ext, ".nia.gz") == 0) || (strcmp(ext, ".nii.gz") == 0) )
      return ext;
#endif

   if( g_opts.debug > 1 )
      fprintf(stderr,"** find_file_ext: failed for name '%s'\n", name);

   return NULL;
}

/*----------------------------------------------------------------------*/
/*! return whether the filename ends in ".gz"
*//*--------------------------------------------------------------------*/
int nifti_is_gzfile(const char* fname)
{
  /* return true if the filename ends with .gz */
  if (fname == NULL) { return 0; }
#ifdef HAVE_ZLIB
  { /* just so len doesn't generate compile warning */
     int len;
     len = strlen(fname);
     if (len < 3) return 0;  /* so we don't search before the name */
     if (strcmp(fname + strlen(fname) - 3,".gz")==0) { return 1; }
  }
#endif
  return 0;
}


/*----------------------------------------------------------------------*/
/*! duplicate the filename, while clearing any extension

    This allocates memory for basename which should eventually be freed.
*//*--------------------------------------------------------------------*/
char * nifti_makebasename(const char* fname)
{
   char *basename, *ext;

   basename=nifti_strdup(fname);

   ext = nifti_find_file_extension(basename);
   if ( ext ) *ext = '\0';  /* clear out extension */
   
   return basename;  /* in either case */
}

/*----------------------------------------------------------------------*/
/*! set nifti's global debug level, for status reporting

    - 0    : quiet, nothing is printed to the terminal, but errors
    - 1    : normal execution (the default)
    - 2, 3 : more details
*//*--------------------------------------------------------------------*/
void nifti_set_debug_level( int level )
{
    g_opts.debug = level;
}

/*----------------------------------------------------------------------*/
/*! set nifti's global skip_blank_ext flag            5 Sep 2006 [rickr]

    explicitly set to 0 or 1
*//*--------------------------------------------------------------------*/
void nifti_set_skip_blank_ext( int skip )
{
    g_opts.skip_blank_ext = skip ? 1 : 0;
}

/*----------------------------------------------------------------------*/
/*! check current directory for existing header file

    \return filename of header on success and NULL if no appropriate file
            could be found 

    NB: it allocates memory for hdrname which should be freed
        when no longer required
*//*-------------------------------------------------------------------*/
char * nifti_findhdrname(const char* fname)
{
   char *basename, *hdrname, *ext;
   char  elist[2][5] = { ".hdr", ".nii" };
   int   efirst;

   /**- check input file(s) for sanity */
   if( !nifti_validfilename(fname) ) return NULL;

   basename = nifti_makebasename(fname);
   if( !basename ) return NULL;   /* only on string alloc failure */

   /**- return filename if it has a valid extension and exists
         (except if it is an .img file (and maybe .gz)) */
   ext = nifti_find_file_extension(fname);
   if ( ext && nifti_fileexists(fname) ) { 
     if ( strncmp(ext,".img",4) != 0 ){
        hdrname = nifti_strdup(fname); 
        free(basename);
        return hdrname; 
     }
   }

   /* So the requested name is a basename, contains .img, or does not exist. */
   /* In any case, use basename. */

   /**- if .img, look for .hdr, .hdr.gz, .nii, .nii.gz, in that order */
   /**- else,    look for .nii, .nii.gz, .hdr, .hdr.gz, in that order */

   /* if we get more extension choices, this could be a loop */

   if ( ext && strncmp(ext,".img",4) == 0 ) efirst = 0;
   else                                     efirst = 1;

   hdrname = (char *)calloc(sizeof(char),strlen(basename)+8);
   if( !hdrname ){
      fprintf(stderr,"** nifti_findhdrname: failed to alloc hdrname\n");
      free(basename);
      return NULL;
   }

   strcpy(hdrname,basename);
   strcat(hdrname,elist[efirst]);
   if (nifti_fileexists(hdrname)) { free(basename); return hdrname; }
#ifdef HAVE_ZLIB
   strcat(hdrname,".gz"); 
   if (nifti_fileexists(hdrname)) { free(basename); return hdrname; }
#endif

   /* okay, try the other possibility */

   efirst = 1 - efirst;

   strcpy(hdrname,basename);
   strcat(hdrname,elist[efirst]);
   if (nifti_fileexists(hdrname)) { free(basename); return hdrname; }
#ifdef HAVE_ZLIB
   strcat(hdrname,".gz"); 
   if (nifti_fileexists(hdrname)) { free(basename); return hdrname; }
#endif

   /**- if nothing has been found, return NULL */
   free(basename); 
   free(hdrname); 
   return NULL;
}


/*------------------------------------------------------------------------*/
/*! check current directory for existing image file

    \param fname filename to check for
    \nifti_type  nifti_type for dataset - this determines whether to
                 first check for ".nii" or ".img" (since both may exist)

    \return filename of data/img file on success and NULL if no appropriate
            file could be found

    NB: it allocates memory for the image filename, which should be freed
        when no longer required
*//*---------------------------------------------------------------------*/
char * nifti_findimgname(const char* fname , int nifti_type)
{
   char *basename, *imgname, ext[2][5] = { ".nii", ".img" };
   int  first;  /* first extension to use */

   /* check input file(s) for sanity */

   if( !nifti_validfilename(fname) ) return NULL;

   basename =  nifti_makebasename(fname);
   imgname = (char *)calloc(sizeof(char),strlen(basename)+8);
   if( !imgname ){
      fprintf(stderr,"** nifti_findimgname: failed to alloc imgname\n");
      free(basename);
      return NULL;
   }

   /**- test for .nii and .img (don't assume input type from image type) */
   /**- if nifti_type = 1, check for .nii first, else .img first         */

   /* if we get 3 or more extensions, can make a loop here... */

   if (nifti_type == NIFTI_FTYPE_NIFTI1_1) first = 0;   /* should match .nii */
   else                                    first = 1;   /* should match .img */

   strcpy(imgname,basename);
   strcat(imgname,ext[first]);
   if (nifti_fileexists(imgname)) { free(basename); return imgname; }
#ifdef HAVE_ZLIB  /* then also check for .gz */
   strcat(imgname,".gz");
   if (nifti_fileexists(imgname)) { free(basename); return imgname; }
#endif

   /* failed to find image file with expected extension, try the other */

   strcpy(imgname,basename);
   strcat(imgname,ext[1-first]);  /* can do this with only 2 choices */
   if (nifti_fileexists(imgname)) { free(basename); return imgname; }
#ifdef HAVE_ZLIB  /* then also check for .gz */
   strcat(imgname,".gz");
   if (nifti_fileexists(imgname)) { free(basename); return imgname; }
#endif

   /**- if nothing has been found, return NULL */
   free(basename); 
   free(imgname); 
   return NULL;
}


/*----------------------------------------------------------------------*/
/*! creates a filename for storing the header, based on nifti_type

   \param   prefix      - this will be copied before the suffix is added
   \param   nifti_type  - determines the extension, unless one is in prefix
   \param   check       - check for existence (fail condition)
   \param   comp        - add .gz for compressed name

   Note that if prefix provides a file suffix, nifti_type is not used.
  
   NB: this allocates memory which should be freed

   \sa nifti_set_filenames
*//*-------------------------------------------------------------------*/
char * nifti_makehdrname(const char * prefix, int nifti_type, int check,
                         int comp)
{
   char * iname, * ext;

   if( !nifti_validfilename(prefix) ) return NULL;

   /* add space for extension, optional ".gz", and null char */
   iname = (char *)calloc(sizeof(char),strlen(prefix)+8);
   if( !iname ){ fprintf(stderr,"** small malloc failure!\n"); return NULL; }
   strcpy(iname, prefix);

   /* use any valid extension */
   if( (ext = nifti_find_file_extension(iname)) != NULL ){
      if( strncmp(ext,".img",4) == 0 )
         memcpy(ext,".hdr",4);   /* then convert img name to hdr */
   }
   /* otherwise, make one up an */
   else if( nifti_type == NIFTI_FTYPE_NIFTI1_1 ) strcat(iname, ".nii");
   else if( nifti_type == NIFTI_FTYPE_ASCII )    strcat(iname, ".nia");
   else                                          strcat(iname, ".hdr");

#ifdef HAVE_ZLIB  /* then also check for .gz */
   if( comp && (!ext || !strstr(iname,".gz")) ) strcat(iname,".gz");
#endif

   /* check for existence failure */
   if( check && nifti_fileexists(iname) ){
      fprintf(stderr,"** failure: header file '%s' already exists\n",iname);
      free(iname);
      return NULL;
   }

   if(g_opts.debug > 2) fprintf(stderr,"+d made header filename '%s'\n", iname);

   return iname;
}
   

/*----------------------------------------------------------------------*/
/*! creates a filename for storing the image, based on nifti_type

   \param   prefix      - this will be copied before the suffix is added
   \param   nifti_type  - determines the extension, unless provided by prefix
   \param   check       - check for existence (fail condition)
   \param   comp        - add .gz for compressed name
  
   Note that if prefix provides a file suffix, nifti_type is not used.
  
   NB: it allocates memory which should be freed

   \sa nifti_set_filenames
*//*-------------------------------------------------------------------*/
char * nifti_makeimgname(const char * prefix, int nifti_type, int check,
                         int comp)
{
   char * iname, * ext;

   if( !nifti_validfilename(prefix) ) return NULL;

   /* add space for extension, optional ".gz", and null char */
   iname = (char *)calloc(sizeof(char),strlen(prefix)+8);
   if( !iname ){ fprintf(stderr,"** small malloc failure!\n"); return NULL; }
   strcpy(iname, prefix);

   /* use any valid extension */
   if( (ext = nifti_find_file_extension(iname)) != NULL ){
      if( strncmp(ext,".hdr",4) == 0 )
         memcpy(ext,".img",4);   /* then convert hdr name to img */
   }
   /* otherwise, make one up */
   else if( nifti_type == NIFTI_FTYPE_NIFTI1_1 ) strcat(iname, ".nii");
   else if( nifti_type == NIFTI_FTYPE_ASCII )    strcat(iname, ".nia");
   else                                          strcat(iname, ".img");

#ifdef HAVE_ZLIB  /* then also check for .gz */
   if( comp && (!ext || !strstr(iname,".gz")) ) strcat(iname,".gz");
#endif

   /* check for existence failure */
   if( check && nifti_fileexists(iname) ){
      fprintf(stderr,"** failure: image file '%s' already exists\n",iname);
      free(iname);
      return NULL;
   }

   if( g_opts.debug > 2 ) fprintf(stderr,"+d made image filename '%s'\n",iname);

   return iname;
}


/*----------------------------------------------------------------------*/
/*! create and set new filenames, based on prefix and image type

   \param nim            pointer to nifti_image in which to set filenames
   \param prefix         (required) prefix for output filenames
   \param check          check for previous existence of filename
                         (existence is an error condition)
   \param set_byte_order flag to set nim->byteorder here
                         (this is probably a logical place to do so)

   \return 0 on successful update
 
   \warning this will free() any existing names and create new ones

   \sa nifti_makeimgname, nifti_makehdrname, nifti_type_and_names_match
*//*--------------------------------------------------------------------*/
int nifti_set_filenames( nifti_image * nim, const char * prefix, int check,
                         int set_byte_order )
{
   int comp = nifti_is_gzfile(prefix);

   if( !nim || !prefix ){
      fprintf(stderr,"** nifti_set_filenames, bad params %p, %p\n",
              (void *)nim,prefix);
      return -1;
   }

   if( g_opts.debug > 1 )
      fprintf(stderr,"+d modifying output filenames using prefix %s\n", prefix);

   if( nim->fname ) free(nim->fname);
   if( nim->iname ) free(nim->iname);
   nim->fname = nifti_makehdrname(prefix, nim->nifti_type, check, comp);
   nim->iname = nifti_makeimgname(prefix, nim->nifti_type, check, comp);
   if( !nim->fname || !nim->iname ){
      LNI_FERR("nifti_set_filename","failed to set prefix for",prefix);
      return -1;
   }

   if( set_byte_order ) nim->byteorder = nifti_short_order() ;

   if( nifti_set_type_from_names(nim) < 0 )
      return -1;

   if( g_opts.debug > 2 )
      fprintf(stderr,"+d have new filenames %s and %s\n",nim->fname,nim->iname);

   return 0;
}


/*--------------------------------------------------------------------------*/
/*! check whether nifti_type matches fname and iname for the nifti_image

    - if type 0 or 2, expect .hdr/.img pair
    - if type 1, expect .nii (and names must match)

    \param nim       given nifti_image
    \param show_warn if set, print a warning message for any mis-match

    \return
        -   1 if the values seem to match
        -   0 if there is a mis-match
        -  -1 if there is not sufficient information to create file(s)

    \sa NIFTI_FTYPE_* codes in nifti1_io.h
    \sa nifti_set_type_from_names, is_valid_nifti_type
*//*------------------------------------------------------------------------*/
int nifti_type_and_names_match( nifti_image * nim, int show_warn )
{
   char func[] = "nifti_type_and_names_match";
   char * ext_h, * ext_i;  /* header and image filename extensions */
   int  errs = 0;          /* error counter */

   /* sanity checks */
   if( !nim ){
      if( show_warn ) fprintf(stderr,"** %s: missing nifti_image\n", func);
      return -1;
   }
   if( !nim->fname ){
      if( show_warn ) fprintf(stderr,"** %s: missing header filename\n", func);
      errs++;
   }
   if( !nim->iname ){
      if( show_warn ) fprintf(stderr,"** %s: missing image filename\n", func);
      errs++;
   }
   if( !is_valid_nifti_type(nim->nifti_type) ){
      if( show_warn )
         fprintf(stderr,"** %s: bad nifti_type %d\n", func, nim->nifti_type);
      errs++;
   }

   if( errs ) return -1;   /* then do not proceed */

   /* get pointers to extensions */
   ext_h = nifti_find_file_extension( nim->fname );
   ext_i = nifti_find_file_extension( nim->iname );

   /* check for filename extensions */
   if( !ext_h ){
      if( show_warn )
         fprintf(stderr,"-d missing NIFTI extension in header filename, %s\n",
                 nim->fname);
      errs++;
   }
   if( !ext_i ){
      if( show_warn )
         fprintf(stderr,"-d missing NIFTI extension in image filename, %s\n",
                 nim->iname);
      errs++;
   }

   if( errs ) return 0;   /* do not proceed, but this is just a mis-match */

   /* general tests */
   if( nim->nifti_type == NIFTI_FTYPE_NIFTI1_1 ){  /* .nii */
      if( strncmp(ext_h,".nii",4) ) {
         if( show_warn )
            fprintf(stderr,
            "-d NIFTI_FTYPE 1, but no .nii extension in header filename, %s\n",
            nim->fname);
         errs++;
      }
      if( strncmp(ext_i,".nii",4) ) {
         if( show_warn )
            fprintf(stderr,
            "-d NIFTI_FTYPE 1, but no .nii extension in image filename, %s\n",
            nim->iname);
         errs++;
      }
      if( strcmp(nim->fname, nim->iname) != 0 ){
         if( show_warn )
            fprintf(stderr,
            "-d NIFTI_FTYPE 1, but header and image filenames differ: %s, %s\n",
            nim->fname, nim->iname);
         errs++;
      }
   }
   else if( (nim->nifti_type == NIFTI_FTYPE_NIFTI1_2) || /* .hdr/.img */
            (nim->nifti_type == NIFTI_FTYPE_ANALYZE) )
   {
      if( strncmp(ext_h,".hdr",4) != 0 ){
         if( show_warn )
            fprintf(stderr,"-d no '.hdr' extension, but NIFTI type is %d, %s\n",
                    nim->nifti_type, nim->fname);
         errs++;
      }
      if( strncmp(ext_i,".img",4) != 0 ){
         if( show_warn )
            fprintf(stderr,"-d no '.img' extension, but NIFTI type is %d, %s\n",
                    nim->nifti_type, nim->iname);
         errs++;
      }
   }
   /* ignore any other nifti_type */

   return 1;
}


/*--------------------------------------------------------------------------*/
/*! check whether the given type is on the "approved" list

    The code is valid if it is non-negative, and does not exceed
    NIFTI_MAX_FTYPE.

    \return 1 if nifti_type is valid, 0 otherwise
    \sa NIFTI_FTYPE_* codes in nifti1_io.h
*//*------------------------------------------------------------------------*/
int is_valid_nifti_type( int nifti_type )
{
   if( nifti_type >= NIFTI_FTYPE_ANALYZE &&   /* smallest type, 0 */
       nifti_type <= NIFTI_MAX_FTYPE )
      return 1;
   return 0;
}


/*--------------------------------------------------------------------------*/
/*! check whether the given type is on the "approved" list

    The type is explicitly checked against the NIFTI_TYPE_* list
    in nifti1.h.

    \return 1 if dtype is valid, 0 otherwise
    \sa NIFTI_TYPE_* codes in nifti1.h
*//*------------------------------------------------------------------------*/
int nifti_is_valid_datatype( int dtype )
{
   if( dtype == NIFTI_TYPE_UINT8        ||
       dtype == NIFTI_TYPE_INT16        ||
       dtype == NIFTI_TYPE_INT32        ||
       dtype == NIFTI_TYPE_FLOAT32      ||
       dtype == NIFTI_TYPE_COMPLEX64    ||
       dtype == NIFTI_TYPE_FLOAT64      ||
       dtype == NIFTI_TYPE_RGB24        ||
       dtype == NIFTI_TYPE_INT8         ||
       dtype == NIFTI_TYPE_UINT16       ||
       dtype == NIFTI_TYPE_UINT32       ||
       dtype == NIFTI_TYPE_INT64        ||
       dtype == NIFTI_TYPE_UINT64       ||
       dtype == NIFTI_TYPE_FLOAT128     ||
       dtype == NIFTI_TYPE_COMPLEX128   ||
       dtype == NIFTI_TYPE_COMPLEX256 ) return 1;
   return 0;
}


/*--------------------------------------------------------------------------*/
/*! set the nifti_type field based on fname and iname

    Note that nifti_type is changed only when it does not match
    the filenames.

    \return 0 on success, -1 on error

    \sa is_valid_nifti_type, nifti_type_and_names_match
*//*------------------------------------------------------------------------*/
int nifti_set_type_from_names( nifti_image * nim )
{
   /* error checking first */
   if( !nim ){ fprintf(stderr,"** NSTFN: no nifti_image\n");  return -1; }

   if( !nim->fname || !nim->iname ){
      fprintf(stderr,"** NSTFN: missing filename(s) fname @ %p, iname @ %p\n",
              nim->fname, nim->iname);
      return -1;
   }

   if( ! nifti_validfilename      ( nim->fname ) ||
       ! nifti_validfilename      ( nim->iname ) ||
       ! nifti_find_file_extension( nim->fname ) ||
       ! nifti_find_file_extension( nim->iname )
     ) {
      fprintf(stderr,"** NSTFN: invalid filename(s) fname='%s', iname='%s'\n",
              nim->fname, nim->iname);
      return -1;
   }

   if( g_opts.debug > 2 )
      fprintf(stderr,"-d verify nifti_type from filenames: %d",nim->nifti_type);

   /* not too picky here, do what must be done, and then verify */
   if( strcmp(nim->fname, nim->iname) == 0 )          /* one file, type 1 */
      nim->nifti_type = NIFTI_FTYPE_NIFTI1_1;
   else if( nim->nifti_type == NIFTI_FTYPE_NIFTI1_1 ) /* cannot be type 1 */
      nim->nifti_type = NIFTI_FTYPE_NIFTI1_2;

   if( g_opts.debug > 2 ) fprintf(stderr," -> %d\n",nim->nifti_type);

   if( g_opts.debug > 1 )  /* warn user about anything strange */
      nifti_type_and_names_match(nim, 1);

   if( is_valid_nifti_type(nim->nifti_type) ) return 0;  /* success! */

   fprintf(stderr,"** NSTFN: bad nifti_type %d, for '%s' and '%s'\n",
           nim->nifti_type, nim->fname, nim->iname);

   return -1;
}


/*--------------------------------------------------------------------------*/
/*! Determine if this is a NIFTI-formatted file.

   <pre>
   \return  0 if file looks like ANALYZE 7.5 [checks sizeof_hdr field == 348]
            1 if file marked as NIFTI (header+data in 1 file)
            2 if file marked as NIFTI (header+data in 2 files)
           -1 if it can't tell, file doesn't exist, etc.
   </pre>
*//*------------------------------------------------------------------------*/
int is_nifti_file( const char *hname )
{
   struct nifti_1_header nhdr ;
   znzFile fp ;
   int ii ;
   char *tmpname;

   /* bad input name? */

   if( !nifti_validfilename(hname) ) return -1 ;

   /* open file */

   tmpname = nifti_findhdrname(hname);
   if( tmpname == NULL ){
      if( g_opts.debug > 0 )
         fprintf(stderr,"** no header file found for '%s'\n",hname);
      return -1;
   }
   fp = znzopen( tmpname , "rb" , nifti_is_gzfile(tmpname) ) ;
   free(tmpname);
   if (znz_isnull(fp))                      return -1 ;  /* bad open? */

   /* read header, close file */

   ii = znzread( &nhdr , 1 , sizeof(nhdr) , fp ) ;
   znzclose( fp ) ;
   if( ii < (int) sizeof(nhdr) )               return -1 ;  /* bad read? */

   /* check for NIFTI-ness */

   if( NIFTI_VERSION(nhdr) != 0 ){
     return ( NIFTI_ONEFILE(nhdr) ) ? 1 : 2 ;
   }

   /* check for ANALYZE-ness (sizeof_hdr field == 348) */

   ii = nhdr.sizeof_hdr ;
   if( ii == (int)sizeof(nhdr) ) return 0 ;  /* matches */

   /* try byte-swapping header */

   swap_4(ii) ;
   if( ii == (int)sizeof(nhdr) ) return 0 ;  /* matches */

   return -1 ;                          /* not good */
}

static int print_hex_vals( const char * data, int nbytes, FILE * fp )
{
   int c;

   if ( !data || nbytes < 1 || !fp ) return -1;

   fputs("0x", fp);
   for ( c = 0; c < nbytes; c++ )
      fprintf(fp, " %x", data[c]);

   return 0;
}

/*----------------------------------------------------------------------*/
/*! display the contents of the nifti_1_header (send to stdout)

   \param info if non-NULL, print this character string
   \param hp   pointer to nifti_1_header
*//*--------------------------------------------------------------------*/
int disp_nifti_1_header( const char * info, const nifti_1_header * hp )
{
   int c;

   fputs( "-------------------------------------------------------\n", stdout );
   if ( info )  fputs( info, stdout );
   if ( !hp  ){ fputs(" ** no nifti_1_header to display!\n",stdout); return 1; }

   fprintf(stdout," nifti_1_header :\n"
           "    sizeof_hdr     = %d\n"
           "    data_type[10]  = ", hp->sizeof_hdr);
   print_hex_vals(hp->data_type, 10, stdout);
   fprintf(stdout, "\n"
           "    db_name[18]    = ");
   print_hex_vals(hp->db_name, 18, stdout);
   fprintf(stdout, "\n"
           "    extents        = %d\n"
           "    session_error  = %d\n"
           "    regular        = 0x%x\n"
           "    dim_info       = 0x%x\n",
      hp->extents, hp->session_error, hp->regular, hp->dim_info );
   fprintf(stdout, "    dim[8]         =");
   for ( c = 0; c < 8; c++ ) fprintf(stdout," %d", hp->dim[c]);
   fprintf(stdout, "\n"
           "    intent_p1      = %f\n"
           "    intent_p2      = %f\n"
           "    intent_p3      = %f\n"
           "    intent_code    = %d\n"
           "    datatype       = %d\n"
           "    bitpix         = %d\n"
           "    slice_start    = %d\n"
           "    pixdim[8]      =",
           hp->intent_p1, hp->intent_p2, hp->intent_p3, hp->intent_code,
           hp->datatype, hp->bitpix, hp->slice_start);
   /* break pixdim over 2 lines */
   for ( c = 0; c < 4; c++ ) fprintf(stdout," %f", hp->pixdim[c]);
   fprintf(stdout, "\n                    ");
   for ( c = 4; c < 8; c++ ) fprintf(stdout," %f", hp->pixdim[c]);
   fprintf(stdout, "\n"
           "    vox_offset     = %f\n"
           "    scl_slope      = %f\n"
           "    scl_inter      = %f\n"
           "    slice_end      = %d\n"
           "    slice_code     = %d\n"
           "    xyzt_units     = 0x%x\n"
           "    cal_max        = %f\n"
           "    cal_min        = %f\n"
           "    slice_duration = %f\n"
           "    toffset        = %f\n"
           "    glmax          = %d\n"
           "    glmin          = %d\n",
           hp->vox_offset, hp->scl_slope, hp->scl_inter, hp->slice_end,
           hp->slice_code, hp->xyzt_units, hp->cal_max, hp->cal_min,
           hp->slice_duration, hp->toffset, hp->glmax, hp->glmin);
   fprintf(stdout,
           "    descrip        = '%.80s'\n"
           "    aux_file       = '%.24s'\n"
           "    qform_code     = %d\n"
           "    sform_code     = %d\n"
           "    quatern_b      = %f\n"
           "    quatern_c      = %f\n"
           "    quatern_d      = %f\n"
           "    qoffset_x      = %f\n"
           "    qoffset_y      = %f\n"
           "    qoffset_z      = %f\n"
           "    srow_x[4]      = %f, %f, %f, %f\n"
           "    srow_y[4]      = %f, %f, %f, %f\n"
           "    srow_z[4]      = %f, %f, %f, %f\n"
           "    intent_name    = '%-.16s'\n"
           "    magic          = '%-.4s'\n",
           hp->descrip, hp->aux_file, hp->qform_code, hp->sform_code,
           hp->quatern_b, hp->quatern_c, hp->quatern_d,
           hp->qoffset_x, hp->qoffset_y, hp->qoffset_z,
           hp->srow_x[0], hp->srow_x[1], hp->srow_x[2], hp->srow_x[3],
           hp->srow_y[0], hp->srow_y[1], hp->srow_y[2], hp->srow_y[3],
           hp->srow_z[0], hp->srow_z[1], hp->srow_z[2], hp->srow_z[3],
           hp->intent_name, hp->magic);
   fputs( "-------------------------------------------------------\n", stdout );
   fflush(stdout);

   return 0;
}


#undef  ERREX
#define ERREX(msg)                                           \
 do{ fprintf(stderr,"** ERROR: nifti_convert_nhdr2nim: %s\n", (msg) ) ;  \
     return NULL ; } while(0)

/*----------------------------------------------------------------------*/
/*! convert a nifti_1_header into a nift1_image
  
   \return an allocated nifti_image, or NULL on failure
*//*--------------------------------------------------------------------*/
nifti_image* nifti_convert_nhdr2nim(struct nifti_1_header nhdr,
                                    const char * fname)
{
   int   ii , doswap , ioff ;
   int   is_nifti , is_onefile ;
   nifti_image *nim;

   nim = (nifti_image *) calloc( 1 , sizeof(nifti_image) ) ;
   if( !nim ) ERREX("failed to allocate nifti image");

   /* be explicit with pointers */
   nim->fname = NULL;
   nim->iname = NULL;
   nim->data = NULL;

   /**- check if we must swap bytes */
  
   doswap = need_nhdr_swap(nhdr.dim[0], nhdr.sizeof_hdr); /* swap data flag */

   if( doswap < 0 ){
      if( doswap == -1 ) ERREX("bad dim[0]") ;
      ERREX("bad sizeof_hdr") ;  /* else */
   }

   /**- determine if this is a NIFTI-1 compliant header */

   is_nifti = NIFTI_VERSION(nhdr) ;
   /*
    * before swapping header, record the Analyze75 orient code
    */
   if(!is_nifti)
     {
     /**- in analyze75, the orient code is at the same address as
      *   qform_code, but it's just one byte
      *   the qform_code will be zero, at which point you can check
      *   analyze75_orient if you care to.
      */
     unsigned char c = *((char *)(&nhdr.qform_code));
     nim->analyze75_orient = (analyze_75_orient_code)c;
     }
   if( doswap ) {
      if ( g_opts.debug > 3 ) disp_nifti_1_header("-d ni1 pre-swap: ", &nhdr);
      swap_nifti_header( &nhdr , is_nifti ) ;
   }

   if ( g_opts.debug > 2 ) disp_nifti_1_header("-d nhdr2nim : ", &nhdr);

   if( nhdr.datatype == DT_BINARY ||
       nhdr.datatype == DT_UNKNOWN  )    ERREX("bad datatype") ;

   if( nhdr.dim[1] <= 0 )                ERREX("bad dim[1]") ;

   /* fix bad dim[] values in the defined dimension range */
   for( ii=2 ; ii <= nhdr.dim[0] ; ii++ )
     if( nhdr.dim[ii] <= 0 ) nhdr.dim[ii] = 1 ;

   /* fix any remaining bad dim[] values, so garbage does not propagate */
   /* (only values 0 or 1 seem rational, otherwise set to arbirary 1)   */
   for( ii=nhdr.dim[0]+1 ; ii <= 7 ; ii++ )
     if( nhdr.dim[ii] != 1 && nhdr.dim[ii] != 0) nhdr.dim[ii] = 1 ;

#if 0  /* rely on dim[0], do not attempt to modify it   16 Nov 2005 [rickr] */

   /**- get number of dimensions (ignoring dim[0] now) */
   for( ii=7 ; ii >= 2 ; ii-- )            /* loop backwards until we  */
     if( nhdr.dim[ii] > 1 ) break ;        /* find a dim bigger than 1 */
   ndim = ii ;
#endif

   /**- set bad grid spacings to 1.0 */

   for( ii=1 ; ii <= nhdr.dim[0] ; ii++ ){
     if( nhdr.pixdim[ii] == 0.0         ||
         !IS_GOOD_FLOAT(nhdr.pixdim[ii])  ) nhdr.pixdim[ii] = 1.0 ;
   }

  is_onefile = is_nifti && NIFTI_ONEFILE(nhdr) ;

  if( is_nifti ) nim->nifti_type = (is_onefile) ? NIFTI_FTYPE_NIFTI1_1
                                                : NIFTI_FTYPE_NIFTI1_2 ;
  else           nim->nifti_type = NIFTI_FTYPE_ANALYZE ;
  
  ii = nifti_short_order() ;
  if( doswap )   nim->byteorder = REVERSE_ORDER(ii) ;
  else           nim->byteorder = ii ;
  

  /**- set dimensions of data array */
  
  nim->ndim = nim->dim[0] = nhdr.dim[0];
  nim->nx   = nim->dim[1] = nhdr.dim[1];
  nim->ny   = nim->dim[2] = nhdr.dim[2];
  nim->nz   = nim->dim[3] = nhdr.dim[3];
  nim->nt   = nim->dim[4] = nhdr.dim[4];
  nim->nu   = nim->dim[5] = nhdr.dim[5];
  nim->nv   = nim->dim[6] = nhdr.dim[6];
  nim->nw   = nim->dim[7] = nhdr.dim[7];

  for( ii=1, nim->nvox=1; ii <= nhdr.dim[0]; ii++ )
     nim->nvox *= nhdr.dim[ii];
  
  /**- set the type of data in voxels and how many bytes per voxel */
  
  nim->datatype = nhdr.datatype ;
  
  nifti_datatype_sizes( nim->datatype , &(nim->nbyper) , &(nim->swapsize) ) ;
  if( nim->nbyper == 0 ){ free(nim); ERREX("bad datatype"); }
  
  /**- set the grid spacings */
  
  nim->dx = nim->pixdim[1] = nhdr.pixdim[1] ;
  nim->dy = nim->pixdim[2] = nhdr.pixdim[2] ;
  nim->dz = nim->pixdim[3] = nhdr.pixdim[3] ;
  nim->dt = nim->pixdim[4] = nhdr.pixdim[4] ;
  nim->du = nim->pixdim[5] = nhdr.pixdim[5] ;
  nim->dv = nim->pixdim[6] = nhdr.pixdim[6] ;
  nim->dw = nim->pixdim[7] = nhdr.pixdim[7] ;
  
  /**- compute qto_xyz transformation from pixel indexes (i,j,k) to (x,y,z) */
  
  if( !is_nifti || nhdr.qform_code <= 0 ){
    /**- if not nifti or qform_code <= 0, use grid spacing for qto_xyz */
    
    nim->qto_xyz.m[0][0] = nim->dx ;  /* grid spacings */
    nim->qto_xyz.m[1][1] = nim->dy ;  /* along diagonal */
    nim->qto_xyz.m[2][2] = nim->dz ;
    
    /* off diagonal is zero */
    
    nim->qto_xyz.m[0][1]=nim->qto_xyz.m[0][2]=nim->qto_xyz.m[0][3] = 0.0;
    nim->qto_xyz.m[1][0]=nim->qto_xyz.m[1][2]=nim->qto_xyz.m[1][3] = 0.0;
    nim->qto_xyz.m[2][0]=nim->qto_xyz.m[2][1]=nim->qto_xyz.m[2][3] = 0.0;
    
    /* last row is always [ 0 0 0 1 ] */
    
    nim->qto_xyz.m[3][0]=nim->qto_xyz.m[3][1]=nim->qto_xyz.m[3][2] = 0.0;
    nim->qto_xyz.m[3][3]= 1.0 ;
    
    nim->qform_code = NIFTI_XFORM_UNKNOWN ;
    
    if( g_opts.debug > 1 ) fprintf(stderr,"-d no qform provided\n");
  } else {
    /**- else NIFTI: use the quaternion-specified transformation */
    
    nim->quatern_b = FIXED_FLOAT( nhdr.quatern_b ) ;
    nim->quatern_c = FIXED_FLOAT( nhdr.quatern_c ) ;
    nim->quatern_d = FIXED_FLOAT( nhdr.quatern_d ) ;
    
    nim->qoffset_x = FIXED_FLOAT(nhdr.qoffset_x) ;
    nim->qoffset_y = FIXED_FLOAT(nhdr.qoffset_y) ;
    nim->qoffset_z = FIXED_FLOAT(nhdr.qoffset_z) ;
    
    nim->qfac = (nhdr.pixdim[0] < 0.0) ? -1.0 : 1.0 ;  /* left-handedness? */
    
    nim->qto_xyz = nifti_quatern_to_mat44(
                      nim->quatern_b, nim->quatern_c, nim->quatern_d,
                      nim->qoffset_x, nim->qoffset_y, nim->qoffset_z,
                      nim->dx       , nim->dy       , nim->dz       ,
                      nim->qfac                                      ) ;

    nim->qform_code = nhdr.qform_code ;

    if( g_opts.debug > 1 )
       nifti_disp_matrix_orient("-d qform orientations:\n", nim->qto_xyz);
  }
  
  /**- load inverse transformation (x,y,z) -> (i,j,k) */
  
  nim->qto_ijk = nifti_mat44_inverse( nim->qto_xyz ) ;
  
  /**- load sto_xyz affine transformation, if present */
  
  if( !is_nifti || nhdr.sform_code <= 0 ){
    /**- if not nifti or sform_code <= 0, then no sto transformation */
    
    nim->sform_code = NIFTI_XFORM_UNKNOWN ;

    if( g_opts.debug > 1 ) fprintf(stderr,"-d no sform provided\n");
    
  } else {
    /**- else set the sto transformation from srow_*[] */
    
    nim->sto_xyz.m[0][0] = nhdr.srow_x[0] ;
    nim->sto_xyz.m[0][1] = nhdr.srow_x[1] ;
    nim->sto_xyz.m[0][2] = nhdr.srow_x[2] ;
    nim->sto_xyz.m[0][3] = nhdr.srow_x[3] ;
    
    nim->sto_xyz.m[1][0] = nhdr.srow_y[0] ;
    nim->sto_xyz.m[1][1] = nhdr.srow_y[1] ;
    nim->sto_xyz.m[1][2] = nhdr.srow_y[2] ;
    nim->sto_xyz.m[1][3] = nhdr.srow_y[3] ;
    
    nim->sto_xyz.m[2][0] = nhdr.srow_z[0] ;
    nim->sto_xyz.m[2][1] = nhdr.srow_z[1] ;
    nim->sto_xyz.m[2][2] = nhdr.srow_z[2] ;
    nim->sto_xyz.m[2][3] = nhdr.srow_z[3] ;
    
    /* last row is always [ 0 0 0 1 ] */
    
    nim->sto_xyz.m[3][0]=nim->sto_xyz.m[3][1]=nim->sto_xyz.m[3][2] = 0.0;
    nim->sto_xyz.m[3][3]= 1.0 ;
    
    nim->sto_ijk = nifti_mat44_inverse( nim->sto_xyz ) ;
    
    nim->sform_code = nhdr.sform_code ;

    if( g_opts.debug > 1 )
       nifti_disp_matrix_orient("-d sform orientations:\n", nim->sto_xyz);
  }
  
  /**- set miscellaneous NIFTI stuff */
  
  if( is_nifti ){
    nim->scl_slope   = FIXED_FLOAT( nhdr.scl_slope ) ;
    nim->scl_inter   = FIXED_FLOAT( nhdr.scl_inter ) ;
    
    nim->intent_code = nhdr.intent_code ;
    
    nim->intent_p1 = FIXED_FLOAT( nhdr.intent_p1 ) ;
    nim->intent_p2 = FIXED_FLOAT( nhdr.intent_p2 ) ;
    nim->intent_p3 = FIXED_FLOAT( nhdr.intent_p3 ) ;
    
    nim->toffset   = FIXED_FLOAT( nhdr.toffset ) ;
    
    memcpy(nim->intent_name,nhdr.intent_name,15); nim->intent_name[15] = '\0';
    
    nim->xyz_units  = XYZT_TO_SPACE(nhdr.xyzt_units) ;
    nim->time_units = XYZT_TO_TIME (nhdr.xyzt_units) ;
    
    nim->freq_dim  = DIM_INFO_TO_FREQ_DIM ( nhdr.dim_info ) ;
    nim->phase_dim = DIM_INFO_TO_PHASE_DIM( nhdr.dim_info ) ;
    nim->slice_dim = DIM_INFO_TO_SLICE_DIM( nhdr.dim_info ) ;
    
    nim->slice_code     = nhdr.slice_code  ;
    nim->slice_start    = nhdr.slice_start ;
    nim->slice_end      = nhdr.slice_end   ;
    nim->slice_duration = FIXED_FLOAT(nhdr.slice_duration) ;
  }
  
  /**- set Miscellaneous ANALYZE stuff */
  
  nim->cal_min = FIXED_FLOAT(nhdr.cal_min) ;
  nim->cal_max = FIXED_FLOAT(nhdr.cal_max) ;
  
  memcpy(nim->descrip ,nhdr.descrip ,79) ; nim->descrip [79] = '\0' ;
  memcpy(nim->aux_file,nhdr.aux_file,23) ; nim->aux_file[23] = '\0' ;
  
   /**- set ioff from vox_offset (but at least sizeof(header)) */

   is_onefile = is_nifti && NIFTI_ONEFILE(nhdr) ;

   if( is_onefile ){
     ioff = (int)nhdr.vox_offset ;
     if( ioff < (int) sizeof(nhdr) ) ioff = (int) sizeof(nhdr) ;
   } else {
     ioff = (int)nhdr.vox_offset ;
   }
   nim->iname_offset = ioff ;


   /**- deal with file names if set */
   if (fname!=NULL) {
       nifti_set_filenames(nim,fname,0,0);
       if (nim->iname==NULL)  { ERREX("bad filename"); }
   } else { 
     nim->fname = NULL;  
     nim->iname = NULL; 
   }

   /* clear extension fields */
   nim->num_ext = 0;
   nim->ext_list = NULL;

   return nim;
}

#undef  ERREX
#define ERREX(msg)                                           \
 do{ fprintf(stderr,"** ERROR: nifti_image_open(%s): %s\n",  \
             (hname != NULL) ? hname : "(null)" , (msg) ) ;  \
     return fptr ; } while(0)

/***************************************************************
 * nifti_image_open
 ***************************************************************/
/*! znzFile nifti_image_open( char *hname, char *opts , nifti_image **nim)
    \brief Read in NIFTI-1 or ANALYZE-7.5 file (pair) header information into a nifti_image struct.

    - The image data is not read from disk (it may be read later using
        nifti_image_load(), for example).
    - The image data will be stored in whatever data format the
        input data is; no scaling will be applied.
    - DT_BINARY data is not supported.
    - nifti_image_free() can be used to delete the returned struct,
        when you are done with it.

    \param hname filename of dataset .hdr or .nii file
    \param opts  options string for opening the header file
    \param nim   pointer to pointer to nifti_image struct
                 (this routine allocates the nifti_image struct)
    \return file pointer (gzippable) to the file with the image data,
                 ready for reading.
        <br>NULL if something fails badly.
    \sa nifti_image_load, nifti_image_free
 */
znzFile nifti_image_open(const char * hname, char * opts, nifti_image ** nim)
{
  znzFile fptr=NULL;
  /* open the hdr and reading it in, but do not load the data  */
  *nim = nifti_image_read(hname,0);
  /* open the image file, ready for reading (compressed works for all reads) */
  if( ((*nim) == NULL)      || ((*nim)->iname == NULL) ||
      ((*nim)->nbyper <= 0) || ((*nim)->nvox <= 0)       )
     ERREX("bad header info") ;
                                                                                
  /* open image data file */
  fptr = znzopen( (*nim)->iname, opts, nifti_is_gzfile((*nim)->iname) );
  if( znz_isnull(fptr) ) ERREX("Can't open data file") ;
                                                                                
  return fptr;
}


/*----------------------------------------------------------------------*/
/*! return an allocated and filled nifti_1_header struct

    Read the binary header from disk, and swap bytes if necessary.

    \return an allocated nifti_1_header struct, or NULL on failure

    \param hname   name of file containing header
    \param swapped if not NULL, return whether header bytes were swapped
    \param check   flag to check for invalid nifti_1_header

    \warning ASCII header type is not supported

    \sa nifti_image_read, nifti_image_free, nifti_image_read_bricks
*//*--------------------------------------------------------------------*/
nifti_1_header * nifti_read_header(const char * hname, int * swapped, int check)
{
   nifti_1_header   nhdr, * hptr;
   znzFile          fp;
   int              bytes, lswap;
   char           * hfile;
   char             fname[] = { "nifti_read_header" };
   
   /* determine file name to use for header */
   hfile = nifti_findhdrname(hname);
   if( hfile == NULL ){
      if( g_opts.debug > 0 )
         LNI_FERR(fname,"failed to find header file for", hname);
      return NULL;
   } else if( g_opts.debug > 1 )
      fprintf(stderr,"-d %s: found header filename '%s'\n",fname,hfile);

   fp = znzopen( hfile, "rb", nifti_is_gzfile(hfile) );
   if( znz_isnull(fp) ){
      if( g_opts.debug > 0 ) LNI_FERR(fname,"failed to open header file",hfile);
      free(hfile);
      return NULL;
   }

   free(hfile);  /* done with filename */

   if( has_ascii_header(fp) == 1 ){
      znzclose( fp );
      if( g_opts.debug > 0 )
         LNI_FERR(fname,"ASCII header type not supported",hname);
      return NULL;
   }

   /* read the binary header */ 
   bytes = znzread( &nhdr, 1, sizeof(nhdr), fp );
   znzclose( fp );                      /* we are done with the file now */

   if( bytes < (int)sizeof(nhdr) ){
      if( g_opts.debug > 0 ){
         LNI_FERR(fname,"bad binary header read for file", hname);
         fprintf(stderr,"  - read %d of %d bytes\n",bytes, (int)sizeof(nhdr));
      }
      return NULL;
   }

   /* now just decide on byte swapping */
   lswap = need_nhdr_swap(nhdr.dim[0], nhdr.sizeof_hdr); /* swap data flag */
   if( check && lswap < 0 ){
      LNI_FERR(fname,"bad nifti_1_header for file", hname);
      return NULL;
   }

   if( lswap ) {
      if ( g_opts.debug > 3 ) disp_nifti_1_header("-d nhdr pre-swap: ", &nhdr);
      swap_nifti_header( &nhdr , NIFTI_VERSION(nhdr) ) ;
   }

   if ( g_opts.debug > 2 ) disp_nifti_1_header("-d nhdr post-swap: ", &nhdr);

   if ( check && ! nifti_hdr_looks_good(&nhdr) ){
      LNI_FERR(fname,"nifti_1_header looks bad for file", hname);
      return NULL;
   }

   /* all looks good, so allocate memory for and return the header */
   hptr = (nifti_1_header *)malloc(sizeof(nifti_1_header));
   if( ! hptr ){
      fprintf(stderr,"** nifti_read_hdr: failed to alloc nifti_1_header\n");
      return NULL;
   }

   if( swapped ) *swapped = lswap;  /* only if they care <sniff!> */

   memcpy(hptr, &nhdr, sizeof(nifti_1_header));

   return hptr;
}


/*----------------------------------------------------------------------*/
/*! decide if this nifti_1_header structure looks reasonable

   Check dim[0], dim[1], sizeof_hdr, and datatype.
   Check magic string for "n+1".
   Maybe more tests will follow.

   \return 1 if the header seems valid, 0 otherwise

   \sa nifti_nim_is_valid, valid_nifti_extensions
*//*--------------------------------------------------------------------*/
int nifti_hdr_looks_good(const nifti_1_header * hdr)
{
   int nbyper, swapsize;
   int c, errs = 0;

   /* check dim[0] and sizeof_hdr */
   if( need_nhdr_swap(hdr->dim[0], hdr->sizeof_hdr) < 0 ){
      if( g_opts.debug > 0 )
         fprintf(stderr,"** bad nhdr fields: dim0, sizeof_hdr = %d, %d\n",
                 hdr->dim[0], hdr->sizeof_hdr);
      errs++;
   }

   /* check the valid dimension sizes (maybe dim[0] is bad) */
   for( c = 1; c <= hdr->dim[0] && c <= 7; c++ )
      if( hdr->dim[c] <= 0 ){
         if( g_opts.debug > 0 )
            fprintf(stderr,"** bad nhdr field: dim[%d] = %d\n",c,hdr->dim[c]);
         errs++;
      }

   /* check the magic string */
   if( (hdr->magic[0] != 'n')                           ||
       (hdr->magic[1] != 'i' && hdr->magic[1] != '+')   ||
       (hdr->magic[2] != '1')                           ||
       (hdr->magic[3] != '\0') )
   {
      if( g_opts.debug > 0 )
         fprintf(stderr,
            "** bad nhdr field: magic = '%.4s', should be \"n+1\" or \"ni1\"\n"
            "   (in hex) magic = 0x%02x%02x%02x%02x\n"
            "        should be = 0x6e2b3100  or  0x6e693100\n",
            hdr->magic, hdr->magic[0], hdr->magic[1],
            hdr->magic[2], hdr->magic[3]);
      errs++;
   }

   /* check the datatype */
   if( hdr->datatype == DT_BINARY || hdr->datatype == DT_UNKNOWN ){
      if( g_opts.debug > 0 )
         fprintf(stderr,"** bad nhdr field: datatype = %d\n",hdr->datatype);
      errs++;
   }

   nifti_datatype_sizes(hdr->datatype, &nbyper, &swapsize);
   if( nbyper == 0 ){
      if( g_opts.debug > 0 )
         fprintf(stderr,"** bad nhdr field: datatype = %d\n",hdr->datatype);
      errs++;
   }

   if( errs ) return 0;  /* problems */

   if( g_opts.debug > 2 ) fprintf(stderr,"-d nifti header looks good\n");

   return 1;   /* looks good */
}


/*----------------------------------------------------------------------
 * check whether byte swapping is needed
 * 
 * dim[0] should be in [0,7], and sizeof_hdr should be accurate
 *
 * \returns  > 0 : needs swap
 *             0 : does not need swap
 *           < 0 : error condition
 *----------------------------------------------------------------------*/
static int need_nhdr_swap( short dim0, int hdrsize )
{
   short d0    = dim0;     /* so we won't have to swap them on the stack */
   int   hsize = hdrsize;

   if( d0 != 0 ){     /* then use it for the check */
      if( d0 > 0 && d0 <= 7 ) return 0;

      nifti_swap_2bytes(1, &d0);        /* swap? */
      if( d0 > 0 && d0 <= 7 ) return 1;

      if( g_opts.debug > 1 ){
         fprintf(stderr,"** bad swapped d0 = %d, unswapped = ", d0);
         nifti_swap_2bytes(1, &d0);        /* swap? */
         fprintf(stderr,"%d\n", d0);
      }

      return -1;        /* bad, naughty d0 */
   } 

   /* dim[0] == 0 should not happen, but could, so try hdrsize */
   if( hsize == sizeof(nifti_1_header) ) return 0;

   nifti_swap_4bytes(1, &hsize);     /* swap? */
   if( hsize == sizeof(nifti_1_header) ) return 1;

   if( g_opts.debug > 1 ){
      fprintf(stderr,"** bad swapped hsize = %d, unswapped = ", hsize);
      nifti_swap_4bytes(1, &hsize);        /* swap? */
      fprintf(stderr,"%d\n", hsize);
   }

   return -2;     /* bad, naughty hsize */
}


/* use macro LNI_FILE_ERROR instead of ERREX()
#undef  ERREX
#define ERREX(msg)                                           \
 do{ fprintf(stderr,"** ERROR: nifti_image_read(%s): %s\n",  \
             (hname != NULL) ? hname : "(null)" , (msg) ) ;  \
     return NULL ; } while(0)
*/


/***************************************************************
 * nifti_image_read
 ***************************************************************/
/*! \brief Read a nifti header and optionally the data, creating a nifti_image.

        - The data buffer will be byteswapped if necessary.
        - The data buffer will not be scaled.
        - The data buffer is allocated with calloc().

    \param hname filename of the nifti dataset
    \param read_data Flag, true=read data blob, false=don't read blob.
    \return A pointer to the nifti_image data structure.

    \sa nifti_image_free, nifti_free_extensions, nifti_image_read_bricks
*/
nifti_image *nifti_image_read( const char *hname , int read_data )
{
   struct nifti_1_header  nhdr ;
   nifti_image           *nim ;
   znzFile                fp ;
   int                    rv, ii , filesize, remaining;
   char                   fname[] = { "nifti_image_read" };
   char                  *hfile=NULL;

   if( g_opts.debug > 1 ){
      fprintf(stderr,"-d image_read from '%s', read_data = %d",hname,read_data);
#ifdef HAVE_ZLIB
      fprintf(stderr,", HAVE_ZLIB = 1\n");
#else
      fprintf(stderr,", HAVE_ZLIB = 0\n");
#endif
   }

   /**- determine filename to use for header */
   hfile = nifti_findhdrname(hname);
   if( hfile == NULL ){
      if(g_opts.debug > 0)
         LNI_FERR(fname,"failed to find header file for", hname);
      return NULL;  /* check return */
   } else if( g_opts.debug > 1 )
      fprintf(stderr,"-d %s: found header filename '%s'\n",fname,hfile);

   if( nifti_is_gzfile(hfile) ) filesize = -1;  /* unknown */
   else                         filesize = nifti_get_filesize(hfile);

   fp = znzopen(hfile, "rb", nifti_is_gzfile(hfile));
   if( znz_isnull(fp) ){
      if( g_opts.debug > 0 ) LNI_FERR(fname,"failed to open header file",hfile);
      free(hfile);
      return NULL;
   }

   rv = has_ascii_header( fp );
   if( rv < 0 ){
      if( g_opts.debug > 0 ) LNI_FERR(fname,"short header read",hfile);
      znzclose( fp );
      free(hfile);
      return NULL;
   }
   else if ( rv == 1 )  /* process special file type */
      return nifti_read_ascii_image( fp, hfile, filesize, read_data );

   /* else, just process normally */

   /**- read binary header */

   ii = znzread( &nhdr , 1 , sizeof(nhdr) , fp ) ;       /* read the thing */

   /* keep file open so we can check for exts. after nifti_convert_nhdr2nim() */

   if( ii < (unsigned int) sizeof(nhdr) ){
      if( g_opts.debug > 0 ){
         LNI_FERR(fname,"bad binary header read for file", hfile);
         fprintf(stderr,"  - read %d of %d bytes\n",ii, (int)sizeof(nhdr));
      }
      znzclose(fp) ;
      free(hfile);
      return NULL;
   }

   /* create output image struct and set it up */

   /**- convert all nhdr fields to nifti_image fields */
   nim = nifti_convert_nhdr2nim(nhdr,hfile);

   if( nim == NULL ){
      znzclose( fp ) ;                                   /* close the file */
      if( g_opts.debug > 0 )
         LNI_FERR(fname,"cannot create nifti image from header",hfile);
      free(hfile); /* had to save this for debug message */
      return NULL;
   }

   if( g_opts.debug > 3 ){
      fprintf(stderr,"+d nifti_image_read(), have nifti image:\n");
      if( g_opts.debug > 2 ) nifti_image_infodump(nim);
   }

   /**- check for extensions (any errors here means no extensions) */
   if( NIFTI_ONEFILE(nhdr) ) remaining = nim->iname_offset - sizeof(nhdr);
   else                      remaining = filesize - sizeof(nhdr);

   (void)nifti_read_extensions(nim, fp, remaining);

   znzclose( fp ) ;                                      /* close the file */
   free(hfile);

   /**- read the data if desired, then bug out */
   if( read_data ){
      if( nifti_image_load( nim ) < 0 ){
         nifti_image_free(nim);          /* take ball, go home. */
         return NULL;
      }
   }
   else nim->data = NULL ;

   return nim ;
}


/*----------------------------------------------------------------------
 * has_ascii_header  - see if the NIFTI header is an ASCII format
 *
 * If the file starts with the ASCII string "<nifti_image", then
 * process the dataset as a type-3 .nia file.
 *
 * return:  -1 on error, 1 if true, or 0 if false
 *
 * NOTE: this is NOT part of the NIFTI-1 standard
 *----------------------------------------------------------------------*/
static int has_ascii_header( znzFile fp )
{
   char  buf[16];
   int   nread;

   if( znz_isnull(fp) ) return 0;

   nread = znzread( buf, 1, 12, fp );
   buf[12] = '\0';

   if( nread < 12 ) return -1;

   znzrewind(fp);  /* move back to the beginning, and check */

   if( strcmp(buf, "<nifti_image") == 0 ) return 1;

   return 0;
}


/*----------------------------------------------------------------------*/
/*! nifti_read_ascii_image  - process as a type-3 .nia image file

   return NULL on failure
  
   NOTE: this is NOT part of the NIFTI-1 standard
*//*--------------------------------------------------------------------*/
nifti_image * nifti_read_ascii_image(znzFile fp, char *fname, int flen,
                                     int read_data)
{
   nifti_image * nim;
   int           slen, txt_size, remain, rv = 0;
   char        * sbuf, lfunc[25] = { "nifti_read_ascii_image" };

   if( nifti_is_gzfile(fname) ){
      LNI_FERR(lfunc, "compressed file with negative offset", fname);
      free(fname);  znzclose(fp);  return NULL;
   }
   slen = flen;  /* slen will be our buffer length */

   if( g_opts.debug > 1 )
      fprintf(stderr,"-d %s: have ASCII NIFTI file of size %d\n",fname,slen);

   if( slen > 65530 ) slen = 65530 ;
   sbuf = (char *)calloc(sizeof(char),slen+1) ;
   if( !sbuf ){
      fprintf(stderr,"** %s: failed to alloc %d bytes for sbuf",lfunc,65530);
      free(fname);  znzclose(fp);  return NULL;
   }
   znzread( sbuf , 1 , slen , fp ) ;
   nim = nifti_image_from_ascii( sbuf, &txt_size ) ; free( sbuf ) ;
   if( nim == NULL ){
      LNI_FERR(lfunc,"failed nifti_image_from_ascii()",fname);
      free(fname);  znzclose(fp);  return NULL;
   }
   nim->nifti_type = NIFTI_FTYPE_ASCII ;

   /* compute remaining space for extensions */
   remain = flen - txt_size - (int)nifti_get_volsize(nim);
   if( remain > 4 ){
      /* read extensions (reposition file pointer, first) */
      znzseek(fp, txt_size, SEEK_SET);
      (void) nifti_read_extensions(nim, fp, remain);
   }

   free(fname);
   znzclose( fp ) ;

   nim->iname_offset = -1 ;  /* check from the end of the file */

   if( read_data ) rv = nifti_image_load( nim ) ;
   else            nim->data = NULL ;

   /* check for nifti_image_load() failure, maybe bail out */
   if( read_data && rv != 0 ){
      if( g_opts.debug > 1 )
         fprintf(stderr,"-d failed image_load, free nifti image struct\n");
      free(nim);
      return NULL;
   }

   return nim ;
}


/*----------------------------------------------------------------------
 * Read the extensions into the nifti_image struct   08 Dec 2004 [rickr]
 *
 * This function is called just after the header struct is read in, and
 * it is assumed the file pointer has not moved.  The value in remain
 * is assumed to be accurate, reflecting the bytes of space for potential
 * extensions.
 *
 * return the number of extensions read in, or < 0 on error
 *----------------------------------------------------------------------*/
static int nifti_read_extensions( nifti_image *nim, znzFile fp, int remain )
{
   nifti1_extender    extdr;      /* defines extension existence  */
   nifti1_extension   extn;       /* single extension to process  */
   nifti1_extension * Elist;      /* list of processed extensions */
   int                posn, count;

   if( !nim || znz_isnull(fp) ) {
      if( g_opts.debug > 0 )
         fprintf(stderr,"** nifti_read_extensions: bad inputs (%p,%p)\n",
                 (void *)nim, (void *)fp);
      return -1;
   }

   posn = znztell(fp);
 
   if( (posn != sizeof(nifti_1_header)) &&
       (nim->nifti_type != NIFTI_FTYPE_ASCII) )
      fprintf(stderr,"** WARNING: posn not header size (%d, %d)\n",
              posn, (int)sizeof(nifti_1_header));

   if( g_opts.debug > 2 )
      fprintf(stderr,"-d nre: posn = %d, offset = %d, type = %d, remain = %d\n",
              posn, nim->iname_offset, nim->nifti_type, remain);

   if( remain < 16 ){
      if( g_opts.debug > 2 ){
         if( g_opts.skip_blank_ext )
            fprintf(stderr,"-d no extender in '%s' is okay, as "
                           "skip_blank_ext is set\n",nim->fname);
         else
            fprintf(stderr,"-d remain=%d, no space for extensions\n",remain);
      }
      return 0;
   }

   count = znzread( extdr.extension, 1, 4, fp ); /* get extender */

   if( count < 4 ){
      if( g_opts.debug > 1 )
         fprintf(stderr,"-d file '%s' is too short for an extender\n",
                 nim->fname);
      return 0;
   }

   if( extdr.extension[0] != 1 ){
      if( g_opts.debug > 2 )
         fprintf(stderr,"-d extender[0] (%d) shows no extensions for '%s'\n",
                 extdr.extension[0], nim->fname);
      return 0;
   }

   remain -= 4;
   if( g_opts.debug > 2 )
      fprintf(stderr,"-d found valid 4-byte extender, remain = %d\n", remain);

   /* so we expect extensions, but have no idea of how many there may be */

   count = 0;
   Elist = NULL;
   while (nifti_read_next_extension(&extn, nim, remain, fp) > 0)
   {
      if( nifti_add_exten_to_list(&extn, &Elist, count+1) < 0 ){
         if( g_opts.debug > 0 )
            fprintf(stderr,"** failed adding ext %d to list\n", count);
         return -1;
      }

      /* we have a new extension */
      if( g_opts.debug > 1 ){
         fprintf(stderr,"+d found extension #%d, code = 0x%x, size = %d\n",
                 count, extn.ecode, extn.esize);
         if( extn.ecode == NIFTI_ECODE_AFNI && g_opts.debug > 2 ) /* ~XML */
            fprintf(stderr,"   AFNI extension: %.*s\n",
                    extn.esize-8,extn.edata);
         else if( extn.ecode == NIFTI_ECODE_COMMENT && g_opts.debug > 2 )
            fprintf(stderr,"   COMMENT extension: %.*s\n",        /* TEXT */
                    extn.esize-8,extn.edata);
      }
      remain -= extn.esize;
      count++;
   }

   if( g_opts.debug > 2 ) fprintf(stderr,"+d found %d extension(s)\n", count);

   nim->num_ext = count;
   nim->ext_list = Elist;

   return count;
}


/*----------------------------------------------------------------------*/
/*! nifti_add_extension - add an extension, with a copy of the data

   Add an extension to the nim->ext_list array.
   Fill this extension with a copy of the data, noting the
       length and extension code.

   \param nim    - nifti_image to add extension to
   \param data   - raw extension data
   \param length - length of raw extension data
   \param ecode  - extension code

   \sa extension codes NIFTI_ECODE_* in nifti1_io.h
   \sa nifti_free_extensions, valid_nifti_extensions, nifti_copy_extensions

   \return 0 on success, -1 on error (and free the entire list)
*//*--------------------------------------------------------------------*/
int nifti_add_extension(nifti_image *nim, const char * data, int len, int ecode)
{
   nifti1_extension ext;

   /* error are printed in functions */
   if( nifti_fill_extension(&ext, data, len, ecode) )                 return -1;
   if( nifti_add_exten_to_list(&ext, &nim->ext_list, nim->num_ext+1)) return -1;

   nim->num_ext++;  /* success, so increment */

   return 0;
}


/*----------------------------------------------------------------------*/
/* nifti_add_exten_to_list     - add a new nifti1_extension to the list

   We will append via "malloc, copy and free", because on an error,
   the list will revert to the previous one (sorry realloc(), only
   quality dolphins get to become part of St@rk!st brand tunafish).

   return 0 on success, -1 on error (and free the entire list)
*//*--------------------------------------------------------------------*/
static int nifti_add_exten_to_list( nifti1_extension *  new_ext,
                                    nifti1_extension ** list, int new_length )
{
   nifti1_extension * tmplist;
  
   tmplist = *list;
   *list = (nifti1_extension *)malloc(new_length * sizeof(nifti1_extension));

   /* check for failure first */
   if( ! *list ){
      fprintf(stderr,"** failed to alloc %d extension structs (%d bytes)\n",
              new_length, new_length*(int)sizeof(nifti1_extension));
      if( !tmplist ) return -1;  /* no old list to lose */

      *list = tmplist;  /* reset list to old one */
      return -1;
   }

   /* if an old list exists, copy the pointers and free the list */
   if( tmplist ){
      memcpy(*list, tmplist, (new_length-1)*sizeof(nifti1_extension));
      free(tmplist);
   }

   /* for some reason, I just don't like struct copy... */
   (*list)[new_length-1].esize = new_ext->esize;
   (*list)[new_length-1].ecode = new_ext->ecode;
   (*list)[new_length-1].edata = new_ext->edata;

   if( g_opts.debug > 2 )
      fprintf(stderr,"+d allocated and appended extension #%d to list\n",
              new_length);

   return 0;
}


/*----------------------------------------------------------------------*/
/* nifti_fill_extension  - given data and length, fill an extension struct

   Allocate memory for data, copy data, set the size and code.

   return 0 on success, -1 on error (and free the entire list)
*//*--------------------------------------------------------------------*/
static int nifti_fill_extension( nifti1_extension *ext, const char * data,
                                int len, int ecode)
{
   int esize;

   if( !ext || !data || len < 0 ){
      fprintf(stderr,"** fill_ext: bad params (%p,%p,%d)\n",
              (void *)ext, data, len);
      return -1;
   } else if( ! nifti_is_valid_ecode(ecode) ){
      fprintf(stderr,"** fill_ext: invalid ecode %d\n", ecode);
      return -1;
   }

   /* compute esize, first : len+8, and take ceiling up to a mult of 16 */
   esize = len+8;
   if( esize & 0xf ) esize = (esize + 0xf) & ~0xf;
   ext->esize = esize;

   /* allocate esize-8 (maybe more than len), using calloc for fill */
   ext->edata = (char *)calloc(esize-8, sizeof(char));
   if( !ext->edata ){
      fprintf(stderr,"** NFE: failed to alloc %d bytes for extension\n",len);
      return -1;
   }

   memcpy(ext->edata, data, len);  /* copy the data, using len */
   ext->ecode = ecode;             /* set the ecode */

   if( g_opts.debug > 2 )
      fprintf(stderr,"+d alloc %d bytes for ext len %d, ecode %d, esize %d\n",
              esize-8, len, ecode, esize);

   return 0;
}


/*----------------------------------------------------------------------
 * nifti_read_next_extension  - read a single extension from the file
 *
 * return (>= 0 is okay):
 *
 *     success      : esize
 *     no extension : 0
 *     error        : -1
 *----------------------------------------------------------------------*/
static int nifti_read_next_extension( nifti1_extension * nex, nifti_image *nim,
                                      int remain, znzFile fp )
{
   int swap = nim->byteorder != nifti_short_order();
   int count, size, code;

   /* first clear nex */
   nex->esize = nex->ecode = 0;
   nex->edata = NULL;

   if( remain < 16 ){
      if( g_opts.debug > 2 )
         fprintf(stderr,"-d only %d bytes remain, so no extension\n", remain);
      return 0;
   }

   /* must start with 4-byte size and code */
   count = znzread( &size, 4, 1, fp );
   if( count == 1 ) count += znzread( &code, 4, 1, fp );

   if( count != 2 ){
      if( g_opts.debug > 2 )
         fprintf(stderr,"-d current extension read failed\n");
      znzseek(fp, -4*count, SEEK_CUR); /* back up past any read */
      return 0;                        /* no extension, no error condition */
   }

   if( swap ){
      if( g_opts.debug > 2 )
         fprintf(stderr,"-d pre-swap exts: code %d, size %d\n", code, size);

      nifti_swap_4bytes(1, &size);
      nifti_swap_4bytes(1, &code);
   }

   if( g_opts.debug > 2 )
      fprintf(stderr,"-d potential extension: code %d, size %d\n", code, size);

   if( !nifti_check_extension(nim, size, code, remain) ){
      if( znzseek(fp, -8, SEEK_CUR) < 0 ){      /* back up past any read */
         fprintf(stderr,"** failure to back out of extension read!\n");
         return -1;
      }
      return 0;
   }

   /* now get the actual data */
   nex->esize = size;
   nex->ecode = code;

   size -= 8;  /* subtract space for size and code in extension */
   nex->edata = (char *)malloc(size * sizeof(char));
   if( !nex->edata ){
      fprintf(stderr,"** failed to allocate %d bytes for extension\n",size);
      return -1;
   }

   count = znzread(nex->edata, 1, size, fp);
   if( count < size ){
      if( g_opts.debug > 0 )
         fprintf(stderr,"-d read only %d (of %d) bytes for extension\n",
                 count, size);
      free(nex->edata);
      nex->edata = NULL;
      return -1;
   }

   /* success! */
   if( g_opts.debug > 2 )
      fprintf(stderr,"+d successfully read extension, code %d, size %d\n",
              nex->ecode, nex->esize);

   return nex->esize;
}


/*----------------------------------------------------------------------*/
/*! for each extension, check code, size and data pointer
*//*--------------------------------------------------------------------*/
int valid_nifti_extensions(const nifti_image * nim)
{
   nifti1_extension * ext;
   int                c, errs;

   if( nim->num_ext <= 0 || nim->ext_list == NULL ){
      if( g_opts.debug > 2 ) fprintf(stderr,"-d empty extension list\n");
      return 0;
   }

   /* for each extension, check code, size and data pointer */
   ext = nim->ext_list;
   errs = 0;
   for ( c = 0; c < nim->num_ext; c++ ){
      if( ! nifti_is_valid_ecode(ext->ecode) ) {
         if( g_opts.debug > 1 )
            fprintf(stderr,"-d ext %d, invalid code %d\n", c, ext->ecode);
         errs++;
      }

      if( ext->esize <= 0 ){
         if( g_opts.debug > 1 )
            fprintf(stderr,"-d ext %d, bad size = %d\n", c, ext->esize);
         errs++;
      } else if( ext->esize & 0xf ){
         if( g_opts.debug > 1 )
            fprintf(stderr,"-d ext %d, size %d not multiple of 16\n",
                    c, ext->esize);
         errs++;
      }

      if( ext->edata == NULL ){
         if( g_opts.debug > 1 ) fprintf(stderr,"-d ext %d, missing data\n", c);
         errs++;
      }

      ext++;
   }

   if( errs > 0 ){
      if( g_opts.debug > 0 )
         fprintf(stderr,"-d had %d extension errors, none will be written\n",
                 errs);
      return 0;
   }

   /* if we're here, we're good */
   return 1;
}


/*----------------------------------------------------------------------*/
/*! check whether the extension code is valid

    \return 1 if valid, 0 otherwise
*//*--------------------------------------------------------------------*/
int nifti_is_valid_ecode( int ecode )
{
   if( ecode < NIFTI_ECODE_IGNORE  ||   /* minimum code number (0) */
       ecode > NIFTI_MAX_ECODE     ||   /* maximum code number     */
       ecode & 1 )                      /* cannot be odd           */
      return 0;

   return 1;
}


/*----------------------------------------------------------------------
 * check for valid size and code, as well as can be done
 *----------------------------------------------------------------------*/
static int nifti_check_extension(nifti_image *nim, int size, int code, int rem)
{
   /* check for bad code before bad size */
   if( ! nifti_is_valid_ecode(code) ) {
      if( g_opts.debug > 2 )
         fprintf(stderr,"-d invalid extension code %d\n",code);
      return 0;
   }

   if( size < 16 ){
      if( g_opts.debug > 2 )
         fprintf(stderr,"-d ext size %d, no extension\n",size);
      return 0;
   }

   if( size > rem ){
      if( g_opts.debug > 2 )
         fprintf(stderr,"-d ext size %d, space %d, no extension\n", size, rem);
      return 0;
   }

   if( size & 0xf ){
      if( g_opts.debug > 2 )
         fprintf(stderr,"-d nifti extension size %d not multiple of 16\n",size);
      return 0;
   }

   if( nim->nifti_type == NIFTI_FTYPE_ASCII && size > LNI_MAX_NIA_EXT_LEN ){
      if( g_opts.debug > 2 )
         fprintf(stderr,"-d NVE, bad nifti_type 3 size %d\n", size);
      return 0;
   }

   return 1;
}


/*----------------------------------------------------------------------
 * nifti_image_load_prep  - prepare to read data
 *
 * Check nifti_image fields, open the file and seek to the appropriate
 * offset for reading.
 *
 * return NULL on failure
 *----------------------------------------------------------------------*/
static znzFile nifti_image_load_prep( nifti_image *nim )
{
   /* set up data space, open data file and seek, then call nifti_read_buffer */
   size_t ntot , ii , ioff;
   znzFile fp;
   char    fname[] = { "nifti_image_load_prep" };

   /**- perform sanity checks */
   if( nim == NULL      || nim->iname == NULL ||
       nim->nbyper <= 0 || nim->nvox <= 0       )
   {
      if ( g_opts.debug > 0 ){
         if( !nim ) fprintf(stderr,"** ERROR: N_image_load: no nifti image\n");
         else fprintf(stderr,"** ERROR: N_image_load: bad params (%p,%d,%u)\n",
                      nim->iname, nim->nbyper, (unsigned)nim->nvox);
      }
      return NULL;
   }

   ntot = nifti_get_volsize(nim) ; /* total bytes to read */

   /**- open image data file */

   fp = znzopen(nim->iname, "rb", nifti_is_gzfile(nim->iname));
   if( znz_isnull(fp) ){
      if( g_opts.debug > 0 ) LNI_FERR(fname,"cannot open data file",nim->iname);
      return NULL;
   }

   /**- get image offset: a negative offset means to figure from end of file */
   if( nim->iname_offset < 0 ){
     if( nifti_is_gzfile(nim->iname) ){
        if( g_opts.debug > 0 )
           LNI_FERR(fname,"negative offset for compressed file",nim->iname);
        znzclose(fp);
        return NULL;
     }
     ii = nifti_get_filesize( nim->iname ) ;
     if( ii <= 0 ){
        if( g_opts.debug > 0 ) LNI_FERR(fname,"empty data file",nim->iname);
        znzclose(fp);
        return NULL;
     }
     ioff = (ii > ntot) ? ii-ntot : 0 ;
   } else {                              /* non-negative offset   */
     ioff = nim->iname_offset ;          /* means use it directly */
   }

   /**- seek to the appropriate read position */
   if( znzseek(fp , ioff , SEEK_SET) < 0 ){
      fprintf(stderr,"** could not seek to offset %u in file '%s'\n",
              (unsigned)ioff, nim->iname);
      znzclose(fp);
      return NULL;
   }

   /**- and return the File pointer */
   return fp;
}


/*----------------------------------------------------------------------
 * nifti_image_load
 *----------------------------------------------------------------------*/
/*! \fn int nifti_image_load( nifti_image *nim )
    \brief Load the image blob into a previously initialized nifti_image.

        - If not yet set, the data buffer is allocated with calloc().
        - The data buffer will be byteswapped if necessary.
        - The data buffer will not be scaled.

    This function is used to read the image from disk.  It should be used
    after a function such as nifti_image_read(), so that the nifti_image
    structure is already initialized.

    \param  nim pointer to a nifti_image (previously initialized)
    \return 0 on success, -1 on failure
    \sa     nifti_image_read, nifti_image_free, nifti_image_unload
*/
int nifti_image_load( nifti_image *nim )
{
   /* set up data space, open data file and seek, then call nifti_read_buffer */
   size_t ntot , ii ;
   znzFile fp ;

   /**- open the file and position the FILE pointer */
   fp = nifti_image_load_prep( nim );

   if( fp == NULL ){
      if( g_opts.debug > 0 )
         fprintf(stderr,"** nifti_image_load, failed load_prep\n");
      return -1;
   }

   ntot = nifti_get_volsize(nim);

   /**- if the data pointer is not yet set, get memory space for the image */

   if( nim->data == NULL )
   {
     nim->data = (void *)calloc(1,ntot) ;  /* create image memory */
     if( nim->data == NULL ){
        if( g_opts.debug > 0 )
           fprintf(stderr,"** failed to alloc %d bytes for image data\n",
                   (int)ntot);
        znzclose(fp);
        return -1;
     }
   }

   /**- now that everything is set up, do the reading */
   ii = nifti_read_buffer(fp,nim->data,ntot,nim);
   if( ii < ntot ){
      znzclose(fp) ;
      free(nim->data) ;
      nim->data = NULL ;
      return -1 ;  /* errors were printed in nifti_read_buffer() */
   }

   /**- close the file */
   znzclose( fp ) ;

   return 0 ;
}


/* 30 Nov 2004 [rickr]
#undef  ERREX
#define ERREX(msg)                                               \
 do{ fprintf(stderr,"** ERROR: nifti_read_buffer: %s\n",(msg)) ;  \
     return 0; } while(0)
*/

/*----------------------------------------------------------------------*/
/*! read ntot bytes of data from an open file and byte swaps if necessary

   note that nifti_image is required for information on datatype, bsize
   (for any needed byte swapping), etc.

   This function does not allocate memory, so dataptr must be valid.
*//*--------------------------------------------------------------------*/
size_t nifti_read_buffer(znzFile fp, void* dataptr, size_t ntot, 
                                nifti_image *nim)
{
  size_t ii;

  if( dataptr == NULL ){
     if( g_opts.debug > 0 )
        fprintf(stderr,"** ERROR: nifti_read_buffer: NULL dataptr\n");
     return -1;
  }

  ii = znzread( dataptr , 1 , ntot , fp ) ;             /* data input */
  
  /* if read was short, fail */
  if( ii < ntot ){ 
    if( g_opts.debug > 0 )
       fprintf(stderr,"++ WARNING: nifti_read_buffer(%s):\n"
               "   data bytes needed = %u\n"
               "   data bytes input  = %u\n"
               "   number missing    = %u (set to 0)\n",
               nim->iname , (unsigned int)ntot ,
               (unsigned int)ii , (unsigned int)(ntot-ii) ) ;
    /* memset( (char *)(dataptr)+ii , 0 , ntot-ii ) ;  now failure [rickr] */
    return -1 ;
  }

  if( g_opts.debug > 2 )
    fprintf(stderr,"+d nifti_read_buffer: read %u bytes\n", (unsigned)ii);
  
  /* byte swap array if needed */
  
  if( nim->swapsize > 1 && nim->byteorder != nifti_short_order() )
    nifti_swap_Nbytes( ntot / nim->swapsize , nim->swapsize , dataptr ) ;
  
#ifdef isfinite
{
  /* check input float arrays for goodness, and fix bad floats */
  int fix_count = 0 ;
  
  switch( nim->datatype ){
    
    case NIFTI_TYPE_FLOAT32:
    case NIFTI_TYPE_COMPLEX64:{
        register float *far = (float *)dataptr ; register size_t jj,nj ;
        nj = ntot / sizeof(float) ;
        for( jj=0 ; jj < nj ; jj++ )   /* count fixes 30 Nov 2004 [rickr] */
           if( !IS_GOOD_FLOAT(far[jj]) ){
              far[jj] = 0 ;
              fix_count++ ;
           }
      }
      break ;
    
    case NIFTI_TYPE_FLOAT64:
    case NIFTI_TYPE_COMPLEX128:{
        register double *far = (double *)dataptr ; register size_t jj,nj ;
        nj = ntot / sizeof(double) ;
        for( jj=0 ; jj < nj ; jj++ )   /* count fixes 30 Nov 2004 [rickr] */
           if( !IS_GOOD_FLOAT(far[jj]) ){
              far[jj] = 0 ;
              fix_count++ ;
           }
      }
      break ;
    
  }

  if( g_opts.debug > 1 )
     fprintf(stderr,"+d in image, %d bad floats were set to 0\n", fix_count);
}
#endif
  
  return ii;
}

/*--------------------------------------------------------------------------*/
/*! Unload the data in a nifti_image struct, but keep the metadata.
*//*------------------------------------------------------------------------*/
void nifti_image_unload( nifti_image *nim )
{
   if( nim != NULL && nim->data != NULL ){
     free(nim->data) ; nim->data = NULL ;
   }
   return ;
}

/*--------------------------------------------------------------------------*/
/*! free 'everything' about a nifti_image struct (including the passed struct)

    free (only fields which are not NULL):
      - fname and iname
      - data
      - any ext_list[i].edata
      - ext_list
      - nim
*//*------------------------------------------------------------------------*/
void nifti_image_free( nifti_image *nim )
{
   if( nim == NULL ) return ;
   if( nim->fname != NULL ) free(nim->fname) ;
   if( nim->iname != NULL ) free(nim->iname) ;
   if( nim->data  != NULL ) free(nim->data ) ;
   (void)nifti_free_extensions( nim ) ;
   free(nim) ; return ;
}


/*--------------------------------------------------------------------------*/
/*! free the nifti extensions

    - If any edata pointer is set in the extension list, free() it.
    - Free ext_list, if it is set.
    - Clear num_ext and ext_list from nim.

    \return 0 on success, -1 on error

    \sa nifti_add_extension, nifti_copy_extensions
*//*------------------------------------------------------------------------*/
int nifti_free_extensions( nifti_image *nim )
{
   int c ;
   if( nim == NULL ) return -1;
   if( nim->num_ext > 0 && nim->ext_list ){
      for( c = 0; c < nim->num_ext; c++ )
         if ( nim->ext_list[c].edata ) free(nim->ext_list[c].edata);
      free(nim->ext_list);
   }
   /* or if it is inconsistent, warn the user (if we are not in quiet mode) */
   else if ( (nim->num_ext > 0 || nim->ext_list != NULL) && (g_opts.debug > 0) )
      fprintf(stderr,"** warning: nifti extension num/ptr mismatch (%d,%p)\n",
              nim->num_ext, (void *)nim->ext_list);

   if( g_opts.debug > 2 )
      fprintf(stderr,"+d free'd %d extension(s)\n", nim->num_ext);

   nim->num_ext = 0;
   nim->ext_list = NULL;

   return 0;
}


/*--------------------------------------------------------------------------*/
/*! Print to stdout some info about a nifti_image struct.
*//*------------------------------------------------------------------------*/
void nifti_image_infodump( const nifti_image *nim )
{
   char *str = nifti_image_to_ascii( nim ) ;
   /* stdout -> stderr   2 Dec 2004 [rickr] */
   if( str != NULL ){ fputs(str,stderr) ; free(str) ; }
   return ;
}


/*--------------------------------------------------------------------------
 * nifti_write_buffer just check for a null znzFile and call znzwrite
 *--------------------------------------------------------------------------*/
/*! \fn size_t nifti_write_buffer(znzFile fp, void *buffer, size_t numbytes)
    \brief write numbytes of buffer to file, fp

    \param fp           File pointer (from znzopen) to gzippable nifti datafile
    \param buffer       data buffer to be written
    \param numbytes     number of bytes in buffer to write
    \return number of bytes successfully written
*/
size_t nifti_write_buffer(znzFile fp, const void *buffer, size_t numbytes)
{
   /* Write all the image data at once (no swapping here) */
   size_t ss;
   if (znz_isnull(fp)){
      fprintf(stderr,"** ERROR: nifti_write_buffer: null file pointer\n");
      return 0;
   }
   ss = znzwrite( (void*)buffer , 1 , numbytes , fp ) ;
   return ss;
}


/*----------------------------------------------------------------------*/
/*! write the nifti_image data to file (from nim->data or from NBL)

   If NBL is not NULL, write the data from that structure.  Otherwise,
   write it out from nim->data.  No swapping is done here.

   \param  fp  : File pointer
   \param  nim : nifti_image corresponding to the data
   \param  NBL : optional source of write data (if NULL use nim->data)
  
   \return 0 on success, -1 on failure

   Note: the nifti_image byte_order is set as that of the current CPU.
         This is because such a conversion was made to the data upon
         reading, while byte_order was not set (so the programs would
         know what format the data was on disk).  Effectively, since
         byte_order should match what is on disk, it should bet set to
         that of the current CPU whenever new filenames are assigned.
*//*--------------------------------------------------------------------*/
int nifti_write_all_data(znzFile fp, nifti_image * nim,
                         const nifti_brick_list * NBL)
{
   size_t ss;
   int    bnum;

   if( !NBL ){ /* just write one buffer and get out of here */
      if( nim->data == NULL ){
         fprintf(stderr,"** NWAD: no image data to write\n");
         return -1;
      }

      ss = nifti_write_buffer(fp,nim->data,nim->nbyper * nim->nvox);
      if (ss < nim->nbyper * nim->nvox){
         fprintf(stderr,
            "** ERROR: NWAD: wrote only %u of %u bytes to file\n",
            (unsigned)ss, (unsigned)(nim->nbyper * nim->nvox));
         return -1;
      }

      if( g_opts.debug > 1 )
         fprintf(stderr,"+d wrote single image of %u bytes\n", (unsigned)ss);
   } else {
      if( ! NBL->bricks || NBL->nbricks <= 0 || NBL->bsize <= 0 ){
         fprintf(stderr,"** NWAD: no brick data to write (%p,%d,%u)\n",
                 (void *)NBL->bricks, NBL->nbricks, (unsigned)NBL->bsize);
         return -1;
      }

      for( bnum = 0; bnum < NBL->nbricks; bnum++ ){
         ss = nifti_write_buffer(fp, NBL->bricks[bnum], NBL->bsize);
         if( ss < NBL->bsize ){
            fprintf(stderr,
              "** NWAD ERROR: wrote %u of %u bytes of brick %d of %d to file",
               (unsigned)ss, (unsigned)NBL->bsize, bnum+1, NBL->nbricks);
            return -1;
         }
      }
      if( g_opts.debug > 1 )
         fprintf(stderr,"+d wrote image of %d brick(s), each of %u bytes\n",
                 NBL->nbricks, (unsigned int)NBL->bsize);
   }

   /* mark as being in this CPU byte order */
   nim->byteorder = nifti_short_order() ;

   return 0;
}

/* return number of extensions written, or -1 on error */
static int nifti_write_extensions(znzFile fp, nifti_image *nim)
{
   nifti1_extension * list;
   char               extdr[4] = { 0, 0, 0, 0 };
   int                c, size, ok = 1;

   if( znz_isnull(fp) || !nim || nim->num_ext < 0 ){
      if( g_opts.debug > 0 )
         fprintf(stderr,"** nifti_write_extensions, bad params\n");
      return -1;
   }

   /* if no extensions and user requests it, skip extender */
   if( g_opts.skip_blank_ext && (nim->num_ext == 0 || ! nim->ext_list ) ){
      if( g_opts.debug > 1 )
         fprintf(stderr,"-d no exts and skip_blank_ext set, "
                        "so skipping 4-byte extender\n");
      return 0;
   }

   /* if invalid extension list, clear num_ext */
   if( ! valid_nifti_extensions(nim) ) nim->num_ext = 0;

   /* write out extender block */
   if( nim->num_ext > 0 ) extdr[0] = 1;
   if( nifti_write_buffer(fp, extdr, 4) != 4 ){
      fprintf(stderr,"** failed to write extender\n");
      return -1;
   }

   list = nim->ext_list;
   for ( c = 0; c < nim->num_ext; c++ ){
      size = nifti_write_buffer(fp, &list->esize, sizeof(int));
      ok = (size == (int)sizeof(int));
      if( ok ){
         size = nifti_write_buffer(fp, &list->ecode, sizeof(int));
         ok = (size == (int)sizeof(int));
      }
      if( ok ){
         size = nifti_write_buffer(fp, list->edata, list->esize - 8);
         ok = (size == list->esize - 8);
      }

      if( !ok ){
         fprintf(stderr,"** failed while writing extension #%d\n",c);
         return -1;
      }

      list++;
   }

   if( g_opts.debug > 1 )
      fprintf(stderr,"+d wrote out %d extension(s)\n", nim->num_ext);

   return nim->num_ext;
}


/*----------------------------------------------------------------------*/
/*! basic initialization of a nifti_image struct (to a 1x1x1 image)
*//*--------------------------------------------------------------------*/
nifti_image* nifti_simple_init_nim(void)
{
  nifti_image *nim;
  struct nifti_1_header nhdr;
  int nbyper, swapsize;
  
   memset(&nhdr,0,sizeof(nhdr)) ;  /* zero out header, to be safe */
  
   nhdr.sizeof_hdr = sizeof(nhdr) ;
   nhdr.regular    = 'r' ;           /* for some stupid reason */

   nhdr.dim[0] = 3 ;
   nhdr.dim[1] = 1 ; nhdr.dim[2] = 1 ; nhdr.dim[3] = 1 ;
   nhdr.dim[4] = 0 ;

   nhdr.pixdim[0] = 0.0 ;
   nhdr.pixdim[1] = 1.0 ; nhdr.pixdim[2] = 1.0 ;
   nhdr.pixdim[3] = 1.0 ;

   nhdr.datatype = DT_FLOAT32 ;
   nifti_datatype_sizes( nhdr.datatype , &nbyper, &swapsize );
   nhdr.bitpix   = 8 * nbyper ;

   strcpy(nhdr.magic, "n+1");  /* init to single file */

   nim = nifti_convert_nhdr2nim(nhdr,NULL);
   nim->fname = NULL;
   nim->iname = NULL;
   return nim;
}


/*----------------------------------------------------------------------*/
/*! basic initialization of a nifti_1_header struct (with given dimensions)

   Return an allocated nifti_1_header struct, based on the given
   dimensions and datatype.
 
   \param arg_dims  : optional dim[8] array (default {3,1,1,1,0,0,0,0})
   \param arg_dtype : optional datatype (default DT_FLOAT32)

   \return pointer to allocated nifti_1_header struct
*//*--------------------------------------------------------------------*/
nifti_1_header * nifti_make_new_header(const int arg_dims[], int arg_dtype)
{
   nifti_1_header * nhdr;
   const int        default_dims[8] = { 3, 1, 1, 1, 0, 0, 0, 0 };
   const int      * dim;  /* either passed or default dims  */
   int              dtype; /* either passed or default dtype */
   int              c, nbyper, swapsize;
  
   /* if arg_dims is passed, apply it */
   if( arg_dims ) dim = arg_dims;
   else           dim = default_dims;

   /* validate dim: if there is any problem, apply default_dims */
   if( dim[0] < 1 || dim[0] > 7 ) {
      fprintf(stderr,"** nifti_simple_hdr_with_dims: bad dim[0]=%d\n",dim[0]);
      dim = default_dims;
   } else {
      for( c = 1; c <= dim[0]; c++ )
         if( dim[c] < 1 )
         {
            fprintf(stderr,
                "** nifti_simple_hdr_with_dims: bad dim[%d]=%d\n",c,dim[c]);
            dim = default_dims;
            break;
         }
   }

   /* validate dtype, too */
   dtype = arg_dtype;
   if( ! nifti_is_valid_datatype(dtype) ) {
      fprintf(stderr,"** nifti_simple_hdr_with_dims: bad dtype %d\n",dtype);
      dtype = DT_FLOAT32;
   }

   /* now populate the header struct */

   if( g_opts.debug > 1 )
      fprintf(stderr,"+d nifti_make_new_header, dim[0] = %d, datatype = %d\n",
              dim[0], dtype);

   nhdr = (nifti_1_header *)calloc(1,sizeof(nifti_1_header));
   if( !nhdr ){
      fprintf(stderr,"** nifti_make_new_header: failed to alloc hdr\n");
      return NULL;
   }
  
   nhdr->sizeof_hdr = sizeof(nifti_1_header) ;
   nhdr->regular    = 'r' ;           /* for some stupid reason */

   /* init dim and pixdim */
   nhdr->dim[0] = dim[0] ;
   nhdr->pixdim[0] = 0.0;
   for( c = 1; c <= dim[0]; c++ ) {
      nhdr->dim[c] = dim[c];
      nhdr->pixdim[c] = 1.0;
   }

   nhdr->datatype = dtype ;
   nifti_datatype_sizes( nhdr->datatype , &nbyper, &swapsize );
   nhdr->bitpix   = 8 * nbyper ;

   strcpy(nhdr->magic, "n+1");  /* init to single file */

   return nhdr;
}


/*----------------------------------------------------------------------*/
/*! basic creation of a nifti_image struct

   Create a nifti_image from the given dimensions and data type.
   Optinally, allocate zero-filled data.
  
   \param dims      : optional dim[8]   (default {3,1,1,1,0,0,0,0})
   \param datatype  : optional datatype (default DT_FLOAT32)
   \param data_fill : if flag is set, allocate zero-filled data for image

   \return pointer to allocated nifti_image struct
*//*--------------------------------------------------------------------*/
nifti_image * nifti_make_new_nim(const int dims[], int datatype, int data_fill)
{
   nifti_image    * nim;
   nifti_1_header * nhdr;

   nhdr = nifti_make_new_header(dims, datatype);
   if( !nhdr ) return NULL;  /* error already printed */

   nim = nifti_convert_nhdr2nim(*nhdr,NULL);
   free(nhdr);               /* in any case, we are done with this */
   if( !nim ){
      fprintf(stderr,"** NMNN: nifti_convert_nhdr2nim failure\n");
      return NULL;
   }

   if( g_opts.debug > 1 )
      fprintf(stderr,"+d nifti_make_new_nim, data_fill = %d\n",data_fill);

   if( data_fill ) {
      nim->data = calloc(nim->nvox, nim->nbyper);

      /* if we cannot allocate data, take ball and go home */
      if( !nim->data ) {
         fprintf(stderr,"** NMNN: failed to alloc %u bytes for data\n",
                 (unsigned)nim->nvox);
         nifti_image_free(nim);
         nim = NULL;
      }
   }

   return nim;
}


/*----------------------------------------------------------------------*/
/*! convert a nifti_image structure to a nifti_1_header struct

    No allocation is done, this should be used via structure copy.
    As in:
    <pre>
    nifti_1_header my_header;
    my_header = nifti_convert_nim2nhdr(my_nim_pointer);
    </pre>
*//*--------------------------------------------------------------------*/
struct nifti_1_header nifti_convert_nim2nhdr(const nifti_image * nim)
{
   struct nifti_1_header nhdr;

   memset(&nhdr,0,sizeof(nhdr)) ;  /* zero out header, to be safe */


   /**- load the ANALYZE-7.5 generic parts of the header struct */

   nhdr.sizeof_hdr = sizeof(nhdr) ;
   nhdr.regular    = 'r' ;             /* for some stupid reason */

   nhdr.dim[0] = nim->ndim ;
   nhdr.dim[1] = nim->nx ; nhdr.dim[2] = nim->ny ; nhdr.dim[3] = nim->nz ;
   nhdr.dim[4] = nim->nt ; nhdr.dim[5] = nim->nu ; nhdr.dim[6] = nim->nv ;
   nhdr.dim[7] = nim->nw ;

   nhdr.pixdim[0] = 0.0 ;
   nhdr.pixdim[1] = nim->dx ; nhdr.pixdim[2] = nim->dy ;
   nhdr.pixdim[3] = nim->dz ; nhdr.pixdim[4] = nim->dt ;
   nhdr.pixdim[5] = nim->du ; nhdr.pixdim[6] = nim->dv ;
   nhdr.pixdim[7] = nim->dw ;

   nhdr.datatype = nim->datatype ;
   nhdr.bitpix   = 8 * nim->nbyper ;

   if( nim->cal_max > nim->cal_min ){
     nhdr.cal_max = nim->cal_max ;
     nhdr.cal_min = nim->cal_min ;
   }

   if( nim->scl_slope != 0.0 ){
     nhdr.scl_slope = nim->scl_slope ;
     nhdr.scl_inter = nim->scl_inter ;
   }

   if( nim->descrip[0] != '\0' ){
     memcpy(nhdr.descrip ,nim->descrip ,79) ; nhdr.descrip[79] = '\0' ;
   }
   if( nim->aux_file[0] != '\0' ){
     memcpy(nhdr.aux_file ,nim->aux_file ,23) ; nhdr.aux_file[23] = '\0' ;
   }

   /**- Load NIFTI specific stuff into the header */

   if( nim->nifti_type > NIFTI_FTYPE_ANALYZE ){ /* then not ANALYZE */

     if( nim->nifti_type == NIFTI_FTYPE_NIFTI1_1 ) strcpy(nhdr.magic,"n+1") ;
     else                                          strcpy(nhdr.magic,"ni1") ;

     nhdr.pixdim[1] = fabs(nhdr.pixdim[1]) ; nhdr.pixdim[2] = fabs(nhdr.pixdim[2]) ;
     nhdr.pixdim[3] = fabs(nhdr.pixdim[3]) ; nhdr.pixdim[4] = fabs(nhdr.pixdim[4]) ;
     nhdr.pixdim[5] = fabs(nhdr.pixdim[5]) ; nhdr.pixdim[6] = fabs(nhdr.pixdim[6]) ;
     nhdr.pixdim[7] = fabs(nhdr.pixdim[7]) ;

     nhdr.intent_code = nim->intent_code ;
     nhdr.intent_p1   = nim->intent_p1 ;
     nhdr.intent_p2   = nim->intent_p2 ;
     nhdr.intent_p3   = nim->intent_p3 ;
     if( nim->intent_name[0] != '\0' ){
       memcpy(nhdr.intent_name,nim->intent_name,15) ;
       nhdr.intent_name[15] = '\0' ;
     }

     nhdr.vox_offset  = (float) nim->iname_offset ;
     nhdr.xyzt_units  = SPACE_TIME_TO_XYZT( nim->xyz_units, nim->time_units ) ;
     nhdr.toffset     = nim->toffset ;

     if( nim->qform_code > 0 ){
       nhdr.qform_code = nim->qform_code ;
       nhdr.quatern_b  = nim->quatern_b ;
       nhdr.quatern_c  = nim->quatern_c ;
       nhdr.quatern_d  = nim->quatern_d ;
       nhdr.qoffset_x  = nim->qoffset_x ;
       nhdr.qoffset_y  = nim->qoffset_y ;
       nhdr.qoffset_z  = nim->qoffset_z ;
       nhdr.pixdim[0]  = (nim->qfac >= 0.0) ? 1.0 : -1.0 ;
     }

     if( nim->sform_code > 0 ){
       nhdr.sform_code = nim->sform_code ;
       nhdr.srow_x[0]  = nim->sto_xyz.m[0][0] ;
       nhdr.srow_x[1]  = nim->sto_xyz.m[0][1] ;
       nhdr.srow_x[2]  = nim->sto_xyz.m[0][2] ;
       nhdr.srow_x[3]  = nim->sto_xyz.m[0][3] ;
       nhdr.srow_y[0]  = nim->sto_xyz.m[1][0] ;
       nhdr.srow_y[1]  = nim->sto_xyz.m[1][1] ;
       nhdr.srow_y[2]  = nim->sto_xyz.m[1][2] ;
       nhdr.srow_y[3]  = nim->sto_xyz.m[1][3] ;
       nhdr.srow_z[0]  = nim->sto_xyz.m[2][0] ;
       nhdr.srow_z[1]  = nim->sto_xyz.m[2][1] ;
       nhdr.srow_z[2]  = nim->sto_xyz.m[2][2] ;
       nhdr.srow_z[3]  = nim->sto_xyz.m[2][3] ;
     }

     nhdr.dim_info = FPS_INTO_DIM_INFO( nim->freq_dim ,
                                        nim->phase_dim , nim->slice_dim ) ;
     nhdr.slice_code     = nim->slice_code ;
     nhdr.slice_start    = nim->slice_start ;
     nhdr.slice_end      = nim->slice_end ;
     nhdr.slice_duration = nim->slice_duration ;
   }

   return nhdr;
}


/*----------------------------------------------------------------------*/
/*! \fn int nifti_copy_extensions(nifti_image * nim_dest, nifti_image * nim_src)
    \brief copy the nifti1_extension list from src to dest

    Duplicate the list of nifti1_extensions.  The dest structure must
    be clear of extensions.
    \return 0 on success, -1 on failure

    \sa nifti_add_extension, nifti_free_extensions
*/
int nifti_copy_extensions(nifti_image * nim_dest, const nifti_image * nim_src)
{
   char   * data;
   size_t   bytes;
   int      c, size, old_size;

   if( nim_dest->num_ext > 0 || nim_dest->ext_list != NULL ){
      fprintf(stderr,"** will not copy extensions over existing ones\n");
      return -1;
   }

   if( g_opts.debug > 1 )
      fprintf(stderr,"+d duplicating %d extension(s)\n", nim_src->num_ext);

   if( nim_src->num_ext <= 0 ) return 0;

   bytes = nim_src->num_ext * sizeof(nifti1_extension);  /* I'm lazy */
   nim_dest->ext_list = (nifti1_extension *)malloc(bytes);
   if( !nim_dest->ext_list ){
      fprintf(stderr,"** failed to allocate %d nifti1_extension structs\n",
              nim_src->num_ext);
      return -1;
   }

   /* copy the extension data */
   nim_dest->num_ext = 0;
   for( c = 0; c < nim_src->num_ext; c++ ){
      size = old_size = nim_src->ext_list[c].esize;
      if( size & 0xf ) size = (size + 0xf) & ~0xf; /* make multiple of 16 */
      if( g_opts.debug > 2 )
         fprintf(stderr,"+d dup'ing ext #%d of size %d (from size %d)\n",
                 c, size, old_size);
      /* data length is size-8, as esize includes space for esize and ecode */
      data = (char *)calloc(size-8,sizeof(char));      /* maybe size > old */
      if( !data ){
         fprintf(stderr,"** failed to alloc %d bytes for extention\n", size);
         if( c == 0 ) { free(nim_dest->ext_list); nim_dest->ext_list = NULL; }
         /* otherwise, keep what we have (a.o.t. deleting them all) */
         return -1;
      }
      /* finally, fill the new structure */
      nim_dest->ext_list[c].esize = size;
      nim_dest->ext_list[c].ecode = nim_src->ext_list[c].ecode;
      nim_dest->ext_list[c].edata = data;
      memcpy(data, nim_src->ext_list[c].edata, old_size-8);

      nim_dest->num_ext++;
   }

   return 0;
}


/*----------------------------------------------------------------------*/
/*! compute the total size of all extensions

    \return the total of all esize fields

    Note that each esize includes 4 bytes for ecode, 4 bytes for esize,
    and the bytes used for the data.  Each esize also needs to be a
    multiple of 16, so it may be greater than the sum of its 3 parts.
*//*--------------------------------------------------------------------*/
int nifti_extension_size(nifti_image *nim)
{
   int c, size = 0;

   if( !nim || nim->num_ext <= 0 ) return 0;

   if( g_opts.debug > 2 ) fprintf(stderr,"-d ext sizes:");

   for ( c = 0; c < nim->num_ext; c++ ){
      size += nim->ext_list[c].esize;
      if( g_opts.debug > 2 ) fprintf(stderr,"  %d",nim->ext_list[c].esize);
   }

   if( g_opts.debug > 2 ) fprintf(stderr," (total = %d)\n",size);

   return size;
}


/*----------------------------------------------------------------------*/
/*! set the nifti_image iname_offset field, based on nifti_type

    - if writing to 2 files, set offset to 0
    - if writing to a single NIFTI-1 file, set the offset to
         352 + total extension size, then align to 16-byte boundary
    - if writing an ASCII header, set offset to -1
*//*--------------------------------------------------------------------*/
void nifti_set_iname_offset(nifti_image *nim)
{
   int offset;

   switch( nim->nifti_type ){

     default:  /* writing into 2 files */
       /* we only write files with 0 offset in the 2 file format */
       nim->iname_offset = 0 ;
     break ;

     /* NIFTI-1 single binary file - always update */
     case NIFTI_FTYPE_NIFTI1_1:
       offset = nifti_extension_size(nim)+sizeof(struct nifti_1_header)+4;
       /* be sure offset is aligned to a 16 byte boundary */
       if ( ( offset % 16 ) != 0 )  offset = ((offset + 0xf) & ~0xf);
       if( nim->iname_offset != offset ){
          if( g_opts.debug > 1 )
             fprintf(stderr,"+d changing offset from %d to %d\n",
                  nim->iname_offset, offset);
          nim->iname_offset = offset;
       }
     break ;

     /* non-standard case: NIFTI-1 ASCII header + binary data (single file) */
     case NIFTI_FTYPE_ASCII:
       nim->iname_offset = -1 ;             /* compute offset from filesize */
     break ;
   }
}


/*----------------------------------------------------------------------*/
/*! write the nifti_image dataset to disk, optionally including data

   This is just a front-end for nifti_image_write_hdr_img2.

   \param nim        nifti_image to write to disk
   \param write_data write options (see nifti_image_write_hdr_img2)
   \param opts       file open options ("wb" from nifti_image_write)

   \sa nifti_image_write, nifti_image_write_hdr_img2, nifti_image_free,
       nifti_set_filenames
*//*--------------------------------------------------------------------*/
znzFile nifti_image_write_hdr_img( nifti_image *nim , int write_data , 
                                          const char* opts )
{
  return nifti_image_write_hdr_img2(nim,write_data,opts,NULL,NULL);
}


#undef  ERREX
#define ERREX(msg)                                                \
 do{ fprintf(stderr,"** ERROR: nifti_image_write_hdr_img: %s\n",(msg)) ;  \
     return fp ; } while(0)


/* ----------------------------------------------------------------------*/
/*! This writes the header (and optionally the image data) to file
 *
 * If the image data file is left open it returns a valid znzFile handle.
 * It also uses imgfile as the open image file is not null, and modifies
 * it inside.
 * 
 * \param nim        nifti_image to write to disk
 * \param write_opts flags whether to write data and/or close file (see below)
 * \param opts       file-open options, probably "wb" from nifti_image_write()
 * \param imgfile    optional open znzFile struct, for writing image data
                     (may be NULL)
 * \param NBL        optional nifti_brick_list, containing the image data
                     (may be NULL)
 * 
 * Values for write_opts mode are based on two binary flags
 * ( 0/1 for no-write/write data, and 0/2 for close/leave-open files ) :
 *    -   0 = do not write data and close (do not open data file)
 *    -   1 = write data        and close
 *    -   2 = do not write data and leave data file open 
 *    -   3 = write data        and leave data file open 
 *
 * \sa nifti_image_write, nifti_image_write_hdr_img, nifti_image_free,
 *     nifti_set_filenames
*//*---------------------------------------------------------------------*/
znzFile nifti_image_write_hdr_img2(nifti_image *nim, int write_opts, 
               const char * opts, znzFile imgfile, const nifti_brick_list * NBL)
{
   struct nifti_1_header nhdr ;
   znzFile               fp=NULL;
   size_t                ss ;
   int                   write_data, leave_open;
   char                  func[] = { "nifti_image_write_hdr_img2" };

   write_data = write_opts & 1;  /* just separate the bits now */
   leave_open = write_opts & 2;

   if( ! nim                              ) ERREX("NULL input") ;
   if( ! nifti_validfilename(nim->fname)  ) ERREX("bad fname input") ;
   if( write_data && ! nim->data && ! NBL ) ERREX("no image data") ;

   nifti_set_iname_offset(nim);

   if( g_opts.debug > 1 ){
      fprintf(stderr,"-d writing nifti file '%s'...\n", nim->fname);
      if( g_opts.debug > 2 )
         fprintf(stderr,"-d nifti type %d, offset %d\n",
                 nim->nifti_type, nim->iname_offset);
   }

   if( nim->nifti_type == NIFTI_FTYPE_ASCII )   /* non-standard case */
      return nifti_write_ascii_image(nim,NBL,opts,write_data,leave_open);

   nhdr = nifti_convert_nim2nhdr(nim);    /* create the nifti1_header struct */

   /* if writing to 2 files, make sure iname is set and different from fname */
   if( nim->nifti_type != NIFTI_FTYPE_NIFTI1_1 ){
       if( nim->iname && strcmp(nim->iname,nim->fname) == 0 ){
         free(nim->iname) ; nim->iname = NULL ;
       }
       if( nim->iname == NULL ){ /* then make a new one */
         nim->iname = nifti_makeimgname(nim->fname,nim->nifti_type,0,0);
         if( nim->iname == NULL ) return NULL;  
       }
   }

   /* if we have an imgfile and will write the header there, use it */
   if( ! znz_isnull(imgfile) && nim->nifti_type == NIFTI_FTYPE_NIFTI1_1 ){
      if( g_opts.debug > 2 ) fprintf(stderr,"+d using passed file for hdr\n");
      fp = imgfile;
   }
   else {
      if( g_opts.debug > 2 )
         fprintf(stderr,"+d opening output file %s [%s]\n",nim->fname,opts);
      fp = znzopen( nim->fname , opts , nifti_is_gzfile(nim->fname) ) ;
      if( znz_isnull(fp) ){
         LNI_FERR(func,"cannot open output file",nim->fname);
         return fp;
      }
   }

   /* write the header and extensions */

   ss = znzwrite(&nhdr , 1 , sizeof(nhdr) , fp); /* write header */
   if( ss < sizeof(nhdr) ){
      LNI_FERR(func,"bad header write to output file",nim->fname);
      znzclose(fp); return fp;
   }

   /* partial file exists, and errors have been printed, so ignore return */
   if( nim->nifti_type != NIFTI_FTYPE_ANALYZE )
      (void)nifti_write_extensions(fp,nim);

   /* if the header is all we want, we are done */
   if( ! write_data && ! leave_open ){
      if( g_opts.debug > 2 ) fprintf(stderr,"-d header is all we want: done\n");
      znzclose(fp); return(fp);
   }

   if( nim->nifti_type != NIFTI_FTYPE_NIFTI1_1 ){ /* get a new file pointer */
      znzclose(fp);         /* first, close header file */
      if( ! znz_isnull(imgfile) ){
         if(g_opts.debug > 2) fprintf(stderr,"+d using passed file for img\n");
         fp = imgfile;
      }
      else {
         if( g_opts.debug > 2 )
            fprintf(stderr,"+d opening img file '%s'\n", nim->iname);
         fp = znzopen( nim->iname , opts , nifti_is_gzfile(nim->iname) ) ;
         if( znz_isnull(fp) ) ERREX("cannot open image file") ;
      }
   }

   znzseek(fp, nim->iname_offset, SEEK_SET);  /* in any case, seek to offset */

   if( write_data ) nifti_write_all_data(fp,nim,NBL);
   if( ! leave_open ) znzclose(fp);

   return fp;
}


/*----------------------------------------------------------------------*/
/*! write a nifti_image to disk in ASCII format
*//*--------------------------------------------------------------------*/
znzFile nifti_write_ascii_image(nifti_image *nim, const nifti_brick_list * NBL,
                              const char *opts, int write_data, int leave_open)
{
   znzFile   fp;
   char    * hstr;
                                                                                
   hstr = nifti_image_to_ascii( nim ) ;  /* get header in ASCII form */
   if( ! hstr ){ fprintf(stderr,"** failed image_to_ascii()\n"); return NULL; }
                                                                                
   fp = znzopen( nim->fname , opts , nifti_is_gzfile(nim->fname) ) ;
   if( znz_isnull(fp) ){
      free(hstr);
      fprintf(stderr,"** failed to open '%s' for ascii write\n",nim->fname);
      return fp;
   }
                                                                                
   znzputs(hstr,fp);                                               /* header */
   nifti_write_extensions(fp,nim);                             /* extensions */
                                                                                
   if ( write_data   ) { nifti_write_all_data(fp,nim,NBL); }         /* data */
   if ( ! leave_open ) { znzclose(fp); }
   free(hstr);
   return fp;  /* returned but may be closed */
}


/*--------------------------------------------------------------------------*/
/*! Write a nifti_image to disk.

   Since data is properly byte-swapped upon reading, it is assumed
   to be in the byte-order of the current CPU at write time.  Thus,
   nim->byte_order should match that of the current CPU.  Note that
   the nifti_set_filenames() function takes the flag, set_byte_order.

   The following fields of nim affect how the output appears:
    - nifti_type = 0 ==> ANALYZE-7.5 format file pair will be written
    - nifti_type = 1 ==> NIFTI-1 format single file will be written
                         (data offset will be 352+extensions)
    - nifti_type = 2 ==> NIFTI_1 format file pair will be written
    - nifti_type = 3 ==> NIFTI_1 ASCII single file will be written
    - fname is the name of the output file (header or header+data)
    - if a file pair is being written, iname is the name of the data file
    - existing files WILL be overwritten with extreme prejudice
    - if qform_code > 0, the quatern_*, qoffset_*, and qfac fields determine
      the qform output, NOT the qto_xyz matrix; if you want to compute these
      fields from the qto_xyz matrix, you can use the utility function
      nifti_mat44_to_quatern()

   \sa nifti_image_write_bricks, nifti_image_free, nifti_set_filenames,
       nifti_image_write_hdr_img
*//*------------------------------------------------------------------------*/
void nifti_image_write( nifti_image *nim )
{
   znzFile fp = nifti_image_write_hdr_img(nim,1,"wb");
   if( fp ){
      if( g_opts.debug > 2 ) fprintf(stderr,"-d niw: done with znzFile\n");
      free(fp);
   }
   if( g_opts.debug > 1 ) fprintf(stderr,"-d nifti_image_write: done\n");
}


/*----------------------------------------------------------------------*/
/*! similar to nifti_image_write, but data is in NBL struct, not nim->data

   \sa nifti_image_write, nifti_image_free, nifti_set_filenames, nifti_free_NBL
*//*--------------------------------------------------------------------*/
void nifti_image_write_bricks( nifti_image *nim, const nifti_brick_list * NBL )
{
   znzFile fp = nifti_image_write_hdr_img2(nim,1,"wb",NULL,NBL);
   if( fp ){
      if( g_opts.debug > 2 ) fprintf(stderr,"-d niwb: done with znzFile\n");
      free(fp);
   }
   if( g_opts.debug > 1 ) fprintf(stderr,"-d niwb: done writing bricks\n");
}


/*----------------------------------------------------------------------*/
/*! copy the nifti_image structure, without data

    Duplicate the structure, including fname, iname and extensions.
    Leave the data pointer as NULL.
*//*--------------------------------------------------------------------*/
nifti_image * nifti_copy_nim_info(const nifti_image * src)
{
  nifti_image *dest;
  dest = (nifti_image *)calloc(1,sizeof(nifti_image));
  if( !dest ){
     fprintf(stderr,"** NCNI: failed to alloc nifti_image\n");
     return NULL;
  }
  memcpy(dest, src, sizeof(nifti_image));
  if( src->fname ) dest->fname = nifti_strdup(src->fname);
  if( src->iname ) dest->iname = nifti_strdup(src->iname);
  dest->num_ext = 0;
  dest->ext_list = NULL;
  /* errors will be printed in NCE(), continue in either case */
  (void)nifti_copy_extensions(dest, src);

  dest->data = NULL;

  return dest;
}


/*------------------------------------------------------------------------*/
/* Un-escape a C string in place -- that is, convert XML escape sequences
   back into their characters.  (This can be done in place since the
   replacement is always smaller than the input.)  Escapes recognized are:
     -  &lt;   ->  <
     -  &gt;   ->  >
     -  &quot; ->  "
     -  &apos; ->  '
     -  &amp;  ->  &
   Also replace CR LF pair (Microsoft), or CR alone (Macintosh) with
   LF (Unix), per the XML standard.
   Return value is number of replacements made (if you care).
--------------------------------------------------------------------------*/

#undef  CR
#undef  LF
#define CR 0x0D
#define LF 0x0A

static int unescape_string( char *str )
{
   int ii,jj , nn,ll ;

   if( str == NULL ) return 0 ;                /* no string? */
   ll = strlen(str) ; if( ll == 0 ) return 0 ;

   /* scan for escapes: &something; */

   for( ii=jj=nn=0 ; ii<ll ; ii++,jj++ ){ /* scan at ii; results go in at jj */

     if( str[ii] == '&' ){  /* start of escape? */

             if( ii+3 < ll        &&   /* &lt; */
                 str[ii+1] == 'l' &&
                 str[ii+2] == 't' &&
                 str[ii+3] == ';'   ){ str[jj] = '<' ; ii += 3 ; nn++ ; }

        else if( ii+3 < ll        &&   /* &gt; */
                 str[ii+1] == 'g' &&
                 str[ii+2] == 't' &&
                 str[ii+3] == ';'   ){ str[jj] = '>' ; ii += 3 ; nn++ ; }

        else if( ii+5 < ll        &&   /* &quot; */
                 str[ii+1] == 'q' &&
                 str[ii+2] == 'u' &&
                 str[ii+3] == 'o' &&
                 str[ii+4] == 't' &&
                 str[ii+5] == ';'   ){ str[jj] = '"' ; ii += 5 ; nn++ ; }

        else if( ii+5 < ll        &&   /* &apos; */
                 str[ii+1] == 'a' &&
                 str[ii+2] == 'p' &&
                 str[ii+3] == 'o' &&
                 str[ii+4] == 's' &&
                 str[ii+5] == ';'   ){ str[jj] = '\'' ; ii += 5 ; nn++ ; }

        else if( ii+4 < ll        &&  /* &amp; */
                 str[ii+1] == 'a' &&
                 str[ii+2] == 'm' &&
                 str[ii+3] == 'p' &&
                 str[ii+4] == ';'   ){ str[jj] = '&' ; ii += 4 ; nn++ ; }

        /* although the comments above don't mention it,
           we also look for XML style numeric escapes
           of the forms &#32; (decimal) and &#xfd; (hex) */

        else if( ii+3 < ll        &&
                 str[ii+1] == '#' &&
                 isdigit(str[ii+2]) ){   /* &#dec; */

           unsigned int val='?' ; int kk=ii+3 ;
           while( kk < ll && kk != ';' ) kk++ ;
           sscanf( str+ii+2 , "%u" , &val ) ;
           str[jj] = (char) val ; ii = kk ; nn++ ;
        }

        else if( ii+4 < ll        &&
                 str[ii+1] == '#' &&
                 str[ii+2] == 'x' &&
                 isxdigit(str[ii+3]) ){   /* &#hex; */

           unsigned int val='?' ; int kk=ii+4 ;
           while( kk < ll && kk != ';' ) kk++ ;
           sscanf( str+ii+3 , "%x" , &val ) ;
           str[jj] = (char) val ; ii = kk ; nn++ ;
        }

        /* didn't start a recognized escape, so just copy as normal */

        else if( jj < ii ){ str[jj] = str[ii] ; }

     } else if( str[ii] == CR ) {  /* is a carriage return */

        if( str[ii+1] == LF ){ str[jj] = LF ; ii++ ; nn++ ; }  /* CR LF */
        else                 { str[jj] = LF ;      ; nn++ ; }  /* CR only */

     } else { /* is a normal character, just copy to output */

             if( jj < ii ){ str[jj] = str[ii] ; }
     }

     /* at this point, ii=index of last character used up in scan
                       jj=index of last character written to (jj <= ii) */
   }

   if( jj < ll ) str[jj] = '\0' ; /* end string properly */

   return nn ;
}

/*------------------------------------------------------------------------*/
/* Quotize (and escapize) one string, returning a new string.
   Approximately speaking, this is the inverse of unescape_string().
   The result should be free()-ed when you are done with it.
--------------------------------------------------------------------------*/

static char *escapize_string( const char * str )
{
   int ii,jj , lstr,lout ;
   char *out ;

   if( str == NULL || (lstr=strlen(str)) == 0 ){      /* 0 length */
     out = nifti_strdup("''") ; return out ;                /* string?? */
   }

   lout = 4 ;                      /* initialize length of output */
   for( ii=0 ; ii < lstr ; ii++ ){ /* count characters for output */
     switch( str[ii] ){
       case '&':  lout += 5 ; break ;  /* replace '&' with "&amp;" */

       case '<':
       case '>':  lout += 4 ; break ;  /* replace '<' with "&lt;" */

       case '"' :
       case '\'': lout += 6 ; break ;  /* replace '"' with "&quot;" */

       case CR:
       case LF:   lout += 6 ; break ;  /* replace CR with "&#x0d;"
                                                  LF with "&#x0a;" */

       default: lout++ ; break ;      /* copy all other chars */
     }
   }
   out = (char *)calloc(1,lout) ;     /* allocate output string */
   if( !out ){
      fprintf(stderr,"** escapize_string: failed to alloc %d bytes\n",lout);
      return NULL;
   }
   out[0] = '\'' ;                    /* opening quote mark */
   for( ii=0,jj=1 ; ii < lstr ; ii++ ){
      switch( str[ii] ){
         default: out[jj++] = str[ii] ; break ;  /* normal characters */

         case '&':  memcpy(out+jj,"&amp;",5)  ; jj+=5 ; break ;

         case '<':  memcpy(out+jj,"&lt;",4)   ; jj+=4 ; break ;
         case '>':  memcpy(out+jj,"&gt;",4)   ; jj+=4 ; break ;

         case '"' : memcpy(out+jj,"&quot;",6) ; jj+=6 ; break ;

         case '\'': memcpy(out+jj,"&apos;",6) ; jj+=6 ; break ;

         case CR:   memcpy(out+jj,"&#x0d;",6) ; jj+=6 ; break ;
         case LF:   memcpy(out+jj,"&#x0a;",6) ; jj+=6 ; break ;
      }
   }
   out[jj++] = '\''  ;  /* closing quote mark */
   out[jj]   = '\0' ;  /* terminate the string */
   return out ;
}

/*---------------------------------------------------------------------------*/
/*! Dump the information in a NIFTI image header to an XML-ish ASCII string
   that can later be converted back into a NIFTI header in
   nifti_image_from_ascii().

   The resulting string can be free()-ed when you are done with it.
*//*-------------------------------------------------------------------------*/
char *nifti_image_to_ascii( const nifti_image *nim )
{
   char *buf , *ebuf ; int nbuf ;

   if( nim == NULL ) return NULL ;   /* stupid caller */

   buf = (char *)calloc(1,65534); nbuf = 0; /* longer than needed, to be safe */
   if( !buf ){
      fprintf(stderr,"** NITA: failed to alloc %d bytes\n",65534);
      return NULL;
   }

   sprintf( buf , "<nifti_image\n" ) ;   /* XML-ish opener */

   sprintf( buf+strlen(buf) , "  nifti_type = '%s'\n" ,
              (nim->nifti_type == NIFTI_FTYPE_NIFTI1_1) ? "NIFTI-1+"
             :(nim->nifti_type == NIFTI_FTYPE_NIFTI1_2) ? "NIFTI-1"
             :(nim->nifti_type == NIFTI_FTYPE_ASCII   ) ? "NIFTI-1A"
             :                         "ANALYZE-7.5" ) ;

   /** Strings that we don't control (filenames, etc.) that might
       contain "weird" characters (like quotes) are "escaped":
       - A few special characters are replaced by XML-style escapes, using
         the function escapize_string().
       - On input, function unescape_string() reverses this process.
       - The result is that the NIFTI ASCII-format header is XML-compliant. */

   ebuf = escapize_string(nim->fname) ;
   sprintf( buf+strlen(buf) , "  header_filename = %s\n",ebuf); free(ebuf);

   ebuf = escapize_string(nim->iname) ;
   sprintf( buf+strlen(buf) , "  image_filename = %s\n", ebuf); free(ebuf);

   sprintf( buf+strlen(buf) , "  image_offset = '%d'\n" , nim->iname_offset );

                       sprintf( buf+strlen(buf), "  ndim = '%d'\n", nim->ndim);
                       sprintf( buf+strlen(buf), "  nx = '%d'\n",   nim->nx  );
   if( nim->ndim > 1 ) sprintf( buf+strlen(buf), "  ny = '%d'\n",   nim->ny  );
   if( nim->ndim > 2 ) sprintf( buf+strlen(buf), "  nz = '%d'\n",   nim->nz  );
   if( nim->ndim > 3 ) sprintf( buf+strlen(buf), "  nt = '%d'\n",   nim->nt  );
   if( nim->ndim > 4 ) sprintf( buf+strlen(buf), "  nu = '%d'\n",   nim->nu  );
   if( nim->ndim > 5 ) sprintf( buf+strlen(buf), "  nv = '%d'\n",   nim->nv  );
   if( nim->ndim > 6 ) sprintf( buf+strlen(buf), "  nw = '%d'\n",   nim->nw  );
                       sprintf( buf+strlen(buf), "  dx = '%g'\n",   nim->dx  );
   if( nim->ndim > 1 ) sprintf( buf+strlen(buf), "  dy = '%g'\n",   nim->dy  );
   if( nim->ndim > 2 ) sprintf( buf+strlen(buf), "  dz = '%g'\n",   nim->dz  );
   if( nim->ndim > 3 ) sprintf( buf+strlen(buf), "  dt = '%g'\n",   nim->dt  );
   if( nim->ndim > 4 ) sprintf( buf+strlen(buf), "  du = '%g'\n",   nim->du  );
   if( nim->ndim > 5 ) sprintf( buf+strlen(buf), "  dv = '%g'\n",   nim->dv  );
   if( nim->ndim > 6 ) sprintf( buf+strlen(buf), "  dw = '%g'\n",   nim->dw  );

   sprintf( buf+strlen(buf) , "  datatype = '%d'\n" , nim->datatype ) ;
   sprintf( buf+strlen(buf) , "  datatype_name = '%s'\n" ,
                              nifti_datatype_string(nim->datatype) ) ;

   sprintf( buf+strlen(buf) , "  nvox = '%u'\n" , (unsigned)nim->nvox ) ;
   sprintf( buf+strlen(buf) , "  nbyper = '%d'\n" , nim->nbyper ) ;

   sprintf( buf+strlen(buf) , "  byteorder = '%s'\n" ,
            (nim->byteorder==MSB_FIRST) ? "MSB_FIRST" : "LSB_FIRST" ) ;

   if( nim->cal_min < nim->cal_max ){
     sprintf( buf+strlen(buf) , "  cal_min = '%g'\n", nim->cal_min ) ;
     sprintf( buf+strlen(buf) , "  cal_max = '%g'\n", nim->cal_max ) ;
   }

   if( nim->scl_slope != 0.0 ){
     sprintf( buf+strlen(buf) , "  scl_slope = '%g'\n" , nim->scl_slope ) ;
     sprintf( buf+strlen(buf) , "  scl_inter = '%g'\n" , nim->scl_inter ) ;
   }

   if( nim->intent_code > 0 ){
     sprintf( buf+strlen(buf) , "  intent_code = '%d'\n", nim->intent_code ) ;
     sprintf( buf+strlen(buf) , "  intent_code_name = '%s'\n" ,
                                nifti_intent_string(nim->intent_code) ) ;
     sprintf( buf+strlen(buf) , "  intent_p1 = '%g'\n" , nim->intent_p1 ) ;
     sprintf( buf+strlen(buf) , "  intent_p2 = '%g'\n" , nim->intent_p2 ) ;
     sprintf( buf+strlen(buf) , "  intent_p3 = '%g'\n" , nim->intent_p3 ) ;

     if( nim->intent_name[0] != '\0' ){
       ebuf = escapize_string(nim->intent_name) ;
       sprintf( buf+strlen(buf) , "  intent_name = %s\n",ebuf) ;
       free(ebuf) ;
     }
   }

   if( nim->toffset != 0.0 )
     sprintf( buf+strlen(buf) , "  toffset = '%g'\n",nim->toffset ) ;

   if( nim->xyz_units > 0 )
     sprintf( buf+strlen(buf) ,
              "  xyz_units = '%d'\n"
              "  xyz_units_name = '%s'\n" ,
              nim->xyz_units , nifti_units_string(nim->xyz_units) ) ;

   if( nim->time_units > 0 )
     sprintf( buf+strlen(buf) ,
              "  time_units = '%d'\n"
              "  time_units_name = '%s'\n" ,
              nim->time_units , nifti_units_string(nim->time_units) ) ;

   if( nim->freq_dim > 0 )
     sprintf( buf+strlen(buf) , "  freq_dim = '%d'\n",nim->freq_dim ) ;
   if( nim->phase_dim > 0 )
     sprintf( buf+strlen(buf) , "  phase_dim = '%d'\n",nim->phase_dim ) ;
   if( nim->slice_dim > 0 )
     sprintf( buf+strlen(buf) , "  slice_dim = '%d'\n",nim->slice_dim ) ;
   if( nim->slice_code > 0 )
     sprintf( buf+strlen(buf) ,
              "  slice_code = '%d'\n"
              "  slice_code_name = '%s'\n" ,
              nim->slice_code , nifti_slice_string(nim->slice_code) ) ;
   if( nim->slice_start >= 0 && nim->slice_end > nim->slice_start )
     sprintf( buf+strlen(buf) ,
              "  slice_start = '%d'\n"
              "  slice_end = '%d'\n"  , nim->slice_start , nim->slice_end ) ;
   if( nim->slice_duration != 0.0 )
     sprintf( buf+strlen(buf) , "  slice_duration = '%g'\n",
              nim->slice_duration ) ;

   if( nim->descrip[0] != '\0' ){
     ebuf = escapize_string(nim->descrip) ;
     sprintf( buf+strlen(buf) , "  descrip = %s\n",ebuf) ;
     free(ebuf) ;
   }

   if( nim->aux_file[0] != '\0' ){
     ebuf = escapize_string(nim->aux_file) ;
     sprintf( buf+strlen(buf) , "  aux_file = %s\n",ebuf) ;
     free(ebuf) ;
   }

   if( nim->qform_code > 0 ){
     int i,j,k ;

     sprintf( buf+strlen(buf) ,
              "  qform_code = '%d'\n"
              "  qform_code_name = '%s'\n"
     "  qto_xyz_matrix = '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g'\n" ,
         nim->qform_code      , nifti_xform_string(nim->qform_code) ,
         nim->qto_xyz.m[0][0] , nim->qto_xyz.m[0][1] ,
         nim->qto_xyz.m[0][2] , nim->qto_xyz.m[0][3] ,
         nim->qto_xyz.m[1][0] , nim->qto_xyz.m[1][1] ,
         nim->qto_xyz.m[1][2] , nim->qto_xyz.m[1][3] ,
         nim->qto_xyz.m[2][0] , nim->qto_xyz.m[2][1] ,
         nim->qto_xyz.m[2][2] , nim->qto_xyz.m[2][3] ,
         nim->qto_xyz.m[3][0] , nim->qto_xyz.m[3][1] ,
         nim->qto_xyz.m[3][2] , nim->qto_xyz.m[3][3]  ) ;

     sprintf( buf+strlen(buf) ,
     "  qto_ijk_matrix = '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g'\n" ,
         nim->qto_ijk.m[0][0] , nim->qto_ijk.m[0][1] ,
         nim->qto_ijk.m[0][2] , nim->qto_ijk.m[0][3] ,
         nim->qto_ijk.m[1][0] , nim->qto_ijk.m[1][1] ,
         nim->qto_ijk.m[1][2] , nim->qto_ijk.m[1][3] ,
         nim->qto_ijk.m[2][0] , nim->qto_ijk.m[2][1] ,
         nim->qto_ijk.m[2][2] , nim->qto_ijk.m[2][3] ,
         nim->qto_ijk.m[3][0] , nim->qto_ijk.m[3][1] ,
         nim->qto_ijk.m[3][2] , nim->qto_ijk.m[3][3]  ) ;

     sprintf( buf+strlen(buf) ,
              "  quatern_b = '%g'\n"
              "  quatern_c = '%g'\n"
              "  quatern_d = '%g'\n"
              "  qoffset_x = '%g'\n"
              "  qoffset_y = '%g'\n"
              "  qoffset_z = '%g'\n"
              "  qfac = '%g'\n" ,
         nim->quatern_b , nim->quatern_c , nim->quatern_d ,
         nim->qoffset_x , nim->qoffset_y , nim->qoffset_z , nim->qfac ) ;

     nifti_mat44_to_orientation( nim->qto_xyz , &i,&j,&k ) ;
     if( i > 0 && j > 0 && k > 0 )
       sprintf( buf+strlen(buf) ,
                "  qform_i_orientation = '%s'\n"
                "  qform_j_orientation = '%s'\n"
                "  qform_k_orientation = '%s'\n" ,
                nifti_orientation_string(i) ,
                nifti_orientation_string(j) ,
                nifti_orientation_string(k)  ) ;
   }

   if( nim->sform_code > 0 ){
     int i,j,k ;

     sprintf( buf+strlen(buf) ,
              "  sform_code = '%d'\n"
              "  sform_code_name = '%s'\n"
     "  sto_xyz_matrix = '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g'\n" ,
         nim->sform_code      , nifti_xform_string(nim->sform_code) ,
         nim->sto_xyz.m[0][0] , nim->sto_xyz.m[0][1] ,
         nim->sto_xyz.m[0][2] , nim->sto_xyz.m[0][3] ,
         nim->sto_xyz.m[1][0] , nim->sto_xyz.m[1][1] ,
         nim->sto_xyz.m[1][2] , nim->sto_xyz.m[1][3] ,
         nim->sto_xyz.m[2][0] , nim->sto_xyz.m[2][1] ,
         nim->sto_xyz.m[2][2] , nim->sto_xyz.m[2][3] ,
         nim->sto_xyz.m[3][0] , nim->sto_xyz.m[3][1] ,
         nim->sto_xyz.m[3][2] , nim->sto_xyz.m[3][3]  ) ;

     sprintf( buf+strlen(buf) ,
     "  sto_ijk matrix = '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g'\n" ,
         nim->sto_ijk.m[0][0] , nim->sto_ijk.m[0][1] ,
         nim->sto_ijk.m[0][2] , nim->sto_ijk.m[0][3] ,
         nim->sto_ijk.m[1][0] , nim->sto_ijk.m[1][1] ,
         nim->sto_ijk.m[1][2] , nim->sto_ijk.m[1][3] ,
         nim->sto_ijk.m[2][0] , nim->sto_ijk.m[2][1] ,
         nim->sto_ijk.m[2][2] , nim->sto_ijk.m[2][3] ,
         nim->sto_ijk.m[3][0] , nim->sto_ijk.m[3][1] ,
         nim->sto_ijk.m[3][2] , nim->sto_ijk.m[3][3]  ) ;

     nifti_mat44_to_orientation( nim->sto_xyz , &i,&j,&k ) ;
     if( i > 0 && j > 0 && k > 0 )
       sprintf( buf+strlen(buf) ,
                "  sform_i_orientation = '%s'\n"
                "  sform_j_orientation = '%s'\n"
                "  sform_k_orientation = '%s'\n" ,
                nifti_orientation_string(i) ,
                nifti_orientation_string(j) ,
                nifti_orientation_string(k)  ) ;
   }

   sprintf( buf+strlen(buf) , "  num_ext = '%d'\n", nim->num_ext ) ;

   sprintf( buf+strlen(buf) , "/>\n" ) ;   /* XML-ish closer */

   nbuf = strlen(buf) ;
   buf  = (char *)realloc((void *)buf, nbuf+1); /* cut back to proper length */
   if( !buf ) fprintf(stderr,"** NITA: failed to realloc %d bytes\n",nbuf+1);
   return buf ;
}

/*---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*! get the byte order for this CPU

    - LSB_FIRST means least significant byte, first (little endian)
    - MSB_FIRST means most significant byte, first (big endian)
*//*--------------------------------------------------------------------*/
int nifti_short_order(void)   /* determine this CPU's byte order */
{
   union { unsigned char bb[2] ;
           short         ss    ; } fred ;

   fred.bb[0] = 1 ; fred.bb[1] = 0 ;

   return (fred.ss == 1) ? LSB_FIRST : MSB_FIRST ;
}

/*---------------------------------------------------------------------------*/

#undef  QQNUM
#undef  QNUM
#undef  QSTR

/* macro to check lhs string against "n1"; if it matches,
   interpret rhs string as a number, and put it into nim->"n2" */

#define QQNUM(n1,n2) if( strcmp(lhs,#n1)==0 ) nim->n2=strtod(rhs,NULL)

/* same, but where "n1" == "n2" */

#define QNUM(nam)    QQNUM(nam,nam)

/* macro to check lhs string against "nam"; if it matches,
   put rhs string into nim->"nam" string, with max length = "ml" */

#define QSTR(nam,ml) if( strcmp(lhs,#nam) == 0 )                           \
                       strncpy(nim->nam,rhs,ml), nim->nam[ml]='\0'

/*---------------------------------------------------------------------------*/
/*! Take an XML-ish ASCII string and create a NIFTI image header to match.

    NULL is returned if enough information isn't present in the input string.
    - The image data can later be loaded with nifti_image_load().
    - The struct returned here can be liberated with nifti_image_free().
    - Not a lot of error checking is done here to make sure that the
      input values are reasonable!
*//*-------------------------------------------------------------------------*/
nifti_image *nifti_image_from_ascii( const char *str, int * bytes_read )
{
   char lhs[1024] , rhs[1024] ;
   int ii , spos, nn , slen ;
   nifti_image *nim ;              /* will be output */

   if( str == NULL || *str == '\0' ) return NULL ;  /* bad input!? */

   /* scan for opening string */

   spos = 0 ; slen = strlen(str) ;
   ii = sscanf( str+spos , "%1023s%n" , lhs , &nn ) ; spos += nn ;
   if( ii == 0 || strcmp(lhs,"<nifti_image") != 0 ) return NULL ;

   /* create empty image struct */

   nim = (nifti_image *) calloc( 1 , sizeof(nifti_image) ) ;
   if( !nim ){
      fprintf(stderr,"** NIFA: failed to alloc nifti_image\n");
      return NULL;
   }

   nim->nx = nim->ny = nim->nz = nim->nt
           = nim->nu = nim->nv = nim->nw = 1 ;
   nim->dx = nim->dy = nim->dz = nim->dt
           = nim->du = nim->dv = nim->dw = nim->qfac = 1.0 ;

   nim->byteorder = nifti_short_order() ;

   /* starting at str[spos], scan for "equations" of the form
         lhs = 'rhs'
      and assign rhs values into the struct component named by lhs */

   while(1){

     while( isspace(str[spos]) ) spos++ ;  /* skip whitespace */
     if( str[spos] == '\0' ) break ;       /* end of string? */

     /* get lhs string */

     ii = sscanf( str+spos , "%1023s%n" , lhs , &nn ) ; spos += nn ;
     if( ii == 0 || strcmp(lhs,"/>") == 0 ) break ;  /* end of input? */

     /* skip whitespace and the '=' marker */

     while( isspace(str[spos]) || str[spos] == '=' ) spos++ ;
     if( str[spos] == '\0' ) break ;       /* end of string? */

     /* if next character is a quote ', copy everything up to next '
        otherwise, copy everything up to next nonblank              */

     if( str[spos] == '\'' ){
        ii = spos+1 ;
        while( str[ii] != '\0' && str[ii] != '\'' ) ii++ ;
        nn = ii-spos-1 ; if( nn > 1023 ) nn = 1023 ;
        memcpy(rhs,str+spos+1,nn) ; rhs[nn] = '\0' ;
        spos = (str[ii] == '\'') ? ii+1 : ii ;
     } else {
        ii = sscanf( str+spos , "%1023s%n" , rhs , &nn ) ; spos += nn ;
        if( ii == 0 ) break ;  /* nothing found? */
     }
     unescape_string(rhs) ;  /* remove any XML escape sequences */

     /* Now can do the assignment, based on lhs string.
        Start with special cases that don't fit the QNUM/QSTR macros. */

     if( strcmp(lhs,"nifti_type") == 0 ){
            if( strcmp(rhs,"ANALYZE-7.5") == 0 )
               nim->nifti_type = NIFTI_FTYPE_ANALYZE ;
       else if( strcmp(rhs,"NIFTI-1+")    == 0 )
               nim->nifti_type = NIFTI_FTYPE_NIFTI1_1 ;
       else if( strcmp(rhs,"NIFTI-1")     == 0 )
               nim->nifti_type = NIFTI_FTYPE_NIFTI1_2 ;
       else if( strcmp(rhs,"NIFTI-1A")    == 0 )
               nim->nifti_type = NIFTI_FTYPE_ASCII ;
     }
     else if( strcmp(lhs,"header_filename") == 0 ){
       nim->fname = nifti_strdup(rhs) ;
     }
     else if( strcmp(lhs,"image_filename") == 0 ){
       nim->iname = nifti_strdup(rhs) ;
     }
     else if( strcmp(lhs,"sto_xyz_matrix") == 0 ){
       sscanf( rhs , "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f" ,
               &(nim->sto_xyz.m[0][0]) , &(nim->sto_xyz.m[0][1]) ,
               &(nim->sto_xyz.m[0][2]) , &(nim->sto_xyz.m[0][3]) ,
               &(nim->sto_xyz.m[1][0]) , &(nim->sto_xyz.m[1][1]) ,
               &(nim->sto_xyz.m[1][2]) , &(nim->sto_xyz.m[1][3]) ,
               &(nim->sto_xyz.m[2][0]) , &(nim->sto_xyz.m[2][1]) ,
               &(nim->sto_xyz.m[2][2]) , &(nim->sto_xyz.m[2][3]) ,
               &(nim->sto_xyz.m[3][0]) , &(nim->sto_xyz.m[3][1]) ,
               &(nim->sto_xyz.m[3][2]) , &(nim->sto_xyz.m[3][3])  ) ;
     }
     else if( strcmp(lhs,"byteorder") == 0 ){
       if( strcmp(rhs,"MSB_FIRST") == 0 ) nim->byteorder = MSB_FIRST ;
       if( strcmp(rhs,"LSB_FIRST") == 0 ) nim->byteorder = LSB_FIRST ;
     }
     else QQNUM(image_offset,iname_offset) ;
     else QNUM(datatype) ;
     else QNUM(ndim) ;
     else QNUM(nx) ;
     else QNUM(ny) ;
     else QNUM(nz) ;
     else QNUM(nt) ;
     else QNUM(nu) ;
     else QNUM(nv) ;
     else QNUM(nw) ;
     else QNUM(dx) ;
     else QNUM(dy) ;
     else QNUM(dz) ;
     else QNUM(dt) ;
     else QNUM(du) ;
     else QNUM(dv) ;
     else QNUM(dw) ;
     else QNUM(cal_min) ;
     else QNUM(cal_max) ;
     else QNUM(scl_slope) ;
     else QNUM(scl_inter) ;
     else QNUM(intent_code) ;
     else QNUM(intent_p1) ;
     else QNUM(intent_p2) ;
     else QNUM(intent_p3) ;
     else QSTR(intent_name,15) ;
     else QNUM(toffset) ;
     else QNUM(xyz_units) ;
     else QNUM(time_units) ;
     else QSTR(descrip,79) ;
     else QSTR(aux_file,23) ;
     else QNUM(qform_code) ;
     else QNUM(quatern_b) ;
     else QNUM(quatern_c) ;
     else QNUM(quatern_d) ;
     else QNUM(qoffset_x) ;
     else QNUM(qoffset_y) ;
     else QNUM(qoffset_z) ;
     else QNUM(qfac) ;
     else QNUM(sform_code) ;
     else QNUM(freq_dim) ;
     else QNUM(phase_dim) ;
     else QNUM(slice_dim) ;
     else QNUM(slice_code) ;
     else QNUM(slice_start) ;
     else QNUM(slice_end) ;
     else QNUM(slice_duration) ;
     else QNUM(num_ext) ;

   } /* end of while loop */

   if( bytes_read ) *bytes_read = spos+1;         /* "process" last '\n' */

   /* do miscellaneous checking and cleanup */

   if( nim->ndim <= 0 ){ nifti_image_free(nim); return NULL; } /* bad! */

   nifti_datatype_sizes( nim->datatype, &(nim->nbyper), &(nim->swapsize) );
   if( nim->nbyper == 0 ){ nifti_image_free(nim); return NULL; } /* bad! */

   nim->dim[0] = nim->ndim ;
   nim->dim[1] = nim->nx ; nim->pixdim[1] = nim->dx ;
   nim->dim[2] = nim->ny ; nim->pixdim[2] = nim->dy ;
   nim->dim[3] = nim->nz ; nim->pixdim[3] = nim->dz ;
   nim->dim[4] = nim->nt ; nim->pixdim[4] = nim->dt ;
   nim->dim[5] = nim->nu ; nim->pixdim[5] = nim->du ;
   nim->dim[6] = nim->nv ; nim->pixdim[6] = nim->dv ;
   nim->dim[7] = nim->nw ; nim->pixdim[7] = nim->dw ;

   nim->nvox = (size_t)nim->nx * nim->ny * nim->nz
                     * nim->nt * nim->nu * nim->nv * nim->nw ;

   if( nim->qform_code > 0 )
     nim->qto_xyz = nifti_quatern_to_mat44(
                      nim->quatern_b, nim->quatern_c, nim->quatern_d,
                      nim->qoffset_x, nim->qoffset_y, nim->qoffset_z,
                      nim->dx       , nim->dy       , nim->dz       ,
                      nim->qfac                                      ) ;
   else
     nim->qto_xyz = nifti_quatern_to_mat44(
                      0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,
                      nim->dx , nim->dy , nim->dz , 0.0 ) ;


   nim->qto_ijk = nifti_mat44_inverse( nim->qto_xyz ) ;

   if( nim->sform_code > 0 )
     nim->sto_ijk = nifti_mat44_inverse( nim->sto_xyz ) ;

   return nim ;
}


/*---------------------------------------------------------------------------*/
/*! validate the nifti_image

    \return 1 if the structure seems valid, otherwise 0

    \sa nifti_nim_has_valid_dims, nifti_hdr_looks_good
*//*-------------------------------------------------------------------------*/
int nifti_nim_is_valid(nifti_image * nim, int complain)
{
   int errs = 0;

   if( !nim ){
      fprintf(stderr,"** is_valid_nim: nim is NULL\n");
      return 0;
   }

   if( g_opts.debug > 2 ) fprintf(stderr,"-d nim_is_valid check...\n");

   /**- check that dim[] matches the individual values ndim, nx, ny, ... */
   if( ! nifti_nim_has_valid_dims(nim,complain) ){
      if( !complain ) return 0;
      errs++;
   }

   /* might check nbyper, pixdim, q/sforms, swapsize, nifti_type, ... */

   /**- be explicit in return of 0 or 1 */
   if( errs > 0 ) return 0;
   else           return 1;
}

/*---------------------------------------------------------------------------*/
/*! validate nifti dimensions

    \return 1 if valid, 0 if not

    \sa nifti_nim_is_valid, nifti_hdr_looks_good

    rely on dim[] as the master
*//*-------------------------------------------------------------------------*/
int nifti_nim_has_valid_dims(nifti_image * nim, int complain)
{
   size_t prod;
   int    c, errs = 0;

   /**- start with dim[0]: failure here is considered terminal */
   if( nim->dim[0] <= 0 || nim->dim[0] > 7 ){
      errs++;
      if( complain )
         fprintf(stderr,"** NVd: dim[0] (%d) out of range [1,7]\n",nim->dim[0]);
      return 0;
   }

   /**- check whether ndim equals dim[0] */
   if( nim->ndim != nim->dim[0] ){
      errs++;
      if( ! complain ) return 0;
      fprintf(stderr,"** NVd: ndim != dim[0] (%d,%d)\n",nim->ndim,nim->dim[0]);
   }

   /**- compare each dim[i] to the proper nx, ny, ... */
   if( ( (nim->dim[0] >= 1) && (nim->dim[1] != nim->nx) ) ||
       ( (nim->dim[0] >= 2) && (nim->dim[2] != nim->ny) ) ||
       ( (nim->dim[0] >= 3) && (nim->dim[3] != nim->nz) ) ||
       ( (nim->dim[0] >= 4) && (nim->dim[4] != nim->nt) ) ||
       ( (nim->dim[0] >= 5) && (nim->dim[5] != nim->nu) ) ||
       ( (nim->dim[0] >= 6) && (nim->dim[6] != nim->nv) ) ||
       ( (nim->dim[0] >= 7) && (nim->dim[7] != nim->nw) )   ){
      errs++;
      if( !complain ) return 0;
      fprintf(stderr,"** NVd mismatch: dims    = %d,%d,%d,%d,%d,%d,%d\n"
                     "                 nxyz... = %d,%d,%d,%d,%d,%d,%d\n",
                     nim->dim[1], nim->dim[2], nim->dim[3],
                     nim->dim[4], nim->dim[5], nim->dim[6], nim->dim[7],
                     nim->nx, nim->ny, nim->nz,
                     nim->nt, nim->nu, nim->nv, nim->nw );
   }

   if( g_opts.debug > 2 ){
      fprintf(stderr,"-d check dim[%d] =", nim->dim[0]);
      for( c = 0; c < 7; c++ ) fprintf(stderr," %d", nim->dim[c]);
      fputc('\n', stderr);
   }

   /**- check the dimensions, and that their product matches nvox */
   prod = 1;
   for( c = 1; c <= nim->dim[0]; c++ ){
      if( nim->dim[c] > 0)
         prod *= nim->dim[c];
      else if( nim->dim[c] <= 0 ){
         if( !complain ) return 0;
         fprintf(stderr,"** NVd: dim[%d] (=%d) <= 0\n",c, nim->dim[c]);
         errs++;
      }
   }
   if( prod != nim->nvox ){
      if( ! complain ) return 0;
      fprintf(stderr,"** NVd: nvox does not match %d-dim product (%u, %u)\n",
              nim->dim[0], (unsigned)nim->nvox, (unsigned)prod);
      errs++;
   }

   /**- if debug, warn about any remaining dim that is neither 0, nor 1 */
   /*   (values in dims above dim[0] are undefined, as reminded by Cinly
         Ooi and Alle Meije Wink)                   16 Nov 2005 [rickr] */
   if( g_opts.debug > 1 )
      for( c = nim->dim[0]+1; c <= 7; c++ )
         if( nim->dim[c] != 0 && nim->dim[c] != 1 )
            fprintf(stderr,"** NVd warning: dim[%d] = %d, but ndim = %d\n",
                    c, nim->dim[c], nim->dim[0]);

   if( g_opts.debug > 2 )
      fprintf(stderr,"-d nim_has_valid_dims check, errs = %d\n", errs);

   /**- return invalid or valid */
   if( errs > 0 ) return 0;
   else           return 1;
}


/*---------------------------------------------------------------------------*/
/*! read a nifti image, collapsed across dimensions according to dims[8]  <pre>

    This function may be used to read parts of a nifti dataset, such as
    the time series for a single voxel, or perhaps a slice.  It is similar
    to nifti_image_load(), though the passed 'data' parameter is used for
    returning the image, not nim->data.

    \param nim  given nifti_image struct, corresponding to the data file
    \param dims given list of dimensions (see below)
    \param data pointer to data pointer (if *data is NULL, data will be
                allocated, otherwise not)

    Here, dims is an array of 8 ints, similar to nim->dim[8].  While dims[0]
    is unused at this point, the other indices specify which dimensions to
    collapse (and at which index), and which not to collapse.  If dims[i] is
    set to -1, then that entire dimension will be read in, from index 0 to
    index (nim->dim[i] - 1).  If dims[i] >= 0, then only that index will be
    read in (so dims[i] must also be < nim->dim[i]).

    Example: given  nim->dim[8] = { 4, 64, 64, 21, 80, 1, 1, 1 } (4-D dataset)

      if dims[8] = { 0,  5,  4, 17, -1, -1, -1, -1 }
         -> read time series for voxel i,j,k = 5,4,17

      if dims[8] = { 0, -1, -1, -1, 17, -1, -1, -1 }
         -> read single volume at time point 17

    Example: given  nim->dim[8] = { 6, 64, 64, 21, 80, 4, 3, 1 } (6-D dataset)

      if dims[8] = { 0, 5, 4, 17, -1, 2, 1, 0 }
         -> read time series for the voxel i,j,k = 5,4,17, and dim 5,6 = 2,1

      if dims[8] = { 0, 5, 4, -1, -1, 0, 0, 0 }
         -> read time series for slice at i,j = 5,4, and dim 5,6,7 = 0,0,0
            (note that dims[7] is not relevant, but must be 0 or -1)

    If *data is NULL, then *data will be set as a pointer to new memory,
    allocated here for the resulting collapsed image data.

      e.g. { int    dims[8] = { 0,  5,  4, 17, -1, -1, -1, -1 };
             void * data    = NULL;
             ret_val = nifti_read_collapsed_image(nim, dims, &data);
             if( ret_val > 0 ){
                process_time_series(data);
                if( data != NULL ) free(data);
             }
           }

    NOTE: If *data is not NULL, then it will be assumed that it points to
          valid memory, sufficient to hold the results.  This is done for
          speed and possibly repeated calls to this function.

      e.g. { int    dims[8] = { 0,  -1, -1, -1, -1, -1, -1, -1 };
             void * data    = NULL;
             for( zslice = 0; zslice < nzslices; zslice++ ){
                dims[3] = zslice;
                ret_val = nifti_read_collapsed_image(nim, dims, &data);
                if( ret_val > 0 ) process_slice(zslice, data);
             }
             if( data != NULL ) free(data);
           }

    \return
        -  the total number of bytes read, or < 0 on failure
        -  the read and byte-swapped data, in 'data'            </pre>

    \sa nifti_image_read, nifti_image_free, nifti_image_read_bricks
        nifti_image_load
*//*-------------------------------------------------------------------------*/
int nifti_read_collapsed_image( nifti_image * nim, const int dims [8],
                                void ** data )
{
   znzFile fp;
   int     pivots[8], prods[8], nprods; /* sizes are bounded by dims[], so 8 */
   int     c, bytes;

   /** - check pointers for sanity */
   if( !nim || !dims || !data ){
      fprintf(stderr,"** nifti_RCI: bad params %p, %p, %p\n",
              (void *)nim, (void *)dims, (void *)data);
      return -1;
   }

   if( g_opts.debug > 2 ){
      fprintf(stderr,"-d read_collapsed_image:\n        dims =");
      for(c = 0; c < 8; c++) fprintf(stderr," %3d", dims[c]);
      fprintf(stderr,"\n   nim->dims =");
      for(c = 0; c < 8; c++) fprintf(stderr," %3d", nim->dim[c]);
      fputc('\n', stderr);
   }

   /** - verify that dim[] makes sense */
   if( ! nifti_nim_is_valid(nim, g_opts.debug > 0) ){
      fprintf(stderr,"** invalid nim (file is '%s')\n", nim->fname );
      return -1;
   }

   /** - verify that dims[] makes sense for this dataset */
   for( c = 1; c <= nim->dim[0]; c++ ){
      if( dims[c] >= nim->dim[c] ){
         fprintf(stderr,"** nifti_RCI: dims[%d] >= nim->dim[%d] (%d,%d)\n",
                 c, c, dims[c], nim->dim[c]);
         return -1;
      }
   }

   /** - prepare pivot list - pivots are fixed indices */
   if( make_pivot_list(nim, dims, pivots, prods, &nprods) < 0 ) return -1;

   bytes = rci_alloc_mem(data, prods, nprods, nim->nbyper);
   if( bytes < 0 ) return -1;

   /** - open the image file for reading at the appropriate offset */
   fp = nifti_image_load_prep( nim );
   if( ! fp ){ free(*data);  *data = NULL;  return -1; }     /* failure */

   /** - call the recursive reading function, passing nim, the pivot info,
         location to store memory, and file pointer and position */
   c = rci_read_data(nim, pivots,prods,nprods,dims,
                     (char *)*data, fp, znztell(fp));
   znzclose(fp);   /* in any case, close the file */
   if( c < 0 ){ free(*data);  *data = NULL;  return -1; }    /* failure */

   if( g_opts.debug > 1 )
      fprintf(stderr,"+d read %d bytes of collapsed image from %s\n",
              bytes, nim->fname);

   return bytes;
}


/* local function to find strides per dimension. assumes 7D size and
** stride array.
*/
static void
compute_strides(int *strides,const int *size,int nbyper)
{
  int i;
  strides[0] = nbyper;
  for(i = 1; i < 7; i++)
    {
    strides[i] = size[i-1] * strides[i-1];
    }
}

/*---------------------------------------------------------------------------*/
/*! read an arbitrary subregion from a nifti image

    This function may be used to read a single arbitary subregion of any
    rectangular size from a nifti dataset, such as a small 5x5x5 subregion
    around the center of a 3D image.

    \param nim  given nifti_image struct, corresponding to the data file
    \param start_index the index location of the first voxel that will be returned
    \param region_size the size of the subregion to be returned
    \param data pointer to data pointer (if *data is NULL, data will be
                allocated, otherwise not)

    Example: given  nim->dim[8] = {3, 64, 64, 64, 1, 1, 1, 1 } (3-D dataset)

      if start_index[7] = { 29,  29, 29, 0, 0, 0, 0 } and
         region_size[7] = {  5,   5,  5, 1, 1, 1, 1 }
         -> read 5x5x5 region starting with the first voxel location at (29,29,29)

    NOTE: If *data is not NULL, then it will be assumed that it points to
          valid memory, sufficient to hold the results.  This is done for
          speed and possibly repeated calls to this function.
    \return
        -  the total number of bytes read, or < 0 on failure
        -  the read and byte-swapped data, in 'data'            </pre>

    \sa nifti_image_read, nifti_image_free, nifti_image_read_bricks
        nifti_image_load, nifti_read_collapsed_image
*//*-------------------------------------------------------------------------*/
int nifti_read_subregion_image( nifti_image * nim,
                                int *start_index,
                                int *region_size,
                                void ** data )
{
  znzFile fp;                   /* file to read */
  int i,j,k,l,m,n;              /* indices for dims */
  long int bytes = 0;           /* total # bytes read */
  int total_alloc_size;         /* size of buffer allocation */
  char *readptr;                /* where in *data to read next */
  int strides[7];               /* strides between dimensions */
  int collapsed_dims[8];        /* for read_collapsed_image */
  int *image_size;              /* pointer to dimensions in header */
  long int initial_offset;
  long int offset;              /* seek offset for reading current row */

  /* probably ignored, but set to ndim for consistency*/
  collapsed_dims[0] = nim->ndim;

  /* build a dims array for collapsed image read */
  for(i = 0; i < nim->ndim; i++)
    {
    /* if you take the whole extent in this dimension */
    if(start_index[i] == 0 &&
       region_size[i] == nim->dim[i+1])
      {
      collapsed_dims[i+1] = -1;
      }
    /* if you specify a single element in this dimension */
    else if(region_size[i] == 1)
      {
      collapsed_dims[i+1] = start_index[i];
      }
    else
      {
      collapsed_dims[i+1] = -2; /* sentinel value */
      }
    }
  /* fill out end of collapsed_dims */
  for(i = nim->ndim ; i < 7; i++)
    {
    collapsed_dims[i+1] = -1;
    }

  /* check to see whether collapsed read is possible */
  for(i = 1; i <= nim->ndim; i++)
    {
    if(collapsed_dims[i] == -2)
      {
      break;
      }
    }

  /* if you get through all the dimensions without hitting 
  ** a subrange of size > 1, a collapsed read is possible
  */
  if(i > nim->ndim)
    {
    return nifti_read_collapsed_image(nim, collapsed_dims, data);
    }

  /* point past first element of dim, which holds nim->ndim */
  image_size = &(nim->dim[1]);

  /* check region sizes for sanity */
  for(i = 0; i < nim->ndim; i++)
    {
    if(start_index[i]  + region_size[i] > image_size[i])
      {
      if(g_opts.debug > 1)
        {
        fprintf(stderr,"region doesn't fit within image size\n");
        }
      return -1;
      }
    }

  /* get the file open */
  fp = nifti_image_load_prep( nim );
  /* the current offset is just past the nifti header, save
   * location so that SEEK_SET can be used below
   */
  initial_offset = znztell(fp);  
  /* get strides*/
  compute_strides(strides,image_size,nim->nbyper);

  total_alloc_size = nim->nbyper; /* size of pixel */

  /* find alloc size */
  for(i = 0; i < nim->ndim; i++)
    {
    total_alloc_size *= region_size[i];
    }
  /* allocate buffer, if necessary */
  if(*data == 0)
    {
    *data = (void *)malloc(total_alloc_size);
    }

  if(*data == 0)
    {
    if(g_opts.debug > 1)
      {
      fprintf(stderr,"allocation of %d bytes failed\n",total_alloc_size);
      return -1;
      }
    }

  /* point to start of data buffer as char * */
  readptr = *((char **)data);
  {
  /* can't assume that start_index and region_size have any more than
  ** nim->ndim elements so make local copies, filled out to seven elements
  */
  int si[7], rs[7];
  for(i = 0; i < nim->ndim; i++)
    {
    si[i] = start_index[i];
    rs[i] = region_size[i];
    }
  for(i = nim->ndim; i < 7; i++)
    {
    si[i] = 0;
    rs[i] = 1;
    }
  /* loop through subregion and read a row at a time */
  for(i = si[6]; i < (si[6] + rs[6]); i++)
    {
    for(j = si[5]; j < (si[5] + rs[5]); j++)
      {
      for(k = si[4]; k < (si[4] + rs[4]); k++)
        {
        for(l = si[3]; l < (si[3] + rs[3]); l++)
          {
          for(m = si[2]; m < (si[2] + rs[2]); m++)
            {
            for(n = si[1]; n < (si[1] + rs[1]); n++)
              {
              int nread,read_amount;
              offset = initial_offset +
                (i * strides[6]) +
                (j * strides[5]) +
                (k * strides[4]) +
                (l * strides[3]) +
                (m * strides[2]) +
                (n * strides[1]) +
                (si[0] * strides[0]);
              znzseek(fp, offset, SEEK_SET); /* seek to current row */
              read_amount = rs[0] * nim->nbyper; /* read a row of the subregion*/
              nread = nifti_read_buffer(fp, readptr, read_amount, nim);
              if(nread != read_amount)
                {
                if(g_opts.debug > 1)
                  {
                  fprintf(stderr,"read of %d bytes failed\n",read_amount);
                  return -1;
                  }
                }
              bytes += nread;
              readptr += read_amount;
              }
            }
          }
        }
      }
    }
  }
  return bytes;
}


/* read the data from the file pointed to by fp

   - this a recursive function, so start with the base case
   - data is now (char *) for easy incrementing

   return 0 on success, < 0 on failure
*/
static int rci_read_data(nifti_image * nim, int * pivots, int * prods,
      int nprods, const int dims[], char * data, znzFile fp, size_t base_offset)
{
   size_t sublen, offset, read_size;
   int    c;

   /* bad check first - base_offset may not have been checked */
   if( nprods <= 0 ){
      fprintf(stderr,"** rci_read_data, bad prods, %d\n", nprods);
      return -1;
   }

   /* base case: actually read the data */
   if( nprods == 1 ){
      size_t nread, bytes;

      /* make sure things look good here */
      if( *pivots != 0 ){
         fprintf(stderr,"** rciRD: final pivot == %d!\n", *pivots);
         return -1;
      }

      /* so just seek and read (prods[0] * nbyper) bytes from the file */
      znzseek(fp, base_offset, SEEK_SET);
      bytes = (size_t)prods[0] * nim->nbyper;
      nread = nifti_read_buffer(fp, data, bytes, nim);
      if( nread != bytes ){
         fprintf(stderr,"** rciRD: read only %u of %u bytes from '%s'\n",
                 (unsigned)nread, (unsigned)bytes, nim->fname);
         return -1;
      } else if( g_opts.debug > 3 )
         fprintf(stderr,"+d successful read of %u bytes at offset %u\n",
                 (unsigned)bytes, (unsigned)base_offset);

      return 0;  /* done with base case - return success */
   }

   /* not the base case, so do a set of reduced reads */

   /* compute size of sub-brick: all dimensions below pivot */
   for( c = 1, sublen = 1; c < *pivots; c++ ) sublen *= nim->dim[c];

   /* compute number of values to read, i.e. remaining prods */
   for( c = 1, read_size = 1; c < nprods; c++ ) read_size *= prods[c];
   read_size *= nim->nbyper;  /* and multiply by bytes per voxel */

   /* now repeatedly compute offsets, and recursively read */
   for( c = 0; c < prods[0]; c++ ){
      /* offset is (c * sub-block size (including pivot dim))   */
      /*         + (dims[] index into pivot sub-block)          */
      /* the unneeded multiplication is to make this more clear */
      offset = (size_t)c * sublen * nim->dim[*pivots] +
               (size_t)sublen * dims[*pivots];
      offset *= nim->nbyper;

      if( g_opts.debug > 3 )
         fprintf(stderr,"-d reading %u bytes, foff %u + %u, doff %u\n",
                 (unsigned)read_size, (unsigned)base_offset, (unsigned)offset,
                 (unsigned)(c*read_size));

      /* now read the next level down, adding this offset */
      if( rci_read_data(nim, pivots+1, prods+1, nprods-1, dims,
                    data + c * read_size, fp, base_offset + offset) < 0 )
         return -1;
   }

   return 0;
}


/* allocate memory for all collapsed image data

   If *data is already set, do not allocate, but still calculate
   size for debug report.

   return total size on success, and < 0 on failure
*/
static int rci_alloc_mem(void ** data, int prods[8], int nprods, int nbyper )
{
   int size, index;

   if( nbyper < 0 || nprods < 1 || nprods > 8 ){
      fprintf(stderr,"** rci_am: bad params, %d, %d\n", nbyper, nprods);
      return -1;
   }

   for( index = 0, size = 1; index < nprods; index++ )
       size *= prods[index];

   size *= nbyper;

   if( ! *data ){   /* then allocate what is needed */
      if( g_opts.debug > 1 )
         fprintf(stderr,"+d alloc %d (= %d x %d) bytes for collapsed image\n",
                 size, size/nbyper, nbyper);

      *data = malloc(size);   /* actually allocate the memory */
      if( ! *data ){
         fprintf(stderr,"** rci_am: failed to alloc %d bytes for data\n", size);
         return -1;
      }
   } else if( g_opts.debug > 1 )
      fprintf(stderr,"-d rci_am: *data already set, need %d (%d x %d) bytes\n",
              size, size/nbyper, nbyper);

   return size;
}


/* prepare a pivot list for reading

   The pivot points are the indices into dims where the calling function
   wants to collapse a dimension.  The last pivot should always be zero
   (note that we have space for that in the lists).
*/
static int make_pivot_list(nifti_image * nim, const int dims[], int pivots[],
                                              int prods[], int * nprods )
{
   int len, index;

   len = 0;
   index = nim->dim[0];
   while( index > 0 ){
      prods[len] = 1;
      while( index > 0 && (nim->dim[index] == 1 || dims[index] == -1) ){
         prods[len] *= nim->dim[index];
         index--;
      }
      pivots[len] = index;
      len++;
      index--;  /* fine, let it drop out at -1 */
   }

   /* make sure to include 0 as a pivot (instead of just 1, if it is) */
   if( pivots[len-1] != 0 ){
      pivots[len] = 0;
      prods[len] = 1;
      len++;
   }

   *nprods = len;

   if( g_opts.debug > 2 ){
      fprintf(stderr,"+d pivot list created, pivots :");
      for(index = 0; index < len; index++) fprintf(stderr," %d", pivots[index]);
      fprintf(stderr,", prods :");
      for(index = 0; index < len; index++) fprintf(stderr," %d", prods[index]);
      fputc('\n',stderr);
   }

   return 0;
}


#undef ISEND
#define ISEND(c) ( (c)==']' || (c)=='}' || (c)=='\0' )

/*---------------------------------------------------------------------*/
/*! Get an integer list in the range 0..(nvals-1), from the
   character string str.  If we call the output pointer fred,
   then fred[0] = number of integers in the list (> 0), and
        fred[i] = i-th integer in the list for i=1..fred[0].
   If on return, fred == NULL or fred[0] == 0, then something is
   wrong, and the caller must deal with that.

   Syntax of input string:
     - initial '{' or '[' is skipped, if present
     - ends when '}' or ']' or end of string is found
     - contains entries separated by commas
     - entries have one of these forms:
       - a single number
       - a dollar sign '$', which means nvals-1
       - a sequence of consecutive numbers in the form "a..b" or
         "a-b", where "a" and "b" are single numbers (or '$')
       - a sequence of evenly spaced numbers in the form
         "a..b(c)" or "a-b(c)", where "c" encodes the step
     - Example:  "[2,7..4,3..9(2)]" decodes to the list
         2 7 6 5 4 3 5 7 9
     - entries should be in the range 0..nvals-1

   (borrowed, with permission, from thd_intlist.c)
*//*-------------------------------------------------------------------*/
int * nifti_get_intlist( int nvals , const char * str )
{
   int *subv = NULL ;
   int ii , ipos , nout , slen ;
   int ibot,itop,istep , nused ;
   char *cpt ;

   /* Meaningless input? */
   if( nvals < 1 ) return NULL ;

   /* No selection list? */
   if( str == NULL || str[0] == '\0' ) return NULL ;

   /* skip initial '[' or '{' */
   subv    = (int *) malloc( sizeof(int) * 2 ) ;
   subv[0] = nout = 0 ;

   ipos = 0 ;
   if( str[ipos] == '[' || str[ipos] == '{' ) ipos++ ;

   if( g_opts.debug > 1 )
      fprintf(stderr,"-d making int_list (vals = %d) from '%s'\n", nvals, str);

   /**- for each sub-selector until end of input... */

   slen = strlen(str) ;
   while( ipos < slen && !ISEND(str[ipos]) ){

      while( isspace(str[ipos]) ) ipos++ ;   /* skip blanks */
      if( ISEND(str[ipos]) ) break ;         /* done */

      /**- get starting value */

      if( str[ipos] == '$' ){  /* special case */
         ibot = nvals-1 ; ipos++ ;
      } else {                 /* decode an integer */
         ibot = strtol( str+ipos , &cpt , 10 ) ;
         if( ibot < 0 ){
           fprintf(stderr,"** ERROR: list index %d is out of range 0..%d\n",
                   ibot,nvals-1) ;
           free(subv) ; return NULL ;
         }
         if( ibot >= nvals ){
           fprintf(stderr,"** ERROR: list index %d is out of range 0..%d\n",
                   ibot,nvals-1) ;
           free(subv) ; return NULL ;
         }
         nused = (cpt-(str+ipos)) ;
         if( ibot == 0 && nused == 0 ){
           fprintf(stderr,"** ERROR: list syntax error '%s'\n",str+ipos) ;
           free(subv) ; return NULL ;
         }
         ipos += nused ;
      }

      while( isspace(str[ipos]) ) ipos++ ;   /* skip blanks */

      /**- if that's it for this sub-selector, add one value to list */

      if( str[ipos] == ',' || ISEND(str[ipos]) ){
         nout++ ;
         subv = (int *) realloc( (char *)subv , sizeof(int) * (nout+1) ) ;
         subv[0]    = nout ;
         subv[nout] = ibot ;
         if( ISEND(str[ipos]) ) break ; /* done */
         ipos++ ; continue ;            /* re-start loop at next sub-selector */
      }

      /**- otherwise, must have '..' or '-' as next inputs */

      if( str[ipos] == '-' ){
         ipos++ ;
      } else if( str[ipos] == '.' && str[ipos+1] == '.' ){
         ipos++ ; ipos++ ;
      } else {
         fprintf(stderr,"** ERROR: index list syntax is bad: '%s'\n",
                 str+ipos) ;
         free(subv) ; return NULL ;
      }

      /**- get ending value for loop now */

      if( str[ipos] == '$' ){  /* special case */
         itop = nvals-1 ; ipos++ ;
      } else {                 /* decode an integer */
         itop = strtol( str+ipos , &cpt , 10 ) ;
         if( itop < 0 ){
           fprintf(stderr,"** ERROR: index %d is out of range 0..%d\n",
                   itop,nvals-1) ;
           free(subv) ; return NULL ;
         }
         if( itop >= nvals ){
           fprintf(stderr,"** ERROR: index %d is out of range 0..%d\n",
                   itop,nvals-1) ;
           free(subv) ; return NULL ;
         }
         nused = (cpt-(str+ipos)) ;
         if( itop == 0 && nused == 0 ){
           fprintf(stderr,"** ERROR: index list syntax error '%s'\n",str+ipos) ;
           free(subv) ; return NULL ;
         }
         ipos += nused ;
      }

      /**- set default loop step */

      istep = (ibot <= itop) ? 1 : -1 ;

      while( isspace(str[ipos]) ) ipos++ ;                  /* skip blanks */

      /**- check if we have a non-default loop step */

      if( str[ipos] == '(' ){  /* decode an integer */
         ipos++ ;
         istep = strtol( str+ipos , &cpt , 10 ) ;
         if( istep == 0 ){
           fprintf(stderr,"** ERROR: index loop step is 0!\n") ;
           free(subv) ; return NULL ;
         }
         nused = (cpt-(str+ipos)) ;
         ipos += nused ;
         if( str[ipos] == ')' ) ipos++ ;
         if( (ibot-itop)*istep > 0 ){
           fprintf(stderr,"** WARNING: index list '%d..%d(%d)' means nothing\n",
                   ibot,itop,istep ) ;
         }
      }

      /**- add values to output */

      for( ii=ibot ; (ii-itop)*istep <= 0 ; ii += istep ){
         nout++ ;
         subv = (int *) realloc( (char *)subv , sizeof(int) * (nout+1) ) ;
         subv[0]    = nout ;
         subv[nout] = ii ;
      }

      /**- check if we have a comma to skip over */

      while( isspace(str[ipos]) ) ipos++ ;                  /* skip blanks */
      if( str[ipos] == ',' ) ipos++ ;                       /* skip commas */

   }  /* end of loop through selector string */

   if( g_opts.debug > 1 ) {
      fprintf(stderr,"+d int_list (vals = %d): ", subv[0]);
      for( ii = 1; ii <= subv[0]; ii++ ) fprintf(stderr,"%d ", subv[ii]);
      fputc('\n',stderr);
   }

   if( subv[0] == 0 ){ free(subv); subv = NULL; }
   return subv ;
}
