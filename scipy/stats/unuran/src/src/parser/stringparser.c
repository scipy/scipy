/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      stringparser.c.                                              *
 *                                                                           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *                                                                           *
 *   The function unur_str2gen() takes a character string as its argument.   *
 *   This string is parsed in the following ways and the information         *
 *   obtained is used to make a generator object. (It returns NULL when      *
 *   it is not possible to make a generator object either due to a syntax    *
 *   error in the given string, or due to invalid data.                      *
 *                                                                           *
 *   First the given string is copied into a working string. Moreover letter *
 *   letters are transformed into their lower case counterparts and spaces   *
 *   are eliminated (i.e. the working string does not contain any white      *
 *   space characters any more.) Then the working string is splitted into    *
 *   several blocks, separated by ampersands '&'. Each block consists of     *
 *   <key>=<value> pairs, separated by semicolons ';'. The first key in      *
 *   each block is used to indicate the corresponding block.                 *
 *                                                                           *
 *   We have three different blocks (yet):                                   *
 *      distr  ... contains the definition of the distribution               *
 *      method ... contains the description of the transformation methdod    *
 *      urng   ... contains the definition of the uniform RNG                *
 *                 (currently it only supports URNGs from the PRNG library)  *
 *                                                                           *
 *   The 'distr' block must be the very first block and is obligatory.       *
 *   For that reason the keyword "distr" is optional and can be omitted      *
 *   (together with the '=' character). Moreover it is almost ignored while  *
 *   parsing the string. To avoid some possible confussion it only has to    *
 *   start with the letter 'd' (if it is given at all). E.g. "dnormal" can   *
 *   be used as keyword.                                                     *
 *   The value of the "distr" key is used to get a distribution object       *
 *   either via a unur_distr_<value>() call for a standard distribution      *
 *   or the unur_distr_<value>_new() calls are used to get an object for a   *
 *   generic distribution. However not all generic distributions are         *
 *   supported yet (see the list %UNSUPPORTED_DISTR_TYPES in file            *
 *   make_stringparser.pl).                                                  *
 *   The parameters for the standard distribution are given as a list of     *
 *   numbers enclosed in parenthesis. There must not be any character (other *
 *   than white space) between the name of the standard distribution and     *
 *   this list.                                                              *
 *                                                                           *
 *   All the other blocks are optional and can be arranged in arbitrary      *
 *   order.                                                                  *
 *   The value of the "method" key is the name of the transformation method  *
 *   and is used to execute the unur_<value>_new() call.                     *
 *                                                                           *
 *   The value of the "urng" key is passed to the PRNG interface.            *
 *   However, it only works when using the PRNG library is enabled           *
 *   --with-urng-prng flag and libprng is linked by the executable.          *
 *                                                                           *
 *   In each block consecuting <key>=<value> pairs, separated by semicolons  *
 *   ';', are used to set parameters. The name of the parameter is given as  *
 *   key, the argument list is given as value. The arguments of the          *
 *   corresponding set calls are given as tokens, separated by commata ','.  *
 *   There are two types of tokens:                                          *
 *      single tokens, that represent numbers, and                           *
 *      list, i.e. a list of single tokens, separated by commata ',',        *
 *            enclosed in parenthesis.                                       *
 *   The value (including the character '=') can be omitted when no argument *
 *   are required.                                                           *
 *   The key is deduced from the UNU.RAN set calls, such that it uses only   *
 *   the part of the set call beyond "..._set_".                             *
 *   Not all set commands are supported yet (see the output to STDERR of     *
 *   the make_stringparser.pl script).                                       *
 *   As a rule of thumb only those set calls are supported that only have    *
 *   numbers or arrays of numbers as arguments (besides the pointer to the   *
 *   distribution or parameter object).                                      *
 *   (See the code of make_stringparser.pl for details.)                     *
 *                                                                           *
 *   There are additional transformations before executing the necessary     *
 *   set calls (listed by argument types)                                    *
 *      int:      'true' and 'on' are transformed to 1,                      *
 *                'false' and 'off' are transformed to 0.                    *
 *                a missing argument is transformed to 1.                    *
 *      int, int: Instead of two integers, a list with at least two          *
 *                numbers can be given.                                      *
 *                'inf[...]' is transformed to INT_MAX,                      *
 *                '-inf[...]' is transformed to INT_MIN.                     *
 *      unsigned: string is interpreted as hexadecimal number.               *
 *      double:   'inf[...]' is transformed to INFINITY,                     *
 *                '-inf[...]' is transformed to -INFINITY.                   *
 *      double, double: Instead of two doubles, a list with at least two     *
 *                numbers can be given.                                      *
 *                'inf[...]' is transformed to INFINITY.                     *
 *      int, double*: If the first argument is missing, it is replaced       *
 *                by the size of the array.                                  *
 *                If the second argument is missing, the NULL pointer is     *
 *                used instead an array as argument.                         *
 *      double*, int: Only for distributions!                                *
 *                If the second argument 'int' is missing, it is replaced    *
 *                by the size of the array.                                  *
 *   Important: there is no support for 'infinity' in lists!                 * 
 *   (See also the source below for more details.)                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   REMARK:                                                                 *
 *   The stringparser always uses the default debugging flag.                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

#include <stdarg.h>
#include <ctype.h>
#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen.h>
#include "parser_source.h"
#include "parser.h"

#include <urng/urng.h>

#include <distr/distr_source.h>
#include <distr/distr.h>
#include <distr/cemp.h>
#include <distr/cont.h>
#include <distr/corder.h>
#include <distr/cvemp.h>
#include <distr/discr.h>

#include <distributions/unur_distributions.h>

#include <methods/arou.h>
#include <methods/ars.h>
#include <methods/auto.h>
#include <methods/cstd.h>
#include <methods/dari.h>
#include <methods/dau.h>
#include <methods/dgt.h>
#include <methods/dsrou.h>
#include <methods/dss.h>
#include <methods/dstd.h>
#include <methods/empk.h>
#include <methods/empl.h>
#include <methods/gibbs.h>
#include <methods/hinv.h>
#include <methods/hist.h>
#include <methods/hitro.h>
#include <methods/hrb.h>
#include <methods/hrd.h>
#include <methods/hri.h>
#include <methods/itdr.h>
#include <methods/mcorr.h>
#include <methods/mvstd.h>
#include <methods/mvtdr.h>
#include <methods/ninv.h>
#include <methods/norta.h>
#include <methods/nrou.h>
#include <methods/pinv.h>
#include <methods/srou.h>
#include <methods/ssr.h>
#include <methods/tabl.h>
#include <methods/tdr.h>
#include <methods/unif.h>
#include <methods/utdr.h>
#include <methods/vempk.h>
#include <methods/vnrou.h>

#if defined(UNUR_URNG_UNURAN) && defined(UNURAN_HAS_PRNG)
#include <uniform/urng_prng.h>
#endif

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "STRING"       /* (pseudo) type of generator                 */

/*---------------------------------------------------------------------------*/

static struct unur_distr *_unur_str_distr( char *str_distr );
/*---------------------------------------------------------------------------*/
/* get distribution object for given distribution.                           */
/*---------------------------------------------------------------------------*/

static struct unur_distr *_unur_str_distr_new( char *distribution );
/*---------------------------------------------------------------------------*/
/* get new distribution object.                                              */
/*---------------------------------------------------------------------------*/

static struct unur_distr *_unur_str_distr_make_os( UNUR_DISTR *distr, 
						   const char *key,
						   char *type_args, char **args );
/*---------------------------------------------------------------------------*/
/* Make distribution object for order statistics of given distribution.      */
/*---------------------------------------------------------------------------*/

static int _unur_str_distr_set( UNUR_DISTR **ptr_distr, const char *key, char *value );
/*---------------------------------------------------------------------------*/
/* set parameters for distribution.                                          */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Set parameter in given distribution object.                               */
/*---------------------------------------------------------------------------*/
typedef int distr_set_i( UNUR_DISTR *distr, int i );
typedef int distr_set_ii( UNUR_DISTR *distr, int i1, int i2 );
typedef int distr_set_d( UNUR_DISTR *distr, double d );
typedef int distr_set_dd( UNUR_DISTR *distr, double d1, double d2 );
typedef int distr_set_Di( UNUR_DISTR *distr, const double *array, int size );
typedef int distr_set_C( UNUR_DISTR *distr, const char *string );
/*---------------------------------------------------------------------------*/
/* types of set calls                                                        */
/*---------------------------------------------------------------------------*/

static int _unur_str_distr_set_i( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				  distr_set_i set );
/*---------------------------------------------------------------------------*/
/* integer.                                                                  */
/*---------------------------------------------------------------------------*/

static int _unur_str_distr_set_ii( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				   distr_set_ii set );
/*---------------------------------------------------------------------------*/
/* 2 integers.                                                               */
/*---------------------------------------------------------------------------*/

static int _unur_str_distr_set_d( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				  distr_set_d set );
/*---------------------------------------------------------------------------*/
/* double.                                                                   */
/*---------------------------------------------------------------------------*/

static int _unur_str_distr_set_dd( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				   distr_set_dd set );
/*---------------------------------------------------------------------------*/
/* 2 doubles.                                                                */
/*---------------------------------------------------------------------------*/

static int _unur_str_distr_set_Di( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				   distr_set_Di set );
/*---------------------------------------------------------------------------*/
/* array of doubles & integer.                                               */
/*---------------------------------------------------------------------------*/

static int _unur_str_distr_set_C( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				  distr_set_C set );
/*---------------------------------------------------------------------------*/
/* constant string.                                                          */
/*---------------------------------------------------------------------------*/


static struct unur_par *_unur_str_par( char *str_method, const UNUR_DISTR *distr,
				       struct unur_slist *mlist );
/*---------------------------------------------------------------------------*/
/* get parameter object for distribution and method.                         */
/*---------------------------------------------------------------------------*/

static struct unur_par *_unur_str_par_new( const char *method, const UNUR_DISTR *distr );
/*---------------------------------------------------------------------------*/
/* get new parameter object for method.                                      */
/*---------------------------------------------------------------------------*/

static int _unur_str_par_set( UNUR_PAR *par, const char *key, char *value,
			      struct unur_slist *mlist );
/*---------------------------------------------------------------------------*/
/* set parameters for method.                                                */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Set parameter in given parameter object.                                 */
/*---------------------------------------------------------------------------*/
typedef int par_set_void( UNUR_PAR *par );
typedef int par_set_i( UNUR_PAR *par, int i );
typedef int par_set_ii( UNUR_PAR *par, int i1, int i2 );
typedef int par_set_u( UNUR_PAR *par, unsigned u );
typedef int par_set_d( UNUR_PAR *par, double d );
typedef int par_set_dd( UNUR_PAR *par, double d1, double d2 );
typedef int par_set_iD( UNUR_PAR *par, int size, const double *array );
typedef int par_set_Di( UNUR_PAR *par, const double *array, int size );
/*---------------------------------------------------------------------------*/
/* types of set calls                                                        */
/*---------------------------------------------------------------------------*/

static int _unur_str_par_set_void( UNUR_PAR *par, const char *key, char *type_args, char **args,
				   par_set_void set );
/*---------------------------------------------------------------------------*/
/* void (no argument required).                                              */
/*---------------------------------------------------------------------------*/

static int _unur_str_par_set_i( UNUR_PAR *par, const char *key, char *type_args, char **args,
				par_set_i set );
/*---------------------------------------------------------------------------*/
/* integer.                                                                  */
/*---------------------------------------------------------------------------*/

static int _unur_str_par_set_ii( UNUR_PAR *par, const char *key, char *type_args, char **args,
				   par_set_ii set );
/*---------------------------------------------------------------------------*/
/* 2 integers.                                                               */
/*---------------------------------------------------------------------------*/

static int _unur_str_par_set_u( UNUR_PAR *par, const char *key, char *type_args, char **args,
				par_set_u set );
/*---------------------------------------------------------------------------*/
/* unsigned integer.                                                         */
/*---------------------------------------------------------------------------*/

static int _unur_str_par_set_d( UNUR_PAR *par, const char *key, char *type_args, char **args,
				par_set_d set );
/*---------------------------------------------------------------------------*/
/* double.                                                                   */
/*---------------------------------------------------------------------------*/

static int _unur_str_par_set_dd( UNUR_PAR *par, const char *key, char *type_args, char **args,
				 par_set_dd set );

/*---------------------------------------------------------------------------*/
/* 2 doubles.                                                                */
/*---------------------------------------------------------------------------*/

static int _unur_str_par_set_iD( UNUR_PAR *par, const char *key, char *type_args, char **args,
				 par_set_iD set, struct unur_slist *mlist );
/*---------------------------------------------------------------------------*/
/* integer & array of doubles.                                               */
/*---------------------------------------------------------------------------*/

static int _unur_str_par_set_Di( UNUR_PAR *par, const char *key, char *type_args, char **args,
				 par_set_Di set, struct unur_slist *mlist );
/*---------------------------------------------------------------------------*/
/* array of doubles & integer.                                               */
/*---------------------------------------------------------------------------*/


static UNUR_URNG *_unur_str2urng( char *str_urng );
/*---------------------------------------------------------------------------*/
/* get uniform RNG.                                                          */
/*---------------------------------------------------------------------------*/


static int _unur_str_set_args( char *value, char *type_args, char **args, int max_args );
/*---------------------------------------------------------------------------*/
/* parse argument string for set call.                                       */
/*---------------------------------------------------------------------------*/

static int _unur_parse_ilist( char *liststr, int **iarray );
/*---------------------------------------------------------------------------*/
/* Process a comma separated list of integers.                               */
/*---------------------------------------------------------------------------*/

static int _unur_parse_dlist( char *liststr, double **darray );
/*---------------------------------------------------------------------------*/
/* Process a comma separated list of numbers.                                */
/*---------------------------------------------------------------------------*/

static int _unur_atoi ( const char *str );
/*---------------------------------------------------------------------------*/
/* Convert string into its integer representation.                           */
/*---------------------------------------------------------------------------*/

static unsigned _unur_atou ( const char *str );
/*---------------------------------------------------------------------------*/
/* Convert string into its unsigned representation (using hexadecimal).      */
/*---------------------------------------------------------------------------*/

static double _unur_atod ( const char *str );
/*---------------------------------------------------------------------------*/
/* Convert string into its double value.                                     */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_str_debug_string( int level, const char *key, const char *value );
/*---------------------------------------------------------------------------*/
/* print key & value info into LOG file.                                     */
/*---------------------------------------------------------------------------*/

static void _unur_str_debug_distr( int level, const char *name, double *params, int n_params );
/*---------------------------------------------------------------------------*/
/* write info about distribution into LOG file.                              */
/*---------------------------------------------------------------------------*/

static void _unur_str_debug_set( int level, const char *key, const char *type, ... );

/*---------------------------------------------------------------------------*/
/* write info about set command into LOG file.                               */
/*---------------------------------------------------------------------------*/

#endif
/*---------------------------------------------------------------------------*/


static void _unur_str_error_unknown( const char *file, int line, const char *key, const char *type );
/*---------------------------------------------------------------------------*/
/* print error message: unknown keyword.                                     */
/*---------------------------------------------------------------------------*/

static void _unur_str_error_invalid( const char *file, int line, const char *key, const char *type );
/*---------------------------------------------------------------------------*/
/* print error message: invalid data for keyword.                            */
/*---------------------------------------------------------------------------*/

static void _unur_str_error_args( const char *file, int line, const char *key );
/*---------------------------------------------------------------------------*/
/* print error message: invalid argument string for set call.                */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define _unur_error_unknown(key,what) \
   do { \
     _unur_str_error_unknown( __FILE__,__LINE__, (key), (what) ); \
   } while (0)

/* invalid data for key */
#define _unur_error_invalid(key,what) \
   do { \
     _unur_str_error_invalid( __FILE__,__LINE__, (key), (what) ); \
   } while (0)

/* invalid argument string for set calls */
#define _unur_error_args(key) \
   do { \
     _unur_str_error_args( __FILE__,__LINE__, (key) ); \
   } while (0)


/*---------------------------------------------------------------------------*/
/* constants                                                                 */

#define MAX_SET_ARGS  (10)  
/* maximal number of arguments for set calls (must be at least 2) */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Load lists of ..._new and ..._set calls                                **/
/*****************************************************************************/

#include "stringparser_lists.ch"

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_gen *
unur_str2gen (const char *string)
     /*----------------------------------------------------------------------*/
     /* get generator object for distribution and method                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   string ... string that contains description of generator           */
     /*                                                                      */
     /* return:                                                              */
     /*   generator object                                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_DISTR *distr = NULL;       /* distribution object */
  UNUR_PAR   *par   = NULL;       /* parameter object */
  UNUR_GEN   *gen   = NULL;       /* generator object    */
  UNUR_URNG  *urng  = NULL;       /* uniform random number generator */

  char *str_distr   = NULL;       /* string for distribution */
  char *str_method  = NULL;       /* string for method       */
  char *str_urng    = NULL;       /* string for uniform RNG  */

  char *str = NULL;               /* pointer to working string */
  char *token;

  struct unur_slist *mlist;       /* list of allocated memory blocks */

  /* check arguments */
  _unur_check_NULL( GENTYPE,string,NULL );

  /* (empty) list of allocated memory blocks */ 
  mlist = _unur_slist_new();

  /* Prepare string for processing:            */
  /*   Remove all white spaces,                */
  /*   convert to lower case letters,          */
  /*   copy into working string                */
  str = _unur_parser_prepare_string( string );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"[input]",str);
#endif

  /* split string into blocks separated by ampersands '&' */

  /* the first block must contain the distribution */
  str_distr = strtok(str, "&");

  /* new get all the other tokens */
  for ( token  = strtok(NULL, "&"); 
        token != NULL;     
        token  = strtok(NULL, "&") ) {
    
    /* get type of block */
    if (!strncmp(token,"method=",(size_t)7)) {
      str_method = token;
    }
    else if (!strncmp(token,"urng=",(size_t)5)) {
      str_urng = token;
    }
    else {
      _unur_error_unknown(token,"category");
      _unur_slist_free( mlist );
      if (str) free(str);
      return NULL;
    }
  }
  
  /* make distribution object */
  distr = _unur_str_distr(str_distr);
  if ( distr == NULL ) {
    /* error */
    _unur_slist_free( mlist );
    if (str) free(str);
    return NULL;
  }

  /* get parameter object */
  if ( str_method != NULL )
    /* method is provided in string */
    par = _unur_str_par(str_method, distr, mlist);
  else
    /* otherwise use default method */
    par = unur_auto_new(distr);

  /* make generator object */
  gen = unur_init(par);

  /* destroy distribution object */
  unur_distr_free(distr);

  /* set uniform random number generator -- if provided */
  if ( str_urng != NULL )
    if (gen != NULL) {
      /* we need a valid generator object */
      if ((urng = _unur_str2urng(str_urng)) != NULL )
	unur_chg_urng(gen, urng);
    }

  /* free allocated memory blocks */
  _unur_slist_free(mlist);
  if (str) free(str);

  /* check result: nothing to do */
  /*  if ( gen == NULL ) { */
  /*    error */
  /*    return NULL; */
  /*  } */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"",NULL);
#endif

  /* return pointer to generator object */
  return gen;

} /* end of unur_str2gen() */

/*---------------------------------------------------------------------------*/

struct unur_par *
_unur_str2par (const struct unur_distr *distr, const char *string, struct unur_slist **mlist )
     /*----------------------------------------------------------------------*/
     /* get parameter object for distribution and method                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   string ... string that contains description of method              */
     /*   mlist  ... pointer to simple list for allocated memory blocks      */
     /*                                                                      */
     /* return:                                                              */
     /*   parameter object                                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_PAR *par = NULL;           /* parameter object */
  char *str = NULL;               /* pointer to working string */

  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );
  _unur_check_NULL( GENTYPE,string,NULL );

  /* (empty) list of allocated memory blocks */ 
  *mlist = _unur_slist_new();

  /* Prepare string for processing:            */
  /*   Remove all white spaces,                */
  /*   convert to lower case letters,          */
  /*   copy into working string                */
  str = _unur_parser_prepare_string( string );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"[input]",str);
#endif

  /* get parameter object */
  par = _unur_str_par(str, distr, *mlist);

  /* free allocated memory blocks */
  if (str) free(str);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag)
    if ( par != NULL )
      _unur_str_debug_string(0,"",NULL);
#endif

  /* return pointer to parameter object */
  return par;

} /* end of _unur_str2par() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_str2distr (const char *string)
     /*----------------------------------------------------------------------*/
     /* get distribution object for given distribution                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   string ... string that contains description of distribution        */
     /*                                                                      */
     /* return:                                                              */
     /*   distribution object                                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_DISTR *distr = NULL;       /* distribution object */
  char *str = NULL;               /* pointer to working string */

  /* check arguments */
  _unur_check_NULL( GENTYPE,string,NULL );

  /* Prepare string for processing:            */
  /*   Remove all white spaces,                */
  /*   convert to lower case letters,          */
  /*   copy into working string                */
  str = _unur_parser_prepare_string( string );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"[input]",str);
#endif

  /* make distribution object */
  distr = _unur_str_distr(str);

  /* free allocated memory blocks */
  if (str) free(str);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag)
    if ( distr != NULL )
      _unur_str_debug_string(0,"",NULL);
#endif

  /* return pointer to distribution object */
  return distr;

} /* end of unur_str2distr() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
unur_makegen_ssu( const char *distrstr, const char *methodstr, UNUR_URNG *urng )
     /*----------------------------------------------------------------------*/
     /* get generator object for distribution and method                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distrstr  ... distribution described by a string                   */
     /*   methodstr ... chosen method described by a string                  */
     /*   urng      ... pointer to uniform random number generator           */
     /*                 (NULL = default generator)                           */ 
     /*                                                                      */
     /* return:                                                              */
     /*   generator object                                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_DISTR *distr = NULL;       /* distribution object */
  UNUR_PAR   *par   = NULL;       /* parameter object */
  UNUR_GEN   *gen   = NULL;       /* generator object    */

  char *str_distr   = NULL;       /* string for distribution */
  char *str_method  = NULL;       /* string for method       */

  struct unur_slist *mlist;       /* list of allocated memory blocks */

  /* check arguments */
  _unur_check_NULL( GENTYPE, distrstr, NULL );

  /* (empty) list of allocated memory blocks */
  mlist = _unur_slist_new();

  /* Prepare string for processing:            */
  /*   Remove all white spaces,                */
  /*   convert to lower case letters,          */
  /*   copy into working string                */
  str_distr = _unur_parser_prepare_string( distrstr );
  str_method = (methodstr) 
    ? _unur_parser_prepare_string( methodstr )
    : NULL;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag) {
    _unur_str_debug_string(0,"[input-distr]",str_distr);
    _unur_str_debug_string(0,"[input-method]",
			   (str_method) ? str_method : "(NULL)" );
  }
#endif

  do {
    /* make distribution object */
    distr = _unur_str_distr(str_distr);
    if (distr == NULL) break;

    /* get parameter object */
    if ( str_method != NULL && strlen(str_method)>0 )
      /* method is provided in string */
      par = _unur_str_par(str_method, distr, mlist);
    else
      /* otherwise use default method */
      par = unur_auto_new(distr);
    if (par == NULL) break;
    
    /* make generator object */
    gen = unur_init(par);
    if (gen == NULL) break;

    /* set uniform RNG */
    if (urng != NULL) 
      unur_chg_urng(gen, urng);

  } while (0);

  /* destroy distribution object */
  unur_distr_free(distr);

  /* free allocated memory blocks */
  _unur_slist_free(mlist);
  if (str_distr) free(str_distr);
  if (str_method) free(str_method);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"",NULL);
#endif

  /* return pointer to generator object */
  return gen;

} /* end of unur_makegen_ssu() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
unur_makegen_dsu( const struct unur_distr *distr, const char *methodstr, UNUR_URNG *urng )
     /*----------------------------------------------------------------------*/
     /* get generator object for distribution and method                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr     ... distribution object                                  */
     /*   methodstr ... chosen method described by a string                  */
     /*   urng      ... pointer to uniform random number generator           */
     /*                 (NULL = default generator)                           */ 
     /*                                                                      */
     /* return:                                                              */
     /*   generator object                                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_PAR   *par   = NULL;       /* parameter object */
  UNUR_GEN   *gen   = NULL;       /* generator object    */

  char *str_method  = NULL;       /* string for method       */

  struct unur_slist *mlist;       /* list of allocated memory blocks */

  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );

  /* (empty) list of allocated memory blocks */
  mlist = _unur_slist_new();

  /* Prepare string for processing:            */
  /*   Remove all white spaces,                */
  /*   convert to lower case letters,          */
  /*   copy into working string                */
  str_method = (methodstr)
    ? _unur_parser_prepare_string( methodstr )
    : NULL;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag) {
    _unur_str_debug_string(0,"[input-distr]","(distribution object)");
    _unur_str_debug_string(0,"[input-method]",
			   (str_method) ? str_method : "(NULL)" );
  }
#endif

  do {
    /* get parameter object */
    if ( str_method != NULL && strlen(str_method)>0 )
      /* method is provided in string */
      par = _unur_str_par(str_method, distr, mlist);
    else
      /* otherwise use default method */
      par = unur_auto_new(distr);
    if (par == NULL) break;
    
    /* make generator object */
    gen = unur_init(par);
    if (gen == NULL) break;

    /* set uniform RNG */
    if (urng != NULL)
      unur_chg_urng(gen, urng);

  } while (0);

  /* free allocated memory blocks */
  _unur_slist_free(mlist);
  if (str_method) free(str_method);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"",NULL);
#endif

  /* return pointer to generator object */
  return gen;

} /* end of unur_makegen_dsu() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Distributions                                                          **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_str_distr( char *str_distr )
     /*----------------------------------------------------------------------*/
     /* get distribution object for given distribution                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   str_distr ... string that contains description of distribution     */
     /*                                                                      */
     /* return:                                                              */
     /*   distribution object                                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *distr = NULL;

  char *token;             /* pointer to token in string */
  char *next;              /* pointer to next token in string */
  char *key, *value;       /* the key and its value */


  /* tokenize the string using ';' as separator and handle the tokens */
  for ( token = next = str_distr;
	next != NULL && *token != '\0';
	token = next ) {

    /* find end of token (i.e. the pointer to next token) */
    next = strchr(token,';');
    if (next != NULL) {
      /* next separator found */
      *next = '\0';     /* terminate token string    */
      next++;           /* set pointer to next token */
    }

    /* token points to a key */
    key = token;

    /* get the value from the pair key=value and terminate key string */
    value = strchr(key, '=');
    if (value != NULL) {
      *value = '\0';    /* terminate key string        */
      value++;          /* set pointer to value string */
    }

    /* the first key is treated differently: get new distribution object */
    if (key == str_distr) {
      /* the keyword "distr" is optional */
      if (value == NULL) {
	/* no "distribution" keyword given, so the key points to the value */
	value = key;
      }
      else {
	if (*key != 'd') {
	  _unur_error(GENTYPE,UNUR_ERR_STR_SYNTAX,"key for distribution does not start with 'd'"); 
	  _unur_distr_free(distr);  /* remark: added for ROOT to make coverity integrity manager happy */ 
	  return NULL;
	}
      }

      /* get new distribution object */
      if (distr != NULL) {
	_unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
	_unur_distr_free(distr);
      }

      distr = _unur_str_distr_new(value);
      if (distr == NULL) {
	/* error */
	/* _unur_error(GENTYPE,UNUR_ERR_STR,"setting distribution failed");  */
	return NULL;
      }
    }

    /* set parameters for distribution */
    else {
      if (_unur_str_distr_set(&distr, key, value)!=UNUR_SUCCESS ) {
	/* error */
	_unur_distr_free(distr);
	/* _unur_error(GENTYPE,UNUR_ERR_STR,"setting distribution failed");  */
	return NULL;
      }
    }
  }

  return distr;

} /* end of _unur_str_distr() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_str_distr_make_os (UNUR_DISTR *distr, const char *key, char *type_args, char **args)
     /*----------------------------------------------------------------------*/
     /* Make distribution object for order statistics of given distribution. */
     /* The given distribution `distr' is destroyed.                         */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "ii"   (int,int)                                               */
     /*                                                                      */
     /* input: "tt"  (exactly two single tokens)                             */
     /*        "L"   (exactly one list)                                      */
     /*                                                                      */
     /* action:                                                              */
     /*    if "tt": convert both to int.                                     */
     /*             execute set call.                                        */
     /*    if "L":  convert list entries to int.                             */
     /*             execute set call with first two entries.                 */
     /*                (error if list has less than two entries)             */
     /*    otherwise: error                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr     ... pointer to distribution object                       */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object for order statistics                */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  int *iarray = NULL;                 /* pointer to working array            */
  struct unur_distr *os = NULL;       /* pointer to order statistics object  */

  /* either exactly two arguments are given */
  if ( !strcmp(type_args, "tt") ) {
    iarray = _unur_xmalloc( 2*sizeof(double) );
    iarray[0] = _unur_atoi( args[0] );
    iarray[1] = _unur_atoi( args[1] );
  }

  /* or one list with 2 entries is given */
  else if ( !strcmp(type_args, "L") ) {
    if ( _unur_parse_ilist( args[0], &iarray ) < 2 ) {
      free (iarray);
      iarray = NULL;
    }
  }

  if (iarray == NULL ) {
    /* data list empty (incorrect syntax) */
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return NULL;
  }

  /* now make object to order statistics */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"ii",iarray[0],iarray[1] );
#endif

  os =  unur_distr_corder_new( distr, iarray[0], iarray[1] );

  /* free memory */
  _unur_distr_free(distr);
  free (iarray);
  
  /* return order statistics (or NULL we could not make such a thing) */
  return os;

} /* end of _unur_str_distr_make_os() */

/*---------------------------------------------------------------------------*/

int
_unur_str_distr_set_i (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_i set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given distribution object.                          */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "i"   (int)                                                    */
     /*                                                                      */
     /* input: "t"  (exactly one single token)                               */
     /*        void (no argument given)                                      */
     /*                                                                      */
     /* action:                                                              */
     /*    convert token to int.                                             */
     /*      (if no argument is given, use 1.)                               */
     /*    execute set call.                                                 */
     /*    error if input != "t"                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr     ... pointer to distribution object                       */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  int iarg;

  if ( !strcmp(type_args, "t") ) {
    iarg = _unur_atoi( args[0] );
  }
  else if ( strlen(type_args) == 0 ) {
    iarg = 1;
  }
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
  
#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"i",iarg);
#endif

  return set(distr,iarg);

} /* end of _unur_str_distr_set_i() */

/*---------------------------------------------------------------------------*/

int
_unur_str_distr_set_ii (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_ii set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given distribution object.                          */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "ii"   (int,int)                                               */
     /*                                                                      */
     /* input: "tt"  (exactly two single tokens)                             */
     /*        "L"   (exactly one list)                                      */
     /*                                                                      */
     /* action:                                                              */
     /*    if "tt": convert both to int.                                     */
     /*             execute set call.                                        */
     /*    if "L":  convert list entries to int.                             */
     /*             execute set call with first two entries.                 */
     /*                (error if list has less than two entries)             */
     /*    otherwise: error                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr     ... pointer to distribution object                       */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  int *iarray = NULL;
  int iarg[2];
  int result;

  /* either exactly two arguments are given */
  if ( !strcmp(type_args, "tt") ) {
    iarg[0] = _unur_atoi( args[0] );
    iarg[1] = _unur_atoi( args[1] );

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"ii",iarg[0],iarg[1]);
#endif

    return set( distr,iarg[0],iarg[1] );
  }

  /* or one list with 2 entries is given */
  else if ( !strcmp(type_args, "L") ) {
    if ( _unur_parse_ilist( args[0], &iarray ) < 2 ) {
      _unur_error_args(key);
      free (iarray);
#ifdef UNUR_ENABLE_LOGGING
      /* write info into LOG file */
      if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	_unur_str_debug_set(2,key,"!");
#endif
      return UNUR_ERR_STR_INVALID;
    }
    result = set( distr,iarray[0],iarray[1] );

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"ii",iarray[0],iarray[1] );
#endif

    free (iarray);
    return result;
  }
  
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }

} /* end of _unur_str_distr_set_ii() */

/*---------------------------------------------------------------------------*/

int
_unur_str_distr_set_d (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_d set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given distribution object.                          */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "d"   (double)                                                 */
     /*                                                                      */
     /* input: "t"  (exactly one single token)                               */
     /*                                                                      */
     /* action:                                                              */
     /*    convert token to double.                                          */
     /*    execute set call.                                                 */
     /*    error if input != "t"                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr     ... pointer to distribution object                       */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  double darg;

  if ( strcmp(type_args, "t") ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
  
  darg = _unur_atod( args[0] );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"d",darg);
#endif

  return set(distr,darg);

} /* end of _unur_str_distr_set_d() */

/*---------------------------------------------------------------------------*/

int
_unur_str_distr_set_dd (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_dd set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given distribution object.                          */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "dd"   (double,double)                                         */
     /*                                                                      */
     /* input: "tt"  (exactly two single tokens)                             */
     /*        "L"   (exactly one list)                                      */
     /*                                                                      */
     /* action:                                                              */
     /*    if "tt": convert both to double.                                  */
     /*             execute set call.                                        */
     /*    if "L":  convert list entries to double.                          */
     /*             execute set call with first two entries.                 */
     /*                (error if list has less than two entries)             */
     /*    otherwise: error                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr     ... pointer to distribution object                       */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  double *darray = NULL;
  double darg[2];
  int result;

  /* either exactly two arguments are given */
  if ( !strcmp(type_args, "tt") ) {
    darg[0] = _unur_atod( args[0] );
    darg[1] = _unur_atod( args[1] );

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"dd",darg[0],darg[1]);
#endif

    return set( distr,darg[0],darg[1] );
  }

  /* or one list with 2 entries is given */
  else if ( !strcmp(type_args, "L") ) {
    if ( _unur_parse_dlist( args[0], &darray ) < 2 ) {
      _unur_error_args(key);
      free (darray);
#ifdef UNUR_ENABLE_LOGGING
      /* write info into LOG file */
      if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	_unur_str_debug_set(2,key,"!");
#endif
      return UNUR_ERR_STR_INVALID;
    }
    result = set( distr,darray[0],darray[1] );

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"dd",darray[0],darray[1] );
#endif

    free (darray);
    return result;
  }
  
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }

} /* end of _unur_str_distr_set_dd() */

/*---------------------------------------------------------------------------*/

int
_unur_str_distr_set_Di (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_Di set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given distribution object.                          */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "Di"   (array of doubles, int)                                 */
     /*                                                                      */
     /* input: "Lt"  (list & single token)                                   */
     /*        "L"   (exactly one list)                                      */
     /*                                                                      */
     /* action:                                                              */
     /*    if "L":  convert list to array of doubles.                        */
     /*             execute set with size if array as second argument.       */
     /*    if "Lt": convert list to array of doubles.                        */
     /*             convert second token to int.                             */
     /*             set int to min of int and size of array.                 */
     /*             execute set call with int and array.                     */
     /*    otherwise: error                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr     ... pointer to distribution object                       */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  int result;
  int t_size;
  int size = -1;
  double *darray = NULL;

  /* either an integer and a list are given */
  if ( !strcmp(type_args, "Lt") ) {
    t_size = _unur_atoi( args[1] );
    size = _unur_parse_dlist( args[0], &darray );
    if (size > t_size)
      size = t_size;
  }
    
  /* or only a list is given */
  else if ( !strcmp(type_args, "L") ) {
    size = _unur_parse_dlist( args[0], &darray );
  }

  /* or there is a syntax error */
  if ( !(size>0) ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    result = UNUR_ERR_STR_INVALID;
  }

  else {
    /* execute set command */
    result = set( distr,darray,size );
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"Di",darray,size,size);
#endif
  }

  if (darray) free (darray);

  return result;

} /* end of _unur_str_distr_set_Di() */

/*---------------------------------------------------------------------------*/

int
_unur_str_distr_set_C (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_C set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given distribution object.                          */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "C"   (string = array of char)                                 */
     /*                                                                      */
     /* input: "s"  (exactly one string)                                     */
     /*                                                                      */
     /* action:                                                              */
     /*    execute set call.                                                 */
     /*    error if input != "s"                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr     ... pointer to distribution object                       */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   's'   for single string                            */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  char *string;

  if ( strcmp(type_args, "s") ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }

  /* get string */
  string = args[0];
  
#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"C",string);
#endif

  return set(distr,string);

} /* end of _unur_str_distr_set_C() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Methods                                                                **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_par *
_unur_str_par( char *str_method, const UNUR_DISTR *distr, struct unur_slist *mlist )
     /*----------------------------------------------------------------------*/
     /* get parameter object for distribution and method                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   str_method ... string that contains method description             */
     /*   distr      ... pointer to distribution object                      */
     /*   mlist      ... list of allocated memory blocks                     */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par = NULL;     /* pointer to parameter object */

  char *token;             /* pointer to token in string */
  char *next;              /* pointer to next token in string */
  char *key, *value;       /* the key and its value */

  /* tokenize the string using ';' as separator and handle the tokens */
  for ( token = next = str_method;
	next != NULL && *token != '\0';
	token = next ) {

    /* find end of token (i.e. the pointer to next token) */
    next = strchr(token,';');
    if (next != NULL) {
      /* next separator found */
      *next = '\0';     /* terminate token string    */
      next++;           /* set pointer to next token */
    }

    /* token points to a key */
    key = token;

    /* get the value from the pair key=value and terminate key string */
    value = strchr(key, '=');
    if (value != NULL) {
      *value = '\0';    /* terminate key string        */
      value++;          /* set pointer to value string */
    }

    /* the first key is treated differently: get new parameter object */
    if (key == str_method) {
      /* the keyword "method" is optional */
      if (value == NULL) {
	/* no "method" keyword given, so the key points to the value */
	value = key;
      }
      else {
	if (*key != 'm') {
	  _unur_error(GENTYPE,UNUR_ERR_STR_SYNTAX,"key for method does not start with 'm'"); 
	  return NULL;
	}
      }
      
      /* get new parameter object */
      par = _unur_str_par_new(value,distr);
      if (par == NULL) {
	/* error */
	_unur_error(GENTYPE,UNUR_ERR_STR,"setting method failed"); 
	return NULL;
      }
    }

    /* set parameters for method */
    else {
      if ( !_unur_str_par_set(par, key, value, mlist) ) {
	; /* error: */
      }
    }
  }

  /* end */
  return par;

} /* end of _unur_str_par() */

/*---------------------------------------------------------------------------*/

int
_unur_str_par_set_void (UNUR_PAR *par, const char *key, 
			char *type_args, char **args ATTRIBUTE__UNUSED, par_set_void set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given parameter object.                             */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: void   (no arguments required)                                 */
     /*                                                                      */
     /* input: void  (no arguments should be given)                          */
     /*                                                                      */
     /* action:                                                              */
     /*    execute set call.                                                 */
     /*    warning if any argument is given.                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter object                          */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"v");
#endif

  if (*type_args != '\0')
    _unur_error_args(key);

  return set(par);

} /* end of _unur_str_par_set_void() */

/*---------------------------------------------------------------------------*/

int
_unur_str_par_set_i (UNUR_PAR *par, const char *key, char *type_args, char **args, par_set_i set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given parameter object.                             */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "i"   (int)                                                    */
     /*                                                                      */
     /* input: "t"   (exactly one single token)                              */
     /*         void (no argument given)                                     */
     /*                                                                      */
     /* action:                                                              */
     /*    convert token to int.                                             */
     /*      (if no argument is given, use 1.)                               */
     /*    execute set call.                                                 */
     /*    error if input != "t"                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter object                          */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  int iarg;

  if ( !strcmp(type_args, "t") ) {
    iarg = _unur_atoi( args[0] );
  }
  else if ( strlen(type_args) == 0 ) {
    iarg = 1;
  }
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
  
#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"i",iarg);
#endif

  return set(par,iarg);

} /* end of _unur_str_par_set_i() */

/*---------------------------------------------------------------------------*/

int
_unur_str_par_set_ii (UNUR_PAR *par, const char *key, char *type_args, char **args, par_set_ii set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given parameter object.                             */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "ii"   (int,int)                                               */
     /*                                                                      */
     /* input: "tt"  (exactly two single tokens)                             */
     /*        "L"   (exactly one list)                                      */
     /*                                                                      */
     /* action:                                                              */
     /*    if "tt": convert both to int.                                     */
     /*             execute set call.                                        */
     /*    if "L":  convert list entries to int.                             */
     /*             execute set call with first two entries.                 */
     /*                (error if list has less than two entries)             */
     /*    otherwise: error                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter object                          */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  int *iarray = NULL;
  int iarg[2];
  int result;

  /* either exactly two arguments are given */
  if ( !strcmp(type_args, "tt") ) {
    iarg[0] = _unur_atoi( args[0] );
    iarg[1] = _unur_atoi( args[1] );

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"ii",iarg[0],iarg[1]);
#endif

    return set( par,iarg[0],iarg[1] );
  }

  /* or one list with 2 entries is given */
  else if ( !strcmp(type_args, "L") ) {
    if ( _unur_parse_ilist( args[0], &iarray ) < 2 ) {
      _unur_error_args(key);
      free (iarray);
#ifdef UNUR_ENABLE_LOGGING
      /* write info into LOG file */
      if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	_unur_str_debug_set(2,key,"!");
#endif
      return UNUR_ERR_STR_INVALID;
    }
    result = set( par,iarray[0],iarray[1] );

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"ii",iarray[0],iarray[1] );
#endif

    free (iarray);
    return result;
  }
  
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }

} /* end of _unur_str_par_set_ii() */

/*---------------------------------------------------------------------------*/

int
_unur_str_par_set_u (UNUR_PAR *par, const char *key, char *type_args, char **args, par_set_u set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given parameter object.                             */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "u"   (unsigned)                                               */
     /*                                                                      */
     /* input: "t"  (exactly one single token)                               */
     /*                                                                      */
     /* action:                                                              */
     /*    convert token to unsigned.                                        */
     /*    execute set call.                                                 */
     /*    error if input != "t"                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter object                          */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  unsigned uarg;

  if ( strcmp(type_args, "t") ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
  
  uarg = _unur_atou( args[0] );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"u",uarg);
#endif

  return set(par,uarg);

} /* end of _unur_str_par_set_u() */

/*---------------------------------------------------------------------------*/

int
_unur_str_par_set_d (UNUR_PAR *par, const char *key, char *type_args, char **args, par_set_d set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given parameter object.                             */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "d"   (double)                                                 */
     /*                                                                      */
     /* input: "t"  (exactly one single token)                               */
     /*                                                                      */
     /* action:                                                              */
     /*    convert token to double.                                          */
     /*    execute set call.                                                 */
     /*    error if input != "t"                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter object                          */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  double darg;

  if ( strcmp(type_args, "t") ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
  
  darg = _unur_atod( args[0] );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"d",darg);
#endif

  return set(par,darg);

} /* end of _unur_str_par_set_d() */

/*---------------------------------------------------------------------------*/

int
_unur_str_par_set_dd (UNUR_PAR *par, const char *key, char *type_args, char **args, par_set_dd set)
     /*----------------------------------------------------------------------*/
     /* Set parameter in given parameter object.                             */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "dd"   (double,double)                                         */
     /*                                                                      */
     /* input: "tt"  (exactly two single tokens)                             */
     /*        "L"   (exactly one list)                                      */
     /*                                                                      */
     /* action:                                                              */
     /*    if "tt": convert both to double.                                  */
     /*             execute set call.                                        */
     /*    if "L":  convert list entries to double.                          */
     /*             execute set call with first two entries.                 */
     /*                (error if list has less than two entries)             */
     /*    otherwise: error                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter object                          */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  double *darray = NULL;
  double darg[2];
  int result;

  /* either exactly two arguments are given */
  if ( !strcmp(type_args, "tt") ) {
    darg[0] = _unur_atod( args[0] );
    darg[1] = _unur_atod( args[1] );

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"dd",darg[0],darg[1]);
#endif

    return set( par,darg[0],darg[1] );
  }

  /* or one list with 2 entries is given */
  else if ( !strcmp(type_args, "L") ) {
    if ( _unur_parse_dlist( args[0], &darray ) < 2 ) {
      _unur_error_args(key);
      free (darray);
#ifdef UNUR_ENABLE_LOGGING
      /* write info into LOG file */
      if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	_unur_str_debug_set(2,key,"!");
#endif
      return UNUR_ERR_STR_INVALID;
    }
    result = set( par,darray[0],darray[1] );

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"dd",darray[0],darray[1] );
#endif

    free (darray);
    return result;
  }
  
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }

} /* end of _unur_str_par_set_dd() */

/*---------------------------------------------------------------------------*/

int
_unur_str_par_set_iD (UNUR_PAR *par, const char *key, char *type_args, char **args,
		      par_set_iD set, struct unur_slist *mlist )
     /*----------------------------------------------------------------------*/
     /* Set parameter in given parameter object.                             */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "iD"   (int, array of doubles)                                 */
     /*                                                                      */
     /* input: "tL"  (single token & list)                                   */
     /*        "L"   (exactly one list)                                      */
     /*        "t"   (single token)                                          */
     /*                                                                      */
     /* action:                                                              */
     /*    if "L":  convert list to array of doubles.                        */
     /*             execute set with size if array as first argument.        */
     /*    if "t":  convert token to int.                                    */ 
     /*             execute set call with NULL as its second argument.       */
     /*    if "tL": convert list to array of doubles.                        */
     /*             convert first token to int.                              */
     /*             set int to min of int and size of array.                 */
     /*             execute set call with int and array.                     */
     /*    otherwise: error                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter object                          */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*   mlist     ... list of allocated memory blocks                      */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  int result;
  int t_size;
  int size = -1;
  double *darray = NULL;

  /* either an integer and a list are given */
  if ( !strcmp(type_args, "tL") ) {
    t_size = _unur_atoi( args[0] );
    size = _unur_parse_dlist( args[1], &darray );
    if (size > 0) {
      if (size > t_size)  size = t_size;
    }
    else {
      /* empty list */
      size = t_size;
      if (darray) free (darray);
      darray = NULL;
    }
  }

  /* or only an integer is given */
  else if ( !strcmp(type_args, "t") ) {
    size = _unur_atoi( args[0] );
    darray = NULL;
  }
  
  /* or only a list is given */
  else if ( !strcmp(type_args, "L") ) {
    size = _unur_parse_dlist( args[0], &darray );
  }
  
  /* or there is a syntax error */
  if ( !(size>0) ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    result = UNUR_ERR_STR_INVALID;
  }
  
  else {
    /* execute set command */
    result = set( par,size,darray );
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"iD",size,darray,size);
#endif
  }

  /* store pointer to allocated array */
  if (darray) 
    _unur_slist_append( mlist, darray );

  return result;

} /* end of _unur_str_par_set_iD() */

/*---------------------------------------------------------------------------*/

int
_unur_str_par_set_Di (UNUR_PAR *par, const char *key, char *type_args, char **args,
		      par_set_Di set, struct unur_slist *mlist )
     /*----------------------------------------------------------------------*/
     /* Set parameter in given parameter object.                             */
     /* The list of arguments and their types are given as string.           */
     /*                                                                      */
     /* type: "Di"   (array of doubles, int)                                 */
     /*                                                                      */
     /* input: "Lt"  (list & single token & list)                            */
     /*                                                                      */
     /* action:                                                              */
     /*    convert list to array of doubles.                                 */
     /*    convert token to int.                                             */
     /*    execute set call with array and int.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter object                          */
     /*   key       ... name of parameter to be set                          */
     /*   type_args ... contains types of arguments. Possible values:        */
     /*                   't'   for single token                             */
     /*                   'L'   for a list                                   */
     /*   args      ... array of pointers to argument strings                */
     /*   set       ... pointer to UNU.RAN function for set call             */
     /*   mlist     ... list of allocated memory blocks                      */
     /*                                                                      */
     /* return:                                                              */
     /*   return code from set call, i.e.                                    */
     /*      UNUR_SUCCESS ... on success                                     */
     /*      error code   ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  int result;
  int t_size;
  int size;
  double *darray = NULL;

  if ( !strcmp(type_args, "Lt") ) {
    t_size = _unur_atoi( args[1] );
    size = _unur_parse_dlist( args[0], &darray );
    if (size > 0) {
      result = set( par,darray,t_size );
#ifdef UNUR_ENABLE_LOGGING
      /* write info into LOG file */
      if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	_unur_str_debug_set(2,key,"Di",darray,size,size);
#endif

      /* store pointer to allocated array */
      if (darray) 
	_unur_slist_append( mlist, darray );

      return result;
    }
  }

  /* there is a syntax error */
  _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"!");
#endif

  return UNUR_ERR_STR_INVALID;
} /* end of _unur_str_par_set_Di() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Uniform RNG                                                            **/
/*****************************************************************************/

UNUR_URNG *
_unur_str2urng( char *str_urng ATTRIBUTE__UNUSED)
     /*----------------------------------------------------------------------*/
     /* get uniform RNG                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   str_urng ... string that contains description of uniform RNG       */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to uniform RNG                                             */
     /*                                                                      */
     /* The string must not contain any other keywords than "urng".          */
     /* the argument is passed unchanged to the PRNG package.                */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
/*---------------------------------------------------------------------------*/
#if defined(UNUR_URNG_UNURAN) && defined(UNURAN_HAS_PRNG)
/*---------------------------------------------------------------------------*/

  UNUR_URNG *urng = NULL;

  char *token;             /* pointer to token in string */
  char *next;              /* pointer to next token in string */
  char *key, *value;       /* the key and its value */

  /* tokenize the string using ';' as separator and handle the tokens */
  for ( token = next = str_urng;
	next != NULL && *token != '\0';
	token = next ) {

    /* find end of token (i.e. the pointer to next token) */
    next = strchr(token,';');
    if (next != NULL) {
      /* next separator found */
      *next = '\0';     /* terminate token string    */
      next++;           /* set pointer to next token */
    }
    
    /* token points to a key */
    key = token;

    /* get the value from the pair key=value and terminate key string */
    value = strchr(key, '=');
    if (value != NULL) {
      *value = '\0';    /* terminate key string        */
      value++;          /* set pointer to value string */
    }
    
    /* get (default) parameter object for method */
    if ( !strcmp( key, "urng") ) {
      if (key == str_urng) {
	/* "urng" must be the first key.                                  */
	/* (on the other hand, the first keyword always must be "urng"    */
	/* since otherwise we had not called this subroutine.             */
	/* hence, we only check, whether the keyword "urng" is unique.)   */
	urng = unur_urng_prng_new(value);
	if (urng == NULL) {
	  /* error */
	  _unur_error(GENTYPE,UNUR_ERR_STR,"setting URNG failed"); 
	  return NULL;
	}
	else {
#ifdef UNUR_ENABLE_LOGGING
	  /* write info into LOG file */
	  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	    _unur_str_debug_string(1,"urng",value);
#endif
	  ;
	}
      }
      else {
	if (urng) unur_urng_free (urng);
	_unur_error(GENTYPE,UNUR_ERR_STR_SYNTAX,"urng must be first key"); 
	return NULL;
      }
    }

    /* there are no other parameters */
    else {
      /* error: unknown parameter */
      if (urng) unur_urng_free (urng);
      _unur_error_unknown(key,"parameter for uniform RNG");
      return NULL;
    }
  }

  /*return pointer to uniform RNG */
  return urng;

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/

  _unur_error(GENTYPE,UNUR_ERR_STR,"setting URNG requires PRNG library enabled"); 
  return NULL;

/*---------------------------------------------------------------------------*/
#endif   /*  end defined(UNUR_URNG_UNURAN) && defined(UNURAN_HAS_PRNG) */
/*---------------------------------------------------------------------------*/
} /* end of _unur_str2urng() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Misc                                                                   **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_str_set_args( char *value, char *type_args, char **args, int max_args )
     /*----------------------------------------------------------------------*/
     /* parse argument string for set call                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   value     ... string that contains list of arguments for set call  */
     /*   type_args ... array containing info about arguments                */
     /*   args      ... pointers to argument strings                         */
     /*   max_args  ... maximal number of argments                           */
     /*                                                                      */
     /* tokenize string with arguments and store the pointer to each         */
     /* argument in 'args'. Evaluate the type of each argument and store it  */
     /* in 'type_args'. We have the following types:                         */
     /*    't' ... a single token, or                                        */
     /*    'L' ... a list, indicated by (...)                                */
     /*                                                                      */
     /* return:                                                              */
     /*   number of parsed arguments                                         */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                           */
     /*----------------------------------------------------------------------*/
{
  char *token;             /* pointer to token in string */
  char *next;              /* pointer to next token in string */
  int n_args;              /* number of arguments in set call */

  /* initialize pointers */
  n_args = 0;
  *type_args = '\0';
  *args = NULL;

  /* we have two types of arguments:
     't' ... a single token, or
     'L' ... a list, indicated by (...)
  */

  /* split value string into separate arguments */
  if (value && *value != '\0') {
    /* there is at least one argument */

    for ( token = next = value;
	  next != NULL && *token != '\0' && n_args < max_args;
	  token = next, n_args++ ) {
      
      /* get type of argument */

      if ( *token == '(' ) {
	/* argument is a list */
	type_args[n_args] = 'L';
	/* skip to next character */
	token++;
	/* store pointer to argument */
	args[n_args] = token;
	/* find end of list using ')' as terminating character */
	next = strchr(token,')');
	if (next != NULL) {
	  /* end of list found */
	  *next = '\0';      /* terminate token string    */
	  token = ++next;    /* set pointer to first character after ')' */
	  /* first character after ')' must be ',' */
	  if ( !(*token == ',' || *token == '\0') ) {
	    _unur_error(GENTYPE,UNUR_ERR_STR_SYNTAX,") not followed by comma"); 
	    return -1;
	  }
	}
	else {
	  /* no ')' found */
	  token = NULL;
	}
      }

      else if ( *token == '"' ) {
	/* argument is a string */
	type_args[n_args] = 's';
	/* skip to next character */
	token++;
	/* store pointer to argument */
	args[n_args] = token;

	/* find end of string using '"' as terminating character */
	next = strchr(token,'"');
	if (next != NULL) {
	  /* end of list found */
	  *next = '\0';      /* terminate token string    */
	  token = ++next;    /* set pointer to first character after ')' */
	  /* first character after ')' must be ',' */
	  if ( !(*token == ',' || *token == '\0') ) {
	    _unur_error(GENTYPE,UNUR_ERR_STR_SYNTAX,"closing '\"' not followed by comma"); 
	    return -1;
	  }
	}
	else {
	  /* no '"' found */
	  token = NULL;
	}
      }

      else {
	/* argument is a sigle item */
	type_args[n_args] = 't';
	/* store pointer to argument */
	args[n_args] = token;
      }

      /* find end of argument string using ',' as separator */
      if (token) {
	next = strchr(token,',');
	if (next != NULL) {
	  /* next separator found */
	  *next = '\0';     /* terminate token string    */
	  next++;           /* set pointer to next token */
	}
      }
      else {
	next = NULL;
      }
    }

    /* terminate string */
    type_args[n_args] = '\0';

    /* If there are two many arguments, we simply assume a syntax error */
    /* (This case should not happen.)                                   */
    if (n_args >= max_args) { 
      _unur_error( GENTYPE, UNUR_ERR_STR_SYNTAX, "too many arguments");
      return -1;
    }
  }

  /* return number of parsed arguments */
  return n_args;

} /* end of _unur_str_set_args() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Convert to numbers                                                     **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_parse_ilist( char *liststr, int **iarray )
     /*----------------------------------------------------------------------*/
     /* Process a comma separated list of integers, convert it to ints and   */
     /* store it an array. The list (string) is terminated by a closing      */
     /* parenthesis ')' or '\0'.                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   liststr ... string containing list of numbers                      */
     /*   iarray  ... pointer to adress of double array for storing result   */
     /*                                                                      */
     /* return:                                                              */
     /*   number of elements extracted                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   as a side effect an int array of appropriate size is allocated.    */
     /*----------------------------------------------------------------------*/
{
  int *iarr = NULL;     /* pointer to int array */
  int n_iarray = 0;     /* number of elements in processed list */
  int n_alloc = 0;      /* number allocated elements */
  char *token;

  /* return 0 if no list is given */
  if (liststr == NULL) {
    *iarray = iarr;
    return 0;
  }

  /* discard possible leading ',' or '(' characters in front of the numbers */
  while ( *liststr == ',' || *liststr == '(' ){
    liststr++; /* next char */
  }

  /* new get integers */
  for ( token  = strtok(liststr, ",)"); 
        token != NULL;     
        token  = strtok(NULL, ",)") ) {

    /* check if there is enough space in iarr */
    if (n_iarray >= n_alloc) {
      /* no --> allocate memory */
      n_alloc += 100;
      iarr = _unur_xrealloc( iarr, n_alloc * sizeof(int) );
    }
    
    /* get int */
    iarr[n_iarray++] = _unur_atoi(token);

  } /* end while -- all elements of list read */

  /* store pointer to double array */
  *iarray = iarr;

  /* return number of elements in array */
  return n_iarray;

} /* end of _unur_parse_ilist() */

/*---------------------------------------------------------------------------*/

int 
_unur_parse_dlist( char *liststr, double **darray )
     /*----------------------------------------------------------------------*/
     /* Process a comma separated list of numbers (doubles), convert it to   */
     /* doubles and store it an array. The list (string) is terminated by    */
     /* a closing parenthesis ')' or '\0'.                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   liststr ... string containing list of numbers                      */
     /*   darray  ... pointer to adress of double array for storing result   */
     /*                                                                      */
     /* return:                                                              */
     /*   number of elements extracted                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   as a side effect a double array of appropriate size is allocated.  */
     /*----------------------------------------------------------------------*/
{
  double *darr = NULL;  /* pointer to double array */
  int n_darray = 0;     /* number of elements in processed list */
  int n_alloc = 0;      /* number allocated elements */
  char *token;          /* pointer to token in string */
  char *next;           /* pointer to next token in string */

  /* return 0 if no list is given */
  if (liststr == NULL) {
    *darray = NULL;
    return 0;
  }

  /* extract doubles from string and write them to array of doubles.  */
  /* end of list is indicated by right bracket ')' or '\0'.           */  

  /* pointer to first token in string */
  token = liststr;

  /* remove (optional) leading parenthesis '('                        */
  while (*token != '\0' && *token == '(')
    ++token;

  /* tokenize the string using ',' as separator and handle the tokens */
  for ( next = token;
	next != NULL && *token != '\0' &&*token != ')';
	token = next ) {

    /* find end of token (i.e. the pointer to next token) */
    next = strchr(token,',');
    if (next != NULL) {
      /* next separator found */
      *next = '\0';     /* terminate token string    */
      next++;           /* set pointer to next token */
    }

   /* check if there is enough space in darr */
    if (n_darray >= n_alloc) {
      /* no --> allocate memory */
      n_alloc += 100;
      darr = _unur_xrealloc( darr, n_alloc * sizeof(double) );
    }

    /* extract double and write to array and increase counter */
    darr[n_darray] = _unur_atod(token);
    n_darray++;

  }

  /* store pointer to double array */
  *darray = darr;

  /* return number of elements in array */
  return n_darray;

} /* end of _unur_parse_dlist() */

/*---------------------------------------------------------------------------*/

int 
_unur_atoi ( const char *str )
     /*----------------------------------------------------------------------*/
     /* Converts the string into its integer representation.                 */
     /* 'true' and 'on' are converted to 1, 'false' and 'off' to 0.          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   str      ... pointer to string                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   integer                                                            */
     /*----------------------------------------------------------------------*/
{
  if ( !strcmp(str,"true") || !strcmp(str,"on") )
    return 1;

  else if ( !strcmp(str,"false") || !strcmp(str,"off") )
    return 0;

  else if ( !strncmp(str,"inf",(size_t)3) )
    return INT_MAX;

  else if ( !strncmp(str,"-inf",(size_t)4) )
    return INT_MIN;

  else
    return atoi(str);

} /* end of _unur_atoi() */

/*---------------------------------------------------------------------------*/

unsigned
_unur_atou ( const char *str )
     /*----------------------------------------------------------------------*/
     /* Converts the string into its unsigned representation.                */
     /* strings are interprated as hexadecimal numbers.                      */
     /* 'true' and 'on' are converted to 1u, 'false' and 'off' to 0u.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   str      ... pointer to string                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   unsigned                                                           */
     /*----------------------------------------------------------------------*/
{
  char *tail;

  if ( !strcmp(str,"true") || !strcmp(str,"on") )
    return 1u;

  else if ( !strcmp(str,"false") || !strcmp(str,"off") )
    return 0u;

  else
    return ((unsigned) strtoul(str, &tail, 16));
} /* end of _unur_atou() */

/*---------------------------------------------------------------------------*/

double
_unur_atod ( const char *str )
     /*----------------------------------------------------------------------*/
     /* Converts the string into its double value.                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   str      ... pointer to string                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   double                                                             */
     /*----------------------------------------------------------------------*/
{
  if ( !strncmp(str,"inf",(size_t)3) )
    return INFINITY;

  else if ( !strncmp(str,"-inf",(size_t)4) )
    return -INFINITY;

  else
    return atof(str);
} /* end of _unur_atod() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void 
_unur_str_debug_string( int level, const char *key, const char *value  )
     /*----------------------------------------------------------------------*/
     /* print key & value info into LOG file                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   level ... level of indentation                                     */
     /*   key   ... pointer to key string                                    */
     /*   value ... pointer to value string                                  */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  
  LOG = unur_get_stream();

/*    fprintf(LOG,"%s: String Interface for UNU.RAN\n",GENTYPE); */
  fprintf(LOG,"%s: ",GENTYPE);
  for (; level>0; level--) 
    fprintf(LOG,"\t");
  fprintf(LOG,"%s",key);
  if (value)
    fprintf(LOG,": %s",value);
  fprintf(LOG,"\n");

/*    fprintf(LOG,"%s:\n",GENTYPE); */

} /* end of _unur_str_debug_string() */

/*---------------------------------------------------------------------------*/

void 
_unur_str_debug_distr( int level, const char *name, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into LOG file                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   level    ... level of indentation                                  */
     /*   name     ... name of distribution                                  */
     /*   params   ... array of parameters                                   */
     /*   n_params ... number of parameters                                  */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;
  
  LOG = unur_get_stream();

  fprintf(LOG,"%s: ",GENTYPE);
  for (; level>0; level--) 
    fprintf(LOG,"\t");

  fprintf(LOG,"distribution = %s (",name);
  if (n_params >0) {
    fprintf(LOG,"%g",params[0]);
    for (i=1; i<n_params; i++)
      fprintf(LOG,", %g",params[i]);
  }
  fprintf(LOG,")\n");

} /* end of _unur_str_debug_distr() */

/*---------------------------------------------------------------------------*/

void
_unur_str_debug_set( int level, const char *key, const char *type, ... )
     /*----------------------------------------------------------------------*/
     /* write info about set command into LOG file                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   level ... level of indentation                                     */
     /*   key   ... pointer to key string                                    */
     /*   type  ... list of types of argments                                */
     /*   ...   ... variable list of arguments                               */
     /*----------------------------------------------------------------------*/
{
  va_list ap;
  FILE *LOG;
  
  va_start(ap, type);

  LOG = unur_get_stream();

  /* print key name */
  fprintf(LOG,"%s: ",GENTYPE);
  for (; level>0; level--) 
    fprintf(LOG,"\t");
  fprintf(LOG,"%s: ",key);

  while (1) {

    switch (*type) {
    case 'v':
      fprintf(LOG," -none-");
      break;
    case 'd':
      fprintf(LOG," %g",va_arg(ap,double));
      break;
    case 'i':
      fprintf(LOG," %d",va_arg(ap,int));
      break;
    case 'u':
      fprintf(LOG," %x",va_arg(ap,unsigned int));
      break;
    case 'C':
      fprintf(LOG," %s",va_arg(ap,char *));
      break;
    case 'D': {
      int i,size;
      double *darray;
      darray = va_arg(ap,double *);
      size = va_arg(ap,int);
      if (size > 0 && darray != NULL) {
	fprintf(LOG," (%g",darray[0]);
	for (i=1; i<size; i++)
	  fprintf(LOG,",%g",darray[i]);
	fprintf(LOG,")");
      }
      else
	fprintf(LOG," (empty)");
      break;
    }
    case '!':
    default:
      fprintf(LOG," syntax error");
      break;
    }

    if ( *(++type) != '\0' )
      /* skip. there is a next argument to be processed */
      fprintf(LOG,",");
    else
      /* done */
      break;
  }

  /* dd, iD, Di */

  fprintf(LOG,"\n");
  fflush(LOG);   /* in case of a segmentation fault */

  va_end(ap);

} /* end of _unur_str_debug_set() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

void
_unur_str_error_unknown( const char *file, int line, const char *key, const char *type )
     /*----------------------------------------------------------------------*/
     /* print error message: unknown keyword.                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   file ... file name (inserted by __FILE__)                          */
     /*   line ... line number in source file (inserted by __LINE__)         */
     /*   key  ... keyword                                                   */
     /*   type ... type of keyword                                           */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *reason = _unur_string_new();
  _unur_string_append( reason, "unknown %s: '%s'", type, key );
  _unur_error_x( GENTYPE, file, line, "error", UNUR_ERR_STR_UNKNOWN,reason->text);
  _unur_string_free( reason );
} /* end of _unur_str_error_unknown() */

/*---------------------------------------------------------------------------*/

void
_unur_str_error_invalid( const char *file, int line, const char *key, const char *type )
     /*----------------------------------------------------------------------*/
     /* print error message: invalid data for keyword.                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   file ... file name (inserted by __FILE__)                          */
     /*   line ... line number in source file (inserted by __LINE__)         */
     /*   key  ... keyword                                                   */
     /*   type ... type of keyword                                           */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *reason = _unur_string_new();
  _unur_string_append( reason, "invalid data for %s '%s'", type, key );
  _unur_error_x( GENTYPE, file, line, "error", UNUR_ERR_STR_INVALID,reason->text);
  _unur_string_free( reason );
} /* end of _unur_str_error_invalid() */

/*---------------------------------------------------------------------------*/

void
_unur_str_error_args( const char *file, int line, const char *key )
     /*----------------------------------------------------------------------*/
     /* print error message: invalid argument string for set call.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   file ... file name (inserted by __FILE__)                          */
     /*   line ... line number in source file (inserted by __LINE__)         */
     /*   key  ... keyword                                                   */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *reason = _unur_string_new();
  _unur_string_append( reason, "invalid argument string for '%s'", key );
  _unur_error_x( GENTYPE, file, line, "error", UNUR_ERR_STR_INVALID,reason->text);
  _unur_string_free( reason );
} /* end of _unur_str_error_args() */

/*---------------------------------------------------------------------------*/
