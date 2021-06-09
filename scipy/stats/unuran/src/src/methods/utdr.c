/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      utdr.h                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection with three points of contact   *
 *              T(x) = -1/sqrt(x)     (T_c with c= -1/2)                     *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF and mode of a T-concave distribution                       *
 *      produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:  pointer to the density, mode of the density                  *
 *                                                                           *
 *   PARAMETERS:                                                             *
 *      double  il, ir       ... left and right boundary of domain           *
 *                               (default: +/- INFINITY)                     *
 *      double  pdf_area     ... area below PDF (need not be 1)              *
 *                               (default: 1.)                               *
 *      double *pdf_param    ... parameters of PDF                           *
 *                               (default: NULL)                             *
 *      int     n_pdf_param  ... number of parameters of PDF                 *
 *                               (default: 0)                                *
 *      double  c_factor     ... constant for choosing the constr. points    *
 *                               (default: 0.664)                            *
 *      double  delta_factor ... constant for approx. first derivative       *
 *                               (default: 0.00001)                          *
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
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] Hoermann W. (1995): A rejection technique for sampling from         *
 *       T-concave distributions, ACM TOMS 21, p. 182-193                    *
 *       (see Algorithm UTDR)                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * ..... beschreibung ....                                                   *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "utdr.h"
#include "utdr_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define UTDR_VARFLAG_VERIFY     0x01u   /* flag for verifying mode           */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define UTDR_DEBUG_REINIT    0x00000010u   /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define UTDR_SET_CPFACTOR       0x001u
#define UTDR_SET_DELTA          0x002u
#define UTDR_SET_PDFMODE        0x004u   /* PDF at mode is set               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "UTDR"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_utdr_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_utdr_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_utdr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_utdr_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_utdr_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_utdr_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_utdr_sample( struct unur_gen *generator );
static double _unur_utdr_sample_check( struct unur_gen *generator );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_utdr_hat( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute hat and squeezes.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_utdr_debug_init( const struct unur_gen *gen, 
				   double ttly, double ttlys, double ttry, double ttrys,
				   double cfac, int setupok, double c );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_utdr_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_utdr_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_utdr_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /* call to PDF         */

/*---------------------------------------------------------------------------*/

#define _unur_utdr_getSAMPLE(gen) \
   ( ((gen)->variant & UTDR_VARFLAG_VERIFY) \
     ? _unur_utdr_sample_check : _unur_utdr_sample )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_utdr_new( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*                                                                      */
     /* comment:                                                             */
     /*   if the area below the PDF is not close to 1 it is necessary to     */
     /*   set pdf_area to an approximate value of its area (+/- 30 % is ok). */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); 
    return NULL;
  }

  /**TODOWH:
     Kann man eigentlich einen bound fuer die area auch
     im Distribution-objekt eingeben????? 
     so dass es klar ist, dass da nur ein bound ist?

     JL: ich habe ich SSR die moeglichkeit eingebaut, dass nur eine
     obere schranke eingegeben wurde. dass einzige problem dabei ist, 
     dass dann der squeeze nicht mehr verwendet werden kann. 
     man kann aber die verwendung des squeezes (der eigentlich nur bei
     sehr teuren dichten etwas bringt, mittels 
     int unur_ssr_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
     abschalten.

     waere das ein brauchbarer weg fuer Dich?
     
     alternativ koennten wir auch ein unur_distr_set_approx_pdfarea()
     einbauen. die unterscheidung zwischen exactem wert und naehrung kann 
     man mittels set-flag machen. (irgendwie ist mir bis jetzt nie die idee
     gekommen, dass so zu machen.) ist das so besser?
  **/

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_utdr_par) );
  COOKIE_SET(par,CK_UTDR_PAR);

  /* copy input */
  par->distr       = distr;   /* pointer to distribution object              */

  /* set default values */
  PAR->c_factor     = 0.664; 
          /* optimal value for the normal distribution, which is good for 
	     all bell-shaped densities. The minimax approach for that 
	     transformation has c_factor=2. */
  PAR->delta_factor = 0.00001;
          /* constant for choosing delta to replace the tangent.
	     default should not be changed if used doubles have at least
	     10 decimal digits precision. */

  PAR->fm        = -1.;                /* PDF at mode (unknown)               */
  PAR->hm        = -1.;                /* square of PDF at mode (unknown)     */

  par->method   = UNUR_METH_UTDR;     /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_utdr_init;

  return par;

} /* end of unur_utdr_new() */

/*****************************************************************************/

int 
unur_utdr_set_pdfatmode( UNUR_PAR *par, double fmode )
     /*----------------------------------------------------------------------*/
     /* Set PDF at mode. if set, the PDF at the mode is never changed.       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   fmode ... PDF at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, UTDR );

  /* check new parameter for generator */
  if (fmode <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode)");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->fm = fmode;             /* PDF at mode */
  PAR->hm = -1./sqrt(fmode);   /* transformed PDF at mode */

  /* changelog */
  par->set |= UTDR_SET_PDFMODE;

  return UNUR_SUCCESS;

} /* end of unur_utdr_set_pdfatmode() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_set_cpfactor( struct unur_par *par, double cp_factor )
     /*----------------------------------------------------------------------*/
     /* set factor for position of left and right construction point         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   cp_factor ... factor                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, UTDR );

  /* check new parameter for generator */
  /** TODO: welche werte fuer c sind zulaessig / sinnvoll ? 
  zulaessig ist jedes c>0, man koennte damit falsche Flaechenangaben kompensieren.
  Wenn area genau bekannt ist, ist ein c > 2 (2 ist der minimax approach) so weit
  ich weiss nie sinnvoll. Ich denke aber, das sollte man besser nicht prinzipiell
  verbieten, hoechstens eine warnung.**/
  if (cp_factor <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp-factor <= 0");
    return UNUR_ERR_PAR_SET;
  }

  if (cp_factor > 2.1)
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp-factor > 2 not recommended. skip");

  /* store date */
  PAR->c_factor = cp_factor;

  /* changelog */
  par->set |= UTDR_SET_CPFACTOR;

  return UNUR_SUCCESS;

} /* end of unur_utdr_set_cpfactor() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_set_deltafactor( struct unur_par *par, double delta )
     /*----------------------------------------------------------------------*/
     /* set factor for replacing tangents by secants                         */
     /* (which then are move above the transformed density)                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   delta ... delta-factor                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, UTDR );

  /* check new parameter for generator */
  /** TODO: welche werte fuer delta sind zulaessig / sinnvoll ? 
      von einem rein numerischen Gesichtspunkt kann man sagen: Delta_factor
      soll nicht kleiner sein als sqrt(epsilon) (mit epsilon meine ich die
      kleinste Zahl die einen Unterschied zwischen 1 und 1+epsilon erkennen
      laesst). Dann ist garantiert, dass wir durch Ausloeschung nicht mehr
      als die haelfte der mantisse verlieren
      Von Algorithmus ueberlegungen her sollte delta>=xl*deltafactor sicher
      nicht groesser sein als ein Bruchteil (wohl hoechstens 1/10) zwischen
      mode und Design-punkt, weil sonst die Anahmewahrscheinlichkeit schlechter
      wird. Das ergibt sicher probleme fuer Dichten wo der mode sehr gross
      und der skale-parameter im Vergleich dazu klein ist.
      Ich habe den Algorithmus so geaendert, dass er, falls es dieses
      Problem gibt, eine (roundoff-problem) Warnung ausspuckt, und als
      Delta 1/100*|designpoint-mode| nimmt Empirisch habe ich dafuer
      festgestellt, dass damit die Flaeche unterm Hut Kaum mehr als 1%
      zunimmt und das ist ja nicht schlimm**/
  if (delta <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"delta <= 0");
    return UNUR_ERR_PAR_SET;
  }
  if (delta > 0.1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"delta must be small");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->delta_factor = delta;

  /* changelog */
  par->set |= UTDR_SET_DELTA;

  return UNUR_SUCCESS;

} /* end of unur_utdr_set_delta() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_set_verify( struct unur_par *par, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, UTDR );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | UTDR_VARFLAG_VERIFY) : (par->variant & (~UTDR_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_utdr_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_chg_verify( struct unur_gen *gen, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, UTDR, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;

  if (verify)
    /* turn verify mode on */
    gen->variant |= UTDR_VARFLAG_VERIFY;
  else
    /* turn verify mode off */
    gen->variant &= ~UTDR_VARFLAG_VERIFY;

  SAMPLE = _unur_utdr_getSAMPLE(gen);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_utdr_chg_verify() */

/*****************************************************************************/

int
unur_utdr_chg_pdfatmode( struct unur_gen *gen, double fmode )
     /*----------------------------------------------------------------------*/
     /* change value of PDF at mode                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   fmode ... PDF at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, UTDR, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (fmode <= 0.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"PDF(mode)");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  GEN->fm = fmode;             /* PDF at mode */
  GEN->hm = -1./sqrt(fmode);   /* transformed PDF at mode */

  /* changelog */
  gen->set |= UTDR_SET_PDFMODE;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_utdr_chg_pdfatmode() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_utdr_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_UTDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_UTDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_utdr_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* check parameters */
  if (_unur_utdr_check_par(gen) != UNUR_SUCCESS) {
    _unur_utdr_free(gen); return NULL;
  }

  /* create hat and squeeze (setup procedure) */
  if ( _unur_utdr_hat(gen)!=UNUR_SUCCESS ) {
    _unur_utdr_free(gen); return NULL;
  }

  /* hat successfully created */
  return gen;

} /* end of _unur_utdr_init() */

/*---------------------------------------------------------------------------*/

int
_unur_utdr_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* re-initialize (existing) generator.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int rcode;

  /* check parameters */
  if ( (rcode = _unur_utdr_check_par(gen)) != UNUR_SUCCESS)
    return rcode;

  /* update left and right boundary for algorithm */
  GEN->il = DISTR.BD_LEFT;
  GEN->ir = DISTR.BD_RIGHT;

  /* (re)set sampling routine */
  SAMPLE = _unur_utdr_getSAMPLE(gen);

  /* compute universal bounding rectangle */
  return _unur_utdr_hat( gen );
} /* end of _unur_utdr_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_utdr_create( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* allocate memory for generator                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to (empty) generator object with default settings          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_UTDR_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_utdr_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_UTDR_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_utdr_getSAMPLE(gen);
  gen->destroy = _unur_utdr_free;
  gen->clone = _unur_utdr_clone;
  gen->reinit = _unur_utdr_reinit;

  /* copy some parameters into generator object */
  GEN->il = DISTR.BD_LEFT;           /* left boundary of domain               */
  GEN->ir = DISTR.BD_RIGHT;          /* right boundary of domain              */
  GEN->fm = PAR->fm;                  /* PDF at mode                           */
  GEN->hm = PAR->hm;                  /* square root of PDF at mode            */
  GEN->c_factor = PAR->c_factor;
  GEN->delta_factor = PAR->delta_factor;

  /* initialize parameters */
  /** TODO !!! **/
  /** ist das wirklich so noetig ?? **/
  GEN->vollc = 0.; 
  GEN->volcompl = 0.; 
  GEN->voll = 0.; 
  GEN->al = 0.; 
  GEN->ar = 0.; 
  GEN->col = 0.; 
  GEN->cor = 0.; 
  GEN->sal = 0.; 
  GEN->sar = 0.; 
  GEN->bl = 0.; 
  GEN->br = 0.; 
  GEN->ttlx = 0.; 
  GEN->ttrx = 0.; 
  GEN->brblvolc = 0.; 
  GEN->drar = 0.; 
  GEN->dlal = 0.; 
  GEN->ooar2 = 0.; 
  GEN->ooal2 = 0.;
  /* constants of the hat and for generation*/

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_utdr_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_utdr_create() */

/*---------------------------------------------------------------------------*/

int
_unur_utdr_check_par( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* check parameters of given distribution and method                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{

  /* check for required data: mode */
  if (!(gen->distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode: try finding it (numerically)");
    if (unur_distr_cont_upd_mode(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode");
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }

  /* check for required data: area */
  if (!(gen->distr->set & UNUR_DISTR_SET_PDFAREA)) {
    if (unur_distr_cont_upd_pdfarea(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"area below PDF");
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }

  /* mode must be in domain */
  if ( (DISTR.mode < DISTR.BD_LEFT) ||
       (DISTR.mode > DISTR.BD_RIGHT) ) {
    /* there is something wrong.
       assume: user has change domain without changing mode.
       but then, she probably has not updated area and is to large */
    _unur_warning(GENTYPE,UNUR_ERR_GEN_DATA,"area and/or CDF at mode");
    DISTR.mode = _unur_max(DISTR.mode,DISTR.BD_LEFT);
    DISTR.mode = _unur_min(DISTR.mode,DISTR.BD_RIGHT);
  }

  return UNUR_SUCCESS;
} /* end of _unur_utdr_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_utdr_clone( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generator object                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
#define CLONE  ((struct unur_utdr_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_UTDR_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_utdr_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_utdr_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 

  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_UTDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_UTDR_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_utdr_free() */

/*****************************************************************************/

double
_unur_utdr_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double u,v,x,help,linx;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_UTDR_GEN,INFINITY);

  while (1) {
    /*2*/
    u = _unur_call_urng(gen->urng) * GEN->volcompl;
    /*2.1*/
    if (u <= GEN->voll) {
      u = GEN->voll-u; /*added to ensure inversion for the hat-generation*/
      x = -GEN->dlal+GEN->ooal2/(u-GEN->col);
      help = GEN->al*(u-GEN->col);
      linx = help*help;
    }
    else {
      if (u <= GEN->vollc) {
	x = (u-GEN->voll) * GEN->brblvolc + GEN->bl;
	linx = GEN->fm;
      }
      else {
	x = - GEN->drar - GEN->ooar2 / (u-GEN->vollc - GEN->cor);
	help = GEN->ar * (u-GEN->vollc - GEN->cor);
	linx = help*help;
      }
    }
    /*2.2*/
    v = _unur_call_urng(gen->urng) * linx;
    /*2.3*/
    if (x<DISTR.mode) {
      if (x >= GEN->ttlx) {
	help = GEN->hm - (DISTR.mode - x) * GEN->sal;
	if (v * help * help <= 1.) return x;
      } 
    }
    else {
      if (x <= GEN->ttrx) {
	help = GEN->hm - (DISTR.mode - x) * GEN->sar;
	if (v * help * help <= 1.) return x; 
      }
    }
    if (v <= PDF(x)) return x; 
  }

} /* end of _unur_utdr_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_utdr_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double u,v,x,help,linx,pdfx,squeezex;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_UTDR_GEN,INFINITY);
  
  while (1) {
    /*2*/
    u = _unur_call_urng(gen->urng) * GEN->volcompl;
    /*2.1*/
    if (u <= GEN->voll) {
      u = GEN->voll-u; /* added to ensure inversion for the hat-generation */
      x = -GEN->dlal+GEN->ooal2/(u-GEN->col);
      help = GEN->al*(u-GEN->col);
      linx = help*help;
    }
    else {
      if (u <= GEN->vollc) {
	x = (u-GEN->voll) * GEN->brblvolc + GEN->bl;
	linx = GEN->fm;
      }
      else {
	x = - GEN->drar - GEN->ooar2 / (u-GEN->vollc - GEN->cor);
	help = GEN->ar * (u-GEN->vollc - GEN->cor);
	linx = help*help;
      }
    }
    /*2.2*/
    v = _unur_call_urng(gen->urng) * linx;
    /*2.3*/
    squeezex=0.;
    if (x<DISTR.mode) {
      if (x >= GEN->ttlx) {
	help = GEN->hm - (DISTR.mode - x) * GEN->sal;
        squeezex=1./(help*help);
	/*        if (v * help * help <= 1.) return x;*/
      } 
    }
    else {
      if (x <= GEN->ttrx) {
	help = GEN->hm - (DISTR.mode - x) * GEN->sar;
        squeezex=1./(help*help);
	/*        if (v * help * help <= 1.) return x; */
      }
    }
    /*evaluate density-function*/
    pdfx=PDF(x);
    
    /* verify hat function */
    if(_unur_FP_less(linx,pdfx))
      { _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
      _unur_log_printf(gen->genid,__FILE__,__LINE__,"x %e PDF(x) %e hat(x) %e squeeze(x) %e", \
		       x,pdfx,linx,squeezex ); 
      }
    /* verify squeeze function */
    if(_unur_FP_less(pdfx,squeezex))
      { _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) < squeeze(x)");
      _unur_log_printf(gen->genid,__FILE__,__LINE__,"x %e PDF(x) %e hat(x) %e squeeze(x) %e", \
		       x,pdfx,linx,squeezex ); 
      }
    if (v <= PDF(x)) return x;
  }
  
} /* end of _unur_utdr_sample_check() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/** TODO gibts da schon eine UNURAN constante? ich hab sie nicht gefunden!! **/
#define SMALL_VAL 1.e-50
/** JL: was soll die leisten?
    wenn sie gebraucht wird, wuerde ich sie mit einer entsprenden beschreibung
    in unuran_config.h stellen.
**/

int
_unur_utdr_hat( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute hat and squeeze                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{ 
  double fm;

  int setupok=1;
  double c,cfac,volc,volr,ttly,ttlys,ttry,ttrys,dl,dr,delta,delta1,delta2,pdfx;

  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen, CK_UTDR_GEN, UNUR_ERR_COOKIE );

  /* compute PDF at mode (if not given by user) */
  if (!(gen->set & UTDR_SET_PDFMODE)) {
    fm = PDF(DISTR.mode);
    /* fm must be positive */
    if (fm <= 0.) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"PDF(mode) <= 0.");
      return UNUR_ERR_GEN_DATA;
    }
    /* step 1.0 of algorithm UTDR */
    GEN->fm = fm;           /* PDF at mode  */
    GEN->hm = -1/sqrt(fm);  /* transformed PDF at mode  */
  }

  /** TODO: inititialisieren notwendig ?? **/

  /** ich habe diese variablen initialisiert, weil ich auf einem
      alpha rechner eine flaoting exception bei der ersten verwendung
      von ttry oder ttrys bekommen habe.
      das ist keine IEEE 764 architektur.
      es stand dann ttry irgend ein bitmuster, dass keiner regulaeren
      zahl entsprochen hat und daher die FPE versucht hat
      (in eine if abfrage !!)
      das zeigt, das da ttry tatsaechlich nicht initialiesiert verwendet 
      wurde. das sollte man sauberer programmieren.
  **/

  ttry = 0.;
  ttrys = 0.;
  ttly = 0.;
  ttlys = 0.;
  dl = 0.;
  dr = 0.;
  volr = 0.;

  /* start of the set-up procedure */

  /* step 1.0 of algorithm UTDR */
  /* see above or in unur_utdr_set_pdfatmode() */

  do {

    /* 1.1 */
    cfac = (setupok) ? GEN->c_factor : 2.;     /* gibt es hier nur zwei varianten ?? ja*/
    c = cfac * DISTR.area/GEN->fm;
    setupok=1;         

    GEN->ttlx = DISTR.mode - c;
    GEN->ttrx = DISTR.mode + c;

    /* 1.2 */
    /** TODO: kann man das nicht loeschen ?? **/
    if (/*GEN->il > -INFINITY &&*/ GEN->ttlx < GEN->il) { 
      /* this is the case of no left tail*/
      GEN->bl = GEN->il;
      GEN->al = 0.;
      GEN->voll = 0.;
      if (GEN->il < DISTR.mode) {
	/* if the left domain border is left of the mode we set ttlx
	   only to use it for the squeeze*/
        GEN->ttlx = DISTR.mode + (GEN->il - DISTR.mode) * 0.6;
        pdfx=PDF(GEN->ttlx);
        if (pdfx > SMALL_VAL)
          GEN->sal = (GEN->hm + 1./sqrt(pdfx)) / (DISTR.mode - GEN->ttlx);
        else 
	  GEN->ttlx = DISTR.mode;
	/* pdfx is too small. We set ttlx=mode that no squeeze is used */
      }  
    }
    else {
     ttlys = PDF(GEN->ttlx);
     if (ttlys < SMALL_VAL) { 
       /* in this case we cut off the left tail*/
       GEN->il = GEN->ttlx;
       GEN->bl = GEN->il;
       GEN->al = 0.;
       GEN->voll = 0.;
       GEN->ttlx=DISTR.mode;
       /* pdfx is too small. We set ttlx=mode that no squeeze is used */
     }
     else {
       ttlys = -1./sqrt(ttlys);
       GEN->sal =  (GEN->hm - ttlys) / (DISTR.mode - GEN->ttlx);

       /* delta1> 0 as -ttlys>0 und sal >0; */
       delta2 = ( GEN->sal > 0. ) ? -ttlys/GEN->sal : -ttlys;
       delta1 = fabs(GEN->ttlx);
       delta = GEN->delta_factor * ((delta1<=delta2) ? delta2 : delta1);
       if (delta > c * 0.01) {
	 delta = UNUR_SQRT_DBL_EPSILON * ((delta1<=delta2) ? delta2 : delta1);
	 /* Using sqrt(DBL_EPSILON) as delta_factor guarantees that not more than
	    half of the precision will be lost when computing the slope GEN->al */
	 if (delta > c * 0.01) {
	   delta = c * 0.01;
	   /* delta is forced to be c * 0.01 although this can 
	      result in numerical inaccuracies when computing 
	      the slope GEN->al. Therefore a warning is
	      issued */
	   _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,
			 "Delta larger than c/100!!, perhaps you can use a mode closer to 0 to remove this problem?");
         }
       }
       
       ttly = -1./sqrt(PDF(GEN->ttlx+delta));
       GEN->al = (ttly-ttlys)/delta;

       if (GEN->al <= 0.) 
	 /* setupok==0 means that this setup is quitted and set-up is restarted
	    with different value for cfac */
	 setupok = 0; 
       else {
	 GEN->bl = GEN->ttlx + (GEN->hm - ttly)/GEN->al;
	 dl = ttly - GEN->al * GEN->ttlx;
	 GEN->voll = -1./(GEN->al * GEN->hm);
	 GEN->col = GEN->voll;
	 if (GEN->il > -INFINITY)
	   GEN->voll += 1./(GEN->al * (GEN->al * GEN->il + dl));
       }
     }
    }

    /* 1.3 */
    if(setupok) {
      if (/*GEN->ir < INFINITY &&*/ GEN->ttrx > GEN->ir) {
	/* this is the case of no right tail */
        GEN->br = GEN->ir;
        GEN->ar = 0.;
        volr = 0.;
        if (GEN->ir > DISTR.mode) {
	  /* if the right domain border is right of the mode we set ttrx
	     only to use it for the squeeze */
          GEN->ttrx = DISTR.mode + (GEN->ir - DISTR.mode) * 0.6;
          pdfx = PDF(GEN->ttrx);
          if (pdfx > SMALL_VAL)
            GEN->sar = (GEN->hm + 1./sqrt(PDF(GEN->ttrx))) / (DISTR.mode - GEN->ttrx);
          else 
	    GEN->ttrx = DISTR.mode;
	  /* pdfx is too small. We set ttrx=mode that no squeeze is used */
        } 
      }
      else {
        ttrys = PDF(GEN->ttrx);
        if (ttrys < SMALL_VAL){
	  /* in this case we cut off the right tail */
          GEN->ir = GEN->ttrx;
          GEN->br = GEN->ir;
          GEN->ar = 0.;
          volr = 0.;
          GEN->ttrx = DISTR.mode;
	  /* pdfx is too small. We set ttrx=mode that no squeeze is used */
	}
	else {
	  ttrys= -1./sqrt(ttrys);
	  /* see 1.2. for explanations */
	  GEN->sar = (GEN->hm - ttrys) / (DISTR.mode - GEN->ttrx);
	  /* delta is positive, da ttrys<0 und sar <0 */
	  delta2 = (GEN->sar<0.) ? ttrys/GEN->sar : -ttrys;
	  delta1 = fabs(GEN->ttrx);
	  delta = GEN->delta_factor * ((delta1<=delta2) ? delta2 : delta1);
	  if (delta > c*0.01) { 
	    delta = UNUR_SQRT_DBL_EPSILON * ((delta1<=delta2) ? delta2 : delta1);
	    /* Using sqrt(DBL_EPSILON) as delta_factor guarantees that not more than
	       have of the precision will be lost when computing the slope GEN->al */
	    if (delta > c*0.01) {
	      delta=c*0.01;
	      /* delta is forced to be c*0.01 allthough this can result in numerical
		 inaccuracies when computing the slope GEN->al. Therefore a warning is
		 issued*/
	      _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,
			    "Delta larger than c/100!!, perhaps you can use a mode closer to 0 to remove this problem?");
	    }
	  }
	  
	  ttry = -1./sqrt(PDF(GEN->ttrx-delta));
	  GEN->ar = (ttrys - ttry)/delta;
	  if (GEN->ar >= 0.) 
	    /* setupok==0 means that this setup is quitted and set-up is 
	       restarted with different value for cfac */
	    setupok = 0;
	  else { 
	    GEN->br = GEN->ttrx + (GEN->hm - ttry) / GEN->ar;
	    dr = ttry - GEN->ar * GEN->ttrx;
	    volr = 1./(GEN->ar * GEN->hm);
	    GEN->cor = volr;
	    if (GEN->ir<INFINITY)
	      volr -= 1./(GEN->ar * (GEN->ar * GEN->ir + dr));
	  }
	}
      }
    }

    /* 1.4 */
    if(setupok) {
      volc = (GEN->br - GEN->bl) * GEN->fm;
      GEN->vollc = GEN->voll + volc;
      GEN->volcompl = GEN->vollc + volr;
      if (volc>0.) 
        GEN->brblvolc = (GEN->br - GEN->bl)/volc;
      if (!_unur_iszero(GEN->ar)) {
        GEN->drar = dr/GEN->ar;
        GEN->ooar2 = 1./(GEN->ar*GEN->ar);
      }
      if (!_unur_iszero(GEN->al)) {
        GEN->dlal = dl/GEN->al;
        GEN->ooal2 = 1./(GEN->al*GEN->al);
      }
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_utdr_debug_init(gen,ttly,ttlys,ttry,ttrys,cfac,setupok,c);
#endif

    if (!_unur_isfsame(cfac,2.)) {
      if(setupok)
        if (GEN->volcompl > 4. * DISTR.area || GEN->volcompl < 0.5 * DISTR.area)
        setupok=0;
    }
    else { 
      if (setupok==0 || GEN->volcompl > 8. * DISTR.area || GEN->volcompl < 0.5 * DISTR.area) {
        _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"; Area below hat too large or zero!! possible reasons: PDF, mode or area below PDF wrong;  density not T-concave\n");
        return 0;
      }
    }

  } while (!setupok);

  return UNUR_SUCCESS;

} /* end of _unur_utdr_hat() */

#undef SMALL_VAL

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_utdr_debug_init( const struct unur_gen *gen,
		       double ttly, double ttlys, double ttry, double ttrys,
		       double cfac, int setupok, double c )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_UTDR_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = utdr(transformed density rejection with 3 points of contact)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_utdr_sample",gen->genid);
  if (gen->variant & UTDR_VARFLAG_VERIFY)
    fprintf(LOG,"_check()\n");
  else
    fprintf(LOG,"()\n");
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: Data for hat and squeeze:\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s:\tc_factor=%e delta_factor=%e real c=%e\n",gen->genid,GEN->c_factor,GEN->delta_factor,c);
  fprintf(LOG,"%s:\ttlx=%e bl=%e mode=%e\n",gen->genid,GEN->ttlx,GEN->bl,DISTR.mode);
  fprintf(LOG,"%s:\tbr=%e trx=%e\n",gen->genid,GEN->br,GEN->ttrx);
  fprintf(LOG,"%s:\ttly=%e tlys=%e al=%e \n",gen->genid,ttly,ttlys,GEN->al);
  fprintf(LOG,"%s:\ttry=%e trys=%e ar=%e \n",gen->genid,ttry,ttrys,GEN->ar);
  fprintf(LOG,"%s:\tcfac=%e setupok=%d volcompl=%e pdf_area=%e\n",gen->genid,cfac,setupok,GEN->volcompl,DISTR.area);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);

} /* end of _unur_utdr_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_utdr_info( struct unur_gen *gen, int help )
     /*----------------------------------------------------------------------*/
     /* create character string that contains information about the          */
     /* given generator object.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   help ... whether to print additional comments                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   mode      = %g   %s\n", unur_distr_cont_get_mode(distr),
		      (distr->set & UNUR_DISTR_SET_MODE_APPROX) ? "[numeric.]" : "");
  _unur_string_append(info,"   area(PDF) = %g\n", DISTR.area);
  _unur_string_append(info,"\n");

  /* method */
  _unur_string_append(info,"method: UTDR (Universal Transformed Density Rejection -- 3 point method)\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   rejection constant = %.2f  [approx.]\n",
		      unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize));
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   deltafactor = %g  %s\n", GEN->delta_factor,
			(gen->set & UTDR_SET_DELTA) ? "" : "[default]");

    if (gen->set & UTDR_SET_PDFMODE)
      _unur_string_append(info,"   pdfatmode = %g\n", GEN->fm);
    
    if (gen->set & UTDR_SET_CPFACTOR)
      _unur_string_append(info,"   cpfactor = %g\n", GEN->c_factor);

    if (gen->variant & UTDR_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");

    _unur_string_append(info,"\n");
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_utdr_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
