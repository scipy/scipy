/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mvtdr_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method MVTDR                              *
 *         (Multi-Variate Transformed Density Rejection)                     *
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


/*---------------------------------------------------------------------------*/
/* Information for constructing the generator                                */

struct unur_mvtdr_par { 

  int max_cones;                  /* maximum number of cones (at least 2^(N+T_STEPS_MIN) */
  int steps_min;                  /* minimum number of triangulation steps   */
  double bound_splitting;         /* bound for splitting cones               */

/* #if MODE == 1 */
/*   double mode_to_boundary;        /\* move mode to boundary if |mode - boundary| / length < MODE_TO_BOUNDARY *\/ */
/* #endif */

};


/*---------------------------------------------------------------------------*/

typedef struct s_vertex           /* -- data for one vertex ---------------- */
{
  struct s_vertex *next;          /* pointer to next vertex in list          */
  int index;                      /* index of vertex                         */
  double *coord;                  /* coordinates of spanning vector(norm = 1), "vertex" */
  double norm;                    /* norm of vertex                          */
} VERTEX;


typedef struct s_cone             /* -- a cone ----------------------------- */
{
  struct s_cone *next;            /* pointer to next cone in list            */
  int level;                      /* level of triangulation                  */
  VERTEX **v;                     /* list of vertices of the cone            */
  double *center;                 /* barycenter of cone                      */
  double logdetf;                 /* log determinant -log((dim-1)!) for cone */
  double alpha;                   /* parameter alpha for hat function        */
  double beta;                    /* parameter alpha for hat function        */
  double *gv;                     /* <g,v> for all vertices v                */
  double logai;                   /* (log of) coefficient for marginal density */
  double tp;                      /* coordinate of touching point            */
  double Hi;                      /* volume under hat in cone                */
  double Hsum;                    /* accumulated sum of volumes              */
  double Tfp;                     /* value of transformed density at touching point */
  double height;                  /* height of pyramid                       */
} CONE;

typedef struct s_edge_table       /* -- hash table for edges --------------- */
{
  int  index[2];                  /* index of incident vertices              */
  VERTEX *vertex;                 /* index of corresponding vertex (=barycenter) */
  struct s_edge_table *next;      /* next entry in list                      */
} E_TABLE;


typedef struct s_tp_arg           /* -- argument for tp function ----------- */
{
  double t;                       /* touching point                          */
  double logH;                    /* log of volume below hat                 */
  CONE *c;                        /* parameters                              */
  UNUR_GEN *gen;                  /* pointer to MVTDR generator object       */
  int status;                     /* status of cone parameters               */
} TP_ARG;


enum {                            /* -- possible status of cone paramters -- */
  MVTDR_CONE_OK      = 0x000,     /* cone is ready to use                    */
  MVTDR_CONE_DOMAIN  = 0x001,     /* touching point out of support of PDF    */
  MVTDR_CONE_INVALID = 0x002      /* parameters invalid (Hi not finite?)     */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_mvtdr_gen { 
  int  dim;                       /* dimension of distribution               */
  int  has_domain;                /* whether the domain of distribution has given domain */
  double max_gamma;               /* upper bound for gamma variaties         */

  const double *center;           /* center of distribution                  */

  CONE *cone;                     /* root of list of cones                   */
  CONE *last_cone;                /* pointer to last cone in list            */
  int n_cone;                     /* number of cones                         */
  int max_cones;                  /* maximum number of cones                 */
  double bound_splitting;         /* bound for splitting cones               */

  VERTEX *vertex;                 /* root of list of vertices                */
  VERTEX *last_vertex;            /* pointer to last vertex in list          */
  int n_vertex;                   /* number of vertices                      */

  E_TABLE **etable;               /* pointer to edge table                   */
  int etable_size;                /* size of edge table                      */

  CONE **guide;                   /* pointer to guide table                  */
  int guide_size;                 /* size of guide table                     */

  double *S;                      /* working array for storing point on simples */
  double *g;                      /* working array for vector g (direction of sweeping plane) */
  double *tp_coord;               /* working array for storing coordinates of touching point of hat */
  double *tp_mcoord;              /* working array for storing coordinates of touching point of hat moved into center */
  double *tp_Tgrad;               /* working array for storing gradient of transformed density at tp */

  double Htot;                    /* total volume below hat                  */
  int steps_min;                  /* minimum number of triangulation steps   */
  int n_steps;                    /* (highest) number of triangulation steps */

  double pdfcenter;               /* PDF at center                           */
};

/*---------------------------------------------------------------------------*/
