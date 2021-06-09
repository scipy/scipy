/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_eval.c                                           *
 *                                                                           *
 *   Evaluate function tree for given argument x.                            *
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

/*****************************************************************************/
/** API                                                                     **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double
_unur_fstr_eval_tree (const struct ftreenode *root, double x)
     /*----------------------------------------------------------------------*/
     /* Evaluate function tree at x                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   root ... pointer to root of function tree                          */
     /*   x    ... argument for which function should be evaluated           */
     /*                                                                      */
     /* return:                                                              */
     /*   result of computation                                              */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{  
  /* check arguments */
  CHECK_NULL(root,INFINITY);   COOKIE_CHECK(root,CK_FSTR_TNODE,INFINITY);
  return _unur_fstr_eval_node( root, x );
} /* end of _unur_fstr_eval_tree() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Routines for evaluating nodes of the function tree                      **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#define CHECK_INF(x)    if(_unur_FP_is_infinity((x))) return INFINITY;
#define CHECK_INFS(l,r) do { CHECK_INF((l)); CHECK_INF((r)); } while(0)

#define NODE_ARGS  double l ATTRIBUTE__UNUSED, double r ATTRIBUTE__UNUSED
/*---------------------------------------------------------------------------*/

double v_dummy  (NODE_ARGS) { return 0.; }
double v_const  (NODE_ARGS) { return 0.; }  /* nothing to do, value in node */

double v_less   (NODE_ARGS) { return (double)(l <  r); }
double v_equal  (NODE_ARGS) { return (double)(_unur_FP_same(l,r)); }
double v_greater(NODE_ARGS) { return (double)(l >  r); }
double v_less_or(NODE_ARGS) { return (double)(l <= r); }
double v_unequal(NODE_ARGS) { return (double)(!_unur_FP_same(l,r)); }
double v_grtr_or(NODE_ARGS) { return (double)(l >= r); }

double v_plus   (NODE_ARGS) { return (l + r); }
double v_minus  (NODE_ARGS) { return (l - r); }
double v_mul    (NODE_ARGS) { return (l * r); }
double v_div    (NODE_ARGS) { return (l / r); }
double v_power  (NODE_ARGS) { return pow(l,r); }

double v_mod    (NODE_ARGS) { return (double)((int)l % (int)r); }
double v_exp    (NODE_ARGS) { return exp(r); }
double v_log    (NODE_ARGS) { return (r<=0.) ? INFINITY : log(r); }
double v_sin    (NODE_ARGS) { CHECK_INF(r); return sin(r); }
double v_cos    (NODE_ARGS) { CHECK_INF(r); return cos(r); }
double v_tan    (NODE_ARGS) { CHECK_INF(r); return tan(r); }
double v_sec    (NODE_ARGS) { double cosr; CHECK_INF(r); cosr=cos(r); 
                                               return _unur_iszero(cosr) ? INFINITY : 1./cosr; }
double v_sqrt   (NODE_ARGS) { return (r<0.) ? INFINITY : sqrt(r); }
double v_abs    (NODE_ARGS) { return fabs(r); }
double v_sgn    (NODE_ARGS) { return ((r<0.) ? -1. : ((r>0.) ? 1. : 0.)); }

/*---------------------------------------------------------------------------*/
#undef CHECK_INF
#undef CHECK_INFS
#undef NODE_ARGS
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Evaluate function                                                       **/
/*****************************************************************************/

double
_unur_fstr_eval_node (const struct ftreenode *node, double x)
     /*----------------------------------------------------------------------*/
     /* Evaluate function tree starting from `node' at x                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   root ... pointer to node in function tree                          */
     /*   x    ... argument for which function should be evaluated           */
     /*                                                                      */
     /* return:                                                              */
     /*   result of computation                                              */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  double val_l, val_r;

  /* check arguments */
  CHECK_NULL(node,INFINITY);   COOKIE_CHECK(node,CK_FSTR_TNODE,INFINITY);

  switch (node->type) {
  case S_UCONST:
  case S_SCONST:
    /* node contains constant */
    return node->val;

  case S_UIDENT:
    /* variable */
    return x;

  default:
    /* use evaluation function */
    /* compute values at leaves */
    val_l = (node->left)  ? _unur_fstr_eval_node(node->left, x) : 0. ;
    val_r = (node->right) ? _unur_fstr_eval_node(node->right,x) : 0. ;
    return (*symbol[node->token].vcalc)(val_l,val_r);
  }
} /* end of _unur_fstr_eval_node() */

/*---------------------------------------------------------------------------*/
