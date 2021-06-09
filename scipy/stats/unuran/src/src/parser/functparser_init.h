/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_init.c                                           *
 *                                                                           *
 *   Init and destroy function tree.                                         *
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

struct ftreenode *
_unur_fstr2tree (const char *functstr)
     /*----------------------------------------------------------------------*/
     /* Compute funtion tree from string.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   functstr ... string containing function definition                 */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to root of function tree                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  return _unur_fstr_2_tree( functstr, FALSE );
} /* end of _unur_fstr2tree() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_fstr2tree_DefFunct (const char *functstr)
     /*----------------------------------------------------------------------*/
     /* Compute funtion tree from string.                                    */
     /* (Same as _unur_fstr2tree() but string must start with "f(x)=".       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   functstr ... string containing function definition                 */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to root of function tree (function term only!)             */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  return _unur_fstr_2_tree( functstr, TRUE );
} /* end of _unur_fstr2tree_DefFunct() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_fstr_dup_tree (const struct ftreenode *root)
     /*----------------------------------------------------------------------*/
     /* Duplicate function tree rooted at root                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   root ... pointer to root of function tree                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to duplicated tree                                         */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *dup;

  if (root==NULL) return NULL;

  /* check arguments */
  COOKIE_CHECK(root,CK_FSTR_TNODE,NULL);

  dup = _unur_xmalloc(sizeof(struct ftreenode));
  memcpy(dup,root,sizeof(struct ftreenode));
  if (root->left)  dup->left  = _unur_fstr_dup_tree(root->left);
  if (root->right) dup->right = _unur_fstr_dup_tree(root->right);

  return dup;

} /* end of _unur_fstr_dup_tree() */

/*---------------------------------------------------------------------------*/

void
_unur_fstr_free (struct ftreenode *root)  
     /*----------------------------------------------------------------------*/
     /* Destroy function tree.                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   root ... pointer to root of function tree                          */
     /*----------------------------------------------------------------------*/
{ 
  if( root != NULL ) {
    /* check arguments */
    COOKIE_CHECK(root,CK_FSTR_TNODE,RETURN_VOID);

    if (root->left)  _unur_fstr_free(root->left);
    if (root->right) _unur_fstr_free(root->right);
    free(root); 
  } 
} /* end of _unur_fstr_free() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Auxilliary routines                                                     **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_fstr_2_tree (const char *functstr, int withDefFunct)
     /*----------------------------------------------------------------------*/
     /* Compute funtion tree from string.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   functstr     ... string containing function definition             */
     /*   withDefFunct ... whether string is assumed to start with "f(x)="   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to root of function tree                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct parser_data *pdata;
  struct ftreenode *root;

  /* check arguments */
  _unur_check_NULL( GENTYPE,functstr,NULL );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag)
    _unur_fstr_debug_input(functstr);
#endif

  /* initialize parser */
  pdata = _unur_fstr_parser_init(functstr);

  /* check for errors (eg. empty string) */
  if (pdata == NULL)
    return NULL;

  /* tokenize function string */
  _unur_fstr_tokenize(pdata);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag)
    _unur_fstr_debug_token(pdata);
#endif

  /* exit in case of error */
  if (pdata->perrno) {
    _unur_fstr_parser_free(pdata);
    return NULL;
  }

  /* parse list of token */
  if (withDefFunct) {
    struct ftreenode *tmp = _unur_FunctDefinition(pdata);
    root = tmp->right;
    /* clear left subtree */
    _unur_fstr_free(tmp->left);
    free(tmp);
  }
  else {
    root = _unur_Expression(pdata);
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if ((_unur_default_debugflag & UNUR_DEBUG_SETUP) && root)
    _unur_fstr_debug_tree(pdata,root);
#endif

  /* check for possible errors */
  if (pdata->tno < pdata->n_tokens && !pdata->perrno)
    _unur_fstr_error_parse(pdata,ERR_UNFINISHED,__LINE__); 
  if (pdata->perrno) {
    _unur_fstr_parser_free(pdata);
    _unur_fstr_free(root);
    return NULL;
  }

  /* free working space */
  _unur_fstr_parser_free(pdata);

  /* return pointer to function tree */
  return root; 
} /* end of _unur_fstr_2_tree() */

/*---------------------------------------------------------------------------*/
