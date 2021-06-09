/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_parser.c                                         *
 *                                                                           *
 *   Parse tokenized function string and construct function tree             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold             *
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
/** Parser                                                                  **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_FunctDefinition (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /* Get user defined function.                                           */
     /*                                                                      */
     /*  FunctDefinition ::= DefFunctDesignator '=' Expression               */
     /*                                                                      */
     /*                    '='                                               */
     /*                   /   \                                              */
     /*  DefFunctDesignator     Expression                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token; 

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* left hand side: DefFunctDesignator */
  left = _unur_DefFunctDesignator(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(left); return NULL;
  }

  /* next token must be "=" sign */
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       strcmp(symb,"=") != 0 ) {
    _unur_fstr_free(left);
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_EQUAL,__LINE__); 
  }

  /* right hand side: function term */
  right = _unur_Expression(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(left); _unur_fstr_free(right);
    return NULL;
  }

  /* store function in node */
  node = _unur_fstr_create_node(symb,0.,token,left,right); 

  /* return pointer to function */
  return node; 
} /* end of _unur_FunctDefinition() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_DefFunctDesignator (struct parser_data *pdata) 
     /*----------------------------------------------------------------------*/
     /* Definition of user defined function.                                 */
     /*                                                                      */
     /*  DefFunctDesignator ::= Identifier '(' DefParameterlist ')'          */
     /*                                                                      */
     /*       Identifier                                                     */
     /*      /          \                                                    */
     /*  NULL            DefParameterlist                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *params; 
  char             *fsymb, *symb;
  int              n_params; 
  int              funct, token;

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* get function identifier */
  if ( _unur_fstr_next_token(pdata,&funct,&fsymb) != UNUR_SUCCESS ||
       symbol[funct].type != S_UFUNCT )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_FUNCT,__LINE__); 

  /* read opening parenthesis '(' */
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       symb[0] != '(' )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_OPEN_P,__LINE__); 

  /* read the parameter list */
  params = _unur_DefParameterlist(pdata,&n_params);
  if (pdata->perrno) {
    _unur_fstr_free(params);
    return NULL;
  }

  /* read closing parenthesis ')' */
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       symb[0] != ')' ) {
    _unur_fstr_free(params);
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_CLOSE_P,__LINE__); 
  }

  /* store function header in node */
  node = _unur_fstr_create_node(fsymb,0.,funct,NULL,params); 

  /* return pointer to function desigatior */
  return node;
} /* end of _unur_DefFunctDesignator() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_DefParameterlist(struct parser_data *pdata, int *n_params) 
     /*----------------------------------------------------------------------*/
     /* Parameter list for user defined function.                            */
     /*                                                                      */
     /*  DefParameterlist ::= '(' Identifier [ ',' Identifier ] ')'          */
     /*                                                                      */
     /*       Identifier                                 ','                 */
     /*      /          \      or:                      /   \                */
     /*  NULL            NULL      more identifiers tree     Identifier      */
     /*                                                     /          \     */
     /*                                                 NULL            NULL */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*   n_params  ... detected number of parameters for defined function   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* read user defined identifier, i.e. a variable */ 
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       symbol[token].type != S_UIDENT )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_VAR,__LINE__);

  /* make node for first parameter of function and set */
  /* counter for parameters to 1                       */
  node = _unur_fstr_create_node(symb,0.,token,NULL,NULL); 
  *n_params = 1;

  /* scan token list while we find a list separator `,' */
  while ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
	  symb[0] == ',' ) {

    /* get next variable */
    if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
	 symbol[token].type != S_UIDENT ) {
      _unur_fstr_free(node);
      return _unur_fstr_error_parse(pdata,ERR_EXPECT_VAR,__LINE__);
    }

    /* old node becomes left node of `,' node */
    left = node; 

    /* make node for next variable (becomes right node) */
    /* and update counter for parameters                */
    right = _unur_fstr_create_node(symb,0.,token,NULL,NULL); 
    (*n_params)++; 

    /* make node for `,' separator */
    node = _unur_fstr_create_node(",",0.,s_comma,left,right); 
  }

  /* set token pointer to first element after parameter list */
  --(pdata->tno);

  /* return pointer to parameter list */
  return node; 
} /* end of _unur_DefParameterlist() */ 

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_Expression (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /* Expression ::= SimpleExpression [ RelationOperator SimpleExpression ] */ 
     /*                                                                      */
     /*                                    RelationOperator                  */
     /* SimpleExpression  or:             /                \                 */
     /*                   SimpleExpression                  SimpleExpression */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* read simple expression from function string */
  left = _unur_SimpleExpression(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(left); return NULL; 
  }

  /* get next token */
  if ( _unur_fstr_next_token(pdata,&token,&symb) == UNUR_SUCCESS &&
       symbol[token].type == S_REL_OP ) {
    /* relation operator  --> read r.h.s.*/
    right = _unur_SimpleExpression(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(left); _unur_fstr_free(right);
      return NULL;
    } 
    /* create a new node */
    node = _unur_fstr_create_node(symb,0.,token,left,right); 
  }

  else {
    /* we only have a simple expression */
    /* hence we set token pointer to get same token again */
    --(pdata->tno);
    node = left; 
  } 

  /* return pointer to Expression */
  return node; 
} /* end of _unur_Expression() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_SimpleExpression (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  SimpleExpression ::= STerm { AddingOperator Term }                  */
     /*                                                                      */
     /*                                AddingOperator                        */
     /*  STerm  or:                   /              \                       */
     /*                more terms tree                Term                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* get next Term in string */
  node = _unur_STerm(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(node); return NULL;
  }

  /* get next token */
  while ( _unur_fstr_next_token(pdata,&token,&symb) == UNUR_SUCCESS &&
	  symbol[token].type == S_ADD_OP) {
    /* get Term after adding operator and        */
    /* use node a right node of new operator node */
    left = node; 

    right = _unur_Term(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(left); _unur_fstr_free(right);
      return NULL;
    }

    node = _unur_fstr_create_node(symb,0.,token,left,right); 
  }
  
  /* set token pointer to get same token again */
  --(pdata->tno);

  /* return pointer */
  return node; 
} /* end of _unur_SimpleExpression() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_STerm (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  STerm ::= [ '+' | '-' ] Term                                        */
     /*                                                                      */
     /*                        '-'                                           */
     /*  Term  or:            /   \                                          */
     /*                    '0'     Term                                      */
     /*                   /   \                                              */
     /*               NULL     NULL                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb; 
  int              token;

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* get next token */
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS)
    /* there is no next token */
    return NULL;

  if ( symb[0] == '-' ) {
    /* term with negative sign                 */
    /* "-" is interpreted as binary operator   */
    /* thus "0" is added in front of it         */
    left = _unur_fstr_create_node(NULL,0.,s_uconst,NULL,NULL); 
    right = _unur_Term(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(left); _unur_fstr_free(right);
      return NULL; 
    }

    node = _unur_fstr_create_node(symb,0.,token,left,right); 
  }

  else {
    /* term with positive sign or without any sign */
    if( symb[0] != '+' ) {
      /* set pointer to previous token again */
      --(pdata->tno);
    }
    node = _unur_Term(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(node); return NULL; 
    }
  }

  /* return pointer to term */
  return node; 
} /* end of _unur_STerm() */ 

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_Term (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  Term ::= Factor [ MultiplyingOperator Factor ]                      */
     /*                                                                      */
     /*                                   MultiplyingOperator                */
     /*  Factor  or:                     /                   \               */
     /*                 more factors tree                     Factor         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* get next factor of multiplication */
  node = _unur_Factor(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(node); return NULL;
  }

  /* get next token */
  while ( _unur_fstr_next_token(pdata,&token,&symb) == UNUR_SUCCESS &&
	  symbol[token].type == S_MUL_OP ) {
    /* get Factor after multiplication operator and */
    /* use node a right node of new operator node   */
     left = node; 

     right = _unur_Factor(pdata);
     if (pdata->perrno) {
       _unur_fstr_free(left); _unur_fstr_free(right);
       return NULL;
     }

     node = _unur_fstr_create_node(symb,0.,token,left,right); 
  }

  /* set token pointer to get same token again */
  --(pdata->tno);

  /* return pointer to term */
  return node; 
} /* end of _unur_Term() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_Factor (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  Factor ::= Base [ '^' Exponent ]                                    */
     /*                                                                      */
     /*                          '^'                                         */
     /*  Bas_Exp  or:           /   \                                        */
     /*                  Bas_Exp     Bas_Exp                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* get base of factor */
  left = _unur_Bas_Exp(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(left); return NULL;
  }
  
  /* get next token */
  if ( _unur_fstr_next_token(pdata,&token,&symb) == UNUR_SUCCESS &&
       symb[0] == '^' ) {
    /* get exponent of factor */
    right = _unur_Bas_Exp(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(left); _unur_fstr_free(right);
      return NULL;
    }

    /* and create node for '^' operator */
    node = _unur_fstr_create_node(symb,0.,token,left,right); 
  }

  else {
    /* no exponent                         */
    /* set pointer to previous token again */
    --(pdata->tno);
    /* and return base */
    node = left; 
  } 

  /* return pointer to factor */
  return node; 
} /* end of _unur_Factor() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_Bas_Exp (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  Base ::= Exponent                                                   */
     /*  Exponent ::= UnsignedConstant | Identifier | FunctDesignator        */ 
     /*               | '(' Expression ')'                                   */ 
     /*                                                                      */
     /*       UnsignedConstant                Identifier                     */
     /*      /                \      or     /          \       or            */
     /*  NULL                  NULL      NULL            NULL                */
     /*                                                                      */
     /*  FunctDesignator    or    Expression                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node; 
  char             *symb;
  int              token;

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* get next token */
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS)
    return _unur_fstr_error_parse(pdata,7,__LINE__); 

  /* constant or and variable */
  if( symbol[token].type==S_UCONST ||
      symbol[token].type==S_UIDENT ||
      symbol[token].type==S_SCONST ) {
    /* make a new end node */
    node = _unur_fstr_create_node(symb,0.,token,NULL,NULL); 
  }

  /* system function */
  else if( symbol[token].type == S_SFUNCT ) {
    /* set pointer to previous token again */
    --(pdata->tno);
    /* and get function */
    node = _unur_FunctDesignator(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(node); return NULL;
    }
  }
  
  else if( symb[0] == '(' ) {
    /* opening parenthesis  --> read expression in side parenthesis */
    node = _unur_Expression(pdata); 
    if (pdata->perrno) {
      _unur_fstr_free(node); return NULL;
    }
    
    /* next symbol must be closing parenthesis */
    if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
	 symb[0] != ')' )
      return _unur_fstr_error_parse(pdata,ERR_EXPECT_CLOSE_P,__LINE__);
  }
  
  else {
    /* unkown symbol */
    --(pdata->tno);
    return _unur_fstr_error_parse(pdata,ERR_UNKNOWN_SYMBOL,__LINE__);
  } 

  /* return pointer to base or exponent of an expression */
  return node; 
} /* end of _unur_Bas_Exp() */ 

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_FunctDesignator (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  FunctDesignator ::= FuncIdentifier '(' ActualParameterlist ')'      */
     /*                                                                      */
     /*       Identifier                                                     */
     /*      /          \                                                    */
     /*  NULL            ActualParameterlist                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *params; 
  char             *fsymb, *symb;
  int              funct, token;
  int              n_params; 

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* get function identifier for system function */
  if ( _unur_fstr_next_token(pdata,&funct,&fsymb) != UNUR_SUCCESS ||
       symbol[funct].type != S_SFUNCT )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_FUNCT,__LINE__);

  /* get number of parameters for this function */
  n_params = symbol[funct].info;

  /* read opening parenthesis '(' */
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       symb[0] != '(' )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_OPEN_P,__LINE__);

  /* read the parameter list */
  params = _unur_ActualParameterlist(pdata,n_params);
  if (pdata->perrno) {
    _unur_fstr_free(params); return NULL;
  }

  /* read closing parenthesis ')' */
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       symb[0] != ')' ) {
    _unur_fstr_free(params);
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_CLOSE_P,__LINE__);
  }

  /* store function in new node */
  node = _unur_fstr_create_node(fsymb,0.,funct,NULL,params); 

  /* return pointer to function */
  return node; 
} /* end of _unur_FunctDesignator() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_ActualParameterlist (struct parser_data *pdata, int n_params) 
     /*----------------------------------------------------------------------*/
     /*  ActualParameterlist ::= ActualParameter [ ',' ActualParameter ]     */
     /*                                                                      */
     /*                                           ','                        */
     /*  Expression  or:                         /   \                       */
     /*                     more expressions tree     Expression             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*   n_params  ... number of expected parameters in list                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;
  int              c_params;   /* counter for parameters */

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* read first parameter from string ...  */
  node = _unur_Expression(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(node); return NULL;
  }

  /* .. and set counter for parameters to 1 */
  c_params = 1; 

  /* scan token list while we find a list separator `,' */
  while ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
	  symb[0] == ',' ) {

    /* update counter for parameters */
    c_params++; 
    if (c_params > n_params) {
      _unur_fstr_free(node);
      return _unur_fstr_error_parse(pdata,ERR_INVALID_N_PARAMS,__LINE__);
    }

    /* old node becomes left node of `,' node */
    left = node; 

    /* make node for next variable (becomes right node) */
    right = _unur_Expression(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(left); _unur_fstr_free(right);
      return NULL;
    }

    /* make node for `,' separator */
    node = _unur_fstr_create_node(",",0.,s_comma,left,right); 
  }

  /* set token pointer to first element after parameter list */
  --(pdata->tno);

  /* check number of parameters */
  if (c_params < n_params) {
    _unur_fstr_free(node);
    return _unur_fstr_error_parse(pdata,ERR_INVALID_N_PARAMS,__LINE__);
  }

  /* return pointer to parameter list */
  return node; 
} /* end of _unur_ActualParameterlist() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Simplify tree                                                           **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_fstr_simplification (const char *symb, int token,
			   struct ftreenode *left, struct ftreenode *right) 
     /*----------------------------------------------------------------------*/
     /* Try to simpify nodes                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   symb   ... symbol string                                           */
     /*   token  ... token for which new node should be created              */
     /*   left   ... pointer to left node                                    */
     /*   right  ... pointer to right node                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   Pointer to simplified node                                         */
     /*   NULL if no simplication is possible                                */
     /*                                                                      */
     /* Side effects:                                                        */
     /*   'left' or 'right' may be freed when the node is simplified.        */
     /*   This cannot happen if NULL is returned.                            */
     /*----------------------------------------------------------------------*/
{
  /* some booleans */
  int l_const = left  && (left->type  == S_SCONST || left->type  == S_UCONST); 
  int r_const = right && (right->type == S_SCONST || right->type == S_UCONST); 
  int l_0     = (l_const && _unur_iszero(left->val));
  int l_1     = (l_const && _unur_isone(left->val));
  int r_0     = (r_const && _unur_iszero(right->val));
  int r_1     = (r_const && _unur_isone(right->val));
  int and;

  char s = symb[0];

  /*             Exp                    Exp
   *            /   \                  /   \ 
   *        NULL    ','       ==>     X     Y
   *               /   \ 
   *              X     Y
   */ 
  if ( left == NULL && right && right->symbol[0] == ',' ) {
    COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    right->token  = token;
    right->symbol = symbol[token].name; 
    right->type   = symbol[token].type;
    return right; 
   }

  /*          Operator            Operator
   *            /   \      or      /   \        ==>     Const (result of computation)
   *       Const     Const     NULL     Const
   */ 
  if ( (l_const || left==NULL) && r_const && s!=',') { 
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    /* compute new value */
    right->val   = ( (left) 
		     ? (*symbol[token].vcalc)(left->val,right->val)
		     : (*symbol[token].vcalc)(0.,right->val) );
    right->token = s_uconst;
    right->type  = S_UCONST;
    right->left  = NULL; 
    right->right = NULL;
    _unur_fstr_free(left);
    return right; 
  } 

  /*          '+'               '*'
   *         /   \      or     /   \        ==>     X
   *        0     X           1     X
   */ 
  if ( (l_0 && s=='+' ) || (l_1 && s=='*') ) { 
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(left);
    return right;
  } 

  /*          '+'               '-'  
   *         /   \      or     /   \      or
   *        X     0           X     0
   *
   *          '*'               '/'         
   *         /   \      or     /   \      or
   *        X     1           X     1       
   *
   *          '^'
   *         /   \                        ==>     X
   *        X     1
   */ 
  if ( (r_0 && (s=='+' || s=='-')) ||
       (r_1 && (s=='*' || s=='/' || s=='^')) ) {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(right);
    return left;
  }

  /*          '*'               '*'         
   *         /   \      or     /   \      or
   *        0     X           X     0       
   *
   *          '/'               '^'         
   *         /   \      or     /   \      or
   *        0     X           0     X       
   *
   *         "and"             "and"         
   *         /   \      or     /   \      ==>     0
   *        0     X           X     0       
   */
  and = (strcmp(symb,"and")==0);
  if ( l_0 && (s=='*' || s=='/' || s=='^' || and) ) {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(right);
    return left;
  }
  if (r_0 && (s=='*' || and ) ) {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(left);
    return right;
  }

  /*          '^'               '^'
   *         /   \      or     /   \        ==>     1
   *        X     0           1     X
   */ 
  if (r_0 && s=='^') {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(left);
    right->val = 1.;
    return right;
  }
  if (l_1 && s=='^') {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(right);
    return left;
  }

  /*              '/'              
   *          ___/   \___
   *         /           \ 
   *        X             X        ==>     1     
   *       / \           / \
   *   NULL   NULL   NULL   NULL
   */ 
  if ( ( symb[0] == '/' &&
	 left  && left->left==NULL  && left->right==NULL  && 
	 right && right->left==NULL && right->right==NULL &&
	 strcmp(left->symbol,right->symbol)== 0 ) ) {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(left);
    right->token = s_uconst;
    right->symbol= symbol[s_uconst].name; 
    right->val   = 1.;
    right->type  = S_UCONST;
    right->left  = NULL; 
    right->right = NULL;
    return right; 
  }

  /* no simplification */
  return NULL; 
} /* _unur_fstr_simplification() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_reorganize (struct ftreenode *node) 
     /*----------------------------------------------------------------------*/
     /* Try to reorganize tree at node                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   node   ... pointer to node in tree                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   1 if successful                                                    */
     /*   0 otherwise                                                        */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *left, *right, *tmp;
  char symb;
  int l_const, r_const;
  int rl_0, ll_0;

  /* check arguments */
  CHECK_NULL(node,0);  COOKIE_CHECK(node,CK_FSTR_TNODE,0);

  left  = node->left;
  right = node->right;
  symb = node->symbol[0];

  /* some booleans */
  l_const = left  && (left->type  == S_SCONST || left->type  == S_UCONST); 
  r_const = right && (right->type == S_SCONST || right->type == S_UCONST); 
  rl_0 = (right && right->left && right->left->type == S_UCONST && _unur_iszero(right->left->val));
  ll_0 = (left  && left->left  && left->left->type  == S_UCONST && _unur_iszero(left->left->val));

  /*          Operator            Operator
   *            /   \      or      /   \        ==>     Const (result of computation)
   *       Const     Const     NULL     Const
   */ 
  if ( (l_const || left==NULL) && r_const && symb!=',') { 
    /* compute new value */
    node->val   = ( (left) 
		    ? (*symbol[node->token].vcalc)(left->val,right->val)
		    : (*symbol[node->token].vcalc)(0.,right->val) );
    node->token = s_uconst;
    node->type  = S_UCONST;
    node->left  = NULL; 
    node->right = NULL;
    if (left)  free(left);
    if (right) free(right);
    return 1; 
  } 

  /*            '+'                  '-'  
   *           /   \                /   \ 
   *          X    '-'      ==>    X     Y
   *              /   \
   *             0     Y
   */ 
  if ( rl_0 && symb=='+' && right->symbol[0]=='-' ) {
    node->symbol = symbol[s_minus].name; 
    node->token  = s_minus; 
    node->type   = symbol[s_minus].type; 
    node->right  = right->right; 
    free(right->left);
    free(right);
    return 1;
  }

  /*            '-'                  '+'  
   *           /   \                /   \ 
   *          X    '-'      ==>    X     Y
   *              /   \
   *             0     Y
   */ 
  if ( rl_0 && symb=='-' && right->symbol[0]=='-' ) {
    node->symbol = symbol[s_plus].name; 
    node->token  = s_plus; 
    node->type   = symbol[s_plus].type; 
    node->right  = right->right; 
    free(right->left);
    free(right);
    return 1;
  }

  /*            '+'                  '-'  
   *           /   \                /   \ 
   *         '-'    Y     ==>      Y     X
   *        /   \
   *       0     X
   */ 
  if ( ll_0 && symb=='+' && left->symbol[0]=='-' ) {
    node->symbol = symbol[s_minus].name; 
    node->token  = s_minus; 
    node->type   = symbol[s_minus].type; 
    tmp = node->right;
    node->right  = left->right; 
    node->left   = tmp;
    free(left->left);
    free(left);
    return 1;
  }

  /*            '*'                  '-'  
   *           /   \                /   \ 
   *          X    '-'      ==>    0    '*'
   *              /   \                /   \
   *             0     Y              X     Y
   */ 
  if ( rl_0 && symb=='*' && right->symbol[0]=='-' ) {
    node->symbol = symbol[s_minus].name; 
    node->token  = s_minus; 
    node->type   = symbol[s_minus].type; 
    right->symbol= symbol[s_mul].name; 
    right->token = s_mul; 
    right->type  = symbol[s_mul].type; 
    tmp = left;
    node->left   = node->right->left;
    node->right  = tmp;
    return 1;
  }

  /* no reorganization */
  return 0;

} /* end of _unur_fstr_reorganize() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Auxilliary routines                                                     **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_fstr_next_token (struct parser_data *pdata, int *token, char **symb)
     /*----------------------------------------------------------------------*/
     /* Get next token from list.                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata  ... pointer to parser object                                */
     /*   token  ... to store token                                          */
     /*   symb   ... to store symbol for token                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(pdata,UNUR_ERR_NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);

  if (pdata->tno < pdata->n_tokens) {
    /* return token and increment scan position */
    *token = pdata->token[pdata->tno];
    *symb = pdata->tpos[pdata->tno];
    ++(pdata->tno);
    return UNUR_SUCCESS;
  }
  else {
    /* no more tokens */
    ++(pdata->tno);
    return UNUR_ERR_SILENT;
  }

} /* end of _unur_fstr_next_token() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_fstr_create_node (const char *symb, double val, int token, 
			struct ftreenode *left, struct ftreenode *right) 
     /*----------------------------------------------------------------------*/
     /* Create new node.                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   symb   ... symbol string                                           */
     /*   val    ... value of constant (only used if symb is NULL)           */
     /*   token  ... token for which node should be created                  */
     /*   left   ... pointer to left node                                    */
     /*   right  ... pointer to right node                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new node                                                */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node; 

  if ( symb && (node = _unur_fstr_simplification(symb,token,left,right)) ) {
    /* node has been simplified  -->  use left node, remove right node */
  }
  else {
    /* make new node */
    node = _unur_xmalloc(sizeof(struct ftreenode)); 
    COOKIE_SET(node,CK_FSTR_TNODE);

    node->symbol = symbol[token].name; 
    node->token  = token; 
    node->type   = symbol[token].type; 
    node->left   = left; 
    node->right  = right; 

    /* compute and/or store constants in val field */
    switch (symbol[token].type) {
    case S_UCONST:      /* user defined constant, i.e. a number */
      node->val = (symb) ? atof(symb) : val;  break;
    case S_SCONST:      /* system constant */
      node->val = symbol[token].val;  break;
    default:
      node->val = 0.;
    }
  } 

  /* it might be possible to reorganize the tree */
  _unur_fstr_reorganize(node); 
  
  /* return node */
  return node; 
  
} /* end of _unur_fstr_create_node() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Error messages                                                          **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_fstr_error_parse ( struct parser_data *pdata, int perrno, int line )
     /*----------------------------------------------------------------------*/
     /* Print error message when parsing function string                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata  ... pointer to parser object                                */
     /*   perrno ... error number                                            */
     /*   line   ... line (to be inserted by __LINE__)                       */
     /*                                                                      */
     /* return:                                                              */
     /*   NULL                                                               */
     /*----------------------------------------------------------------------*/
{ 
  int i;
  struct unur_string *reason;

  /* check arguments */
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);

  /* set parser error */
  if (!pdata->perrno) pdata->perrno = perrno;

  /* create string for reason of error */
  reason = _unur_string_new();
  _unur_string_append( reason, "%s: ", _unur_fstr_error_code(perrno) );
  for (i=0; i<pdata->tno-1; i++)
    _unur_string_append( reason, "%s ", pdata->tpos[i]);

  if (i<pdata->n_tokens)
    _unur_string_append( reason, " -->%s<--  ", pdata->tpos[i]);
  else
    _unur_string_append( reason, " <--  ");

  for (i++; i<pdata->n_tokens; i++)
    _unur_string_append( reason, "%s ",pdata->tpos[i]);
  
  /* report error */
  _unur_error_x( GENTYPE, __FILE__, line, "error", UNUR_ERR_FSTR_SYNTAX,reason->text);

  /* free working space */
  _unur_string_free( reason );

  return NULL; 

} /* _unur_fstr_error_parse() */ 

/*---------------------------------------------------------------------------*/

const char *
_unur_fstr_error_code ( int perrno )
     /*----------------------------------------------------------------------*/
     /* Print message for error number                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   perrno ... error number                                            */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to message string                                          */
     /*----------------------------------------------------------------------*/
{
  switch (perrno) {
  case ERR_UNFINISHED:
    return "incomplete. not all tokens parsed";
  case ERR_UNKNOWN_SYMBOL:
    return "unknown symbol in function string";
  case ERR_EXPECT_EQUAL:
    return "expected symbol: '='";
  case ERR_EXPECT_OPEN_P:
    return "expected symbol: '('";
  case ERR_EXPECT_CLOSE_P:
    return "expected symbol: ')'";
  case ERR_INVALID_N_PARAMS:
    return "invalid number of parameters for function";
  case ERR_EXPECT_FUNCT:
    return "function (name) expected";
  case ERR_EXPECT_VAR:
    return "user identifier (variable name) expected";
  case ERR_MISSING:
    return "more tokens expected";

  default:
    return "";
  }
} /* end of _unur_fstr_error_code() */

/*---------------------------------------------------------------------------*/
