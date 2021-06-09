/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_scanner.c                                        *
 *                                                                           *
 *   Scan and tokenize function string                                       *
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
/** Initialize and destroy parser object                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct parser_data *
_unur_fstr_parser_init ( const char *fstr )
     /*----------------------------------------------------------------------*/
     /* Create and initialize parser object for given function string.       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   fstr ... string containing function definition                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to created and initialized parser object                   */
     /*                                                                      */
     /* error:                                                               */
     /*    return NULL                                                       */
     /*----------------------------------------------------------------------*/
{
  static int symbols_initialized = FALSE;
  struct parser_data *pdata;

  /* first we have to prepare the table of known symbols */
  if (symbols_initialized == FALSE) {
    _unur_fstr_symbols_init();
    symbols_initialized = TRUE;
  }

  /* allocate memory for parser object */
  pdata = _unur_xmalloc(sizeof(struct parser_data));
  COOKIE_SET(pdata,CK_FSTR_PDATA);

  /* make a working copy of the function string,                */
  /* remove all white spaces and convert to lower case letters. */
  pdata->fstr = _unur_parser_prepare_string(fstr);
  pdata->len_fstr = strlen(pdata->fstr);

  /* check whether string is empty */
  if (pdata->len_fstr <= 0) {
    /* empty, there is nothing to do */
    _unur_error(GENTYPE,UNUR_ERR_STR,"empty string"); 
    free (pdata->fstr);
    free(pdata);
    return NULL;
  }

  /* make arrays to store tokens in string */
  pdata->token = _unur_xmalloc( (pdata->len_fstr+1) * sizeof(int) );
  pdata->tpos  = _unur_xmalloc( (pdata->len_fstr+1) * sizeof(char *) );
  pdata->tstr  = _unur_xmalloc( (2*pdata->len_fstr+1) * sizeof(char) );

  /* initialize arrays */
  pdata->n_tokens = 0;
  memset(pdata->token,0  ,(size_t)pdata->len_fstr);
  memset(pdata->tpos, 0  ,(size_t)pdata->len_fstr);
  memset(pdata->tstr,'\0',(size_t)pdata->len_fstr);

  /* initialize for scanning */
  pdata->scanpos = 0;     /* scan position at beginning */
  pdata->lastpos = -1;
  pdata->perrno = 0;

  /* names of user defined symbols */
  pdata->variable_name = NULL;
  pdata->function_name = NULL;

  /* initialize data for parsing */
  pdata->tno = 0;

  /* return pointer to parser object */
  return pdata;

} /* end of _unur_fstr_parser_init() */

/*---------------------------------------------------------------------------*/

void
_unur_fstr_symbols_init (void)
     /*----------------------------------------------------------------------*/
     /* Prepare table of known symbols for usage.                            */
     /*                                                                      */
     /* parameters: none                                                     */
     /*----------------------------------------------------------------------*/
{
  int i;
  char *s;

  /* none of these positions found */
  _ros_start = 0;
  _nas_start = 0;
  _ans_start = 0;

  /* find marker in list */
  for (i=0; !_end; i++) {
    s = symbol[i].name;
    if (!_ros_start) {
      if ( strcmp(s,"_ROS") == 0) _ros_start = i;
      continue;
    }
    if (!_nas_start) {
      if ( strcmp(s,"_NAS") == 0) _nas_start = i;
      continue;
    }
    if (!_ans_start) {
      if ( strcmp(s,"_ANS") == 0) _ans_start = i;
      continue;
    }
    if (strcmp(s,"_END") == 0) _end = i;
  }

  /* end of region markers for list */
  _ros_end = _nas_start;
  _nas_end = _ans_start;
  _ans_end = _end;

  /* find location of special symbols */
  s_comma = _unur_fstr_find_symbol(",",_nas_start,_nas_end);
  s_minus = _unur_fstr_find_symbol("-",_nas_start,_nas_end);
  s_plus  = _unur_fstr_find_symbol("+",_nas_start,_nas_end);
  s_mul   = _unur_fstr_find_symbol("*",_nas_start,_nas_end);
  s_div   = _unur_fstr_find_symbol("/",_nas_start,_nas_end);
  s_power = _unur_fstr_find_symbol("^",_nas_start,_nas_end);

} /* end of _unur_fstr_symbols_init() */

/*---------------------------------------------------------------------------*/

void
_unur_fstr_parser_free ( struct parser_data *pdata )
     /*----------------------------------------------------------------------*/
     /* Destroy parser object.                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*----------------------------------------------------------------------*/
{
  if (pdata) {
    COOKIE_CHECK(pdata,CK_FSTR_PDATA,RETURN_VOID);
    free(pdata->fstr);
    free(pdata->token);
    free(pdata->tpos);
    free(pdata->tstr);
    free(pdata);
  }

} /* end of _unur_fstr_parser_free() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Tokenize function string                                                **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_fstr_tokenize (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /* Tokenize function string                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   error code of called subroutines                                   */
     /*----------------------------------------------------------------------*/
{
  int token;

  int n_token = 0;                    /* counter for token      */
  char *symb = pdata->tstr;           /* array to store token   */

  /* check arguments */
  CHECK_NULL(pdata,UNUR_ERR_COOKIE);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);

  /* locate token in function string and copy into token string */
  while ((token = _unur_fstr_next_symbol(pdata,symb)) != S_NOSYMBOL) {
    pdata->token[n_token] = token;    /* marker for token */
    pdata->tpos[n_token] = symb;  /* name of token    */
    n_token++;                        /* increment counter for tokens */
    symb += pdata->scanpos - pdata->lastpos + 1;  /* skip to unused space in string */ 
  }

  /* store total number of tokens */
  pdata->n_tokens = n_token;

  /* set token pointer to first token */
  pdata->tno = 0;

  /* return error code */
  return pdata->perrno;

} /* end of _unur_fstr_tokenize() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_next_symbol (struct parser_data *pdata, char *symb)
     /*----------------------------------------------------------------------*/
     /* Get next symbol in function string.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   symb  ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* return:                                                              */
     /*   error code of called subroutines                                   */
     /*----------------------------------------------------------------------*/
{
  int token;
  int errcode = 0;
  char c;

  /* check arguments */
  CHECK_NULL(pdata,S_NOSYMBOL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,S_NOSYMBOL);

  /* store position of pointer */
  pdata->lastpos = pdata->scanpos;

  if (pdata->scanpos >= pdata->len_fstr)
    /* end of string */
    return S_NOSYMBOL;
  
  /* next character in string */
  c = pdata->fstr[pdata->scanpos];

  /* new get next symbol in list */
  if ( (c >= '0' && c <= '9') || c == '.') {
    /* Unsigned Constant */
    _unur_fstr_UnsignedConstant(pdata,symb);
    token = s_uconst;
  }

  else if (c >=  'a' && c <= 'z') {
    /* Identifier */
    _unur_fstr_Identifier(pdata,symb);

    if ( ( (token = _unur_fstr_find_symbol(symb,_ans_start,_ans_end)) == 0 ) && 
	 ( (token = _unur_fstr_find_user_defined(pdata,symb,pdata->fstr[pdata->scanpos])) <= 0 ) )
      errcode = ERR_UNKNOWN_SYMBOL;
  }

  else if ( c == '<' || c == '>' || c == '=' || c == '!' ) {
    /* Relation Operator */
    _unur_fstr_RelationOperator(pdata,symb);

    if ((token = _unur_fstr_find_symbol(symb,_ros_start,_ros_end)) <= 0 )
      errcode = ERR_UNKNOWN_SYMBOL;
  }

  else {
    symb[0] = c; symb[1] = '\0';           /* all other charactors */
    (pdata->scanpos)++;
    if ((token = _unur_fstr_find_symbol(symb,_nas_start,_nas_end)) <= 0 )
      errcode = ERR_UNKNOWN_SYMBOL;
  }

  /* set errorcode */
  pdata->perrno = errcode;

  if (errcode) {
    _unur_fstr_error_scan (pdata,symb,__LINE__);
  }

  return (errcode) ? S_NOSYMBOL : token;

} /* end of _unur_fstr_next_symbol() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_find_symbol (const char *symb, int start, int end)
     /*----------------------------------------------------------------------*/
     /* find symbol in table between position (start+1) and (end-1)          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   symb  ... symbol to look for                                       */
     /*   start ... starting point for searching                             */
     /*   end   ... endpoint for searching                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   location in table                                                  */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* search for symbol in table */
  for (i = start + 1; i < end; i++)
    if (strcmp(symb,symbol[i].name) == 0) break;

  /* return location if symbol found and 0 otherwise */
  return ((i < end ) ? i : 0);

} /* end of _unur_fstr_find_symbol() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_find_user_defined (struct parser_data *pdata, char *symb, int next_char)
     /*----------------------------------------------------------------------*/
     /* Find user defined symbol.                                            */
     /* If there are no user defined symbols yet, store it.                  */
     /* If there are already user defined symbols, and this one is new,      */
     /* make return an errorcode.                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*   symb      ... symbol to look for                                   */
     /*   next_char ... character following given symbol                     */
     /*                                                                      */
     /* return:                                                              */
     /*   location in table                                                  */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(pdata,0);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,0);

  /* we use next_char to distinguish between variables and functions */

  if (next_char == '(') {
    /* symbol is user defined function */
    if (pdata->function_name == NULL) {
      /* new function --> store name */
      pdata->function_name = symb;
      /* return marker for user defined identifier */
      return s_ufunct;
    }
    else
      /* the identifier name must match with variable name */
      return (strcmp(pdata->function_name,symb) == 0) ? s_ufunct : 0;
  }
 
  else {
    /* symbol is user defined identifier (i.e. a variable) */
    if (pdata->variable_name == NULL) {
      /* new variable --> store name */
      pdata->variable_name = symb;
      /* return marker for user defined identifier */
      return s_uident;
    }
    else
      /* the identifier name must match with variable name */
      return (strcmp(pdata->variable_name,symb) == 0) ? s_uident : 0;
  }

} /* end of _unur_fstr_find_user_defined() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Scan function string                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_fstr_UnsignedConstant (struct parser_data *pdata, char *uc)
     /*----------------------------------------------------------------------*/
     /* Get Unsigned Constant                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   uc    ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* Syntax:                                                              */
     /*   UnsignedConstant ::= UnsignedInteger | UnsignedReal                */
     /*   UnsignedInteger  ::= DigitSequence                                 */
     /*   UnsignedReal ::= UnsignedInteger ['.' DigitSequence] ['e' ScaleFactor] */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* store scan position */
  int startpos = pdata->scanpos;

  /* check arguments */
  CHECK_NULL(pdata,UNUR_ERR_NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);

  /* copy digit sequence into uc */
  _unur_fstr_DigitalSequence(pdata,uc);

  if( pdata->fstr[pdata->scanpos] == '.' ) {
    /* decimal point  --> copy point and following digits */
    *(uc + pdata->scanpos - startpos) = '.';
    (pdata->scanpos)++;
    _unur_fstr_DigitalSequence(pdata, uc + pdata->scanpos - startpos);
  }

  if( pdata->fstr[pdata->scanpos] == 'e' ) {
    /* exponent --> copy indicator 'E' and following [sign and] digits */
    *(uc + pdata->scanpos - startpos) = 'e';
    (pdata->scanpos)++;
    _unur_fstr_ScaleFactor(pdata, uc + pdata->scanpos - startpos);
  }

  return UNUR_SUCCESS;
} /* end of _unur_fstr_UnsignedConstant() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_DigitalSequence (struct parser_data *pdata, char *ds)
     /*----------------------------------------------------------------------*/
     /* Get Digital Sequence                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   ds    ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* Syntax:                                                              */
     /*   DigitSequence ::= Digit [ Digit [...] ]                            */
     /*   Digit         ::= '0' | '1' | '2' | ... | '8' | '9'                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(pdata,UNUR_ERR_NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);

  /* copy digit */
  while ( (*ds = pdata->fstr[pdata->scanpos]) >= '0' && *ds <= '9' ) {
     ds++;
     (pdata->scanpos)++;
  }
  /* terminate string */
  *ds = '\0';

  return UNUR_SUCCESS;
} /* end of _unur_fstr_DigitalSequence() */

/*---------------------------------------------------------------------------*/

int 
_unur_fstr_ScaleFactor (struct parser_data *pdata, char *sf)
     /*----------------------------------------------------------------------*/
     /* Get Scale Factor                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   sf    ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* Syntax:                                                              */
     /*   ScaleFactor ::= [Sign] DigitSequence                               */
     /*   Sign        ::= '+' | '-'                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(pdata,UNUR_ERR_NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);

  /* copy sign */
  if ( (sf[0] = pdata->fstr[pdata->scanpos]) == '+' || sf[0] == '-' ) {
     sf++;
     (pdata->scanpos)++;
  }
  /* copy digital sequence (after sign) */
  _unur_fstr_DigitalSequence(pdata,sf);

  return UNUR_SUCCESS;
} /* _unur_fstr_ScaleFactor() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_Identifier (struct parser_data *pdata, char *id)
     /*----------------------------------------------------------------------*/
     /* Get Identifier                                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   id    ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* Syntax:                                                              */
     /*   Identifier ::= Letter [ Letter | Digit [...] ]                     */
     /*   Letter     ::= 'a' | 'b' | ... | 'z' | '_'                         */
     /*   Digit      ::= '0' | '1' | ... | '9'                               */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(pdata,UNUR_ERR_NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);

  /* copy word */
  while ( ((*id = pdata->fstr[pdata->scanpos]) >= 'a' && *id <= 'z')
	  || *id == '_' 
	  || ( *id >= '0' && *id <= '9')) {
    id++;
    (pdata->scanpos)++;
  }
  /* terminate string */
  *id = '\0';

  return UNUR_SUCCESS;
} /* end of _unur_fstr_Identifier() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_RelationOperator (struct parser_data *pdata, char *ro)
     /*----------------------------------------------------------------------*/
     /* Get Relation Operator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   ro    ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* Syntax:                                                              */
     /*   RelationOperator ::= RelationChar [ RelationChar ]                 */ 
     /*   RelationChar     ::= '<' | '>' | '=' | '!'                         */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(pdata,1);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,1);

  /* copy relation operator */
  while ((*ro = pdata->fstr[pdata->scanpos]) == '<' || *ro == '>' || *ro == '=' || *ro == '!' ) {
    ro++;
    (pdata->scanpos)++;
  }
  /* terminate string */
  *ro = '\0';

  return UNUR_SUCCESS;
} /* _unur_fstr_RelationOperator() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Error messages                                                          **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

void
_unur_fstr_error_scan (const struct parser_data *pdata, const char *symb, int line)
     /*----------------------------------------------------------------------*/
     /* Print error message when scanning function string                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   symb  ... pointer to unknown symbol                                */  
     /*----------------------------------------------------------------------*/
{
  struct unur_string *reason;
  char *c;
  
  /* check arguments */
  CHECK_NULL(pdata,RETURN_VOID);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,RETURN_VOID);

  /* create string for reason of error */
  reason = _unur_string_new();

  /* print unknown symbol */
  _unur_string_append( reason, "unknown symbol '%s': ", symb );

  /* print scanned part of function string to error stream */
  for (c=pdata->fstr; c < pdata->fstr+pdata->lastpos; c++) 
    _unur_string_append( reason, "%c", *c );

  /* print remaining part of function string including unknown symbol */
  _unur_string_append( reason, "    --> %s", pdata->fstr + pdata->lastpos);

  /* report error */
  _unur_error_x( GENTYPE, __FILE__, line, "error", UNUR_ERR_FSTR_SYNTAX,reason->text);

  /* free working space */
  _unur_string_free( reason );

} /* end of _unur_fstr_error_scan() */

/*---------------------------------------------------------------------------*/
