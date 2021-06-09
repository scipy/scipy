/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser.c                                                *
 *                                                                           *
 *   Parser function string, evaluate function, print programm code.         *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given a string for a function.                                       *
 *         The string is parser.                                             *
 *         A tree representing the function term is generated.               *
 *         A tree for the derivative of the function is generated.           *
 *         The tree is used to evalute the corresponding function for an x.  *
 *         The source code for a program is produced.                        *
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
 *   [1] Aho, A. and Ullman, J. (1972): The Theory of Parsing, Translating   *
 *       and Compiling. Prentice-Hall.                                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Description:                                                             *
 *                                                                           *
 *  Function parser for use with UNU.RAN. For list of symbols and functions  *
 *  see functparser_symbols.h.                                               *
 *                                                                           *
 *  Important:                                                               *
 *  Blanks are ignored, e.g. "x>-1 and x < 1" is transformed to              *
 *  "x>-1andx<1" and thus the parser aborts at the unknwn symbol "andx".     *
 *  Use parenthesis to solve this problem: (x>-1) and (x<1).                 *
 *                                                                           *
 *  The first unknown symbol is treated as the variable.                     *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Pattern in Backus-Naur notation:                                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *                                                                           *
 *  FunctDefinition ::= DefFunctDesignator '=' Expression                    *
 *                                                                           *
 *                      '='                                                  *
 *                     /   \                                                 *
 *    DefFunctDesignator     Expression                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  DefFunctDesignator ::= Identifier '(' DefParameterlist ')'               *
 *                                                                           *
 *         Identifier                                                        *
 *        /          \                                                       *
 *    NULL            DefParameterlist                                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  DefParameterlist ::= '(' Identifier [ ',' Identifier ] ')'               *
 *                                                                           *
 *                                                                           *
 *         Identifier                                 ','                    *
 *        /          \      or:                      /   \                   *
 *    NULL            NULL      more identifiers tree     Identifier         *
 *                                                     /          \          *
 *                                                 NULL            NULL      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  Expression ::= SimpleExpression [ RelationOperator SimpleExpression ]    * 
 *                                                                           *
 *                                      RelationOperator                     *
 *   SimpleExpression  or:             /                \                    *
 *                     SimpleExpression                  SimpleExpression    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  SimpleExpression ::= STerm { AddingOperator Term }                       *
 *                                                                           *
 *                                  AddingOperator                           *
 *    STerm  or:                   /              \                          *
 *                  more terms tree                Term                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  STerm ::= [ '+' | '-' ] Term                                             *
 *                                                                           *
 *                          '-'                                              *
 *    Term  or:            /   \                                             *
 *                      '0'     Term                                         *
 *                     /   \                                                 *
 *                 NULL     NULL                                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  Term ::= Factor [ MultiplyingOperator Factor ]                           *
 *                                                                           *
 *                                     MultiplyingOperator                   *
 *    Factor  or:                     /                   \                  *
 *                   more factors tree                     Factor            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  Factor ::= Base [ '^' Exponent ]                                         *
 *                                                                           *
 *                            '^'                                            *
 *    Bas_Exp  or:           /   \                                           *
 *                    Bas_Exp     Bas_Exp                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  Base ::= Exponent                                                        *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  Exponent ::= UnsignedConstant | Identifier | FunctDesignator             * 
 *               | '(' Expression ')'                                        *
 *                                                                           *
 *         UnsignedConstant                Identifier                        *
 *        /                \      or     /          \       or               *
 *    NULL                  NULL      NULL            NULL                   *
 *                                                                           *
 *    FunctDesignator  or  Expression                                        *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  FunctDesignator ::= FuncIdentifier '(' ActualParameterlist ')'           *
 *                                                                           *
 *         Identifier                                                        *
 *        /          \                                                       *
 *    NULL            ActualParameterlist                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  ActualParameterlist ::= ActualParameter [ ',' ActualParameter ]          *
 *                                                                           *
 *                                             ','                           *
 *    Expression  or:                         /   \                          *
 *                       more expressions tree     Expression                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *  UnsignedConstant::= UnsignedInteger | UnsignedReal                       *
 *  UnsignedInteger ::= DigitSequence                                        *
 *  UnsignedReal    ::= UnsignedInteger ['.' DigitSequence] ['e' ScaleFactor]*
 *  DigitSequence   ::= Digit [ Digit [...] ]                                *
 *  Digit           ::= '0' | '1' | '2' | ... | '8' | '9'                    *
 *  ScaleFactor     ::= [Sign] DigitSequence                                 *
 *  Sign            ::= '+' | '-'                                            *
 *  Identifier      ::= Letter [ Letter | Digit [...] ]                      *
 *  Letter          ::= 'a' | 'b' | ... | 'z' | '_'                          *
 *  RelationOperator::= RelationChar [ RelationChar ]                        *
 *  RelationChar    ::= '<' | '=' | '>'                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Simplification rules:                                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *             Exp                    Exp                                    *
 *            /   \                  /   \                                   *
 *        NULL    ','       ==>     X     Y                                  *
 *               /   \                                                       *
 *              X     Y                                                      *
 *---------------------------------------------------------------------------*
 *          Operator            Operator                                     *
 *            /   \      or      /   \        ==>  Const (compute)           *
 *       Const     Const     NULL     Const                                  *
 *---------------------------------------------------------------------------*
 *          '+'               '*'                                            *
 *         /   \      or     /   \        ==>     X                          *
 *        0     X           1     X                                          *
 *---------------------------------------------------------------------------*
 *          '+'               '-'                                            *
 *         /   \      or     /   \      or                                   *
 *        X     0           X     0                                          *
 *                                                                           *
 *          '*'               '/'                                            *     
 *         /   \      or     /   \      or                                   *
 *        X     1           X     1                                          *
 *                                                                           *
 *          '^'                                                              *
 *         /   \                        ==>     X                            *
 *        X     1                                                            *
 *---------------------------------------------------------------------------*
 *          '*'               '*'                                            *
 *         /   \      or     /   \      or                                   *
 *        0     X           X     0                                          *
 *                                                                           *
 *          '/'               '^'                                            *
 *         /   \      or     /   \      or                                   *
 *        0     X           0     X                                          *
 *                                                                           *
 *         "and"             "and"                                           *
 *         /   \      or     /   \      ==>     0                            *
 *        0     X           X     0                                          *
 *---------------------------------------------------------------------------*
 *          '^'               '^'                                            *
 *         /   \      or     /   \        ==>     1                          *
 *        X     0           1     X                                          *
 *---------------------------------------------------------------------------*
 *              '/'                                                          *
 *          ___/   \___                                                      *
 *         /           \                                                     *
 *        X             X        ==>     1                                   *
 *       / \           / \                                                   *
 *   NULL   NULL   NULL   NULL                                               *
 *---------------------------------------------------------------------------*
 *            '+'                  '-'                                       *
 *           /   \                /   \                                      *
 *          X    '-'      ==>    X     Y                                     *
 *              /   \                                                        *
 *             0     Y                                                       *
 *---------------------------------------------------------------------------*
 *            '-'                  '+'                                       *
 *           /   \                /   \                                      *
 *          X    '-'      ==>    X     Y                                     *
 *              /   \                                                        *
 *             0     Y                                                       *
 *---------------------------------------------------------------------------*
 *            '+'                  '-'                                       *
 *           /   \                /   \                                      *
 *         '-'    Y     ==>      Y     X                                     *
 *        /   \                                                              *
 *       0     X                                                             *
 *---------------------------------------------------------------------------*
 *            '*'                  '-'                                       *
 *           /   \                /   \                                      *
 *          X    '-'      ==>    0    '*'                                    *
 *              /   \                /   \                                   *
 *             0     Y              X     Y                                  *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <ctype.h>
#include <unur_source.h>
#include "parser.h"
#include "parser_source.h"
#include "functparser_source.h"

/*---------------------------------------------------------------------------*/

#define GENTYPE "FSTRING"      /* (pseudo) type of generator                 */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Structures and definitions                                              **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Structure for storing data while tokizing and parsing the function string */

struct parser_data {
  char  *fstr;          /* pointer to function string                        */
  int   *token;         /* array to store marker for each token in string    */
  char  *tstr;          /* working array for tokenized string                */
  char **tpos;          /* array to store pointers to each token in string   */
  int    tno;           /* pointer to token in list                          */
  int    n_tokens;      /* total number of tokens                            */
  char  *variable_name; /* name of user defined identifier                   */
  char  *function_name; /* name of user defined function                     */
  int    scanpos;       /* pointer to location in string                     */
  int    lastpos;       /* position in string before scanning next token     */
  int    len_fstr;      /* length of function string                         */
  int    perrno;        /* error code                                        */

#ifdef UNUR_COOKIES
  unsigned cookie;      /* magic cookie                                      */
#endif
};

/*---------------------------------------------------------------------------*/
/* Error codes                                                               */

enum {
  SUCCESS = 0,          /* executed succesfully, no errors detected          */
  ERR_UNFINISHED,       /* incomplete. not all tokens parsed                 */
  ERR_UNKNOWN_SYMBOL,   /* unknown symbol in function string                 */
  ERR_EXPECT_EQUAL,     /* expected symbol: '='                              */
  ERR_EXPECT_OPEN_P,    /* expected symbol: '('                              */
  ERR_EXPECT_CLOSE_P,   /* expected symbol: ')'                              */
  ERR_INVALID_N_PARAMS, /* invalid number of parameters for function         */
  ERR_EXPECT_FUNCT,     /* function (name) expected                          */ 
  ERR_EXPECT_VAR,       /* user identifier (variable name) expected          */
  ERR_MISSING           /* more tokens expected                              */ 
};


/*****************************************************************************/
/** List of known symbols                                                   **/
/*****************************************************************************/

#define PARSER
#include "functparser_symbols.h"
#undef PARSER

/*****************************************************************************/
/** Prototypes                                                              **/
/*****************************************************************************/

/**-------------------------------------------------------------------------**/
/** functparser_init.c                                                      **/
/**   Init and destroy function tree.                                       **/
/**-------------------------------------------------------------------------**/

static struct ftreenode *_unur_fstr_2_tree (const char *functstr, int withDefFunct);
/*---------------------------------------------------------------------------*/
/* Compute funtion tree from string.                                         */
/*---------------------------------------------------------------------------*/


/**-------------------------------------------------------------------------**/
/** functparser_scanner.c                                                   **/
/**   Scan and tokenize function string                                     **/
/**-------------------------------------------------------------------------**/

/*---------------------------------------------------------------------------*/
/* Initialize and destroy parser                                             */
/*---------------------------------------------------------------------------*/

static struct parser_data *_unur_fstr_parser_init (const char *fstr);
/*---------------------------------------------------------------------------*/
/* Create and initialize parser object for given function string.            */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_symbols_init (void);
/*---------------------------------------------------------------------------*/
/* Prepare table of known symbols for usage.                                 */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_parser_free (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/* Destroy parser object.                                                    */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Scan function string                                                      */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_tokenize (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/* Tokenize function string.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_next_symbol (struct parser_data *pdata, char *symb);
/*---------------------------------------------------------------------------*/
/* Get next symbol in function string.                                       */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_find_symbol (const char *symb, int start, int end);
/*---------------------------------------------------------------------------*/
/* Find symbol in table between position (start+1) and (end-1).              */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_find_user_defined (struct parser_data *pdata, char *symb, int next_char);
/*---------------------------------------------------------------------------*/
/* Find user defined symbol.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_UnsignedConstant (struct parser_data *pdata, char *uc);
/*---------------------------------------------------------------------------*/
/* Get Unsigned Constant.                                                    */
/*  UnsignedConstant ::= UnsignedInteger | UnsignedReal                      */
/*  UnsignedInteger  ::= DigitSequence                                       */
/*  UnsignedReal ::= UnsignedInteger ['.' DigitSequence] ['e' ScaleFactor]   */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_DigitalSequence (struct parser_data *pdata, char *ds);
/*---------------------------------------------------------------------------*/
/* Get Digital Sequence.                                                     */
/*  DigitSequence ::= Digit [ Digit [...] ]                                  */
/*  Digit         ::= '0' | '1' | '2' | ... | '8' | '9'                      */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_ScaleFactor (struct parser_data *pdata, char *sf);
/*---------------------------------------------------------------------------*/
/* Get Scale Factor.                                                         */
/*  ScaleFactor ::= [Sign] DigitSequence                                     */
/*  Sign        ::= '+' | '-'                                                */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_Identifier (struct parser_data *pdata, char *id);
/*---------------------------------------------------------------------------*/
/* Get Identifier.                                                           */
/*  Identifier ::= Letter [ Letter | Digit [...] ]                           */
/*  Letter     ::= 'a' | 'b' | ... | 'z' | '_'                               */
/*  Digit      ::= '0' | '1' | ... | '9'                                     */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_RelationOperator (struct parser_data *pdata, char *ro);
/*---------------------------------------------------------------------------*/
/* Get Relation Operator.                                                    */
/*  RelationOperator ::= RelationChar [ RelationChar [...] ]                 */ 
/*  RelationChar     ::= '<' | '=' | '>'                                     */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_error_scan (const struct parser_data *pdata, const char *symb, int line);
/*---------------------------------------------------------------------------*/
/* Print error message when scanning function string.                        */
/*---------------------------------------------------------------------------*/


/**-------------------------------------------------------------------------**/
/** functparser_parser.c                                                    **/
/**   Parse tokenized function string and construct function tree           **/
/**-------------------------------------------------------------------------**/

/*---------------------------------------------------------------------------*/
/* Parser.                                                                   */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_FunctDefinition (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/* Get user defined function.                                                */
/*  FunctDefinition ::= DefFunctDesignator '=' Expression                    */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_DefFunctDesignator (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/* Definition of user defined function.                                      */
/*  DefFunctDesignator ::= Identifier '(' DefParameterlist ')'               */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_DefParameterlist (struct parser_data *pdata, int *n_params);
/*---------------------------------------------------------------------------*/
/* Parameter list for user defined function.                                 */
/*  DefParameterlist ::= '(' Identifier [ ',' Identifier ] ')'                */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_Expression (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  Expression ::= SimpleExpression [ RelationOperator SimpleExpression ]    */ 
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_SimpleExpression (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  SimpleExpression ::= STerm { AddingOperator Term }                       */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_STerm (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  STerm ::= [ '+' | '-' ] Term                                             */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_Term (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  Term ::= Factor [ MultiplyingOperator Factor ]                           */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_Factor (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  Factor ::= Base [ '^' Exponent ]                                         */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_Bas_Exp (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  Base ::= Exponent                                                        */
/*  Exponent ::= UnsignedConstant | Identifier | FunctDesignator |           */ 
/*               '(' Expression ')'                                          */ 
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_FunctDesignator (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  FunctDesignator ::= FuncIdentifier '(' ActualParameterlist ')'           */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_ActualParameterlist (struct parser_data *pdata, int n_params);
/*---------------------------------------------------------------------------*/
/*  ActualParameterlist ::= ActualParameter [ ',' ActualParameter ]          */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Simplify tree.                                                            */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_fstr_simplification (const char *symb, int token,
						    struct ftreenode *left,
						    struct ftreenode *right);
/*---------------------------------------------------------------------------*/
/* Try to simpify nodes.                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_reorganize (struct ftreenode *node);
/*---------------------------------------------------------------------------*/
/* Try to reorganize tree at node.                                           */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Auxilliary routines                                                       */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_next_token (struct parser_data *pdata, int *token, char **symb);
/*---------------------------------------------------------------------------*/
/* Get next token from list.                                                 */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_fstr_create_node (const char *symb, double val, int token,
						 struct ftreenode *left,
						 struct ftreenode *right);
/*---------------------------------------------------------------------------*/
/* Create new node.                                                          */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Error messages                                                            */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_fstr_error_parse ( struct parser_data *pdata, int perrno, int line );
/*---------------------------------------------------------------------------*/
/* Print error message when parsing function string.                         */
/*---------------------------------------------------------------------------*/

static const char *_unur_fstr_error_code ( int perrno );
/*---------------------------------------------------------------------------*/
/* Print message for error number.                                           */
/*---------------------------------------------------------------------------*/


/**-------------------------------------------------------------------------**/
/** functparser_eval.c                                                      **/
/**   Evaluate function tree for given argument x.                          **/
/**                                                                         **/
/** (see also functparser_symbols.h)                                        **/
/**-------------------------------------------------------------------------**/

static double _unur_fstr_eval_node (const struct ftreenode *node, double x);
/*---------------------------------------------------------------------------*/
/* Evaluate function tree starting from `node' at x                          */
/*---------------------------------------------------------------------------*/


/**-------------------------------------------------------------------------**/
/** functparser_deriv.c                                                     **/
/**   Compute function tree for derivative.                                 **/
/**                                                                         **/
/** (see also functparser_symbols.h)                                        **/
/**-------------------------------------------------------------------------**/

static void _unur_fstr_error_deriv (const struct ftreenode *node, int line);
/*---------------------------------------------------------------------------*/
/* Print error message for unknown derivative.                               */
/*---------------------------------------------------------------------------*/


/**-------------------------------------------------------------------------**/
/** functparser_stringgen.c                                                 **/
/**   Make string for function given by its tree.                           **/
/**-------------------------------------------------------------------------**/

static int _unur_fstr_node2string ( struct unur_string *output, const struct ftreenode *node,
				    const char *variable, const char *function, int spaces );
/*---------------------------------------------------------------------------*/
/* Produce string from function tree.                                        */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_print ( struct unur_string *output, const char *symb, double number );
/*---------------------------------------------------------------------------*/
/* Print string or number into output string.                                */
/* The number is only printed if symb is the NULL pointer.                   */
/*---------------------------------------------------------------------------*/


/**-------------------------------------------------------------------------**/
/** functparser_debug.c                                                     **/
/**   Debugging tools for function parser.                                  **/
/**-------------------------------------------------------------------------**/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_debug_input ( const char *fstr );
/*---------------------------------------------------------------------------*/
/* Print function string on output stream.                                   */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_debug_token ( const struct parser_data *pdata );
/*---------------------------------------------------------------------------*/
/* Print tokenized string on output stream.                                  */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_debug_tree ( const struct parser_data *pdata,
				    const struct ftreenode *root );
/*---------------------------------------------------------------------------*/
/* Print function tree.                                                      */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_debug_show_tree (const struct parser_data *pdata,
					const struct ftreenode *node,
					int level, int location);
/*---------------------------------------------------------------------------*/
/* Print function tree by recursion.                                         */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_debug_deriv (const struct ftreenode *funct,
				    const struct ftreenode *deriv);
/*---------------------------------------------------------------------------*/
/* Print function and its derivative.                                        */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Include sources                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#include "functparser_init.h"
#include "functparser_scanner.h"
#include "functparser_parser.h"
#include "functparser_eval.h"
#include "functparser_deriv.h"
#include "functparser_stringgen.h"
#include "functparser_debug.h"
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** End                                                                     **/
/*****************************************************************************/
