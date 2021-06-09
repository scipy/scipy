/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_symbols.c                                        *
 *                                                                           *
 *   Table of known symbols for function parser                              *
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
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/*
   #define one of the two marcos 'PARSER' or 'CODEGEN' before #including 
   this header file:

   'PARSER'  ... for parsing and evalution of function string
   'CODEGEN' ... for the code generating routines

   These two macros should be #undef-ined after this file has been #included!
*/

/* #define PARSER */
/* #define CODEGEN */

#if defined(PARSER) && defined(CODEGEN)
#  error \
   You have to #define either 'PARSER' or 'CODEGEN' before including this header file!
#endif

#if !defined(PARSER) && !defined(CODEGEN)
#  error \
   You have to #define either 'PARSER' or 'CODEGEN' before including this header file!
#endif

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Routines for evaluating nodes of the function tree                      **/
/*****************************************************************************/

#ifdef PARSER
static double v_dummy  (double l, double r);
static double v_less   (double l, double r);
static double v_equal  (double l, double r);
static double v_greater(double l, double r);
static double v_less_or(double l, double r);
static double v_unequal(double l, double r);
static double v_grtr_or(double l, double r);
static double v_plus   (double l, double r); 
static double v_minus  (double l, double r); 
static double v_mul    (double l, double r);
static double v_div    (double l, double r);
static double v_mod    (double l, double r);
static double v_power  (double l, double r);
static double v_const  (double l, double r);
static double v_exp    (double l, double r);
static double v_log    (double l, double r);
static double v_sin    (double l, double r); 
static double v_cos    (double l, double r);
static double v_tan    (double l, double r);
static double v_sec    (double l, double r); 
static double v_sqrt   (double l, double r);
static double v_abs    (double l, double r);
static double v_sgn    (double l, double r);
#endif

/*****************************************************************************/
/** Routines for computing derivatives                                      **/
/*****************************************************************************/

#ifdef PARSER
static struct ftreenode *d_error ( const struct ftreenode *node, int *error );
static struct ftreenode *d_const ( const struct ftreenode *node, int *error );
static struct ftreenode *d_var   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_add   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_mul   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_div   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_power ( const struct ftreenode *node, int *error );
static struct ftreenode *d_exp   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_log   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_sin   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_cos   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_tan   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_sec   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_sqrt  ( const struct ftreenode *node, int *error );
static struct ftreenode *d_abs   ( const struct ftreenode *node, int *error );
#endif


/*****************************************************************************/
/** Routines for printing C code                                            **/
/*****************************************************************************/
#ifdef CODEGEN

/* error */
static unsigned C_error  ( struct unur_string *output, const struct ftreenode *node, const char *variable );

/* constants and variables */
static unsigned C_const  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_var    ( struct unur_string *output, const struct ftreenode *node, const char *variable );

/* prefix operators (functions) (eg. exp(x)) */
static unsigned C_prefix_generic ( struct unur_string *output, const struct ftreenode *node, 
				   const char *variable, const char *symbol );
/* operators where C routine is the same as parser name */
static unsigned C_prefix ( struct unur_string *output, const struct ftreenode *node, const char *variable );
/* special operators */
static unsigned C_power  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_sec    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_abs    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_sgn    ( struct unur_string *output, const struct ftreenode *node, const char *variable );

/* infix (binary) operators (eg. x + y) */
static unsigned C_infix_generic  ( struct unur_string *output, const struct ftreenode *node,
				   const char *variable, const char *symbol );
/* operators where C routine is the same as parser name */
static unsigned C_infix  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
/* special operators */
static unsigned C_equal  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_unequal( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_minus  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_mod    ( struct unur_string *output, const struct ftreenode *node, const char *variable );

#endif

/*****************************************************************************/
/** Routines for printing FORTRAN code                                      **/
/*****************************************************************************/
#ifdef CODEGEN

/* error */
static unsigned F_error  ( struct unur_string *output, const struct ftreenode *node, const char *variable );

/* constants and variables */
static unsigned F_const  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_var    ( struct unur_string *output, const struct ftreenode *node, const char *variable );

/* prefix operators (functions) (eg. exp(x)) */
static unsigned F_prefix_generic ( struct unur_string *output, const struct ftreenode *node, 
				   const char *variable, const char *symbol );
/* operators where C routine is the same as parser name */
static unsigned F_prefix ( struct unur_string *output, const struct ftreenode *node, const char *variable );
/* relation operators */
static unsigned F_lt     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_le     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_gt     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_ge     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_eq     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_ne     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
/* special operators */
static unsigned F_sec    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_sgn    ( struct unur_string *output, const struct ftreenode *node, const char *variable );

/* infix (binary) operators (eg. x + y) */
static unsigned F_infix_generic  ( struct unur_string *output, const struct ftreenode *node,
				   const char *variable, const char *symbol );
/* operators where FORTRAN routine is the same as parser name */
static unsigned F_infix  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
/* special operators */
static unsigned F_minus  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_power  ( struct unur_string *output, const struct ftreenode *node, const char *variable );

/* static unsigned F_equal  ( struct unur_string *output, const struct ftreenode *node, const char *variable ); */
/* static unsigned F_unequal( struct unur_string *output, const struct ftreenode *node, const char *variable ); */

#endif

/*****************************************************************************/
/** Routines for printing JAVA code                                         **/
/*****************************************************************************/
#ifdef CODEGEN

/* error */
static unsigned J_error  ( struct unur_string *output, const struct ftreenode *node, const char *variable );

/* constants and variables */
static unsigned J_const  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_var    ( struct unur_string *output, const struct ftreenode *node, const char *variable );

/* prefix operators (functions) (eg. exp(x)) */
static unsigned J_prefix_generic ( struct unur_string *output, const struct ftreenode *node, 
				   const char *variable, const char *symbol );
/* operators where C routine is the same as parser name */
static unsigned J_prefix ( struct unur_string *output, const struct ftreenode *node, const char *variable );
/* relation operators */
static unsigned J_lt     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_le     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_gt     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_ge     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_eq     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_ne     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
/* special operators */
static unsigned J_power  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_sec    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_sgn    ( struct unur_string *output, const struct ftreenode *node, const char *variable );

/* infix (binary) operators (eg. x + y) */
static unsigned J_infix_generic  ( struct unur_string *output, const struct ftreenode *node,
				   const char *variable, const char *symbol );
/* operators where C routine is the same as parser name */
static unsigned J_infix  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
/* special operators */
static unsigned J_minus  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_mod    ( struct unur_string *output, const struct ftreenode *node, const char *variable );

#endif

/*****************************************************************************/
/** List of known symbols                                                   **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Symbol types for tokens                                                   */

enum {
  S_NOSYMBOL = 0,    /* token is not a symbol / is unknown                   */
  S_SFUNCT,          /* system function                                      */
  S_SCONST,          /* system constant                                      */
  S_UIDENT,          /* user defined identifier (variable)                   */
  S_UFUNCT,          /* user defined function                                */
  S_UCONST,          /* user defined constant                                */
  S_REL_OP,          /* relation operator                                    */
  S_ADD_OP,          /* addition operator                                    */
  S_MUL_OP,          /* multiplication operator                              */
  S_HPR_OP,          /* higher priority operator                             */
  S_OTHERS           /* other operator                                       */
};

/*  IMPORTANT: symbol types S_ADD_OP, S_MUL_OP and S_HPR_OP indicate         */
/*  priorities of the corresponding operators.                               */

/*---------------------------------------------------------------------------*/
/* structure for storing known symbols                                       */

#define SYMBLENGTH 10           /* Maximum length for symbols                */

/* structure for storing known symbols                                       */
struct symbols { 
  char   name[SYMBLENGTH];       /* name of symbols (e.g. "sin")             */ 

#ifdef PARSER
  int    type;                   /* type of symbol */
  int    info;                   /* priority or 
				    number of argument for system function   */
  double val;                    /* value of constant (0. otherwise)         */

  double (*vcalc)(double l, double r);
                                 /* function for computing value of node     */
  struct ftreenode *(*dcalc)(const struct ftreenode *node, int *error); 
                                 /* function for computing derivate          */
#endif

#ifdef CODEGEN
  unsigned (*node2C)(struct unur_string *output, const struct ftreenode *node, const char *variable );
                                 /* function for printing C code             */
  unsigned (*node2F)(struct unur_string *output, const struct ftreenode *node, const char *variable );
                                 /* function for printing FORTRAN code       */
  unsigned (*node2J)(struct unur_string *output, const struct ftreenode *node, const char *variable );
                                 /* function for printing JAVA code          */
#endif
};

/*---------------------------------------------------------------------------*/
/* List of known symbols                                                     */

/*  IMPORTANT: Symbols that do NOT start with a digit, a letter or `.'       */
/*  must not have more than ONE single character.                            */
/*  The only exception are relation operators.                               */
/*                                                                           */
/*  IMPORTANT: symbol types S_ADD_OP, S_MUL_OP and S_HPR_OP indicate         */
/*  priorities of the corresponding operators.                               */
/*                                                                           */
/*  In the current implementation system function must not have more than    */
/*  two arguments.                                                           */
/*                                                                           */
/*  Constant should have highest priority to avoid enclosing parenthesis     */
/*  in string output.                                                        */
/*                                                                           */
/*  The first entry in the list must be empty!                               */
/*                                                                           */
/*  There must be FOUR markers in the list (in this ordering!):              */
/*     _ROS   ... start of relation operators                                */
/*     _NAS   ... start of non-alphanumeric symbols (single characters!)     */
/*     _ANS   ... start of alphanumeric symbols (constants and functions)    */
/*     _END   ... end of list                                                */
/*                                                                           */

/* entries needed for parsing and evaluation of terms */
#ifdef PARSER
#  define S(name,type,priority,value,funct,deriv) \
   name,type,priority,value,funct,deriv
#else
#  define S(name,type,priority,value,funct,deriv) \
   name,
#endif

/* code generating functings */
#ifdef CODEGEN
#  define CG(cc,fortran,java)  cc,fortran,java
#else
#  define CG(cc,fortran,java)
#endif

static struct symbols symbol[] = {   
  /*  symbol,	         priority, evaluation routine,     C function,          Java function */
  /*          type,          value,           derivative              F function,             */

  /* void */
  {S (""    , S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},

  /* user defined symbols: */                        
  /*    constant           */                     
  {S ("UCONST",S_UCONST , 9, 0.0 , v_const  , d_const) CG (C_const  , F_const , J_const  )},
  /*    function           */
  {S ("UFUNCT",S_UFUNCT , 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},
  /*    variable           */
  {S ("VAR" , S_UIDENT  , 9, 0.0 , v_dummy  , d_var  ) CG (C_var    , F_var   , J_var    )},

  /* marker for relation operators */
  {S ("_ROS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},

  /* relation operators  */
  {S ("<"   , S_REL_OP  , 1, 0.0 , v_less   , d_const) CG (C_infix  , F_lt    , J_lt     )},
  {S ("="   , S_REL_OP  , 1, 0.0 , v_equal  , d_const) CG (C_equal  , F_eq    , J_eq     )},
  {S ("=="  , S_REL_OP  , 1, 0.0 , v_equal  , d_const) CG (C_equal  , F_eq    , J_eq     )},
  {S (">"   , S_REL_OP  , 1, 0.0 , v_greater, d_const) CG (C_infix  , F_gt    , J_gt     )},
  {S ("<="  , S_REL_OP  , 1, 0.0 , v_less_or, d_const) CG (C_infix  , F_le    , J_le     )},
  {S ("<>"  , S_REL_OP  , 1, 0.0 , v_unequal, d_const) CG (C_unequal, F_ne    , J_ne     )},
  {S ("!="  , S_REL_OP  , 1, 0.0 , v_unequal, d_const) CG (C_unequal, F_ne    , J_ne     )},
  {S (">="  , S_REL_OP  , 1, 0.0 , v_grtr_or, d_const) CG (C_infix  , F_ge    , J_ge     )},

  /* marker for non-alphanumeric symbols */
  {S ("_NAS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},

  /* special symbols */
  {S ("("   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},
  {S (")"   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},
  {S (","   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},

  /* arithmetic operators */
  {S ("+"   , S_ADD_OP  , 2, 0.0 , v_plus   , d_add  ) CG (C_infix  , F_infix , J_infix  )},
  {S ("-"   , S_ADD_OP  , 2, 0.0 , v_minus  , d_add  ) CG (C_minus  , F_minus , J_minus  )},
  {S ("*"   , S_MUL_OP  , 4, 0.0 , v_mul    , d_mul  ) CG (C_infix  , F_infix , J_infix  )},
  {S ("/"   , S_MUL_OP  , 4, 0.0 , v_div    , d_div  ) CG (C_infix  , F_infix , J_infix  )},
  {S ("^"   , S_HPR_OP  , 5, 0.0 , v_power  , d_power) CG (C_power  , F_power , J_power  )},

  /* marker for alphanumeric symbols */
  {S ("_ANS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},

  /* logical operators: removed */

  /* system constants */
  {S ("pi"  , S_SCONST  , 9, M_PI, v_const  , d_const) CG (C_const  , F_const , J_const  )},
  {S ("e"   , S_SCONST  , 9, M_E , v_const  , d_const) CG (C_const  , F_const , J_const  )},

  /* system functions */
  {S ("mod" , S_SFUNCT  , 2, 0.0 , v_mod    , d_const) CG (C_mod    , F_prefix, J_mod    )},
  {S ("exp" , S_SFUNCT  , 1, 0.0 , v_exp    , d_exp  ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("log" , S_SFUNCT  , 1, 0.0 , v_log    , d_log  ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("sin" , S_SFUNCT  , 1, 0.0 , v_sin    , d_sin  ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("cos" , S_SFUNCT  , 1, 0.0 , v_cos    , d_cos  ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("tan" , S_SFUNCT  , 1, 0.0 , v_tan    , d_tan  ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("sec" , S_SFUNCT  , 1, 0.0 , v_sec    , d_sec  ) CG (C_sec    , F_sec   , J_sec    )},
  {S ("sqrt", S_SFUNCT  , 1, 0.0 , v_sqrt   , d_sqrt ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("abs" , S_SFUNCT  , 1, 0.0 , v_abs    , d_abs  ) CG (C_abs    , F_prefix, J_prefix )},
  {S ("sgn" , S_SFUNCT  , 1, 0.0 , v_sgn    , d_const) CG (C_sgn    , F_sgn   , J_sgn    )},

  /* marker for end-of-table */
  {S ("_END", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error )},
};

#undef S
#undef GC

/*---------------------------------------------------------------------------*/

#ifdef PARSER

/*  location of special symbols in table                                     */
/** DO NOT CHANGE THEIR POSITION IN TABLE **/
static int s_uconst = 1;      /* user defined constant                       */
static int s_ufunct = 2;      /* user defined function                       */
static int s_uident = 3;      /* user defined variable                       */

/* location of special symbols */
static int s_comma, s_minus, s_plus, s_mul, s_div, s_power;

/* Marker for different regions in symbol table                              */
static int _ros_start, _ros_end;    /* relation symbols                      */
static int _nas_start, _nas_end;    /* non-alphanumeric symbols              */
static int _ans_start, _ans_end;    /* alphanumeric symbols                  */
static int _end;                    /* end of list                           */

#endif

/*---------------------------------------------------------------------------*/
/* #undef PARSER */
/* #undef CODEGEN */
/*---------------------------------------------------------------------------*/
