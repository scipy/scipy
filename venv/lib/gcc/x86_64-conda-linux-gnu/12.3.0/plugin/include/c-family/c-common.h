/* Definitions for c-common.cc.
   Copyright (C) 1987-2022 Free Software Foundation, Inc.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3, or (at your option) any later
version.

GCC is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with GCC; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.  */

#ifndef GCC_C_COMMON_H
#define GCC_C_COMMON_H

#include "splay-tree.h"
#include "cpplib.h"
#include "alias.h"
#include "tree.h"
#include "fold-const.h"
#include "wide-int-bitmask.h"

/* In order for the format checking to accept the C frontend
   diagnostic framework extensions, you must include this file before
   diagnostic-core.h, not after.  The C front end formats are a subset of those
   for C++, so they are the appropriate set to use in common code;
   cp-tree.h overrides this for C++.  */
#if defined(GCC_DIAGNOSTIC_CORE_H)
#error \
In order for the format checking to accept the C front end diagnostic \
framework extensions, you must include this file before diagnostic-core.h \
never after.
#endif
#ifndef GCC_DIAG_STYLE
#define GCC_DIAG_STYLE __gcc_cdiag__
#endif
#include "diagnostic-core.h"

/* Usage of TREE_LANG_FLAG_?:
   0: IDENTIFIER_MARKED (used by search routines).
      C_MAYBE_CONST_EXPR_INT_OPERANDS (in C_MAYBE_CONST_EXPR, for C)
   1: C_DECLARED_LABEL_FLAG (in LABEL_DECL)
      STATEMENT_LIST_STMT_EXPR (in STATEMENT_LIST)
      C_MAYBE_CONST_EXPR_NON_CONST (in C_MAYBE_CONST_EXPR, for C)
   2: unused
   3: STATEMENT_LIST_HAS_LABEL (in STATEMENT_LIST)
   4: unused
*/

/* Reserved identifiers.  This is the union of all the keywords for C,
   C++, and Objective-C.  All the type modifiers have to be in one
   block at the beginning, because they are used as mask bits.  There
   are 28 type modifiers; if we add many more we will have to redesign
   the mask mechanism.  */

enum rid
{
  /* Modifiers: */
  /* C, in empirical order of frequency.  */
  RID_STATIC = 0,
  RID_UNSIGNED, RID_LONG,    RID_CONST, RID_EXTERN,
  RID_REGISTER, RID_TYPEDEF, RID_SHORT, RID_INLINE,
  RID_VOLATILE, RID_SIGNED,  RID_AUTO,  RID_RESTRICT,
  RID_NORETURN, RID_ATOMIC,

  /* C extensions */
  RID_COMPLEX, RID_THREAD, RID_SAT,

  /* C++ */
  RID_FRIEND, RID_VIRTUAL, RID_EXPLICIT, RID_EXPORT, RID_MUTABLE,

  /* ObjC ("PQ" reserved words - they do not appear after a '@' and
     are keywords only in specific contexts)  */
  RID_IN, RID_OUT, RID_INOUT, RID_BYCOPY, RID_BYREF, RID_ONEWAY,

  /* ObjC ("PATTR" reserved words - they do not appear after a '@' 
     and are keywords only as property attributes)  */
  RID_GETTER, RID_SETTER,
  RID_READONLY, RID_READWRITE,
  RID_ASSIGN, RID_RETAIN, RID_COPY,
  RID_PROPATOMIC, RID_NONATOMIC,

  /* ObjC nullability support keywords that also can appear in the
     property attribute context.  These values should remain contiguous
     with the other property attributes.  */
  RID_NULL_UNSPECIFIED, RID_NULLABLE, RID_NONNULL, RID_NULL_RESETTABLE,

  /* C (reserved and imaginary types not implemented, so any use is a
     syntax error) */
  RID_IMAGINARY,

  /* C */
  RID_INT,     RID_CHAR,   RID_FLOAT,    RID_DOUBLE, RID_VOID,
  RID_ENUM,    RID_STRUCT, RID_UNION,    RID_IF,     RID_ELSE,
  RID_WHILE,   RID_DO,     RID_FOR,      RID_SWITCH, RID_CASE,
  RID_DEFAULT, RID_BREAK,  RID_CONTINUE, RID_RETURN, RID_GOTO,
  RID_SIZEOF,

  /* C extensions */
  RID_ASM,       RID_TYPEOF,   RID_ALIGNOF,  RID_ATTRIBUTE,  RID_VA_ARG,
  RID_EXTENSION, RID_IMAGPART, RID_REALPART, RID_LABEL,      RID_CHOOSE_EXPR,
  RID_TYPES_COMPATIBLE_P,      RID_BUILTIN_COMPLEX,	     RID_BUILTIN_SHUFFLE,
  RID_BUILTIN_SHUFFLEVECTOR,   RID_BUILTIN_CONVERTVECTOR,   RID_BUILTIN_TGMATH,
  RID_BUILTIN_HAS_ATTRIBUTE,   RID_BUILTIN_ASSOC_BARRIER,
  RID_DFLOAT32, RID_DFLOAT64, RID_DFLOAT128,

  /* TS 18661-3 keywords, in the same sequence as the TI_* values.  */
  RID_FLOAT16,
  RID_FLOATN_NX_FIRST = RID_FLOAT16,
  RID_FLOAT32,
  RID_FLOAT64,
  RID_FLOAT128,
  RID_FLOAT32X,
  RID_FLOAT64X,
  RID_FLOAT128X,
#define CASE_RID_FLOATN_NX						\
  case RID_FLOAT16: case RID_FLOAT32: case RID_FLOAT64: case RID_FLOAT128: \
  case RID_FLOAT32X: case RID_FLOAT64X: case RID_FLOAT128X

  RID_FRACT, RID_ACCUM, RID_AUTO_TYPE, RID_BUILTIN_CALL_WITH_STATIC_CHAIN,

  /* "__GIMPLE", for the GIMPLE-parsing extension to the C frontend. */
  RID_GIMPLE,

  /* "__PHI", for parsing PHI function in GIMPLE FE.  */
  RID_PHI,

  /* "__RTL", for the RTL-parsing extension to the C frontend.  */
  RID_RTL,

  /* C11 */
  RID_ALIGNAS, RID_GENERIC,

  /* This means to warn that this is a C++ keyword, and then treat it
     as a normal identifier.  */
  RID_CXX_COMPAT_WARN,

  /* GNU transactional memory extension */
  RID_TRANSACTION_ATOMIC, RID_TRANSACTION_RELAXED, RID_TRANSACTION_CANCEL,

  /* Too many ways of getting the name of a function as a string */
  RID_FUNCTION_NAME, RID_PRETTY_FUNCTION_NAME, RID_C99_FUNCTION_NAME,

  /* C++ (some of these are keywords in Objective-C as well, but only
     if they appear after a '@') */
  RID_BOOL,     RID_WCHAR,    RID_CLASS,
  RID_PUBLIC,   RID_PRIVATE,  RID_PROTECTED,
  RID_TEMPLATE, RID_NULL,     RID_CATCH,
  RID_DELETE,   RID_FALSE,    RID_NAMESPACE,
  RID_NEW,      RID_OFFSETOF, RID_OPERATOR,
  RID_THIS,     RID_THROW,    RID_TRUE,
  RID_TRY,      RID_TYPENAME, RID_TYPEID,
  RID_USING,    RID_CHAR16,   RID_CHAR32,

  /* casts */
  RID_CONSTCAST, RID_DYNCAST, RID_REINTCAST, RID_STATCAST,

  /* C++ extensions */
  RID_ADDRESSOF,               RID_BASES,
  RID_BUILTIN_LAUNDER,         RID_DIRECT_BASES,
  RID_HAS_NOTHROW_ASSIGN,      RID_HAS_NOTHROW_CONSTRUCTOR,
  RID_HAS_NOTHROW_COPY,        RID_HAS_TRIVIAL_ASSIGN,
  RID_HAS_TRIVIAL_CONSTRUCTOR, RID_HAS_TRIVIAL_COPY,
  RID_HAS_TRIVIAL_DESTRUCTOR,  RID_HAS_UNIQUE_OBJ_REPRESENTATIONS,
  RID_HAS_VIRTUAL_DESTRUCTOR,  RID_BUILTIN_BIT_CAST,
  RID_IS_ABSTRACT,             RID_IS_AGGREGATE,
  RID_IS_BASE_OF,              RID_IS_CLASS,
  RID_IS_EMPTY,                RID_IS_ENUM,
  RID_IS_FINAL,                RID_IS_LAYOUT_COMPATIBLE,
  RID_IS_LITERAL_TYPE,
  RID_IS_POINTER_INTERCONVERTIBLE_BASE_OF,
  RID_IS_POD,                  RID_IS_POLYMORPHIC,
  RID_IS_SAME_AS,
  RID_IS_STD_LAYOUT,           RID_IS_TRIVIAL,
  RID_IS_TRIVIALLY_ASSIGNABLE, RID_IS_TRIVIALLY_CONSTRUCTIBLE,
  RID_IS_TRIVIALLY_COPYABLE,
  RID_IS_UNION,                RID_UNDERLYING_TYPE,
  RID_IS_ASSIGNABLE,           RID_IS_CONSTRUCTIBLE,
  RID_IS_NOTHROW_ASSIGNABLE,   RID_IS_NOTHROW_CONSTRUCTIBLE,

  /* C++11 */
  RID_CONSTEXPR, RID_DECLTYPE, RID_NOEXCEPT, RID_NULLPTR, RID_STATIC_ASSERT,

  /* C++20 */
  RID_CONSTINIT, RID_CONSTEVAL,

  /* char8_t */
  RID_CHAR8,

  /* C++ concepts */
  RID_CONCEPT, RID_REQUIRES,

  /* C++ modules.  */
  RID__MODULE, RID__IMPORT, RID__EXPORT, /* Internal tokens.  */

  /* C++ coroutines */
  RID_CO_AWAIT, RID_CO_YIELD, RID_CO_RETURN,

  /* C++ transactional memory.  */
  RID_ATOMIC_NOEXCEPT, RID_ATOMIC_CANCEL, RID_SYNCHRONIZED,

  /* Objective-C ("AT" reserved words - they are only keywords when
     they follow '@')  */
  RID_AT_ENCODE,   RID_AT_END,
  RID_AT_CLASS,    RID_AT_ALIAS,     RID_AT_DEFS,
  RID_AT_PRIVATE,  RID_AT_PROTECTED, RID_AT_PUBLIC,  RID_AT_PACKAGE,
  RID_AT_PROTOCOL, RID_AT_SELECTOR,
  RID_AT_THROW,	   RID_AT_TRY,       RID_AT_CATCH,
  RID_AT_FINALLY,  RID_AT_SYNCHRONIZED, 
  RID_AT_OPTIONAL, RID_AT_REQUIRED, RID_AT_PROPERTY,
  RID_AT_SYNTHESIZE, RID_AT_DYNAMIC,
  RID_AT_INTERFACE,
  RID_AT_IMPLEMENTATION,

  /* Named address support, mapping the keyword to a particular named address
     number.  Named address space 0 is reserved for the generic address.  If
     there are more than 254 named addresses, the addr_space_t type will need
     to be grown from an unsigned char to unsigned short.  */
  RID_ADDR_SPACE_0,		/* generic address */
  RID_ADDR_SPACE_1,
  RID_ADDR_SPACE_2,
  RID_ADDR_SPACE_3,
  RID_ADDR_SPACE_4,
  RID_ADDR_SPACE_5,
  RID_ADDR_SPACE_6,
  RID_ADDR_SPACE_7,
  RID_ADDR_SPACE_8,
  RID_ADDR_SPACE_9,
  RID_ADDR_SPACE_10,
  RID_ADDR_SPACE_11,
  RID_ADDR_SPACE_12,
  RID_ADDR_SPACE_13,
  RID_ADDR_SPACE_14,
  RID_ADDR_SPACE_15,

  RID_FIRST_ADDR_SPACE = RID_ADDR_SPACE_0,
  RID_LAST_ADDR_SPACE = RID_ADDR_SPACE_15,

  /* __intN keywords.  The _N_M here doesn't correspond to the intN
     in the keyword; use the bitsize in int_n_t_data_t[M] for that.
     For example, if int_n_t_data_t[0].bitsize is 13, then RID_INT_N_0
     is for __int13.  */

  /* Note that the range to use is RID_FIRST_INT_N through
     RID_FIRST_INT_N + NUM_INT_N_ENTS - 1 and c-parser.cc has a list of
     all RID_INT_N_* in a case statement.  */

  RID_INT_N_0,
  RID_INT_N_1,
  RID_INT_N_2,
  RID_INT_N_3,

  RID_FIRST_INT_N = RID_INT_N_0,
  RID_LAST_INT_N = RID_INT_N_3,

  RID_MAX,

  RID_FIRST_MODIFIER = RID_STATIC,
  RID_LAST_MODIFIER = RID_ONEWAY,

  RID_FIRST_CXX11 = RID_CONSTEXPR,
  RID_LAST_CXX11 = RID_STATIC_ASSERT,
  RID_FIRST_CXX20 = RID_CONSTINIT,
  RID_LAST_CXX20 = RID_CONSTINIT,
  RID_FIRST_AT = RID_AT_ENCODE,
  RID_LAST_AT = RID_AT_IMPLEMENTATION,
  RID_FIRST_PQ = RID_IN,
  RID_LAST_PQ = RID_ONEWAY,
  RID_FIRST_PATTR = RID_GETTER,
  RID_LAST_PATTR = RID_NULL_RESETTABLE
};

#define OBJC_IS_AT_KEYWORD(rid) \
  ((unsigned int) (rid) >= (unsigned int) RID_FIRST_AT && \
   (unsigned int) (rid) <= (unsigned int) RID_LAST_AT)

#define OBJC_IS_PQ_KEYWORD(rid) \
  ((unsigned int) (rid) >= (unsigned int) RID_FIRST_PQ && \
   (unsigned int) (rid) <= (unsigned int) RID_LAST_PQ)

/* Keywords permitted in an @property attribute context.  */
#define OBJC_IS_PATTR_KEYWORD(rid) \
  ((((unsigned int) (rid) >= (unsigned int) RID_FIRST_PATTR && \
     (unsigned int) (rid) <= (unsigned int) RID_LAST_PATTR)) \
   || rid == RID_CLASS)

/* OBJC_IS_CXX_KEYWORD recognizes the 'CXX_OBJC' keywords (such as
   'class') which are shared in a subtle way between Objective-C and
   C++.  When the lexer is lexing in Objective-C/Objective-C++, if it
   finds '@' followed by one of these identifiers (eg, '@class'), it
   recognizes the whole as an Objective-C keyword.  If the identifier
   is found elsewhere, it follows the rules of the C/C++ language.
 */
#define OBJC_IS_CXX_KEYWORD(rid) \
  (rid == RID_CLASS || rid == RID_SYNCHRONIZED			\
   || rid == RID_PUBLIC || rid == RID_PROTECTED || rid == RID_PRIVATE	\
   || rid == RID_TRY || rid == RID_THROW || rid == RID_CATCH)

/* The elements of `ridpointers' are identifier nodes for the reserved
   type names and storage classes.  It is indexed by a RID_... value.  */
extern GTY ((length ("(int) RID_MAX"))) tree *ridpointers;

/* Standard named or nameless data types of the C compiler.  */

enum c_tree_index
{
    CTI_CHAR8_TYPE,
    CTI_CHAR16_TYPE,
    CTI_CHAR32_TYPE,
    CTI_WCHAR_TYPE,
    CTI_UNDERLYING_WCHAR_TYPE,
    CTI_WINT_TYPE,
    CTI_SIGNED_SIZE_TYPE, /* For format checking only.  */
    CTI_UNSIGNED_PTRDIFF_TYPE, /* For format checking only.  */
    CTI_INTMAX_TYPE,
    CTI_UINTMAX_TYPE,
    CTI_WIDEST_INT_LIT_TYPE,
    CTI_WIDEST_UINT_LIT_TYPE,

    /* Types for <stdint.h>, that may not be defined on all
       targets.  */
    CTI_SIG_ATOMIC_TYPE,
    CTI_INT8_TYPE,
    CTI_INT16_TYPE,
    CTI_INT32_TYPE,
    CTI_INT64_TYPE,
    CTI_UINT8_TYPE,
    CTI_UINT16_TYPE,
    CTI_UINT32_TYPE,
    CTI_UINT64_TYPE,
    CTI_INT_LEAST8_TYPE,
    CTI_INT_LEAST16_TYPE,
    CTI_INT_LEAST32_TYPE,
    CTI_INT_LEAST64_TYPE,
    CTI_UINT_LEAST8_TYPE,
    CTI_UINT_LEAST16_TYPE,
    CTI_UINT_LEAST32_TYPE,
    CTI_UINT_LEAST64_TYPE,
    CTI_INT_FAST8_TYPE,
    CTI_INT_FAST16_TYPE,
    CTI_INT_FAST32_TYPE,
    CTI_INT_FAST64_TYPE,
    CTI_UINT_FAST8_TYPE,
    CTI_UINT_FAST16_TYPE,
    CTI_UINT_FAST32_TYPE,
    CTI_UINT_FAST64_TYPE,
    CTI_INTPTR_TYPE,
    CTI_UINTPTR_TYPE,

    CTI_CHAR_ARRAY_TYPE,
    CTI_CHAR8_ARRAY_TYPE,
    CTI_CHAR16_ARRAY_TYPE,
    CTI_CHAR32_ARRAY_TYPE,
    CTI_WCHAR_ARRAY_TYPE,
    CTI_STRING_TYPE,
    CTI_CONST_STRING_TYPE,

    /* Type for boolean expressions (bool in C++, int in C).  */
    CTI_TRUTHVALUE_TYPE,
    CTI_TRUTHVALUE_TRUE,
    CTI_TRUTHVALUE_FALSE,

    CTI_DEFAULT_FUNCTION_TYPE,

    CTI_NULL,

    /* These are not types, but we have to look them up all the time.  */
    CTI_FUNCTION_NAME_DECL,
    CTI_PRETTY_FUNCTION_NAME_DECL,
    CTI_C99_FUNCTION_NAME_DECL,

    CTI_MODULE_HWM,
    /* Below here entities change during compilation.  */

    CTI_SAVED_FUNCTION_NAME_DECLS,

    CTI_MAX
};

#define C_CPP_HASHNODE(id) \
  (&(((struct c_common_identifier *) (id))->node))
#define C_RID_CODE(id) \
  ((enum rid) (((struct c_common_identifier *) (id))->node.rid_code))
#define C_SET_RID_CODE(id, code) \
  (((struct c_common_identifier *) (id))->node.rid_code = (unsigned char) code)

/* Identifier part common to the C front ends.  Inherits from
   tree_identifier, despite appearances.  */
struct GTY(()) c_common_identifier {
  struct tree_common common;
  struct cpp_hashnode node;
};

/* An entry in the reserved keyword table.  */

struct c_common_resword
{
  const char *const word;
  ENUM_BITFIELD(rid) const rid : 16;
  const unsigned int disable   : 16;
};

/* Mode used to build pointers (VOIDmode means ptr_mode).  */

extern machine_mode c_default_pointer_mode;

/* Extra cpp_ttype values for C++.  */

/* A token type for template-ids.  If a template-id is processed while
   parsing tentatively, it is replaced with a CPP_TEMPLATE_ID token;
   the value of the CPP_TEMPLATE_ID is whatever was returned by
   cp_parser_template_id.  */
#define CPP_TEMPLATE_ID ((enum cpp_ttype) (CPP_KEYWORD + 1))

/* A token type for nested-name-specifiers.  If a
   nested-name-specifier is processed while parsing tentatively, it is
   replaced with a CPP_NESTED_NAME_SPECIFIER token; the value of the
   CPP_NESTED_NAME_SPECIFIER is whatever was returned by
   cp_parser_nested_name_specifier_opt.  */
#define CPP_NESTED_NAME_SPECIFIER ((enum cpp_ttype) (CPP_TEMPLATE_ID + 1))

/* A token type for pre-parsed C++0x decltype.  */
#define CPP_DECLTYPE ((enum cpp_ttype) (CPP_NESTED_NAME_SPECIFIER + 1))

/* A token type for pre-parsed primary-expression (lambda- or statement-).  */
#define CPP_PREPARSED_EXPR ((enum cpp_ttype) (CPP_DECLTYPE + 1))

/* The number of token types, including C++-specific ones.  */
#define N_CP_TTYPES ((int) (CPP_PREPARSED_EXPR + 1))

/* Disable mask.  Keywords are disabled if (reswords[i].disable &
   mask) is _true_.  Thus for keywords which are present in all
   languages the disable field is zero.  */

#define D_CONLY		0x0001	/* C only (not in C++).  */
#define D_CXXONLY	0x0002	/* C++ only (not in C).  */
#define D_C99		0x0004	/* In C, C99 only.  */
#define D_CXX11         0x0008	/* In C++, C++11 only.  */
#define D_EXT		0x0010	/* GCC extension.  */
#define D_EXT89		0x0020	/* GCC extension incorporated in C99.  */
#define D_ASM		0x0040	/* Disabled by -fno-asm.  */
#define D_OBJC		0x0080	/* In Objective C and neither C nor C++.  */
#define D_CXX_OBJC	0x0100	/* In Objective C, and C++, but not C.  */
#define D_CXXWARN	0x0200	/* In C warn with -Wcxx-compat.  */
#define D_CXX_CONCEPTS  0x0400	/* In C++, only with concepts.  */
#define D_TRANSMEM	0X0800	/* C++ transactional memory TS.  */
#define D_CXX_CHAR8_T	0X1000	/* In C++, only with -fchar8_t.  */
#define D_CXX20		0x2000  /* In C++, C++20 only.  */
#define D_CXX_COROUTINES 0x4000  /* In C++, only with coroutines.  */
#define D_CXX_MODULES	0x8000  /* In C++, only with modules.  */

#define D_CXX_CONCEPTS_FLAGS D_CXXONLY | D_CXX_CONCEPTS
#define D_CXX_CHAR8_T_FLAGS D_CXXONLY | D_CXX_CHAR8_T
#define D_CXX_MODULES_FLAGS (D_CXXONLY | D_CXX_MODULES)
#define D_CXX_COROUTINES_FLAGS (D_CXXONLY | D_CXX_COROUTINES)

/* The reserved keyword table.  */
extern const struct c_common_resword c_common_reswords[];

/* The number of items in the reserved keyword table.  */
extern const unsigned int num_c_common_reswords;

#define char8_type_node			c_global_trees[CTI_CHAR8_TYPE]
#define char16_type_node		c_global_trees[CTI_CHAR16_TYPE]
#define char32_type_node		c_global_trees[CTI_CHAR32_TYPE]
#define wchar_type_node			c_global_trees[CTI_WCHAR_TYPE]
#define underlying_wchar_type_node	c_global_trees[CTI_UNDERLYING_WCHAR_TYPE]
#define wint_type_node			c_global_trees[CTI_WINT_TYPE]
#define signed_size_type_node		c_global_trees[CTI_SIGNED_SIZE_TYPE]
#define unsigned_ptrdiff_type_node	c_global_trees[CTI_UNSIGNED_PTRDIFF_TYPE]
#define intmax_type_node		c_global_trees[CTI_INTMAX_TYPE]
#define uintmax_type_node		c_global_trees[CTI_UINTMAX_TYPE]
#define widest_integer_literal_type_node c_global_trees[CTI_WIDEST_INT_LIT_TYPE]
#define widest_unsigned_literal_type_node c_global_trees[CTI_WIDEST_UINT_LIT_TYPE]

#define sig_atomic_type_node		c_global_trees[CTI_SIG_ATOMIC_TYPE]
#define int8_type_node			c_global_trees[CTI_INT8_TYPE]
#define int16_type_node			c_global_trees[CTI_INT16_TYPE]
#define int32_type_node			c_global_trees[CTI_INT32_TYPE]
#define int64_type_node			c_global_trees[CTI_INT64_TYPE]
#define uint8_type_node			c_global_trees[CTI_UINT8_TYPE]
#define c_uint16_type_node		c_global_trees[CTI_UINT16_TYPE]
#define c_uint32_type_node		c_global_trees[CTI_UINT32_TYPE]
#define c_uint64_type_node		c_global_trees[CTI_UINT64_TYPE]
#define int_least8_type_node		c_global_trees[CTI_INT_LEAST8_TYPE]
#define int_least16_type_node		c_global_trees[CTI_INT_LEAST16_TYPE]
#define int_least32_type_node		c_global_trees[CTI_INT_LEAST32_TYPE]
#define int_least64_type_node		c_global_trees[CTI_INT_LEAST64_TYPE]
#define uint_least8_type_node		c_global_trees[CTI_UINT_LEAST8_TYPE]
#define uint_least16_type_node		c_global_trees[CTI_UINT_LEAST16_TYPE]
#define uint_least32_type_node		c_global_trees[CTI_UINT_LEAST32_TYPE]
#define uint_least64_type_node		c_global_trees[CTI_UINT_LEAST64_TYPE]
#define int_fast8_type_node		c_global_trees[CTI_INT_FAST8_TYPE]
#define int_fast16_type_node		c_global_trees[CTI_INT_FAST16_TYPE]
#define int_fast32_type_node		c_global_trees[CTI_INT_FAST32_TYPE]
#define int_fast64_type_node		c_global_trees[CTI_INT_FAST64_TYPE]
#define uint_fast8_type_node		c_global_trees[CTI_UINT_FAST8_TYPE]
#define uint_fast16_type_node		c_global_trees[CTI_UINT_FAST16_TYPE]
#define uint_fast32_type_node		c_global_trees[CTI_UINT_FAST32_TYPE]
#define uint_fast64_type_node		c_global_trees[CTI_UINT_FAST64_TYPE]
#define intptr_type_node		c_global_trees[CTI_INTPTR_TYPE]
#define uintptr_type_node		c_global_trees[CTI_UINTPTR_TYPE]

#define truthvalue_type_node		c_global_trees[CTI_TRUTHVALUE_TYPE]
#define truthvalue_true_node		c_global_trees[CTI_TRUTHVALUE_TRUE]
#define truthvalue_false_node		c_global_trees[CTI_TRUTHVALUE_FALSE]

#define char_array_type_node		c_global_trees[CTI_CHAR_ARRAY_TYPE]
#define char8_array_type_node		c_global_trees[CTI_CHAR8_ARRAY_TYPE]
#define char16_array_type_node		c_global_trees[CTI_CHAR16_ARRAY_TYPE]
#define char32_array_type_node		c_global_trees[CTI_CHAR32_ARRAY_TYPE]
#define wchar_array_type_node		c_global_trees[CTI_WCHAR_ARRAY_TYPE]
#define string_type_node		c_global_trees[CTI_STRING_TYPE]
#define const_string_type_node		c_global_trees[CTI_CONST_STRING_TYPE]

#define default_function_type		c_global_trees[CTI_DEFAULT_FUNCTION_TYPE]

#define function_name_decl_node		c_global_trees[CTI_FUNCTION_NAME_DECL]
#define pretty_function_name_decl_node	c_global_trees[CTI_PRETTY_FUNCTION_NAME_DECL]
#define c99_function_name_decl_node		c_global_trees[CTI_C99_FUNCTION_NAME_DECL]
#define saved_function_name_decls	c_global_trees[CTI_SAVED_FUNCTION_NAME_DECLS]

/* The node for C++ `__null'.  */
#define null_node                       c_global_trees[CTI_NULL]

extern GTY(()) tree c_global_trees[CTI_MAX];

/* Mark which labels are explicitly declared.
   These may be shadowed, and may be referenced from nested functions.  */
#define C_DECLARED_LABEL_FLAG(label) TREE_LANG_FLAG_1 (label)

enum c_language_kind
{
  clk_c		= 0,		/* C90, C94, C99, C11 or C2X */
  clk_objc	= 1,		/* clk_c with ObjC features.  */
  clk_cxx	= 2,		/* ANSI/ISO C++ */
  clk_objcxx	= 3		/* clk_cxx with ObjC features.  */
};

/* To test for a specific language use c_language, defined by each
   front end.  For "ObjC features" or "not C++" use the macros.  */
extern c_language_kind c_language;

#define c_dialect_cxx()		((c_language & clk_cxx) != 0)
#define c_dialect_objc()	((c_language & clk_objc) != 0)

/* The various name of operator that appears in error messages. */
enum ref_operator {
  /* NULL */
  RO_NULL,
  /* array indexing */
  RO_ARRAY_INDEXING,
  /* unary * */
  RO_UNARY_STAR,
  /* -> */
  RO_ARROW,
  /* implicit conversion */
  RO_IMPLICIT_CONVERSION,
  /* ->* */
  RO_ARROW_STAR
};

/* Information about a statement tree.  */

struct GTY(()) stmt_tree_s {
  /* A stack of statement lists being collected.  */
  vec<tree, va_gc> *x_cur_stmt_list;

  /* In C++, Nonzero if we should treat statements as full
     expressions.  In particular, this variable is non-zero if at the
     end of a statement we should destroy any temporaries created
     during that statement.  Similarly, if, at the end of a block, we
     should destroy any local variables in this block.  Normally, this
     variable is nonzero, since those are the normal semantics of
     C++.

     This flag has no effect in C.  */
  int stmts_are_full_exprs_p;
};

typedef struct stmt_tree_s *stmt_tree;

/* Global state pertinent to the current function.  Some C dialects
   extend this structure with additional fields.  */

struct GTY(()) c_language_function {
  /* While we are parsing the function, this contains information
     about the statement-tree that we are building.  */
  struct stmt_tree_s x_stmt_tree;

  /* Vector of locally defined typedefs, for
     -Wunused-local-typedefs.  */
  vec<tree, va_gc> *local_typedefs;
};

#define stmt_list_stack (current_stmt_tree ()->x_cur_stmt_list)

/* When building a statement-tree, this is the current statement list
   being collected.  */
#define cur_stmt_list	(stmt_list_stack->last ())

#define building_stmt_list_p() (stmt_list_stack && !stmt_list_stack->is_empty())

/* Language-specific hooks.  */

/* If non-NULL, this function is called after a precompile header file
   is loaded.  */
extern void (*lang_post_pch_load) (void);

extern void push_file_scope (void);
extern void pop_file_scope (void);
extern stmt_tree current_stmt_tree (void);
extern tree push_stmt_list (void);
extern tree pop_stmt_list (tree);
extern tree add_stmt (tree);
extern void push_cleanup (tree, tree, bool);

extern tree build_modify_expr (location_t, tree, tree, enum tree_code,
			       location_t, tree, tree);
extern tree build_indirect_ref (location_t, tree, ref_operator);

extern bool has_c_linkage (const_tree decl);
extern bool c_decl_implicit (const_tree);

/* Switches common to the C front ends.  */

/* Nonzero means don't output line number information.  */

extern char flag_no_line_commands;

/* Nonzero causes -E output not to be done, but directives such as
   #define that have side effects are still obeyed.  */

extern char flag_no_output;

/* Nonzero means dump macros in some fashion; contains the 'D', 'M',
   'N' or 'U' of the command line switch.  */

extern char flag_dump_macros;

/* Nonzero means pass #include lines through to the output.  */

extern char flag_dump_includes;

/* Nonzero means process PCH files while preprocessing.  */

extern bool flag_pch_preprocess;

/* The file name to which we should write a precompiled header, or
   NULL if no header will be written in this compile.  */

extern const char *pch_file;

/* Nonzero if an ISO standard was selected.  It rejects macros in the
   user's namespace.  */

extern int flag_iso;

/* C/ObjC language option variables.  */


/* Nonzero means allow type mismatches in conditional expressions;
   just make their values `void'.  */

extern int flag_cond_mismatch;

/* Nonzero means enable C89 Amendment 1 features.  */

extern int flag_isoc94;

/* Nonzero means use the ISO C99 (or later) dialect of C.  */

extern int flag_isoc99;

/* Nonzero means use the ISO C11 (or later) dialect of C.  */

extern int flag_isoc11;

/* Nonzero means use the ISO C2X dialect of C.  */

extern int flag_isoc2x;

/* Nonzero means that we have builtin functions, and main is an int.  */

extern int flag_hosted;

/* ObjC language option variables.  */


/* Tells the compiler that this is a special run.  Do not perform any
   compiling, instead we are to test some platform dependent features
   and output a C header file with appropriate definitions.  */

extern int print_struct_values;

/* Tells the compiler what is the constant string class for ObjC.  */

extern const char *constant_string_class_name;


/* C++ language option variables.  */

/* The reference version of the ABI for -Wabi.  */

extern int warn_abi_version;

/* Return TRUE if one of {flag_abi_version,flag_abi_compat_version} is
   less than N and the other is at least N.  */
#define abi_compat_version_crosses(N)		\
  (abi_version_at_least(N)			\
   != (flag_abi_compat_version == 0		\
       || flag_abi_compat_version >= (N)))

/* Return TRUE if one of {flag_abi_version,warn_abi_version} is
   less than N and the other is at least N, for use by -Wabi.  */
#define abi_version_crosses(N)			\
  (abi_version_at_least(N)			\
   != (warn_abi_version == 0			\
       || warn_abi_version >= (N)))

/* The supported C++ dialects.  */

enum cxx_dialect {
  cxx_unset,
  /* C++98 with TC1  */
  cxx98,
  cxx03 = cxx98,
  /* C++11  */
  cxx0x,
  cxx11 = cxx0x,
  /* C++14 */
  cxx14,
  /* C++17 */
  cxx17,
  /* C++20 */
  cxx20,
  /* C++23 */
  cxx23
};

/* The C++ dialect being used. C++98 is the default.  */
extern enum cxx_dialect cxx_dialect;

/* Maximum template instantiation depth.  This limit is rather
   arbitrary, but it exists to limit the time it takes to notice
   excessively recursive template instantiations.  */

extern int max_tinst_depth;

/* Nonzero means that we should not issue warnings about problems that
   occur when the code is executed, because the code being processed
   is not expected to be executed.  This is set during parsing.  This
   is used for cases like sizeof() and "0 ? a : b".  This is a count,
   not a bool, because unexecuted expressions can nest.  */

extern int c_inhibit_evaluation_warnings;

/* Whether lexing has been completed, so subsequent preprocessor
   errors should use the compiler's input_location.  */

extern bool done_lexing;

/* C types are partitioned into three subsets: object, function, and
   incomplete types.  */
#define C_TYPE_OBJECT_P(type) \
  (TREE_CODE (type) != FUNCTION_TYPE && TYPE_SIZE (type))

#define C_TYPE_INCOMPLETE_P(type) \
  (TREE_CODE (type) != FUNCTION_TYPE && TYPE_SIZE (type) == 0)

#define C_TYPE_FUNCTION_P(type) \
  (TREE_CODE (type) == FUNCTION_TYPE)

/* For convenience we define a single macro to identify the class of
   object or incomplete types.  */
#define C_TYPE_OBJECT_OR_INCOMPLETE_P(type) \
  (!C_TYPE_FUNCTION_P (type))

/* Return true if TYPE is a vector type that should be subject to the GNU
   vector extensions (as opposed to a vector type that is used only for
   the purposes of defining target-specific built-in functions).  */

inline bool
gnu_vector_type_p (const_tree type)
{
  return TREE_CODE (type) == VECTOR_TYPE && !TYPE_INDIVISIBLE_P (type);
}

struct visibility_flags
{
  unsigned inpragma : 1;	/* True when in #pragma GCC visibility.  */
  unsigned inlines_hidden : 1;	/* True when -finlineshidden in effect.  */
};

/* These enumerators are possible types of unsafe conversions.  */
enum conversion_safety {
  /* The conversion is safe.  */
  SAFE_CONVERSION = 0,
  /* Another type of conversion with problems.  */
  UNSAFE_OTHER,
  /* Conversion between signed and unsigned integers.  */
  UNSAFE_SIGN,
  /* Conversions that reduce the precision of reals including conversions
     from reals to integers.  */
  UNSAFE_REAL,
  /* Conversions from complex to reals or integers, that discard imaginary
     component.  */
  UNSAFE_IMAGINARY
};

/* Global visibility options.  */
extern struct visibility_flags visibility_options;

/* Attribute table common to the C front ends.  */
extern const struct attribute_spec c_common_attribute_table[];
extern const struct attribute_spec c_common_format_attribute_table[];

/* Pointer to function to lazily generate the VAR_DECL for __FUNCTION__ etc.
   ID is the identifier to use, NAME is the string.
   TYPE_DEP indicates whether it depends on type of the function or not
   (i.e. __PRETTY_FUNCTION__).  */

extern tree (*make_fname_decl) (location_t, tree, int);

/* In c-decl.cc and cp/tree.cc.  FIXME.  */
extern void c_register_addr_space (const char *str, addr_space_t as);

/* In c-common.cc.  */
extern bool in_late_binary_op;
extern const char *c_addr_space_name (addr_space_t as);
extern tree identifier_global_value (tree);
extern tree identifier_global_tag (tree);
extern bool names_builtin_p (const char *);
extern tree c_linkage_bindings (tree);
extern void record_builtin_type (enum rid, const char *, tree);
extern tree build_void_list_node (void);
extern void start_fname_decls (void);
extern void finish_fname_decls (void);
extern const char *fname_as_string (int);
extern tree fname_decl (location_t, unsigned, tree);

extern int check_user_alignment (const_tree, bool, bool);
extern bool check_function_arguments (location_t loc, const_tree, const_tree,
				      int, tree *, vec<location_t> *);
extern void check_function_arguments_recurse (void (*)
					      (void *, tree,
					       unsigned HOST_WIDE_INT),
					      void *, tree,
					      unsigned HOST_WIDE_INT,
					      opt_code);
extern bool check_builtin_function_arguments (location_t, vec<location_t>,
					      tree, tree, int, tree *);
extern void check_function_format (const_tree, tree, int, tree *,
				   vec<location_t> *);
extern bool attribute_fallthrough_p (tree);
extern tree handle_format_attribute (tree *, tree, tree, int, bool *);
extern tree handle_format_arg_attribute (tree *, tree, tree, int, bool *);
extern bool c_common_handle_option (size_t, const char *, HOST_WIDE_INT, int,
				    location_t,
				    const struct cl_option_handlers *);
extern bool default_handle_c_option (size_t, const char *, int);
extern tree c_common_type_for_mode (machine_mode, int);
extern tree c_common_type_for_size (unsigned int, int);
extern tree c_common_fixed_point_type_for_size (unsigned int, unsigned int,
						int, int);
extern tree c_common_unsigned_type (tree);
extern tree c_common_signed_type (tree);
extern tree c_common_signed_or_unsigned_type (int, tree);
extern void c_common_init_ts (void);
extern tree c_build_bitfield_integer_type (unsigned HOST_WIDE_INT, int);
extern enum conversion_safety unsafe_conversion_p (tree, tree, tree, bool);
extern bool decl_with_nonnull_addr_p (const_tree);
extern tree c_fully_fold (tree, bool, bool *, bool = false);
extern tree c_wrap_maybe_const (tree, bool);
extern tree c_common_truthvalue_conversion (location_t, tree);
extern void c_apply_type_quals_to_decl (int, tree);
extern tree c_sizeof_or_alignof_type (location_t, tree, bool, bool, int);
extern tree c_alignof_expr (location_t, tree);
/* Print an error message for invalid operands to arith operation CODE.
   NOP_EXPR is used as a special case (see truthvalue_conversion).  */
extern void binary_op_error (rich_location *, enum tree_code, tree, tree);
extern tree fix_string_type (tree);
extern tree convert_and_check (location_t, tree, tree, bool = false);
extern bool c_determine_visibility (tree);
extern bool vector_types_compatible_elements_p (tree, tree);
extern void mark_valid_location_for_stdc_pragma (bool);
extern bool valid_location_for_stdc_pragma_p (void);
extern void set_float_const_decimal64 (void);
extern void clear_float_const_decimal64 (void);
extern bool float_const_decimal64_p (void);

extern bool keyword_begins_type_specifier (enum rid);
extern bool keyword_is_storage_class_specifier (enum rid);
extern bool keyword_is_type_qualifier (enum rid);
extern bool keyword_is_decl_specifier (enum rid);
extern unsigned max_align_t_align (void);
extern bool cxx_fundamental_alignment_p (unsigned);
extern bool pointer_to_zero_sized_aggr_p (tree);
extern bool bool_promoted_to_int_p (tree);
extern tree fold_for_warn (tree);
extern tree c_common_get_narrower (tree, int *);
extern bool get_attribute_operand (tree, unsigned HOST_WIDE_INT *);
extern void c_common_finalize_early_debug (void);

/* Used by convert_and_check; in front ends.  */
extern tree convert_init (tree, tree);

#define c_sizeof(LOC, T)  c_sizeof_or_alignof_type (LOC, T, true, false, 1)
#define c_alignof(LOC, T) c_sizeof_or_alignof_type (LOC, T, false, false, 1)

/* Subroutine of build_binary_op, used for certain operations.  */
extern tree shorten_binary_op (tree result_type, tree op0, tree op1, bool bitwise);

/* Subroutine of build_binary_op, used for comparison operations.
   See if the operands have both been converted from subword integer types
   and, if so, perhaps change them both back to their original type.  */
extern tree shorten_compare (location_t, tree *, tree *, tree *,
			     enum tree_code *);

extern tree pointer_int_sum (location_t, enum tree_code, tree, tree,
			     bool = true);

/* Add qualifiers to a type, in the fashion for C.  */
extern tree c_build_qualified_type (tree, int, tree = NULL_TREE, size_t = 0);

/* Build tree nodes and builtin functions common to both C and C++ language
   frontends.  */
extern void c_common_nodes_and_builtins (void);

extern void disable_builtin_function (const char *);

extern void set_compound_literal_name (tree decl);

extern tree build_va_arg (location_t, tree, tree);

extern const unsigned int c_family_lang_mask;
extern unsigned int c_common_option_lang_mask (void);
extern void c_common_diagnostics_set_defaults (diagnostic_context *);
extern bool c_common_complain_wrong_lang_p (const struct cl_option *);
extern void c_common_init_options_struct (struct gcc_options *);
extern void c_common_init_options (unsigned int, struct cl_decoded_option *);
extern bool c_common_post_options (const char **);
extern bool c_common_init (void);
extern void c_common_finish (void);
extern void c_common_parse_file (void);
extern FILE *get_dump_info (int, dump_flags_t *);
extern alias_set_type c_common_get_alias_set (tree);
extern void c_register_builtin_type (tree, const char*);
extern bool c_promoting_integer_type_p (const_tree);
extern bool self_promoting_args_p (const_tree);
extern tree strip_pointer_operator (tree);
extern tree strip_pointer_or_array_types (tree);
extern HOST_WIDE_INT c_common_to_target_charset (HOST_WIDE_INT);

/* This is the basic parsing function.  */
extern void c_parse_file (void);

extern void c_parse_final_cleanups (void);

/* These macros provide convenient access to the various _STMT nodes.  */

/* Nonzero if a given STATEMENT_LIST represents the outermost binding
   if a statement expression.  */
#define STATEMENT_LIST_STMT_EXPR(NODE) \
  TREE_LANG_FLAG_1 (STATEMENT_LIST_CHECK (NODE))

/* Nonzero if a label has been added to the statement list.  */
#define STATEMENT_LIST_HAS_LABEL(NODE) \
  TREE_LANG_FLAG_3 (STATEMENT_LIST_CHECK (NODE))

/* C_MAYBE_CONST_EXPR accessors.  */
#define C_MAYBE_CONST_EXPR_PRE(NODE)			\
  TREE_OPERAND (C_MAYBE_CONST_EXPR_CHECK (NODE), 0)
#define C_MAYBE_CONST_EXPR_EXPR(NODE)			\
  TREE_OPERAND (C_MAYBE_CONST_EXPR_CHECK (NODE), 1)
#define C_MAYBE_CONST_EXPR_INT_OPERANDS(NODE)		\
  TREE_LANG_FLAG_0 (C_MAYBE_CONST_EXPR_CHECK (NODE))
#define C_MAYBE_CONST_EXPR_NON_CONST(NODE)		\
  TREE_LANG_FLAG_1 (C_MAYBE_CONST_EXPR_CHECK (NODE))
#define EXPR_INT_CONST_OPERANDS(EXPR)			\
  (INTEGRAL_TYPE_P (TREE_TYPE (EXPR))			\
   && (TREE_CODE (EXPR) == INTEGER_CST			\
       || (TREE_CODE (EXPR) == C_MAYBE_CONST_EXPR	\
	   && C_MAYBE_CONST_EXPR_INT_OPERANDS (EXPR))))

/* In a FIELD_DECL, nonzero if the decl was originally a bitfield.  */
#define DECL_C_BIT_FIELD(NODE) \
  (DECL_LANG_FLAG_4 (FIELD_DECL_CHECK (NODE)) == 1)
#define SET_DECL_C_BIT_FIELD(NODE) \
  (DECL_LANG_FLAG_4 (FIELD_DECL_CHECK (NODE)) = 1)
#define CLEAR_DECL_C_BIT_FIELD(NODE) \
  (DECL_LANG_FLAG_4 (FIELD_DECL_CHECK (NODE)) = 0)

/* True if the decl was an unnamed bitfield.  */
#define DECL_UNNAMED_BIT_FIELD(NODE) \
  (DECL_C_BIT_FIELD (NODE) && !DECL_NAME (NODE))

extern tree do_case (location_t, tree, tree);
extern tree build_stmt (location_t, enum tree_code, ...);
extern tree build_real_imag_expr (location_t, enum tree_code, tree);

/* These functions must be defined by each front-end which implements
   a variant of the C language.  They are used in c-common.cc.  */

extern tree build_unary_op (location_t, enum tree_code, tree, bool);
extern tree build_binary_op (location_t, enum tree_code, tree, tree, bool);
extern tree perform_integral_promotions (tree);

/* These functions must be defined by each front-end which implements
   a variant of the C language.  They are used by port files.  */

extern tree default_conversion (tree);

/* Given two integer or real types, return the type for their sum.
   Given two compatible ANSI C types, returns the merged type.  */

extern tree common_type (tree, tree);

extern tree decl_constant_value (tree);

/* Handle increment and decrement of boolean types.  */
extern tree boolean_increment (enum tree_code, tree);

extern int case_compare (splay_tree_key, splay_tree_key);

extern tree c_add_case_label (location_t, splay_tree, tree, tree, tree);
extern bool c_switch_covers_all_cases_p (splay_tree, tree);
extern bool c_block_may_fallthru (const_tree);

extern tree build_function_call (location_t, tree, tree);

extern tree build_function_call_vec (location_t, vec<location_t>, tree,
				     vec<tree, va_gc> *, vec<tree, va_gc> *,
				     tree = NULL_TREE);

extern tree resolve_overloaded_builtin (location_t, tree, vec<tree, va_gc> *);

extern tree finish_label_address_expr (tree, location_t);

/* Same function prototype, but the C and C++ front ends have
   different implementations.  Used in c-common.cc.  */
extern tree lookup_label (tree);
extern tree lookup_name (tree);
extern bool lvalue_p (const_tree);
extern bool instantiation_dependent_expression_p (tree);

extern bool vector_targets_convertible_p (const_tree t1, const_tree t2);
extern bool vector_types_convertible_p (const_tree t1, const_tree t2, bool emit_lax_note);
extern tree c_build_vec_perm_expr (location_t, tree, tree, tree, bool = true);
extern tree c_build_shufflevector (location_t, tree, tree,
				   const vec<tree> &, bool = true);
extern tree c_build_vec_convert (location_t, tree, location_t, tree, bool = true);

extern void init_c_lex (void);

extern void c_cpp_builtins (cpp_reader *);
extern void c_cpp_builtins_optimize_pragma (cpp_reader *, tree, tree);
extern bool c_cpp_diagnostic (cpp_reader *, enum cpp_diagnostic_level,
			      enum cpp_warning_reason, rich_location *,
			      const char *, va_list *)
     ATTRIBUTE_GCC_DIAG(5,0);
extern int c_common_has_attribute (cpp_reader *, bool);
extern int c_common_has_builtin (cpp_reader *);

extern bool parse_optimize_options (tree, bool);

/* Positive if an implicit `extern "C"' scope has just been entered;
   negative if such a scope has just been exited.  */
extern GTY(()) int pending_lang_change;

/* Information recorded about each file examined during compilation.  */

struct c_fileinfo
{
  int time;	/* Time spent in the file.  */

  /* Flags used only by C++.
     INTERFACE_ONLY nonzero means that we are in an "interface" section
     of the compiler.  INTERFACE_UNKNOWN nonzero means we cannot trust
     the value of INTERFACE_ONLY.  If INTERFACE_UNKNOWN is zero and
     INTERFACE_ONLY is zero, it means that we are responsible for
     exporting definitions that others might need.  */
  short interface_only;
  short interface_unknown;
};

struct c_fileinfo *get_fileinfo (const char *);
extern void dump_time_statistics (void);

extern bool c_dump_tree (void *, tree);

extern void verify_sequence_points (tree);

extern tree fold_offsetof (tree, tree = size_type_node,
			   tree_code ctx = ERROR_MARK);

extern int complete_array_type (tree *, tree, bool);
extern void complete_flexible_array_elts (tree);

extern tree builtin_type_for_size (int, bool);

extern void c_common_mark_addressable_vec (tree);

extern void set_underlying_type (tree);
extern bool user_facing_original_type_p (const_tree);
extern void record_types_used_by_current_var_decl (tree);
extern vec<tree, va_gc> *make_tree_vector (void);
extern void release_tree_vector (vec<tree, va_gc> *);
extern vec<tree, va_gc> *make_tree_vector_single (tree);
extern vec<tree, va_gc> *make_tree_vector_from_list (tree);
extern vec<tree, va_gc> *make_tree_vector_from_ctor (tree);
extern vec<tree, va_gc> *make_tree_vector_copy (const vec<tree, va_gc> *);

/* Used for communication between c_common_type_for_mode and
   c_register_builtin_type.  */
extern GTY(()) tree registered_builtin_types;

/* Read SOURCE_DATE_EPOCH from environment to have a deterministic
   timestamp to replace embedded current dates to get reproducible
   results.  Returns -1 if SOURCE_DATE_EPOCH is not defined.  */
extern time_t cb_get_source_date_epoch (cpp_reader *pfile);

/* The value (as a unix timestamp) corresponds to date
   "Dec 31 9999 23:59:59 UTC", which is the latest date that __DATE__ and
   __TIME__ can store.  */
#define MAX_SOURCE_DATE_EPOCH HOST_WIDE_INT_C (253402300799)

/* Callback for libcpp for offering spelling suggestions for misspelled
   directives.  */
extern const char *cb_get_suggestion (cpp_reader *, const char *,
				      const char *const *);

extern GTY(()) string_concat_db *g_string_concat_db;

class substring_loc;
extern const char *c_get_substring_location (const substring_loc &substr_loc,
					     location_t *out_loc);

/* In c-gimplify.cc.  */
typedef struct bc_state
{
  tree bc_label[2];
} bc_state_t;
extern void save_bc_state (bc_state_t *);
extern void restore_bc_state (bc_state_t *);
extern tree c_genericize_control_stmt (tree *, int *, void *,
				       walk_tree_fn, walk_tree_lh);
extern void c_genericize (tree);
extern int c_gimplify_expr (tree *, gimple_seq *, gimple_seq *);
extern tree c_build_bind_expr (location_t, tree, tree);

/* In c-lex.cc.  */
extern enum cpp_ttype
conflict_marker_get_final_tok_kind (enum cpp_ttype tok1_kind);

/* In c-pch.cc  */
extern void pch_init (void);
extern void pch_cpp_save_state (void);
extern int c_common_valid_pch (cpp_reader *pfile, const char *name, int fd);
extern void c_common_read_pch (cpp_reader *pfile, const char *name, int fd,
			       const char *orig);
extern void c_common_write_pch (void);
extern void c_common_no_more_pch (void);
extern void c_common_pch_pragma (cpp_reader *pfile, const char *);

/* In *-checksum.c */
extern const unsigned char executable_checksum[16];

/* In c-cppbuiltin.cc  */
extern void builtin_define_std (const char *macro);
extern void builtin_define_with_value (const char *, const char *, int);
extern void builtin_define_with_int_value (const char *, HOST_WIDE_INT);
extern void builtin_define_type_sizeof (const char *, tree);
extern void c_stddef_cpp_builtins (void);
extern void fe_file_change (const line_map_ordinary *);
extern void c_parse_error (const char *, enum cpp_ttype, tree, unsigned char,
			   rich_location *richloc);

/* In c-ppoutput.cc  */
extern void init_pp_output (FILE *);
extern void preprocess_file (cpp_reader *);
extern void pp_file_change (const line_map_ordinary *);
extern void pp_dir_change (cpp_reader *, const char *);
extern bool check_missing_format_attribute (tree, tree);

/* In c-omp.cc  */
typedef wide_int_bitmask omp_clause_mask;

#define OMP_CLAUSE_MASK_1 omp_clause_mask (1)

enum c_omp_clause_split
{
  C_OMP_CLAUSE_SPLIT_TARGET = 0,
  C_OMP_CLAUSE_SPLIT_TEAMS,
  C_OMP_CLAUSE_SPLIT_DISTRIBUTE,
  C_OMP_CLAUSE_SPLIT_PARALLEL,
  C_OMP_CLAUSE_SPLIT_FOR,
  C_OMP_CLAUSE_SPLIT_SIMD,
  C_OMP_CLAUSE_SPLIT_COUNT,
  C_OMP_CLAUSE_SPLIT_SECTIONS = C_OMP_CLAUSE_SPLIT_FOR,
  C_OMP_CLAUSE_SPLIT_TASKLOOP = C_OMP_CLAUSE_SPLIT_FOR,
  C_OMP_CLAUSE_SPLIT_LOOP = C_OMP_CLAUSE_SPLIT_FOR,
  C_OMP_CLAUSE_SPLIT_MASKED = C_OMP_CLAUSE_SPLIT_DISTRIBUTE
};

enum c_omp_region_type
{
  C_ORT_OMP			= 1 << 0,
  C_ORT_ACC			= 1 << 1,
  C_ORT_DECLARE_SIMD		= 1 << 2,
  C_ORT_TARGET			= 1 << 3,
  C_ORT_OMP_DECLARE_SIMD	= C_ORT_OMP | C_ORT_DECLARE_SIMD,
  C_ORT_OMP_TARGET		= C_ORT_OMP | C_ORT_TARGET
};

extern tree c_finish_omp_master (location_t, tree);
extern tree c_finish_omp_masked (location_t, tree, tree);
extern tree c_finish_omp_taskgroup (location_t, tree, tree);
extern tree c_finish_omp_critical (location_t, tree, tree, tree);
extern tree c_finish_omp_ordered (location_t, tree, tree);
extern void c_finish_omp_barrier (location_t);
extern tree c_finish_omp_atomic (location_t, enum tree_code, enum tree_code,
				 tree, tree, tree, tree, tree, tree, bool,
				 enum omp_memory_order, bool, bool = false);
extern bool c_omp_depend_t_p (tree);
extern void c_finish_omp_depobj (location_t, tree, enum omp_clause_depend_kind,
				 tree);
extern void c_finish_omp_flush (location_t, int);
extern void c_finish_omp_taskwait (location_t);
extern void c_finish_omp_taskyield (location_t);
extern tree c_finish_omp_for (location_t, enum tree_code, tree, tree, tree,
			      tree, tree, tree, tree, bool);
extern bool c_omp_check_loop_iv (tree, tree, walk_tree_lh);
extern bool c_omp_check_loop_iv_exprs (location_t, enum tree_code, tree, int,
				       tree, tree, tree, walk_tree_lh);
extern tree c_finish_oacc_wait (location_t, tree, tree);
extern tree c_oacc_split_loop_clauses (tree, tree *, bool);
extern void c_omp_split_clauses (location_t, enum tree_code, omp_clause_mask,
				 tree, tree *);
extern tree c_omp_declare_simd_clauses_to_numbers (tree, tree);
extern void c_omp_declare_simd_clauses_to_decls (tree, tree);
extern bool c_omp_predefined_variable (tree);
extern enum omp_clause_default_kind c_omp_predetermined_sharing (tree);
extern enum omp_clause_defaultmap_kind c_omp_predetermined_mapping (tree);
extern tree c_omp_check_context_selector (location_t, tree);
extern void c_omp_mark_declare_variant (location_t, tree, tree);
extern void c_omp_adjust_map_clauses (tree, bool);

enum c_omp_directive_kind {
  C_OMP_DIR_STANDALONE,
  C_OMP_DIR_CONSTRUCT,
  C_OMP_DIR_DECLARATIVE,
  C_OMP_DIR_UTILITY,
  C_OMP_DIR_INFORMATIONAL
};

struct c_omp_directive {
  const char *first, *second, *third;
  unsigned int id;
  enum c_omp_directive_kind kind;
  bool simd;
};

extern const struct c_omp_directive *c_omp_categorize_directive (const char *,
								 const char *,
								 const char *);

/* Return next tree in the chain for chain_next walking of tree nodes.  */
static inline tree
c_tree_chain_next (tree t)
{
  /* TREE_CHAIN of a type is TYPE_STUB_DECL, which is different
     kind of object, never a long chain of nodes.  Prefer
     TYPE_NEXT_VARIANT for types.  */
  if (CODE_CONTAINS_STRUCT (TREE_CODE (t), TS_TYPE_COMMON))
    return TYPE_NEXT_VARIANT (t);
  /* Otherwise, if there is TREE_CHAIN, return it.  */
  if (CODE_CONTAINS_STRUCT (TREE_CODE (t), TS_COMMON))
    return TREE_CHAIN (t);
  return NULL;
}

/* Mask used by tm_stmt_attr.  */
#define TM_STMT_ATTR_OUTER	2
#define TM_STMT_ATTR_ATOMIC	4
#define TM_STMT_ATTR_RELAXED	8

/* Mask used by tm_attr_to_mask and tm_mask_to_attr.  Note that these
   are ordered specifically such that more restrictive attributes are
   at lower bit positions.  This fact is known by the C++ tm attribute
   inheritance code such that least bit extraction (mask & -mask) results
   in the most restrictive attribute.  */
#define TM_ATTR_SAFE			1
#define TM_ATTR_CALLABLE		2
#define TM_ATTR_PURE			4
#define TM_ATTR_IRREVOCABLE		8
#define TM_ATTR_MAY_CANCEL_OUTER	16

/* A suffix-identifier value doublet that represents user-defined literals
   for C++-0x.  */
enum overflow_type {
  OT_UNDERFLOW = -1,
  OT_NONE,
  OT_OVERFLOW
};

struct GTY(()) tree_userdef_literal {
  struct tree_base base;
  tree suffix_id;
  tree value;
  tree num_string;
  enum overflow_type overflow;
};

#define USERDEF_LITERAL_SUFFIX_ID(NODE) \
  (((struct tree_userdef_literal *)USERDEF_LITERAL_CHECK (NODE))->suffix_id)

#define USERDEF_LITERAL_VALUE(NODE) \
  (((struct tree_userdef_literal *)USERDEF_LITERAL_CHECK (NODE))->value)

#define USERDEF_LITERAL_OVERFLOW(NODE) \
  (((struct tree_userdef_literal *)USERDEF_LITERAL_CHECK (NODE))->overflow)

#define USERDEF_LITERAL_NUM_STRING(NODE) \
  (((struct tree_userdef_literal *)USERDEF_LITERAL_CHECK (NODE))->num_string)

#define USERDEF_LITERAL_TYPE(NODE) \
  (TREE_TYPE (USERDEF_LITERAL_VALUE (NODE)))

extern tree build_userdef_literal (tree suffix_id, tree value,
				   enum overflow_type overflow,
				   tree num_string);


/* WHILE_STMT accessors. These give access to the condition of the
   while statement and the body of the while statement, respectively.  */
#define WHILE_COND(NODE)	TREE_OPERAND (WHILE_STMT_CHECK (NODE), 0)
#define WHILE_BODY(NODE)	TREE_OPERAND (WHILE_STMT_CHECK (NODE), 1)

/* DO_STMT accessors. These give access to the condition of the do
   statement and the body of the do statement, respectively.  */
#define DO_COND(NODE)		TREE_OPERAND (DO_STMT_CHECK (NODE), 0)
#define DO_BODY(NODE)		TREE_OPERAND (DO_STMT_CHECK (NODE), 1)

/* FOR_STMT accessors. These give access to the init statement,
   condition, update expression, and body of the for statement,
   respectively.  */
#define FOR_INIT_STMT(NODE)	TREE_OPERAND (FOR_STMT_CHECK (NODE), 0)
#define FOR_COND(NODE)		TREE_OPERAND (FOR_STMT_CHECK (NODE), 1)
#define FOR_EXPR(NODE)		TREE_OPERAND (FOR_STMT_CHECK (NODE), 2)
#define FOR_BODY(NODE)		TREE_OPERAND (FOR_STMT_CHECK (NODE), 3)
#define FOR_SCOPE(NODE)		TREE_OPERAND (FOR_STMT_CHECK (NODE), 4)

#define SWITCH_STMT_COND(NODE)	TREE_OPERAND (SWITCH_STMT_CHECK (NODE), 0)
#define SWITCH_STMT_BODY(NODE)	TREE_OPERAND (SWITCH_STMT_CHECK (NODE), 1)
#define SWITCH_STMT_TYPE(NODE)	TREE_OPERAND (SWITCH_STMT_CHECK (NODE), 2)
#define SWITCH_STMT_SCOPE(NODE)	TREE_OPERAND (SWITCH_STMT_CHECK (NODE), 3)
/* True if there are case labels for all possible values of switch cond, either
   because there is a default: case label or because the case label ranges cover
   all values.  */
#define SWITCH_STMT_ALL_CASES_P(NODE) \
  TREE_LANG_FLAG_0 (SWITCH_STMT_CHECK (NODE))
/* True if the body of a switch stmt contains no BREAK_STMTs.  */
#define SWITCH_STMT_NO_BREAK_P(NODE) \
  TREE_LANG_FLAG_2 (SWITCH_STMT_CHECK (NODE))


/* Nonzero if NODE is the target for genericization of 'break' stmts.  */
#define LABEL_DECL_BREAK(NODE) \
  DECL_LANG_FLAG_0 (LABEL_DECL_CHECK (NODE))

/* Nonzero if NODE is the target for genericization of 'continue' stmts.  */
#define LABEL_DECL_CONTINUE(NODE) \
  DECL_LANG_FLAG_1 (LABEL_DECL_CHECK (NODE))

extern bool convert_vector_to_array_for_subscript (location_t, tree *, tree);

/* Possibe cases of scalar_to_vector conversion.  */
enum stv_conv {
  stv_error,        /* Error occurred.  */
  stv_nothing,      /* Nothing happened.  */
  stv_firstarg,     /* First argument must be expanded.  */
  stv_secondarg     /* Second argument must be expanded.  */
};

extern enum stv_conv scalar_to_vector (location_t loc, enum tree_code code,
				       tree op0, tree op1, bool);

extern tree find_inv_trees (tree *, int *, void *);
extern tree replace_inv_trees (tree *, int *, void *);

extern bool reject_gcc_builtin (const_tree, location_t = UNKNOWN_LOCATION);
extern bool valid_array_size_p (location_t, const_tree, tree, bool = true);
extern void invalid_array_size_error (location_t, cst_size_error,
				      const_tree, const_tree);

/* In c-warn.cc.  */
extern void constant_expression_warning (tree);
extern void constant_expression_error (tree);
extern void overflow_warning (location_t, tree, tree = NULL_TREE);
extern void warn_logical_operator (location_t, enum tree_code, tree,
				   enum tree_code, tree, enum tree_code, tree);
extern void warn_tautological_cmp (const op_location_t &, enum tree_code,
				   tree, tree);
extern void warn_logical_not_parentheses (location_t, enum tree_code, tree,
					  tree);
extern bool warn_if_unused_value (const_tree, location_t, bool = false);
extern bool strict_aliasing_warning (location_t, tree, tree);
extern void sizeof_pointer_memaccess_warning (location_t *, tree,
					      vec<tree, va_gc> *, tree *,
					      bool (*) (tree, tree));
extern void check_main_parameter_types (tree decl);
extern void warnings_for_convert_and_check (location_t, tree, tree, tree);
extern void c_do_switch_warnings (splay_tree, location_t, tree, tree, bool);
extern void warn_for_omitted_condop (location_t, tree);
extern bool warn_for_restrict (unsigned, tree *, unsigned);
extern void warn_for_address_or_pointer_of_packed_member (tree, tree);
extern void warn_parm_array_mismatch (location_t, tree, tree);
extern void maybe_warn_sizeof_array_div (location_t, tree, tree, tree, tree);
extern void do_warn_array_compare (location_t, tree_code, tree, tree);

/* Places where an lvalue, or modifiable lvalue, may be required.
   Used to select diagnostic messages in lvalue_error and
   readonly_error.  */
enum lvalue_use {
  lv_assign,
  lv_increment,
  lv_decrement,
  lv_addressof,
  lv_asm
};

extern void lvalue_error (location_t, enum lvalue_use);
extern void invalid_indirection_error (location_t, tree, ref_operator);
extern void readonly_error (location_t, tree, enum lvalue_use);
extern void warn_array_subscript_with_type_char (location_t, tree);
extern void warn_about_parentheses (location_t,
				    enum tree_code,
				    enum tree_code, tree,
				    enum tree_code, tree);
extern void warn_for_unused_label (tree label);
extern void warn_for_div_by_zero (location_t, tree divisor);
extern void warn_for_memset (location_t, tree, tree, int);
extern void warn_for_sign_compare (location_t,
				   tree orig_op0, tree orig_op1,
				   tree op0, tree op1,
				   tree result_type,
				   enum tree_code resultcode);
extern void do_warn_double_promotion (tree, tree, tree, const char *,
				      location_t);
extern void do_warn_unused_parameter (tree);
extern void record_locally_defined_typedef (tree);
extern void maybe_record_typedef_use (tree);
extern void maybe_warn_unused_local_typedefs (void);
extern void maybe_warn_bool_compare (location_t, enum tree_code, tree, tree);
extern bool maybe_warn_shift_overflow (location_t, tree, tree);
extern void warn_duplicated_cond_add_or_warn (location_t, tree, vec<tree> **);
extern bool diagnose_mismatched_attributes (tree, tree);
extern tree do_warn_duplicated_branches_r (tree *, int *, void *);
extern void warn_for_multistatement_macros (location_t, location_t,
					    location_t, enum rid);

/* In c-attribs.cc.  */
extern bool attribute_takes_identifier_p (const_tree);
extern tree handle_deprecated_attribute (tree *, tree, tree, int, bool *);
extern tree handle_unused_attribute (tree *, tree, tree, int, bool *);
extern tree handle_fallthrough_attribute (tree *, tree, tree, int, bool *);
extern int parse_tm_stmt_attr (tree, int);
extern int tm_attr_to_mask (tree);
extern tree tm_mask_to_attr (int);
extern tree find_tm_attribute (tree);
extern const struct attribute_spec::exclusions attr_cold_hot_exclusions[];
extern const struct attribute_spec::exclusions attr_noreturn_exclusions[];
extern tree handle_noreturn_attribute (tree *, tree, tree, int, bool *);
extern bool has_attribute (location_t, tree, tree, tree (*)(tree));
extern tree build_attr_access_from_parms (tree, bool);

/* In c-format.cc.  */
extern bool valid_format_string_type_p (tree);

/* A bitmap of flags to positional_argument.  */
enum posargflags {
  /* Consider positional attribute argument value zero valid.  */
  POSARG_ZERO = 1,
  /* Consider positional attribute argument value valid if it refers
     to the ellipsis (i.e., beyond the last typed argument).  */
  POSARG_ELLIPSIS = 2
};

extern tree positional_argument (const_tree, const_tree, tree, tree_code,
				 int = 0, int = posargflags ());

extern enum flt_eval_method
excess_precision_mode_join (enum flt_eval_method, enum flt_eval_method);

extern int c_flt_eval_method (bool ts18661_p);
extern void add_no_sanitize_value (tree node, unsigned int flags);

extern void maybe_add_include_fixit (rich_location *, const char *, bool);
extern void maybe_suggest_missing_token_insertion (rich_location *richloc,
						   enum cpp_ttype token_type,
						   location_t prev_token_loc);
extern tree braced_lists_to_strings (tree, tree);

#if CHECKING_P
namespace selftest {
  /* Declarations for specific families of tests within c-family,
     by source file, in alphabetical order.  */
  extern void c_diagnostic_cc_tests (void);
  extern void c_format_cc_tests (void);
  extern void c_indentation_cc_tests (void);
  extern void c_opt_problem_cc_tests (void);
  extern void c_pretty_print_cc_tests (void);
  extern void c_spellcheck_cc_tests (void);

  /* The entrypoint for running all of the above tests.  */
  extern void c_family_tests (void);
} // namespace selftest
#endif /* #if CHECKING_P */

#endif /* ! GCC_C_COMMON_H */
