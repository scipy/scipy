/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: parser.h                                                          *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for parser                                    *
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

=NODE  StringAPI   String Interface
=UP TOP [25]

=DESCRIPTION

   The string interface (string API) provided by the
   unur_str2gen() call is the easiest way to use UNU.RAN. This
   function takes a character string as its argument. The string is
   parsed and the information obtained is used to create a generator
   object. It returns NULL if this fails, either due to a syntax
   error, or due to invalid data. In both cases @code{unur_error} is
   set to the corresponding error codes
   (@pxref{Error_reporting,,Error reporting}).
   Additionally there exists the call unur_str2distr() that only
   produces a distribution object.

   Notice that the string interface does not implement all features of
   the UNU.RAN library. For trickier tasks it might be necessary to use
   the UNU.RAN calls.

   In @ref{Examples}, all examples are given using both the 
   UNU.RAN standard API and this convenient string API.
   The corresponding programm codes are equivalent.

=END

*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_GEN *unur_str2gen( const char *string );
/* 
   Get a generator object for the distribution, method and uniform
   random number generator as described in the given @var{string}.
   See @ref{StringSyntax,,Syntax of String Interface}, for details.
*/

UNUR_DISTR *unur_str2distr( const char *string );
/*
   Get a distribution object for the distribution described in
   @var{string}. 
   See @ref{StringSyntax,,Syntax of String Interface},
   and @ref{StringDistr,,Distribution String},
   for details. However, only the block for the distribution object is
   allowed.
*/

UNUR_GEN *unur_makegen_ssu( const char *distrstr, const char *methodstr, UNUR_URNG *urng );
/* */

UNUR_GEN *unur_makegen_dsu( const UNUR_DISTR *distribution, const char *methodstr, UNUR_URNG *urng );
/* 
   Make a generator object for the distribution, method and uniform
   random number generator. The distribution can be given either as
   string @var{distrstr} or as a distribution object @var{distr}.
   The method must be given as a string @var{methodstr}.
   For the syntax of these strings see 
   @ref{StringSyntax,,Syntax of String Interface}.
   However, the @code{method} keyword is optional for these calls
   and can be omitted. If @var{methodstr} is the empty (blank) string
   or NULL method AUTO is used.
   The uniform random number generator is optional. If 
   @var{urng} is NULL then the default uniform random number generator
   is used.
*/

/*
=EON
*/

/*...........................................................................*/


UNUR_PAR *_unur_str2par( const UNUR_DISTR *distribution, const char *method, struct unur_slist **mlist );
/*
   Get a parameter object for the given distribution object and for the 
   method described in @var{method}. 
   See @ref{StringSyntax,,Syntax of String Interface},
   and @ref{StringMethod,,Method String},
   for details. However, only the block for the method is allowed.

   @var{mlist} is used to store memory blocks that have been allocated
   while parsing and making the parameter object. These blocks
   should (must) be free using _unur_slist_free() after initializing
   the generator object to avoid a memory leak.

   (Not part of manual!!)
*/

/*---------------------------------------------------------------------------*/
/*

=NODE  StringSyntax   Syntax of String Interface
=UP StringAPI [10]

=DESCRIPTION

   The given string holds information about the requested distribution
   and (optional) about the sampling method and the uniform random
   number generator invoked. The interpretation of the string is not
   case-sensitive, all white spaces are ignored.

   The string consists of up to three blocks, separated by ampersands
   @code{&}.

   Each block consists of @code{<key>=<value>} pairs, separated by
   semicolons @code{;}.

   The first key in each block is used to indicate each block.
   We have three different blocks with the following (first) keys:
   @table @code 
   @item distr
      definition of the distribution
      (@pxref{StringDistr,,Distribution String}).

   @item method
      description of the transformation method
      (@pxref{StringMethod,,Method String}).

   @item urng
      uniform random number generation
      (@pxref{StringURNG,,Uniform RNG String}).
   @end table

   The @code{distr} block must be the very first block and is
   obligatory. All the other blocks are optional and can be arranged
   in arbitrary order.

   For details see the following description of each block.

   In the following example
   @smallexample
   distr = normal(3.,0.75); domain = (0,inf) & method = tdr; c = 0 
   @end smallexample
   we have a distribution block for the truncated normal distribution
   with mean 3 and standard deviation 0.75 on domain (0,infinity);
   and block for choosing method TDR with parameter c set to 0.

   @sp 1
   The @code{<key>=<value>} pairs that follow the first (initial) pair
   in each block are used to set parameters.
   The name of the parameter is given by the @code{<key>} string. It is
   deduced from the UNU.RAN set calls by taking the part after
   @code{@dots{}_set_}.
   The @code{<value>} string holds the parameters to be
   set, separated by commata @code{,}.
   There are three types of parameters:
   @table @emph
   @item string @code{"@dots{}"}
      i.e. any sequence of characters enclosed by double quotes 
      @code{"@dots{}"}, 
   @item list @code{(@dots{},@dots{})}
      i.e. list of @emph{numbers}, separated by commata @code{,},
      enclosed in parenthesis @code{(...)}, and
   @item number
      a sequence of characters that is not enclosed by quotes
      @code{"@dots{}"} or parenthesis @code{(...)}. 
      It is interpreted as float or integer depending on the type of 
      the corresponding parameter.
   @end table
   The @code{<value>} string (including the character @code{=}) can be
   omitted when no argument is required.

   At the moment not all @command{set} calls are supported.
   The syntax for the @code{<value>} can be directly derived from the
   corresponding @command{set} calls. To simplify the syntax additional
   shortcuts are possible. The following table lists the parameters for
   the @code{set} calls that are supported by the string interface; the
   entry in parenthesis gives the type of the argument as 
   @code{<value>} string:

   @table @code
   @item int  @i{(number)}:
      The @i{number} is interpreted as an integer.
      @code{true} and @code{on} are transformed to @code{1},
      @code{false} and @code{off} are transformed to @code{0}.
      A missing argument is interpreted as @code{1}.

   @item int, int  @i{(number, number} @r{or} @i{list)}:
      The two numbers or the first two entries in the list are
      interpreted as a integers.
      @code{inf} and @code{-inf} are transformed to @code{INT_MAX} and
      @code{INT_MIN} respectively, i.e. the largest and smallest
      integers that can be represented by the computer.

   @item unsigned @i{(number)}: 
      The @i{number} is interpreted as an unsigned hexadecimal
      integer.

   @item double  @i{(number)}:
      The number is interpreted as a floating point number.
      @code{inf} is transformed to @code{UNUR_INFINITY}.

   @item double, double  @i{(number, number} @r{or} @i{list)}:
      The two numbers or the first two entries in the list are
      interpreted as a floating point numbers.
      @code{inf} is transformed to @code{UNUR_INFINITY}. However using
      @code{inf} in the list might not work for all versions of C. Then it
      is recommended to use two single numbers instead of a list.

   @item int, double*  @i{([number,] list} @r{or} @i{number)}:
      @itemize @minus
      @item
         The list is interpreted as a double array.
	 The (first) number as its length.
	 If it is less than the actual size of the array only the
	 first entries of the array are used.
      @item
         If only the list is given (i.e., if the first number is omitted),
	 the first number is set to the actual size of the array.
      @item
	 If only the number is given (i.e., if the list is omitted), the NULL
	 pointer is used instead an array as argument.
      @end itemize

   @item double*, int  @i{(list [,number])}:
      The list is interpreted as a double array.
      The (second) number as its length. 
      If the length is omitted, it is replaced by the actual size of the
      array. (Only in the @code{distribution} block!)

   @item char*  @i{(string)}:
      The character string is passed as is to the corresponding set
      call.

   @end table

   Notice that missing entries in a list of numbers are interpreted as 
   @code{0}. E.g, a the list @code{(1,,3)} is read as @code{(1,0,3)}, the
   list @code{(1,2,)} as @code{(1,2,0)}.

   The the list of @code{key} strings in
   @ref{KeysDistr,,Keys for Distribution String}, and 
   @ref{KeysMethod,,Keys for Method String}, for further details.

=EON */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*

=NODEX  StringDistr    Distribution String
=UP StringAPI [20]

=DESCRIPTION

   The @code{distr} block must be the very first block and is
   obligatory. For that reason the keyword @code{distr} is optional and
   can be omitted (together with the @code{=} character). 
   Moreover it is ignored while parsing the string. However, to
   avoid some possible confusion it has to start with the
   letter @code{d} (if it is given at all). 

   The value of the @code{distr} key is used to get the distribution
   object, either via a @command{unur_distr_<value>} call for a standard
   distribution via a @command{unur_distr_<value>_new} call to get an
   object of a generic distribution. 
   However not all generic distributions are supported yet.

   The parameters for the standard distribution are given
   as a list. There must not be any character (other than white space)
   between the name of the standard distribution and the opening
   parenthesis @code{(} of this list. E.g., to get a beta distribution,
   use
   @smallexample
      distr = beta(2,4)
   @end smallexample

   To get an object for a discrete distribution with probability
   vector (0.5,0.2,0.3), use 
   @smallexample
      distr = discr; pv = (0.5,0.2,0.3)
   @end smallexample

   It is also possible to set a PDF, PMF, or CDF using a string.
   E.g., to create a continuous distribution with PDF proportional to 
   @code{exp(-sqrt(2+(x-1)^2) + (x-1))} and domain (0,inf) use
   @smallexample
      distr = cont; pdf = "exp(-sqrt(2+(x-1)^2) + (x-1))"
   @end smallexample
   Notice: If this string is used in an unur_str2distr() or
   unur_str2gen() call the double quotes @code{"} must be protected by
   @code{\"}. Alternatively, single quotes may be used instead 
   @smallexample
      distr = cont; pdf = 'exp(-sqrt(2+(x-1)^2) + (x-1))'
   @end smallexample

   For the details of function strings see
   @ref{StringFunct,,Function String}.

=EON
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*

=NODEX  StringMethod   Method String
=UP StringAPI [30]

=DESCRIPTION

   The key @code{method} is obligatory, it must be the first key and its
   value is the name of a method suitable for the choosen standard
   distribution. E.g., if method AROU is chosen, use
   @smallexample
      method = arou
   @end smallexample

   Of course the all following keys dependend on the method choosen at
   first. All corresponding @command{set} calls of UNU.RAN are available
   and the key is the string after the @command{unur_<methodname>_set_}
   part of the command. E.g., UNU.RAN provides the command 
   @command{unur_arou_set_max_sqhratio} to set a parameter of method AROU.
   To call this function via the string-interface, the
   key @code{max_sqhratio} can be used:
   @smallexample
      max_sqhratio = 0.9
   @end smallexample
   Additionally the keyword @code{debug} can be used to set debugging
   flags (see @ref{Debug,,Debugging}, for details).
   
   If this block is omitted, a suitable default method is used. Notice
   however that the default method may change in future versions of
   UNU.RAN.

=EON
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*

=NODEX  StringURNG    Uniform RNG String
=UP StringAPI [40]

=DESCRIPTION

   The value of the @code{urng} key is passed to the PRNG interface (see
   @ifinfo
      @xref{Top,,Overview,prng,PRNG Manual}.
   @end ifinfo
   @ifnotinfo
      @uref{http://statmath.wu.ac.at/prng/manual/,PRNG manual}
   @end ifnotinfo
   for details).
   However it only works when using the PRNG library is enabled, 
   see @ref{Installation} for details. There are no other keys.

   IMPORTANT: UNU.RAN creates a new uniform random number generator for
   the generator object. The pointer to this uniform generator 
   has to be read and saved via a unur_get_urng() call in order to
   clear the memory @emph{before} the UNU.RAN generator object is
   destroyed.

   If this block is omitted the UNU.RAN default generator is used
   (which @emph{must not} be destroyed). 

=EON
*/
/*---------------------------------------------------------------------------*/
