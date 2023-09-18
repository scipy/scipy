/* Map (unsigned int) keys to (source file, line, column) triples.
   Copyright (C) 2001-2022 Free Software Foundation, Inc.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3, or (at your option) any
later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.

 In other words, you are welcome to use, share and improve this program.
 You are forbidden to forbid anyone else to use, share and improve
 what you give them.   Help stamp out software-hoarding!  */

#ifndef LIBCPP_LINE_MAP_H
#define LIBCPP_LINE_MAP_H

#ifndef GTY
#define GTY(x) /* nothing */
#endif

/* Both gcc and emacs number source *lines* starting at 1, but
   they have differing conventions for *columns*.

   GCC uses a 1-based convention for source columns,
   whereas Emacs's M-x column-number-mode uses a 0-based convention.

   For example, an error in the initial, left-hand
   column of source line 3 is reported by GCC as:

      some-file.c:3:1: error: ...etc...

   On navigating to the location of that error in Emacs
   (e.g. via "next-error"),
   the locus is reported in the Mode Line
   (assuming M-x column-number-mode) as:

     some-file.c   10%   (3, 0)

   i.e. "3:1:" in GCC corresponds to "(3, 0)" in Emacs.  */

/* The type of line numbers.  */
typedef unsigned int linenum_type;

/* A type for doing arithmetic on line numbers.  */
typedef long long linenum_arith_t;

/* A function for for use by qsort for comparing line numbers.  */

inline int compare (linenum_type lhs, linenum_type rhs)
{
  /* Avoid truncation issues by using linenum_arith_t for the comparison,
     and only consider the sign of the result.  */
  linenum_arith_t diff = (linenum_arith_t)lhs - (linenum_arith_t)rhs;
  if (diff)
    return diff > 0 ? 1 : -1;
  return 0;
}

/* Reason for creating a new line map with linemap_add.  */
enum lc_reason
{
  LC_ENTER = 0,		/* Begin #include.  */
  LC_LEAVE,		/* Return to including file.  */
  LC_RENAME,		/* Other reason for name change.  */
  LC_RENAME_VERBATIM,	/* Likewise, but "" != stdin.  */
  LC_ENTER_MACRO,	/* Begin macro expansion.  */
  LC_MODULE,		/* A (C++) Module.  */
  /* FIXME: add support for stringize and paste.  */
  LC_HWM /* High Water Mark.  */
};

/* The typedef "location_t" is a key within the location database,
   identifying a source location or macro expansion, along with range
   information, and (optionally) a pointer for use by gcc.

   This key only has meaning in relation to a line_maps instance.  Within
   gcc there is a single line_maps instance: "line_table", declared in
   gcc/input.h and defined in gcc/input.cc.

   The values of the keys are intended to be internal to libcpp,
   but for ease-of-understanding the implementation, they are currently
   assigned as follows:

  Actual     | Value                         | Meaning
  -----------+-------------------------------+-------------------------------
  0x00000000 | UNKNOWN_LOCATION (gcc/input.h)| Unknown/invalid location.
  -----------+-------------------------------+-------------------------------
  0x00000001 | BUILTINS_LOCATION             | The location for declarations
             |   (gcc/input.h)               | in "<built-in>"
  -----------+-------------------------------+-------------------------------
  0x00000002 | RESERVED_LOCATION_COUNT       | The first location to be
             | (also                         | handed out, and the
             |  ordmap[0]->start_location)   | first line in ordmap 0
  -----------+-------------------------------+-------------------------------
             | ordmap[1]->start_location     | First line in ordmap 1
             | ordmap[1]->start_location+32  | First column in that line
             |   (assuming range_bits == 5)  |
             | ordmap[1]->start_location+64  | 2nd column in that line
             | ordmap[1]->start_location+4096| Second line in ordmap 1
             |   (assuming column_bits == 12)
             |
             |   Subsequent lines are offset by (1 << column_bits),
             |   e.g. 4096 for 12 bits, with a column value of 0 representing
             |   "the whole line".
             |
             |   Within a line, the low "range_bits" (typically 5) are used for
             |   storing short ranges, so that there's an offset of
             |     (1 << range_bits) between individual columns within a line,
             |   typically 32.
             |   The low range_bits store the offset of the end point from the
             |   start point, and the start point is found by masking away
             |   the range bits.
             |
             |   For example:
             |      ordmap[1]->start_location+64    "2nd column in that line"
             |   above means a caret at that location, with a range
             |   starting and finishing at the same place (the range bits
             |   are 0), a range of length 1.
             |
             |   By contrast:
             |      ordmap[1]->start_location+68
             |   has range bits 0x4, meaning a caret with a range starting at
             |   that location, but with endpoint 4 columns further on: a range
             |   of length 5.
             |
             |   Ranges that have caret != start, or have an endpoint too
             |   far away to fit in range_bits are instead stored as ad-hoc
             |   locations.  Hence for range_bits == 5 we can compactly store
             |   tokens of length <= 32 without needing to use the ad-hoc
             |   table.
             |
             |   This packing scheme means we effectively have
             |     (column_bits - range_bits)
             |   of bits for the columns, typically (12 - 5) = 7, for 128
             |   columns; longer line widths are accomodated by starting a
             |   new ordmap with a higher column_bits.
             |
             | ordmap[2]->start_location-1   | Final location in ordmap 1
  -----------+-------------------------------+-------------------------------
             | ordmap[2]->start_location     | First line in ordmap 2
             | ordmap[3]->start_location-1   | Final location in ordmap 2
  -----------+-------------------------------+-------------------------------
             |                               | (etc)
  -----------+-------------------------------+-------------------------------
             | ordmap[n-1]->start_location   | First line in final ord map
             |                               | (etc)
             | set->highest_location - 1     | Final location in that ordmap
  -----------+-------------------------------+-------------------------------
             | set->highest_location         | Location of the where the next
             |                               | ordinary linemap would start
  -----------+-------------------------------+-------------------------------
             |                               |
             |                  VVVVVVVVVVVVVVVVVVVVVVVVVVV
             |                  Ordinary maps grow this way
             |
             |                    (unallocated integers)
             |
  0x60000000 | LINE_MAP_MAX_LOCATION_WITH_COLS
             |   Beyond this point, ordinary linemaps have 0 bits per column:
             |   each increment of the value corresponds to a new source line.
             |
  0x70000000 | LINE_MAP_MAX_LOCATION
             |   Beyond the point, we give up on ordinary maps; attempts to
             |   create locations in them lead to UNKNOWN_LOCATION (0).
             |
             |                    (unallocated integers)
             |
             |                   Macro maps grow this way
             |                   ^^^^^^^^^^^^^^^^^^^^^^^^
             |                               |
  -----------+-------------------------------+-------------------------------
             | LINEMAPS_MACRO_LOWEST_LOCATION| Locations within macro maps
             | macromap[m-1]->start_location | Start of last macro map
             |                               |
  -----------+-------------------------------+-------------------------------
             | macromap[m-2]->start_location | Start of penultimate macro map
  -----------+-------------------------------+-------------------------------
             | macromap[1]->start_location   | Start of macro map 1
  -----------+-------------------------------+-------------------------------
             | macromap[0]->start_location   | Start of macro map 0
  0x7fffffff | MAX_LOCATION_T                | Also used as a mask for
             |                               | accessing the ad-hoc data table
  -----------+-------------------------------+-------------------------------
  0x80000000 | Start of ad-hoc values; the lower 31 bits are used as an index
  ...        | into the line_table->location_adhoc_data_map.data array.
  0xffffffff | UINT_MAX                      |
  -----------+-------------------------------+-------------------------------

   Examples of location encoding.

   Packed ranges
   =============

   Consider encoding the location of a token "foo", seen underlined here
   on line 523, within an ordinary line_map that starts at line 500:

                 11111111112
        12345678901234567890
     522
     523   return foo + bar;
                  ^~~
     524

   The location's caret and start are both at line 523, column 11; the
   location's finish is on the same line, at column 13 (an offset of 2
   columns, for length 3).

   Line 523 is offset 23 from the starting line of the ordinary line_map.

   caret == start, and the offset of the finish fits within 5 bits, so
   this can be stored as a packed range.

   This is encoded as:
      ordmap->start
         + (line_offset << ordmap->m_column_and_range_bits)
         + (column << ordmap->m_range_bits)
         + (range_offset);
   i.e. (for line offset 23, column 11, range offset 2):
      ordmap->start
         + (23 << 12)
         + (11 << 5)
         + 2;
   i.e.:
      ordmap->start + 0x17162
   assuming that the line_map uses the default of 7 bits for columns and
   5 bits for packed range (giving 12 bits for m_column_and_range_bits).


   "Pure" locations
   ================

   These are a special case of the above, where
      caret == start == finish
   They are stored as packed ranges with offset == 0.
   For example, the location of the "f" of "foo" could be stored
   as above, but with range offset 0, giving:
      ordmap->start
         + (23 << 12)
         + (11 << 5)
         + 0;
   i.e.:
      ordmap->start + 0x17160


   Unoptimized ranges
   ==================

   Consider encoding the location of the binary expression
   below:

                 11111111112
        12345678901234567890
     522
     523   return foo + bar;
                  ~~~~^~~~~
     524

   The location's caret is at the "+", line 523 column 15, but starts
   earlier, at the "f" of "foo" at column 11.  The finish is at the "r"
   of "bar" at column 19.

   This can't be stored as a packed range since start != caret.
   Hence it is stored as an ad-hoc location e.g. 0x80000003.

   Stripping off the top bit gives us an index into the ad-hoc
   lookaside table:

     line_table->location_adhoc_data_map.data[0x3]

   from which the caret, start and finish can be looked up,
   encoded as "pure" locations:

     start  == ordmap->start + (23 << 12) + (11 << 5)
            == ordmap->start + 0x17160  (as above; the "f" of "foo")

     caret  == ordmap->start + (23 << 12) + (15 << 5)
            == ordmap->start + 0x171e0

     finish == ordmap->start + (23 << 12) + (19 << 5)
            == ordmap->start + 0x17260

   To further see how location_t works in practice, see the
   worked example in libcpp/location-example.txt.  */
typedef unsigned int location_t;

/* Do not track column numbers higher than this one.  As a result, the
   range of column_bits is [12, 18] (or 0 if column numbers are
   disabled).  */
const unsigned int LINE_MAP_MAX_COLUMN_NUMBER = (1U << 12);

/* Do not pack ranges if locations get higher than this.
   If you change this, update:
     gcc.dg/plugin/location-overflow-test-*.c.  */
const location_t LINE_MAP_MAX_LOCATION_WITH_PACKED_RANGES = 0x50000000;

/* Do not track column numbers if locations get higher than this.
   If you change this, update:
     gcc.dg/plugin/location-overflow-test-*.c.  */
const location_t LINE_MAP_MAX_LOCATION_WITH_COLS = 0x60000000;

/* Highest possible source location encoded within an ordinary map.  */
const location_t LINE_MAP_MAX_LOCATION = 0x70000000;

/* A range of source locations.

   Ranges are closed:
   m_start is the first location within the range,
   m_finish is the last location within the range.

   We may need a more compact way to store these, but for now,
   let's do it the simple way, as a pair.  */
struct GTY(()) source_range
{
  location_t m_start;
  location_t m_finish;

  /* We avoid using constructors, since various structs that
     don't yet have constructors will embed instances of
     source_range.  */

  /* Make a source_range from a location_t.  */
  static source_range from_location (location_t loc)
  {
    source_range result;
    result.m_start = loc;
    result.m_finish = loc;
    return result;
  }

  /* Make a source_range from a pair of location_t.  */
  static source_range from_locations (location_t start,
				      location_t finish)
  {
    source_range result;
    result.m_start = start;
    result.m_finish = finish;
    return result;
  }
};

/* Memory allocation function typedef.  Works like xrealloc.  */
typedef void *(*line_map_realloc) (void *, size_t);

/* Memory allocator function that returns the actual allocated size,
   for a given requested allocation.  */
typedef size_t (*line_map_round_alloc_size_func) (size_t);

/* A line_map encodes a sequence of locations.
   There are two kinds of maps. Ordinary maps and macro expansion
   maps, a.k.a macro maps.

   A macro map encodes source locations of tokens that are part of a
   macro replacement-list, at a macro expansion point. E.g, in:

            #define PLUS(A,B) A + B

   No macro map is going to be created there, because we are not at a
   macro expansion point. We are at a macro /definition/ point. So the
   locations of the tokens of the macro replacement-list (i.e, A + B)
   will be locations in an ordinary map, not a macro map.

   On the other hand, if we later do:

        int a = PLUS (1,2);

   The invocation of PLUS here is a macro expansion. So we are at a
   macro expansion point. The preprocessor expands PLUS (1,2) and
   replaces it with the tokens of its replacement-list: 1 + 2. A macro
   map is going to be created to hold (or rather to map, haha ...) the
   locations of the tokens 1, + and 2. The macro map also records the
   location of the expansion point of PLUS. That location is mapped in
   the map that is active right before the location of the invocation
   of PLUS.  */

/* This contains GTY mark-up to support precompiled headers.
   line_map is an abstract class, only derived objects exist.  */
struct GTY((tag ("0"), desc ("MAP_ORDINARY_P (&%h) ? 1 : 2"))) line_map {
  location_t start_location;

  /* Size and alignment is (usually) 4 bytes.  */
};

/* An ordinary line map encodes physical source locations. Those
   physical source locations are called "spelling locations".
   
   Physical source file TO_FILE at line TO_LINE at column 0 is represented
   by the logical START_LOCATION.  TO_LINE+L at column C is represented by
   START_LOCATION+(L*(1<<m_column_and_range_bits))+(C*1<<m_range_bits), as
   long as C<(1<<effective range bits), and the result_location is less than
   the next line_map's start_location.
   (The top line is line 1 and the leftmost column is column 1; line/column 0
   means "entire file/line" or "unknown line/column" or "not applicable".)

   The highest possible source location is MAX_LOCATION_T.  */
struct GTY((tag ("1"))) line_map_ordinary : public line_map {
  /* Base class is 4 bytes.  */

  /* 4 bytes of integers, each 1 byte for easy extraction/insertion.  */

  /* The reason for creation of this line map.  */
  ENUM_BITFIELD (lc_reason) reason : 8;

  /* SYSP is one for a system header, two for a C system header file
     that therefore needs to be extern "C" protected in C++, and zero
     otherwise.  This field isn't really needed now that it's in
     cpp_buffer.  */
  unsigned char sysp;

  /* Number of the low-order location_t bits used for column numbers
     and ranges.  */
  unsigned int m_column_and_range_bits : 8;

  /* Number of the low-order "column" bits used for storing short ranges
     inline, rather than in the ad-hoc table.
     MSB                                                                 LSB
     31                                                                    0
     +-------------------------+-------------------------------------------+
     |                         |<---map->column_and_range_bits (e.g. 12)-->|
     +-------------------------+-----------------------+-------------------+
     |                         | column_and_range_bits | map->range_bits   |
     |                         |   - range_bits        |                   |
     +-------------------------+-----------------------+-------------------+
     | row bits                | effective column bits | short range bits  |
     |                         |    (e.g. 7)           |   (e.g. 5)        |
     +-------------------------+-----------------------+-------------------+ */
  unsigned int m_range_bits : 8;

  /* Pointer alignment boundary on both 32 and 64-bit systems.  */

  const char *to_file;
  linenum_type to_line;

  /* Location from whence this line map was included.  For regular
     #includes, this location will be the last location of a map.  For
     outermost file, this is 0.  For modules it could be anywhere
     within a map.  */
  location_t included_from;

  /* Size is 20 or 24 bytes, no padding  */
};

/* This is the highest possible source location encoded within an
   ordinary or macro map.  */
const location_t MAX_LOCATION_T = 0x7FFFFFFF;

struct cpp_hashnode;

/* A macro line map encodes location of tokens coming from a macro
   expansion.
   
   The offset from START_LOCATION is used to index into
   MACRO_LOCATIONS; this holds the original location of the token.  */
struct GTY((tag ("2"))) line_map_macro : public line_map {
  /* Base is 4 bytes.  */

  /* The number of tokens inside the replacement-list of MACRO.  */
  unsigned int n_tokens;

  /* Pointer alignment boundary.  */

  /* The cpp macro whose expansion gave birth to this macro map.  */
  struct cpp_hashnode *
    GTY ((nested_ptr (union tree_node,
		      "%h ? CPP_HASHNODE (GCC_IDENT_TO_HT_IDENT (%h)) : NULL",
		      "%h ? HT_IDENT_TO_GCC_IDENT (HT_NODE (%h)) : NULL")))
    macro;

  /* This array of location is actually an array of pairs of
     locations. The elements inside it thus look like:

           x0,y0, x1,y1, x2,y2, ...., xn,yn.

     where n == n_tokens;

     Remember that these xI,yI are collected when libcpp is about to
     expand a given macro.

     yI is the location in the macro definition, either of the token
     itself or of a macro parameter that it replaces.

     Imagine this:

	#define PLUS(A, B) A + B  <--- #1

	int a = PLUS (1,2); <--- #2

     There is a macro map for the expansion of PLUS in #2.  PLUS is
     expanded into its expansion-list.  The expansion-list is the
     replacement-list of PLUS where the macro parameters are replaced
     with their arguments.  So the replacement-list of PLUS is made of
     the tokens:

        A, +, B

     and the expansion-list is made of the tokens:

        1, +, 2

     Let's consider the case of token "+".  Its y1 [yI for I == 1] is
     its spelling location in #1.

     y0 (thus for token "1") is the spelling location of A in #1.

     And y2 (of token "2") is the spelling location of B in #1.

     When the token is /not/ an argument for a macro, xI is the same
     location as yI.  Otherwise, xI is the location of the token
     outside this macro expansion.  If this macro was expanded from
     another macro expansion, xI is a virtual location representing
     the token in that macro expansion; otherwise, it is the spelling
     location of the token.

     Note that a virtual location is a location returned by
     linemap_add_macro_token.  It encodes the relevant locations (x,y
     pairs) of that token across the macro expansions from which it
     (the token) might come from.

     In the example above x1 (for token "+") is going to be the same
     as y1.  x0 is the spelling location for the argument token "1",
     and x2 is the spelling location for the argument token "2".  */
  location_t * GTY((atomic)) macro_locations;

  /* This is the location of the expansion point of the current macro
     map.  It's the location of the macro name.  That location is held
     by the map that was current right before the current one. It
     could have been either a macro or an ordinary map, depending on
     if we are in a nested expansion context not.  */
  location_t expansion;

  /* Size is 20 or 32 (4 bytes padding on 64-bit).  */
};

#if CHECKING_P && (GCC_VERSION >= 2007)

/* Assertion macro to be used in line-map code.  */
#define linemap_assert(EXPR)                  \
  do {                                                \
    if (! (EXPR))                             \
      abort ();                                       \
  } while (0)

/* Assert that becomes a conditional expression when checking is disabled at
   compilation time.  Use this for conditions that should not happen but if
   they happen, it is better to handle them gracefully rather than crash
   randomly later.
   Usage:

   if (linemap_assert_fails(EXPR)) handle_error(); */
#define linemap_assert_fails(EXPR) __extension__ \
  ({linemap_assert (EXPR); false;})

#else
/* Include EXPR, so that unused variable warnings do not occur.  */
#define linemap_assert(EXPR) ((void)(0 && (EXPR)))
#define linemap_assert_fails(EXPR) (! (EXPR))
#endif

/* Get whether location LOC is an ordinary location.  */

inline bool
IS_ORDINARY_LOC (location_t loc)
{
  return loc < LINE_MAP_MAX_LOCATION;
}

/* Get whether location LOC is an ad-hoc location.  */

inline bool
IS_ADHOC_LOC (location_t loc)
{
  return loc > MAX_LOCATION_T;
}

/* Categorize line map kinds.  */

inline bool
MAP_ORDINARY_P (const line_map *map)
{
  return IS_ORDINARY_LOC (map->start_location);
}

/* Return TRUE if MAP encodes locations coming from a macro
   replacement-list at macro expansion point.  */
bool
linemap_macro_expansion_map_p (const line_map *);

/* Assert that MAP encodes locations of tokens that are not part of
   the replacement-list of a macro expansion, downcasting from
   line_map * to line_map_ordinary *.  */

inline line_map_ordinary *
linemap_check_ordinary (line_map *map)
{
  linemap_assert (MAP_ORDINARY_P (map));
  return (line_map_ordinary *)map;
}

/* Assert that MAP encodes locations of tokens that are not part of
   the replacement-list of a macro expansion, downcasting from
   const line_map * to const line_map_ordinary *.  */

inline const line_map_ordinary *
linemap_check_ordinary (const line_map *map)
{
  linemap_assert (MAP_ORDINARY_P (map));
  return (const line_map_ordinary *)map;
}

/* Assert that MAP is a macro expansion and downcast to the appropriate
   subclass.  */

inline line_map_macro *linemap_check_macro (line_map *map)
{
  linemap_assert (!MAP_ORDINARY_P (map));
  return (line_map_macro *)map;
}

/* Assert that MAP is a macro expansion and downcast to the appropriate
   subclass.  */

inline const line_map_macro *
linemap_check_macro (const line_map *map)
{
  linemap_assert (!MAP_ORDINARY_P (map));
  return (const line_map_macro *)map;
}

/* Read the start location of MAP.  */

inline location_t
MAP_START_LOCATION (const line_map *map)
{
  return map->start_location;
}

/* Get the starting line number of ordinary map MAP.  */

inline linenum_type
ORDINARY_MAP_STARTING_LINE_NUMBER (const line_map_ordinary *ord_map)
{
  return ord_map->to_line;
}

/* Return a positive value if map encodes locations from a system
   header, 0 otherwise. Returns 1 if ordinary map MAP encodes locations
   in a system header and 2 if it encodes locations in a C system header
   that therefore needs to be extern "C" protected in C++.  */

inline unsigned char
ORDINARY_MAP_IN_SYSTEM_HEADER_P (const line_map_ordinary *ord_map)
{
  return ord_map->sysp;
}

/* TRUE if this line map is for a module (not a source file).  */

inline bool
MAP_MODULE_P (const line_map *map)
{
  return (MAP_ORDINARY_P (map)
	  && linemap_check_ordinary (map)->reason == LC_MODULE);
}

/* Get the filename of ordinary map MAP.  */

inline const char *
ORDINARY_MAP_FILE_NAME (const line_map_ordinary *ord_map)
{
  return ord_map->to_file;
}

/* Get the cpp macro whose expansion gave birth to macro map MAP.  */

inline cpp_hashnode *
MACRO_MAP_MACRO (const line_map_macro *macro_map)
{
  return macro_map->macro;
}

/* Get the number of tokens inside the replacement-list of the macro
   that led to macro map MAP.  */

inline unsigned int
MACRO_MAP_NUM_MACRO_TOKENS (const line_map_macro *macro_map)
{
  return macro_map->n_tokens;
}

/* Get the array of pairs of locations within macro map MAP.
   See the declaration of line_map_macro for more information.  */

inline location_t *
MACRO_MAP_LOCATIONS (const line_map_macro *macro_map)
{
  return macro_map->macro_locations;
}

/* Get the location of the expansion point of the macro map MAP.  */

inline location_t
MACRO_MAP_EXPANSION_POINT_LOCATION (const line_map_macro *macro_map)
{
  return macro_map->expansion;
}

/* The abstraction of a set of location maps. There can be several
   types of location maps. This abstraction contains the attributes
   that are independent from the type of the map.

   Essentially this is just a vector of T_linemap_subclass,
   which can only ever grow in size.  */

struct GTY(()) maps_info_ordinary {
  /* This array contains the "ordinary" line maps, for all
     events other than macro expansion
     (e.g. when a new preprocessing unit starts or ends).  */
  line_map_ordinary * GTY ((length ("%h.used"))) maps;

  /* The total number of allocated maps.  */
  unsigned int allocated;

  /* The number of elements used in maps. This number is smaller
     or equal to ALLOCATED.  */
  unsigned int used;

  mutable unsigned int cache;
};

struct GTY(()) maps_info_macro {
  /* This array contains the macro line maps.
     A macro line map is created whenever a macro expansion occurs.  */
  line_map_macro * GTY ((length ("%h.used"))) maps;

  /* The total number of allocated maps.  */
  unsigned int allocated;

  /* The number of elements used in maps. This number is smaller
     or equal to ALLOCATED.  */
  unsigned int used;

  mutable unsigned int cache;
};

/* Data structure to associate a source_range together with an arbitrary
   data pointer with a source location.  */
struct GTY(()) location_adhoc_data {
  location_t locus;
  source_range src_range;
  void * GTY((skip)) data;
};

struct htab;

/* The following data structure encodes a location with some adhoc data
   and maps it to a new unsigned integer (called an adhoc location)
   that replaces the original location to represent the mapping.

   The new adhoc_loc uses the highest bit as the enabling bit, i.e. if the
   highest bit is 1, then the number is adhoc_loc. Otherwise, it serves as
   the original location. Once identified as the adhoc_loc, the lower 31
   bits of the integer is used to index the location_adhoc_data array,
   in which the locus and associated data is stored.  */

struct GTY(()) location_adhoc_data_map {
  struct htab * GTY((skip)) htab;
  location_t curr_loc;
  unsigned int allocated;
  struct location_adhoc_data GTY((length ("%h.allocated"))) *data;
};

/* A set of chronological line_map structures.  */
class GTY(()) line_maps {
public:

  ~line_maps ();
  
  maps_info_ordinary info_ordinary;

  maps_info_macro info_macro;

  /* Depth of the include stack, including the current file.  */
  unsigned int depth;

  /* If true, prints an include trace a la -H.  */
  bool trace_includes;

  /* True if we've seen a #line or # 44 "file" directive.  */
  bool seen_line_directive;

  /* Highest location_t "given out".  */
  location_t highest_location;

  /* Start of line of highest location_t "given out".  */
  location_t highest_line;

  /* The maximum column number we can quickly allocate.  Higher numbers
     may require allocating a new line_map.  */
  unsigned int max_column_hint;

  /* The allocator to use when resizing 'maps', defaults to xrealloc.  */
  line_map_realloc GTY((callback)) reallocator;

  /* The allocators' function used to know the actual size it
     allocated, for a certain allocation size requested.  */
  line_map_round_alloc_size_func GTY((callback)) round_alloc_size;

  struct location_adhoc_data_map location_adhoc_data_map;

  /* The special location value that is used as spelling location for
     built-in tokens.  */
  location_t builtin_location;

  /* The default value of range_bits in ordinary line maps.  */
  unsigned int default_range_bits;

  unsigned int num_optimized_ranges;
  unsigned int num_unoptimized_ranges;
};

/* Returns the number of allocated maps so far. MAP_KIND shall be TRUE
   if we are interested in macro maps, FALSE otherwise.  */
inline unsigned int
LINEMAPS_ALLOCATED (const line_maps *set, bool map_kind)
{
  if (map_kind)
    return set->info_macro.allocated;
  else
    return set->info_ordinary.allocated;
}

/* As above, but by reference (e.g. as an lvalue).  */

inline unsigned int &
LINEMAPS_ALLOCATED (line_maps *set, bool map_kind)
{
  if (map_kind)
    return set->info_macro.allocated;
  else
    return set->info_ordinary.allocated;
}

/* Returns the number of used maps so far. MAP_KIND shall be TRUE if
   we are interested in macro maps, FALSE otherwise.*/
inline unsigned int
LINEMAPS_USED (const line_maps *set, bool map_kind)
{
  if (map_kind)
    return set->info_macro.used;
  else
    return set->info_ordinary.used;
}

/* As above, but by reference (e.g. as an lvalue).  */

inline unsigned int &
LINEMAPS_USED (line_maps *set, bool map_kind)
{
  if (map_kind)
    return set->info_macro.used;
  else
    return set->info_ordinary.used;
}

/* Returns the index of the last map that was looked up with
   linemap_lookup. MAP_KIND shall be TRUE if we are interested in
   macro maps, FALSE otherwise.  */
inline unsigned int &
LINEMAPS_CACHE (const line_maps *set, bool map_kind)
{
  if (map_kind)
    return set->info_macro.cache;
  else
    return set->info_ordinary.cache;
}

/* Return the map at a given index.  */
inline line_map *
LINEMAPS_MAP_AT (const line_maps *set, bool map_kind, int index)
{
  if (map_kind)
    return &set->info_macro.maps[index];
  else
    return &set->info_ordinary.maps[index];
}

/* Returns the last map used in the line table SET. MAP_KIND
   shall be TRUE if we are interested in macro maps, FALSE
   otherwise.*/
inline line_map *
LINEMAPS_LAST_MAP (const line_maps *set, bool map_kind)
{
  return LINEMAPS_MAP_AT (set, map_kind,
			  LINEMAPS_USED (set, map_kind) - 1);
}

/* Returns the last map that was allocated in the line table SET.
   MAP_KIND shall be TRUE if we are interested in macro maps, FALSE
   otherwise.*/
inline line_map *
LINEMAPS_LAST_ALLOCATED_MAP (const line_maps *set, bool map_kind)
{
  return LINEMAPS_MAP_AT (set, map_kind,
			  LINEMAPS_ALLOCATED (set, map_kind) - 1);
}

/* Returns a pointer to the memory region where ordinary maps are
   allocated in the line table SET.  */
inline line_map_ordinary *
LINEMAPS_ORDINARY_MAPS (const line_maps *set)
{
  return set->info_ordinary.maps;
}

/* Returns the INDEXth ordinary map.  */
inline line_map_ordinary *
LINEMAPS_ORDINARY_MAP_AT (const line_maps *set, int index)
{
  linemap_assert (index >= 0
		  && (unsigned int)index < LINEMAPS_USED (set, false));
  return (line_map_ordinary *)LINEMAPS_MAP_AT (set, false, index);
}

/* Return the number of ordinary maps allocated in the line table
   SET.  */
inline unsigned int
LINEMAPS_ORDINARY_ALLOCATED (const line_maps *set)
{
  return LINEMAPS_ALLOCATED (set, false);
}

/* Return the number of ordinary maps used in the line table SET.  */
inline unsigned int
LINEMAPS_ORDINARY_USED (const line_maps *set)
{
  return LINEMAPS_USED (set, false);
}

/* Return the index of the last ordinary map that was looked up with
   linemap_lookup.  */
inline unsigned int &
LINEMAPS_ORDINARY_CACHE (const line_maps *set)
{
  return LINEMAPS_CACHE (set, false);
}

/* Returns a pointer to the last ordinary map used in the line table
   SET.  */
inline line_map_ordinary *
LINEMAPS_LAST_ORDINARY_MAP (const line_maps *set)
{
  return (line_map_ordinary *)LINEMAPS_LAST_MAP (set, false);
}

/* Returns a pointer to the last ordinary map allocated the line table
   SET.  */
inline line_map_ordinary *
LINEMAPS_LAST_ALLOCATED_ORDINARY_MAP (const line_maps *set)
{
  return (line_map_ordinary *)LINEMAPS_LAST_ALLOCATED_MAP (set, false);
}

/* Returns a pointer to the beginning of the region where macro maps
   are allocated.  */
inline line_map_macro *
LINEMAPS_MACRO_MAPS (const line_maps *set)
{
  return set->info_macro.maps;
}

/* Returns the INDEXth macro map.  */
inline line_map_macro *
LINEMAPS_MACRO_MAP_AT (const line_maps *set, int index)
{
  linemap_assert (index >= 0
		  && (unsigned int)index < LINEMAPS_USED (set, true));
  return (line_map_macro *)LINEMAPS_MAP_AT (set, true, index);
}

/* Returns the number of macro maps that were allocated in the line
   table SET.  */
inline unsigned int
LINEMAPS_MACRO_ALLOCATED (const line_maps *set)
{
  return LINEMAPS_ALLOCATED (set, true);
}

/* Returns the number of macro maps used in the line table SET.  */
inline unsigned int
LINEMAPS_MACRO_USED (const line_maps *set)
{
  return LINEMAPS_USED (set, true);
}

/* Return the index of the last macro map that was looked up with
   linemap_lookup.  */
inline unsigned int &
LINEMAPS_MACRO_CACHE (const line_maps *set)
{
  return LINEMAPS_CACHE (set, true);
}

/* Returns the last macro map used in the line table SET.  */
inline line_map_macro *
LINEMAPS_LAST_MACRO_MAP (const line_maps *set)
{
  return (line_map_macro *)LINEMAPS_LAST_MAP (set, true);
}

/* Returns the lowest location [of a token resulting from macro
   expansion] encoded in this line table.  */
inline location_t
LINEMAPS_MACRO_LOWEST_LOCATION (const line_maps *set)
{
  return LINEMAPS_MACRO_USED (set)
         ? MAP_START_LOCATION (LINEMAPS_LAST_MACRO_MAP (set))
         : MAX_LOCATION_T + 1;
}

/* Returns the last macro map allocated in the line table SET.  */
inline line_map_macro *
LINEMAPS_LAST_ALLOCATED_MACRO_MAP (const line_maps *set)
{
  return (line_map_macro *)LINEMAPS_LAST_ALLOCATED_MAP (set, true);
}

extern location_t get_combined_adhoc_loc (line_maps *, location_t,
					  source_range, void *);
extern void *get_data_from_adhoc_loc (const line_maps *, location_t);
extern location_t get_location_from_adhoc_loc (const line_maps *,
					       location_t);

extern source_range get_range_from_loc (line_maps *set, location_t loc);

/* Get whether location LOC is a "pure" location, or
   whether it is an ad-hoc location, or embeds range information.  */

bool
pure_location_p (line_maps *set, location_t loc);

/* Given location LOC within SET, strip away any packed range information
   or ad-hoc information.  */

extern location_t get_pure_location (line_maps *set, location_t loc);

/* Combine LOC and BLOCK, giving a combined adhoc location.  */

inline location_t
COMBINE_LOCATION_DATA (class line_maps *set,
		       location_t loc,
		       source_range src_range,
		       void *block)
{
  return get_combined_adhoc_loc (set, loc, src_range, block);
}

extern void rebuild_location_adhoc_htab (class line_maps *);

/* Initialize a line map set.  SET is the line map set to initialize
   and BUILTIN_LOCATION is the special location value to be used as
   spelling location for built-in tokens.  This BUILTIN_LOCATION has
   to be strictly less than RESERVED_LOCATION_COUNT.  */
extern void linemap_init (class line_maps *set,
			  location_t builtin_location);

/* Check for and warn about line_maps entered but not exited.  */

extern void linemap_check_files_exited (class line_maps *);

/* Return a location_t for the start (i.e. column==0) of
   (physical) line TO_LINE in the current source file (as in the
   most recent linemap_add).   MAX_COLUMN_HINT is the highest column
   number we expect to use in this line (but it does not change
   the highest_location).  */

extern location_t linemap_line_start
(class line_maps *set, linenum_type to_line,  unsigned int max_column_hint);

/* Allocate a raw block of line maps, zero initialized.  */
extern line_map *line_map_new_raw (line_maps *, bool, unsigned);

/* Add a mapping of logical source line to physical source file and
   line number. This function creates an "ordinary map", which is a
   map that records locations of tokens that are not part of macro
   replacement-lists present at a macro expansion point.

   The text pointed to by TO_FILE must have a lifetime
   at least as long as the lifetime of SET.  An empty
   TO_FILE means standard input.  If reason is LC_LEAVE, and
   TO_FILE is NULL, then TO_FILE, TO_LINE and SYSP are given their
   natural values considering the file we are returning to.

   A call to this function can relocate the previous set of
   maps, so any stored line_map pointers should not be used.  */
extern const line_map *linemap_add
  (class line_maps *, enum lc_reason, unsigned int sysp,
   const char *to_file, linenum_type to_line);

/* Create a macro map.  A macro map encodes source locations of tokens
   that are part of a macro replacement-list, at a macro expansion
   point. See the extensive comments of struct line_map and struct
   line_map_macro, in line-map.h.

   This map shall be created when the macro is expanded. The map
   encodes the source location of the expansion point of the macro as
   well as the "original" source location of each token that is part
   of the macro replacement-list. If a macro is defined but never
   expanded, it has no macro map.  SET is the set of maps the macro
   map should be part of.  MACRO_NODE is the macro which the new macro
   map should encode source locations for.  EXPANSION is the location
   of the expansion point of MACRO. For function-like macros
   invocations, it's best to make it point to the closing parenthesis
   of the macro, rather than the the location of the first character
   of the macro.  NUM_TOKENS is the number of tokens that are part of
   the replacement-list of MACRO.  */
const line_map_macro *linemap_enter_macro (line_maps *, cpp_hashnode *,
					   location_t, unsigned int);

/* Create a source location for a module.  The creator must either do
   this after the TU is tokenized, or deal with saving and restoring
   map state.  */

extern location_t linemap_module_loc
  (line_maps *, location_t from, const char *name);
extern void linemap_module_reparent
  (line_maps *, location_t loc, location_t new_parent);

/* Restore the linemap state such that the map at LWM-1 continues.
   Return start location of the new map.  */
extern unsigned linemap_module_restore
  (line_maps *, unsigned lwm);

/* Given a logical source location, returns the map which the
   corresponding (source file, line, column) triplet can be deduced
   from. Since the set is built chronologically, the logical lines are
   monotonic increasing, and so the list is sorted and we can use a
   binary search. If no line map have been allocated yet, this
   function returns NULL.  */
extern const line_map *linemap_lookup
  (const line_maps *, location_t);

unsigned linemap_lookup_macro_index (const line_maps *, location_t);

/* Returns TRUE if the line table set tracks token locations across
   macro expansion, FALSE otherwise.  */
bool linemap_tracks_macro_expansion_locs_p (class line_maps *);

/* Return the name of the macro associated to MACRO_MAP.  */
const char* linemap_map_get_macro_name (const line_map_macro *);

/* Return a positive value if LOCATION is the locus of a token that is
   located in a system header, O otherwise. It returns 1 if LOCATION
   is the locus of a token that is located in a system header, and 2
   if LOCATION is the locus of a token located in a C system header
   that therefore needs to be extern "C" protected in C++.

   Note that this function returns 1 if LOCATION belongs to a token
   that is part of a macro replacement-list defined in a system
   header, but expanded in a non-system file.  */
int linemap_location_in_system_header_p (class line_maps *,
					 location_t);

/* Return TRUE if LOCATION is a source code location of a token that is part of
   a macro expansion, FALSE otherwise.  */
bool linemap_location_from_macro_expansion_p (const line_maps *,
					      location_t);

/* TRUE if LOCATION is a source code location of a token that is part of the
   definition of a macro, FALSE otherwise.  */
bool linemap_location_from_macro_definition_p (class line_maps *,
					       location_t);

/* With the precondition that LOCATION is the locus of a token that is
   an argument of a function-like macro MACRO_MAP and appears in the
   expansion of MACRO_MAP, return the locus of that argument in the
   context of the caller of MACRO_MAP.  */

extern location_t linemap_macro_map_loc_unwind_toward_spelling
  (line_maps *set, const line_map_macro *macro_map, location_t location);

/* location_t values from 0 to RESERVED_LOCATION_COUNT-1 will
   be reserved for libcpp user as special values, no token from libcpp
   will contain any of those locations.  */
const location_t RESERVED_LOCATION_COUNT = 2;

/* Converts a map and a location_t to source line.  */
inline linenum_type
SOURCE_LINE (const line_map_ordinary *ord_map, location_t loc)
{
  return ((loc - ord_map->start_location)
	  >> ord_map->m_column_and_range_bits) + ord_map->to_line;
}

/* Convert a map and location_t to source column number.  */
inline linenum_type
SOURCE_COLUMN (const line_map_ordinary *ord_map, location_t loc)
{
  return ((loc - ord_map->start_location)
	  & ((1 << ord_map->m_column_and_range_bits) - 1)) >> ord_map->m_range_bits;
}


inline location_t
linemap_included_from (const line_map_ordinary *ord_map)
{
  return ord_map->included_from;
}

/* The linemap containing the included-from location of MAP.  */
const line_map_ordinary *linemap_included_from_linemap
  (line_maps *set, const line_map_ordinary *map);

/* True if the map is at the bottom of the include stack.  */

inline bool
MAIN_FILE_P (const line_map_ordinary *ord_map)
{
  return ord_map->included_from == 0;
}

/* Encode and return a location_t from a column number. The
   source line considered is the last source line used to call
   linemap_line_start, i.e, the last source line which a location was
   encoded from.  */
extern location_t
linemap_position_for_column (class line_maps *, unsigned int);

/* Encode and return a source location from a given line and
   column.  */
location_t
linemap_position_for_line_and_column (line_maps *set,
				      const line_map_ordinary *,
				      linenum_type, unsigned int);

/* Encode and return a location_t starting from location LOC and
   shifting it by OFFSET columns.  This function does not support
   virtual locations.  */
location_t
linemap_position_for_loc_and_offset (class line_maps *set,
				     location_t loc,
				     unsigned int offset);

/* Return the file this map is for.  */
inline const char *
LINEMAP_FILE (const line_map_ordinary *ord_map)
{
  return ord_map->to_file;
}

/* Return the line number this map started encoding location from.  */
inline linenum_type
LINEMAP_LINE (const line_map_ordinary *ord_map)
{
  return ord_map->to_line;
}

/* Return a positive value if map encodes locations from a system
   header, 0 otherwise. Returns 1 if MAP encodes locations in a
   system header and 2 if it encodes locations in a C system header
   that therefore needs to be extern "C" protected in C++.  */
inline unsigned char
LINEMAP_SYSP (const line_map_ordinary *ord_map)
{
  return ord_map->sysp;
}

const struct line_map *first_map_in_common (line_maps *set,
					    location_t loc0,
					    location_t loc1,
					    location_t *res_loc0,
					    location_t *res_loc1);

/* Return a positive value if PRE denotes the location of a token that
   comes before the token of POST, 0 if PRE denotes the location of
   the same token as the token for POST, and a negative value
   otherwise.  */
int linemap_compare_locations (class line_maps *set,
			       location_t   pre,
			       location_t   post);

/* Return TRUE if LOC_A denotes the location a token that comes
   topogically before the token denoted by location LOC_B, or if they
   are equal.  */
inline bool
linemap_location_before_p (class line_maps *set,
			   location_t loc_a,
			   location_t loc_b)
{
  return linemap_compare_locations (set, loc_a, loc_b) >= 0;
}

typedef struct
{
  /* The name of the source file involved.  */
  const char *file;

  /* The line-location in the source file.  */
  int line;

  int column;

  void *data;

  /* In a system header?. */
  bool sysp;
} expanded_location;

class range_label;

/* A hint to diagnostic_show_locus on how to print a source range within a
   rich_location.

   Typically this is SHOW_RANGE_WITH_CARET for the 0th range, and
   SHOW_RANGE_WITHOUT_CARET for subsequent ranges,
   but the Fortran frontend uses SHOW_RANGE_WITH_CARET repeatedly for
   printing things like:

       x = x + y
           1   2
       Error: Shapes for operands at (1) and (2) are not conformable

   where "1" and "2" are notionally carets.  */

enum range_display_kind
{
  /* Show the pertinent source line(s), the caret, and underline(s).  */
  SHOW_RANGE_WITH_CARET,

  /* Show the pertinent source line(s) and underline(s), but don't
     show the caret (just an underline).  */
  SHOW_RANGE_WITHOUT_CARET,

  /* Just show the source lines; don't show the range itself.
     This is for use when displaying some line-insertion fix-it hints (for
     showing the user context on the change, for when it doesn't make sense
     to highlight the first column on the next line).  */
  SHOW_LINES_WITHOUT_RANGE
};

/* A location within a rich_location: a caret&range, with
   the caret potentially flagged for display, and an optional
   label.  */

struct location_range
{
  location_t m_loc;

  enum range_display_kind m_range_display_kind;

  /* If non-NULL, the label for this range.  */
  const range_label *m_label;
};

/* A partially-embedded vec for use within rich_location for storing
   ranges and fix-it hints.

   Elements [0..NUM_EMBEDDED) are allocated within m_embed, after
   that they are within the dynamically-allocated m_extra.

   This allows for static allocation in the common case, whilst
   supporting the rarer case of an arbitrary number of elements.

   Dynamic allocation is not performed unless it's needed.  */

template <typename T, int NUM_EMBEDDED>
class semi_embedded_vec
{
 public:
  semi_embedded_vec ();
  ~semi_embedded_vec ();

  unsigned int count () const { return m_num; }
  T& operator[] (int idx);
  const T& operator[] (int idx) const;

  void push (const T&);
  void truncate (int len);

 private:
  int m_num;
  T m_embedded[NUM_EMBEDDED];
  int m_alloc;
  T *m_extra;
};

/* Constructor for semi_embedded_vec.  In particular, no dynamic allocation
   is done.  */

template <typename T, int NUM_EMBEDDED>
semi_embedded_vec<T, NUM_EMBEDDED>::semi_embedded_vec ()
: m_num (0), m_alloc (0), m_extra (NULL)
{
}

/* semi_embedded_vec's dtor.  Release any dynamically-allocated memory.  */

template <typename T, int NUM_EMBEDDED>
semi_embedded_vec<T, NUM_EMBEDDED>::~semi_embedded_vec ()
{
  XDELETEVEC (m_extra);
}

/* Look up element IDX, mutably.  */

template <typename T, int NUM_EMBEDDED>
T&
semi_embedded_vec<T, NUM_EMBEDDED>::operator[] (int idx)
{
  linemap_assert (idx < m_num);
  if (idx < NUM_EMBEDDED)
    return m_embedded[idx];
  else
    {
      linemap_assert (m_extra != NULL);
      return m_extra[idx - NUM_EMBEDDED];
    }
}

/* Look up element IDX (const).  */

template <typename T, int NUM_EMBEDDED>
const T&
semi_embedded_vec<T, NUM_EMBEDDED>::operator[] (int idx) const
{
  linemap_assert (idx < m_num);
  if (idx < NUM_EMBEDDED)
    return m_embedded[idx];
  else
    {
      linemap_assert (m_extra != NULL);
      return m_extra[idx - NUM_EMBEDDED];
    }
}

/* Append VALUE to the end of the semi_embedded_vec.  */

template <typename T, int NUM_EMBEDDED>
void
semi_embedded_vec<T, NUM_EMBEDDED>::push (const T& value)
{
  int idx = m_num++;
  if (idx < NUM_EMBEDDED)
    m_embedded[idx] = value;
  else
    {
      /* Offset "idx" to be an index within m_extra.  */
      idx -= NUM_EMBEDDED;
      if (NULL == m_extra)
	{
	  linemap_assert (m_alloc == 0);
	  m_alloc = 16;
	  m_extra = XNEWVEC (T, m_alloc);
	}
      else if (idx >= m_alloc)
	{
	  linemap_assert (m_alloc > 0);
	  m_alloc *= 2;
	  m_extra = XRESIZEVEC (T, m_extra, m_alloc);
	}
      linemap_assert (m_extra);
      linemap_assert (idx < m_alloc);
      m_extra[idx] = value;
    }
}

/* Truncate to length LEN.  No deallocation is performed.  */

template <typename T, int NUM_EMBEDDED>
void
semi_embedded_vec<T, NUM_EMBEDDED>::truncate (int len)
{
  linemap_assert (len <= m_num);
  m_num = len;
}

class fixit_hint;
class diagnostic_path;

/* A "rich" source code location, for use when printing diagnostics.
   A rich_location has one or more carets&ranges, where the carets
   are optional.  These are referred to as "ranges" from here.
   Typically the zeroth range has a caret; other ranges sometimes
   have carets.

   The "primary" location of a rich_location is the caret of range 0,
   used for determining the line/column when printing diagnostic
   text, such as:

      some-file.c:3:1: error: ...etc...

   Additional ranges may be added to help the user identify other
   pertinent clauses in a diagnostic.

   Ranges can (optionally) be given labels via class range_label.

   rich_location instances are intended to be allocated on the stack
   when generating diagnostics, and to be short-lived.

   Examples of rich locations
   --------------------------

   Example A
   *********
      int i = "foo";
              ^
   This "rich" location is simply a single range (range 0), with
   caret = start = finish at the given point.

   Example B
   *********
      a = (foo && bar)
          ~~~~~^~~~~~~
   This rich location has a single range (range 0), with the caret
   at the first "&", and the start/finish at the parentheses.
   Compare with example C below.

   Example C
   *********
      a = (foo && bar)
           ~~~ ^~ ~~~
   This rich location has three ranges:
   - Range 0 has its caret and start location at the first "&" and
     end at the second "&.
   - Range 1 has its start and finish at the "f" and "o" of "foo";
     the caret is not flagged for display, but is perhaps at the "f"
     of "foo".
   - Similarly, range 2 has its start and finish at the "b" and "r" of
     "bar"; the caret is not flagged for display, but is perhaps at the
     "b" of "bar".
   Compare with example B above.

   Example D (Fortran frontend)
   ****************************
       x = x + y
           1   2
   This rich location has range 0 at "1", and range 1 at "2".
   Both are flagged for caret display.  Both ranges have start/finish
   equal to their caret point.  The frontend overrides the diagnostic
   context's default caret character for these ranges.

   Example E (range labels)
   ************************
      printf ("arg0: %i  arg1: %s arg2: %i",
                               ^~
                               |
                               const char *
              100, 101, 102);
                   ~~~
                   |
                   int
   This rich location has two ranges:
   - range 0 is at the "%s" with start = caret = "%" and finish at
     the "s".  It has a range_label ("const char *").
   - range 1 has start/finish covering the "101" and is not flagged for
     caret printing.  The caret is at the start of "101", where its
     range_label is printed ("int").

   Fix-it hints
   ------------

   Rich locations can also contain "fix-it hints", giving suggestions
   for the user on how to edit their code to fix a problem.  These
   can be expressed as insertions, replacements, and removals of text.
   The edits by default are relative to the zeroth range within the
   rich_location, but optionally they can be expressed relative to
   other locations (using various overloaded methods of the form
   rich_location::add_fixit_*).

   For example:

   Example F: fix-it hint: insert_before
   *************************************
      ptr = arr[0];
	    ^~~~~~
	    &
   This rich location has a single range (range 0) covering "arr[0]",
   with the caret at the start.  The rich location has a single
   insertion fix-it hint, inserted before range 0, added via
     richloc.add_fixit_insert_before ("&");

   Example G: multiple fix-it hints: insert_before and insert_after
   ****************************************************************
      #define FN(ARG0, ARG1, ARG2) fn(ARG0, ARG1, ARG2)
				      ^~~~  ^~~~  ^~~~
				      (   ) (   ) (   )
   This rich location has three ranges, covering "arg0", "arg1",
   and "arg2", all with caret-printing enabled.
   The rich location has 6 insertion fix-it hints: each arg
   has a pair of insertion fix-it hints, suggesting wrapping
   them with parentheses: one a '(' inserted before,
   the other a ')' inserted after, added via
     richloc.add_fixit_insert_before (LOC, "(");
   and
     richloc.add_fixit_insert_after (LOC, ")");

   Example H: fix-it hint: removal
   *******************************
     struct s {int i};;
		      ^
		      -
   This rich location has a single range at the stray trailing
   semicolon, along with a single removal fix-it hint, covering
   the same range, added via:
     richloc.add_fixit_remove ();

   Example I: fix-it hint: replace
   *******************************
      c = s.colour;
	    ^~~~~~
	    color
   This rich location has a single range (range 0) covering "colour",
   and a single "replace" fix-it hint, covering the same range,
   added via
     richloc.add_fixit_replace ("color");

   Example J: fix-it hint: line insertion
   **************************************

     3 | #include <stddef.h>
     + |+#include <stdio.h>
     4 | int the_next_line;

   This rich location has a single range at line 4 column 1, marked
   with SHOW_LINES_WITHOUT_RANGE (to avoid printing a meaningless caret
   on the "i" of int).  It has a insertion fix-it hint of the string
   "#include <stdio.h>\n".

   Adding a fix-it hint can fail: for example, attempts to insert content
   at the transition between two line maps may fail due to there being no
   location_t value to express the new location.

   Attempts to add a fix-it hint within a macro expansion will fail.

   There is only limited support for newline characters in fix-it hints:
   only hints with newlines which insert an entire new line are permitted,
   inserting at the start of a line, and finishing with a newline
   (with no interior newline characters).  Other attempts to add
   fix-it hints containing newline characters will fail.
   Similarly, attempts to delete or replace a range *affecting* multiple
   lines will fail.

   The rich_location API handles these failures gracefully, so that
   diagnostics can attempt to add fix-it hints without each needing
   extensive checking.

   Fix-it hints within a rich_location are "atomic": if any hints can't
   be applied, none of them will be (tracked by the m_seen_impossible_fixit
   flag), and no fix-its hints will be displayed for that rich_location.
   This implies that diagnostic messages need to be worded in such a way
   that they make sense whether or not the fix-it hints are displayed,
   or that richloc.seen_impossible_fixit_p () should be checked before
   issuing the diagnostics.  */

class rich_location
{
 public:
  /* Constructors.  */

  /* Constructing from a location.  */
  rich_location (line_maps *set, location_t loc,
		 const range_label *label = NULL);

  /* Destructor.  */
  ~rich_location ();

  /* The class manages the memory pointed to by the elements of
     the M_FIXIT_HINTS vector and is not meant to be copied or
     assigned.  */
  rich_location (const rich_location &) = delete;
  void operator= (const rich_location &) = delete;

  /* Accessors.  */
  location_t get_loc () const { return get_loc (0); }
  location_t get_loc (unsigned int idx) const;

  void
  add_range (location_t loc,
	     enum range_display_kind range_display_kind
	       = SHOW_RANGE_WITHOUT_CARET,
	     const range_label *label = NULL);

  void
  set_range (unsigned int idx, location_t loc,
	     enum range_display_kind range_display_kind);

  unsigned int get_num_locations () const { return m_ranges.count (); }

  const location_range *get_range (unsigned int idx) const;
  location_range *get_range (unsigned int idx);

  expanded_location get_expanded_location (unsigned int idx);

  void
  override_column (int column);

  /* Fix-it hints.  */

  /* Methods for adding insertion fix-it hints.  */

  /* Suggest inserting NEW_CONTENT immediately before the primary
     range's start.  */
  void
  add_fixit_insert_before (const char *new_content);

  /* Suggest inserting NEW_CONTENT immediately before the start of WHERE.  */
  void
  add_fixit_insert_before (location_t where,
			   const char *new_content);

  /* Suggest inserting NEW_CONTENT immediately after the end of the primary
     range.  */
  void
  add_fixit_insert_after (const char *new_content);

  /* Suggest inserting NEW_CONTENT immediately after the end of WHERE.  */
  void
  add_fixit_insert_after (location_t where,
			  const char *new_content);

  /* Methods for adding removal fix-it hints.  */

  /* Suggest removing the content covered by range 0.  */
  void
  add_fixit_remove ();

  /* Suggest removing the content covered between the start and finish
     of WHERE.  */
  void
  add_fixit_remove (location_t where);

  /* Suggest removing the content covered by SRC_RANGE.  */
  void
  add_fixit_remove (source_range src_range);

  /* Methods for adding "replace" fix-it hints.  */

  /* Suggest replacing the content covered by range 0 with NEW_CONTENT.  */
  void
  add_fixit_replace (const char *new_content);

  /* Suggest replacing the content between the start and finish of
     WHERE with NEW_CONTENT.  */
  void
  add_fixit_replace (location_t where,
		     const char *new_content);

  /* Suggest replacing the content covered by SRC_RANGE with
     NEW_CONTENT.  */
  void
  add_fixit_replace (source_range src_range,
		     const char *new_content);

  unsigned int get_num_fixit_hints () const { return m_fixit_hints.count (); }
  fixit_hint *get_fixit_hint (int idx) const { return m_fixit_hints[idx]; }
  fixit_hint *get_last_fixit_hint () const;
  bool seen_impossible_fixit_p () const { return m_seen_impossible_fixit; }

  /* Set this if the fix-it hints are not suitable to be
     automatically applied.

     For example, if you are suggesting more than one
     mutually exclusive solution to a problem, then
     it doesn't make sense to apply all of the solutions;
     manual intervention is required.

     If set, then the fix-it hints in the rich_location will
     be printed, but will not be added to generated patches,
     or affect the modified version of the file.  */
  void fixits_cannot_be_auto_applied ()
  {
    m_fixits_cannot_be_auto_applied = true;
  }

  bool fixits_can_be_auto_applied_p () const
  {
    return !m_fixits_cannot_be_auto_applied;
  }

  /* An optional path through the code.  */
  const diagnostic_path *get_path () const { return m_path; }
  void set_path (const diagnostic_path *path) { m_path = path; }

  /* A flag for hinting that the diagnostic involves character encoding
     issues, and thus that it will be helpful to the user if we show some
     representation of how the characters in the pertinent source lines
     are encoded.
     The default is false (i.e. do not escape).
     When set to true, non-ASCII bytes in the pertinent source lines will
     be escaped in a manner controlled by the user-supplied option
     -fdiagnostics-escape-format=, so that the user can better understand
     what's going on with the encoding in their source file.  */
  bool escape_on_output_p () const { return m_escape_on_output; }
  void set_escape_on_output (bool flag) { m_escape_on_output = flag; }

private:
  bool reject_impossible_fixit (location_t where);
  void stop_supporting_fixits ();
  void maybe_add_fixit (location_t start,
			location_t next_loc,
			const char *new_content);

public:
  static const int STATICALLY_ALLOCATED_RANGES = 3;

protected:
  line_maps *m_line_table;
  semi_embedded_vec <location_range, STATICALLY_ALLOCATED_RANGES> m_ranges;

  int m_column_override;

  bool m_have_expanded_location;
  bool m_seen_impossible_fixit;
  bool m_fixits_cannot_be_auto_applied;
  bool m_escape_on_output;

  expanded_location m_expanded_location;

  static const int MAX_STATIC_FIXIT_HINTS = 2;
  semi_embedded_vec <fixit_hint *, MAX_STATIC_FIXIT_HINTS> m_fixit_hints;

  const diagnostic_path *m_path;
};

/* A struct for the result of range_label::get_text: a NUL-terminated buffer
   of localized text, and a flag to determine if the caller should "free" the
   buffer.  */

class label_text
{
public:
  label_text ()
  : m_buffer (NULL), m_caller_owned (false)
  {}

  void maybe_free ()
  {
    if (m_caller_owned)
      free (m_buffer);
  }

  /* Create a label_text instance that borrows BUFFER from a
     longer-lived owner.  */
  static label_text borrow (const char *buffer)
  {
    return label_text (const_cast <char *> (buffer), false);
  }

  /* Create a label_text instance that takes ownership of BUFFER.  */
  static label_text take (char *buffer)
  {
    return label_text (buffer, true);
  }

  /* Take ownership of the buffer, copying if necessary.  */
  char *take_or_copy ()
  {
    if (m_caller_owned)
      return m_buffer;
    else
      return xstrdup (m_buffer);
  }

  char *m_buffer;
  bool m_caller_owned;

private:
  label_text (char *buffer, bool owned)
  : m_buffer (buffer), m_caller_owned (owned)
  {}
};

/* Abstract base class for labelling a range within a rich_location
   (e.g. for labelling expressions with their type).

   Generating the text could require non-trivial work, so this work
   is delayed (via the "get_text" virtual function) until the diagnostic
   printing code "knows" it needs it, thus avoiding doing it e.g. for
   warnings that are filtered by command-line flags.  This virtual
   function also isolates libcpp and the diagnostics subsystem from
   the front-end and middle-end-specific code for generating the text
   for the labels.

   Like the rich_location instances they annotate, range_label instances
   are intended to be allocated on the stack when generating diagnostics,
   and to be short-lived.  */

class range_label
{
 public:
  virtual ~range_label () {}

  /* Get localized text for the label.
     The RANGE_IDX is provided, allowing for range_label instances to be
     shared by multiple ranges if need be (the "flyweight" design pattern).  */
  virtual label_text get_text (unsigned range_idx) const = 0;
};

/* A fix-it hint: a suggested insertion, replacement, or deletion of text.
   We handle these three types of edit with one class, by representing
   them as replacement of a half-open range:
       [start, next_loc)
   Insertions have start == next_loc: "replace" the empty string at the
   start location with the new string.
   Deletions are replacement with the empty string.

   There is only limited support for newline characters in fix-it hints
   as noted above in the comment for class rich_location.
   A fixit_hint instance can have at most one newline character; if
   present, the newline character must be the final character of
   the content (preventing e.g. fix-its that split a pre-existing line).  */

class fixit_hint
{
 public:
  fixit_hint (location_t start,
	      location_t next_loc,
	      const char *new_content);
  ~fixit_hint () { free (m_bytes); }

  bool affects_line_p (const char *file, int line) const;
  location_t get_start_loc () const { return m_start; }
  location_t get_next_loc () const { return m_next_loc; }
  bool maybe_append (location_t start,
		     location_t next_loc,
		     const char *new_content);

  const char *get_string () const { return m_bytes; }
  size_t get_length () const { return m_len; }

  bool insertion_p () const { return m_start == m_next_loc; }

  bool ends_with_newline_p () const;

 private:
  /* We don't use source_range here since, unlike most places,
     this is a half-open/half-closed range:
       [start, next_loc)
     so that we can support insertion via start == next_loc.  */
  location_t m_start;
  location_t m_next_loc;
  char *m_bytes;
  size_t m_len;
};


/* This is enum is used by the function linemap_resolve_location
   below.  The meaning of the values is explained in the comment of
   that function.  */
enum location_resolution_kind
{
  LRK_MACRO_EXPANSION_POINT,
  LRK_SPELLING_LOCATION,
  LRK_MACRO_DEFINITION_LOCATION
};

/* Resolve a virtual location into either a spelling location, an
   expansion point location or a token argument replacement point
   location.  Return the map that encodes the virtual location as well
   as the resolved location.

   If LOC is *NOT* the location of a token resulting from the
   expansion of a macro, then the parameter LRK (which stands for
   Location Resolution Kind) is ignored and the resulting location
   just equals the one given in argument.

   Now if LOC *IS* the location of a token resulting from the
   expansion of a macro, this is what happens.

   * If LRK is set to LRK_MACRO_EXPANSION_POINT
   -------------------------------

   The virtual location is resolved to the first macro expansion point
   that led to this macro expansion.

   * If LRK is set to LRK_SPELLING_LOCATION
   -------------------------------------

   The virtual location is resolved to the locus where the token has
   been spelled in the source.   This can follow through all the macro
   expansions that led to the token.

   * If LRK is set to LRK_MACRO_DEFINITION_LOCATION
   --------------------------------------

   The virtual location is resolved to the locus of the token in the
   context of the macro definition.

   If LOC is the locus of a token that is an argument of a
   function-like macro [replacing a parameter in the replacement list
   of the macro] the virtual location is resolved to the locus of the
   parameter that is replaced, in the context of the definition of the
   macro.

   If LOC is the locus of a token that is not an argument of a
   function-like macro, then the function behaves as if LRK was set to
   LRK_SPELLING_LOCATION.

   If LOC_MAP is not NULL, *LOC_MAP is set to the map encoding the
   returned location.  Note that if the returned location wasn't originally
   encoded by a map, the *MAP is set to NULL.  This can happen if LOC
   resolves to a location reserved for the client code, like
   UNKNOWN_LOCATION or BUILTINS_LOCATION in GCC.  */

location_t linemap_resolve_location (class line_maps *,
				     location_t loc,
				     enum location_resolution_kind lrk,
				     const line_map_ordinary **loc_map);

/* Suppose that LOC is the virtual location of a token coming from the
   expansion of a macro M.  This function then steps up to get the
   location L of the point where M got expanded.  If L is a spelling
   location inside a macro expansion M', then this function returns
   the point where M' was expanded.  LOC_MAP is an output parameter.
   When non-NULL, *LOC_MAP is set to the map of the returned
   location.  */
location_t linemap_unwind_toward_expansion (class line_maps *,
					    location_t loc,
					    const line_map **loc_map);

/* If LOC is the virtual location of a token coming from the expansion
   of a macro M and if its spelling location is reserved (e.g, a
   location for a built-in token), then this function unwinds (using
   linemap_unwind_toward_expansion) the location until a location that
   is not reserved and is not in a system header is reached.  In other
   words, this unwinds the reserved location until a location that is
   in real source code is reached.

   Otherwise, if the spelling location for LOC is not reserved or if
   LOC doesn't come from the expansion of a macro, the function
   returns LOC as is and *MAP is not touched.

   *MAP is set to the map of the returned location if the later is
   different from LOC.  */
location_t linemap_unwind_to_first_non_reserved_loc (class line_maps *,
						     location_t loc,
						     const line_map **map);

/* Expand source code location LOC and return a user readable source
   code location.  LOC must be a spelling (non-virtual) location.  If
   it's a location < RESERVED_LOCATION_COUNT a zeroed expanded source
   location is returned.  */
expanded_location linemap_expand_location (class line_maps *,
					   const line_map *,
					   location_t loc);

/* Statistics about maps allocation and usage as returned by
   linemap_get_statistics.  */
struct linemap_stats
{
  long num_ordinary_maps_allocated;
  long num_ordinary_maps_used;
  long ordinary_maps_allocated_size;
  long ordinary_maps_used_size;
  long num_expanded_macros;
  long num_macro_tokens;
  long num_macro_maps_used;
  long macro_maps_allocated_size;
  long macro_maps_used_size;
  long macro_maps_locations_size;
  long duplicated_macro_maps_locations_size;
  long adhoc_table_size;
  long adhoc_table_entries_used;
};

/* Return the highest location emitted for a given file for which
   there is a line map in SET.  FILE_NAME is the file name to
   consider.  If the function returns TRUE, *LOC is set to the highest
   location emitted for that file.  */
bool linemap_get_file_highest_location (class line_maps * set,
					const char *file_name,
					location_t *loc);

/* Compute and return statistics about the memory consumption of some
   parts of the line table SET.  */
void linemap_get_statistics (line_maps *, struct linemap_stats *);

/* Dump debugging information about source location LOC into the file
   stream STREAM. SET is the line map set LOC comes from.  */
void linemap_dump_location (line_maps *, location_t, FILE *);

/* Dump line map at index IX in line table SET to STREAM.  If STREAM
   is NULL, use stderr.  IS_MACRO is true if the caller wants to
   dump a macro map, false otherwise.  */
void linemap_dump (FILE *, line_maps *, unsigned, bool);

/* Dump line table SET to STREAM.  If STREAM is NULL, stderr is used.
   NUM_ORDINARY specifies how many ordinary maps to dump.  NUM_MACRO
   specifies how many macro maps to dump.  */
void line_table_dump (FILE *, line_maps *, unsigned int, unsigned int);

/* An enum for distinguishing the various parts within a location_t.  */

enum location_aspect
{
  LOCATION_ASPECT_CARET,
  LOCATION_ASPECT_START,
  LOCATION_ASPECT_FINISH
};

/* The rich_location class requires a way to expand location_t instances.
   We would directly use expand_location_to_spelling_point, which is
   implemented in gcc/input.cc, but we also need to use it for rich_location
   within genmatch.cc.
   Hence we require client code of libcpp to implement the following
   symbol.  */
extern expanded_location
linemap_client_expand_location_to_spelling_point (location_t,
						  enum location_aspect);

#endif /* !LIBCPP_LINE_MAP_H  */
