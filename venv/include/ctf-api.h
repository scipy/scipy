/* Public API to libctf.
   Copyright (C) 2019-2023 Free Software Foundation, Inc.

   This file is part of libctf.

   libctf is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 3, or (at your option) any later
   version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; see the file COPYING.  If not see
   <http://www.gnu.org/licenses/>.  */

/* This header file defines the interfaces available from the CTF debugger
   library, libctf.  This API can be used by a debugger to operate on data in
   the Compact ANSI-C Type Format (CTF).  */

#ifndef	_CTF_API_H
#define	_CTF_API_H

#include <sys/types.h>
#include <ctf.h>
#include <zlib.h>

#ifdef	__cplusplus
extern "C"
{
#endif

/* Clients can open one or more CTF containers and obtain a pointer to an
   opaque ctf_dict_t.  Types are identified by an opaque ctf_id_t token.
   They can also open or create read-only archives of CTF containers in a
   ctf_archive_t.

   These opaque definitions allow libctf to evolve without breaking clients.  */

typedef struct ctf_dict ctf_dict_t;
typedef struct ctf_archive_internal ctf_archive_t;
typedef unsigned long ctf_id_t;

/* This opaque definition allows libctf to accept BFD data structures without
   importing all the BFD noise into users' namespaces.  */

struct bfd;

/* If the debugger needs to provide the CTF library with a set of raw buffers
   for use as the CTF data, symbol table, and string table, it can do so by
   filling in ctf_sect_t structures and passing them to ctf_bufopen.

   The contents of this structure must always be in native endianness.  At read
   time, the symbol table endianness is derived from the BFD target (if BFD is
   in use): if a BFD target is not in use, please call ctf_symsect_endianness or
   ctf_arc_symsect_endianness.  */

typedef struct ctf_sect
{
  const char *cts_name;		  /* Section name (if any).  */
  const void *cts_data;		  /* Pointer to section data.  */
  size_t cts_size;		  /* Size of data in bytes.  */
  size_t cts_entsize;		  /* Size of each section entry (symtab only).  */
} ctf_sect_t;

/* A minimal symbol extracted from a linker's internal symbol table
   representation.  The symbol name can be given either via st_name or via a
   strtab offset in st_nameidx, which corresponds to one of the string offsets
   communicated via the ctf_link_add_strtab callback.   */

typedef struct ctf_link_sym
{
  /* The st_name and st_nameidx will not be accessed outside the call to
     ctf_link_shuffle_syms.  If you set st_nameidx to offset zero, make sure
     to set st_nameidx_set as well.  */

  const char *st_name;
  size_t st_nameidx;
  int st_nameidx_set;
  uint32_t st_symidx;
  uint32_t st_shndx;
  uint32_t st_type;
  uint32_t st_value;
} ctf_link_sym_t;

/* Flags applying to this specific link.  */

/* Share all types that are not in conflict.  The default.  */
#define CTF_LINK_SHARE_UNCONFLICTED 0x0

/* Share only types that are used by multiple inputs.  */
#define CTF_LINK_SHARE_DUPLICATED 0x1

/* Do a nondeduplicating link, or otherwise deduplicate "less hard", trading off
   CTF output size for link time.  */
#define CTF_LINK_NONDEDUP 0x2

/* Create empty outputs for all registered CU mappings even if no types are
   emitted into them.  */
#define CTF_LINK_EMPTY_CU_MAPPINGS 0x4

/* Omit the content of the variables section.  */
#define CTF_LINK_OMIT_VARIABLES_SECTION 0x8

/* If *unset*, filter out entries corresponding to linker-reported symbols
   from the variable section, and filter out all entries with no linker-reported
   symbols from the data object and function info sections: if set, do no
   filtering and leave all entries in place.  (This is a negative-sense flag
   because it is rare to want symbols the linker has not reported as present to
   stick around in the symtypetab sections nonetheless: relocatable links are
   the only likely case.)  */
#define CTF_LINK_NO_FILTER_REPORTED_SYMS 0x10

/* Symbolic names for CTF sections.  */

typedef enum ctf_sect_names
  {
   CTF_SECT_HEADER,
   CTF_SECT_LABEL,
   CTF_SECT_OBJT,
   CTF_SECT_OBJTIDX = CTF_SECT_OBJT,
   CTF_SECT_FUNC,
   CTF_SECT_FUNCIDX = CTF_SECT_FUNC,
   CTF_SECT_VAR,
   CTF_SECT_TYPE,
   CTF_SECT_STR
  } ctf_sect_names_t;

/* Encoding information for integers, floating-point values, and certain other
   intrinsics can be obtained by calling ctf_type_encoding, below.  The flags
   field will contain values appropriate for the type defined in <ctf.h>.  */

typedef struct ctf_encoding
{
  uint32_t cte_format;		 /* Data format (CTF_INT_* or CTF_FP_* flags).  */
  uint32_t cte_offset;		 /* Offset of value in bits.  */
  uint32_t cte_bits;		 /* Size of storage in bits.  */
} ctf_encoding_t;

typedef struct ctf_membinfo
{
  ctf_id_t ctm_type;		/* Type of struct or union member.  */
  unsigned long ctm_offset;	/* Offset of member in bits.  */
} ctf_membinfo_t;

typedef struct ctf_arinfo
{
  ctf_id_t ctr_contents;	/* Type of array contents.  */
  ctf_id_t ctr_index;		/* Type of array index.  */
  uint32_t ctr_nelems;		/* Number of elements.  */
} ctf_arinfo_t;

typedef struct ctf_funcinfo
{
  ctf_id_t ctc_return;		/* Function return type.  */
  uint32_t ctc_argc;		/* Number of typed arguments to function.  */
  uint32_t ctc_flags;		/* Function attributes (see below).  */
} ctf_funcinfo_t;

typedef struct ctf_lblinfo
{
  ctf_id_t ctb_type;		/* Last type associated with the label.  */
} ctf_lblinfo_t;

typedef struct ctf_snapshot_id
{
  unsigned long dtd_id;		/* Highest DTD ID at time of snapshot.  */
  unsigned long snapshot_id;	/* Snapshot id at time of snapshot.  */
} ctf_snapshot_id_t;

#define	CTF_FUNC_VARARG	0x1	/* Function arguments end with varargs.  */

/* Functions that return a ctf_id_t use the following value to indicate failure.
   ctf_errno can be used to obtain an error code.  Functions that return
   a straight integral -1 also use ctf_errno.  */
#define	CTF_ERR	((ctf_id_t) -1L)

/* This macro holds information about all the available ctf errors.
   It is used to form both an enum holding all the error constants,
   and also the error strings themselves.  To use, define _CTF_FIRST
   and _CTF_ITEM to expand as you like, then mention the macro name.
   See the enum after this for an example.  */
#define _CTF_ERRORS \
  _CTF_FIRST (ECTF_FMT, "File is not in CTF or ELF format.")	\
  _CTF_ITEM (ECTF_BFDERR, "BFD error.")				\
  _CTF_ITEM (ECTF_CTFVERS, "CTF dict version is too new for libctf.") \
  _CTF_ITEM (ECTF_BFD_AMBIGUOUS, "Ambiguous BFD target.")	\
  _CTF_ITEM (ECTF_SYMTAB, "Symbol table uses invalid entry size.") \
  _CTF_ITEM (ECTF_SYMBAD, "Symbol table data buffer is not valid.") \
  _CTF_ITEM (ECTF_STRBAD, "String table data buffer is not valid.") \
  _CTF_ITEM (ECTF_CORRUPT, "File data structure corruption detected.") \
  _CTF_ITEM (ECTF_NOCTFDATA, "File does not contain CTF data.") \
  _CTF_ITEM (ECTF_NOCTFBUF, "Buffer does not contain CTF data.") \
  _CTF_ITEM (ECTF_NOSYMTAB, "Symbol table information is not available.") \
  _CTF_ITEM (ECTF_NOPARENT, "The parent CTF dictionary is unavailable.") \
  _CTF_ITEM (ECTF_DMODEL, "Data model mismatch.") \
  _CTF_ITEM (ECTF_LINKADDEDLATE, "File added to link too late.") \
  _CTF_ITEM (ECTF_ZALLOC, "Failed to allocate (de)compression buffer.") \
  _CTF_ITEM (ECTF_DECOMPRESS, "Failed to decompress CTF data.") \
  _CTF_ITEM (ECTF_STRTAB, "External string table is not available.") \
  _CTF_ITEM (ECTF_BADNAME, "String name offset is corrupt.") \
  _CTF_ITEM (ECTF_BADID, "Invalid type identifier.") \
  _CTF_ITEM (ECTF_NOTSOU, "Type is not a struct or union.") \
  _CTF_ITEM (ECTF_NOTENUM, "Type is not an enum.") \
  _CTF_ITEM (ECTF_NOTSUE, "Type is not a struct, union, or enum.") \
  _CTF_ITEM (ECTF_NOTINTFP, "Type is not an integer, float, or enum.") \
  _CTF_ITEM (ECTF_NOTARRAY, "Type is not an array.") \
  _CTF_ITEM (ECTF_NOTREF, "Type does not reference another type.") \
  _CTF_ITEM (ECTF_NAMELEN, "Buffer is too small to hold type name.") \
  _CTF_ITEM (ECTF_NOTYPE, "No type found corresponding to name.") \
  _CTF_ITEM (ECTF_SYNTAX, "Syntax error in type name.") \
  _CTF_ITEM (ECTF_NOTFUNC, "Symbol table entry or type is not a function.") \
  _CTF_ITEM (ECTF_NOFUNCDAT, "No function information available for function.") \
  _CTF_ITEM (ECTF_NOTDATA, "Symbol table entry does not refer to a data object.") \
  _CTF_ITEM (ECTF_NOTYPEDAT, "No type information available for symbol.") \
  _CTF_ITEM (ECTF_NOLABEL, "No label found corresponding to name.") \
  _CTF_ITEM (ECTF_NOLABELDATA, "File does not contain any labels.") \
  _CTF_ITEM (ECTF_NOTSUP, "Feature not supported.") \
  _CTF_ITEM (ECTF_NOENUMNAM, "Enum element name not found.") \
  _CTF_ITEM (ECTF_NOMEMBNAM, "Member name not found.") \
  _CTF_ITEM (ECTF_RDONLY, "CTF container is read-only.") \
  _CTF_ITEM (ECTF_DTFULL, "CTF type is full (no more members allowed).") \
  _CTF_ITEM (ECTF_FULL, "CTF container is full.") \
  _CTF_ITEM (ECTF_DUPLICATE, "Duplicate member or variable name.") \
  _CTF_ITEM (ECTF_CONFLICT, "Conflicting type is already defined.") \
  _CTF_ITEM (ECTF_OVERROLLBACK, "Attempt to roll back past a ctf_update.") \
  _CTF_ITEM (ECTF_COMPRESS, "Failed to compress CTF data.") \
  _CTF_ITEM (ECTF_ARCREATE, "Error creating CTF archive.") \
  _CTF_ITEM (ECTF_ARNNAME, "Name not found in CTF archive.") \
  _CTF_ITEM (ECTF_SLICEOVERFLOW, "Overflow of type bitness or offset in slice.") \
  _CTF_ITEM (ECTF_DUMPSECTUNKNOWN, "Unknown section number in dump.") \
  _CTF_ITEM (ECTF_DUMPSECTCHANGED, "Section changed in middle of dump.") \
  _CTF_ITEM (ECTF_NOTYET, "Feature not yet implemented.") \
  _CTF_ITEM (ECTF_INTERNAL, "Internal error: assertion failure.") \
  _CTF_ITEM (ECTF_NONREPRESENTABLE, "Type not representable in CTF.") \
  _CTF_ITEM (ECTF_NEXT_END, "End of iteration.") \
  _CTF_ITEM (ECTF_NEXT_WRONGFUN, "Wrong iteration function called.") \
  _CTF_ITEM (ECTF_NEXT_WRONGFP, "Iteration entity changed in mid-iterate.") \
  _CTF_ITEM (ECTF_FLAGS, "CTF header contains flags unknown to libctf.") \
  _CTF_ITEM (ECTF_NEEDSBFD, "This feature needs a libctf with BFD support.") \
  _CTF_ITEM (ECTF_INCOMPLETE, "Type is not a complete type.") \
  _CTF_ITEM (ECTF_NONAME, "Type name must not be empty.")

#define	ECTF_BASE	1000	/* Base value for libctf errnos.  */

enum
  {
#define _CTF_FIRST(NAME, STR) NAME = ECTF_BASE
#define _CTF_ITEM(NAME, STR) , NAME
_CTF_ERRORS
#undef _CTF_ITEM
#undef _CTF_FIRST
  };

#define ECTF_NERR (ECTF_NONAME - ECTF_BASE + 1) /* Count of CTF errors.  */

/* The CTF data model is inferred to be the caller's data model or the data
   model of the given object, unless ctf_setmodel is explicitly called.  */
#define	CTF_MODEL_ILP32 1	/* Object data model is ILP32.  */
#define	CTF_MODEL_LP64  2	/* Object data model is LP64.  */
#ifdef _LP64
# define CTF_MODEL_NATIVE CTF_MODEL_LP64
#else
# define CTF_MODEL_NATIVE CTF_MODEL_ILP32
#endif

/* Dynamic CTF containers can be created using ctf_create.  The ctf_add_*
   routines can be used to add new definitions to the dynamic container.
   New types are labeled as root or non-root to determine whether they are
   visible at the top-level program scope when subsequently doing a lookup.  */

#define	CTF_ADD_NONROOT	0	/* Type only visible in nested scope.  */
#define	CTF_ADD_ROOT	1	/* Type visible at top-level scope.  */

/* Flags for ctf_member_next.  */

#define CTF_MN_RECURSE 0x1	/* Recurse into unnamed members.  */

/* These typedefs are used to define the signature for callback functions that
   can be used with the iteration and visit functions below.  There is also a
   family of iteration functions that do not require callbacks.  */

typedef int ctf_visit_f (const char *name, ctf_id_t type, unsigned long offset,
			 int depth, void *arg);
typedef int ctf_member_f (const char *name, ctf_id_t membtype,
			  unsigned long offset, void *arg);
typedef int ctf_enum_f (const char *name, int val, void *arg);
typedef int ctf_variable_f (const char *name, ctf_id_t type, void *arg);
typedef int ctf_type_f (ctf_id_t type, void *arg);
typedef int ctf_type_all_f (ctf_id_t type, int flag, void *arg);
typedef int ctf_label_f (const char *name, const ctf_lblinfo_t *info,
			 void *arg);
typedef int ctf_archive_member_f (ctf_dict_t *fp, const char *name, void *arg);
typedef int ctf_archive_raw_member_f (const char *name, const void *content,
				      size_t len, void *arg);
typedef char *ctf_dump_decorate_f (ctf_sect_names_t sect,
				   char *line, void *arg);

typedef struct ctf_dump_state ctf_dump_state_t;

/* Iteration state for the _next functions, and allocators/copiers/freers for
   it.  (None of these are needed for the simple case of iterating to the end:
   the _next function allocate and free the iterators for you.)  */

typedef struct ctf_next ctf_next_t;
extern ctf_next_t *ctf_next_create (void);
extern void ctf_next_destroy (ctf_next_t *);
extern ctf_next_t *ctf_next_copy (ctf_next_t *);

/* Opening.  These mostly return an abstraction over both CTF files and CTF
   archives: so they can be used to open both.  CTF files will appear to be an
   archive with one member named '.ctf'.  The low-level functions
   ctf_simple_open and ctf_bufopen return ctf_dict_t's directly, and cannot
   be used on CTF archives.

   Some of these functions take raw symtab and strtab section content in the
   form of ctf_sect_t structures.  For CTF in ELF files, these should be
   extracted from .dynsym and its associated string table (usually .dynsym)
   whenever the CTF_F_DYNSTR flag is set in the CTF preamble (which it almost
   always will be for linked objects, but not for .o files).  */

extern ctf_archive_t *ctf_bfdopen (struct bfd *, int *);
extern ctf_archive_t *ctf_bfdopen_ctfsect (struct bfd *, const ctf_sect_t *,
					   int *);
extern ctf_archive_t *ctf_fdopen (int fd, const char *filename,
				  const char *target, int *errp);
extern ctf_archive_t *ctf_open (const char *filename,
				const char *target, int *errp);
extern void ctf_close (ctf_archive_t *);
extern ctf_sect_t ctf_getdatasect (const ctf_dict_t *);
extern ctf_sect_t ctf_getsymsect (const ctf_dict_t *);
extern ctf_sect_t ctf_getstrsect (const ctf_dict_t *);
extern void ctf_symsect_endianness (ctf_dict_t *, int little_endian);
extern ctf_archive_t *ctf_get_arc (const ctf_dict_t *);
extern ctf_archive_t *ctf_arc_open (const char *, int *);
extern ctf_archive_t *ctf_arc_bufopen (const ctf_sect_t *,
				       const ctf_sect_t *,
				       const ctf_sect_t *,
				       int *);
extern void ctf_arc_symsect_endianness (ctf_archive_t *, int little_endian);
extern void ctf_arc_close (ctf_archive_t *);
extern ctf_dict_t *ctf_arc_lookup_symbol (ctf_archive_t *,
					  unsigned long symidx,
					  ctf_id_t *, int *errp);
extern ctf_dict_t *ctf_arc_lookup_symbol_name (ctf_archive_t *,
					       const char *name,
					       ctf_id_t *, int *errp);
extern void ctf_arc_flush_caches (ctf_archive_t *);
extern ctf_dict_t *ctf_dict_open (const ctf_archive_t *,
				  const char *, int *);
extern ctf_dict_t *ctf_dict_open_sections (const ctf_archive_t *,
					   const ctf_sect_t *,
					   const ctf_sect_t *,
					   const char *, int *);
extern size_t ctf_archive_count (const ctf_archive_t *);

/* The next functions return or close real CTF files, or write out CTF archives,
   not opaque containers around either.  */

extern ctf_dict_t *ctf_simple_open (const char *, size_t, const char *, size_t,
				   size_t, const char *, size_t, int *);
extern ctf_dict_t *ctf_bufopen (const ctf_sect_t *, const ctf_sect_t *,
				const ctf_sect_t *, int *);
extern void ctf_ref (ctf_dict_t *);
extern void ctf_dict_close (ctf_dict_t *);

extern int ctf_arc_write (const char *, ctf_dict_t **, size_t,
			  const char **, size_t);
extern int ctf_arc_write_fd (int, ctf_dict_t **, size_t, const char **,
			     size_t);

extern const char *ctf_cuname (ctf_dict_t *);
extern int ctf_cuname_set (ctf_dict_t *, const char *);
extern ctf_dict_t *ctf_parent_dict (ctf_dict_t *);
extern const char *ctf_parent_name (ctf_dict_t *);
extern int ctf_parent_name_set (ctf_dict_t *, const char *);
extern int ctf_type_isparent (ctf_dict_t *, ctf_id_t);
extern int ctf_type_ischild (ctf_dict_t *, ctf_id_t);

extern int ctf_import (ctf_dict_t *, ctf_dict_t *);
extern int ctf_setmodel (ctf_dict_t *, int);
extern int ctf_getmodel (ctf_dict_t *);

extern void ctf_setspecific (ctf_dict_t *, void *);
extern void *ctf_getspecific (ctf_dict_t *);

extern int ctf_errno (ctf_dict_t *);
extern const char *ctf_errmsg (int);
extern int ctf_version (int);

extern int ctf_func_info (ctf_dict_t *, unsigned long, ctf_funcinfo_t *);
extern int ctf_func_args (ctf_dict_t *, unsigned long, uint32_t, ctf_id_t *);
extern int ctf_func_type_info (ctf_dict_t *, ctf_id_t, ctf_funcinfo_t *);
extern int ctf_func_type_args (ctf_dict_t *, ctf_id_t, uint32_t, ctf_id_t *);

extern ctf_id_t ctf_lookup_by_name (ctf_dict_t *, const char *);
extern ctf_id_t ctf_lookup_by_symbol (ctf_dict_t *, unsigned long);
extern ctf_id_t ctf_lookup_by_symbol_name (ctf_dict_t *, const char *);
extern ctf_id_t ctf_symbol_next (ctf_dict_t *, ctf_next_t **,
				 const char **name, int functions);
extern ctf_id_t ctf_lookup_variable (ctf_dict_t *, const char *);

extern ctf_id_t ctf_type_resolve (ctf_dict_t *, ctf_id_t);
extern char *ctf_type_aname (ctf_dict_t *, ctf_id_t);
extern char *ctf_type_aname_raw (ctf_dict_t *, ctf_id_t);
extern ssize_t ctf_type_lname (ctf_dict_t *, ctf_id_t, char *, size_t);
extern char *ctf_type_name (ctf_dict_t *, ctf_id_t, char *, size_t);
extern const char *ctf_type_name_raw (ctf_dict_t *, ctf_id_t);
extern ssize_t ctf_type_size (ctf_dict_t *, ctf_id_t);
extern ssize_t ctf_type_align (ctf_dict_t *, ctf_id_t);
extern int ctf_type_kind (ctf_dict_t *, ctf_id_t);
extern int ctf_type_kind_forwarded (ctf_dict_t *, ctf_id_t);
extern ctf_id_t ctf_type_reference (ctf_dict_t *, ctf_id_t);
extern ctf_id_t ctf_type_pointer (ctf_dict_t *, ctf_id_t);
extern int ctf_type_encoding (ctf_dict_t *, ctf_id_t, ctf_encoding_t *);
extern int ctf_type_visit (ctf_dict_t *, ctf_id_t, ctf_visit_f *, void *);
extern int ctf_type_cmp (ctf_dict_t *, ctf_id_t, ctf_dict_t *, ctf_id_t);
extern int ctf_type_compat (ctf_dict_t *, ctf_id_t, ctf_dict_t *, ctf_id_t);

extern int ctf_member_info (ctf_dict_t *, ctf_id_t, const char *,
			    ctf_membinfo_t *);
extern int ctf_array_info (ctf_dict_t *, ctf_id_t, ctf_arinfo_t *);

extern const char *ctf_enum_name (ctf_dict_t *, ctf_id_t, int);
extern int ctf_enum_value (ctf_dict_t *, ctf_id_t, const char *, int *);

extern void ctf_label_set (ctf_dict_t *, const char *);
extern const char *ctf_label_get (ctf_dict_t *);

extern const char *ctf_label_topmost (ctf_dict_t *);
extern int ctf_label_info (ctf_dict_t *, const char *, ctf_lblinfo_t *);

extern int ctf_member_count (ctf_dict_t *, ctf_id_t);
extern int ctf_member_iter (ctf_dict_t *, ctf_id_t, ctf_member_f *, void *);
extern ssize_t ctf_member_next (ctf_dict_t *, ctf_id_t, ctf_next_t **,
				const char **name, ctf_id_t *membtype,
				int flags);
extern int ctf_enum_iter (ctf_dict_t *, ctf_id_t, ctf_enum_f *, void *);
extern const char *ctf_enum_next (ctf_dict_t *, ctf_id_t, ctf_next_t **,
				  int *);
extern int ctf_type_iter (ctf_dict_t *, ctf_type_f *, void *);
extern int ctf_type_iter_all (ctf_dict_t *, ctf_type_all_f *, void *);
extern ctf_id_t ctf_type_next (ctf_dict_t *, ctf_next_t **,
			       int *flag, int want_hidden);
extern int ctf_label_iter (ctf_dict_t *, ctf_label_f *, void *);
extern int ctf_label_next (ctf_dict_t *, ctf_next_t **, const char **); /* TBD */
extern int ctf_variable_iter (ctf_dict_t *, ctf_variable_f *, void *);
extern ctf_id_t ctf_variable_next (ctf_dict_t *, ctf_next_t **,
				   const char **);
extern int ctf_archive_iter (const ctf_archive_t *, ctf_archive_member_f *,
			     void *);
extern ctf_dict_t *ctf_archive_next (const ctf_archive_t *, ctf_next_t **,
				     const char **, int skip_parent, int *errp);

/* This function alone does not currently operate on CTF files masquerading
   as archives, and returns -EINVAL: the raw data is no longer available.  It is
   expected to be used only by archiving tools, in any case, which have no need
   to deal with non-archives at all.  */
extern int ctf_archive_raw_iter (const ctf_archive_t *,
				 ctf_archive_raw_member_f *, void *);
extern char *ctf_dump (ctf_dict_t *, ctf_dump_state_t **state,
		       ctf_sect_names_t sect, ctf_dump_decorate_f *,
		       void *arg);

/* Error-warning reporting: an 'iterator' that returns errors and warnings from
   the error/warning list, in order of emission.  Errors and warnings are popped
   after return: the caller must free the returned error-text pointer.  */
extern char *ctf_errwarning_next (ctf_dict_t *, ctf_next_t **,
				  int *is_warning, int *errp);

extern ctf_id_t ctf_add_array (ctf_dict_t *, uint32_t,
			       const ctf_arinfo_t *);
extern ctf_id_t ctf_add_const (ctf_dict_t *, uint32_t, ctf_id_t);
extern ctf_id_t ctf_add_enum_encoded (ctf_dict_t *, uint32_t, const char *,
				      const ctf_encoding_t *);
extern ctf_id_t ctf_add_enum (ctf_dict_t *, uint32_t, const char *);
extern ctf_id_t ctf_add_float (ctf_dict_t *, uint32_t,
			       const char *, const ctf_encoding_t *);
extern ctf_id_t ctf_add_forward (ctf_dict_t *, uint32_t, const char *,
				 uint32_t);
extern ctf_id_t ctf_add_function (ctf_dict_t *, uint32_t,
				  const ctf_funcinfo_t *, const ctf_id_t *);
extern ctf_id_t ctf_add_integer (ctf_dict_t *, uint32_t, const char *,
				 const ctf_encoding_t *);
extern ctf_id_t ctf_add_slice (ctf_dict_t *, uint32_t, ctf_id_t, const ctf_encoding_t *);
extern ctf_id_t ctf_add_pointer (ctf_dict_t *, uint32_t, ctf_id_t);
extern ctf_id_t ctf_add_type (ctf_dict_t *, ctf_dict_t *, ctf_id_t);
extern ctf_id_t ctf_add_typedef (ctf_dict_t *, uint32_t, const char *,
				 ctf_id_t);
extern ctf_id_t ctf_add_restrict (ctf_dict_t *, uint32_t, ctf_id_t);
extern ctf_id_t ctf_add_struct (ctf_dict_t *, uint32_t, const char *);
extern ctf_id_t ctf_add_union (ctf_dict_t *, uint32_t, const char *);
extern ctf_id_t ctf_add_struct_sized (ctf_dict_t *, uint32_t, const char *,
				      size_t);
extern ctf_id_t ctf_add_union_sized (ctf_dict_t *, uint32_t, const char *,
				     size_t);
extern ctf_id_t ctf_add_unknown (ctf_dict_t *, uint32_t, const char *);
extern ctf_id_t ctf_add_volatile (ctf_dict_t *, uint32_t, ctf_id_t);

extern int ctf_add_enumerator (ctf_dict_t *, ctf_id_t, const char *, int);
extern int ctf_add_member (ctf_dict_t *, ctf_id_t, const char *, ctf_id_t);
extern int ctf_add_member_offset (ctf_dict_t *, ctf_id_t, const char *,
				  ctf_id_t, unsigned long);
extern int ctf_add_member_encoded (ctf_dict_t *, ctf_id_t, const char *,
				   ctf_id_t, unsigned long,
				   const ctf_encoding_t);

extern int ctf_add_variable (ctf_dict_t *, const char *, ctf_id_t);

extern int ctf_add_objt_sym (ctf_dict_t *, const char *, ctf_id_t);
extern int ctf_add_func_sym (ctf_dict_t *, const char *, ctf_id_t);

extern int ctf_set_array (ctf_dict_t *, ctf_id_t, const ctf_arinfo_t *);

extern ctf_dict_t *ctf_create (int *);
extern int ctf_update (ctf_dict_t *);
extern ctf_snapshot_id_t ctf_snapshot (ctf_dict_t *);
extern int ctf_rollback (ctf_dict_t *, ctf_snapshot_id_t);
extern int ctf_discard (ctf_dict_t *);
extern int ctf_write (ctf_dict_t *, int);
extern int ctf_gzwrite (ctf_dict_t *fp, gzFile fd);
extern int ctf_compress_write (ctf_dict_t * fp, int fd);
extern unsigned char *ctf_write_mem (ctf_dict_t *, size_t *, size_t threshold);

extern int ctf_link_add_ctf (ctf_dict_t *, ctf_archive_t *, const char *);
/* The variable filter should return nonzero if a variable should not
   appear in the output.  */
typedef int ctf_link_variable_filter_f (ctf_dict_t *, const char *, ctf_id_t,
					void *);
extern int ctf_link_set_variable_filter (ctf_dict_t *,
					 ctf_link_variable_filter_f *, void *);
extern int ctf_link (ctf_dict_t *, int flags);
typedef const char *ctf_link_strtab_string_f (uint32_t *offset, void *arg);
extern int ctf_link_add_strtab (ctf_dict_t *, ctf_link_strtab_string_f *,
				void *);
extern int ctf_link_add_linker_symbol (ctf_dict_t *, ctf_link_sym_t *);
extern int ctf_link_shuffle_syms (ctf_dict_t *);
extern unsigned char *ctf_link_write (ctf_dict_t *, size_t *size,
				      size_t threshold);

/* Specialist linker functions.  These functions are not used by ld, but can be
   used by other programs making use of the linker machinery for other purposes
   to customize its output.  */
extern int ctf_link_add_cu_mapping (ctf_dict_t *, const char *from,
				    const char *to);
typedef char *ctf_link_memb_name_changer_f (ctf_dict_t *,
					    const char *, void *);
extern void ctf_link_set_memb_name_changer
  (ctf_dict_t *, ctf_link_memb_name_changer_f *, void *);

extern void ctf_setdebug (int debug);
extern int ctf_getdebug (void);

/* Deprecated aliases for existing functions and types.  */

struct ctf_file;
typedef struct ctf_dict ctf_file_t;
extern void ctf_file_close (ctf_file_t *);
extern ctf_dict_t *ctf_parent_file (ctf_dict_t *);
extern ctf_dict_t *ctf_arc_open_by_name (const ctf_archive_t *,
					 const char *, int *);
extern ctf_dict_t *ctf_arc_open_by_name_sections (const ctf_archive_t *,
						  const ctf_sect_t *,
						  const ctf_sect_t *,
						  const char *, int *);

#ifdef	__cplusplus
}
#endif

#endif				/* _CTF_API_H */
