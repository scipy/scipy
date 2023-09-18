/* ctfc.h - Declarations and definitions related to the CTF container.
   Copyright (C) 2019-2022 Free Software Foundation, Inc.

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

/* This file defines the data structures and functions used by the compiler
   to generate the CTF debug info.  The definitions below are compiler internal
   representations and closely reflect the CTF format requirements in <ctf.h>.

   The contents of the CTF container are used eventually for emission of both
   CTF (ctfout.cc) and BTF debug info (btfout.cc), as the two type debug formats
   are close cousins.  */

#ifndef GCC_CTFC_H
#define GCC_CTFC_H 1

#include "config.h"
#include "system.h"
#include "tree.h"
#include "fold-const.h"
#include "dwarf2ctf.h"
#include "ctf.h"
#include "btf.h"

/* Invalid CTF type ID definition.  */

#define CTF_NULL_TYPEID 0

/* Value to start generating the CTF type ID from.  */

#define CTF_INIT_TYPEID 1

/* CTF type ID.  */

typedef uint64_t ctf_id_t;

/* CTF string table element (list node).  */

typedef struct GTY ((chain_next ("%h.cts_next"))) ctf_string
{
  const char * cts_str;		  /* CTF string.  */
  struct ctf_string * cts_next;   /* A list node.  */
} ctf_string_t;

/* Internal representation of CTF string table.  */

typedef struct GTY (()) ctf_strtable
{
  ctf_string_t * ctstab_head;	    /* Head str ptr.  */
  ctf_string_t * ctstab_tail;	    /* Tail.  new str appended to tail.  */
  int ctstab_num;		    /* Number of strings in the table.  */
  size_t ctstab_len;		    /* Size of string table in bytes.  */
  const char * ctstab_estr;	    /* Empty string "".  */
} ctf_strtable_t;

/* Encoding information for integers, floating-point values etc.  The flags
   field will contain values appropriate for the type defined in <ctf.h>.  */

typedef struct GTY (()) ctf_encoding
{
  unsigned int cte_format;  /* Data format (CTF_INT_* or CTF_FP_* flags).  */
  unsigned int cte_offset;  /* Offset of value in bits.  */
  unsigned int cte_bits;    /* Size of storage in bits.  */
} ctf_encoding_t;

/* Array information for CTF generation.  */

typedef struct GTY (()) ctf_arinfo
{
  ctf_id_t ctr_contents;	/* Type of array contents.  */
  ctf_id_t ctr_index;		/* Type of array index.  */
  unsigned int ctr_nelems;	/* Number of elements.  */
} ctf_arinfo_t;

/* Function information for CTF generation.  */

typedef struct GTY (()) ctf_funcinfo
{
  ctf_id_t ctc_return;		/* Function return type.  */
  unsigned int ctc_argc;	/* Number of typed arguments to function.  */
  unsigned int ctc_flags;	/* Function attributes (see below).  */
} ctf_funcinfo_t;

typedef struct GTY (()) ctf_sliceinfo
{
  unsigned int cts_type;	/* Reference CTF type.  */
  unsigned short cts_offset;	/* Offset in bits of the first bit.  */
  unsigned short cts_bits;	/* Size in bits.  */
} ctf_sliceinfo_t;

/* CTF type representation internal to the compiler.  It closely reflects the
   ctf_type_t type node in <ctf.h> except the GTY (()) tags.  */

typedef struct GTY (()) ctf_itype
{
  uint32_t ctti_name;		/* Reference to name in string table.  */
  uint32_t ctti_info;		/* Encoded kind, variant length (see below).  */
  union GTY ((desc ("0")))
  {
    uint32_t GTY ((tag ("0"))) _size;/* Size of entire type in bytes.  */
    uint32_t GTY ((tag ("1"))) _type;/* Reference to another type.  */
  } _u;
  uint32_t ctti_lsizehi;	/* High 32 bits of type size in bytes.  */
  uint32_t ctti_lsizelo;	/* Low 32 bits of type size in bytes.  */
} ctf_itype_t;

#define ctti_size _u._size
#define ctti_type _u._type

/* Function arguments end with varargs.  */

#define CTF_FUNC_VARARG 0x1

/* Struct/union/enum member definition for CTF generation.  */

typedef struct GTY ((chain_next ("%h.dmd_next"))) ctf_dmdef
{
  const char * dmd_name;	/* Name of this member.  */
  ctf_id_t dmd_type;		/* Type of this member (for sou).  */
  uint32_t dmd_name_offset;	/* Offset of the name in str table.  */
  uint64_t dmd_offset;		/* Offset of this member in bits (for sou).  */
  int dmd_value;		/* Value of this member (for enum).  */
  struct ctf_dmdef * dmd_next;	/* A list node.  */
} ctf_dmdef_t;

#define ctf_dmd_list_next(elem) ((ctf_dmdef_t *)((elem)->dmd_next))

/* Function Argument.  */

typedef struct GTY (()) ctf_func_arg
{
  ctf_id_t farg_type;		  /* Type identifier of the argument.  */
  const char * farg_name;	  /* Name of the argument.  */
  uint32_t farg_name_offset;	  /* Offset of the name in str table.  */
  struct ctf_func_arg * farg_next;/* A list node.  */
} ctf_func_arg_t;

#define ctf_farg_list_next(elem) ((ctf_func_arg_t *)((elem)->farg_next))

/* Type definition for CTF generation.  */

struct GTY ((for_user)) ctf_dtdef
{
  dw_die_ref dtd_key;	      /* Type key for hashing.  */
  const char * dtd_name;      /* Name associated with definition (if any).  */
  ctf_id_t dtd_type;	      /* Type identifier for this definition.  */
  ctf_itype_t dtd_data;	      /* Type node.  */
  bool from_global_func; /* Whether this type was added from a global
			    function.  */
  union GTY ((desc ("ctf_dtu_d_union_selector (&%1)")))
  {
    /* struct, union, or enum.  */
    ctf_dmdef_t * GTY ((tag ("CTF_DTU_D_MEMBERS"))) dtu_members;
    /* array.  */
    ctf_arinfo_t GTY ((tag ("CTF_DTU_D_ARRAY"))) dtu_arr;
    /* integer or float.  */
    ctf_encoding_t GTY ((tag ("CTF_DTU_D_ENCODING"))) dtu_enc;
    /* function.  */
    ctf_func_arg_t * GTY ((tag ("CTF_DTU_D_ARGUMENTS"))) dtu_argv;
    /* slice.  */
    ctf_sliceinfo_t GTY ((tag ("CTF_DTU_D_SLICE"))) dtu_slice;
  } dtd_u;
};

typedef struct ctf_dtdef ctf_dtdef_t;

/* Variable definition for CTF generation.  */

struct GTY ((for_user)) ctf_dvdef
{
  dw_die_ref dvd_key;		/* DWARF DIE corresponding to the variable.  */
  const char * dvd_name;	/* Name associated with variable.  */
  uint32_t dvd_name_offset;	/* Offset of the name in str table.  */
  unsigned int dvd_visibility;	/* External visibility.  0=static,1=global.  */
  ctf_id_t dvd_type;		/* Type of variable.  */
};

typedef struct ctf_dvdef ctf_dvdef_t;

typedef ctf_dvdef_t * ctf_dvdef_ref;
typedef ctf_dtdef_t * ctf_dtdef_ref;

/* Location information for CTF Types and CTF Variables.  */

typedef struct GTY (()) ctf_srcloc
{
  const char * ctsloc_file;
  unsigned int ctsloc_line;
  unsigned int ctsloc_col;
} ctf_srcloc_t;

typedef ctf_srcloc_t * ctf_srcloc_ref;

/* Helper enum and api for the GTY machinery to work on union dtu_d.  */

enum ctf_dtu_d_union_enum {
  CTF_DTU_D_MEMBERS,
  CTF_DTU_D_ARRAY,
  CTF_DTU_D_ENCODING,
  CTF_DTU_D_ARGUMENTS,
  CTF_DTU_D_SLICE
};

enum ctf_dtu_d_union_enum
ctf_dtu_d_union_selector (ctf_dtdef_ref);

struct ctfc_dtd_hasher : ggc_ptr_hash <ctf_dtdef_t>
{
  typedef ctf_dtdef_ref compare_type;

  static hashval_t hash (ctf_dtdef_ref);
  static bool equal (ctf_dtdef_ref, ctf_dtdef_ref);
};

inline hashval_t
ctfc_dtd_hasher::hash (ctf_dtdef_ref dtd)
{
  return htab_hash_pointer (dtd->dtd_key);
}

inline bool
ctfc_dtd_hasher::equal (ctf_dtdef_ref dtd, ctf_dtdef_ref dtd2)
{
  return (dtd->dtd_key == dtd2->dtd_key);
}

struct ctfc_dvd_hasher : ggc_ptr_hash <ctf_dvdef_t>
{
  typedef ctf_dvdef_ref compare_type;

  static hashval_t hash (ctf_dvdef_ref);
  static bool equal (ctf_dvdef_ref, ctf_dvdef_ref);
};

inline hashval_t
ctfc_dvd_hasher::hash (ctf_dvdef_ref dvd)
{
  return htab_hash_pointer (dvd->dvd_key);
}

inline bool
ctfc_dvd_hasher::equal (ctf_dvdef_ref dvd, ctf_dvdef_ref dvd2)
{
  return (dvd->dvd_key == dvd2->dvd_key);
}

/* CTF container structure.
   It is the context passed around when generating ctf debug info.  There is
   one container per translation unit.  */

typedef struct GTY (()) ctf_container
{
  /* CTF Preamble.  */
  unsigned short ctfc_magic;
  unsigned char ctfc_version;
  unsigned char ctfc_flags;
  uint32_t ctfc_cuname_offset;

  /* CTF types.  */
  hash_table <ctfc_dtd_hasher> * GTY (()) ctfc_types;
  /* CTF variables.  */
  hash_table <ctfc_dvd_hasher> * GTY (()) ctfc_vars;
  /* CTF variables to be ignored.  */
  hash_table <ctfc_dvd_hasher> * GTY (()) ctfc_ignore_vars;

  /* CTF string table.  */
  ctf_strtable_t ctfc_strtable;
  /* Auxilliary string table.  At this time, used for keeping func arg names
     for BTF.  */
  ctf_strtable_t ctfc_aux_strtable;

  uint64_t ctfc_num_types;
  uint64_t ctfc_num_stypes;
  uint64_t ctfc_num_global_funcs;
  uint64_t ctfc_num_global_objts;

  /* Number of vlen bytes - the variable length portion after ctf_type_t and
     ctf_stype_t in the CTF section.  This is used to calculate the offsets in
     the CTF header.  */
  uint64_t ctfc_num_vlen_bytes;

  /* Next CTF type id to assign.  */
  ctf_id_t ctfc_nextid;

  /* Specify an explicit length of 0 so that the GC marking routines steer
     clear of marking the CTF vars and CTF types twice. These lists below do
     not own the pointed to objects, they simply hold references to them.  */

  /* List of pre-processed CTF Variables.  CTF requires that the variables
     appear in the sorted order of their names.  */
  ctf_dvdef_t ** GTY ((length ("0"))) ctfc_vars_list;
  /* Count of pre-processed CTF Variables in the list.  */
  uint64_t ctfc_vars_list_count;
  /* List of pre-processed CTF types.  CTF requires that a shared type must
     appear before the type that uses it.  For the compiler, this means types
     are emitted in sorted order of their type IDs.  */
  ctf_dtdef_t ** GTY ((length ("0"))) ctfc_types_list;
  /* List of CTF function types for global functions.  The order of global
     function entries in the CTF funcinfo section is undefined by the
     compiler.  */
  ctf_dtdef_t ** GTY ((length ("0"))) ctfc_gfuncs_list;
  /* List of CTF variables at global scope.  The order of global object entries
     in the CTF objinfo section is undefined by the  compiler.  */
  ctf_dvdef_t ** GTY ((length ("0"))) ctfc_gobjts_list;

  /* Following members are for debugging only.  They do not add functional
     value to the task of CTF creation.  These can be cleaned up once CTF
     generation stabilizes.  */

  /* Keep a count of the number of bytes dumped in asm for debugging
     purposes.  */
  uint64_t ctfc_numbytes_asm;
   /* Total length of all strings in CTF.  */
  size_t ctfc_strlen;
  /* Total length of all strings in aux string table.  */
  size_t ctfc_aux_strlen;

} ctf_container_t;

/* Markers for which string table from the CTF container to use.  */

#define CTF_STRTAB 0	    /* CTF string table.  */
#define CTF_AUX_STRTAB 1    /* CTF auxilliary string table.  */

typedef ctf_container_t * ctf_container_ref;

extern GTY (()) ctf_container_ref tu_ctfc;

extern void ctfc_delete_container (ctf_container_ref);

/* If the next ctf type id is still set to the init value, no ctf records to
   report.  */
extern bool ctfc_is_empty_container (ctf_container_ref);

/* Get the total number of CTF types in the container.  */

extern unsigned int ctfc_get_num_ctf_types (ctf_container_ref);

/* Get the total number of CTF variables in the container.  */

extern unsigned int ctfc_get_num_ctf_vars (ctf_container_ref);

/* Get reference to the CTF string table or the CTF auxilliary
   string table.  */

extern ctf_strtable_t * ctfc_get_strtab (ctf_container_ref, int);

/* Get the length of the specified string table in the CTF container.  */

extern size_t ctfc_get_strtab_len (ctf_container_ref, int);

/* Get the number of bytes to represent the variable length portion of all CTF
   types in the CTF container.  */

extern size_t ctfc_get_num_vlen_bytes (ctf_container_ref);

/* The compiler demarcates whether types are visible at top-level scope or not.
   The only example so far of a type not visible at top-level scope is slices.
   CTF_ADD_NONROOT is used to indicate the latter.  */
#define	CTF_ADD_NONROOT	0	/* CTF type only visible in nested scope.  */
#define	CTF_ADD_ROOT	1	/* CTF type visible at top-level scope.  */

/* These APIs allow to initialize and finalize the CTF machinery and
   to add types to the CTF container associated to the current
   translation unit.  Used in dwarf2ctf.cc.  */

extern void ctf_init (void);
extern void ctf_output (const char * filename);
extern void ctf_finalize (void);

extern void btf_output (const char * filename);
extern void btf_init_postprocess (void);
extern void btf_finalize (void);

extern ctf_container_ref ctf_get_tu_ctfc (void);

extern bool ctf_type_exists (ctf_container_ref, dw_die_ref, ctf_id_t *);

extern void ctf_add_cuname (ctf_container_ref, const char *);

extern ctf_dtdef_ref ctf_dtd_lookup (const ctf_container_ref ctfc,
				     dw_die_ref die);
extern ctf_dvdef_ref ctf_dvd_lookup (const ctf_container_ref ctfc,
				     dw_die_ref die);
extern bool ctf_dvd_ignore_lookup (const ctf_container_ref ctfc,
				   dw_die_ref die);

extern const char * ctf_add_string (ctf_container_ref, const char *,
				    uint32_t *, int);

extern ctf_id_t ctf_add_reftype (ctf_container_ref, uint32_t, ctf_id_t,
				 uint32_t, dw_die_ref);
extern ctf_id_t ctf_add_enum (ctf_container_ref, uint32_t, const char *,
			      HOST_WIDE_INT, dw_die_ref);
extern ctf_id_t ctf_add_slice (ctf_container_ref, uint32_t, ctf_id_t,
			       uint32_t, uint32_t, dw_die_ref);
extern ctf_id_t ctf_add_float (ctf_container_ref, uint32_t, const char *,
			       const ctf_encoding_t *, dw_die_ref);
extern ctf_id_t ctf_add_integer (ctf_container_ref, uint32_t, const char *,
				 const ctf_encoding_t *, dw_die_ref);
extern ctf_id_t ctf_add_unknown (ctf_container_ref, uint32_t, const char *,
				 const ctf_encoding_t *, dw_die_ref);
extern ctf_id_t ctf_add_pointer (ctf_container_ref, uint32_t, ctf_id_t,
				 dw_die_ref);
extern ctf_id_t ctf_add_array (ctf_container_ref, uint32_t,
			       const ctf_arinfo_t *, dw_die_ref);
extern ctf_id_t ctf_add_forward (ctf_container_ref, uint32_t, const char *,
				 uint32_t, dw_die_ref);
extern ctf_id_t ctf_add_typedef (ctf_container_ref, uint32_t, const char *,
				 ctf_id_t, dw_die_ref);
extern ctf_id_t ctf_add_function (ctf_container_ref, uint32_t, const char *,
				  const ctf_funcinfo_t *, dw_die_ref, bool);
extern ctf_id_t ctf_add_sou (ctf_container_ref, uint32_t, const char *,
			     uint32_t, size_t, dw_die_ref);

extern int ctf_add_enumerator (ctf_container_ref, ctf_id_t, const char *,
			       HOST_WIDE_INT, dw_die_ref);
extern int ctf_add_member_offset (ctf_container_ref, dw_die_ref, const char *,
				  ctf_id_t, uint64_t);
extern int ctf_add_function_arg (ctf_container_ref, dw_die_ref,
				 const char *, ctf_id_t);
extern int ctf_add_variable (ctf_container_ref, const char *, ctf_id_t,
			     dw_die_ref, unsigned int, dw_die_ref);

extern ctf_id_t ctf_lookup_tree_type (ctf_container_ref, const tree);
extern ctf_id_t get_btf_id (ctf_id_t);

/* CTF section does not emit location information; at this time, location
   information is needed for BTF CO-RE use-cases.  */

extern int ctfc_get_dtd_srcloc (ctf_dtdef_ref, ctf_srcloc_ref);
extern int ctfc_get_dvd_srcloc (ctf_dvdef_ref, ctf_srcloc_ref);

#endif /* GCC_CTFC_H */
