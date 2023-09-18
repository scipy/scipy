/* SFrame format description.
   Copyright (C) 2022-2023 Free Software Foundation, Inc.

   This file is part of libsframe.

   libsframe is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 3, or (at your option) any later
   version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; see the file COPYING.  If not see
   <http://www.gnu.org/licenses/>.  */

#ifndef	_SFRAME_H
#define	_SFRAME_H

#include <sys/types.h>
#include <limits.h>
#include <stdint.h>

#include "ansidecl.h"

#ifdef	__cplusplus
extern "C"
{
#endif

/* SFrame format.

   SFrame format is a simple format to represent the information needed
   for vanilla virtual stack unwinding.  SFrame format keeps track of the
   minimal necessary information needed for stack unwinding:
     - Canonical Frame Address (CFA)
     - Frame Pointer (FP)
     - Return Address (RA)

   The SFrame section itself has the following structure:

       +--------+------------+---------+
       |  file  |  function  | frame   |
       | header | descriptor |  row    |
       |        |   entries  | entries |
       +--------+------------+---------+

   The file header stores a magic number and version information, flags, and
   the byte offset of each of the sections relative to the end of the header
   itself.  The file header also specifies the total number of Function
   Descriptor Entries, Frame Row Entries and length of the FRE sub-section.

   Following the header is a list of Function Descriptor Entries (FDEs).
   This list may be sorted if the flags in the file header indicate it to be
   so.  The sort order, if applicable, is the order of functions in the
   .text.* sections in the resulting binary artifact.  Each Function
   Descriptor Entry specifies the start PC of a function, the size in bytes
   of the function and an offset to its first Frame Row Entry (FRE).  Each FDE
   additionally also specifies the type of FRE it uses to encode the unwind
   information.

   Next, the Frame Row Entry section is a list of variable size records,
   each of which represent SFrame unwind information for a set of PCs.  A
   singular Frame Row Entry is a self-sufficient record with information on
   how to virtually unwind the stack for the applicable set of PCs.

   */


/* SFrame format versions.  */
#define SFRAME_VERSION_1	1
/* SFrame magic number.  */
#define SFRAME_MAGIC		0xdee2
/* Current version of SFrame format.  */
#define SFRAME_VERSION	SFRAME_VERSION_1

/* Various flags for SFrame.  */

/* Function Descriptor Entries are sorted on PC.  */
#define SFRAME_F_FDE_SORTED	0x1
/* Frame-pointer based unwinding.  */
#define SFRAME_F_FRAME_POINTER 0x2

#define SFRAME_CFA_FIXED_FP_INVALID 0
#define SFRAME_CFA_FIXED_RA_INVALID 0

/* Supported ABIs/Arch.  */
#define SFRAME_ABI_AARCH64_ENDIAN_BIG      1 /* AARCH64 big endian.  */
#define SFRAME_ABI_AARCH64_ENDIAN_LITTLE   2 /* AARCH64 little endian.  */
#define SFRAME_ABI_AMD64_ENDIAN_LITTLE     3 /* AMD64 little endian.  */

/* SFrame FRE types.  */
#define SFRAME_FRE_TYPE_ADDR1	0
#define SFRAME_FRE_TYPE_ADDR2	1
#define SFRAME_FRE_TYPE_ADDR4	2

/* SFrame Function Descriptor Entry types.

   The SFrame format has two possible representations for functions.  The
   choice of which type to use is made according to the instruction patterns
   in the relevant program stub.

   An SFrame FDE of type SFRAME_FDE_TYPE_PCINC is an indication
   that the PCs in the FREs should be treated as increments in bytes.  This is
   used for a bulk of the executable code of a program, which contains
   instructions with no specific pattern.

   An SFrame FDE of type SFRAME_FDE_TYPE_PCMASK is an indication
   that the PCs in the FREs should be treated as masks.  This type is useful
   for the cases when a small pattern of instructions in a program stub is
   repeatedly to cover a specific functionality.  Typical usescases are pltN
   entries, trampolines etc.  */

/* Unwinders perform a (PC >= FRE_START_ADDR) to look up a matching FRE.  */
#define SFRAME_FDE_TYPE_PCINC   0
/* Unwinders perform a (PC & FRE_START_ADDR_AS_MASK >= FRE_START_ADDR_AS_MASK)
   to look up a matching FRE.  */
#define SFRAME_FDE_TYPE_PCMASK  1

typedef struct sframe_preamble
{
  uint16_t sfp_magic;	/* Magic number (SFRAME_MAGIC).  */
  uint8_t sfp_version;	/* Data format version number (SFRAME_VERSION).  */
  uint8_t sfp_flags;	/* Flags.  */
} ATTRIBUTE_PACKED sframe_preamble;

typedef struct sframe_header
{
  sframe_preamble sfh_preamble;
  /* Information about the arch (endianness) and ABI.  */
  uint8_t sfh_abi_arch;
  /* Offset for the Frame Pointer (FP) from CFA may be fixed for some
     ABIs (e.g, in AMD64 when -fno-omit-frame-pointer is used).  When fixed,
     this field specifies the fixed stack frame offset and the individual
     FREs do not need to track it.  When not fixed, it is set to
     SFRAME_CFA_FIXED_FP_INVALID, and the individual FREs may provide
     the applicable stack frame offset, if any.  */
  int8_t sfh_cfa_fixed_fp_offset;
  /* Offset for the Return Address from CFA is fixed for some ABIs
     (e.g., AMD64 has it as CFA-8).  When fixed, the header specifies the
     fixed stack frame offset and the individual FREs do not track it.  When
     not fixed, it is set to SFRAME_CFA_FIXED_RA_INVALID, and individual
     FREs provide the applicable stack frame offset, if any.  */
  int8_t sfh_cfa_fixed_ra_offset;
  /* Number of bytes making up the auxilliary header, if any.
     Some ABI/arch, in the future, may use this space for extending the
     information in SFrame header.  Auxilliary header is contained in
     bytes sequentially following the sframe_header.  */
  uint8_t sfh_auxhdr_len;
  /* Number of SFrame FDEs in this SFrame section.  */
  uint32_t sfh_num_fdes;
  /* Number of SFrame Frame Row Entries.  */
  uint32_t sfh_num_fres;
  /* Number of bytes in the SFrame Frame Row Entry section. */
  uint32_t sfh_fre_len;
  /* Offset of SFrame Function Descriptor Entry section.  */
  uint32_t sfh_fdeoff;
  /* Offset of SFrame Frame Row Entry section.  */
  uint32_t sfh_freoff;
} ATTRIBUTE_PACKED sframe_header;

#define SFRAME_V1_HDR_SIZE(sframe_hdr)	\
  ((sizeof (sframe_header) + (sframe_hdr).sfh_auxhdr_len))

/* Two possible keys for executable (instruction) pointers signing.  */
#define SFRAME_AARCH64_PAUTH_KEY_A    0 /* Key A.  */
#define SFRAME_AARCH64_PAUTH_KEY_B    1 /* Key B.  */

typedef struct sframe_func_desc_entry
{
  /* Function start address.  Encoded as a signed offset, relative to the
     beginning of the current FDE.  */
  int32_t sfde_func_start_address;
  /* Size of the function in bytes.  */
  uint32_t sfde_func_size;
  /* Offset of the first SFrame Frame Row Entry of the function, relative to the
     beginning of the SFrame Frame Row Entry sub-section.  */
  uint32_t sfde_func_start_fre_off;
  /* Number of frame row entries for the function.  */
  uint32_t sfde_func_num_fres;
  /* Additional information for deciphering the unwind information for the
     function.
     - 4-bits: Identify the FRE type used for the function.
     - 1-bit: Identify the FDE type of the function - mask or inc.
     - 1-bit: PAC authorization A/B key (aarch64).
     - 2-bits: Unused.
     ------------------------------------------------------------------------
     |     Unused    |  PAC auth A/B key (aarch64) |  FDE type |   FRE type   |
     |               |        Unused (amd64)       |           |              |
     ------------------------------------------------------------------------
     8               6                             5           4              0     */
  uint8_t sfde_func_info;
} ATTRIBUTE_PACKED sframe_func_desc_entry;

/* Macros to compose and decompose function info in FDE.  */

/* Note: Set PAC auth key to SFRAME_AARCH64_PAUTH_KEY_A by default.  */
#define SFRAME_V1_FUNC_INFO(fde_type, fre_enc_type) \
  (((SFRAME_AARCH64_PAUTH_KEY_A & 0x1) << 5) | \
   (((fde_type) & 0x1) << 4) | ((fre_enc_type) & 0xf))

#define SFRAME_V1_FUNC_FRE_TYPE(data)	  ((data) & 0xf)
#define SFRAME_V1_FUNC_FDE_TYPE(data)	  (((data) >> 4) & 0x1)
#define SFRAME_V1_FUNC_PAUTH_KEY(data)	  (((data) >> 5) & 0x1)

/* Set the pauth key as indicated.  */
#define SFRAME_V1_FUNC_INFO_UPDATE_PAUTH_KEY(pauth_key, fde_info) \
  ((((pauth_key) & 0x1) << 5) | ((fde_info) & 0xdf))

/* Size of stack frame offsets in an SFrame Frame Row Entry.  A single
   SFrame FRE has all offsets of the same size.  Offset size may vary
   across frame row entries.  */
#define SFRAME_FRE_OFFSET_1B	  0
#define SFRAME_FRE_OFFSET_2B	  1
#define SFRAME_FRE_OFFSET_4B	  2

/* An SFrame Frame Row Entry can be SP or FP based.  */
#define SFRAME_BASE_REG_FP	0
#define SFRAME_BASE_REG_SP	1

/* The index at which a specific offset is presented in the variable length
   bytes of an FRE.  */
#define SFRAME_FRE_CFA_OFFSET_IDX   0
/* The RA stack offset, if present, will always be at index 1 in the variable
   length bytes of the FRE.  */
#define SFRAME_FRE_RA_OFFSET_IDX    1
/* The FP stack offset may appear at offset 1 or 2, depending on the ABI as RA
   may or may not be tracked.  */
#define SFRAME_FRE_FP_OFFSET_IDX    2

typedef struct sframe_fre_info
{
  /* Information about
     - 1 bit: base reg for CFA
     - 4 bits: Number of offsets (N).  A value of upto 3 is allowed to track
     all three of CFA, FP and RA (fixed implicit order).
     - 2 bits: information about size of the offsets (S) in bytes.
     Valid values are SFRAME_FRE_OFFSET_1B, SFRAME_FRE_OFFSET_2B,
     SFRAME_FRE_OFFSET_4B
     - 1 bit: Mangled RA state bit (aarch64 only).
     ----------------------------------------------------------------------------------
     | Mangled-RA (aarch64) |  Size of offsets   |   Number of offsets    |   base_reg |
     |  Unused (amd64)      |                    |                        |            |
     ----------------------------------------------------------------------------------
     8                     7                    5                        1            0

     */
  uint8_t fre_info;
} sframe_fre_info;

/* Macros to compose and decompose FRE info.  */

/* Note: Set mangled_ra_p to zero by default.  */
#define SFRAME_V1_FRE_INFO(base_reg_id, offset_num, offset_size) \
  (((0 & 0x1) << 7) | (((offset_size) & 0x3) << 5) | \
   (((offset_num) & 0xf) << 1) | ((base_reg_id) & 0x1))

/* Set the mangled_ra_p bit as indicated.  */
#define SFRAME_V1_FRE_INFO_UPDATE_MANGLED_RA_P(mangled_ra_p, fre_info) \
  ((((mangled_ra_p) & 0x1) << 7) | ((fre_info) & 0x7f))

#define SFRAME_V1_FRE_CFA_BASE_REG_ID(data)	  ((data) & 0x1)
#define SFRAME_V1_FRE_OFFSET_COUNT(data)	  (((data) >> 1) & 0xf)
#define SFRAME_V1_FRE_OFFSET_SIZE(data)		  (((data) >> 5) & 0x3)
#define SFRAME_V1_FRE_MANGLED_RA_P(data)	  (((data) >> 7) & 0x1)

/* SFrame Frame Row Entry definitions.

   Used for both AMD64 and AARCH64.

   An SFrame Frame Row Entry is a self-sufficient record containing SFrame
   unwind info for a range of addresses, starting at the specified offset in
   the function.  Each SFrame Frame Row Entry is followed by S*N bytes, where:
     S is the size of the stack frame offset for the FRE, and
     N is the number of stack frame offsets in the FRE

   The offsets are interpreted in order as follows:

    offset1 (interpreted as CFA = BASE_REG + offset1)

    if RA is being tracked
      offset2 (interpreted as RA = CFA + offset2)
      if FP is being tracked
	offset3 (intrepreted as FP = CFA + offset2)
      fi
    else
      if FP is being tracked
	offset2 (intrepreted as FP = CFA + offset2)
      fi
    fi
*/

/* Used when SFRAME_FRE_TYPE_ADDR1 is specified as FRE type.  */
typedef struct sframe_frame_row_entry_addr1
{
  /* Start address of the frame row entry.  Encoded as an 1-byte unsigned
     offset, relative to the start address of the function.  */
  uint8_t sfre_start_address;
  sframe_fre_info sfre_info;
} ATTRIBUTE_PACKED sframe_frame_row_entry_addr1;

/* Upper limit of start address in sframe_frame_row_entry_addr1
   is 0x100 (not inclusive).  */
#define SFRAME_FRE_TYPE_ADDR1_LIMIT   \
  (1ULL << ((SFRAME_FRE_TYPE_ADDR1 + 1) * 8))

/* Used when SFRAME_FRE_TYPE_ADDR2 is specified as FRE type.  */
typedef struct sframe_frame_row_entry_addr2
{
  /* Start address of the frame row entry.  Encoded as an 2-byte unsigned
     offset, relative to the start address of the function.  */
  uint16_t sfre_start_address;
  sframe_fre_info sfre_info;
} ATTRIBUTE_PACKED sframe_frame_row_entry_addr2;

/* Upper limit of start address in sframe_frame_row_entry_addr2
   is 0x10000 (not inclusive).  */
#define SFRAME_FRE_TYPE_ADDR2_LIMIT   \
  (1ULL << ((SFRAME_FRE_TYPE_ADDR2 * 2) * 8))

/* Used when SFRAME_FRE_TYPE_ADDR4 is specified as FRE type.  */
typedef struct sframe_frame_row_entry_addr4
{
  /* Start address of the frame row entry.  Encoded as a 4-byte unsigned
     offset, relative to the start address of the function.  */
  uint32_t sfre_start_address;
  sframe_fre_info sfre_info;
} ATTRIBUTE_PACKED sframe_frame_row_entry_addr4;

/* Upper limit of start address in sframe_frame_row_entry_addr2
   is 0x100000000 (not inclusive).  */
#define SFRAME_FRE_TYPE_ADDR4_LIMIT   \
  (1ULL << ((SFRAME_FRE_TYPE_ADDR4 * 2) * 8))

#ifdef	__cplusplus
}
#endif

#endif				/* _SFRAME_H */
