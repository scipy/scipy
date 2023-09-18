/* Public API to SFrame.

   Copyright (C) 2022-2023 Free Software Foundation, Inc.

   This file is part of libsframe.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

#ifndef	_SFRAME_API_H
#define	_SFRAME_API_H

#include <sframe.h>
#include <stdbool.h>

#ifdef	__cplusplus
extern "C"
{
#endif

typedef struct sframe_decoder_ctx sframe_decoder_ctx;
typedef struct sframe_encoder_ctx sframe_encoder_ctx;

#define MAX_OFFSET_BYTES (SFRAME_FRE_OFFSET_4B * 2 * 3)

/* User interfacing SFrame Row Entry.
   An abstraction provided by libsframe so the consumer is decoupled from
   the binary format representation of the same.

   The members are best ordered such that they are aligned at their natural
   boundaries.  This helps avoid usage of undesirable misaligned memory
   accesses.  See PR libsframe/29856.  */

typedef struct sframe_frame_row_entry
{
  uint32_t fre_start_addr;
  unsigned char fre_offsets[MAX_OFFSET_BYTES];
  unsigned char fre_info;
} sframe_frame_row_entry;

#define SFRAME_ERR ((int) -1)

/* This macro holds information about all the available SFrame
   errors.  It is used to form both an enum holding all the error
   constants, and also the error strings themselves.  To use, define
   _SFRAME_FIRST and _SFRAME_ITEM to expand as you like, then
   mention the macro name.  See the enum after this for an example.  */
#define _SFRAME_ERRORS \
  _SFRAME_FIRST (SFRAME_ERR_VERSION_INVAL, "SFrame version not supported.") \
  _SFRAME_ITEM (SFRAME_ERR_NOMEM, "Out of Memory.") \
  _SFRAME_ITEM (SFRAME_ERR_INVAL, "Corrupt SFrame.") \
  _SFRAME_ITEM (SFRAME_ERR_BUF_INVAL, "Buffer does not contain SFrame data.") \
  _SFRAME_ITEM (SFRAME_ERR_DCTX_INVAL, "Corrupt SFrame decoder.") \
  _SFRAME_ITEM (SFRAME_ERR_ECTX_INVAL, "Corrupt SFrame encoder.") \
  _SFRAME_ITEM (SFRAME_ERR_FDE_INVAL, "Corrput FDE.") \
  _SFRAME_ITEM (SFRAME_ERR_FRE_INVAL, "Corrupt FRE.") \
  _SFRAME_ITEM (SFRAME_ERR_FDE_NOTFOUND,"FDE not found.") \
  _SFRAME_ITEM (SFRAME_ERR_FDE_NOTSORTED, "FDEs not sorted.") \
  _SFRAME_ITEM (SFRAME_ERR_FRE_NOTFOUND,"FRE not found.") \
  _SFRAME_ITEM (SFRAME_ERR_FREOFFSET_NOPRESENT,"FRE offset not present.")

#define	SFRAME_ERR_BASE	2000	/* Base value for libsframe errnos.  */

enum
  {
#define _SFRAME_FIRST(NAME, STR) NAME = SFRAME_ERR_BASE
#define _SFRAME_ITEM(NAME, STR) , NAME
_SFRAME_ERRORS
#undef _SFRAME_ITEM
#undef _SFRAME_FIRST
  };

/* Count of SFrame errors.  */
#define SFRAME_ERR_NERR (SFRAME_ERR_FREOFFSET_NOPRESENT - SFRAME_ERR_BASE + 1)

/* Get the error message string.  */

extern const char *
sframe_errmsg (int error);

/* Create an FDE function info bye given an FRE_TYPE and an FDE_TYPE.  */

extern unsigned char
sframe_fde_create_func_info (unsigned int fre_type, unsigned int fde_type);

/* Gather the FRE type given the function size.  */

extern unsigned int
sframe_calc_fre_type (size_t func_size);

/* The SFrame Decoder.  */

/* Decode the specified SFrame buffer CF_BUF of size CF_SIZE and return the
   new SFrame decoder context.  Sets ERRP for the caller if any error.  */
extern sframe_decoder_ctx *
sframe_decode (const char *cf_buf, size_t cf_size, int *errp);

/* Free the decoder context.  */
extern void
sframe_decoder_free (sframe_decoder_ctx **dctx);

/* Get the size of the SFrame header from the decoder context DCTX.  */
extern unsigned int
sframe_decoder_get_hdr_size (sframe_decoder_ctx *dctx);

/* Get the SFrame's abi/arch info.  */
extern unsigned char
sframe_decoder_get_abi_arch (sframe_decoder_ctx *dctx);

/* Return the number of function descriptor entries in the SFrame decoder
   DCTX.  */
unsigned int
sframe_decoder_get_num_fidx (sframe_decoder_ctx *dctx);

/* Get the fixed FP offset from the decoder context DCTX.  */
extern int8_t
sframe_decoder_get_fixed_fp_offset (sframe_decoder_ctx *dctx);

/* Get the fixed RA offset from the decoder context DCTX.  */
extern int8_t
sframe_decoder_get_fixed_ra_offset (sframe_decoder_ctx *dctx);

/* Find the function descriptor entry which contains the specified address.  */
extern sframe_func_desc_entry *
sframe_get_funcdesc_with_addr (sframe_decoder_ctx *dctx,
			       int32_t addr, int *errp);

/* Find the SFrame Frame Row Entry which contains the PC.  Returns
   SFRAME_ERR if failure.  */

extern int
sframe_find_fre (sframe_decoder_ctx *ctx, int32_t pc,
		 sframe_frame_row_entry *frep);

/* Get the FRE_IDX'th FRE of the function at FUNC_IDX'th function
   index entry in the SFrame decoder CTX.  Returns error code as
   applicable.  */
extern int
sframe_decoder_get_fre (sframe_decoder_ctx *ctx,
			unsigned int func_idx,
			unsigned int fre_idx,
			sframe_frame_row_entry *fre);

/* Get the data (NUM_FRES, FUNC_START_ADDRESS) from the function
   descriptor entry at index I'th in the decoder CTX.  If failed,
   return error code.  */
extern int
sframe_decoder_get_funcdesc (sframe_decoder_ctx *ctx,
			     unsigned int i,
			     uint32_t *num_fres,
			     uint32_t *func_size,
			     int32_t *func_start_address,
			     unsigned char *func_info);

/* SFrame textual dump.  */
extern void
dump_sframe (sframe_decoder_ctx *decoder, uint64_t addr);

/* Get the base reg id from the FRE info.  Sets errp if fails.  */
extern unsigned int
sframe_fre_get_base_reg_id (sframe_frame_row_entry *fre, int *errp);

/* Get the CFA offset from the FRE.  If the offset is invalid, sets errp.  */
extern int32_t
sframe_fre_get_cfa_offset (sframe_decoder_ctx *dtcx,
			   sframe_frame_row_entry *fre, int *errp);

/* Get the FP offset from the FRE.  If the offset is invalid, sets errp.  */
extern int32_t
sframe_fre_get_fp_offset (sframe_decoder_ctx *dctx,
			  sframe_frame_row_entry *fre, int *errp);

/* Get the RA offset from the FRE.  If the offset is invalid, sets errp.  */
extern int32_t
sframe_fre_get_ra_offset (sframe_decoder_ctx *dctx,
			  sframe_frame_row_entry *fre, int *errp);

/* Get whether the RA is mangled.  */

extern bool
sframe_fre_get_ra_mangled_p (sframe_decoder_ctx *dctx,
			     sframe_frame_row_entry *fre, int *errp);

/* The SFrame Encoder.  */

/* Create an encoder context with the given SFrame format version VER, FLAGS
   and ABI information.  Sets errp if failure.  */
extern sframe_encoder_ctx *
sframe_encode (unsigned char ver, unsigned char flags, int abi,
	       int8_t fixed_fp_offset, int8_t fixed_ra_offset, int *errp);

/* Free the encoder context.  */
extern void
sframe_encoder_free (sframe_encoder_ctx **encoder);

/* Get the size of the SFrame header from the encoder ctx ENCODER.  */
extern unsigned int
sframe_encoder_get_hdr_size (sframe_encoder_ctx *encoder);

/* Get the abi/arch info from the SFrame encoder context CTX.  */
extern unsigned char
sframe_encoder_get_abi_arch (sframe_encoder_ctx *encoder);

/* Return the number of function descriptor entries in the SFrame encoder
   ENCODER.  */
extern unsigned int
sframe_encoder_get_num_fidx (sframe_encoder_ctx *encoder);

/* Add an FRE to function at FUNC_IDX'th function descriptor index entry in
   the encoder context.  */
extern int
sframe_encoder_add_fre (sframe_encoder_ctx *encoder,
			unsigned int func_idx,
			sframe_frame_row_entry *frep);

/* Add a new function descriptor entry with START_ADDR, FUNC_SIZE and NUM_FRES
   to the encoder.  */
extern int
sframe_encoder_add_funcdesc (sframe_encoder_ctx *encoder,
			     int32_t start_addr,
			     uint32_t func_size,
			     unsigned char func_info,
			     uint32_t num_fres);

/* Serialize the contents of the encoder and return the buffer.  ENCODED_SIZE
   is updated to the size of the buffer.  Sets ERRP if failure.  */
extern char  *
sframe_encoder_write (sframe_encoder_ctx *encoder,
		      size_t *encoded_size, int *errp);

#ifdef	__cplusplus
}
#endif

#endif				/* _SFRAME_API_H */
