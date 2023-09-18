/* Communication between GCC and libgomp.

   Copyright (C) 2014-2022 Free Software Foundation, Inc.

   Contributed by Mentor Embedded.

   This file is part of the GNU Offloading and Multi Processing Library
   (libgomp).

   Libgomp is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   Libgomp is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
   more details.

   Under Section 7 of GPL version 3, you are granted additional
   permissions described in the GCC Runtime Library Exception, version
   3.1, as published by the Free Software Foundation.

   You should have received a copy of the GNU General Public License and
   a copy of the GCC Runtime Library Exception along with this program;
   see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
   <http://www.gnu.org/licenses/>.  */

#ifndef GOMP_CONSTANTS_H
#define GOMP_CONSTANTS_H 1

/* Memory mapping types.  */

/* One byte.  */
#define GOMP_MAP_LAST			(1 << 8)

#define GOMP_MAP_FLAG_TO		(1 << 0)
#define GOMP_MAP_FLAG_FROM		(1 << 1)
/* Special map kinds, enumerated starting here.  */
#define GOMP_MAP_FLAG_SPECIAL_0		(1 << 2)
#define GOMP_MAP_FLAG_SPECIAL_1		(1 << 3)
#define GOMP_MAP_FLAG_SPECIAL_2		(1 << 4)
#define GOMP_MAP_FLAG_SPECIAL_3		(1 << 5)
#define GOMP_MAP_FLAG_SPECIAL_4		(1 << 6)
#define GOMP_MAP_FLAG_SPECIAL		(GOMP_MAP_FLAG_SPECIAL_1 \
					 | GOMP_MAP_FLAG_SPECIAL_0)
#define GOMP_MAP_DEEP_COPY		(GOMP_MAP_FLAG_SPECIAL_4 \
					 | GOMP_MAP_FLAG_SPECIAL_2)
/* This value indicates the map was created implicitly according to
   OpenMP rules.  */
#define GOMP_MAP_IMPLICIT		(GOMP_MAP_FLAG_SPECIAL_3 \
					 | GOMP_MAP_FLAG_SPECIAL_4)
/* Mask for entire set of special map kind bits.  */
#define GOMP_MAP_FLAG_SPECIAL_BITS	(GOMP_MAP_FLAG_SPECIAL_0 \
					 | GOMP_MAP_FLAG_SPECIAL_1 \
					 | GOMP_MAP_FLAG_SPECIAL_2 \
					 | GOMP_MAP_FLAG_SPECIAL_3 \
					 | GOMP_MAP_FLAG_SPECIAL_4)
/* Flag to force a specific behavior (or else, trigger a run-time error).  */
#define GOMP_MAP_FLAG_FORCE		(1 << 7)

enum gomp_map_kind
  {
    /* If not already present, allocate.  */
    GOMP_MAP_ALLOC =			0,
    /* ..., and copy to device.  */
    GOMP_MAP_TO =			(GOMP_MAP_ALLOC | GOMP_MAP_FLAG_TO),
    /* ..., and copy from device.  */
    GOMP_MAP_FROM =			(GOMP_MAP_ALLOC | GOMP_MAP_FLAG_FROM),
    /* ..., and copy to and from device.  */
    GOMP_MAP_TOFROM =			(GOMP_MAP_TO | GOMP_MAP_FROM),
    /* The following kind is an internal only map kind, used for pointer based
       array sections.  OMP_CLAUSE_SIZE for these is not the pointer size,
       which is implicitly POINTER_SIZE_UNITS, but the bias.  */
    GOMP_MAP_POINTER =			(GOMP_MAP_FLAG_SPECIAL_0 | 0),
    /* Also internal, behaves like GOMP_MAP_TO, but additionally any
       GOMP_MAP_POINTER records consecutive after it which have addresses
       falling into that range will not be ignored if GOMP_MAP_TO_PSET wasn't
       mapped already.
       For OpenACC attach operations (e.g. copyin of struct members),
       GOMP_MAP_TO_PSET is followed by a single GOMP_MAP_ATTACH mapping
       instead.  */
    GOMP_MAP_TO_PSET =			(GOMP_MAP_FLAG_SPECIAL_0 | 1),
    /* Must already be present.  */
    GOMP_MAP_FORCE_PRESENT =		(GOMP_MAP_FLAG_SPECIAL_0 | 2),
    /* Deallocate a mapping, without copying from device.  */
    GOMP_MAP_DELETE =			(GOMP_MAP_FLAG_SPECIAL_0 | 3),
    /* Is a device pointer.  OMP_CLAUSE_SIZE for these is unused; is implicitly
       POINTER_SIZE_UNITS.  */
    GOMP_MAP_FORCE_DEVICEPTR =		(GOMP_MAP_FLAG_SPECIAL_1 | 0),
    /* OpenACC device_resident.  */
    GOMP_MAP_DEVICE_RESIDENT =		(GOMP_MAP_FLAG_SPECIAL_1 | 1),
    /* OpenACC link.  */
    GOMP_MAP_LINK =			(GOMP_MAP_FLAG_SPECIAL_1 | 2),
    /* Use device data if present, fall back to host address otherwise.  */
    GOMP_MAP_IF_PRESENT =		(GOMP_MAP_FLAG_SPECIAL_1 | 3),
    /* Do not map, copy bits for firstprivate instead.  */
    GOMP_MAP_FIRSTPRIVATE =		(GOMP_MAP_FLAG_SPECIAL | 0),
    /* Similarly, but store the value in the pointer rather than
       pointed by the pointer.  */
    GOMP_MAP_FIRSTPRIVATE_INT =		(GOMP_MAP_FLAG_SPECIAL | 1),
    /* Pointer translate host address into device address and copy that
       back to host.  */
    GOMP_MAP_USE_DEVICE_PTR =		(GOMP_MAP_FLAG_SPECIAL | 2),
    /* Allocate a zero length array section.  Prefer next non-zero length
       mapping over previous non-zero length mapping over zero length mapping
       at the address.  If not already mapped, do nothing (and pointer translate
       to NULL).  */
    GOMP_MAP_ZERO_LEN_ARRAY_SECTION = 	(GOMP_MAP_FLAG_SPECIAL | 3),
    /* Allocate.  */
    GOMP_MAP_FORCE_ALLOC =		(GOMP_MAP_FLAG_FORCE | GOMP_MAP_ALLOC),
    /* ..., and copy to device.  */
    GOMP_MAP_FORCE_TO =			(GOMP_MAP_FLAG_FORCE | GOMP_MAP_TO),
    /* ..., and copy from device.  */
    GOMP_MAP_FORCE_FROM =		(GOMP_MAP_FLAG_FORCE | GOMP_MAP_FROM),
    /* ..., and copy to and from device.  */
    GOMP_MAP_FORCE_TOFROM =		(GOMP_MAP_FLAG_FORCE | GOMP_MAP_TOFROM),
    /* Like GOMP_MAP_USE_DEVICE_PTR above, translate a host to a device
       address.  If translation fails because the target is not mapped,
       continue using the host address. */
    GOMP_MAP_USE_DEVICE_PTR_IF_PRESENT = (GOMP_MAP_FLAG_SPECIAL_2 | 0),
    /* If not already present, allocate.  And unconditionally copy to
       device.  */
    GOMP_MAP_ALWAYS_TO =		(GOMP_MAP_FLAG_SPECIAL_2 | GOMP_MAP_TO),
    /* If not already present, allocate.  And unconditionally copy from
       device.  */
    GOMP_MAP_ALWAYS_FROM =		(GOMP_MAP_FLAG_SPECIAL_2
					 | GOMP_MAP_FROM),
    /* If not already present, allocate.  And unconditionally copy to and from
       device.  */
    GOMP_MAP_ALWAYS_TOFROM =		(GOMP_MAP_FLAG_SPECIAL_2
					 | GOMP_MAP_TOFROM),
    /* Map a sparse struct; the address is the base of the structure, alignment
       it's required alignment, and size is the number of adjacent entries
       that belong to the struct.  The adjacent entries should be sorted by
       increasing address, so it is easy to determine lowest needed address
       (address of the first adjacent entry) and highest needed address
       (address of the last adjacent entry plus its size).  */
    GOMP_MAP_STRUCT =			(GOMP_MAP_FLAG_SPECIAL_2
					 | GOMP_MAP_FLAG_SPECIAL | 0),
    /* On a location of a pointer/reference that is assumed to be already mapped
       earlier, store the translated address of the preceeding mapping.
       No refcount is bumped by this, and the store is done unconditionally.  */
    GOMP_MAP_ALWAYS_POINTER =		(GOMP_MAP_FLAG_SPECIAL_2
					 | GOMP_MAP_FLAG_SPECIAL | 1),
    /* Like GOMP_MAP_POINTER, but allow zero-length array section, i.e. set to
       NULL if target is not mapped.  */
    GOMP_MAP_POINTER_TO_ZERO_LENGTH_ARRAY_SECTION
      =					(GOMP_MAP_FLAG_SPECIAL_2
					 | GOMP_MAP_FLAG_SPECIAL | 2),
    /* Forced deallocation of zero length array section.  */
    GOMP_MAP_DELETE_ZERO_LEN_ARRAY_SECTION
      =					(GOMP_MAP_FLAG_SPECIAL_2
					 | GOMP_MAP_FLAG_SPECIAL | 3),
    /* Decrement usage count and deallocate if zero.  */
    GOMP_MAP_RELEASE =			(GOMP_MAP_FLAG_SPECIAL_2
					 | GOMP_MAP_DELETE),
    /* The attach/detach mappings below use the OMP_CLAUSE_SIZE field as a
       bias.  This will typically be zero, except when mapping an array slice
       with a non-zero base.  In that case the bias will indicate the
       (positive) difference between the start of the actual mapped data and
       the "virtual" origin of the array.
       In OpenACC, attach a pointer to a mapped struct field.  */
    GOMP_MAP_ATTACH =			(GOMP_MAP_DEEP_COPY | 0),
    /* In OpenACC, detach a pointer to a mapped struct field.  */
    GOMP_MAP_DETACH =			(GOMP_MAP_DEEP_COPY | 1),
    /* In OpenACC, detach a pointer to a mapped struct field.  */
    GOMP_MAP_FORCE_DETACH =		(GOMP_MAP_DEEP_COPY
					 | GOMP_MAP_FLAG_FORCE | 1),

    /* Like GOMP_MAP_ATTACH, but allow attaching to zero-length array sections
       (i.e. set to NULL when array section is not mapped) Currently only used
       by OpenMP.  */
    GOMP_MAP_ATTACH_ZERO_LENGTH_ARRAY_SECTION
      =					(GOMP_MAP_DEEP_COPY | 2),

    /* Internal to GCC, not used in libgomp.  */
    /* Do not map, but pointer assign a pointer instead.  */
    GOMP_MAP_FIRSTPRIVATE_POINTER =	(GOMP_MAP_LAST | 1),
    /* Do not map, but pointer assign a reference instead.  */
    GOMP_MAP_FIRSTPRIVATE_REFERENCE =	(GOMP_MAP_LAST | 2),
    /* An attach or detach operation.  Rewritten to the appropriate type during
       gimplification, depending on directive (i.e. "enter data" or
       parallel/kernels region vs. "exit data").  */
    GOMP_MAP_ATTACH_DETACH =		(GOMP_MAP_LAST | 3)
  };

#define GOMP_MAP_COPY_TO_P(X) \
  (!((X) & GOMP_MAP_FLAG_SPECIAL) \
   && ((X) & GOMP_MAP_FLAG_TO))

#define GOMP_MAP_COPY_FROM_P(X) \
  (!((X) & GOMP_MAP_FLAG_SPECIAL) \
   && ((X) & GOMP_MAP_FLAG_FROM))

#define GOMP_MAP_ALWAYS_POINTER_P(X) \
  ((X) == GOMP_MAP_ALWAYS_POINTER)

#define GOMP_MAP_POINTER_P(X) \
  ((X) == GOMP_MAP_POINTER \
   || (X) == GOMP_MAP_POINTER_TO_ZERO_LENGTH_ARRAY_SECTION)

#define GOMP_MAP_ALWAYS_TO_P(X) \
  (((X) == GOMP_MAP_ALWAYS_TO) || ((X) == GOMP_MAP_ALWAYS_TOFROM))

#define GOMP_MAP_ALWAYS_FROM_P(X) \
  (((X) == GOMP_MAP_ALWAYS_FROM) || ((X) == GOMP_MAP_ALWAYS_TOFROM))

#define GOMP_MAP_ALWAYS_P(X) \
  (GOMP_MAP_ALWAYS_TO_P (X) || ((X) == GOMP_MAP_ALWAYS_FROM))

#define GOMP_MAP_IMPLICIT_P(X) \
  (((X) & GOMP_MAP_FLAG_SPECIAL_BITS) == GOMP_MAP_IMPLICIT)


/* Asynchronous behavior.  Keep in sync with
   libgomp/{openacc.h,openacc.f90,openacc_lib.h}:acc_async_t.  */

#define GOMP_ASYNC_NOVAL		-1
#define GOMP_ASYNC_SYNC			-2


/* Device codes.  Keep in sync with
   libgomp/{openacc.h,openacc.f90,openacc_lib.h}:acc_device_t as well as
   libgomp/libgomp-plugin.h.  */
#define GOMP_DEVICE_NONE		0
#define GOMP_DEVICE_DEFAULT		1
#define GOMP_DEVICE_HOST		2
/* #define GOMP_DEVICE_HOST_NONSHM	3 removed.  */
#define GOMP_DEVICE_NOT_HOST		4
#define GOMP_DEVICE_NVIDIA_PTX		5
#define GOMP_DEVICE_INTEL_MIC		6
#define GOMP_DEVICE_HSA			7
#define GOMP_DEVICE_GCN			8

#define GOMP_DEVICE_ICV			-1
#define GOMP_DEVICE_HOST_FALLBACK	-2

/* GOMP_task/GOMP_taskloop* flags argument.  */
#define GOMP_TASK_FLAG_UNTIED		(1 << 0)
#define GOMP_TASK_FLAG_FINAL		(1 << 1)
#define GOMP_TASK_FLAG_MERGEABLE	(1 << 2)
#define GOMP_TASK_FLAG_DEPEND		(1 << 3)
#define GOMP_TASK_FLAG_PRIORITY		(1 << 4)
#define GOMP_TASK_FLAG_UP		(1 << 8)
#define GOMP_TASK_FLAG_GRAINSIZE	(1 << 9)
#define GOMP_TASK_FLAG_IF		(1 << 10)
#define GOMP_TASK_FLAG_NOGROUP		(1 << 11)
#define GOMP_TASK_FLAG_REDUCTION	(1 << 12)
#define GOMP_TASK_FLAG_DETACH		(1 << 13)
#define GOMP_TASK_FLAG_STRICT		(1 << 14)

/* GOMP_target{_ext,update_ext,enter_exit_data} flags argument.  */
#define GOMP_TARGET_FLAG_NOWAIT		(1 << 0)
#define GOMP_TARGET_FLAG_EXIT_DATA	(1 << 1)
/* Internal to libgomp.  */
#define GOMP_TARGET_FLAG_UPDATE		(1U << 31)


/* OpenACC construct flags.  */

/* Force host fallback execution.  */
#define GOACC_FLAG_HOST_FALLBACK	(1 << 0)

/* For legacy reasons, in the ABI, the GOACC_FLAGs are encoded as an inverted
   bitmask.  */
#define GOACC_FLAGS_MARSHAL_OP		BIT_NOT_EXPR
#define GOACC_FLAGS_UNMARSHAL(X)	(~(X))


/* Versions of libgomp and device-specific plugins.  GOMP_VERSION
   should be incremented whenever an ABI-incompatible change is introduced
   to the plugin interface defined in libgomp/libgomp.h.  */
#define GOMP_VERSION	1
#define GOMP_VERSION_NVIDIA_PTX 1
#define GOMP_VERSION_INTEL_MIC 0
#define GOMP_VERSION_GCN 2

#define GOMP_VERSION_PACK(LIB, DEV) (((LIB) << 16) | (DEV))
#define GOMP_VERSION_LIB(PACK) (((PACK) >> 16) & 0xffff)
#define GOMP_VERSION_DEV(PACK) ((PACK) & 0xffff)

#define GOMP_DIM_GANG	0
#define GOMP_DIM_WORKER	1
#define GOMP_DIM_VECTOR	2
#define GOMP_DIM_MAX	3
#define GOMP_DIM_MASK(X) (1u << (X))

/* Varadic launch arguments.  End of list is marked by a zero.  */
#define GOMP_LAUNCH_DIM		1  /* Launch dimensions, op = mask */
#define GOMP_LAUNCH_ASYNC	2  /* Async, op = cst val if not MAX  */
#define GOMP_LAUNCH_WAIT	3  /* Waits, op = num waits.  */
#define GOMP_LAUNCH_CODE_SHIFT	28
#define GOMP_LAUNCH_DEVICE_SHIFT 16
#define GOMP_LAUNCH_OP_SHIFT 0
#define GOMP_LAUNCH_PACK(CODE,DEVICE,OP)	\
  (((CODE) << GOMP_LAUNCH_CODE_SHIFT)		\
   | ((DEVICE) << GOMP_LAUNCH_DEVICE_SHIFT)	\
   | ((OP) << GOMP_LAUNCH_OP_SHIFT))
#define GOMP_LAUNCH_CODE(X) (((X) >> GOMP_LAUNCH_CODE_SHIFT) & 0xf)
#define GOMP_LAUNCH_DEVICE(X) (((X) >> GOMP_LAUNCH_DEVICE_SHIFT) & 0xfff)
#define GOMP_LAUNCH_OP(X) (((X) >> GOMP_LAUNCH_OP_SHIFT) & 0xffff)
#define GOMP_LAUNCH_OP_MAX 0xffff

/* Bitmask to apply in order to find out the intended device of a target
   argument.  */
#define GOMP_TARGET_ARG_DEVICE_MASK		((1 << 7) - 1)
/* The target argument is significant for all devices.  */
#define GOMP_TARGET_ARG_DEVICE_ALL		0

/* Flag set when the subsequent element in the device-specific argument
   values.  */
#define GOMP_TARGET_ARG_SUBSEQUENT_PARAM	(1 << 7)

/* Bitmask to apply to a target argument to find out the value identifier.  */
#define GOMP_TARGET_ARG_ID_MASK			(((1 << 8) - 1) << 8)
/* Target argument index of NUM_TEAMS.  */
#define GOMP_TARGET_ARG_NUM_TEAMS		(1 << 8)
/* Target argument index of THREAD_LIMIT.  */
#define GOMP_TARGET_ARG_THREAD_LIMIT		(2 << 8)

/* If the value is directly embeded in target argument, it should be a 16-bit
   at most and shifted by this many bits.  */
#define GOMP_TARGET_ARG_VALUE_SHIFT		16

/* Dependence types in omp_depend_t objects.  */
#define GOMP_DEPEND_IN			1
#define GOMP_DEPEND_OUT			2
#define GOMP_DEPEND_INOUT		3
#define GOMP_DEPEND_MUTEXINOUTSET	4

/* HSA specific data structures.  */

/* Identifiers of device-specific target arguments.  */
#define GOMP_TARGET_ARG_HSA_KERNEL_ATTRIBUTES	(1 << 8)

#endif
