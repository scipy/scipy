/* Copyright (C) 2005-2022 Free Software Foundation, Inc.
   Contributed by Richard Henderson <rth@redhat.com>.

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

#ifndef _OMP_H
#define _OMP_H 1

#if defined(__GNUC__) && _OPENMP >= 201811
# define __GOMP_DEPRECATED_5_0 __attribute__((__deprecated__))
#else
# define __GOMP_DEPRECATED_5_0
#endif

#if defined(__GNUC__) && _OPENMP >= 202011
# define __GOMP_DEPRECATED_5_1 __attribute__((__deprecated__))
#else
# define __GOMP_DEPRECATED_5_1
#endif

#ifndef _LIBGOMP_OMP_LOCK_DEFINED
#define _LIBGOMP_OMP_LOCK_DEFINED 1
/* These two structures get edited by the libgomp build process to 
   reflect the shape of the two types.  Their internals are private
   to the library.  */

typedef struct
{
  unsigned char _x[4] 
    __attribute__((__aligned__(4)));
} omp_lock_t;

typedef struct
{
  unsigned char _x[16] 
    __attribute__((__aligned__(8)));
} omp_nest_lock_t;
#endif

typedef enum omp_sched_t
{
  omp_sched_static = 1,
  omp_sched_dynamic = 2,
  omp_sched_guided = 3,
  omp_sched_auto = 4,
  omp_sched_monotonic = 0x80000000U
} omp_sched_t;

typedef enum omp_proc_bind_t
{
  omp_proc_bind_false = 0,
  omp_proc_bind_true = 1,
  omp_proc_bind_primary = 2,
  omp_proc_bind_master __GOMP_DEPRECATED_5_1
    = omp_proc_bind_primary,
  omp_proc_bind_close = 3,
  omp_proc_bind_spread = 4
} omp_proc_bind_t;

typedef enum omp_sync_hint_t
{
  omp_sync_hint_none = 0,
  omp_lock_hint_none __GOMP_DEPRECATED_5_0 = omp_sync_hint_none,
  omp_sync_hint_uncontended = 1,
  omp_lock_hint_uncontended __GOMP_DEPRECATED_5_0 = omp_sync_hint_uncontended,
  omp_sync_hint_contended = 2,
  omp_lock_hint_contended __GOMP_DEPRECATED_5_0 = omp_sync_hint_contended,
  omp_sync_hint_nonspeculative = 4,
  omp_lock_hint_nonspeculative __GOMP_DEPRECATED_5_0
    = omp_sync_hint_nonspeculative,
  omp_sync_hint_speculative = 8,
  omp_lock_hint_speculative __GOMP_DEPRECATED_5_0 = omp_sync_hint_speculative
} omp_sync_hint_t;

typedef __GOMP_DEPRECATED_5_0 omp_sync_hint_t omp_lock_hint_t;

typedef struct __attribute__((__aligned__ (sizeof (void *)))) omp_depend_t
{
  char __omp_depend_t__[2 * sizeof (void *)];
} omp_depend_t;

typedef enum omp_pause_resource_t
{
  omp_pause_soft = 1,
  omp_pause_hard = 2
} omp_pause_resource_t;

typedef __UINTPTR_TYPE__ omp_uintptr_t;

#if __cplusplus >= 201103L
# define __GOMP_UINTPTR_T_ENUM : omp_uintptr_t
#else
# define __GOMP_UINTPTR_T_ENUM
#endif

typedef enum omp_memspace_handle_t __GOMP_UINTPTR_T_ENUM
{
  omp_default_mem_space = 0,
  omp_large_cap_mem_space = 1,
  omp_const_mem_space = 2,
  omp_high_bw_mem_space = 3,
  omp_low_lat_mem_space = 4,
  __omp_memspace_handle_t_max__ = __UINTPTR_MAX__
} omp_memspace_handle_t;

typedef enum omp_allocator_handle_t __GOMP_UINTPTR_T_ENUM
{
  omp_null_allocator = 0,
  omp_default_mem_alloc = 1,
  omp_large_cap_mem_alloc = 2,
  omp_const_mem_alloc = 3,
  omp_high_bw_mem_alloc = 4,
  omp_low_lat_mem_alloc = 5,
  omp_cgroup_mem_alloc = 6,
  omp_pteam_mem_alloc = 7,
  omp_thread_mem_alloc = 8,
  __omp_allocator_handle_t_max__ = __UINTPTR_MAX__
} omp_allocator_handle_t;

typedef enum omp_alloctrait_key_t
{
  omp_atk_sync_hint = 1,
  omp_atk_alignment = 2,
  omp_atk_access = 3,
  omp_atk_pool_size = 4,
  omp_atk_fallback = 5,
  omp_atk_fb_data = 6,
  omp_atk_pinned = 7,
  omp_atk_partition = 8
} omp_alloctrait_key_t;

typedef enum omp_alloctrait_value_t
{
  omp_atv_default = (__UINTPTR_TYPE__) -1,
  omp_atv_false = 0,
  omp_atv_true = 1,
  omp_atv_contended = 3,
  omp_atv_uncontended = 4,
  omp_atv_serialized = 5,
  omp_atv_sequential __GOMP_DEPRECATED_5_1 = omp_atv_serialized,
  omp_atv_private = 6,
  omp_atv_all = 7,
  omp_atv_thread = 8,
  omp_atv_pteam = 9,
  omp_atv_cgroup = 10,
  omp_atv_default_mem_fb = 11,
  omp_atv_null_fb = 12,
  omp_atv_abort_fb = 13,
  omp_atv_allocator_fb = 14,
  omp_atv_environment = 15,
  omp_atv_nearest = 16,
  omp_atv_blocked = 17,
  omp_atv_interleaved = 18
} omp_alloctrait_value_t;

typedef struct omp_alloctrait_t
{
  omp_alloctrait_key_t key;
  omp_uintptr_t value;
} omp_alloctrait_t;

typedef enum omp_event_handle_t __GOMP_UINTPTR_T_ENUM
{
  __omp_event_handle_t_max__ = __UINTPTR_MAX__
} omp_event_handle_t;

#ifdef __cplusplus
extern "C" {
# define __GOMP_NOTHROW throw ()
# define __GOMP_DEFAULT_NULL_ALLOCATOR = omp_null_allocator
#else
# define __GOMP_NOTHROW __attribute__((__nothrow__))
# define __GOMP_DEFAULT_NULL_ALLOCATOR
#endif

extern void omp_set_num_threads (int) __GOMP_NOTHROW;
extern int omp_get_num_threads (void) __GOMP_NOTHROW;
extern int omp_get_max_threads (void) __GOMP_NOTHROW;
extern int omp_get_thread_num (void) __GOMP_NOTHROW;
extern int omp_get_num_procs (void) __GOMP_NOTHROW;

extern int omp_in_parallel (void) __GOMP_NOTHROW;

extern void omp_set_dynamic (int) __GOMP_NOTHROW;
extern int omp_get_dynamic (void) __GOMP_NOTHROW;

extern void omp_set_nested (int) __GOMP_NOTHROW __GOMP_DEPRECATED_5_0;
extern int omp_get_nested (void) __GOMP_NOTHROW __GOMP_DEPRECATED_5_0;

extern void omp_init_lock (omp_lock_t *) __GOMP_NOTHROW;
extern void omp_init_lock_with_hint (omp_lock_t *, omp_sync_hint_t)
  __GOMP_NOTHROW;
extern void omp_destroy_lock (omp_lock_t *) __GOMP_NOTHROW;
extern void omp_set_lock (omp_lock_t *) __GOMP_NOTHROW;
extern void omp_unset_lock (omp_lock_t *) __GOMP_NOTHROW;
extern int omp_test_lock (omp_lock_t *) __GOMP_NOTHROW;

extern void omp_init_nest_lock (omp_nest_lock_t *) __GOMP_NOTHROW;
extern void omp_init_nest_lock_with_hint (omp_nest_lock_t *, omp_sync_hint_t)
  __GOMP_NOTHROW;
extern void omp_destroy_nest_lock (omp_nest_lock_t *) __GOMP_NOTHROW;
extern void omp_set_nest_lock (omp_nest_lock_t *) __GOMP_NOTHROW;
extern void omp_unset_nest_lock (omp_nest_lock_t *) __GOMP_NOTHROW;
extern int omp_test_nest_lock (omp_nest_lock_t *) __GOMP_NOTHROW;

extern double omp_get_wtime (void) __GOMP_NOTHROW;
extern double omp_get_wtick (void) __GOMP_NOTHROW;

extern void omp_set_schedule (omp_sched_t, int) __GOMP_NOTHROW;
extern void omp_get_schedule (omp_sched_t *, int *) __GOMP_NOTHROW;
extern int omp_get_thread_limit (void) __GOMP_NOTHROW;
extern void omp_set_max_active_levels (int) __GOMP_NOTHROW;
extern int omp_get_max_active_levels (void) __GOMP_NOTHROW;
extern int omp_get_supported_active_levels (void) __GOMP_NOTHROW;
extern int omp_get_level (void) __GOMP_NOTHROW;
extern int omp_get_ancestor_thread_num (int) __GOMP_NOTHROW;
extern int omp_get_team_size (int) __GOMP_NOTHROW;
extern int omp_get_active_level (void) __GOMP_NOTHROW;

extern int omp_in_final (void) __GOMP_NOTHROW;

extern int omp_get_cancellation (void) __GOMP_NOTHROW;
extern omp_proc_bind_t omp_get_proc_bind (void) __GOMP_NOTHROW;
extern int omp_get_num_places (void) __GOMP_NOTHROW;
extern int omp_get_place_num_procs (int) __GOMP_NOTHROW;
extern void omp_get_place_proc_ids (int, int *) __GOMP_NOTHROW;
extern int omp_get_place_num (void) __GOMP_NOTHROW;
extern int omp_get_partition_num_places (void) __GOMP_NOTHROW;
extern void omp_get_partition_place_nums (int *) __GOMP_NOTHROW;

extern void omp_set_default_device (int) __GOMP_NOTHROW;
extern int omp_get_default_device (void) __GOMP_NOTHROW;
extern int omp_get_num_devices (void) __GOMP_NOTHROW;
extern int omp_get_device_num (void) __GOMP_NOTHROW;
extern int omp_get_num_teams (void) __GOMP_NOTHROW;
extern int omp_get_team_num (void) __GOMP_NOTHROW;

extern int omp_is_initial_device (void) __GOMP_NOTHROW;
extern int omp_get_initial_device (void) __GOMP_NOTHROW;
extern int omp_get_max_task_priority (void) __GOMP_NOTHROW;

extern void omp_fulfill_event (omp_event_handle_t) __GOMP_NOTHROW;

extern void omp_set_num_teams (int) __GOMP_NOTHROW;
extern int omp_get_max_teams (void) __GOMP_NOTHROW;
extern void omp_set_teams_thread_limit (int) __GOMP_NOTHROW;
extern int omp_get_teams_thread_limit (void) __GOMP_NOTHROW;

extern void *omp_target_alloc (__SIZE_TYPE__, int) __GOMP_NOTHROW;
extern void omp_target_free (void *, int) __GOMP_NOTHROW;
extern int omp_target_is_present (const void *, int) __GOMP_NOTHROW;
extern int omp_target_memcpy (void *, const void *, __SIZE_TYPE__,
			      __SIZE_TYPE__, __SIZE_TYPE__, int, int)
  __GOMP_NOTHROW;
extern int omp_target_memcpy_rect (void *, const void *, __SIZE_TYPE__, int,
				   const __SIZE_TYPE__ *,
				   const __SIZE_TYPE__ *,
				   const __SIZE_TYPE__ *,
				   const __SIZE_TYPE__ *,
				   const __SIZE_TYPE__ *, int, int)
  __GOMP_NOTHROW;
extern int omp_target_associate_ptr (const void *, const void *, __SIZE_TYPE__,
				     __SIZE_TYPE__, int) __GOMP_NOTHROW;
extern int omp_target_disassociate_ptr (const void *, int) __GOMP_NOTHROW;

extern void omp_set_affinity_format (const char *) __GOMP_NOTHROW;
extern __SIZE_TYPE__ omp_get_affinity_format (char *, __SIZE_TYPE__)
  __GOMP_NOTHROW;
extern void omp_display_affinity (const char *) __GOMP_NOTHROW;
extern __SIZE_TYPE__ omp_capture_affinity (char *, __SIZE_TYPE__, const char *)
  __GOMP_NOTHROW;

extern int omp_pause_resource (omp_pause_resource_t, int) __GOMP_NOTHROW;
extern int omp_pause_resource_all (omp_pause_resource_t) __GOMP_NOTHROW;

extern omp_allocator_handle_t omp_init_allocator (omp_memspace_handle_t,
						  int,
						  const omp_alloctrait_t [])
  __GOMP_NOTHROW;
extern void omp_destroy_allocator (omp_allocator_handle_t) __GOMP_NOTHROW;
extern void omp_set_default_allocator (omp_allocator_handle_t) __GOMP_NOTHROW;
extern omp_allocator_handle_t omp_get_default_allocator (void) __GOMP_NOTHROW;
extern void omp_free (void *,
		      omp_allocator_handle_t __GOMP_DEFAULT_NULL_ALLOCATOR)
  __GOMP_NOTHROW;
extern void *omp_alloc (__SIZE_TYPE__,
			omp_allocator_handle_t __GOMP_DEFAULT_NULL_ALLOCATOR)
  __GOMP_NOTHROW __attribute__((__malloc__, __malloc__ (omp_free),
				__alloc_size__ (1)));
extern void *omp_aligned_alloc (__SIZE_TYPE__, __SIZE_TYPE__,
				omp_allocator_handle_t
				__GOMP_DEFAULT_NULL_ALLOCATOR)
  __GOMP_NOTHROW __attribute__((__malloc__, __malloc__ (omp_free),
				__alloc_size__ (2), __alloc_align__ (1)));
extern void *omp_calloc (__SIZE_TYPE__, __SIZE_TYPE__,
			 omp_allocator_handle_t __GOMP_DEFAULT_NULL_ALLOCATOR)
  __GOMP_NOTHROW __attribute__((__malloc__, __malloc__ (omp_free),
				__alloc_size__ (1, 2)));
extern void *omp_aligned_calloc (__SIZE_TYPE__, __SIZE_TYPE__, __SIZE_TYPE__,
				 omp_allocator_handle_t
				 __GOMP_DEFAULT_NULL_ALLOCATOR)
  __GOMP_NOTHROW __attribute__((__malloc__, __malloc__ (omp_free),
				__alloc_size__ (2, 3), __alloc_align__ (1)));
extern void *omp_realloc (void *, __SIZE_TYPE__,
			  omp_allocator_handle_t __GOMP_DEFAULT_NULL_ALLOCATOR,
			  omp_allocator_handle_t __GOMP_DEFAULT_NULL_ALLOCATOR)
  __GOMP_NOTHROW __attribute__((__malloc__ (omp_free), __alloc_size__ (2)));

extern void omp_display_env (int) __GOMP_NOTHROW;

#ifdef __cplusplus
}
#endif

#endif /* _OMP_H */
