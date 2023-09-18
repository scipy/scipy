!  Copyright (C) 2005-2022 Free Software Foundation, Inc.
!  Contributed by Jakub Jelinek <jakub@redhat.com>.

!  This file is part of the GNU Offloading and Multi Processing Library
!  (libgomp).

!  Libgomp is free software; you can redistribute it and/or modify it
!  under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3, or (at your option)
!  any later version.

!  Libgomp is distributed in the hope that it will be useful, but WITHOUT ANY
!  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
!  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
!  more details.

!  Under Section 7 of GPL version 3, you are granted additional
!  permissions described in the GCC Runtime Library Exception, version
!  3.1, as published by the Free Software Foundation.

!  You should have received a copy of the GNU General Public License and
!  a copy of the GCC Runtime Library Exception along with this program;
!  see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
!  <http://www.gnu.org/licenses/>.

      integer omp_lock_kind, omp_nest_lock_kind, openmp_version
      parameter (omp_lock_kind = 4)
      parameter (omp_nest_lock_kind = 8)
      integer omp_depend_kind
      parameter (omp_depend_kind = 16)
      integer omp_sched_kind
      parameter (omp_sched_kind = 4)
      integer (omp_sched_kind) omp_sched_static, omp_sched_dynamic
      integer (omp_sched_kind) omp_sched_guided, omp_sched_auto
      parameter (omp_sched_static = 1)
      parameter (omp_sched_dynamic = 2)
      parameter (omp_sched_guided = 3)
      parameter (omp_sched_auto = 4)
      integer omp_proc_bind_kind
      parameter (omp_proc_bind_kind = 4)
      integer (omp_proc_bind_kind) omp_proc_bind_false
      integer (omp_proc_bind_kind) omp_proc_bind_true
      integer (omp_proc_bind_kind) omp_proc_bind_primary
      integer (omp_proc_bind_kind) omp_proc_bind_master
      integer (omp_proc_bind_kind) omp_proc_bind_close
      integer (omp_proc_bind_kind) omp_proc_bind_spread
      parameter (omp_proc_bind_false = 0)
      parameter (omp_proc_bind_true = 1)
      parameter (omp_proc_bind_primary = 2)
      parameter (omp_proc_bind_master = 2)
      parameter (omp_proc_bind_close = 3)
      parameter (omp_proc_bind_spread = 4)
      integer omp_sync_hint_kind
      integer omp_lock_hint_kind
      parameter (omp_sync_hint_kind = 4)
      parameter (omp_lock_hint_kind = omp_sync_hint_kind)
      integer (omp_sync_hint_kind) omp_sync_hint_none
      integer (omp_lock_hint_kind) omp_lock_hint_none
      integer (omp_sync_hint_kind) omp_sync_hint_uncontended
      integer (omp_lock_hint_kind) omp_lock_hint_uncontended
      integer (omp_sync_hint_kind) omp_sync_hint_contended
      integer (omp_sync_hint_kind) omp_lock_hint_contended
      integer (omp_lock_hint_kind) omp_sync_hint_nonspeculative
      integer (omp_lock_hint_kind) omp_lock_hint_nonspeculative
      integer (omp_sync_hint_kind) omp_sync_hint_speculative
      integer (omp_lock_hint_kind) omp_lock_hint_speculative
      parameter (omp_sync_hint_none = 0)
      parameter (omp_lock_hint_none = 0)
      parameter (omp_sync_hint_uncontended = 1)
      parameter (omp_lock_hint_uncontended = 1)
      parameter (omp_sync_hint_contended = 2)
      parameter (omp_lock_hint_contended = 2)
      parameter (omp_sync_hint_nonspeculative = 4)
      parameter (omp_lock_hint_nonspeculative = 4)
      parameter (omp_sync_hint_speculative = 8)
      parameter (omp_lock_hint_speculative = 8)
      parameter (openmp_version = 201511)
      integer omp_pause_resource_kind
      parameter (omp_pause_resource_kind = 4)
      integer (omp_pause_resource_kind) omp_pause_soft
      integer (omp_pause_resource_kind) omp_pause_hard
      parameter (omp_pause_soft = 1)
      parameter (omp_pause_hard = 2)

      integer omp_allocator_handle_kind, omp_alloctrait_key_kind
      integer omp_alloctrait_val_kind, omp_memspace_handle_kind
      integer omp_event_handle_kind
      parameter (omp_allocator_handle_kind = 8)
      parameter (omp_alloctrait_key_kind = 4)
      parameter (omp_alloctrait_val_kind = 8)
      parameter (omp_memspace_handle_kind = 8)
      parameter (omp_event_handle_kind = 8)
      integer (omp_alloctrait_key_kind) omp_atk_sync_hint
      integer (omp_alloctrait_key_kind) omp_atk_alignment
      integer (omp_alloctrait_key_kind) omp_atk_access
      integer (omp_alloctrait_key_kind) omp_atk_pool_size
      integer (omp_alloctrait_key_kind) omp_atk_fallback
      integer (omp_alloctrait_key_kind) omp_atk_fb_data
      integer (omp_alloctrait_key_kind) omp_atk_pinned
      integer (omp_alloctrait_key_kind) omp_atk_partition
      parameter (omp_atk_sync_hint = 1)
      parameter (omp_atk_alignment = 2)
      parameter (omp_atk_access = 3)
      parameter (omp_atk_pool_size = 4)
      parameter (omp_atk_fallback = 5)
      parameter (omp_atk_fb_data = 6)
      parameter (omp_atk_pinned = 7)
      parameter (omp_atk_partition = 8)
      integer (omp_alloctrait_val_kind) omp_atv_false
      integer (omp_alloctrait_val_kind) omp_atv_true
      integer (omp_alloctrait_val_kind) omp_atv_default
      integer (omp_alloctrait_val_kind) omp_atv_contended
      integer (omp_alloctrait_val_kind) omp_atv_uncontended
      integer (omp_alloctrait_val_kind) omp_atv_serialized
      integer (omp_alloctrait_val_kind) omp_atv_sequential
      integer (omp_alloctrait_val_kind) omp_atv_private
      integer (omp_alloctrait_val_kind) omp_atv_all
      integer (omp_alloctrait_val_kind) omp_atv_thread
      integer (omp_alloctrait_val_kind) omp_atv_pteam
      integer (omp_alloctrait_val_kind) omp_atv_cgroup
      integer (omp_alloctrait_val_kind) omp_atv_default_mem_fb
      integer (omp_alloctrait_val_kind) omp_atv_null_fb
      integer (omp_alloctrait_val_kind) omp_atv_abort_fb
      integer (omp_alloctrait_val_kind) omp_atv_allocator_fb
      integer (omp_alloctrait_val_kind) omp_atv_environment
      integer (omp_alloctrait_val_kind) omp_atv_nearest
      integer (omp_alloctrait_val_kind) omp_atv_blocked
      integer (omp_alloctrait_val_kind) omp_atv_interleaved
      parameter (omp_atv_default = -1)
      parameter (omp_atv_false = 0)
      parameter (omp_atv_true = 1)
      parameter (omp_atv_contended = 3)
      parameter (omp_atv_uncontended = 4)
      parameter (omp_atv_serialized = 5)
      parameter (omp_atv_sequential = omp_atv_serialized)
      parameter (omp_atv_private = 6)
      parameter (omp_atv_all = 7)
      parameter (omp_atv_thread = 8)
      parameter (omp_atv_pteam = 9)
      parameter (omp_atv_cgroup = 10)
      parameter (omp_atv_default_mem_fb = 11)
      parameter (omp_atv_null_fb = 12)
      parameter (omp_atv_abort_fb = 13)
      parameter (omp_atv_allocator_fb = 14)
      parameter (omp_atv_environment = 15)
      parameter (omp_atv_nearest = 16)
      parameter (omp_atv_blocked = 17)
      parameter (omp_atv_interleaved = 18)
      integer (omp_allocator_handle_kind) omp_null_allocator
      integer (omp_allocator_handle_kind) omp_default_mem_alloc
      integer (omp_allocator_handle_kind) omp_large_cap_mem_alloc
      integer (omp_allocator_handle_kind) omp_const_mem_alloc
      integer (omp_allocator_handle_kind) omp_high_bw_mem_alloc
      integer (omp_allocator_handle_kind) omp_low_lat_mem_alloc
      integer (omp_allocator_handle_kind) omp_cgroup_mem_alloc
      integer (omp_allocator_handle_kind) omp_pteam_mem_alloc
      integer (omp_allocator_handle_kind) omp_thread_mem_alloc
      parameter (omp_null_allocator = 0)
      parameter (omp_default_mem_alloc = 1)
      parameter (omp_large_cap_mem_alloc = 2)
      parameter (omp_const_mem_alloc = 3)
      parameter (omp_high_bw_mem_alloc = 4)
      parameter (omp_low_lat_mem_alloc = 5)
      parameter (omp_cgroup_mem_alloc = 6)
      parameter (omp_pteam_mem_alloc = 7)
      parameter (omp_thread_mem_alloc = 8)
      integer (omp_memspace_handle_kind) omp_default_mem_space
      integer (omp_memspace_handle_kind) omp_large_cap_mem_space
      integer (omp_memspace_handle_kind) omp_const_mem_space
      integer (omp_memspace_handle_kind) omp_high_bw_mem_space
      integer (omp_memspace_handle_kind) omp_low_lat_mem_space
      parameter (omp_default_mem_space = 0)
      parameter (omp_large_cap_mem_space = 1)
      parameter (omp_const_mem_space = 2)
      parameter (omp_high_bw_mem_space = 3)
      parameter (omp_low_lat_mem_space = 4)

      type omp_alloctrait
        integer (omp_alloctrait_key_kind) key
        integer (omp_alloctrait_val_kind) value
      end type omp_alloctrait

      external omp_init_lock, omp_init_nest_lock
      external omp_init_lock_with_hint
      external omp_init_nest_lock_with_hint
      external omp_destroy_lock, omp_destroy_nest_lock
      external omp_set_lock, omp_set_nest_lock
      external omp_unset_lock, omp_unset_nest_lock
      external omp_set_dynamic, omp_set_nested
      external omp_set_num_threads

      external omp_get_dynamic, omp_get_nested
      logical(4) omp_get_dynamic, omp_get_nested
      external omp_test_lock, omp_in_parallel
      logical(4) omp_test_lock, omp_in_parallel

      external omp_get_max_threads, omp_get_num_procs
      integer(4) omp_get_max_threads, omp_get_num_procs
      external omp_get_num_threads, omp_get_thread_num
      integer(4) omp_get_num_threads, omp_get_thread_num
      external omp_test_nest_lock
      integer(4) omp_test_nest_lock

      external omp_get_wtick, omp_get_wtime
      double precision omp_get_wtick, omp_get_wtime

      external omp_set_schedule, omp_get_schedule
      external omp_get_thread_limit, omp_set_max_active_levels
      external omp_get_max_active_levels, omp_get_level
      external omp_get_ancestor_thread_num, omp_get_team_size
      external omp_get_active_level
      external omp_get_supported_active_levels
      integer(4) omp_get_thread_limit, omp_get_max_active_levels
      integer(4) omp_get_level, omp_get_ancestor_thread_num
      integer(4) omp_get_team_size, omp_get_active_level
      integer(4) omp_get_supported_active_levels

      external omp_in_final
      logical(4) omp_in_final

      external omp_get_cancellation
      logical(4) omp_get_cancellation

      external omp_get_proc_bind
      integer(omp_proc_bind_kind) omp_get_proc_bind

      integer(4) omp_get_num_places
      external omp_get_num_places
      integer(4) omp_get_place_num_procs
      external omp_get_place_num_procs
      external omp_get_place_proc_ids
      integer(4) omp_get_place_num
      external omp_get_place_num
      integer(4) omp_get_partition_num_places
      external omp_get_partition_num_places
      external omp_get_partition_place_nums

      external omp_set_default_device, omp_get_default_device
      external omp_get_num_devices, omp_get_num_teams
      external omp_get_team_num
      integer(4) omp_get_default_device, omp_get_num_devices
      integer(4) omp_get_num_teams, omp_get_team_num

      external omp_is_initial_device
      logical(4) omp_is_initial_device
      external omp_get_initial_device
      integer(4) omp_get_initial_device

      external omp_get_device_num
      integer(4) omp_get_device_num

      external omp_get_max_task_priority
      integer(4) omp_get_max_task_priority

      external omp_set_num_teams, omp_set_teams_thread_limit
      external omp_get_max_teams, omp_get_teams_thread_limit
      integer(4) omp_get_max_teams, omp_get_teams_thread_limit

      external omp_fulfill_event

      external omp_set_affinity_format, omp_get_affinity_format
      external omp_display_affinity, omp_capture_affinity
      integer(4) omp_get_affinity_format
      integer(4) omp_capture_affinity

      external omp_pause_resource, omp_pause_resource_all
      integer(4) omp_pause_resource
      integer(4) omp_pause_resource_all

      external omp_init_allocator
      integer (omp_allocator_handle_kind) omp_init_allocator
      external omp_destroy_allocator
      external omp_set_default_allocator
      external omp_get_default_allocator
      integer (omp_allocator_handle_kind) omp_get_default_allocator

      external omp_display_env

      interface
        function omp_alloc (size, allocator) bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
          use, intrinsic :: omp_lib_kinds
          type(c_ptr) :: omp_alloc
          integer(c_size_t), value :: size
          integer(omp_allocator_handle_kind), value :: allocator
        end function omp_alloc
      end interface

      interface
        function omp_aligned_alloc (alignment, size, allocator) bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
          use, intrinsic :: omp_lib_kinds
          type(c_ptr) :: omp_aligned_alloc
          integer(c_size_t), value :: alignment, size
          integer(omp_allocator_handle_kind), value :: allocator
        end function omp_aligned_alloc
      end interface

      interface
        subroutine omp_free(ptr, allocator) bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr
          use, intrinsic :: omp_lib_kinds
          type(c_ptr), value :: ptr
          integer(omp_allocator_handle_kind), value :: allocator
        end subroutine omp_free
      end interface

      interface
        function omp_calloc (nmemb, size, allocator) bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
          use, intrinsic :: omp_lib_kinds
          type(c_ptr) :: omp_calloc
          integer(c_size_t), value :: nmemb, size
          integer(omp_allocator_handle_kind), value :: allocator
        end function omp_calloc
      end interface

      interface
        function omp_aligned_calloc (alignment, nmemb, size, allocator)   &
     &      bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
          use, intrinsic :: omp_lib_kinds
          type(c_ptr) :: omp_aligned_calloc
          integer(c_size_t), value :: alignment, nmemb, size
          integer(omp_allocator_handle_kind), value :: allocator
        end function omp_aligned_calloc
      end interface

      interface
        function omp_realloc (ptr, size, allocator, free_allocator)      &
     &      bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
          use, intrinsic :: omp_lib_kinds
          type(c_ptr) :: omp_realloc
          type(c_ptr), value :: ptr
          integer(c_size_t), value :: size
          integer(omp_allocator_handle_kind), value :: allocator
          integer(omp_allocator_handle_kind), value :: free_allocator
        end function omp_realloc
      end interface

      interface
        function omp_target_alloc (size, device_num) bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t, c_int
          type(c_ptr) :: omp_target_alloc
          integer(c_size_t), value :: size
          integer(c_int), value :: device_num
        end function omp_target_alloc
      end interface

      interface
        subroutine omp_target_free (device_ptr, device_num) bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_int
          type(c_ptr), value :: device_ptr
          integer(c_int), value :: device_num
        end subroutine omp_target_free
      end interface

      interface
        function omp_target_is_present (ptr, device_num) bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_int
          integer(c_int) :: omp_target_is_present
          type(c_ptr), value :: ptr
          integer(c_int), value :: device_num
        end function omp_target_is_present
      end interface

      interface
        function omp_target_memcpy (dst, src, length, dst_offset,          &
     &                              src_offset, dst_device_num,            &
     &                              src_device_num) bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_size_t
          integer(c_int) :: omp_target_memcpy
          type(c_ptr), value :: dst, src
          integer(c_size_t), value :: length, dst_offset, src_offset
          integer(c_int), value :: dst_device_num, src_device_num
        end function omp_target_memcpy
      end interface

      interface
        function omp_target_memcpy_rect (dst,src,element_size, num_dims,   &
     &                                   volume, dst_offsets,              &
     &                                   src_offsets, dst_dimensions,      &
     &                                   src_dimensions, dst_device_num,   &
     &                                   src_device_num) bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_size_t
          integer(c_int) :: omp_target_memcpy_rect
          type(c_ptr), value :: dst, src
          integer(c_size_t), value :: element_size
          integer(c_int), value :: num_dims
          integer(c_int), value :: dst_device_num, src_device_num
          integer(c_size_t), intent(in) :: volume(*), dst_offsets(*)
          integer(c_size_t), intent(in) :: src_offsets(*)
          integer(c_size_t), intent(in) :: dst_dimensions(*)
          integer(c_size_t), intent(in) :: src_dimensions(*)
        end function omp_target_memcpy_rect
      end interface

      interface
        function omp_target_associate_ptr (host_ptr, device_ptr, size,     &
     &                                     device_offset, device_num)      &
     &      bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t, c_int
          integer(c_int) :: omp_target_associate_ptr
          type(c_ptr), value :: host_ptr, device_ptr
          integer(c_size_t), value :: size, device_offset
          integer(c_int), value :: device_num
        end function omp_target_associate_ptr
      end interface

      interface
        function omp_target_disassociate_ptr (ptr, device_num) bind(c)
          use, intrinsic :: iso_c_binding, only : c_ptr, c_int
          integer(c_int) :: omp_target_disassociate_ptr
          type(c_ptr), value :: ptr
          integer(c_int), value :: device_num
        end function omp_target_disassociate_ptr
      end interface
