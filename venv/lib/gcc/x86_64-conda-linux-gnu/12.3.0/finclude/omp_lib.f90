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

      module omp_lib_kinds
        use iso_c_binding, only: c_int, c_intptr_t
        implicit none
        private :: c_int, c_intptr_t
        integer, parameter :: omp_lock_kind = 4
        integer, parameter :: omp_nest_lock_kind = 8
        integer, parameter :: omp_sched_kind = 4
        integer, parameter :: omp_proc_bind_kind = 4
        integer, parameter :: omp_sync_hint_kind = 4
        integer, parameter :: omp_lock_hint_kind = omp_sync_hint_kind
        integer, parameter :: omp_pause_resource_kind = 4
        integer, parameter :: omp_allocator_handle_kind = c_intptr_t
        integer, parameter :: omp_alloctrait_key_kind = c_int
        integer, parameter :: omp_alloctrait_val_kind = c_intptr_t
        integer, parameter :: omp_memspace_handle_kind = c_intptr_t
        integer, parameter :: omp_depend_kind = 16
        integer, parameter :: omp_event_handle_kind = c_intptr_t
        integer (omp_sched_kind), parameter :: omp_sched_static = 1
        integer (omp_sched_kind), parameter :: omp_sched_dynamic = 2
        integer (omp_sched_kind), parameter :: omp_sched_guided = 3
        integer (omp_sched_kind), parameter :: omp_sched_auto = 4
        integer (omp_proc_bind_kind), &
                 parameter :: omp_proc_bind_false = 0
        integer (omp_proc_bind_kind), &
                 parameter :: omp_proc_bind_true = 1
        integer (omp_proc_bind_kind), &
                 parameter :: omp_proc_bind_primary = 2
        integer (omp_proc_bind_kind), &
                 parameter :: omp_proc_bind_master = 2
        integer (omp_proc_bind_kind), &
                 parameter :: omp_proc_bind_close = 3
        integer (omp_proc_bind_kind), &
                 parameter :: omp_proc_bind_spread = 4
        integer (omp_lock_hint_kind), &
                 parameter :: omp_sync_hint_none = 0
        integer (omp_lock_hint_kind), &
                 parameter :: omp_lock_hint_none = omp_sync_hint_none
        integer (omp_lock_hint_kind), &
                 parameter :: omp_sync_hint_uncontended = 1
        integer (omp_lock_hint_kind), &
                 parameter :: omp_lock_hint_uncontended &
                 = omp_sync_hint_uncontended
        integer (omp_lock_hint_kind), &
                 parameter :: omp_sync_hint_contended = 2
        integer (omp_lock_hint_kind), &
                 parameter :: omp_lock_hint_contended &
                 = omp_sync_hint_contended
        integer (omp_lock_hint_kind), &
                 parameter :: omp_sync_hint_nonspeculative = 4
        integer (omp_lock_hint_kind), &
                 parameter :: omp_lock_hint_nonspeculative &
                 = omp_sync_hint_nonspeculative
        integer (omp_lock_hint_kind), &
                 parameter :: omp_sync_hint_speculative = 8
        integer (omp_lock_hint_kind), &
                 parameter :: omp_lock_hint_speculative &
                 = omp_sync_hint_speculative
        integer (kind=omp_pause_resource_kind), &
                 parameter :: omp_pause_soft = 1
        integer (kind=omp_pause_resource_kind), &
                 parameter :: omp_pause_hard = 2
        integer (kind=omp_alloctrait_key_kind), &
                 parameter :: omp_atk_sync_hint = 1
        integer (kind=omp_alloctrait_key_kind), &
                 parameter :: omp_atk_alignment = 2
        integer (kind=omp_alloctrait_key_kind), &
                 parameter :: omp_atk_access = 3
        integer (kind=omp_alloctrait_key_kind), &
                 parameter :: omp_atk_pool_size = 4
        integer (kind=omp_alloctrait_key_kind), &
                 parameter :: omp_atk_fallback = 5
        integer (kind=omp_alloctrait_key_kind), &
                 parameter :: omp_atk_fb_data = 6
        integer (kind=omp_alloctrait_key_kind), &
                 parameter :: omp_atk_pinned = 7
        integer (kind=omp_alloctrait_key_kind), &
                 parameter :: omp_atk_partition = 8
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_default = -1
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_false = 0
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_true = 1
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_contended = 3
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_uncontended = 4
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_serialized = 5
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_sequential = omp_atv_serialized
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_private = 6
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_all = 7
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_thread = 8
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_pteam = 9
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_cgroup = 10
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_default_mem_fb = 11
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_null_fb = 12
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_abort_fb = 13
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_allocator_fb = 14
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_environment = 15
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_nearest = 16
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_blocked = 17
        integer (kind=omp_alloctrait_val_kind), &
                 parameter :: omp_atv_interleaved = 18
        integer (kind=omp_allocator_handle_kind), &
                 parameter :: omp_null_allocator = 0
        integer (kind=omp_allocator_handle_kind), &
                 parameter :: omp_default_mem_alloc = 1
        integer (kind=omp_allocator_handle_kind), &
                 parameter :: omp_large_cap_mem_alloc = 2
        integer (kind=omp_allocator_handle_kind), &
                 parameter :: omp_const_mem_alloc = 3
        integer (kind=omp_allocator_handle_kind), &
                 parameter :: omp_high_bw_mem_alloc = 4
        integer (kind=omp_allocator_handle_kind), &
                 parameter :: omp_low_lat_mem_alloc = 5
        integer (kind=omp_allocator_handle_kind), &
                 parameter :: omp_cgroup_mem_alloc = 6
        integer (kind=omp_allocator_handle_kind), &
                 parameter :: omp_pteam_mem_alloc = 7
        integer (kind=omp_allocator_handle_kind), &
                 parameter :: omp_thread_mem_alloc = 8
        integer (omp_memspace_handle_kind), &
                 parameter :: omp_default_mem_space = 0
        integer (omp_memspace_handle_kind), &
                 parameter :: omp_large_cap_mem_space = 1
        integer (omp_memspace_handle_kind), &
                 parameter :: omp_const_mem_space = 2
        integer (omp_memspace_handle_kind), &
                 parameter :: omp_high_bw_mem_space = 3
        integer (omp_memspace_handle_kind), &
                 parameter :: omp_low_lat_mem_space = 4

        type omp_alloctrait
          integer (kind=omp_alloctrait_key_kind) key
          integer (kind=omp_alloctrait_val_kind) value
        end type omp_alloctrait
      end module

      module omp_lib
        use omp_lib_kinds
        implicit none
        integer, parameter :: openmp_version = 201511

        interface
          subroutine omp_init_lock (svar)
            use omp_lib_kinds
            integer (omp_lock_kind), intent (out) :: svar
          end subroutine omp_init_lock
        end interface

        interface
          subroutine omp_init_lock_with_hint (svar, hint)
            use omp_lib_kinds
            integer (omp_lock_kind), intent (out) :: svar
            integer (omp_lock_hint_kind), intent (in) :: hint
          end subroutine omp_init_lock_with_hint
        end interface

        interface
          subroutine omp_init_nest_lock (nvar)
            use omp_lib_kinds
            integer (omp_nest_lock_kind), intent (out) :: nvar
          end subroutine omp_init_nest_lock
        end interface

        interface
          subroutine omp_init_nest_lock_with_hint (nvar, hint)
            use omp_lib_kinds
            integer (omp_nest_lock_kind), intent (out) :: nvar
            integer (omp_lock_hint_kind), intent (in) :: hint
          end subroutine omp_init_nest_lock_with_hint
        end interface

        interface
          subroutine omp_destroy_lock (svar)
            use omp_lib_kinds
            integer (omp_lock_kind), intent (inout) :: svar
          end subroutine omp_destroy_lock
        end interface

        interface
          subroutine omp_destroy_nest_lock (nvar)
            use omp_lib_kinds
            integer (omp_nest_lock_kind), intent (inout) :: nvar
          end subroutine omp_destroy_nest_lock
        end interface

        interface
          subroutine omp_set_lock (svar)
            use omp_lib_kinds
            integer (omp_lock_kind), intent (inout) :: svar
          end subroutine omp_set_lock
        end interface

        interface
          subroutine omp_set_nest_lock (nvar)
            use omp_lib_kinds
            integer (omp_nest_lock_kind), intent (inout) :: nvar
          end subroutine omp_set_nest_lock
        end interface

        interface
          subroutine omp_unset_lock (svar)
            use omp_lib_kinds
            integer (omp_lock_kind), intent (inout) :: svar
          end subroutine omp_unset_lock
        end interface

        interface
          subroutine omp_unset_nest_lock (nvar)
            use omp_lib_kinds
            integer (omp_nest_lock_kind), intent (inout) :: nvar
          end subroutine omp_unset_nest_lock
        end interface

        interface omp_set_dynamic
          subroutine omp_set_dynamic (dynamic_threads)
            logical (4), intent (in) :: dynamic_threads
          end subroutine omp_set_dynamic
          subroutine omp_set_dynamic_8 (dynamic_threads)
            logical (8), intent (in) :: dynamic_threads
          end subroutine omp_set_dynamic_8
        end interface

        interface omp_set_nested
          subroutine omp_set_nested (nested)
            logical (4), intent (in) :: nested
          end subroutine omp_set_nested
          subroutine omp_set_nested_8 (nested)
            logical (8), intent (in) :: nested
          end subroutine omp_set_nested_8
        end interface

        interface omp_set_num_threads
          subroutine omp_set_num_threads (num_threads)
            integer (4), intent (in) :: num_threads
          end subroutine omp_set_num_threads
          subroutine omp_set_num_threads_8 (num_threads)
            integer (8), intent (in) :: num_threads
          end subroutine omp_set_num_threads_8
        end interface

        interface
          function omp_get_dynamic ()
            logical (4) :: omp_get_dynamic
          end function omp_get_dynamic
        end interface

        interface
          function omp_get_nested ()
            logical (4) :: omp_get_nested
          end function omp_get_nested
        end interface

        interface
          function omp_in_parallel ()
            logical (4) :: omp_in_parallel
          end function omp_in_parallel
        end interface

        interface
          function omp_test_lock (svar)
            use omp_lib_kinds
            logical (4) :: omp_test_lock
            integer (omp_lock_kind), intent (inout) :: svar
          end function omp_test_lock
        end interface

        interface
          function omp_get_max_threads ()
            integer (4) :: omp_get_max_threads
          end function omp_get_max_threads
        end interface

        interface
          function omp_get_num_procs ()
            integer (4) :: omp_get_num_procs
          end function omp_get_num_procs
        end interface

        interface
          function omp_get_num_threads ()
            integer (4) :: omp_get_num_threads
          end function omp_get_num_threads
        end interface

        interface
          function omp_get_thread_num ()
            integer (4) :: omp_get_thread_num
          end function omp_get_thread_num
        end interface

        interface
          function omp_test_nest_lock (nvar)
            use omp_lib_kinds
            integer (4) :: omp_test_nest_lock
            integer (omp_nest_lock_kind), intent (inout) :: nvar
          end function omp_test_nest_lock
        end interface

        interface
          function omp_get_wtick ()
            double precision :: omp_get_wtick
          end function omp_get_wtick
        end interface

        interface
          function omp_get_wtime ()
            double precision :: omp_get_wtime
          end function omp_get_wtime
        end interface

        interface omp_set_schedule
          subroutine omp_set_schedule (kind, chunk_size)
            use omp_lib_kinds
            integer (omp_sched_kind), intent (in) :: kind
            integer (4), intent (in) :: chunk_size
          end subroutine omp_set_schedule
          subroutine omp_set_schedule_8 (kind, chunk_size)
            use omp_lib_kinds
            integer (omp_sched_kind), intent (in) :: kind
            integer (8), intent (in) :: chunk_size
          end subroutine omp_set_schedule_8
         end interface

        interface omp_get_schedule
          subroutine omp_get_schedule (kind, chunk_size)
            use omp_lib_kinds
            integer (omp_sched_kind), intent (out) :: kind
            integer (4), intent (out) :: chunk_size
          end subroutine omp_get_schedule
          subroutine omp_get_schedule_8 (kind, chunk_size)
            use omp_lib_kinds
            integer (omp_sched_kind), intent (out) :: kind
            integer (8), intent (out) :: chunk_size
          end subroutine omp_get_schedule_8
         end interface

        interface
          function omp_get_thread_limit ()
            integer (4) :: omp_get_thread_limit
          end function omp_get_thread_limit
        end interface

        interface omp_set_max_active_levels
          subroutine omp_set_max_active_levels (max_levels)
            integer (4), intent (in) :: max_levels
          end subroutine omp_set_max_active_levels
          subroutine omp_set_max_active_levels_8 (max_levels)
            integer (8), intent (in) :: max_levels
          end subroutine omp_set_max_active_levels_8
        end interface

        interface
          function omp_get_max_active_levels ()
            integer (4) :: omp_get_max_active_levels
          end function omp_get_max_active_levels
        end interface

        interface
          function omp_get_supported_active_levels ()
            integer (4) :: omp_get_supported_active_levels
          end function omp_get_supported_active_levels
        end interface

        interface
          function omp_get_level ()
            integer (4) :: omp_get_level
          end function omp_get_level
        end interface

        interface omp_get_ancestor_thread_num
          function omp_get_ancestor_thread_num (level)
            integer (4), intent (in) :: level
            integer (4) :: omp_get_ancestor_thread_num
          end function omp_get_ancestor_thread_num
          function omp_get_ancestor_thread_num_8 (level)
            integer (8), intent (in) :: level
            integer (4) :: omp_get_ancestor_thread_num_8
          end function omp_get_ancestor_thread_num_8
        end interface

        interface omp_get_team_size
          function omp_get_team_size (level)
            integer (4), intent (in) :: level
            integer (4) :: omp_get_team_size
          end function omp_get_team_size
          function omp_get_team_size_8 (level)
            integer (8), intent (in) :: level
            integer (4) :: omp_get_team_size_8
          end function omp_get_team_size_8
        end interface

        interface
          function omp_get_active_level ()
            integer (4) :: omp_get_active_level
          end function omp_get_active_level
        end interface

        interface
          function omp_in_final ()
            logical (4) :: omp_in_final
          end function omp_in_final
        end interface

        interface
          function omp_get_cancellation ()
            logical (4) :: omp_get_cancellation
          end function omp_get_cancellation
        end interface

        interface
          function omp_get_proc_bind ()
            use omp_lib_kinds
            integer (omp_proc_bind_kind) :: omp_get_proc_bind
          end function omp_get_proc_bind
        end interface

        interface
          function omp_get_num_places ()
            integer (4) :: omp_get_num_places
          end function omp_get_num_places
        end interface

        interface omp_get_place_num_procs
          function omp_get_place_num_procs (place_num)
            integer (4), intent(in) :: place_num
            integer (4) :: omp_get_place_num_procs
          end function omp_get_place_num_procs

          function omp_get_place_num_procs_8 (place_num)
            integer (8), intent(in) :: place_num
            integer (4) :: omp_get_place_num_procs_8
          end function omp_get_place_num_procs_8
        end interface

        interface omp_get_place_proc_ids
          subroutine omp_get_place_proc_ids (place_num, ids)
            integer (4), intent(in) :: place_num
            integer (4), intent(out) :: ids(*)
          end subroutine omp_get_place_proc_ids

          subroutine omp_get_place_proc_ids_8 (place_num, ids)
            integer (8), intent(in) :: place_num
            integer (8), intent(out) :: ids(*)
          end subroutine omp_get_place_proc_ids_8
        end interface

        interface
          function omp_get_place_num ()
            integer (4) :: omp_get_place_num
          end function omp_get_place_num
        end interface

        interface
          function omp_get_partition_num_places ()
            integer (4) :: omp_get_partition_num_places
          end function omp_get_partition_num_places
        end interface

        interface omp_get_partition_place_nums
          subroutine omp_get_partition_place_nums (place_nums)
            integer (4), intent(out) :: place_nums(*)
          end subroutine omp_get_partition_place_nums

          subroutine omp_get_partition_place_nums_8 (place_nums)
            integer (8), intent(out) :: place_nums(*)
          end subroutine omp_get_partition_place_nums_8
        end interface

        interface omp_set_default_device
          subroutine omp_set_default_device (device_num)
            integer (4), intent (in) :: device_num
          end subroutine omp_set_default_device
          subroutine omp_set_default_device_8 (device_num)
            integer (8), intent (in) :: device_num
          end subroutine omp_set_default_device_8
        end interface

        interface
          function omp_get_default_device ()
            integer (4) :: omp_get_default_device
          end function omp_get_default_device
        end interface

        interface
          function omp_get_num_devices ()
            integer (4) :: omp_get_num_devices
          end function omp_get_num_devices
        end interface

        interface
          function omp_get_num_teams ()
            integer (4) :: omp_get_num_teams
          end function omp_get_num_teams
        end interface

        interface
          function omp_get_team_num ()
            integer (4) :: omp_get_team_num
          end function omp_get_team_num
        end interface

        interface
          function omp_is_initial_device ()
            logical (4) :: omp_is_initial_device
          end function omp_is_initial_device
        end interface

        interface
          function omp_get_initial_device ()
            integer (4) :: omp_get_initial_device
          end function omp_get_initial_device
        end interface

        interface
          function omp_get_device_num ()
            integer (4) :: omp_get_device_num
          end function omp_get_device_num
        end interface

        interface
          function omp_get_max_task_priority ()
            integer (4) :: omp_get_max_task_priority
          end function omp_get_max_task_priority
        end interface

        interface omp_set_num_teams
          subroutine omp_set_num_teams (num_teams)
            integer (4), intent (in) :: num_teams
          end subroutine omp_set_num_teams
          subroutine omp_set_num_teams_8 (num_teams)
            integer (8), intent (in) :: num_teams
          end subroutine omp_set_num_teams_8
        end interface

        interface
          function omp_get_max_teams ()
            integer (4) :: omp_get_max_teams
          end function omp_get_max_teams
        end interface

        interface omp_set_teams_thread_limit
          subroutine omp_set_teams_thread_limit (thread_limit)
            integer (4), intent (in) :: thread_limit
          end subroutine omp_set_teams_thread_limit
          subroutine omp_set_teams_thread_limit_8 (thread_limit)
            integer (8), intent (in) :: thread_limit
          end subroutine omp_set_teams_thread_limit_8
        end interface

        interface
          function omp_get_teams_thread_limit ()
            integer (4) :: omp_get_teams_thread_limit
          end function omp_get_teams_thread_limit
        end interface

        interface
          subroutine omp_fulfill_event (event)
            use omp_lib_kinds
            integer (kind=omp_event_handle_kind), &
              value, intent(in) :: event
          end subroutine omp_fulfill_event
        end interface

        interface
          subroutine omp_set_affinity_format (format)
            character(len=*), intent(in) :: format
          end subroutine omp_set_affinity_format
        end interface

        interface
          function omp_get_affinity_format (buffer)
            integer (4) :: omp_get_affinity_format
            character(len=*), intent(out) :: buffer
          end function omp_get_affinity_format
        end interface

        interface
          subroutine omp_display_affinity (format)
            character(len=*), intent(in) :: format
          end subroutine omp_display_affinity
        end interface

        interface
          function omp_capture_affinity (buffer, format)
            integer (4) :: omp_capture_affinity
            character(len=*), intent(out) :: buffer
            character(len=*), intent(in) :: format
          end function omp_capture_affinity
        end interface

        interface
          function omp_pause_resource (kind, device_num)
            use omp_lib_kinds
            integer (4) :: omp_pause_resource
            integer (kind=omp_pause_resource_kind), &
              intent(in) :: kind
            integer (4) :: device_num
          end function
        end interface

        interface
          function omp_pause_resource_all (kind)
            use omp_lib_kinds
            integer (4) :: omp_pause_resource_all
            integer (kind=omp_pause_resource_kind), &
              intent(in) :: kind
          end function
        end interface

        interface omp_init_allocator
          function omp_init_allocator (memspace, ntraits, traits)
            use omp_lib_kinds
            integer (kind=omp_allocator_handle_kind) omp_init_allocator
            integer (kind=omp_memspace_handle_kind), &
              intent(in) :: memspace
            integer (4), intent(in) :: ntraits
            type (omp_alloctrait), intent(in) :: traits(*)
          end function
          function omp_init_allocator_8 (memspace, ntraits, traits)
            use omp_lib_kinds
            integer (kind=omp_allocator_handle_kind) omp_init_allocator_8
            integer (kind=omp_memspace_handle_kind), &
              intent(in) :: memspace
            integer (8), intent(in) :: ntraits
            type (omp_alloctrait), intent(in) :: traits(*)
          end function
        end interface

        interface
          subroutine omp_destroy_allocator (allocator)
            use omp_lib_kinds
            integer (kind=omp_allocator_handle_kind), &
              intent(in) :: allocator
          end subroutine
        end interface

        interface
          subroutine omp_set_default_allocator (allocator)
            use omp_lib_kinds
            integer (kind=omp_allocator_handle_kind), &
              intent(in) :: allocator
          end subroutine
        end interface

        interface
          function omp_get_default_allocator ()
            use omp_lib_kinds
            integer (kind=omp_allocator_handle_kind) &
               omp_get_default_allocator
          end function
        end interface

        interface omp_display_env
          subroutine omp_display_env (verbose)
            logical (4),intent (in) :: verbose
          end subroutine omp_display_env
          subroutine omp_display_env_8 (verbose)
            logical (8),intent (in) :: verbose
          end subroutine omp_display_env_8
        end interface

        interface
          function omp_alloc (size, allocator) bind(c)
            use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
            import :: omp_allocator_handle_kind
            type(c_ptr) :: omp_alloc
            integer(c_size_t), value :: size
            integer(omp_allocator_handle_kind), value :: allocator
          end function omp_alloc
        end interface

        interface
          function omp_aligned_alloc (alignment, size, allocator) bind(c)
            use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
            import :: omp_allocator_handle_kind
            type(c_ptr) :: omp_aligned_alloc
            integer(c_size_t), value :: alignment, size
            integer(omp_allocator_handle_kind), value :: allocator
          end function omp_aligned_alloc
        end interface

        interface
          subroutine omp_free(ptr, allocator) bind(c)
            use, intrinsic :: iso_c_binding, only : c_ptr
            import :: omp_allocator_handle_kind
            type(c_ptr), value :: ptr
            integer(omp_allocator_handle_kind), value :: allocator
          end subroutine omp_free
        end interface

        interface
          function omp_calloc (nmemb, size, allocator) bind(c)
            use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
            import :: omp_allocator_handle_kind
            type(c_ptr) :: omp_calloc
            integer(c_size_t), value :: nmemb, size
            integer(omp_allocator_handle_kind), value :: allocator
          end function omp_calloc
        end interface

        interface
          function omp_aligned_calloc (alignment, nmemb, size, allocator) bind(c)
            use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
            import :: omp_allocator_handle_kind
            type(c_ptr) :: omp_aligned_calloc
            integer(c_size_t), value :: alignment, nmemb, size
            integer(omp_allocator_handle_kind), value :: allocator
          end function omp_aligned_calloc
        end interface

        interface
          function omp_realloc (ptr, size, allocator, free_allocator) bind(c)
            use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
            import :: omp_allocator_handle_kind
            type(c_ptr) :: omp_realloc
            type(c_ptr), value :: ptr
            integer(c_size_t), value :: size
            integer(omp_allocator_handle_kind), value :: allocator, free_allocator
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
          function omp_target_memcpy (dst, src, length, dst_offset, &
                                      src_offset, dst_device_num, &
                                      src_device_num) bind(c)
            use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_size_t
            integer(c_int) :: omp_target_memcpy
            type(c_ptr), value :: dst, src
            integer(c_size_t), value :: length, dst_offset, src_offset
            integer(c_int), value :: dst_device_num, src_device_num
          end function omp_target_memcpy
        end interface

        interface
          function omp_target_memcpy_rect (dst,src,element_size, num_dims, &
                                           volume, dst_offsets, src_offsets, &
                                           dst_dimensions, src_dimensions, &
                                           dst_device_num, src_device_num) &
              bind(c)
            use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_size_t
            integer(c_int) :: omp_target_memcpy_rect
            type(c_ptr), value :: dst, src
            integer(c_size_t), value :: element_size
            integer(c_int), value :: num_dims, dst_device_num, src_device_num
            integer(c_size_t), intent(in) :: volume(*), dst_offsets(*),  &
                                             src_offsets(*), dst_dimensions(*), &
                                             src_dimensions(*)
          end function omp_target_memcpy_rect
        end interface

        interface
          function omp_target_associate_ptr (host_ptr, device_ptr, size, &
                                             device_offset, device_num) bind(c)
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

#if _OPENMP >= 201811
!GCC$ ATTRIBUTES DEPRECATED :: omp_get_nested, omp_set_nested
#endif

#if _OPENMP >= 202011
!GCC$ ATTRIBUTES DEPRECATED :: omp_proc_bind_master, omp_atv_sequential
#endif

      end module omp_lib
