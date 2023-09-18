!  OpenACC Runtime Library Definitions.

!  Copyright (C) 2014-2022 Free Software Foundation, Inc.

!  Contributed by Tobias Burnus <burnus@net-b.de>
!              and Mentor Embedded.

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

! Keep in sync with config/accel/openacc.f90 and openacc_lib.h.

module openacc_kinds
  use iso_fortran_env, only: int32
  implicit none

  public
  private :: int32

  ! When adding items, also update 'public' setting in 'module openacc' below.

  integer, parameter :: acc_device_kind = int32

  ! Keep in sync with include/gomp-constants.h.
  integer (acc_device_kind), parameter :: acc_device_current = -1
  integer (acc_device_kind), parameter :: acc_device_none = 0
  integer (acc_device_kind), parameter :: acc_device_default = 1
  integer (acc_device_kind), parameter :: acc_device_host = 2
  ! integer (acc_device_kind), parameter :: acc_device_host_nonshm = 3 removed.
  integer (acc_device_kind), parameter :: acc_device_not_host = 4
  integer (acc_device_kind), parameter :: acc_device_nvidia = 5
  integer (acc_device_kind), parameter :: acc_device_radeon = 8

  integer, parameter :: acc_device_property_kind = int32
  ! OpenACC 2.6/2.7/3.0 used acc_device_property; in a spec update the
  ! missing '_kind' was added for consistency.  For backward compatibility, keep:
  integer, parameter :: acc_device_property = acc_device_property_kind

  ! Keep in sync with 'libgomp/libgomp-plugin.h:goacc_property'.
  integer (acc_device_property_kind), parameter :: acc_property_memory = 1
  integer (acc_device_property_kind), parameter :: acc_property_free_memory = 2
  integer (acc_device_property_kind), parameter :: acc_property_name = int(Z'10001')
  integer (acc_device_property_kind), parameter :: acc_property_vendor = int(Z'10002')
  integer (acc_device_property_kind), parameter :: acc_property_driver = int(Z'10003')

  integer, parameter :: acc_handle_kind = int32

  ! Keep in sync with include/gomp-constants.h.
  integer (acc_handle_kind), parameter :: acc_async_noval = -1
  integer (acc_handle_kind), parameter :: acc_async_sync = -2
end module openacc_kinds

module openacc_internal
  use openacc_kinds
  implicit none

  interface
    function acc_get_num_devices_h (devicetype)
      import
      integer acc_get_num_devices_h
      integer (acc_device_kind) devicetype
    end function

    subroutine acc_set_device_type_h (devicetype)
      import
      integer (acc_device_kind) devicetype
    end subroutine

    function acc_get_device_type_h ()
      import
      integer (acc_device_kind) acc_get_device_type_h
    end function

    subroutine acc_set_device_num_h (devicenum, devicetype)
      import
      integer devicenum
      integer (acc_device_kind) devicetype
    end subroutine

    function acc_get_device_num_h (devicetype)
      import
      integer acc_get_device_num_h
      integer (acc_device_kind) devicetype
    end function

    function acc_get_property_h (devicenum, devicetype, property)
      use iso_c_binding, only: c_size_t
      import
      implicit none (type, external)
      integer (c_size_t) :: acc_get_property_h
      integer, value :: devicenum
      integer (acc_device_kind), value :: devicetype
      integer (acc_device_property_kind), value :: property
    end function

    subroutine acc_get_property_string_h (devicenum, devicetype, property, string)
      import
      implicit none (type, external)
      integer, value :: devicenum
      integer (acc_device_kind), value :: devicetype
      integer (acc_device_property_kind), value :: property
      character (*) :: string
    end subroutine

    function acc_async_test_h (arg)
      logical acc_async_test_h
      integer arg
    end function

    function acc_async_test_all_h ()
      logical acc_async_test_all_h
    end function

    subroutine acc_wait_h (arg)
      integer arg
    end subroutine

    subroutine acc_wait_async_h (arg, async)
      integer arg, async
    end subroutine

    subroutine acc_wait_all_h ()
    end subroutine

    subroutine acc_wait_all_async_h (async)
      integer async
    end subroutine

    subroutine acc_init_h (devicetype)
      import
      integer (acc_device_kind) devicetype
    end subroutine

    subroutine acc_shutdown_h (devicetype)
      import
      integer (acc_device_kind) devicetype
    end subroutine

    function acc_on_device_h (devicetype)
      import
      integer (acc_device_kind) devicetype
      logical acc_on_device_h
    end function

    subroutine acc_copyin_32_h (a, len)
      use iso_c_binding, only: c_int32_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
    end subroutine

    subroutine acc_copyin_64_h (a, len)
      use iso_c_binding, only: c_int64_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
    end subroutine

    subroutine acc_copyin_array_h (a)
      type (*), dimension (..), contiguous :: a
    end subroutine

    subroutine acc_present_or_copyin_32_h (a, len)
      use iso_c_binding, only: c_int32_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
    end subroutine

    subroutine acc_present_or_copyin_64_h (a, len)
      use iso_c_binding, only: c_int64_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
    end subroutine

    subroutine acc_present_or_copyin_array_h (a)
      type (*), dimension (..), contiguous :: a
    end subroutine

    subroutine acc_create_32_h (a, len)
      use iso_c_binding, only: c_int32_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
    end subroutine

    subroutine acc_create_64_h (a, len)
      use iso_c_binding, only: c_int64_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
    end subroutine

    subroutine acc_create_array_h (a)
      type (*), dimension (..), contiguous :: a
    end subroutine

    subroutine acc_present_or_create_32_h (a, len)
      use iso_c_binding, only: c_int32_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
    end subroutine

    subroutine acc_present_or_create_64_h (a, len)
      use iso_c_binding, only: c_int64_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
    end subroutine

    subroutine acc_present_or_create_array_h (a)
      type (*), dimension (..), contiguous :: a
    end subroutine

    subroutine acc_copyout_32_h (a, len)
      use iso_c_binding, only: c_int32_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
    end subroutine

    subroutine acc_copyout_64_h (a, len)
      use iso_c_binding, only: c_int64_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
    end subroutine

    subroutine acc_copyout_array_h (a)
      type (*), dimension (..), contiguous :: a
    end subroutine

    subroutine acc_copyout_finalize_32_h (a, len)
      use iso_c_binding, only: c_int32_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
    end subroutine

    subroutine acc_copyout_finalize_64_h (a, len)
      use iso_c_binding, only: c_int64_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
    end subroutine

    subroutine acc_copyout_finalize_array_h (a)
      type (*), dimension (..), contiguous :: a
    end subroutine

    subroutine acc_delete_32_h (a, len)
      use iso_c_binding, only: c_int32_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
    end subroutine

    subroutine acc_delete_64_h (a, len)
      use iso_c_binding, only: c_int64_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
    end subroutine

    subroutine acc_delete_array_h (a)
      type (*), dimension (..), contiguous :: a
    end subroutine

    subroutine acc_delete_finalize_32_h (a, len)
      use iso_c_binding, only: c_int32_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
    end subroutine

    subroutine acc_delete_finalize_64_h (a, len)
      use iso_c_binding, only: c_int64_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
    end subroutine

    subroutine acc_delete_finalize_array_h (a)
      type (*), dimension (..), contiguous :: a
    end subroutine

    subroutine acc_update_device_32_h (a, len)
      use iso_c_binding, only: c_int32_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
    end subroutine

    subroutine acc_update_device_64_h (a, len)
      use iso_c_binding, only: c_int64_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
    end subroutine

    subroutine acc_update_device_array_h (a)
      type (*), dimension (..), contiguous :: a
    end subroutine

    subroutine acc_update_self_32_h (a, len)
      use iso_c_binding, only: c_int32_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
    end subroutine

    subroutine acc_update_self_64_h (a, len)
      use iso_c_binding, only: c_int64_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
    end subroutine

    subroutine acc_update_self_array_h (a)
      type (*), dimension (..), contiguous :: a
    end subroutine

    function acc_is_present_32_h (a, len)
      use iso_c_binding, only: c_int32_t
      logical acc_is_present_32_h
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
    end function

    function acc_is_present_64_h (a, len)
      use iso_c_binding, only: c_int64_t
      logical acc_is_present_64_h
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
    end function

    function acc_is_present_array_h (a)
      logical acc_is_present_array_h
      type (*), dimension (..), contiguous :: a
    end function

    subroutine acc_copyin_async_32_h (a, len, async)
      use iso_c_binding, only: c_int32_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_copyin_async_64_h (a, len, async)
      use iso_c_binding, only: c_int64_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_copyin_async_array_h (a, async)
      use openacc_kinds, only: acc_handle_kind
      type (*), dimension (..), contiguous :: a
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_create_async_32_h (a, len, async)
      use iso_c_binding, only: c_int32_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_create_async_64_h (a, len, async)
      use iso_c_binding, only: c_int64_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_create_async_array_h (a, async)
      use openacc_kinds, only: acc_handle_kind
      type (*), dimension (..), contiguous :: a
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_copyout_async_32_h (a, len, async)
      use iso_c_binding, only: c_int32_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_copyout_async_64_h (a, len, async)
      use iso_c_binding, only: c_int64_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_copyout_async_array_h (a, async)
      use openacc_kinds, only: acc_handle_kind
      type (*), dimension (..), contiguous :: a
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_delete_async_32_h (a, len, async)
      use iso_c_binding, only: c_int32_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_delete_async_64_h (a, len, async)
      use iso_c_binding, only: c_int64_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_delete_async_array_h (a, async)
      use openacc_kinds, only: acc_handle_kind
      type (*), dimension (..), contiguous :: a
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_update_device_async_32_h (a, len, async)
      use iso_c_binding, only: c_int32_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_update_device_async_64_h (a, len, async)
      use iso_c_binding, only: c_int64_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_update_device_async_array_h (a, async)
      use openacc_kinds, only: acc_handle_kind
      type (*), dimension (..), contiguous :: a
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_update_self_async_32_h (a, len, async)
      use iso_c_binding, only: c_int32_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int32_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_update_self_async_64_h (a, len, async)
      use iso_c_binding, only: c_int64_t
      use openacc_kinds, only: acc_handle_kind
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_int64_t) len
      integer (acc_handle_kind) async
    end subroutine

    subroutine acc_update_self_async_array_h (a, async)
      use openacc_kinds, only: acc_handle_kind
      type (*), dimension (..), contiguous :: a
      integer (acc_handle_kind) async
    end subroutine
  end interface

  interface
    function acc_get_num_devices_l (devicetype) &
        bind (C, name = "acc_get_num_devices")
      use iso_c_binding, only: c_int
      integer (c_int) :: acc_get_num_devices_l
      integer (c_int), value :: devicetype
    end function

    subroutine acc_set_device_type_l (devicetype) &
        bind (C, name = "acc_set_device_type")
      use iso_c_binding, only: c_int
      integer (c_int), value :: devicetype
    end subroutine

    function acc_get_device_type_l () &
        bind (C, name = "acc_get_device_type")
      use iso_c_binding, only: c_int
      integer (c_int) :: acc_get_device_type_l
    end function

    subroutine acc_set_device_num_l (devicenum, devicetype) &
        bind (C, name = "acc_set_device_num")
      use iso_c_binding, only: c_int
      integer (c_int), value :: devicenum, devicetype
    end subroutine

    function acc_get_device_num_l (devicetype) &
        bind (C, name = "acc_get_device_num")
      use iso_c_binding, only: c_int
      integer (c_int) :: acc_get_device_num_l
      integer (c_int), value :: devicetype
    end function

    function acc_get_property_l (devicenum, devicetype, property) &
        bind (C, name = "acc_get_property")
      use iso_c_binding, only: c_int, c_size_t
      implicit none (type, external)
      integer (c_size_t) :: acc_get_property_l
      integer (c_int), value :: devicenum
      integer (c_int), value :: devicetype
      integer (c_int), value :: property
    end function

    function acc_get_property_string_l (devicenum, devicetype, property) &
        bind (C, name = "acc_get_property_string")
      use iso_c_binding, only: c_int, c_ptr
      implicit none (type, external)
      type (c_ptr) :: acc_get_property_string_l
      integer (c_int), value :: devicenum
      integer (c_int), value :: devicetype
      integer (c_int), value :: property
    end function

    function acc_async_test_l (a) &
        bind (C, name = "acc_async_test")
      use iso_c_binding, only: c_int
      integer (c_int) :: acc_async_test_l
      integer (c_int), value :: a
    end function

    function acc_async_test_all_l () &
        bind (C, name = "acc_async_test_all")
      use iso_c_binding, only: c_int
      integer (c_int) :: acc_async_test_all_l
    end function

    subroutine acc_wait_l (a) &
        bind (C, name = "acc_wait")
      use iso_c_binding, only: c_int
      integer (c_int), value :: a
    end subroutine

    subroutine acc_wait_async_l (arg, async) &
        bind (C, name = "acc_wait_async")
      use iso_c_binding, only: c_int
      integer (c_int), value :: arg, async
    end subroutine

    subroutine acc_wait_all_l () &
        bind (C, name = "acc_wait_all")
      use iso_c_binding, only: c_int
    end subroutine

    subroutine acc_wait_all_async_l (async) &
        bind (C, name = "acc_wait_all_async")
      use iso_c_binding, only: c_int
      integer (c_int), value :: async
    end subroutine

    subroutine acc_init_l (devicetype) &
        bind (C, name = "acc_init")
      use iso_c_binding, only: c_int
      integer (c_int), value :: devicetype
    end subroutine

    subroutine acc_shutdown_l (devicetype) &
        bind (C, name = "acc_shutdown")
      use iso_c_binding, only: c_int
      integer (c_int), value :: devicetype
    end subroutine

    function acc_on_device_l (devicetype) &
        bind (C, name = "acc_on_device")
      use iso_c_binding, only: c_int
      integer (c_int) :: acc_on_device_l
      integer (c_int), value :: devicetype
    end function

    subroutine acc_copyin_l (a, len) &
        bind (C, name = "acc_copyin")
      use iso_c_binding, only: c_size_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
    end subroutine

    subroutine acc_present_or_copyin_l (a, len) &
        bind (C, name = "acc_present_or_copyin")
      use iso_c_binding, only: c_size_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
    end subroutine

    subroutine acc_create_l (a, len) &
        bind (C, name = "acc_create")
      use iso_c_binding, only: c_size_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
    end subroutine

    subroutine acc_present_or_create_l (a, len) &
        bind (C, name = "acc_present_or_create")
      use iso_c_binding, only: c_size_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
    end subroutine

    subroutine acc_copyout_l (a, len) &
        bind (C, name = "acc_copyout")
      use iso_c_binding, only: c_size_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
    end subroutine

    subroutine acc_copyout_finalize_l (a, len) &
        bind (C, name = "acc_copyout_finalize")
      use iso_c_binding, only: c_size_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
    end subroutine

    subroutine acc_delete_l (a, len) &
        bind (C, name = "acc_delete")
      use iso_c_binding, only: c_size_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
    end subroutine

    subroutine acc_delete_finalize_l (a, len) &
        bind (C, name = "acc_delete_finalize")
      use iso_c_binding, only: c_size_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
    end subroutine

    subroutine acc_update_device_l (a, len) &
        bind (C, name = "acc_update_device")
      use iso_c_binding, only: c_size_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
    end subroutine

    subroutine acc_update_self_l (a, len) &
        bind (C, name = "acc_update_self")
      use iso_c_binding, only: c_size_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
    end subroutine

    function acc_is_present_l (a, len) &
        bind (C, name = "acc_is_present")
      use iso_c_binding, only: c_int32_t, c_size_t
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      integer (c_int32_t) :: acc_is_present_l
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
    end function

    subroutine acc_copyin_async_l (a, len, async) &
        bind (C, name = "acc_copyin_async")
      use iso_c_binding, only: c_size_t, c_int
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
      integer (c_int), value :: async
    end subroutine

    subroutine acc_create_async_l (a, len, async) &
        bind (C, name = "acc_create_async")
      use iso_c_binding, only: c_size_t, c_int
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
      integer (c_int), value :: async
    end subroutine

    subroutine acc_copyout_async_l (a, len, async) &
        bind (C, name = "acc_copyout_async")
      use iso_c_binding, only: c_size_t, c_int
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
      integer (c_int), value :: async
    end subroutine

    subroutine acc_delete_async_l (a, len, async) &
        bind (C, name = "acc_delete_async")
      use iso_c_binding, only: c_size_t, c_int
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
      integer (c_int), value :: async
    end subroutine

    subroutine acc_update_device_async_l (a, len, async) &
        bind (C, name = "acc_update_device_async")
      use iso_c_binding, only: c_size_t, c_int
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
      integer (c_int), value :: async
    end subroutine

    subroutine acc_update_self_async_l (a, len, async) &
        bind (C, name = "acc_update_self_async")
      use iso_c_binding, only: c_size_t, c_int
      !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
      type (*), dimension (*) :: a
      integer (c_size_t), value :: len
      integer (c_int), value :: async
    end subroutine
  end interface
end module openacc_internal

module openacc
  use openacc_kinds
  use openacc_internal
  implicit none

  private

  ! From openacc_kinds
  public :: acc_device_kind
  public :: acc_device_none, acc_device_default, acc_device_host
  public :: acc_device_not_host, acc_device_nvidia, acc_device_radeon

  public :: acc_device_property_kind, acc_device_property
  public :: acc_property_memory, acc_property_free_memory
  public :: acc_property_name, acc_property_vendor, acc_property_driver

  public :: acc_handle_kind
  public :: acc_async_noval, acc_async_sync

  public :: openacc_version

  public :: acc_get_num_devices, acc_set_device_type, acc_get_device_type
  public :: acc_set_device_num, acc_get_device_num
  public :: acc_get_property, acc_get_property_string
  public :: acc_async_test, acc_async_test_all
  public :: acc_wait, acc_async_wait, acc_wait_async
  public :: acc_wait_all, acc_async_wait_all, acc_wait_all_async
  public :: acc_init, acc_shutdown, acc_on_device
  public :: acc_copyin, acc_present_or_copyin, acc_pcopyin, acc_create
  public :: acc_present_or_create, acc_pcreate, acc_copyout, acc_delete
  public :: acc_update_device, acc_update_self, acc_is_present
  public :: acc_copyin_async, acc_create_async, acc_copyout_async
  public :: acc_delete_async, acc_update_device_async, acc_update_self_async
  public :: acc_copyout_finalize, acc_delete_finalize

  integer, parameter :: openacc_version = 201711

  interface acc_get_num_devices
    procedure :: acc_get_num_devices_h
  end interface

  interface acc_set_device_type
    procedure :: acc_set_device_type_h
  end interface

  interface acc_get_device_type
    procedure :: acc_get_device_type_h
  end interface

  interface acc_set_device_num
    procedure :: acc_set_device_num_h
  end interface

  interface acc_get_device_num
    procedure :: acc_get_device_num_h
  end interface

  interface acc_get_property
    procedure :: acc_get_property_h
  end interface

  interface acc_get_property_string
    procedure :: acc_get_property_string_h
  end interface

  interface acc_async_test
    procedure :: acc_async_test_h
  end interface

  interface acc_async_test_all
    procedure :: acc_async_test_all_h
  end interface

  interface acc_wait
    procedure :: acc_wait_h
  end interface

  ! acc_async_wait is an OpenACC 1.0 compatibility name for acc_wait.
  interface acc_async_wait
    procedure :: acc_wait_h
  end interface

  interface acc_wait_async
    procedure :: acc_wait_async_h
  end interface

  interface acc_wait_all
    procedure :: acc_wait_all_h
  end interface

  ! acc_async_wait_all is an OpenACC 1.0 compatibility name for acc_wait_all.
  interface acc_async_wait_all
    procedure :: acc_wait_all_h
  end interface

  interface acc_wait_all_async
    procedure :: acc_wait_all_async_h
  end interface

  interface acc_init
    procedure :: acc_init_h
  end interface

  interface acc_shutdown
    procedure :: acc_shutdown_h
  end interface

  interface acc_on_device
    procedure :: acc_on_device_h
  end interface

  ! acc_malloc: Only available in C/C++
  ! acc_free: Only available in C/C++

  ! As vendor extension, the following code supports both 32bit and 64bit
  ! arguments for "size"; the OpenACC standard only permits default-kind
  ! integers, which are of kind 4 (i.e. 32 bits).
  ! Additionally, the two-argument version also takes arrays as argument.
  ! and the one argument version also scalars. Note that the code assumes
  ! that the arrays are contiguous.

  interface acc_copyin
    procedure :: acc_copyin_32_h
    procedure :: acc_copyin_64_h
    procedure :: acc_copyin_array_h
  end interface

  interface acc_present_or_copyin
    procedure :: acc_present_or_copyin_32_h
    procedure :: acc_present_or_copyin_64_h
    procedure :: acc_present_or_copyin_array_h
  end interface

  interface acc_pcopyin
    procedure :: acc_present_or_copyin_32_h
    procedure :: acc_present_or_copyin_64_h
    procedure :: acc_present_or_copyin_array_h
  end interface

  interface acc_create
    procedure :: acc_create_32_h
    procedure :: acc_create_64_h
    procedure :: acc_create_array_h
  end interface

  interface acc_present_or_create
    procedure :: acc_present_or_create_32_h
    procedure :: acc_present_or_create_64_h
    procedure :: acc_present_or_create_array_h
  end interface

  interface acc_pcreate
    procedure :: acc_present_or_create_32_h
    procedure :: acc_present_or_create_64_h
    procedure :: acc_present_or_create_array_h
  end interface

  interface acc_copyout
    procedure :: acc_copyout_32_h
    procedure :: acc_copyout_64_h
    procedure :: acc_copyout_array_h
  end interface

  interface acc_copyout_finalize
    procedure :: acc_copyout_finalize_32_h
    procedure :: acc_copyout_finalize_64_h
    procedure :: acc_copyout_finalize_array_h
  end interface

  interface acc_delete
    procedure :: acc_delete_32_h
    procedure :: acc_delete_64_h
    procedure :: acc_delete_array_h
  end interface

  interface acc_delete_finalize
    procedure :: acc_delete_finalize_32_h
    procedure :: acc_delete_finalize_64_h
    procedure :: acc_delete_finalize_array_h
  end interface

  interface acc_update_device
    procedure :: acc_update_device_32_h
    procedure :: acc_update_device_64_h
    procedure :: acc_update_device_array_h
  end interface

  interface acc_update_self
    procedure :: acc_update_self_32_h
    procedure :: acc_update_self_64_h
    procedure :: acc_update_self_array_h
  end interface

  ! acc_map_data: Only available in C/C++
  ! acc_unmap_data: Only available in C/C++
  ! acc_deviceptr: Only available in C/C++
  ! acc_hostptr: Only available in C/C++

  interface acc_is_present
    procedure :: acc_is_present_32_h
    procedure :: acc_is_present_64_h
    procedure :: acc_is_present_array_h
  end interface

  ! acc_memcpy_to_device: Only available in C/C++
  ! acc_memcpy_from_device: Only available in C/C++

  interface acc_copyin_async
    procedure :: acc_copyin_async_32_h
    procedure :: acc_copyin_async_64_h
    procedure :: acc_copyin_async_array_h
  end interface

  interface acc_create_async
    procedure :: acc_create_async_32_h
    procedure :: acc_create_async_64_h
    procedure :: acc_create_async_array_h
  end interface

  interface acc_copyout_async
    procedure :: acc_copyout_async_32_h
    procedure :: acc_copyout_async_64_h
    procedure :: acc_copyout_async_array_h
  end interface

  interface acc_delete_async
    procedure :: acc_delete_async_32_h
    procedure :: acc_delete_async_64_h
    procedure :: acc_delete_async_array_h
  end interface

  interface acc_update_device_async
    procedure :: acc_update_device_async_32_h
    procedure :: acc_update_device_async_64_h
    procedure :: acc_update_device_async_array_h
  end interface

  interface acc_update_self_async
    procedure :: acc_update_self_async_32_h
    procedure :: acc_update_self_async_64_h
    procedure :: acc_update_self_async_array_h
  end interface

end module openacc

function acc_get_num_devices_h (devicetype)
  use openacc_internal, only: acc_get_num_devices_l
  use openacc_kinds
  integer acc_get_num_devices_h
  integer (acc_device_kind) devicetype
  acc_get_num_devices_h = acc_get_num_devices_l (devicetype)
end function

subroutine acc_set_device_type_h (devicetype)
  use openacc_internal, only: acc_set_device_type_l
  use openacc_kinds
  integer (acc_device_kind) devicetype
  call acc_set_device_type_l (devicetype)
end subroutine

function acc_get_device_type_h ()
  use openacc_internal, only: acc_get_device_type_l
  use openacc_kinds
  integer (acc_device_kind) acc_get_device_type_h
  acc_get_device_type_h = acc_get_device_type_l ()
end function

subroutine acc_set_device_num_h (devicenum, devicetype)
  use openacc_internal, only: acc_set_device_num_l
  use openacc_kinds
  integer devicenum
  integer (acc_device_kind) devicetype
  call acc_set_device_num_l (devicenum, devicetype)
end subroutine

function acc_get_device_num_h (devicetype)
  use openacc_internal, only: acc_get_device_num_l
  use openacc_kinds
  integer acc_get_device_num_h
  integer (acc_device_kind) devicetype
  acc_get_device_num_h = acc_get_device_num_l (devicetype)
end function

function acc_get_property_h (devicenum, devicetype, property)
  use iso_c_binding, only: c_size_t
  use openacc_internal, only: acc_get_property_l
  use openacc_kinds
  implicit none (type, external)
  integer (c_size_t) :: acc_get_property_h
  integer, value :: devicenum
  integer (acc_device_kind), value :: devicetype
  integer (acc_device_property_kind), value :: property
  acc_get_property_h = acc_get_property_l (devicenum, devicetype, property)
end function

subroutine acc_get_property_string_h (devicenum, devicetype, property, string)
  use iso_c_binding, only: c_char, c_size_t, c_ptr, c_f_pointer, c_associated
  use openacc_internal, only: acc_get_property_string_l
  use openacc_kinds
  implicit none (type, external)
  integer, value :: devicenum
  integer (acc_device_kind), value :: devicetype
  integer (acc_device_property_kind), value :: property
  character (*) :: string

  type (c_ptr) :: cptr
  integer(c_size_t) :: clen, slen, i
  character (kind=c_char, len=1), pointer, contiguous :: sptr (:)

  interface
     function strlen (s) bind (C, name = "strlen")
       use iso_c_binding, only: c_ptr, c_size_t
       type (c_ptr), intent(in), value :: s
       integer (c_size_t) :: strlen
     end function strlen
  end interface

  cptr = acc_get_property_string_l (devicenum, devicetype, property)
  string = ""
  if (.not. c_associated (cptr)) then
     return
  end if

  clen = strlen (cptr)
  call c_f_pointer (cptr, sptr, [clen])

  slen = min (clen, len (string, kind=c_size_t))
  do i = 1, slen
    string (i:i) = sptr (i)
  end do
end subroutine

function acc_async_test_h (arg)
  use openacc_internal, only: acc_async_test_l
  logical acc_async_test_h
  integer arg
  acc_async_test_h = acc_async_test_l (arg) /= 0
end function

function acc_async_test_all_h ()
  use openacc_internal, only: acc_async_test_all_l
  logical acc_async_test_all_h
  acc_async_test_all_h = acc_async_test_all_l () /= 0
end function

subroutine acc_wait_h (arg)
  use openacc_internal, only: acc_wait_l
  integer arg
  call acc_wait_l (arg)
end subroutine

subroutine acc_wait_async_h (arg, async)
  use openacc_internal, only: acc_wait_async_l
  integer arg, async
  call acc_wait_async_l (arg, async)
end subroutine

subroutine acc_wait_all_h ()
  use openacc_internal, only: acc_wait_all_l
  call acc_wait_all_l ()
end subroutine

subroutine acc_wait_all_async_h (async)
  use openacc_internal, only: acc_wait_all_async_l
  integer async
  call acc_wait_all_async_l (async)
end subroutine

subroutine acc_init_h (devicetype)
  use openacc_internal, only: acc_init_l
  use openacc_kinds
  integer (acc_device_kind) devicetype
  call acc_init_l (devicetype)
end subroutine

subroutine acc_shutdown_h (devicetype)
  use openacc_internal, only: acc_shutdown_l
  use openacc_kinds
  integer (acc_device_kind) devicetype
  call acc_shutdown_l (devicetype)
end subroutine

function acc_on_device_h (devicetype)
  use openacc_internal, only: acc_on_device_l
  use openacc_kinds
  integer (acc_device_kind) devicetype
  logical acc_on_device_h
  acc_on_device_h = acc_on_device_l (devicetype) /= 0
end function

subroutine acc_copyin_32_h (a, len)
  use iso_c_binding, only: c_int32_t, c_size_t
  use openacc_internal, only: acc_copyin_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  call acc_copyin_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_copyin_64_h (a, len)
  use iso_c_binding, only: c_int64_t, c_size_t
  use openacc_internal, only: acc_copyin_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  call acc_copyin_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_copyin_array_h (a)
  use openacc_internal, only: acc_copyin_l
  type (*), dimension (..), contiguous :: a
  call acc_copyin_l (a, sizeof (a))
end subroutine

subroutine acc_present_or_copyin_32_h (a, len)
  use iso_c_binding, only: c_int32_t, c_size_t
  use openacc_internal, only: acc_present_or_copyin_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  call acc_present_or_copyin_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_present_or_copyin_64_h (a, len)
  use iso_c_binding, only: c_int64_t, c_size_t
  use openacc_internal, only: acc_present_or_copyin_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  call acc_present_or_copyin_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_present_or_copyin_array_h (a)
  use openacc_internal, only: acc_present_or_copyin_l
  type (*), dimension (..), contiguous :: a
  call acc_present_or_copyin_l (a, sizeof (a))
end subroutine

subroutine acc_create_32_h (a, len)
  use iso_c_binding, only: c_int32_t, c_size_t
  use openacc_internal, only: acc_create_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  call acc_create_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_create_64_h (a, len)
  use iso_c_binding, only: c_int64_t, c_size_t
  use openacc_internal, only: acc_create_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  call acc_create_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_create_array_h (a)
  use openacc_internal, only: acc_create_l
  type (*), dimension (..), contiguous :: a
  call acc_create_l (a, sizeof (a))
end subroutine

subroutine acc_present_or_create_32_h (a, len)
  use iso_c_binding, only: c_int32_t, c_size_t
  use openacc_internal, only: acc_present_or_create_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  call acc_present_or_create_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_present_or_create_64_h (a, len)
  use iso_c_binding, only: c_int64_t, c_size_t
  use openacc_internal, only: acc_present_or_create_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  call acc_present_or_create_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_present_or_create_array_h (a)
  use openacc_internal, only: acc_present_or_create_l
  type (*), dimension (..), contiguous :: a
  call acc_present_or_create_l (a, sizeof (a))
end subroutine

subroutine acc_copyout_32_h (a, len)
  use iso_c_binding, only: c_int32_t, c_size_t
  use openacc_internal, only: acc_copyout_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  call acc_copyout_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_copyout_64_h (a, len)
  use iso_c_binding, only: c_int64_t, c_size_t
  use openacc_internal, only: acc_copyout_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  call acc_copyout_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_copyout_array_h (a)
  use openacc_internal, only: acc_copyout_l
  type (*), dimension (..), contiguous :: a
  call acc_copyout_l (a, sizeof (a))
end subroutine

subroutine acc_copyout_finalize_32_h (a, len)
  use iso_c_binding, only: c_int32_t, c_size_t
  use openacc_internal, only: acc_copyout_finalize_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  call acc_copyout_finalize_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_copyout_finalize_64_h (a, len)
  use iso_c_binding, only: c_int64_t, c_size_t
  use openacc_internal, only: acc_copyout_finalize_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  call acc_copyout_finalize_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_copyout_finalize_array_h (a)
  use openacc_internal, only: acc_copyout_finalize_l
  type (*), dimension (..), contiguous :: a
  call acc_copyout_finalize_l (a, sizeof (a))
end subroutine

subroutine acc_delete_32_h (a, len)
  use iso_c_binding, only: c_int32_t, c_size_t
  use openacc_internal, only: acc_delete_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  call acc_delete_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_delete_64_h (a, len)
  use iso_c_binding, only: c_int64_t, c_size_t
  use openacc_internal, only: acc_delete_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  call acc_delete_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_delete_array_h (a)
  use openacc_internal, only: acc_delete_l
  type (*), dimension (..), contiguous :: a
  call acc_delete_l (a, sizeof (a))
end subroutine

subroutine acc_delete_finalize_32_h (a, len)
  use iso_c_binding, only: c_int32_t, c_size_t
  use openacc_internal, only: acc_delete_finalize_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  call acc_delete_finalize_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_delete_finalize_64_h (a, len)
  use iso_c_binding, only: c_int64_t, c_size_t
  use openacc_internal, only: acc_delete_finalize_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  call acc_delete_finalize_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_delete_finalize_array_h (a)
  use openacc_internal, only: acc_delete_finalize_l
  type (*), dimension (..), contiguous :: a
  call acc_delete_finalize_l (a, sizeof (a))
end subroutine

subroutine acc_update_device_32_h (a, len)
  use iso_c_binding, only: c_int32_t, c_size_t
  use openacc_internal, only: acc_update_device_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  call acc_update_device_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_update_device_64_h (a, len)
  use iso_c_binding, only: c_int64_t, c_size_t
  use openacc_internal, only: acc_update_device_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  call acc_update_device_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_update_device_array_h (a)
  use openacc_internal, only: acc_update_device_l
  type (*), dimension (..), contiguous :: a
  call acc_update_device_l (a, sizeof (a))
end subroutine

subroutine acc_update_self_32_h (a, len)
  use iso_c_binding, only: c_int32_t, c_size_t
  use openacc_internal, only: acc_update_self_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  call acc_update_self_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_update_self_64_h (a, len)
  use iso_c_binding, only: c_int64_t, c_size_t
  use openacc_internal, only: acc_update_self_l
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  call acc_update_self_l (a, int (len, kind = c_size_t))
end subroutine

subroutine acc_update_self_array_h (a)
  use openacc_internal, only: acc_update_self_l
  type (*), dimension (..), contiguous :: a
  call acc_update_self_l (a, sizeof (a))
end subroutine

function acc_is_present_32_h (a, len)
  use iso_c_binding, only: c_int32_t, c_size_t
  use openacc_internal, only: acc_is_present_l
  logical acc_is_present_32_h
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  acc_is_present_32_h = acc_is_present_l (a, int (len, kind = c_size_t)) /= 0
end function

function acc_is_present_64_h (a, len)
  use iso_c_binding, only: c_int64_t, c_size_t
  use openacc_internal, only: acc_is_present_l
  logical acc_is_present_64_h
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  acc_is_present_64_h = acc_is_present_l (a, int (len, kind = c_size_t)) /= 0
end function

function acc_is_present_array_h (a)
  use openacc_internal, only: acc_is_present_l
  logical acc_is_present_array_h
  type (*), dimension (..), contiguous :: a
  acc_is_present_array_h = acc_is_present_l (a, sizeof (a)) /= 0
end function

subroutine acc_copyin_async_32_h (a, len, async)
  use iso_c_binding, only: c_int32_t, c_size_t, c_int
  use openacc_internal, only: acc_copyin_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  integer (acc_handle_kind) async
  call acc_copyin_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_copyin_async_64_h (a, len, async)
  use iso_c_binding, only: c_int64_t, c_size_t, c_int
  use openacc_internal, only: acc_copyin_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  integer (acc_handle_kind) async
  call acc_copyin_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_copyin_async_array_h (a, async)
  use iso_c_binding, only: c_int
  use openacc_internal, only: acc_copyin_async_l
  use openacc_kinds, only: acc_handle_kind
  type (*), dimension (..), contiguous :: a
  integer (acc_handle_kind) async
  call acc_copyin_async_l (a, sizeof (a), int (async, kind = c_int))
end subroutine

subroutine acc_create_async_32_h (a, len, async)
  use iso_c_binding, only: c_int32_t, c_size_t, c_int
  use openacc_internal, only: acc_create_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  integer (acc_handle_kind) async
  call acc_create_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_create_async_64_h (a, len, async)
  use iso_c_binding, only: c_int64_t, c_size_t, c_int
  use openacc_internal, only: acc_create_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  integer (acc_handle_kind) async
  call acc_create_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_create_async_array_h (a, async)
  use iso_c_binding, only: c_int
  use openacc_internal, only: acc_create_async_l
  use openacc_kinds, only: acc_handle_kind
  type (*), dimension (..), contiguous :: a
  integer (acc_handle_kind) async
  call acc_create_async_l (a, sizeof (a), int (async, kind = c_int))
end subroutine

subroutine acc_copyout_async_32_h (a, len, async)
  use iso_c_binding, only: c_int32_t, c_size_t, c_int
  use openacc_internal, only: acc_copyout_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  integer (acc_handle_kind) async
  call acc_copyout_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_copyout_async_64_h (a, len, async)
  use iso_c_binding, only: c_int64_t, c_size_t, c_int
  use openacc_internal, only: acc_copyout_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  integer (acc_handle_kind) async
  call acc_copyout_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_copyout_async_array_h (a, async)
  use iso_c_binding, only: c_int
  use openacc_internal, only: acc_copyout_async_l
  use openacc_kinds, only: acc_handle_kind
  type (*), dimension (..), contiguous :: a
  integer (acc_handle_kind) async
  call acc_copyout_async_l (a, sizeof (a), int (async, kind = c_int))
end subroutine

subroutine acc_delete_async_32_h (a, len, async)
  use iso_c_binding, only: c_int32_t, c_size_t, c_int
  use openacc_internal, only: acc_delete_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  integer (acc_handle_kind) async
  call acc_delete_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_delete_async_64_h (a, len, async)
  use iso_c_binding, only: c_int64_t, c_size_t, c_int
  use openacc_internal, only: acc_delete_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  integer (acc_handle_kind) async
  call acc_delete_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_delete_async_array_h (a, async)
  use iso_c_binding, only: c_int
  use openacc_internal, only: acc_delete_async_l
  use openacc_kinds, only: acc_handle_kind
  type (*), dimension (..), contiguous :: a
  integer (acc_handle_kind) async
  call acc_delete_async_l (a, sizeof (a), int (async, kind = c_int))
end subroutine

subroutine acc_update_device_async_32_h (a, len, async)
  use iso_c_binding, only: c_int32_t, c_size_t, c_int
  use openacc_internal, only: acc_update_device_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  integer (acc_handle_kind) async
  call acc_update_device_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_update_device_async_64_h (a, len, async)
  use iso_c_binding, only: c_int64_t, c_size_t, c_int
  use openacc_internal, only: acc_update_device_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  integer (acc_handle_kind) async
  call acc_update_device_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_update_device_async_array_h (a, async)
  use iso_c_binding, only: c_int
  use openacc_internal, only: acc_update_device_async_l
  use openacc_kinds, only: acc_handle_kind
  type (*), dimension (..), contiguous :: a
  integer (acc_handle_kind) async
  call acc_update_device_async_l (a, sizeof (a), int (async, kind = c_int))
end subroutine

subroutine acc_update_self_async_32_h (a, len, async)
  use iso_c_binding, only: c_int32_t, c_size_t, c_int
  use openacc_internal, only: acc_update_self_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int32_t) len
  integer (acc_handle_kind) async
  call acc_update_self_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_update_self_async_64_h (a, len, async)
  use iso_c_binding, only: c_int64_t, c_size_t, c_int
  use openacc_internal, only: acc_update_self_async_l
  use openacc_kinds, only: acc_handle_kind
  !GCC$ ATTRIBUTES NO_ARG_CHECK :: a
  type (*), dimension (*) :: a
  integer (c_int64_t) len
  integer (acc_handle_kind) async
  call acc_update_self_async_l (a, int (len, kind = c_size_t), int (async, kind = c_int))
end subroutine

subroutine acc_update_self_async_array_h (a, async)
  use iso_c_binding, only: c_int
  use openacc_internal, only: acc_update_self_async_l
  use openacc_kinds, only: acc_handle_kind
  type (*), dimension (..), contiguous :: a
  integer (acc_handle_kind) async
  call acc_update_self_async_l (a, sizeof (a), int (async, kind = c_int))
end subroutine
