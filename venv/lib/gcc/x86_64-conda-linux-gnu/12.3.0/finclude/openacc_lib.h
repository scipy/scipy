!  OpenACC Runtime Library Definitions.			-*- mode: fortran -*-

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

! NOTE: Due to the use of dimension (..), the code only works when compiled
! with -std=f2008ts/gnu/legacy but not with other standard settings.
! Alternatively, the user can use the module version, which permits
! compilation with -std=f95.

! Keep in sync with openacc.f90 and config/accel/openacc.f90.

      integer, parameter :: acc_device_kind = 4

!     Keep in sync with include/gomp-constants.h.
      integer (acc_device_kind), parameter :: acc_device_current = -1
      integer (acc_device_kind), parameter :: acc_device_none = 0
      integer (acc_device_kind), parameter :: acc_device_default = 1
      integer (acc_device_kind), parameter :: acc_device_host = 2
!     integer (acc_device_kind), parameter :: acc_device_host_nonshm = 3
!     removed.
      integer (acc_device_kind), parameter :: acc_device_not_host = 4
      integer (acc_device_kind), parameter :: acc_device_nvidia = 5
      integer (acc_device_kind), parameter :: acc_device_radeon = 8

      integer, parameter :: acc_device_property_kind = 4
! OpenACC 2.6/2.7/3.0 used acc_device_property; in a spec update the
! missing '_kind' was added for consistency.  For backward compatibility, keep:
      integer, parameter :: acc_device_property                         &
     &    = acc_device_property_kind

      integer (acc_device_property_kind), parameter ::                  &
     &    acc_property_memory = 1
      integer (acc_device_property_kind), parameter ::                  &
     &    acc_property_free_memory = 2
      integer (acc_device_property_kind), parameter ::                  &
     &    acc_property_name = int(Z'10001')
      integer (acc_device_property_kind), parameter ::                  &
     &    acc_property_vendor = int(Z'10002')
      integer (acc_device_property_kind), parameter ::                  &
     &    acc_property_driver = int(Z'10003')

      integer, parameter :: acc_handle_kind = 4

!     Keep in sync with include/gomp-constants.h.
      integer (acc_handle_kind), parameter :: acc_async_noval = -1
      integer (acc_handle_kind), parameter :: acc_async_sync = -2

      integer, parameter :: openacc_version = 201711

      interface acc_get_num_devices
        function acc_get_num_devices_h (devicetype)
          import acc_device_kind
          integer acc_get_num_devices_h
          integer (acc_device_kind) devicetype
        end function
      end interface

      interface acc_set_device_type
        subroutine acc_set_device_type_h (devicetype)
          import acc_device_kind
          integer (acc_device_kind) devicetype
        end subroutine
      end interface

      interface acc_get_device_type
        function acc_get_device_type_h ()
          import acc_device_kind
          integer (acc_device_kind) acc_get_device_type_h
        end function
      end interface

      interface acc_set_device_num
        subroutine acc_set_device_num_h (devicenum, devicetype)
          import acc_device_kind
          integer devicenum
          integer (acc_device_kind) devicetype
        end subroutine
      end interface

      interface acc_get_device_num
        function acc_get_device_num_h (devicetype)
          import acc_device_kind
          integer acc_get_device_num_h
          integer (acc_device_kind) devicetype
        end function
      end interface

      interface acc_get_property
        function acc_get_property_h (devicenum, devicetype,             &
     &                               property)
          use iso_c_binding, only: c_size_t
          import acc_device_kind, acc_device_property_kind
          implicit none (type, external)
          integer (c_size_t) :: acc_get_property_h
          integer, value :: devicenum
          integer (acc_device_kind), value :: devicetype
          integer (acc_device_property_kind), value :: property
        end function
      end interface

      interface acc_get_property_string
        subroutine acc_get_property_string_h (devicenum, devicetype,    &
     &                                        property, string)
          import acc_device_kind, acc_device_property_kind
          implicit none (type, external)
          integer, value :: devicenum
          integer (acc_device_kind), value :: devicetype
          integer (acc_device_property_kind), value :: property
          character (*) :: string
        end subroutine
      end interface

      interface acc_async_test
        function acc_async_test_h (arg)
          logical acc_async_test_h
          integer arg
        end function
      end interface

      interface acc_async_test_all
        function acc_async_test_all_h ()
          logical acc_async_test_all_h
        end function
      end interface

      interface acc_wait
        subroutine acc_wait_h (arg)
          integer arg
        end subroutine
      end interface

!     acc_async_wait is an OpenACC 1.0 compatibility name for acc_wait.
      interface acc_async_wait
        procedure :: acc_wait_h
      end interface

      interface acc_wait_async
        subroutine acc_wait_async_h (arg, async)
          integer arg, async
        end subroutine
      end interface

      interface acc_wait_all
        subroutine acc_wait_all_h ()
        end subroutine
      end interface

!     acc_async_wait_all is an OpenACC 1.0 compatibility name for
!     acc_wait_all.
      interface acc_async_wait_all
        procedure :: acc_wait_all_h
      end interface

      interface acc_wait_all_async
        subroutine acc_wait_all_async_h (async)
          integer async
        end subroutine
      end interface

      interface acc_init
        subroutine acc_init_h (devicetype)
          import acc_device_kind
          integer (acc_device_kind) devicetype
        end subroutine
      end interface

      interface acc_shutdown
        subroutine acc_shutdown_h (devicetype)
          import acc_device_kind
          integer (acc_device_kind) devicetype
        end subroutine
      end interface

      interface acc_on_device
        function acc_on_device_h (devicetype)
          import acc_device_kind
          logical acc_on_device_h
          integer (acc_device_kind) devicetype
        end function
      end interface

      ! acc_malloc: Only available in C/C++
      ! acc_free: Only available in C/C++

      interface acc_copyin
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
      end interface

      interface acc_present_or_copyin
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
      end interface

      interface acc_pcopyin
        procedure :: acc_present_or_copyin_32_h
        procedure :: acc_present_or_copyin_64_h
        procedure :: acc_present_or_copyin_array_h
      end interface

      interface acc_create
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
      end interface

      interface acc_present_or_create
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
      end interface

      interface acc_pcreate
        procedure :: acc_present_or_create_32_h
        procedure :: acc_present_or_create_64_h
        procedure :: acc_present_or_create_array_h
      end interface

      interface acc_copyout
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
      end interface

      interface acc_copyout_finalize
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
      end interface

      interface acc_delete
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
      end interface

      interface acc_delete_finalize
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
      end interface

      interface acc_update_device
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
      end interface

      interface acc_update_self
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
      end interface

      ! acc_map_data: Only available in C/C++
      ! acc_unmap_data: Only available in C/C++
      ! acc_deviceptr: Only available in C/C++
      ! acc_hostptr: Only available in C/C++

      interface acc_is_present
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
      end interface

      ! acc_memcpy_to_device: Only available in C/C++
      ! acc_memcpy_from_device: Only available in C/C++

      interface acc_copyin_async
        subroutine acc_copyin_async_32_h (a, len, async)
          use iso_c_binding, only: c_int32_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int32_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_copyin_async_64_h (a, len, async)
          use iso_c_binding, only: c_int64_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int64_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_copyin_async_array_h (a, async_)
          import acc_handle_kind
          type (*), dimension (..), contiguous :: a
          integer (acc_handle_kind) async_
        end subroutine
      end interface

      interface acc_create_async
        subroutine acc_create_async_32_h (a, len, async)
          use iso_c_binding, only: c_int32_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int32_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_create_async_64_h (a, len, async)
          use iso_c_binding, only: c_int64_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int64_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_create_async_array_h (a, async_)
          import acc_handle_kind
          type (*), dimension (..), contiguous :: a
          integer (acc_handle_kind) async_
        end subroutine
      end interface

      interface acc_copyout_async
        subroutine acc_copyout_async_32_h (a, len, async)
          use iso_c_binding, only: c_int32_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int32_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_copyout_async_64_h (a, len, async)
          use iso_c_binding, only: c_int64_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int64_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_copyout_async_array_h (a, async_)
          import acc_handle_kind
          type (*), dimension (..), contiguous :: a
          integer (acc_handle_kind) async_
        end subroutine
      end interface

      interface acc_delete_async
        subroutine acc_delete_async_32_h (a, len, async)
          use iso_c_binding, only: c_int32_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int32_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_delete_async_64_h (a, len, async)
          use iso_c_binding, only: c_int64_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int64_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_delete_async_array_h (a, async_)
          import acc_handle_kind
          type (*), dimension (..), contiguous :: a
          integer (acc_handle_kind) async_
        end subroutine
      end interface

      interface acc_update_device_async
        subroutine acc_update_device_async_32_h (a, len, async)
          use iso_c_binding, only: c_int32_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int32_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_update_device_async_64_h (a, len, async)
          use iso_c_binding, only: c_int64_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int64_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_update_device_async_array_h (a, async_)
          import acc_handle_kind
          type (*), dimension (..), contiguous :: a
          integer (acc_handle_kind) async_
        end subroutine
      end interface

      interface acc_update_self_async
        subroutine acc_update_self_async_32_h (a, len, async)
          use iso_c_binding, only: c_int32_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int32_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_update_self_async_64_h (a, len, async)
          use iso_c_binding, only: c_int64_t
          import acc_handle_kind
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
          type (*), dimension (*) :: a
          integer (c_int64_t) len
          integer (acc_handle_kind) async
        end subroutine

        subroutine acc_update_self_async_array_h (a, async_)
          import acc_handle_kind
          type (*), dimension (..), contiguous :: a
          integer (acc_handle_kind) async_
        end subroutine
      end interface
