/*
 *  This file is part of the MR utility library.
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/** \file ducc0/infra/system.h
 *  Helper functions for determining system resources
 *
 *  \note This functionality accesses the Linux /proc file system and will not
 *        work on other platforms.
 *  \copyright Copyright (C) 2019-2021 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef DUCC0_SYSTEM_H
#define DUCC0_SYSTEM_H

#include <string>
#include <cstddef>

namespace ducc0 {

namespace detail_system {

std::size_t getProcessInfo(const std::string &quantity);
std::size_t getMemInfo(const std::string &quantity);
std::size_t usable_memory();

}

using detail_system::getProcessInfo;
using detail_system::getMemInfo;
using detail_system::usable_memory;

}

#endif
