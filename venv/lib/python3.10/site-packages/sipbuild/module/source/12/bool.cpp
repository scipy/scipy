// This contains all the C++ code that is needed by the sip module.
//
// Copyright (c) 2015 Riverbank Computing Limited <info@riverbankcomputing.com>
//
// This file is part of SIP.
//
// This copy of SIP is licensed for use under the terms of the SIP License
// Agreement.  See the file LICENSE for more details.
//
// This copy of SIP may also used under the terms of the GNU General Public
// License v2 or v3 as published by the Free Software Foundation which can be
// found in the files LICENSE-GPL2 and LICENSE-GPL3 included in this package.
//
// SIP is supplied WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


// Set a C++ bool for the main C implementation of the module.
extern "C" void sipSetBool(void *ptr, int val)
{
    *reinterpret_cast<bool *>(ptr) = !!val;
}
