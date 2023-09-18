#
# Copyright (C) 2015 Ford Motor Company
# Contact: https://www.qt.io/licensing/
#
# This file is part of the QtRemoteObjects module of the Qt Toolkit.
#
# $QT_BEGIN_LICENSE:LGPL$
# Commercial License Usage
# Licensees holding valid commercial Qt licenses may use this file in
# accordance with the commercial license agreement provided with the
# Software or, alternatively, in accordance with the terms contained in
# a written agreement between you and The Qt Company. For licensing terms
# and conditions see https://www.qt.io/terms-conditions. For further
# information use the contact form at https://www.qt.io/contact-us.
#
# GNU Lesser General Public License Usage
# Alternatively, this file may be used under the terms of the GNU Lesser
# General Public License version 3 as published by the Free Software
# Foundation and appearing in the file LICENSE.LGPL3 included in the
# packaging of this file. Please review the following information to
# ensure the GNU Lesser General Public License version 3 requirements
# will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
#
# GNU General Public License Usage
# Alternatively, this file may be used under the terms of the GNU
# General Public License version 2.0 or (at your option) the GNU General
# Public license version 3 or any later version approved by the KDE Free
# Qt Foundation. The licenses are as published by the Free Software
# Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
# included in the packaging of this file. Please review the following
# information to ensure the GNU General Public License requirements will
# be met: https://www.gnu.org/licenses/gpl-2.0.html and
# https://www.gnu.org/licenses/gpl-3.0.html.
#
# $QT_END_LICENSE$

if (NOT TARGET Qt5::repc)
    add_executable(Qt5::repc IMPORTED)

    set(imported_location "${_qt5RemoteObjects_install_prefix}/bin/repc")
    _qt5_RemoteObjects_check_file_exists(${imported_location})

    set_target_properties(Qt5::repc PROPERTIES
        IMPORTED_LOCATION ${imported_location}
    )
    get_target_property(Qt5RemoteObjects_REPC_EXECUTABLE Qt5::repc LOCATION)
endif()

# Create versionless tool targets.
foreach(__qt_tool repc)
    if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND NOT TARGET Qt::${__qt_tool}
       AND TARGET Qt5::${__qt_tool})
        add_executable(Qt::${__qt_tool} IMPORTED)
        get_target_property(__qt_imported_location Qt5::${__qt_tool} IMPORTED_LOCATION)
        set_target_properties(Qt::${__qt_tool}
                              PROPERTIES IMPORTED_LOCATION "${__qt_imported_location}")
    endif()
endforeach()
